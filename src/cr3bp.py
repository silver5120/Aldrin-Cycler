"""
Earth-Moon Circular Restricted Three-Body Problem (CR3BP) Dynamics & Solvers.

All internal computations use canonical dimensionless units:
    - Length unit L* = mean Earth-Moon distance (384,400 km)
    - Mass parameter mu = mu_M / (mu_E + mu_M)
    - Time unit T* = 1 / n, so 1 sidereal lunar month = 2*pi T*
    - Velocity unit V* = L* / T*

In the barycentric rotating (synodic) frame:
    - Earth is located at (-mu, 0)
    - Moon is located at (1 - mu, 0)
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


# Physical constants (JPL / IAU standard values)
MU_EARTH_KM3_S2 = 398600.4418
MU_MOON_KM3_S2 = 4902.8001
L_STAR_KM = 384400.0
R_EARTH_KM = 6378.137
R_MOON_KM = 1737.4

MU = MU_MOON_KM3_S2 / (MU_EARTH_KM3_S2 + MU_MOON_KM3_S2)  # ~0.0121505856
N_RAD_S = np.sqrt((MU_EARTH_KM3_S2 + MU_MOON_KM3_S2) / (L_STAR_KM ** 3))
T_STAR_S = 1.0 / N_RAD_S
T_STAR_DAYS = T_STAR_S / 86400.0
V_STAR_KM_S = L_STAR_KM / T_STAR_S
LUNAR_MONTH_DAYS = 2.0 * np.pi * T_STAR_DAYS  # ~27.32166 days

R_EARTH_ND = R_EARTH_KM / L_STAR_KM
R_MOON_ND = R_MOON_KM / L_STAR_KM


def pseudo_potential(x: np.ndarray | float, y: np.ndarray | float, mu: float = MU):
    """Effective potential Omega(x, y) in the planar CR3BP rotating frame."""
    r1 = np.sqrt((x + mu) ** 2 + y ** 2)
    r2 = np.sqrt((x - 1.0 + mu) ** 2 + y ** 2)
    return 0.5 * (x ** 2 + y ** 2) + (1.0 - mu) / r1 + mu / r2


def jacobi_constant(state: np.ndarray, mu: float = MU) -> np.ndarray | float:
    """
    Compute the Jacobi integral C_J = 2*Omega(x, y) - (vx^2 + vy^2).
    Accepts shape (4,) or (N, 4).
    """
    state = np.asarray(state, dtype=float)
    if state.ndim == 1:
        x, y, vx, vy = state
        return float(2.0 * pseudo_potential(x, y, mu) - (vx ** 2 + vy ** 2))
    x, y, vx, vy = state[:, 0], state[:, 1], state[:, 2], state[:, 3]
    return 2.0 * pseudo_potential(x, y, mu) - (vx ** 2 + vy ** 2)


def cr3bp_eom(_t: float, state: np.ndarray, mu: float = MU) -> np.ndarray:
    """Planar CR3BP equations of motion for state [x, y, vx, vy]."""
    x, y, vx, vy = state
    dx1 = x + mu
    dx2 = x - 1.0 + mu
    r1_3 = (dx1 * dx1 + y * y) ** 1.5
    r2_3 = (dx2 * dx2 + y * y) ** 1.5

    ax = 2.0 * vy + x - (1.0 - mu) * dx1 / r1_3 - mu * dx2 / r2_3
    ay = -2.0 * vx + y - (1.0 - mu) * y / r1_3 - mu * y / r2_3
    return np.array([vx, vy, ax, ay], dtype=float)


def cr3bp_eom_stm(_t: float, aug_state: np.ndarray, mu: float = MU) -> np.ndarray:
    """
    Planar CR3BP equations of motion augmented with the 4x4 State Transition Matrix (STM).
    aug_state has length 20: [x, y, vx, vy, phi_00, ..., phi_33].
    """
    x, y, vx, vy = aug_state[:4]
    phi = aug_state[4:].reshape((4, 4))

    dx1 = x + mu
    dx2 = x - 1.0 + mu
    r1_2 = dx1 * dx1 + y * y
    r2_2 = dx2 * dx2 + y * y
    r1_3 = r1_2 ** 1.5
    r2_3 = r2_2 ** 1.5
    r1_5 = r1_2 ** 2.5
    r2_5 = r2_2 ** 2.5

    ax = 2.0 * vy + x - (1.0 - mu) * dx1 / r1_3 - mu * dx2 / r2_3
    ay = -2.0 * vx + y - (1.0 - mu) * y / r1_3 - mu * y / r2_3

    # Second partial derivatives of Omega
    Uxx = (
        1.0
        - (1.0 - mu) / r1_3
        - mu / r2_3
        + 3.0 * (1.0 - mu) * dx1 * dx1 / r1_5
        + 3.0 * mu * dx2 * dx2 / r2_5
    )
    Uyy = (
        1.0
        - (1.0 - mu) / r1_3
        - mu / r2_3
        + 3.0 * (1.0 - mu) * y * y / r1_5
        + 3.0 * mu * y * y / r2_5
    )
    Uxy = (
        3.0 * (1.0 - mu) * dx1 * y / r1_5
        + 3.0 * mu * dx2 * y / r2_5
    )

    A = np.array(
        [
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
            [Uxx, Uxy, 0.0, 2.0],
            [Uxy, Uyy, -2.0, 0.0],
        ],
        dtype=float,
    )

    dphi = A @ phi
    out = np.empty(20, dtype=float)
    out[:4] = [vx, vy, ax, ay]
    out[4:] = dphi.ravel()
    return out


def compute_lagrange_points(mu: float = MU) -> dict[str, tuple[float, float]]:
    """Compute exact planar coordinates (x, y) of L1..L5 in the rotating frame."""

    def dOmega_dx(x: float) -> float:
        return (
            x
            - (1.0 - mu) * (x + mu) / abs(x + mu) ** 3
            - mu * (x - 1.0 + mu) / abs(x - 1.0 + mu) ** 3
        )

    l1_x = brentq(dOmega_dx, -mu + 1e-6, 1.0 - mu - 1e-6, xtol=1e-14)
    l2_x = brentq(dOmega_dx, 1.0 - mu + 1e-6, 2.0, xtol=1e-14)
    l3_x = brentq(dOmega_dx, -2.0, -mu - 1e-6, xtol=1e-14)
    l4 = (0.5 - mu, np.sqrt(3.0) / 2.0)
    l5 = (0.5 - mu, -np.sqrt(3.0) / 2.0)

    return {
        "L1": (float(l1_x), 0.0),
        "L2": (float(l2_x), 0.0),
        "L3": (float(l3_x), 0.0),
        "L4": (float(l4[0]), float(l4[1])),
        "L5": (float(l5[0]), float(l5[1])),
    }


def propagate(
    state0: np.ndarray,
    t_span: tuple[float, float],
    t_eval: np.ndarray | None = None,
    events=None,
    with_stm: bool = False,
    rtol: float = 1e-12,
    atol: float = 1e-13,
    mu: float = MU,
):
    """High-accuracy DOP853 numerical integration of the planar CR3BP."""
    if with_stm:
        y0 = np.concatenate([np.asarray(state0, dtype=float), np.eye(4).ravel()])
        fun = lambda t, y: cr3bp_eom_stm(t, y, mu)
    else:
        y0 = np.asarray(state0, dtype=float)
        fun = lambda t, y: cr3bp_eom(t, y, mu)

    return solve_ivp(
        fun,
        t_span,
        y0,
        method="DOP853",
        t_eval=t_eval,
        events=events,
        rtol=rtol,
        atol=atol,
    )


def rotating_to_eci(t: np.ndarray, states_rot: np.ndarray, mu: float = MU) -> dict[str, np.ndarray]:
    """
    Transform rotating barycentric states (N, 4) into Earth-Centered Inertial (ECI) frame.
    Returns positions (x_eci, y_eci) and velocities (vx_eci, vy_eci) in dimensionless units,
    plus Moon ECI position (xm_eci, ym_eci).
    """
    t = np.asarray(t, dtype=float)
    x = states_rot[:, 0]
    y = states_rot[:, 1]
    vx = states_rot[:, 2]
    vy = states_rot[:, 3]

    c = np.cos(t)
    s = np.sin(t)

    # Shift origin from barycenter to Earth: x_e = x + mu, y_e = y
    xe = x + mu
    ye = y

    x_eci = xe * c - ye * s
    y_eci = xe * s + ye * c

    # Inertial velocity relative to Earth: v_I = R(t) * (v_rot + omega x r_e), omega = (0,0,1)
    vxe_rot = vx - ye
    vye_rot = vy + xe
    vx_eci = vxe_rot * c - vye_rot * s
    vy_eci = vxe_rot * s + vye_rot * c

    return {
        "x_eci": x_eci,
        "y_eci": y_eci,
        "vx_eci": vx_eci,
        "vy_eci": vy_eci,
        "xm_eci": c,
        "ym_eci": s,
    }


def differential_correct_symmetric(
    x0: float,
    vy0: float,
    t_half_guess: float,
    crossing_index: int = 1,
    fix_x0: bool = True,
    max_iter: int = 30,
    tol: float = 1e-11,
    mu: float = MU,
) -> dict:
    """
    Single-shooting differential corrector for x-axis symmetric periodic orbits in CR3BP.

    Starts at [x0, 0, 0, vy0] and enforces a perpendicular x-axis crossing
    (y = 0, vx = 0) at the `crossing_index`-th crossing of y = 0.
    By CR3BP symmetry, the full period is T = 2 * t_half.
    """
    x_curr = float(x0)
    vy_curr = float(vy0)

    def y_cross(_t, y_aug):
        return y_aug[1]

    y_cross.terminal = False
    y_cross.direction = 0

    for iteration in range(max_iter):
        state0 = np.array([x_curr, 0.0, 0.0, vy_curr], dtype=float)
        t_max = max(t_half_guess * 1.8, 1.0)
        sol = propagate(state0, (0.0, t_max), events=y_cross, with_stm=True, mu=mu)

        # Filter out t ~ 0
        valid_mask = sol.t_events[0] > 1e-3
        t_hits = sol.t_events[0][valid_mask]
        y_hits = sol.y_events[0][valid_mask]

        if len(t_hits) < crossing_index:
            raise RuntimeError(
                f"Only found {len(t_hits)} x-axis crossings (needed {crossing_index}) "
                f"at iter {iteration}, state0={state0}."
            )

        t_half = float(t_hits[crossing_index - 1])
        hit = y_hits[crossing_index - 1]
        state_half = hit[:4]
        phi = hit[4:].reshape((4, 4))
        vx_f = state_half[2]
        vy_f = state_half[3]

        if abs(vx_f) < tol:
            # Propagate full period with STM to compute monodromy matrix & eigenvalues
            full_sol = propagate(state0, (0.0, 2.0 * t_half), with_stm=True, mu=mu)
            final_state = full_sol.y[:4, -1]
            monodromy = full_sol.y[4:, -1].reshape((4, 4))
            eigvals = np.linalg.eigvals(monodromy)
            stability_index = 0.5 * float(np.max(np.abs(eigvals) + 1.0 / np.abs(eigvals)))
            return {
                "state0": state0,
                "half_period": t_half,
                "period": 2.0 * t_half,
                "state_half": state_half,
                "closure_error": float(np.linalg.norm(final_state - state0)),
                "jacobi": jacobi_constant(state0, mu),
                "monodromy": monodromy,
                "eigenvalues": eigvals,
                "stability_index": stability_index,
                "iterations": iteration,
            }

        # State derivative at half-period
        deriv_f = cr3bp_eom(t_half, state_half, mu)
        ax_f = deriv_f[2]

        # Variation of vx(t_half) holding y(t_half) = 0:
        # delta_vx = (phi[2, j] - (ax_f / vy_f) * phi[1, j]) * delta_q_j
        col = 3 if fix_x0 else 0
        dvx_dq = phi[2, col] - (ax_f / vy_f) * phi[1, col]
        delta_q = -vx_f / dvx_dq
        # Damp large jumps for numerical robustness
        max_step = 0.05
        if abs(delta_q) > max_step:
            delta_q = np.sign(delta_q) * max_step

        if fix_x0:
            vy_curr += delta_q
        else:
            x_curr += delta_q
        t_half_guess = t_half

    raise RuntimeError(f"Differential corrector did not converge after {max_iter} iterations.")
