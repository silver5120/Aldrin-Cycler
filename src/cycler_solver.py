"""
Numerical Solver for Earth-Moon Aldrin Cyclers & Astrodynamics Trade Study.

Computes two exact periodic Earth-Moon cycler trajectories in the planar CR3BP:
  1. Aldrin-Genova "Shamrock" (C-2-R) Cycler:
     - 3:1 lunar-resonant circumlunar free-return cycler.
     - Ballistic figure-8 leg (~6.8 days) passing Earth perigee (3,000 km) and
       lunar farside perilune (3,000 km), joined to two sub-lunar phasing loops
       (~9.5 days each) via two tiny symmetric tangential perigee burns (~39 m/s total).
  2. Pure Ballistic Arenstorf 3-Lobe Periodic Cycler (Delta-V = 0 m/s):
     - Solved via STM differential correction in the planar CR3BP.
     - Demonstrates the fundamental ballistic CR3BP trade-off: zero deterministic
       station-keeping Delta-V, but poor Earth proximity (>170,000 km min altitude)
       and strong linear instability (large Floquet stability index).

Also computes taxi rendezvous Delta-V budgets and a Tsiolkovsky rocket-equation
cumulative IMLEO / propellant mass trade study comparing Cycler + Taxi vs. Direct.
"""

from __future__ import annotations

import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

from cr3bp import (
    L_STAR_KM,
    LUNAR_MONTH_DAYS,
    MU,
    MU_EARTH_KM3_S2,
    MU_MOON_KM3_S2,
    R_EARTH_KM,
    R_MOON_KM,
    T_STAR_DAYS,
    V_STAR_KM_S,
    compute_lagrange_points,
    differential_correct_symmetric,
    jacobi_constant,
    propagate,
    rotating_to_eci,
)


def apply_eci_tangential_burn(state: np.ndarray, dv_nd: float, mu: float = MU) -> np.ndarray:
    """Apply an instantaneous tangential Delta-V along the Earth-centered inertial velocity."""
    x, y, vx, vy = state
    xe, ye = x + mu, y
    vxe = vx - ye
    vye = vy + xe
    speed_eci = np.hypot(vxe, vye)
    ux, uy = vxe / speed_eci, vye / speed_eci
    return np.array([x, y, vx + dv_nd * ux, vy + dv_nd * uy], dtype=float)


def solve_shamrock_cycler(
    perilune_alt_km: float = 3000.0,
    perigee_alt_km: float = 3000.0,
    samples_per_leg: int = 400,
) -> dict:
    """
    Solve for the exact symmetric Aldrin-Genova Shamrock (C-2-R) Earth-Moon cycler.

    Trajectory structure (1 full cycle ~25.9 days):
      - Segment 0 (Outbound Figure-8): Outbound Earth perigee -> Lunar farside perilune
      - Segment 1 (Inbound Figure-8):  Lunar farside perilune -> Inbound Earth perigee
      - Burn 1 (Perigee Retrograde):   Tiny tangential burn (-dv1) into sub-lunar phasing
      - Segment 2 (Phasing Loop 1):    Inbound Earth perigee -> Middle x-axis perigee (vx=0)
      - Segment 3 (Phasing Loop 2):    Middle perigee -> Outbound Earth perigee (mirror image)
      - Burn 2 (Perigee Prograde):     Matching tangential burn (+dv1) to re-enter Figure-8
    """
    x_moon_far = (1.0 - MU) + (R_MOON_KM + perilune_alt_km) / L_STAR_KM

    def earth_perigee_event(_t, y):
        return (y[0] + MU) * y[2] + y[1] * y[3]

    earth_perigee_event.direction = 1
    earth_perigee_event.terminal = True

    def get_fig8_perigee(vy0: float):
        sol = propagate([x_moon_far, 0.0, 0.0, vy0], (0.0, 2.0), events=earth_perigee_event)
        if len(sol.t_events[0]) == 0:
            raise RuntimeError(f"No Earth perigee found for vy0={vy0}")
        t_p = float(sol.t_events[0][0])
        s_p = sol.y_events[0][0].copy()
        h_e = float(np.hypot(s_p[0] + MU, s_p[1]) * L_STAR_KM - R_EARTH_KM)
        return t_p, s_p, h_e

    # 1. Solve for lunar farside speed vy0 that yields exact target Earth perigee altitude
    vy0_fig8 = brentq(
        lambda v: get_fig8_perigee(v)[2] - perigee_alt_km,
        -1.655,
        -1.705,
        xtol=1e-13,
    )
    t_half_fig8, s_perigee_in, h_e_check = get_fig8_perigee(vy0_fig8)
    s_perilune = np.array([x_moon_far, 0.0, 0.0, vy0_fig8], dtype=float)
    # Mirror-symmetric outbound Earth perigee state
    s_perigee_out = np.array(
        [s_perigee_in[0], -s_perigee_in[1], -s_perigee_in[2], s_perigee_in[3]],
        dtype=float,
    )

    # 2. Solve for exact retrograde tangential burn dv1 (m/s) that hits y=0 perpendicularly (vx=0)
    #    at the midpoint between the two sub-lunar phasing loops.
    def y_cross(_t, y):
        return y[1]

    y_cross.direction = 0

    def phasing_midpoint(dv_ms: float):
        dv_nd = (dv_ms / 1000.0) / V_STAR_KM_S
        s_after_b1 = apply_eci_tangential_burn(s_perigee_in, dv_nd)
        sol = propagate(s_after_b1, (0.0, 3.5), events=y_cross)
        # First x-axis crossing near Earth after 1 phasing orbit (t > 1.0, r1 < 0.15)
        for t_h, y_h in zip(sol.t_events[0], sol.y_events[0]):
            if t_h > 1.0 and np.hypot(y_h[0] + MU, y_h[1]) < 0.15:
                return float(t_h), y_h.copy(), s_after_b1
        raise RuntimeError(f"No phasing midpoint crossing found for dv={dv_ms} m/s")

    dv1_ms = brentq(lambda dv: phasing_midpoint(dv)[1][2], -30.0, -10.0, xtol=1e-11)
    t_phase_loop, s_mid_perigee, s_after_b1 = phasing_midpoint(dv1_ms)
    dv_burn1_ms = abs(dv1_ms)
    dv_burn2_ms = dv_burn1_ms
    dv_total_ms = dv_burn1_ms + dv_burn2_ms

    # 3. Propagate the complete closed cycle starting from Outbound Earth Perigee (t = 0)
    #    Leg A: Outbound + Inbound Figure-8 (0 -> 2*t_half_fig8)
    t_fig8_total = 2.0 * t_half_fig8
    t_grid_fig8 = np.linspace(0.0, t_fig8_total, samples_per_leg * 2, endpoint=False)
    sol_fig8 = propagate(s_perigee_out, (0.0, t_fig8_total), t_eval=t_grid_fig8)

    # Verify exact arrival at s_perigee_in
    sol_fig8_end = propagate(s_perigee_out, (0.0, t_fig8_total))
    fig8_err = float(np.linalg.norm(sol_fig8_end.y[:, -1] - s_perigee_in))

    #    Leg B: Two Phasing Loops (t_fig8_total -> t_fig8_total + 2*t_phase_loop)
    t_phase_total = 2.0 * t_phase_loop
    t_grid_phase = np.linspace(0.0, t_phase_total, samples_per_leg * 3, endpoint=True)
    sol_phase = propagate(s_after_b1, (0.0, t_phase_total), t_eval=t_grid_phase)

    # Verify exact closure after Burn 2
    s_before_b2 = sol_phase.y[:, -1]
    s_after_b2 = apply_eci_tangential_burn(s_before_b2, +(dv_burn2_ms / 1000.0) / V_STAR_KM_S)
    closure_err = float(np.linalg.norm(s_after_b2 - s_perigee_out))

    # Stitch full cycle timeline
    t_all = np.concatenate([t_grid_fig8, t_fig8_total + t_grid_phase])
    states_all = np.vstack([sol_fig8.y.T, sol_phase.y.T])
    phase_labels = np.concatenate(
        [
            np.where(t_grid_fig8 < t_half_fig8, 0, 1),  # 0: Outbound, 1: Inbound
            np.where(t_grid_phase < t_phase_loop, 2, 3),  # 2: Phasing 1, 3: Phasing 2
        ]
    )

    # Compute inertial frame and physical metrics
    eci = rotating_to_eci(t_all, states_all)
    r_earth_km = np.hypot(states_all[:, 0] + MU, states_all[:, 1]) * L_STAR_KM
    r_moon_km = np.hypot(states_all[:, 0] - (1.0 - MU), states_all[:, 1]) * L_STAR_KM
    alt_earth_km = r_earth_km - R_EARTH_KM
    alt_moon_km = r_moon_km - R_MOON_KM
    v_eci_km_s = np.hypot(eci["vx_eci"], eci["vy_eci"]) * V_STAR_KM_S
    v_rot_km_s = np.hypot(states_all[:, 2], states_all[:, 3]) * V_STAR_KM_S
    cj_arr = jacobi_constant(states_all)

    # Selenocentric velocity at perilune
    # v_moon_rel_inertial = v_rot + omega x r_moon_rel
    xm_rel = s_perilune[0] - (1.0 - MU)
    ym_rel = s_perilune[1]
    vxm_moon = s_perilune[2] - ym_rel
    vym_moon = s_perilune[3] + xm_rel
    v_perilune_moon_km_s = float(np.hypot(vxm_moon, vym_moon) * V_STAR_KM_S)
    r_perilune_km = R_MOON_KM + perilune_alt_km
    v_inf_moon_km_s = float(
        np.sqrt(max(0.0, v_perilune_moon_km_s ** 2 - 2.0 * MU_MOON_KM3_S2 / r_perilune_km))
    )

    # Geocentric ECI velocity at perigee (on figure-8 leg)
    xe_p = s_perigee_out[0] + MU
    ye_p = s_perigee_out[1]
    v_perigee_eci_km_s = float(
        np.hypot(s_perigee_out[2] - ye_p, s_perigee_out[3] + xe_p) * V_STAR_KM_S
    )

    period_nd = t_fig8_total + t_phase_total
    mid_perigee_alt_km = float(np.hypot(s_mid_perigee[0] + MU, s_mid_perigee[1]) * L_STAR_KM - R_EARTH_KM)

    # Compute phasing apogee altitude
    phase_mask = phase_labels >= 2
    phasing_apogee_km = float(np.max(alt_earth_km[phase_mask]))

    # Multi-cycle ECI propagation (3 full cycles = ~78 days) for inertial visualization
    t_multi = []
    states_multi = []
    for k in range(3):
        t_offset = k * period_nd
        if k == 0:
            t_multi.append(t_all)
            states_multi.append(states_all)
        else:
            t_multi.append(t_offset + t_all[1:])
            states_multi.append(states_all[1:])
    t_multi = np.concatenate(t_multi)
    states_multi = np.vstack(states_multi)
    eci_multi = rotating_to_eci(t_multi, states_multi)

    return {
        "name": "Aldrin-Genova Shamrock (C-2-R)",
        "t_nd": t_all,
        "t_days": t_all * T_STAR_DAYS,
        "states_rot": states_all,
        "phase_labels": phase_labels,
        "eci": eci,
        "t_multi_days": t_multi * T_STAR_DAYS,
        "eci_multi": eci_multi,
        "alt_earth_km": alt_earth_km,
        "alt_moon_km": alt_moon_km,
        "v_eci_km_s": v_eci_km_s,
        "v_rot_km_s": v_rot_km_s,
        "jacobi_series": cj_arr,
        "metrics": {
            "period_days": float(period_nd * T_STAR_DAYS),
            "period_lunar_months": float(period_nd / (2.0 * np.pi)),
            "earth_to_moon_transit_days": float(t_half_fig8 * T_STAR_DAYS),
            "fig8_leg_days": float(t_fig8_total * T_STAR_DAYS),
            "phasing_loop_days": float(t_phase_loop * T_STAR_DAYS),
            "perigee_alt_km": float(h_e_check),
            "mid_perigee_alt_km": mid_perigee_alt_km,
            "perilune_alt_km": float(perilune_alt_km),
            "phasing_apogee_alt_km": phasing_apogee_km,
            "v_perigee_eci_km_s": v_perigee_eci_km_s,
            "v_perilune_moon_km_s": v_perilune_moon_km_s,
            "v_inf_moon_km_s": v_inf_moon_km_s,
            "dv_burn1_m_s": float(dv_burn1_ms),
            "dv_burn2_m_s": float(dv_burn2_ms),
            "dv_total_per_cycle_m_s": float(dv_total_ms),
            "dv_annual_m_s": float(dv_total_ms * (365.25 / (period_nd * T_STAR_DAYS))),
            "jacobi_fig8": float(jacobi_constant(s_perilune)),
            "jacobi_phasing": float(jacobi_constant(s_after_b1)),
            "fig8_symmetry_error": fig8_err,
            "closure_error": closure_err,
            "burn1_time_days": float(t_fig8_total * T_STAR_DAYS),
            "burn2_time_days": float(period_nd * T_STAR_DAYS),
            "perilune_time_days": float(t_half_fig8 * T_STAR_DAYS),
            "state_perigee_out": s_perigee_out.tolist(),
            "state_perilune": s_perilune.tolist(),
        },
    }


def solve_ballistic_arenstorf(samples: int = 2000) -> dict:
    """
    Solve the classic 3-lobe ballistic Arenstorf periodic orbit in the Earth-Moon CR3BP.
    Requires Delta-V = 0 m/s deterministically, illustrating the ballistic CR3BP trade-off.
    """
    res = differential_correct_symmetric(
        x0=0.994, vy0=-2.0015851, t_half_guess=8.53, crossing_index=3, fix_x0=True
    )
    t_eval = np.linspace(0.0, res["period"], samples)
    sol = propagate(res["state0"], (0.0, res["period"]), t_eval=t_eval)
    states = sol.y.T
    eci = rotating_to_eci(t_eval, states)

    r_earth_km = np.hypot(states[:, 0] + MU, states[:, 1]) * L_STAR_KM
    r_moon_km = np.hypot(states[:, 0] - (1.0 - MU), states[:, 1]) * L_STAR_KM
    alt_earth_km = r_earth_km - R_EARTH_KM
    alt_moon_km = r_moon_km - R_MOON_KM
    v_eci_km_s = np.hypot(eci["vx_eci"], eci["vy_eci"]) * V_STAR_KM_S
    cj_arr = jacobi_constant(states)

    return {
        "name": "Pure Ballistic Arenstorf (3-Lobe)",
        "t_nd": t_eval,
        "t_days": t_eval * T_STAR_DAYS,
        "states_rot": states,
        "eci": eci,
        "alt_earth_km": alt_earth_km,
        "alt_moon_km": alt_moon_km,
        "v_eci_km_s": v_eci_km_s,
        "jacobi_series": cj_arr,
        "metrics": {
            "period_days": float(res["period"] * T_STAR_DAYS),
            "period_lunar_months": float(res["period"] / (2.0 * np.pi)),
            "min_earth_alt_km": float(np.min(alt_earth_km)),
            "min_moon_alt_km": float(np.min(alt_moon_km)),
            "dv_total_per_cycle_m_s": 0.0,
            "jacobi": float(res["jacobi"]),
            "jacobi_drift": float(np.max(np.abs(cj_arr - res["jacobi"]))),
            "closure_error": float(res["closure_error"]),
            "stability_index": float(res["stability_index"]),
            "max_floquet_multiplier": float(np.max(np.abs(res["eigenvalues"]))),
            "state0": res["state0"].tolist(),
        },
    }


def compute_mission_trade_study(shamrock_metrics: dict) -> dict:
    """
    Compute Earth/Moon taxi Delta-V budgets and cumulative IMLEO mass comparison
    (Aldrin Cycler + Light Taxi vs. Direct Heavy Transfer) via the Tsiolkovsky equation.
    """
    # 1. Earth Taxi: 400 km circular LEO -> rendezvous with cycler at Earth perigee
    r_leo = R_EARTH_KM + 400.0
    r_p_cyc = R_EARTH_KM + shamrock_metrics["perigee_alt_km"]
    v_leo = np.sqrt(MU_EARTH_KM3_S2 / r_leo)
    a_trans = 0.5 * (r_leo + r_p_cyc)
    v_trans_p = np.sqrt(MU_EARTH_KM3_S2 * (2.0 / r_leo - 1.0 / a_trans))
    v_trans_a = np.sqrt(MU_EARTH_KM3_S2 * (2.0 / r_p_cyc - 1.0 / a_trans))
    dv_leo_depart = v_trans_p - v_leo
    dv_perigee_match = shamrock_metrics["v_perigee_eci_km_s"] - v_trans_a
    dv_earth_taxi_km_s = float(dv_leo_depart + dv_perigee_match)

    # Also compute direct TLI from 400 km LEO to 384,400 km lunar transfer ellipse (reference)
    a_tli = 0.5 * (r_leo + L_STAR_KM)
    dv_direct_tli_km_s = float(np.sqrt(MU_EARTH_KM3_S2 * (2.0 / r_leo - 1.0 / a_tli)) - v_leo)

    # 2. Lunar Taxi: Hyperbolic drop from cycler v_inf -> 100 km circular LLO (and return)
    r_llo = R_MOON_KM + 100.0
    v_llo = np.sqrt(MU_MOON_KM3_S2 / r_llo)
    v_inf_l = shamrock_metrics["v_inf_moon_km_s"]
    v_hyp_llo = np.sqrt(v_inf_l ** 2 + 2.0 * MU_MOON_KM3_S2 / r_llo)
    dv_lunar_capture_km_s = float(v_hyp_llo - v_llo)
    dv_lunar_depart_km_s = dv_lunar_capture_km_s  # symmetric return to cycler
    dv_earth_entry_trim_km_s = 0.05  # direct atmospheric entry (aerobraking/splashdown)

    # 3. Tsiolkovsky Rocket Equation comparison (LOX/CH4 or LOX/LH2 storable/cryo mix)
    g0 = 9.80665e-3  # km/s^2
    isp_earth_stage = 420.0  # s (cryogenic upper stage for TLI / taxi boost)
    isp_inspace = 360.0      # s (methalox / storable in-space propulsion)
    ve_earth = isp_earth_stage * g0
    ve_inspace = isp_inspace * g0

    # Dry masses (metric tons) for a 4-crew radiation-shielded transit architecture:
    # - Heavy Transit Habitat (radiation shielding, ECLSS, centrifuge/gym, galley): 35 t
    # - Light Crew Taxi Capsule (short-duration docking + reentry capsule): 10 t
    m_hab_t = 35.0
    m_taxi_t = 10.0
    m_consumables_per_flight_t = 2.0

    # Direct mission (each flight carries the 35 t heavy hab + 10 t reentry capsule = 45 t dry
    # through TLI, LOI into LLO, and TEI back to Earth reentry):
    # Backward mass staging:
    #   TEI burn (LLO -> Earth return) on (m_taxi + m_hab_returned_or_reused):
    #   Suppose direct carries 35 t hab + 10 t capsule = 45 t to LLO and returns 45 t?
    #   Even generously assuming direct only returns the 10 t capsule + 20 t hab (leaving 15 t or
    #   using a 25 t combined Orion+Hab): let's use an apples-to-apples 35 t comfortable transit
    #   habitat + 10 t capsule (45 t total dry) vs Cycler (35 t stays on cycler forever; only
    #   10 t taxi + 2 t consumables = 12 t makes the LEO->Cycler and Cycler<->LLO burns!).
    # Direct round-trip IMLEO per flight (conservative: 35 t hab + 10 t capsule = 45 t to LLO,
    # 10 t capsule + 15 t reusable hab core = 25 t returns via TEI; or full 45 t if reusable):
    # Let's model apples-to-apples: BOTH systems provide a 35 t shielded habitat during transit
    # and deliver a 10 t crew capsule to LLO and back to Earth!
    m_direct_post_tei = m_taxi_t + 0.0  # Even if direct discards the hab in LLO (minimum possible TEI!)
    # Wait: if direct discards the 35 t hab every time, it launches 35 t of new hardware every flight!
    # Either way, Direct must push (m_hab + m_taxi + m_cons + m_tei_prop) through LOI and TLI.
    m_direct_before_tei = (m_taxi_t + 0.5 * m_consumables_per_flight_t) * np.exp(
        dv_lunar_depart_km_s / ve_inspace
    )
    m_direct_before_loi = (m_hab_t + m_direct_before_tei + 0.5 * m_consumables_per_flight_t) * np.exp(
        dv_lunar_capture_km_s / ve_inspace
    )
    imleo_direct_per_flight_t = float(m_direct_before_loi * np.exp(dv_direct_tli_km_s / ve_earth))

    # Cycler architecture:
    #   (a) One-time upfront IMLEO to inject the 35 t shielded Cycler Habitat into the Shamrock orbit
    imleo_cycler_setup_t = float(m_hab_t * np.exp(dv_earth_taxi_km_s / ve_earth))
    #   (b) Per-flight taxi IMLEO (only 10 t taxi + 2 t consumables):
    #       Outbound: taxi departs LEO -> docks with cycler -> undocks at Moon -> LOI into LLO
    #       Return:   taxi departs LLO -> docks with cycler -> undocks at Earth -> direct reentry
    m_taxi_before_tei = (m_taxi_t + 0.5 * m_consumables_per_flight_t) * np.exp(
        (dv_lunar_depart_km_s + dv_earth_entry_trim_km_s) / ve_inspace
    )
    m_taxi_before_loi = (m_taxi_before_tei + 0.5 * m_consumables_per_flight_t) * np.exp(
        dv_lunar_capture_km_s / ve_inspace
    )
    #   Plus cycler station-keeping propellant (~39 m/s per cycle, 1 cycle per flight)
    m_sk_prop_per_cycle_t = float(
        m_hab_t * (np.exp((shamrock_metrics["dv_total_per_cycle_m_s"] / 1000.0) / ve_inspace) - 1.0)
    )
    imleo_cycler_per_flight_t = float(
        (m_taxi_before_loi + m_sk_prop_per_cycle_t) * np.exp(dv_earth_taxi_km_s / ve_earth)
    )

    # Also compute a leaner Apollo/Orion-style direct capsule (18 t dry: 10 t capsule + 8 t service module,
    # viable specifically because Earth-Moon transit is only ~3.4 days, unlike 6-month Earth-Mars!)
    m_apollo_dry = 18.0
    m_apollo_pre_tei = (10.0 + 0.5 * m_consumables_per_flight_t) * np.exp(
        dv_lunar_depart_km_s / ve_inspace
    )
    m_apollo_pre_loi = (
        (m_apollo_dry - 10.0) + m_apollo_pre_tei + 0.5 * m_consumables_per_flight_t
    ) * np.exp(dv_lunar_capture_km_s / ve_inspace)
    imleo_apollo_lean_t = float(m_apollo_pre_loi * np.exp(dv_direct_tli_km_s / ve_earth))

    flights = list(range(1, 26))
    cum_direct_t = [round(imleo_direct_per_flight_t * n, 1) for n in flights]
    cum_apollo_lean_t = [round(imleo_apollo_lean_t * n, 1) for n in flights]
    cum_cycler_t = [
        round(imleo_cycler_setup_t + imleo_cycler_per_flight_t * n, 1) for n in flights
    ]

    # Exact fractional break-even flight counts
    breakeven_vs_heavy = (
        imleo_cycler_setup_t / (imleo_direct_per_flight_t - imleo_cycler_per_flight_t)
        if imleo_direct_per_flight_t > imleo_cycler_per_flight_t
        else float("inf")
    )
    breakeven_vs_lean = (
        imleo_cycler_setup_t / (imleo_apollo_lean_t - imleo_cycler_per_flight_t)
        if imleo_apollo_lean_t > imleo_cycler_per_flight_t
        else float("inf")
    )

    return {
        "dv_earth_taxi_km_s": dv_earth_taxi_km_s,
        "dv_leo_depart_km_s": float(dv_leo_depart),
        "dv_perigee_match_km_s": float(dv_perigee_match),
        "dv_direct_tli_km_s": dv_direct_tli_km_s,
        "dv_lunar_capture_km_s": dv_lunar_capture_km_s,
        "dv_lunar_depart_km_s": dv_lunar_depart_km_s,
        "isp_earth_stage_s": isp_earth_stage,
        "isp_inspace_s": isp_inspace,
        "m_hab_t": m_hab_t,
        "m_taxi_t": m_taxi_t,
        "m_apollo_dry_t": m_apollo_dry,
        "m_sk_prop_per_cycle_t": m_sk_prop_per_cycle_t,
        "imleo_cycler_setup_t": imleo_cycler_setup_t,
        "imleo_cycler_per_flight_t": imleo_cycler_per_flight_t,
        "imleo_direct_per_flight_t": imleo_direct_per_flight_t,
        "imleo_apollo_lean_t": imleo_apollo_lean_t,
        "breakeven_vs_heavy_flights": float(breakeven_vs_heavy),
        "breakeven_vs_lean_flights": float(breakeven_vs_lean),
        "flights": flights,
        "cum_direct_imleo_t": cum_direct_t,
        "cum_apollo_lean_imleo_t": cum_apollo_lean_t,
        "cum_cycler_imleo_t": cum_cycler_t,
        "savings_vs_heavy_at_20_pct": float(
            100.0 * (1.0 - cum_cycler_t[19] / cum_direct_t[19])
        ),
        "savings_vs_lean_at_20_pct": float(
            100.0 * (1.0 - cum_cycler_t[19] / cum_apollo_lean_t[19])
        ),
    }


def export_solution(out_path: Path) -> dict:
    """Solve all orbits + trade study and save compact JSON for web/plotting."""
    shamrock = solve_shamrock_cycler(perilune_alt_km=3000.0, perigee_alt_km=3000.0)
    arenstorf = solve_ballistic_arenstorf(samples=1600)
    trade = compute_mission_trade_study(shamrock["metrics"])
    lagrange = compute_lagrange_points()

    # Downsample smoothly for lightweight JSON (still super crisp for HTML5 Canvas!)
    def decimate(arr, step):
        return np.round(arr[::step], 6).tolist()

    step_s = 2  # ~1000 points per Shamrock cycle
    step_a = 2  # ~800 points per Arenstorf cycle
    step_m = 4  # ~1500 points for 3-cycle ECI trail

    payload = {
        "constants": {
            "mu": MU,
            "L_star_km": L_STAR_KM,
            "T_star_days": T_STAR_DAYS,
            "V_star_km_s": V_STAR_KM_S,
            "R_earth_km": R_EARTH_KM,
            "R_moon_km": R_MOON_KM,
            "lunar_month_days": LUNAR_MONTH_DAYS,
        },
        "lagrange_points": lagrange,
        "shamrock": {
            "metrics": shamrock["metrics"],
            "t_days": decimate(shamrock["t_days"], step_s),
            "t_nd": decimate(shamrock["t_nd"], step_s),
            "x_rot": decimate(shamrock["states_rot"][:, 0], step_s),
            "y_rot": decimate(shamrock["states_rot"][:, 1], step_s),
            "x_eci": decimate(shamrock["eci"]["x_eci"], step_s),
            "y_eci": decimate(shamrock["eci"]["y_eci"], step_s),
            "xm_eci": decimate(shamrock["eci"]["xm_eci"], step_s),
            "ym_eci": decimate(shamrock["eci"]["ym_eci"], step_s),
            "alt_earth_km": np.round(shamrock["alt_earth_km"][::step_s], 1).tolist(),
            "alt_moon_km": np.round(shamrock["alt_moon_km"][::step_s], 1).tolist(),
            "v_eci_km_s": np.round(shamrock["v_eci_km_s"][::step_s], 3).tolist(),
            "jacobi": np.round(shamrock["jacobi_series"][::step_s], 5).tolist(),
            "phase": shamrock["phase_labels"][::step_s].tolist(),
            "multi_x_eci": decimate(shamrock["eci_multi"]["x_eci"], step_m),
            "multi_y_eci": decimate(shamrock["eci_multi"]["y_eci"], step_m),
        },
        "arenstorf": {
            "metrics": arenstorf["metrics"],
            "t_days": decimate(arenstorf["t_days"], step_a),
            "t_nd": decimate(arenstorf["t_nd"], step_a),
            "x_rot": decimate(arenstorf["states_rot"][:, 0], step_a),
            "y_rot": decimate(arenstorf["states_rot"][:, 1], step_a),
            "x_eci": decimate(arenstorf["eci"]["x_eci"], step_a),
            "y_eci": decimate(arenstorf["eci"]["y_eci"], step_a),
            "xm_eci": decimate(arenstorf["eci"]["xm_eci"], step_a),
            "ym_eci": decimate(arenstorf["eci"]["ym_eci"], step_a),
            "alt_earth_km": np.round(arenstorf["alt_earth_km"][::step_a], 1).tolist(),
            "alt_moon_km": np.round(arenstorf["alt_moon_km"][::step_a], 1).tolist(),
            "v_eci_km_s": np.round(arenstorf["v_eci_km_s"][::step_a], 3).tolist(),
            "jacobi": np.round(arenstorf["jacobi_series"][::step_a], 5).tolist(),
        },
        "trade_study": trade,
    }

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return {"shamrock": shamrock, "arenstorf": arenstorf, "trade": trade, "lagrange": lagrange}


if __name__ == "__main__":
    root = Path(__file__).resolve().parent.parent
    result = export_solution(root / "data" / "cycler_solution.json")

    sm = result["shamrock"]["metrics"]
    am = result["arenstorf"]["metrics"]
    tr = result["trade"]

    print("=== ALDRIN-GENOVA SHAMROCK (C-2-R) EARTH-MOON CYCLER ===")
    print(f"  Total Cycle Period      : {sm['period_days']:.3f} days ({sm['period_lunar_months']:.3f} lunar months)")
    print(f"  Earth -> Moon Transit   : {sm['earth_to_moon_transit_days']:.3f} days (Figure-8 leg: {sm['fig8_leg_days']:.3f} d)")
    print(f"  Phasing Loop (x2)       : {sm['phasing_loop_days']:.3f} days each (apogee {sm['phasing_apogee_alt_km']:,.0f} km)")
    print(f"  Earth Perigee Altitude  : {sm['perigee_alt_km']:.1f} km (mid perigee: {sm['mid_perigee_alt_km']:.1f} km)")
    print(f"  Lunar Perilune Altitude : {sm['perilune_alt_km']:.1f} km (v_inf = {sm['v_inf_moon_km_s']:.3f} km/s)")
    print(f"  Perigee ECI Speed       : {sm['v_perigee_eci_km_s']:.3f} km/s")
    print(f"  Station-Keeping Delta-V : {sm['dv_burn1_m_s']:.2f} + {sm['dv_burn2_m_s']:.2f} = {sm['dv_total_per_cycle_m_s']:.2f} m/s per cycle ({sm['dv_annual_m_s']:.1f} m/s/yr)")
    print(f"  Jacobi Constant C_J     : {sm['jacobi_fig8']:.5f} (Fig-8) / {sm['jacobi_phasing']:.5f} (Phasing)")
    print(f"  Closure Error           : {sm['closure_error']:.2e} (dimensionless)")

    print("\n=== PURE BALLISTIC ARENSTORF (3-LOBE CR3BP, DV = 0) ===")
    print(f"  Total Cycle Period      : {am['period_days']:.3f} days ({am['period_lunar_months']:.3f} lunar months)")
    print(f"  Min Earth Altitude      : {am['min_earth_alt_km']:,.1f} km")
    print(f"  Min Moon Altitude       : {am['min_moon_alt_km']:,.1f} km")
    print(f"  Jacobi Constant C_J     : {am['jacobi']:.5f} (drift: {am['jacobi_drift']:.2e})")
    print(f"  Closure Error           : {am['closure_error']:.2e}")
    print(f"  Floquet Stability Index : {am['stability_index']:.2f} (max |lambda| = {am['max_floquet_multiplier']:.2f})")

    print("\n=== MISSION & TAXI DELTA-V TRADE STUDY ===")
    print(f"  Earth Taxi (400km LEO -> Cycler Perigee): {tr['dv_earth_taxi_km_s']:.3f} km/s (vs Direct TLI {tr['dv_direct_tli_km_s']:.3f} km/s)")
    print(f"  Lunar Taxi (Cycler <-> 100km LLO)       : {tr['dv_lunar_capture_km_s']:.3f} km/s each way")
    print(f"  Cycler Upfront IMLEO (35 t shielded hab): {tr['imleo_cycler_setup_t']:.1f} t")
    print(f"  Per-Flight IMLEO (Cycler + 10 t Taxi)   : {tr['imleo_cycler_per_flight_t']:.1f} t (incl. {tr['m_sk_prop_per_cycle_t']:.2f} t SK prop)")
    print(f"  Per-Flight IMLEO (Direct 35 t Hab+Taxi) : {tr['imleo_direct_per_flight_t']:.1f} t")
    print(f"  Per-Flight IMLEO (Lean Apollo 18 t)     : {tr['imleo_apollo_lean_t']:.1f} t")
    print(f"  Break-Even vs Heavy Hab Direct (45 t)   : {tr['breakeven_vs_heavy_flights']:.2f} flights")
    print(f"  Break-Even vs Lean Capsule Direct (18 t): {tr['breakeven_vs_lean_flights']:.2f} flights")
    print(f"  IMLEO Savings at 20 Flights (vs Heavy)  : {tr['savings_vs_heavy_at_20_pct']:.1f}%")
    print(f"  IMLEO Savings at 20 Flights (vs Lean)   : {tr['savings_vs_lean_at_20_pct']:.1f}%")
