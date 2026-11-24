"""Automated physics and numerical accuracy tests for Aldrin-Cycler."""

import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from cr3bp import (
    MU,
    compute_lagrange_points,
    cr3bp_eom,
    jacobi_constant,
    propagate,
)
from cycler_solver import (
    compute_mission_trade_study,
    solve_ballistic_arenstorf,
    solve_shamrock_cycler,
)


def test_lagrange_points_are_equilibria():
    lp = compute_lagrange_points(MU)
    assert set(lp.keys()) == {"L1", "L2", "L3", "L4", "L5"}
    for name, (lx, ly) in lp.items():
        deriv = cr3bp_eom(0.0, np.array([lx, ly, 0.0, 0.0]), MU)
        accel_norm = np.hypot(deriv[2], deriv[3])
        assert accel_norm < 1e-11, f"{name} acceleration not zero: {accel_norm}"


def test_jacobi_constant_conservation():
    state0 = np.array([0.85, 0.0, 0.0, 0.35])
    cj0 = jacobi_constant(state0)
    sol = propagate(state0, (0.0, 2.0), t_eval=np.linspace(0.0, 2.0, 200))
    cj_series = jacobi_constant(sol.y.T)
    max_drift = float(np.max(np.abs(cj_series - cj0)))
    assert max_drift < 1e-11, f"Jacobi drift too large: {max_drift}"


def test_shamrock_cycler_closure_and_physics():
    sh = solve_shamrock_cycler(perilune_alt_km=3000.0, perigee_alt_km=3000.0, samples_per_leg=120)
    m = sh["metrics"]
    assert abs(m["perigee_alt_km"] - 3000.0) < 0.1
    assert abs(m["perilune_alt_km"] - 3000.0) < 0.1
    assert 24.0 < m["period_days"] < 28.0
    assert 3.0 < m["earth_to_moon_transit_days"] < 4.0
    assert 15.0 < m["dv_total_per_cycle_m_s"] < 65.0
    assert m["closure_error"] < 1e-6
    assert m["fig8_symmetry_error"] < 1e-6


def test_ballistic_arenstorf_periodicity():
    ar = solve_ballistic_arenstorf(samples=400)
    m = ar["metrics"]
    assert m["dv_total_per_cycle_m_s"] == 0.0
    assert m["closure_error"] < 1e-8
    assert m["jacobi_drift"] < 1e-10
    assert m["min_moon_alt_km"] < 1000.0
    assert m["min_earth_alt_km"] > 100000.0  # Illustrates ballistic Earth-proximity constraint
    assert m["stability_index"] > 10.0       # Strongly unstable


def test_trade_study_breakeven_ordering():
    sh = solve_shamrock_cycler(perilune_alt_km=3000.0, perigee_alt_km=3000.0, samples_per_leg=60)
    tr = compute_mission_trade_study(sh["metrics"])
    assert tr["dv_earth_taxi_km_s"] > tr["dv_direct_tli_km_s"]  # Oberth penalty at 3000 km perigee
    assert tr["imleo_cycler_per_flight_t"] < tr["imleo_apollo_lean_t"] < tr["imleo_direct_per_flight_t"]
    assert tr["breakeven_vs_heavy_flights"] < tr["breakeven_vs_lean_flights"]
