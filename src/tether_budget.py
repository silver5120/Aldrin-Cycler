"""Earth-Moon momentum-exchange tether (skyhook / rotovator) budgets.

Prints:
  1. Tether-to-payload mass ratio vs tip speed for real fibers (Moravec/Hoyt
     tapered-tether formula: M_t/M_p = sqrt(pi) * x * exp(x^2) * erf(x), x = V/Vc,
     Vc = sqrt(2 * sigma / (SF * rho))).
  2. Full lunar rotovator kinematics (catch/throw speeds, catch angle, g-load,
     tip motion near touchdown).
  3. Elliptical-orbit LEO tether sizing (bottom tip matches LEO, top tip matches TLI).
  4. One-way LEO <-> lunar-surface delta-v and propellant, rockets vs tethers.

All numbers are two-body / impulsive approximations for trade-study use.
"""
import math

MU_E = 398600.4418   # km^3/s^2
MU_M = 4902.8        # km^3/s^2
R_E = 6371.0         # km
R_M = 1737.4         # km
D_EM = 384400.0      # km
G0 = 9.80665e-3      # km/s^2
ISP_CHEM = 370.0     # s, methalox

FIBERS = {  # name: (tensile strength Pa, density kg/m^3)
    "Zylon (PBO)": (5.8e9, 1560.0),
    "Dyneema SK99": (4.1e9, 975.0),
    "Kevlar 49": (3.6e9, 1440.0),
}


def char_velocity(sigma, rho, sf):
    """Fiber characteristic velocity Vc in km/s."""
    return math.sqrt(2.0 * sigma / (sf * rho)) / 1000.0


def mass_ratio(v_tip, vc):
    """Ideal tapered rotating-tether mass / payload mass."""
    x = v_tip / vc
    return math.sqrt(math.pi) * x * math.exp(x * x) * math.erf(x)


def prop_per_kg(dv):
    """kg of propellant per kg of mass delivered after a burn of dv (km/s)."""
    return math.exp(dv / (ISP_CHEM * G0)) - 1.0


def main():
    print("=== 1. Tether mass ratio (tether/payload) ===")
    speeds = (1.0, 1.54, 1.61, 2.4, 3.3, 3.8)
    for sf in (2, 3):
        print(f"Safety factor {sf}")
        for name, (s, rho) in FIBERS.items():
            vc = char_velocity(s, rho, sf)
            row = "  ".join(f"{v:.2f}:{mass_ratio(v, vc):8.1f}" for v in speeds)
            print(f"  {name:13s} Vc={vc:.2f} km/s | {row}")

    print("\n=== 2. Full lunar rotovator (150 km arm, center at 150 km) ===")
    arm, h = 150.0, 150.0
    r = R_M + h
    v_orb = math.sqrt(MU_M / r)
    v_tip = v_orb
    top = v_orb + v_tip
    g_load = v_tip ** 2 / arm / G0
    print(f"  orbit {v_orb:.3f} km/s | bottom tip 0 | top tip {top:.3f} km/s | {g_load:.2f} g")
    # Earth-transfer arrival (Hohmann-like, v_inf ~0.83 km/s)
    r_leo = R_E + 400.0
    v_leo = math.sqrt(MU_E / r_leo)
    v_tli = math.sqrt(MU_E * (2 / r_leo - 2 / (r_leo + D_EM)))
    v_inf = 0.83
    v_arr = math.sqrt(v_inf ** 2 + 2 * MU_M / (R_M + 300.0))
    ang = 2 * math.degrees(math.acos(v_arr / top))
    print(f"  arrival at 300 km: {v_arr:.3f} km/s -> catch {ang:.1f} deg before top "
          f"(tip flight-path angle {ang / 2:.1f} deg)")
    a_up = v_tip ** 2 / arm * 1000 - MU_M / r ** 2 * 1000  # m/s^2
    print(f"  near touchdown the tip decelerates/accelerates vertically at {a_up:.1f} m/s^2")
    for t in (0.5, 1, 2, 5, 10):
        print(f"    t=+/-{t:>4}s  height {0.5 * a_up * t * t:7.1f} m  vz {a_up * t:6.1f} m/s")

    print("\n=== 3. Elliptical-orbit LEO tether ===")
    v_peri = (v_tli + v_leo) / 2
    tip = (v_tli - v_leo) / 2
    a = 1 / (2 / r_leo - v_peri ** 2 / MU_E)
    print(f"  perigee 400 km at {v_peri:.3f} km/s, tip {tip:.3f} km/s, "
          f"apogee alt {2 * a - r_leo - R_E:,.0f} km")
    print(f"  bottom tip {v_peri - tip:.3f} (= LEO {v_leo:.3f}) | top tip {v_peri + tip:.3f} (= TLI {v_tli:.3f})")

    print("\n=== 4. One-way rocket delta-v, LEO <-> lunar surface ===")
    r_llo = R_M + 100.0
    loi = math.sqrt(v_inf ** 2 + 2 * MU_M / r_llo) - math.sqrt(MU_M / r_llo)
    descent = 1.87  # LLO -> surface incl. gravity losses
    tli = v_tli - v_leo
    cases = {
        "Rockets only (propulsive LEO capture)": tli + loi + descent,
        "Rockets, aerobrake at Earth (return leg)": descent + loi + 0.1,
        "Lunar rotovator only": tli + 0.1,
        "Lunar rotovator + aerobrake (return leg)": 0.15,
        "Rotovator + elliptical LEO tether": 0.2,
    }
    print(f"  TLI {tli:.3f} | LOI {loi:.3f} | descent {descent:.2f} km/s")
    for k, dv in cases.items():
        print(f"  {k:42s} {dv:5.2f} km/s  -> {prop_per_kg(dv):5.2f} kg prop per kg delivered")


if __name__ == "__main__":
    main()
