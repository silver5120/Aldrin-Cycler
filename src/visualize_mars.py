"""
Visualization & Animation Generator for the Earth-Mars Aldrin Cycler.

Generates:
  1. assets/mars_cycler_analysis.png        - 4-panel quantitative astrodynamics & entry/stability dashboard
  2. assets/mars_cycler_animation.gif       - Dual-frame (Heliocentric + Sun-Earth Rotating) twin-cycler animation
"""

from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Circle

BG_DARK = "#0b1020"
PANEL_BG = "#11182c"
GRID_COL = "#23304d"
TEXT_MAIN = "#e6edf7"
TEXT_MUTED = "#94a3b8"

C_FAST = "#38bdf8"       # cyan (145 d outbound transit)
C_RETURN = "#f59e0b"     # amber (635 d deep-space return)
C_INBOUND = "#34d399"    # emerald (inbound down-escalator cycler)
C_EARTH = "#3b82f6"
C_MARS = "#ef4444"
C_SUN = "#fbbf24"
C_PURPLE = "#a78bfa"
C_ROSE = "#fb7185"


def _style_ax(ax, title: str, xlabel: str = "", ylabel: str = ""):
    ax.set_facecolor(PANEL_BG)
    ax.set_title(title, color=TEXT_MAIN, fontsize=11.2, fontweight="bold", pad=9)
    if xlabel:
        ax.set_xlabel(xlabel, color=TEXT_MUTED, fontsize=9.3)
    if ylabel:
        ax.set_ylabel(ylabel, color=TEXT_MUTED, fontsize=9.3)
    ax.tick_params(colors=TEXT_MUTED, labelsize=8.4)
    for spine in ax.spines.values():
        spine.set_color(GRID_COL)
        spine.set_linewidth(0.9)
    ax.grid(True, color=GRID_COL, alpha=0.45, linewidth=0.6)


def generate_mars_analysis_figure(data: dict, out_png: Path):
    """Create the 4-panel high-resolution Earth-Mars Aldrin Cycler engineering figure."""
    id_m = data["ideal"]
    cmp = data["campaign_15yr"]
    trj = data["trajectory"]
    ts = trj["mars_trade_study"]

    fig, axes = plt.subplots(2, 2, figsize=(14.8, 11.4), facecolor=BG_DARK)
    fig.subplots_adjust(left=0.065, right=0.96, top=0.90, bottom=0.075, wspace=0.25, hspace=0.30)

    fig.suptitle(
        "Earth–Mars Aldrin Cycler — Heliocentric Precession, Taxi Δv, Mars Entry & 15-Year Stability",
        color=TEXT_MAIN,
        fontsize=15.0,
        fontweight="bold",
        y=0.965,
    )
    fig.text(
        0.5,
        0.932,
        f"a = {id_m['a_au']:.3f} AU, e = {id_m['e']:.3f}  |  Fast Transit: {id_m['tof_fast_days']:.1f} d  |  "
        f"LEO Taxi Δv: {id_m['dv_leo_taxi_ideal_km_s']:.2f} km/s  |  "
        f"Mars Entry: {id_m['v_ei_mars_cycler_km_s']:.2f} km/s  |  "
        f"15-Yr Avg Correction: {cmp['avg_dv_corr_per_cycle_m_s']:.0f} m/s/cycle",
        ha="center",
        color="#7dd3fc",
        fontsize=10.0,
    )

    # ------------------------------------------------------------------
    # Panel 1: Heliocentric 3-Cycle Precession (6.4 Years)
    # ------------------------------------------------------------------
    ax = axes[0, 0]
    _style_ax(
        ax,
        "1. Heliocentric Precession (+48.7°/Cycle via Earth Gravity Assist)",
        "X (AU)",
        "Y (AU)",
    )
    th = np.linspace(0, 2 * np.pi, 300)
    ax.plot(np.cos(th), np.sin(th), color=C_EARTH, alpha=0.45, linewidth=1.2, linestyle="--", label="Earth Orbit (1.00 AU)")
    ax.plot(1.5237 * np.cos(th), 1.5237 * np.sin(th), color=C_MARS, alpha=0.45, linewidth=1.2, linestyle="--", label="Mars Orbit (1.52 AU)")

    ox = np.array(trj["out_x"])
    oy = np.array(trj["out_y"])
    ph = np.array(trj["out_phase"])
    pts_per_cyc = len(ox) // trj["n_cycles_anim"]

    # Plot Cycles 2 & 3 muted to show +48.7 deg rosette precession
    ax.plot(ox[pts_per_cyc:], oy[pts_per_cyc:], color="#64748b", linewidth=1.25, alpha=0.75, label="Cycles 2–3 (+48.7° Precession)")

    # Highlight Cycle 1 fast outbound vs return loop
    ox1, oy1, ph1 = ox[:pts_per_cyc], oy[:pts_per_cyc], ph[:pts_per_cyc]
    fast_idx = np.where(ph1 == 0)[0]
    ret_idx = np.where(ph1 == 1)[0]
    ax.plot(ox1[fast_idx], oy1[fast_idx], color=C_FAST, linewidth=2.5, label=f"Fast Earth→Mars Leg ({id_m['tof_fast_days']:.0f} d)")
    ax.plot(ox1[ret_idx], oy1[ret_idx], color=C_RETURN, linewidth=1.9, alpha=0.9, label=f"Aphelion Return Leg ({id_m['tof_long_days']:.0f} d, 2.22 AU)")

    # Earth departure & Mars arrival markers for Cycle 1
    m_arr = fast_idx[-1]
    ax.scatter([trj["earth_x"][0]], [trj["earth_y"][0]], color=C_EARTH, s=80, zorder=8, edgecolors="white", lw=1.0)
    ax.scatter([trj["mars_x"][m_arr]], [trj["mars_y"][m_arr]], color=C_MARS, s=75, zorder=8, edgecolors="white", lw=1.0)
    ax.text(trj["earth_x"][0] + 0.12, trj["earth_y"][0] - 0.18, "Earth GA\n(Day 0)", color="#93c5fd", fontsize=7.8, fontweight="bold")
    ax.text(trj["mars_x"][m_arr] - 0.15, trj["mars_y"][m_arr] + 0.16, "Mars Encounter\n(Day 145)", color="#fca5a5", fontsize=7.8, fontweight="bold", ha="center")

    ax.add_patch(Circle((0, 0), 0.09, color=C_SUN, zorder=9))
    ax.text(0, -0.22, "Sun", color=C_SUN, fontsize=8.0, ha="center", fontweight="bold")

    ax.set_aspect("equal")
    ax.set_xlim(-2.45, 2.45)
    ax.set_ylim(-2.45, 2.45)
    ax.legend(loc="lower left", fontsize=7.3, facecolor="#0f172a", edgecolor=GRID_COL, labelcolor=TEXT_MAIN, framealpha=0.92)

    # ------------------------------------------------------------------
    # Panel 2: Taxi LEO Delta-V & Mars Entry / Capture Comparison
    # ------------------------------------------------------------------
    ax = axes[0, 1]
    _style_ax(
        ax,
        "2. Propulsion & Entry Demands: Aldrin Cycler vs. Hohmann",
        "",
        "Velocity / Δv (km/s)",
    )
    categories = [
        "LEO (400 km)\nTaxi Escape Δv",
        "Mars Entry\nSpeed v_EI",
        "Mars Propulsive\n400 km LMO Δv",
        "Mars Aerocapture\nTrim Δv",
    ]
    cyc_vals = [
        id_m["dv_leo_taxi_ideal_km_s"],
        id_m["v_ei_mars_cycler_km_s"],
        id_m["dv_mars_propulsive_lmo_km_s"],
        id_m["dv_mars_aerocapture_circ_km_s"],
    ]
    hoh_vals = [
        id_m["dv_leo_hohmann_km_s"],
        id_m["v_ei_mars_hohmann_km_s"],
        id_m["dv_mars_propulsive_lmo_hoh_km_s"],
        0.14,
    ]
    x = np.arange(len(categories))
    w = 0.35
    b1 = ax.bar(x - w / 2, cyc_vals, width=w, color=C_FAST, edgecolor="#0f172a", label="Aldrin Cycler (145 d)")
    b2 = ax.bar(x + w / 2, hoh_vals, width=w, color=C_RETURN, edgecolor="#0f172a", alpha=0.88, label="Hohmann Transfer (259 d)")

    for bar in b1:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + 0.18, f"{h:.2f}", ha="center", va="bottom", color="#bae6fd", fontsize=8.4, fontweight="bold")
    for bar in b2:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + 0.18, f"{h:.2f}", ha="center", va="bottom", color="#fde68a", fontsize=8.4, fontweight="bold")

    ax.annotate(
        f"3.76× Kinetic Energy\n7.3× Conv. Heating\n(Aerocapture mandatory!)",
        xy=(1 - w / 2, cyc_vals[1]),
        xytext=(1.55, 9.6),
        color="#fca5a5",
        fontsize=8.0,
        fontweight="bold",
        arrowprops=dict(arrowstyle="->", color=C_ROSE, lw=1.2),
        bbox=dict(boxstyle="round,pad=0.3", facecolor="#1e293b", edgecolor=C_ROSE, alpha=0.92),
    )

    ax.set_xticks(x)
    ax.set_xticklabels(categories, color=TEXT_MAIN, fontsize=8.4)
    ax.set_ylim(0, 12.4)
    ax.legend(loc="upper right", fontsize=8.0, facecolor="#0f172a", edgecolor=GRID_COL, labelcolor=TEXT_MAIN, framealpha=0.92)

    # ------------------------------------------------------------------
    # Panel 3: 15-Year Eccentric Mars Correction Burns (7 Synodic Cycles)
    # ------------------------------------------------------------------
    ax = axes[1, 0]
    _style_ax(
        ax,
        "3. 15-Year Eccentric Mars (e = 0.0934) Correction Burns per Cycle",
        "Synodic Cycle (Year Start)",
        "Correction Δv per Cycle (m/s)",
    )
    cycs = cmp["cycles"]
    cx = np.arange(len(cycs))
    dv_e = np.array([c["dv_earth_flyby_m_s"] for c in cycs])
    dv_m = np.array([c["dv_mars_corr_m_s"] for c in cycs])
    dv_i = np.array([c["dv_inclination_m_s"] for c in cycs])
    dv_t = np.array([c["dv_nav_tcm_m_s"] for c in cycs])
    r_m = np.array([c["mars_radius_au"] for c in cycs])

    ax.bar(cx, dv_e, width=0.55, color=C_FAST, label="Earth Powered Flyby / Turn Deficit")
    ax.bar(cx, dv_m, bottom=dv_e, width=0.55, color=C_MARS, label="Mars GA Trim")
    ax.bar(cx, dv_i, bottom=dv_e + dv_m, width=0.55, color=C_PURPLE, label="Mars Inclination (1.85°) Trim")
    ax.bar(cx, dv_t, bottom=dv_e + dv_m + dv_i, width=0.55, color=C_INBOUND, label="Nav TCMs (~35 m/s)")

    ax2 = ax.twinx()
    ax2.set_ylabel("Mars Distance at Encounter (AU)", color=C_RETURN, fontsize=9.2)
    ax2.tick_params(axis="y", colors=C_RETURN, labelsize=8.4)
    ax2.spines["right"].set_color(GRID_COL)
    ax2.spines["top"].set_visible(False)
    l_rm, = ax2.plot(cx, r_m, color=C_RETURN, marker="o", linewidth=2.0, markersize=5.5, label="Mars r_M (AU)")
    ax2.set_ylim(1.28, 1.76)

    ax.set_xticks(cx)
    ax.set_xticklabels([f"C{c['cycle']}\n(Yr {c['year_start']:.1f})" for c in cycs], fontsize=8.0)
    ax.set_ylim(0, max(dv_e + dv_m + dv_i + dv_t) * 1.28)
    h1, lb1 = ax.get_legend_handles_labels()
    ax.legend(h1 + [l_rm], lb1 + ["Mars Distance (AU)"], loc="upper right", fontsize=7.3, facecolor="#0f172a", edgecolor=GRID_COL, labelcolor=TEXT_MAIN, framealpha=0.92)

    # ------------------------------------------------------------------
    # Panel 4: Earth-Mars Cumulative IMLEO Trade Study
    # ------------------------------------------------------------------
    ax = axes[1, 1]
    _style_ax(
        ax,
        "4. Why Mars Cyclers Pay Off Faster Than Lunar Cyclers",
        "Number of Crewed Mars Expeditions (1 every 2.14 yrs)",
        "Cumulative Initial Mass in LEO (metric tons)",
    )
    fl = np.array(ts["flights"])
    c_dir = np.array(ts["cum_direct_t"])
    c_cyc = np.array(ts["cum_cycler_t"])

    ax.plot(fl, c_dir, color=C_ROSE, linewidth=2.4, marker="o", markersize=4.5,
            label=f"Direct Mars (50 t Hab + 15 t Lander = {ts['imleo_direct_per_flight_t']:.0f} t/flight)")
    ax.plot(fl, c_cyc, color=C_INBOUND, linewidth=2.5, marker="^", markersize=5.0,
            label=f"Aldrin Cycler ({ts['imleo_cycler_setup_t']:.0f} t setup + {ts['imleo_cycler_per_flight_t']:.0f} t/flight taxi)")

    be = ts["breakeven_flights"]
    y_be = ts["imleo_direct_per_flight_t"] * be
    ax.axvline(be, color="#fbbf24", linestyle=":", linewidth=1.3)
    ax.scatter([be], [y_be], color="#fbbf24", s=75, zorder=8, edgecolors="white")
    ax.annotate(
        f"Break-even: ~{be:.1f} flights!\n(50 t radiation-shielded hab\nis mandatory for 145 d transit;\n{ts['savings_at_10_flights_pct']:.0f}% IMLEO saved by Flight 10)",
        xy=(be, y_be),
        xytext=(be + 1.8, y_be * 2.15),
        color="#fde68a",
        fontsize=8.2,
        fontweight="bold",
        arrowprops=dict(arrowstyle="->", color="#fbbf24", lw=1.2),
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#1e293b", edgecolor="#fbbf24", alpha=0.92),
    )

    ax.set_xlim(1, 12)
    ax.set_ylim(0, max(c_dir) * 1.08)
    ax.legend(loc="upper left", fontsize=7.6, facecolor="#0f172a", edgecolor=GRID_COL, labelcolor=TEXT_MAIN, framealpha=0.92)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180, facecolor=BG_DARK)
    plt.close(fig)
    print(f"Saved {out_png}")


def generate_mars_animation_gif(data: dict, out_gif: Path, n_frames: int = 140, fps: int = 18):
    """Create a smooth dual-frame (Heliocentric + Sun-Earth Rotating) Earth-Mars Twin Cycler GIF."""
    id_m = data["ideal"]
    trj = data["trajectory"]

    t_days = np.array(trj["t_days"])
    t_total = t_days[-1]
    t_u = np.linspace(0.0, t_total, n_frames, endpoint=False)

    ex = np.interp(t_u, t_days, trj["earth_x"])
    ey = np.interp(t_u, t_days, trj["earth_y"])
    mx = np.interp(t_u, t_days, trj["mars_x"])
    my = np.interp(t_u, t_days, trj["mars_y"])
    ox = np.interp(t_u, t_days, trj["out_x"])
    oy = np.interp(t_u, t_days, trj["out_y"])
    orx = np.interp(t_u, t_days, trj["out_rot_x"])
    ory = np.interp(t_u, t_days, trj["out_rot_y"])
    ix = np.interp(t_u, t_days, trj["in_x"])
    iy = np.interp(t_u, t_days, trj["in_y"])
    irx = np.interp(t_u, t_days, trj["in_rot_x"])
    iry = np.interp(t_u, t_days, trj["in_rot_y"])
    ospd = np.interp(t_u, t_days, trj["out_speed_km_s"])
    ode = np.interp(t_u, t_days, trj["out_dist_earth_au"])
    odm = np.interp(t_u, t_days, trj["out_dist_mars_au"])

    T_syn = id_m["T_synodic_days"]
    tof_fast = id_m["tof_fast_days"]

    fig, (ax_hel, ax_rot) = plt.subplots(1, 2, figsize=(11.8, 5.9), facecolor=BG_DARK)
    fig.subplots_adjust(left=0.055, right=0.975, top=0.82, bottom=0.14, wspace=0.14)

    fig.suptitle(
        "Earth–Mars Twin Aldrin Cyclers (3 Synodic Cycles / 6.4 Years)",
        color=TEXT_MAIN,
        fontsize=13.8,
        fontweight="bold",
        y=0.955,
    )
    hud_txt = fig.text(0.5, 0.885, "", ha="center", color=C_FAST, fontsize=10.0, family="monospace")
    footer_txt = fig.text(0.5, 0.042, "", ha="center", color=TEXT_MUTED, fontsize=9.2, family="monospace")

    # Left: Heliocentric Inertial Frame
    _style_ax(ax_hel, "Heliocentric Inertial Frame (AU)", "X (AU)", "Y (AU)")
    th = np.linspace(0, 2 * np.pi, 240)
    ax_hel.plot(np.cos(th), np.sin(th), color=C_EARTH, alpha=0.35, linestyle=":", linewidth=1.2)
    ax_hel.plot(1.5237 * np.cos(th), 1.5237 * np.sin(th), color=C_MARS, alpha=0.35, linestyle=":", linewidth=1.2)
    ax_hel.plot(trj["out_x"], trj["out_y"], color="#475569", alpha=0.32, linewidth=1.1)
    ax_hel.add_patch(Circle((0, 0), 0.08, color=C_SUN, zorder=6))
    ax_hel.set_aspect("equal")
    ax_hel.set_xlim(-2.4, 2.4)
    ax_hel.set_ylim(-2.4, 2.4)

    # Right: Sun-Earth Rotating Frame (Earth fixed at (1, 0))
    _style_ax(ax_rot, "Sun–Earth Rotating Frame (Earth Fixed)", "x_rot (AU)", "y_rot (AU)")
    ax_rot.plot(np.cos(th), np.sin(th), color=C_EARTH, alpha=0.25, linestyle=":", linewidth=1.0)
    ax_rot.plot(1.5237 * np.cos(th), 1.5237 * np.sin(th), color=C_MARS, alpha=0.35, linestyle=":", linewidth=1.2)
    ax_rot.plot(trj["out_rot_x"], trj["out_rot_y"], color="#475569", alpha=0.32, linewidth=1.1)
    ax_rot.add_patch(Circle((0, 0), 0.08, color=C_SUN, zorder=6))
    ax_rot.add_patch(Circle((1.0, 0.0), 0.065, color=C_EARTH, zorder=7))
    ax_rot.text(1.0, -0.22, "Earth", color="#93c5fd", fontsize=8.0, ha="center", fontweight="bold")
    ax_rot.set_aspect("equal")
    ax_rot.set_xlim(-2.4, 2.4)
    ax_rot.set_ylim(-2.4, 2.4)

    # Animated artists
    out_trail_h, = ax_hel.plot([], [], color=C_FAST, linewidth=2.5, zorder=7)
    out_sc_h, = ax_hel.plot([], [], marker="o", markersize=8.0, color="white", markeredgecolor=C_FAST, markeredgewidth=2.0, zorder=9)
    in_sc_h, = ax_hel.plot([], [], marker="s", markersize=6.5, color=C_INBOUND, markeredgecolor="white", markeredgewidth=1.0, zorder=8)
    earth_h, = ax_hel.plot([], [], marker="o", markersize=9.5, color=C_EARTH, zorder=8)
    mars_h, = ax_hel.plot([], [], marker="o", markersize=8.0, color=C_MARS, zorder=8)
    earth_lbl_h = ax_hel.text(0, 0, "Earth", color="#93c5fd", fontsize=7.6, ha="center", fontweight="bold")
    mars_lbl_h = ax_hel.text(0, 0, "Mars", color="#fca5a5", fontsize=7.6, ha="center", fontweight="bold")

    out_trail_r, = ax_rot.plot([], [], color=C_FAST, linewidth=2.5, zorder=7)
    out_sc_r, = ax_rot.plot([], [], marker="o", markersize=8.0, color="white", markeredgecolor=C_FAST, markeredgewidth=2.0, zorder=9)
    in_sc_r, = ax_rot.plot([], [], marker="s", markersize=6.5, color=C_INBOUND, markeredgecolor="white", markeredgewidth=1.0, zorder=8)
    mars_r, = ax_rot.plot([], [], marker="o", markersize=8.0, color=C_MARS, zorder=8)

    trail_len = 24

    def update(frame: int):
        t_cur = t_u[frame]
        cyc_num = int(t_cur // T_syn) + 1
        dt_cyc = t_cur % T_syn
        is_fast = dt_cyc <= tof_fast
        col = C_FAST if is_fast else C_RETURN

        i0 = max(0, frame - trail_len)
        sl = slice(i0, frame + 1)

        out_trail_h.set_data(ox[sl], oy[sl])
        out_trail_h.set_color(col)
        out_sc_h.set_data([ox[frame]], [oy[frame]])
        out_sc_h.set_markeredgecolor(col)
        in_sc_h.set_data([ix[frame]], [iy[frame]])
        earth_h.set_data([ex[frame]], [ey[frame]])
        mars_h.set_data([mx[frame]], [my[frame]])
        earth_lbl_h.set_position((ex[frame] * 1.18, ey[frame] * 1.18))
        mars_lbl_h.set_position((mx[frame] * 1.13, my[frame] * 1.13))

        # Mars in Sun-Earth rotating frame
        th_e = np.arctan2(ey[frame], ex[frame])
        mrx = mx[frame] * np.cos(-th_e) - my[frame] * np.sin(-th_e)
        mry = mx[frame] * np.sin(-th_e) + my[frame] * np.cos(-th_e)

        out_trail_r.set_data(orx[sl], ory[sl])
        out_trail_r.set_color(col)
        out_sc_r.set_data([orx[frame]], [ory[frame]])
        out_sc_r.set_markeredgecolor(col)
        in_sc_r.set_data([irx[frame]], [iry[frame]])
        mars_r.set_data([mrx], [mry])

        tag = "OUTBOUND CREW SPRINT (Earth → Mars, 145 d)" if is_fast else "APHELION RETURN LOOP (Uncrewed, 2.22 AU)"
        if dt_cyc < 22.0 or dt_cyc > T_syn - 22.0:
            tag = "★ EARTH GRAVITY ASSIST (+48.7° Precession) & TAXI DOCKING"
            col_hud = "#93c5fd"
        elif abs(dt_cyc - tof_fast) < 18.0:
            tag = "★ MARS ENCOUNTER (v_EI = 10.9 km/s — Taxi Aerocapture!)"
            col_hud = "#fca5a5"
        else:
            col_hud = col

        hud_txt.set_text(f"Year {t_cur / 365.25:4.2f}  |  Cycle {cyc_num} (Day {dt_cyc:5.1f})  |  {tag}")
        hud_txt.set_color(col_hud)
        footer_txt.set_text(
            f"Helio Speed: {ospd[frame]:5.2f} km/s   |   Dist to Earth: {ode[frame]:5.2f} AU   |   "
            f"Dist to Mars: {odm[frame]:5.2f} AU   |   ■ Green = Inbound Cycler"
        )
        return out_trail_h, out_sc_h, in_sc_h, earth_h, mars_h, earth_lbl_h, mars_lbl_h, out_trail_r, out_sc_r, in_sc_r, mars_r, hud_txt, footer_txt

    anim = FuncAnimation(fig, update, frames=n_frames, interval=1000 // fps, blit=False)
    out_gif.parent.mkdir(parents=True, exist_ok=True)
    anim.save(out_gif, writer=PillowWriter(fps=fps), dpi=95)
    plt.close(fig)
    print(f"Saved {out_gif}")


def main():
    root = Path(__file__).resolve().parent.parent
    data_path = root / "data" / "mars_cycler_solution.json"
    if not data_path.exists():
        from mars_cycler import export_mars_solution
        export_mars_solution(data_path)
    data = json.loads(data_path.read_text(encoding="utf-8"))
    generate_mars_analysis_figure(data, root / "assets" / "mars_cycler_analysis.png")
    generate_mars_animation_gif(data, root / "assets" / "mars_cycler_animation.gif")


if __name__ == "__main__":
    main()
