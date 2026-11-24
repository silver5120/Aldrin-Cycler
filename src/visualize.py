"""
Visualization & Animation Generator for the Earth-Moon Aldrin Cycler.

Generates:
  1. assets/cycler_analysis.png        - 4-panel quantitative astrodynamics summary
  2. assets/aldrin_cycler_animation.gif - Dual-frame (Rotating + Inertial) orbit animation
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

from cr3bp import L_STAR_KM, MU, R_EARTH_KM, R_MOON_KM, pseudo_potential


BG_DARK = "#0b1020"
PANEL_BG = "#11182c"
GRID_COL = "#23304d"
TEXT_MAIN = "#e6edf7"
TEXT_MUTED = "#94a3b8"

C_OUTBOUND = "#38bdf8"   # cyan
C_INBOUND = "#34d399"    # emerald
C_PHASE1 = "#f59e0b"     # amber
C_PHASE2 = "#a78bfa"     # violet
C_ARENSTORF = "#fb7185"  # rose
C_EARTH = "#3b82f6"
C_MOON = "#cbd5e1"
C_BURN = "#ef4444"


def _style_ax(ax, title: str, xlabel: str = "", ylabel: str = ""):
    ax.set_facecolor(PANEL_BG)
    ax.set_title(title, color=TEXT_MAIN, fontsize=11.5, fontweight="bold", pad=9)
    if xlabel:
        ax.set_xlabel(xlabel, color=TEXT_MUTED, fontsize=9.5)
    if ylabel:
        ax.set_ylabel(ylabel, color=TEXT_MUTED, fontsize=9.5)
    ax.tick_params(colors=TEXT_MUTED, labelsize=8.5)
    for spine in ax.spines.values():
        spine.set_color(GRID_COL)
        spine.set_linewidth(0.9)
    ax.grid(True, color=GRID_COL, alpha=0.45, linewidth=0.6)


def generate_analysis_figure(data: dict, out_png: Path):
    """Create the 4-panel high-resolution astrodynamics summary plot."""
    sh = data["shamrock"]
    ar = data["arenstorf"]
    tr = data["trade_study"]
    lp = data["lagrange_points"]
    sm = sh["metrics"]
    am = ar["metrics"]

    fig, axes = plt.subplots(2, 2, figsize=(14.5, 11.2), facecolor=BG_DARK)
    fig.subplots_adjust(left=0.065, right=0.96, top=0.90, bottom=0.07, wspace=0.25, hspace=0.28)

    fig.suptitle(
        "Earth–Moon Aldrin–Genova Shamrock (C-2-R) Cycler — CR3BP Numerical Analysis",
        color=TEXT_MAIN,
        fontsize=15.5,
        fontweight="bold",
        y=0.965,
    )
    fig.text(
        0.5,
        0.932,
        f"Cycle Period: {sm['period_days']:.2f} d  |  Earth→Moon Sprint: {sm['earth_to_moon_transit_days']:.2f} d  |  "
        f"Station-Keeping Δv: {sm['dv_total_per_cycle_m_s']:.1f} m/s/cycle  |  "
        f"Perigee / Perilune: {sm['perigee_alt_km']:,.0f} / {sm['perilune_alt_km']:,.0f} km",
        ha="center",
        color="#7dd3fc",
        fontsize=10.2,
    )

    # ------------------------------------------------------------------
    # Panel 1: Earth-Moon Rotating (Synodic) Frame
    # ------------------------------------------------------------------
    ax = axes[0, 0]
    _style_ax(
        ax,
        "1. Earth–Moon Rotating (Synodic) Frame",
        "x / L*  (1 L* = 384,400 km)",
        "y / L*",
    )

    # Subtle CR3BP effective potential contours
    gx = np.linspace(-1.35, 1.35, 220)
    gy = np.linspace(-1.25, 1.25, 220)
    GX, GY = np.meshgrid(gx, gy)
    CJ_grid = 2.0 * pseudo_potential(GX, GY, MU)
    ax.contour(
        GX,
        GY,
        CJ_grid,
        levels=[2.95, 3.05, 3.188],
        colors=["#334155"],
        linewidths=0.65,
        linestyles=":",
        alpha=0.65,
    )

    # Pure Ballistic Arenstorf comparison
    ax.plot(
        ar["x_rot"],
        ar["y_rot"],
        color=C_ARENSTORF,
        linestyle="--",
        linewidth=1.2,
        alpha=0.65,
        label=f"Ballistic Arenstorf (Δv=0, min hE={am['min_earth_alt_km']/1e3:.0f}k km)",
    )

    # Shamrock segments color-coded
    xr = np.array(sh["x_rot"])
    yr = np.array(sh["y_rot"])
    ph = np.array(sh["phase"])
    seg_meta = [
        (0, C_OUTBOUND, f"Outbound Fig-8 ({sm['earth_to_moon_transit_days']:.1f} d)"),
        (1, C_INBOUND, f"Inbound Fig-8 ({sm['earth_to_moon_transit_days']:.1f} d)"),
        (2, C_PHASE1, f"Phasing Loop 1 ({sm['phasing_loop_days']:.1f} d)"),
        (3, C_PHASE2, f"Phasing Loop 2 ({sm['phasing_loop_days']:.1f} d)"),
    ]
    for code, col, lbl in seg_meta:
        idx = np.where(ph == code)[0]
        if len(idx) > 1:
            end_i = min(idx[-1] + 2, len(xr))
            ax.plot(xr[idx[0]:end_i], yr[idx[0]:end_i], color=col, linewidth=2.1, label=lbl)

    # Burn markers at perigee
    s_p_out = sm["state_perigee_out"]
    ax.scatter(
        [s_p_out[0], s_p_out[0]],
        [s_p_out[1], -s_p_out[1]],
        color=C_BURN,
        s=52,
        zorder=8,
        marker="D",
        edgecolors="white",
        linewidths=0.8,
        label=f"Perigee Burns (±{sm['dv_burn1_m_s']:.1f} m/s)",
    )

    # Lagrange points
    for name, (lx, ly) in lp.items():
        ax.scatter([lx], [ly], color="#64748b", marker="+", s=40, zorder=5)
        dy = 0.06 if ly >= 0 else -0.10
        ax.text(lx, ly + dy, name, color="#94a3b8", fontsize=7.8, ha="center")

    # Earth & Moon
    ax.add_patch(Circle((-MU, 0.0), 0.042, color=C_EARTH, zorder=7))
    ax.add_patch(Circle((1.0 - MU, 0.0), 0.024, color=C_MOON, zorder=7))
    ax.text(-MU - 0.04, -0.11, "Earth", color="#93c5fd", fontsize=8.5, ha="center", fontweight="bold")
    ax.text(1.0 - MU, -0.10, "Moon", color="#e2e8f0", fontsize=8.5, ha="center", fontweight="bold")

    ax.set_aspect("equal")
    ax.set_xlim(-1.25, 1.28)
    ax.set_ylim(-1.12, 1.12)
    ax.legend(
        loc="lower left",
        fontsize=7.3,
        facecolor="#0f172a",
        edgecolor=GRID_COL,
        labelcolor=TEXT_MAIN,
        framealpha=0.9,
    )

    # ------------------------------------------------------------------
    # Panel 2: Earth-Centered Inertial (ECI) Frame (3 Cycles)
    # ------------------------------------------------------------------
    ax = axes[0, 1]
    _style_ax(
        ax,
        "2. Earth-Centered Inertial (ECI) Frame — 3 Cycles (78 Days)",
        "X_ECI / L*",
        "Y_ECI / L*",
    )

    theta = np.linspace(0, 2 * np.pi, 300)
    ax.plot(np.cos(theta), np.sin(theta), color="#64748b", linestyle=":", linewidth=1.2, label="Lunar Orbit (384,400 km)")

    # Plot 3 cycles lightly, then highlight cycle 1 brightly
    ax.plot(
        sh["multi_x_eci"],
        sh["multi_y_eci"],
        color="#475569",
        linewidth=1.1,
        alpha=0.75,
        label="Cycles 2–3 (Inertial Precession)",
    )
    xe = np.array(sh["x_eci"])
    ye = np.array(sh["y_eci"])
    for code, col, _ in seg_meta:
        idx = np.where(ph == code)[0]
        if len(idx) > 1:
            end_i = min(idx[-1] + 2, len(xe))
            lbl = "Cycle 1 (Fig-8 + 2 Phasing)" if code == 0 else None
            ax.plot(xe[idx[0]:end_i], ye[idx[0]:end_i], color=col, linewidth=2.0, label=lbl)

    # Mark Moon position at Perilune encounter of Cycle 1
    p_idx = int(np.argmin(sh["alt_moon_km"]))
    ax.scatter(
        [sh["xm_eci"][p_idx]],
        [sh["ym_eci"][p_idx]],
        color=C_MOON,
        s=95,
        zorder=8,
        edgecolors="#38bdf8",
        linewidths=1.5,
        label=f"Lunar Encounter (Day {sh['t_days'][p_idx]:.1f})",
    )
    ax.add_patch(Circle((0.0, 0.0), 0.045, color=C_EARTH, zorder=8))
    ax.text(0.0, -0.11, "Earth", color="#93c5fd", fontsize=8.5, ha="center", fontweight="bold")

    ax.set_aspect("equal")
    ax.set_xlim(-1.25, 1.25)
    ax.set_ylim(-1.25, 1.25)
    ax.legend(
        loc="lower left",
        fontsize=7.6,
        facecolor="#0f172a",
        edgecolor=GRID_COL,
        labelcolor=TEXT_MAIN,
        framealpha=0.9,
    )

    # ------------------------------------------------------------------
    # Panel 3: Altitude & Inertial Speed Telemetry (1 Cycle)
    # ------------------------------------------------------------------
    ax = axes[1, 0]
    _style_ax(
        ax,
        "3. Altitude & Inertial Speed Over One Cycle (25.96 Days)",
        "Mission Elapsed Time (days)",
        "Altitude (1,000 km)",
    )
    t_d = np.array(sh["t_days"])
    alt_e = np.array(sh["alt_earth_km"]) / 1000.0
    alt_m = np.array(sh["alt_moon_km"]) / 1000.0
    v_eci = np.array(sh["v_eci_km_s"])

    # Shaded phase backgrounds
    ax.axvspan(0, sm["fig8_leg_days"], color=C_OUTBOUND, alpha=0.09)
    ax.axvspan(sm["fig8_leg_days"], sm["period_days"], color=C_PHASE1, alpha=0.06)

    l1, = ax.plot(t_d, alt_e, color=C_OUTBOUND, linewidth=2.0, label="Earth Altitude (10³ km)")
    l2, = ax.plot(t_d, alt_m, color=C_MOON, linewidth=1.6, linestyle="--", label="Moon Altitude (10³ km)")

    ax2 = ax.twinx()
    ax2.set_ylabel("ECI Speed (km/s)", color=C_PHASE1, fontsize=9.5)
    ax2.tick_params(axis="y", colors=C_PHASE1, labelsize=8.5)
    ax2.spines["right"].set_color(GRID_COL)
    ax2.spines["top"].set_visible(False)
    l3, = ax2.plot(t_d, v_eci, color=C_PHASE1, linewidth=1.5, alpha=0.9, label="ECI Speed (km/s)")

    # Key event annotations
    ax.axvline(sm["perilune_time_days"], color=C_INBOUND, linestyle=":", linewidth=1.1)
    ax.axvline(sm["burn1_time_days"], color=C_BURN, linestyle=":", linewidth=1.2)
    ax.text(
        sm["perilune_time_days"] + 0.35,
        420,
        f"Perilune\n{sm['perilune_alt_km']:,.0f} km\n(Day {sm['perilune_time_days']:.1f})",
        color=C_INBOUND,
        fontsize=7.8,
        va="top",
        bbox=dict(boxstyle="round,pad=0.25", facecolor=PANEL_BG, edgecolor=C_INBOUND, alpha=0.88, lw=0.7),
    )
    ax.text(
        sm["burn1_time_days"] + 0.45,
        420,
        f"Burn 1: −{sm['dv_burn1_m_s']:.1f} m/s\n(Oberth Phasing)",
        color="#fca5a5",
        fontsize=7.8,
        va="top",
        bbox=dict(boxstyle="round,pad=0.25", facecolor=PANEL_BG, edgecolor=C_BURN, alpha=0.88, lw=0.7),
    )

    ax.set_xlim(0, sm["period_days"])
    ax.set_ylim(-10, 460)
    ax2.set_ylim(0, 10.8)
    lines = [l1, l2, l3]
    ax.legend(
        lines,
        [l.get_label() for l in lines],
        loc="upper right",
        fontsize=7.6,
        facecolor="#0f172a",
        edgecolor=GRID_COL,
        labelcolor=TEXT_MAIN,
        framealpha=0.9,
    )

    # ------------------------------------------------------------------
    # Panel 4: Tsiolkovsky Cumulative IMLEO Trade Study
    # ------------------------------------------------------------------
    ax = axes[1, 1]
    _style_ax(
        ax,
        "4. Cumulative IMLEO: Why We Don't Fly This Yet",
        "Number of Crewed Lunar Round-Trips",
        "Cumulative Initial Mass in LEO (metric tons)",
    )
    fl = np.array(tr["flights"])
    c_heavy = np.array(tr["cum_direct_imleo_t"])
    c_lean = np.array(tr["cum_apollo_lean_imleo_t"])
    c_cyc = np.array(tr["cum_cycler_imleo_t"])

    ax.plot(
        fl,
        c_heavy,
        color=C_ARENSTORF,
        linewidth=2.2,
        marker="o",
        markersize=3.5,
        label=f"Direct w/ 35 t Shielded Hab ({tr['imleo_direct_per_flight_t']:.0f} t/flight)",
    )
    ax.plot(
        fl,
        c_lean,
        color=C_PHASE1,
        linewidth=2.2,
        linestyle="--",
        marker="s",
        markersize=3.5,
        label=f"Lean Apollo/Orion 18 t Capsule ({tr['imleo_apollo_lean_t']:.0f} t/flight)",
    )
    ax.plot(
        fl,
        c_cyc,
        color=C_INBOUND,
        linewidth=2.5,
        marker="^",
        markersize=3.8,
        label=f"Aldrin Cycler ({tr['imleo_cycler_setup_t']:.0f} t setup + {tr['imleo_cycler_per_flight_t']:.0f} t/flight)",
    )

    # Shade break-even region vs lean Apollo/Orion
    be_lean = tr["breakeven_vs_lean_flights"]
    if np.isfinite(be_lean) and be_lean < 25:
        y_be = tr["imleo_apollo_lean_t"] * be_lean
        ax.axvline(be_lean, color="#fbbf24", linestyle=":", linewidth=1.3)
        ax.scatter([be_lean], [y_be], color="#fbbf24", s=70, zorder=8, edgecolors="white")
        ax.annotate(
            f"Break-even vs Lean Capsule:\n~{be_lean:.1f} flights\n(>10 yrs at Artemis cadence!)",
            xy=(be_lean, y_be),
            xytext=(min(be_lean + 2.2, 13.5), y_be * 0.52),
            color="#fde68a",
            fontsize=8.0,
            fontweight="bold",
            arrowprops=dict(arrowstyle="->", color="#fbbf24", lw=1.2),
            bbox=dict(boxstyle="round,pad=0.35", facecolor="#1e293b", edgecolor="#fbbf24", alpha=0.92),
        )

    ax.set_xlim(1, 25)
    ax.set_ylim(0, max(c_heavy) * 1.06)
    ax.legend(
        loc="upper left",
        fontsize=7.6,
        facecolor="#0f172a",
        edgecolor=GRID_COL,
        labelcolor=TEXT_MAIN,
        framealpha=0.92,
    )

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180, facecolor=BG_DARK)
    plt.close(fig)
    print(f"Saved {out_png}")


def generate_animation_gif(data: dict, out_gif: Path, n_frames: int = 180, fps: int = 20):
    """Create a smooth dual-frame (Rotating + ECI) animated GIF of the Shamrock Cycler."""
    sh = data["shamrock"]
    sm = sh["metrics"]
    lp = data["lagrange_points"]

    t_days = np.array(sh["t_days"])
    x_rot = np.array(sh["x_rot"])
    y_rot = np.array(sh["y_rot"])
    x_eci = np.array(sh["x_eci"])
    y_eci = np.array(sh["y_eci"])
    xm_eci = np.array(sh["xm_eci"])
    ym_eci = np.array(sh["ym_eci"])
    alt_e = np.array(sh["alt_earth_km"])
    alt_m = np.array(sh["alt_moon_km"])
    v_eci = np.array(sh["v_eci_km_s"])
    phase = np.array(sh["phase"])

    # Resample uniformly in mission elapsed time so orbital speed looks physically authentic!
    t_uniform = np.linspace(0.0, t_days[-1], n_frames, endpoint=False)
    xr_f = np.interp(t_uniform, t_days, x_rot)
    yr_f = np.interp(t_uniform, t_days, y_rot)
    xe_f = np.interp(t_uniform, t_days, x_eci)
    ye_f = np.interp(t_uniform, t_days, y_eci)
    xme_f = np.interp(t_uniform, t_days, xm_eci)
    yme_f = np.interp(t_uniform, t_days, ym_eci)
    alte_f = np.interp(t_uniform, t_days, alt_e)
    altm_f = np.interp(t_uniform, t_days, alt_m)
    ve_f = np.interp(t_uniform, t_days, v_eci)
    ph_f = np.round(np.interp(t_uniform, t_days, phase)).astype(int)

    phase_names = {
        0: ("OUTBOUND FIGURE-8 (Earth → Moon)", C_OUTBOUND),
        1: ("INBOUND FIGURE-8 (Moon → Earth)", C_INBOUND),
        2: ("PHASING LOOP 1 (Sub-Lunar)", C_PHASE1),
        3: ("PHASING LOOP 2 (Sub-Lunar)", C_PHASE2),
    }

    fig, (ax_rot, ax_eci) = plt.subplots(1, 2, figsize=(11.6, 5.8), facecolor=BG_DARK)
    fig.subplots_adjust(left=0.055, right=0.975, top=0.82, bottom=0.14, wspace=0.14)

    title_txt = fig.suptitle(
        "Aldrin–Genova Shamrock (C-2-R) Earth–Moon Cycler",
        color=TEXT_MAIN,
        fontsize=14,
        fontweight="bold",
        y=0.955,
    )
    hud_txt = fig.text(
        0.5,
        0.885,
        "",
        ha="center",
        color="#7dd3fc",
        fontsize=10.2,
        family="monospace",
    )
    footer_txt = fig.text(
        0.5,
        0.042,
        "",
        ha="center",
        color=TEXT_MUTED,
        fontsize=9.4,
        family="monospace",
    )

    # Static background - Rotating frame
    _style_ax(ax_rot, "Earth–Moon Rotating (Synodic) Frame", "x / L*", "y / L*")
    for code, col in [(0, C_OUTBOUND), (1, C_INBOUND), (2, C_PHASE1), (3, C_PHASE2)]:
        idx = np.where(phase == code)[0]
        if len(idx) > 1:
            end_i = min(idx[-1] + 2, len(x_rot))
            ax_rot.plot(x_rot[idx[0]:end_i], y_rot[idx[0]:end_i], color=col, alpha=0.28, linewidth=1.6)
    for name, (lx, ly) in lp.items():
        ax_rot.scatter([lx], [ly], color="#475569", marker="+", s=26, zorder=4)
        ax_rot.text(lx, ly + (0.06 if ly >= 0 else -0.095), name, color="#64748b", fontsize=7.2, ha="center")
    ax_rot.add_patch(Circle((-MU, 0.0), 0.044, color=C_EARTH, zorder=6))
    ax_rot.add_patch(Circle((1.0 - MU, 0.0), 0.025, color=C_MOON, zorder=6))
    ax_rot.text(-MU, -0.115, "Earth", color="#93c5fd", fontsize=8.2, ha="center", fontweight="bold")
    ax_rot.text(1.0 - MU, -0.105, "Moon", color="#e2e8f0", fontsize=8.2, ha="center", fontweight="bold")
    ax_rot.set_aspect("equal")
    ax_rot.set_xlim(-1.22, 1.25)
    ax_rot.set_ylim(-1.12, 1.12)

    # Static background - ECI frame
    _style_ax(ax_eci, "Earth-Centered Inertial (ECI) Frame", "X_ECI / L*", "Y_ECI / L*")
    th = np.linspace(0, 2 * np.pi, 240)
    ax_eci.plot(np.cos(th), np.sin(th), color="#475569", linestyle=":", linewidth=1.2)
    ax_eci.plot(x_eci, y_eci, color="#64748b", alpha=0.25, linewidth=1.4)
    ax_eci.add_patch(Circle((0.0, 0.0), 0.045, color=C_EARTH, zorder=6))
    ax_eci.text(0.0, -0.115, "Earth", color="#93c5fd", fontsize=8.2, ha="center", fontweight="bold")
    ax_eci.set_aspect("equal")
    ax_eci.set_xlim(-1.22, 1.22)
    ax_eci.set_ylim(-1.18, 1.18)

    # Animated artists
    trail_rot, = ax_rot.plot([], [], color=C_OUTBOUND, linewidth=2.6, zorder=7)
    sc_rot, = ax_rot.plot([], [], marker="o", markersize=8.5, color="white", markeredgecolor=C_OUTBOUND, markeredgewidth=2.0, zorder=9)

    trail_eci, = ax_eci.plot([], [], color=C_OUTBOUND, linewidth=2.6, zorder=7)
    sc_eci, = ax_eci.plot([], [], marker="o", markersize=8.5, color="white", markeredgecolor=C_OUTBOUND, markeredgewidth=2.0, zorder=9)
    moon_eci, = ax_eci.plot([], [], marker="o", markersize=11.0, color=C_MOON, markeredgecolor="#94a3b8", zorder=8)
    moon_lbl = ax_eci.text(0, 0, "Moon", color="#e2e8f0", fontsize=8.0, ha="center", fontweight="bold")

    trail_len = 28

    def update(frame: int):
        idx_slice = [(frame - k) % n_frames for k in range(trail_len, -1, -1)]
        col = phase_names[ph_f[frame]][1]
        p_label = phase_names[ph_f[frame]][0]

        # Break trail across cycle seam in ECI to avoid jump line
        t_cur = t_uniform[frame]
        if frame < trail_len:
            eci_idx = list(range(0, frame + 1))
        else:
            eci_idx = list(range(frame - trail_len, frame + 1))

        trail_rot.set_data(xr_f[idx_slice], yr_f[idx_slice])
        trail_rot.set_color(col)
        sc_rot.set_data([xr_f[frame]], [yr_f[frame]])
        sc_rot.set_markeredgecolor(col)

        trail_eci.set_data(xe_f[eci_idx], ye_f[eci_idx])
        trail_eci.set_color(col)
        sc_eci.set_data([xe_f[frame]], [ye_f[frame]])
        sc_eci.set_markeredgecolor(col)

        moon_eci.set_data([xme_f[frame]], [yme_f[frame]])
        moon_lbl.set_position((xme_f[frame] * 1.13, yme_f[frame] * 1.13))

        burn_tag = ""
        if abs(t_cur - sm["burn1_time_days"]) < 0.45:
            burn_tag = "  [★ BURN 1: -18.7 m/s]"
        elif t_cur < 0.35 or abs(t_cur - sm["period_days"]) < 0.45:
            burn_tag = "  [★ BURN 2: +18.7 m/s]"
        elif abs(t_cur - sm["perilune_time_days"]) < 0.35:
            burn_tag = "  [★ LUNAR PERILUNE: 3,000 km]"

        hud_txt.set_text(f"Day {t_cur:5.2f} / {sm['period_days']:.2f}   |   {p_label}{burn_tag}")
        hud_txt.set_color("#fca5a5" if "BURN" in burn_tag else col)
        footer_txt.set_text(
            f"Earth Alt: {alte_f[frame]:9,.0f} km   |   Moon Alt: {altm_f[frame]:9,.0f} km   |   "
            f"ECI Speed: {ve_f[frame]:5.2f} km/s"
        )
        return trail_rot, sc_rot, trail_eci, sc_eci, moon_eci, moon_lbl, hud_txt, footer_txt

    anim = FuncAnimation(fig, update, frames=n_frames, interval=1000 // fps, blit=False)
    out_gif.parent.mkdir(parents=True, exist_ok=True)
    anim.save(out_gif, writer=PillowWriter(fps=fps), dpi=110)
    plt.close(fig)
    print(f"Saved {out_gif}")


def main():
    root = Path(__file__).resolve().parent.parent
    data_path = root / "data" / "cycler_solution.json"
    if not data_path.exists():
        from cycler_solver import export_solution
        export_solution(data_path)
    data = json.loads(data_path.read_text(encoding="utf-8"))
    generate_analysis_figure(data, root / "assets" / "cycler_analysis.png")
    generate_animation_gif(data, root / "assets" / "aldrin_cycler_animation.gif")


if __name__ == "__main__":
    main()
