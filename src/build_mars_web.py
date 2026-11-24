"""
Build a self-contained interactive HTML5 Earth-Mars Aldrin Cycler simulator (web/mars.html)
embedding the exact numerical solution from data/mars_cycler_solution.json.
"""

from __future__ import annotations

import json
from pathlib import Path


MARS_HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>Earth–Mars Aldrin Cycler Simulation</title>
<style>
  * { box-sizing: border-box; }
  html, body {
    margin: 0;
    padding: 0;
    font-family: var(--font-sans, system-ui, -apple-system, Segoe UI, Roboto, sans-serif);
    color: var(--foreground, #e6edf7);
  }
  .wrap {
    width: 100%;
    display: flex;
    flex-direction: column;
    gap: 12px;
  }
  .toolbar {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    justify-content: space-between;
    gap: 8px;
    padding: 12px 16px;
    border-radius: var(--radius, 10px);
    background: var(--card, #11182c);
    border: 1px solid var(--border, #23304d);
  }
  .btn-group {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 6px;
  }
  button {
    font-family: inherit;
    font-size: 12px;
    font-weight: 600;
    padding: 6px 11px;
    border-radius: calc(var(--radius, 10px) - 3px);
    border: 1px solid var(--border, #2a3a5c);
    background: var(--secondary, #1a243d);
    color: var(--secondary-foreground, var(--foreground, #e6edf7));
    cursor: pointer;
    transition: opacity 0.15s ease;
  }
  button:hover { opacity: 0.9; }
  button.active {
    background: var(--primary, #38bdf8);
    color: var(--primary-foreground, #09111f);
    border-color: var(--primary, #38bdf8);
  }
  .scrubber-row {
    display: flex;
    align-items: center;
    gap: 10px;
    width: 100%;
    padding-top: 6px;
  }
  .scrubber-row input[type="range"] {
    flex: 1;
    accent-color: var(--primary, #38bdf8);
    cursor: pointer;
  }
  .mono {
    font-family: var(--font-mono, ui-monospace, SFMono-Regular, Menlo, Consolas, monospace);
    font-variant-numeric: tabular-nums;
  }
  .telemetry {
    display: grid;
    grid-template-columns: repeat(4, minmax(0, 1fr));
    gap: 8px;
  }
  .stat {
    padding: 10px 12px;
    border-radius: var(--radius, 10px);
    background: var(--card, #11182c);
    border: 1px solid var(--border, #23304d);
  }
  .stat-label {
    font-size: 10.5px;
    text-transform: uppercase;
    letter-spacing: 0.04em;
    color: var(--muted-foreground, #94a3b8);
    margin-bottom: 3px;
  }
  .stat-val {
    font-size: 14px;
    font-weight: 700;
    color: var(--foreground, #e6edf7);
  }
  .phase-pill {
    display: inline-flex;
    align-items: center;
    gap: 6px;
    font-size: 11.5px;
    font-weight: 600;
    padding: 4px 10px;
    border-radius: 999px;
    background: var(--secondary, #1a243d);
    border: 1px solid var(--border, #23304d);
  }
  .dot {
    width: 8px;
    height: 8px;
    border-radius: 50%;
    display: inline-block;
  }
  .views {
    display: grid;
    grid-template-columns: repeat(2, minmax(0, 1fr));
    gap: 12px;
  }
  .panel {
    padding: 14px 16px 12px 16px;
    border-radius: var(--radius, 10px);
    background: var(--card, #11182c);
    border: 1px solid var(--border, #23304d);
    display: flex;
    flex-direction: column;
    gap: 6px;
  }
  .panel-head {
    display: flex;
    align-items: baseline;
    justify-content: space-between;
    gap: 8px;
  }
  .panel-title {
    font-size: 13px;
    font-weight: 700;
  }
  .panel-sub {
    font-size: 11px;
    color: var(--muted-foreground, #94a3b8);
  }
  canvas.orbit-canvas {
    width: 100%;
    height: 300px;
    display: block;
    border-radius: calc(var(--radius, 10px) - 4px);
    background: var(--background, #0b1020);
  }
  canvas.chart-canvas {
    width: 100%;
    height: 205px;
    display: block;
  }
  .legend {
    display: flex;
    flex-wrap: wrap;
    gap: 9px;
    font-size: 10.8px;
    color: var(--muted-foreground, #94a3b8);
    padding-top: 2px;
  }
  .legend-item {
    display: inline-flex;
    align-items: center;
    gap: 5px;
  }
  .swatch {
    width: 10px;
    height: 3px;
    border-radius: 2px;
    display: inline-block;
  }
  @media (max-width: 600px) {
    .views { grid-template-columns: 1fr; }
    .telemetry { grid-template-columns: repeat(2, minmax(0, 1fr)); }
    canvas.orbit-canvas { height: 260px; }
  }
</style>
</head>
<body>
<div class="wrap">
  <div class="toolbar">
    <div class="btn-group">
      <button id="btnPlay" class="active">Pause</button>
      <button id="btnTwin" class="active">Twin Cyclers (Up + Down)</button>
      <button class="spd" data-speed="0.5">0.5×</button>
      <button class="spd active" data-speed="1">1×</button>
      <button class="spd" data-speed="2.5">2.5×</button>
    </div>
    <span id="phaseBadge" class="phase-pill"><span class="dot" id="phaseDot"></span><span id="phaseText">Fast Earth → Mars (145 d)</span></span>
    <div class="scrubber-row">
      <span class="mono" style="font-size:11.5px; color: var(--muted-foreground, #94a3b8);" id="timeReadout">Yr 0.00 · Cyc 1 (Day 0)</span>
      <input id="timeSlider" type="range" min="0" max="1000" value="0" />
    </div>
  </div>

  <div class="telemetry">
    <div class="stat">
      <div class="stat-label">LEO Taxi Escape Δv</div>
      <div class="stat-val mono" id="valLeoDv">4.92 km/s</div>
    </div>
    <div class="stat">
      <div class="stat-label">Mars Entry v_EI</div>
      <div class="stat-val mono" id="valMarsEi">10.87 km/s</div>
    </div>
    <div class="stat">
      <div class="stat-label">Helio Speed</div>
      <div class="stat-val mono" id="valSpeed">34.85 km/s</div>
    </div>
    <div class="stat">
      <div class="stat-label">Dist Earth / Mars</div>
      <div class="stat-val mono" id="valDist">0.00 / 0.89 AU</div>
    </div>
  </div>

  <div class="views">
    <div class="panel">
      <div class="panel-head">
        <span class="panel-title">Heliocentric Frame (3 Cycles)</span>
        <span class="panel-sub">+48.7° apsidal precession / Earth GA</span>
      </div>
      <canvas id="helioCanvas" class="orbit-canvas"></canvas>
      <div class="legend">
        <span class="legend-item"><span class="swatch" style="background:var(--chart-1,#38bdf8)"></span>Outbound Sprint (145 d)</span>
        <span class="legend-item"><span class="swatch" style="background:var(--warning,#f59e0b)"></span>Return Loop (635 d)</span>
        <span class="legend-item"><span class="swatch" style="background:var(--success,#34d399)"></span>Inbound Cycler</span>
      </div>
    </div>

    <div class="panel">
      <div class="panel-head">
        <span class="panel-title">Sun–Earth Rotating Frame</span>
        <span class="panel-sub">Earth pinned at (1, 0) AU</span>
      </div>
      <canvas id="rotCanvas" class="orbit-canvas"></canvas>
      <div class="legend">
        <span class="legend-item"><span class="swatch" style="background:#3b82f6"></span>Earth (1.00 AU)</span>
        <span class="legend-item"><span class="swatch" style="background:#ef4444"></span>Mars (1.52 AU)</span>
        <span class="legend-item"><span class="swatch" style="background:#ffffff"></span>Cycler + Crew Taxi</span>
      </div>
    </div>
  </div>

  <div class="views">
    <div class="panel">
      <div class="panel-head">
        <span class="panel-title">Taxi Δv &amp; Mars Entry vs. Hohmann</span>
        <span class="panel-sub">km/s · Aerocapture mandatory</span>
      </div>
      <canvas id="dvCanvas" class="chart-canvas"></canvas>
      <div class="legend">
        <span class="legend-item"><span class="swatch" style="background:var(--chart-1,#38bdf8)"></span>Aldrin Cycler (145 d)</span>
        <span class="legend-item"><span class="swatch" style="background:var(--warning,#f59e0b)"></span>Hohmann (259 d)</span>
      </div>
    </div>

    <div class="panel">
      <div class="panel-head">
        <span class="panel-title">15-Yr Eccentric Mars Correction Δv</span>
        <span class="panel-sub" id="corrSub">e_M = 0.0934 · 7 synodic cycles</span>
      </div>
      <canvas id="corrCanvas" class="chart-canvas"></canvas>
      <div class="legend">
        <span class="legend-item"><span class="swatch" style="background:var(--chart-1,#38bdf8)"></span>Earth Flyby/Turn</span>
        <span class="legend-item"><span class="swatch" style="background:var(--destructive,#ef4444)"></span>Mars Trim</span>
        <span class="legend-item"><span class="swatch" style="background:var(--chart-4,#a78bfa)"></span>Inclination + TCM</span>
      </div>
    </div>
  </div>
</div>

<script>
const DATA = __MARS_DATA_JSON__;

(function () {
  const cs = getComputedStyle(document.documentElement);
  const cssVar = (name, fb) => (cs.getPropertyValue(name) || "").trim() || fb;

  const COL = {
    bg: cssVar("--background", "#0b1020"),
    fg: cssVar("--foreground", "#e6edf7"),
    muted: cssVar("--muted-foreground", "#94a3b8"),
    fast: cssVar("--chart-1", "#38bdf8"),
    ret: cssVar("--warning", "#f59e0b"),
    inb: cssVar("--success", "#34d399"),
    purple: cssVar("--chart-4", "#a78bfa"),
    mars: "#ef4444",
    earth: "#3b82f6",
    sun: "#fbbf24",
  };

  const TRJ = DATA.trajectory;
  const IDEAL = DATA.ideal;
  const CMP = DATA.campaign_15yr;
  const T_SYN = IDEAL.T_synodic_days;
  const TOF_FAST = IDEAL.tof_fast_days;
  const T_MAX = TRJ.t_days[TRJ.t_days.length - 1];

  let playing = true;
  let showTwin = true;
  let speed = 1.0;
  let frac = 0.0;
  let lastTs = null;

  const helioCanvas = document.getElementById("helioCanvas");
  const rotCanvas = document.getElementById("rotCanvas");
  const dvCanvas = document.getElementById("dvCanvas");
  const corrCanvas = document.getElementById("corrCanvas");

  const slider = document.getElementById("timeSlider");
  const timeReadout = document.getElementById("timeReadout");
  const valLeoDv = document.getElementById("valLeoDv");
  const valMarsEi = document.getElementById("valMarsEi");
  const valSpeed = document.getElementById("valSpeed");
  const valDist = document.getElementById("valDist");
  const phaseDot = document.getElementById("phaseDot");
  const phaseText = document.getElementById("phaseText");
  const btnPlay = document.getElementById("btnPlay");
  const btnTwin = document.getElementById("btnTwin");

  valLeoDv.textContent = `${IDEAL.dv_leo_taxi_ideal_km_s.toFixed(2)} km/s`;
  valMarsEi.textContent = `${IDEAL.v_ei_mars_cycler_km_s.toFixed(2)} km/s`;
  document.getElementById("corrSub").textContent =
    `Avg ${Math.round(CMP.avg_dv_corr_per_cycle_m_s)} m/s/cyc (e = 0.0934)`;

  function resizeCanvas(canvas) {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const rect = canvas.getBoundingClientRect();
    const w = Math.max(180, Math.round(rect.width));
    const h = Math.max(160, Math.round(rect.height));
    if (canvas.width !== w * dpr || canvas.height !== h * dpr) {
      canvas.width = w * dpr;
      canvas.height = h * dpr;
    }
    const ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx, w, h };
  }

  function sampleAt(tTarget) {
    const ts = TRJ.t_days;
    const n = ts.length;
    if (tTarget <= ts[0]) return { i: 0, a: 0 };
    if (tTarget >= ts[n - 1]) return { i: n - 2, a: 1 };
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (ts[mid] <= tTarget) lo = mid; else hi = mid;
    }
    return { i: lo, a: (tTarget - ts[lo]) / ((ts[hi] - ts[lo]) || 1e-9) };
  }

  const lerp = (arr, i, a) => arr[i] + (arr[i + 1] - arr[i]) * a;

  function drawSpaceBg(ctx, w, h, lim) {
    ctx.fillStyle = COL.bg;
    ctx.fillRect(0, 0, w, h);
    const scale = Math.min(w, h) / (2.0 * lim);
    const cx = w * 0.5, cy = h * 0.5;
    ctx.strokeStyle = "rgba(148,163,184,0.12)";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, cy); ctx.lineTo(w, cy);
    ctx.moveTo(cx, 0); ctx.lineTo(cx, h);
    ctx.stroke();

    // Earth & Mars orbit rings
    ctx.setLineDash([3, 4]);
    ctx.strokeStyle = "rgba(59,130,246,0.35)";
    ctx.beginPath(); ctx.arc(cx, cy, 1.0 * scale, 0, Math.PI * 2); ctx.stroke();
    ctx.strokeStyle = "rgba(239,68,68,0.35)";
    ctx.beginPath(); ctx.arc(cx, cy, 1.5237 * scale, 0, Math.PI * 2); ctx.stroke();
    ctx.setLineDash([]);

    // Sun
    ctx.fillStyle = COL.sun;
    ctx.beginPath(); ctx.arc(cx, cy, 6.5, 0, Math.PI * 2); ctx.fill();
    return { cx, cy, scale };
  }

  function drawHelio(s, tNow) {
    const { ctx, w, h } = resizeCanvas(helioCanvas);
    const { cx, cy, scale } = drawSpaceBg(ctx, w, h, 2.42);
    const toXY = (x, y) => [cx + x * scale, cy - y * scale];

    // Full 3-cycle precession backdrop
    ctx.strokeStyle = "rgba(148,163,184,0.22)";
    ctx.lineWidth = 1.25;
    ctx.beginPath();
    for (let k = 0; k < TRJ.out_x.length; k++) {
      const [px, py] = toXY(TRJ.out_x[k], TRJ.out_y[k]);
      if (k === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    }
    ctx.stroke();

    // Active fading trail for Outbound Cycler
    const dtCyc = tNow % T_SYN;
    const isFast = dtCyc <= TOF_FAST;
    const activeCol = isFast ? COL.fast : COL.ret;
    const trailDays = 175;
    const steps = 42;
    ctx.lineWidth = 2.6;
    ctx.lineCap = "round";
    for (let k = 1; k <= steps; k++) {
      const t0 = tNow - trailDays * (k / steps);
      const t1 = tNow - trailDays * ((k - 1) / steps);
      if (t0 < 0) continue;
      const p0 = sampleAt(t0), p1 = sampleAt(t1);
      const [x0, y0] = toXY(lerp(TRJ.out_x, p0.i, p0.a), lerp(TRJ.out_y, p0.i, p0.a));
      const [x1, y1] = toXY(lerp(TRJ.out_x, p1.i, p1.a), lerp(TRJ.out_y, p1.i, p1.a));
      ctx.strokeStyle = activeCol;
      ctx.globalAlpha = 0.9 * (1 - k / steps);
      ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
    }
    ctx.globalAlpha = 1.0;

    // Inbound cycler trail & marker
    if (showTwin) {
      for (let k = 1; k <= 24; k++) {
        const t0 = tNow - 110 * (k / 24);
        const t1 = tNow - 110 * ((k - 1) / 24);
        if (t0 < 0) continue;
        const p0 = sampleAt(t0), p1 = sampleAt(t1);
        const [x0, y0] = toXY(lerp(TRJ.in_x, p0.i, p0.a), lerp(TRJ.in_y, p0.i, p0.a));
        const [x1, y1] = toXY(lerp(TRJ.in_x, p1.i, p1.a), lerp(TRJ.in_y, p1.i, p1.a));
        ctx.strokeStyle = COL.inb;
        ctx.globalAlpha = 0.6 * (1 - k / 24);
        ctx.lineWidth = 1.8;
        ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
      }
      ctx.globalAlpha = 1.0;
      const [ix, iy] = toXY(lerp(TRJ.in_x, s.i, s.a), lerp(TRJ.in_y, s.i, s.a));
      ctx.fillStyle = COL.inb;
      ctx.strokeStyle = "#ffffff";
      ctx.lineWidth = 1.2;
      ctx.fillRect(ix - 3.8, iy - 3.8, 7.6, 7.6);
      ctx.strokeRect(ix - 3.8, iy - 3.8, 7.6, 7.6);
    }

    // Earth & Mars
    const [ex, ey] = toXY(lerp(TRJ.earth_x, s.i, s.a), lerp(TRJ.earth_y, s.i, s.a));
    const [mx, my] = toXY(lerp(TRJ.mars_x, s.i, s.a), lerp(TRJ.mars_y, s.i, s.a));
    ctx.fillStyle = COL.earth;
    ctx.beginPath(); ctx.arc(ex, ey, 5.8, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = COL.mars;
    ctx.beginPath(); ctx.arc(mx, my, 5.0, 0, Math.PI * 2); ctx.fill();

    ctx.font = "bold 10px system-ui, sans-serif";
    ctx.textAlign = "center";
    ctx.fillStyle = "#93c5fd"; ctx.fillText("Earth", ex, ey - 9);
    ctx.fillStyle = "#fca5a5"; ctx.fillText("Mars", mx, my - 9);

    // Outbound Cycler (+ Taxi glow during 145 d crew sprint)
    const [ox, oy] = toXY(lerp(TRJ.out_x, s.i, s.a), lerp(TRJ.out_y, s.i, s.a));
    if (isFast) {
      ctx.strokeStyle = COL.fast;
      ctx.globalAlpha = 0.38;
      ctx.lineWidth = 4;
      ctx.beginPath(); ctx.arc(ox, oy, 9.5, 0, Math.PI * 2); ctx.stroke();
      ctx.globalAlpha = 1.0;
    }
    ctx.fillStyle = "#ffffff";
    ctx.strokeStyle = activeCol;
    ctx.lineWidth = 2.2;
    ctx.beginPath(); ctx.arc(ox, oy, 4.8, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
  }

  function drawRot(s, tNow) {
    const { ctx, w, h } = resizeCanvas(rotCanvas);
    const { cx, cy, scale } = drawSpaceBg(ctx, w, h, 2.42);
    const toXY = (x, y) => [cx + x * scale, cy - y * scale];

    // Full synodic loop backdrop
    ctx.strokeStyle = "rgba(148,163,184,0.22)";
    ctx.lineWidth = 1.25;
    ctx.beginPath();
    for (let k = 0; k < TRJ.out_rot_x.length; k++) {
      const [px, py] = toXY(TRJ.out_rot_x[k], TRJ.out_rot_y[k]);
      if (k === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    }
    ctx.stroke();

    const dtCyc = tNow % T_SYN;
    const isFast = dtCyc <= TOF_FAST;
    const activeCol = isFast ? COL.fast : COL.ret;
    const trailDays = 175;
    const steps = 42;
    ctx.lineWidth = 2.6;
    ctx.lineCap = "round";
    for (let k = 1; k <= steps; k++) {
      const t0 = tNow - trailDays * (k / steps);
      const t1 = tNow - trailDays * ((k - 1) / steps);
      if (t0 < 0) continue;
      const p0 = sampleAt(t0), p1 = sampleAt(t1);
      const [x0, y0] = toXY(lerp(TRJ.out_rot_x, p0.i, p0.a), lerp(TRJ.out_rot_y, p0.i, p0.a));
      const [x1, y1] = toXY(lerp(TRJ.out_rot_x, p1.i, p1.a), lerp(TRJ.out_rot_y, p1.i, p1.a));
      ctx.strokeStyle = activeCol;
      ctx.globalAlpha = 0.9 * (1 - k / steps);
      ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
    }
    ctx.globalAlpha = 1.0;

    if (showTwin) {
      const [ix, iy] = toXY(lerp(TRJ.in_rot_x, s.i, s.a), lerp(TRJ.in_rot_y, s.i, s.a));
      ctx.fillStyle = COL.inb;
      ctx.strokeStyle = "#ffffff";
      ctx.lineWidth = 1.2;
      ctx.fillRect(ix - 3.8, iy - 3.8, 7.6, 7.6);
      ctx.strokeRect(ix - 3.8, iy - 3.8, 7.6, 7.6);
    }

    // Fixed Earth at (1, 0)
    const [ex, ey] = toXY(1.0, 0.0);
    ctx.fillStyle = COL.earth;
    ctx.beginPath(); ctx.arc(ex, ey, 6.2, 0, Math.PI * 2); ctx.fill();
    ctx.font = "bold 10px system-ui, sans-serif";
    ctx.textAlign = "center";
    ctx.fillStyle = "#93c5fd"; ctx.fillText("Earth", ex, ey + 16);

    // Mars rotated into Sun-Earth frame
    const exH = lerp(TRJ.earth_x, s.i, s.a), eyH = lerp(TRJ.earth_y, s.i, s.a);
    const mxH = lerp(TRJ.mars_x, s.i, s.a), myH = lerp(TRJ.mars_y, s.i, s.a);
    const thE = Math.atan2(eyH, exH);
    const mrx = mxH * Math.cos(-thE) - myH * Math.sin(-thE);
    const mry = mxH * Math.sin(-thE) + myH * Math.cos(-thE);
    const [mx, my] = toXY(mrx, mry);
    ctx.fillStyle = COL.mars;
    ctx.beginPath(); ctx.arc(mx, my, 5.0, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = "#fca5a5"; ctx.fillText("Mars", mx, my - 9);

    const [ox, oy] = toXY(lerp(TRJ.out_rot_x, s.i, s.a), lerp(TRJ.out_rot_y, s.i, s.a));
    ctx.fillStyle = "#ffffff";
    ctx.strokeStyle = activeCol;
    ctx.lineWidth = 2.2;
    ctx.beginPath(); ctx.arc(ox, oy, 4.8, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
  }

  function drawDvChart() {
    const { ctx, w, h } = resizeCanvas(dvCanvas);
    const pad = { l: 38, r: 10, t: 16, b: 36 };
    const pw = w - pad.l - pad.r, ph = h - pad.t - pad.b;
    ctx.clearRect(0, 0, w, h);

    const items = [
      { lbl: "LEO Taxi", cyc: IDEAL.dv_leo_taxi_ideal_km_s, hoh: IDEAL.dv_leo_hohmann_km_s },
      { lbl: "Mars v_EI", cyc: IDEAL.v_ei_mars_cycler_km_s, hoh: IDEAL.v_ei_mars_hohmann_km_s },
      { lbl: "Prop. LMO", cyc: IDEAL.dv_mars_propulsive_lmo_km_s, hoh: IDEAL.dv_mars_propulsive_lmo_hoh_km_s },
      { lbl: "Aero Trim", cyc: IDEAL.dv_mars_aerocapture_circ_km_s, hoh: 0.14 },
    ];
    const maxY = 12.0;
    const yMap = v => pad.t + ph - (v / maxY) * ph;

    ctx.strokeStyle = "rgba(148,163,184,0.15)";
    ctx.fillStyle = COL.muted;
    ctx.font = "10px ui-monospace, monospace";
    for (let g = 0; g <= 4; g++) {
      const val = (maxY * g) / 4;
      const y = yMap(val);
      ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(w - pad.r, y); ctx.stroke();
      ctx.textAlign = "right";
      ctx.fillText(val.toFixed(0), pad.l - 5, y + 3.5);
    }

    const slot = pw / items.length;
    const bw = Math.min(26, slot * 0.34);
    items.forEach((it, idx) => {
      const cx = pad.l + slot * (idx + 0.5);
      // Cycler bar
      const y1 = yMap(it.cyc);
      ctx.fillStyle = COL.fast;
      ctx.fillRect(cx - bw - 1.5, y1, bw, pad.t + ph - y1);
      ctx.fillStyle = COL.fg;
      ctx.textAlign = "center";
      ctx.font = "bold 9.5px ui-monospace, monospace";
      ctx.fillText(it.cyc.toFixed(1), cx - bw * 0.5 - 1.5, y1 - 4);

      // Hohmann bar
      const y2 = yMap(it.hoh);
      ctx.fillStyle = COL.ret;
      ctx.fillRect(cx + 1.5, y2, bw, pad.t + ph - y2);
      ctx.fillStyle = COL.muted;
      ctx.font = "9.5px ui-monospace, monospace";
      ctx.fillText(it.hoh.toFixed(1), cx + bw * 0.5 + 1.5, y2 - 4);

      ctx.fillStyle = COL.fg;
      ctx.font = "10px system-ui, sans-serif";
      ctx.fillText(it.lbl, cx, h - 14);
    });
  }

  function drawCorrChart() {
    const { ctx, w, h } = resizeCanvas(corrCanvas);
    const pad = { l: 42, r: 10, t: 16, b: 36 };
    const pw = w - pad.l - pad.r, ph = h - pad.t - pad.b;
    ctx.clearRect(0, 0, w, h);

    const cycs = CMP.cycles;
    const maxY = Math.ceil(Math.max(...cycs.map(c => c.dv_total_corr_m_s)) / 200) * 200;
    const yMap = v => pad.t + ph - (v / maxY) * ph;

    ctx.strokeStyle = "rgba(148,163,184,0.15)";
    ctx.fillStyle = COL.muted;
    ctx.font = "10px ui-monospace, monospace";
    for (let g = 0; g <= 4; g++) {
      const val = (maxY * g) / 4;
      const y = yMap(val);
      ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(w - pad.r, y); ctx.stroke();
      ctx.textAlign = "right";
      ctx.fillText(Math.round(val), pad.l - 5, y + 3.5);
    }

    const slot = pw / cycs.length;
    const bw = Math.min(28, slot * 0.58);
    cycs.forEach((c, idx) => {
      const cx = pad.l + slot * (idx + 0.5);
      const x0 = cx - bw * 0.5;
      let acc = 0;
      const parts = [
        [c.dv_earth_flyby_m_s, COL.fast],
        [c.dv_mars_corr_m_s, COL.mars],
        [c.dv_inclination_m_s + c.dv_nav_tcm_m_s, COL.purple],
      ];
      parts.forEach(([val, col]) => {
        const yTop = yMap(acc + val);
        const yBot = yMap(acc);
        ctx.fillStyle = col;
        ctx.fillRect(x0, yTop, bw, Math.max(0, yBot - yTop));
        acc += val;
      });
      ctx.fillStyle = COL.fg;
      ctx.font = "bold 9px ui-monospace, monospace";
      ctx.textAlign = "center";
      ctx.fillText(Math.round(c.dv_total_corr_m_s), cx, yMap(acc) - 4);
      ctx.fillStyle = COL.muted;
      ctx.font = "10px system-ui, sans-serif";
      ctx.fillText(`C${c.cycle}`, cx, h - 18);
      ctx.font = "9px ui-monospace, monospace";
      ctx.fillText(`${c.mars_radius_au.toFixed(2)}`, cx, h - 6);
    });
  }

  function renderFrame() {
    const tNow = frac * T_MAX;
    const s = sampleAt(tNow);
    const cyc = Math.floor(tNow / T_SYN) + 1;
    const dtCyc = tNow % T_SYN;
    const isFast = dtCyc <= TOF_FAST;

    timeReadout.textContent = `Yr ${(tNow / 365.25).toFixed(2)} · Cyc ${cyc} (Day ${Math.round(dtCyc)})`;
    valSpeed.textContent = `${lerp(TRJ.out_speed_km_s, s.i, s.a).toFixed(2)} km/s`;
    valDist.textContent = `${lerp(TRJ.out_dist_earth_au, s.i, s.a).toFixed(2)} / ${lerp(TRJ.out_dist_mars_au, s.i, s.a).toFixed(2)} AU`;
    slider.value = Math.round(frac * 1000);

    let pLbl = isFast ? "Fast Earth → Mars (145 d, Crewed)" : "Deep-Space Return (635 d, 2.22 AU)";
    let pCol = isFast ? COL.fast : COL.ret;
    if (dtCyc < 20 || dtCyc > T_SYN - 20) {
      pLbl = "★ Earth Gravity Assist & Taxi Docking";
      pCol = COL.earth;
    } else if (Math.abs(dtCyc - TOF_FAST) < 16) {
      pLbl = "★ Mars Encounter (10.9 km/s Aerocapture!)";
      pCol = COL.mars;
    }
    phaseDot.style.background = pCol;
    phaseText.textContent = pLbl;

    drawHelio(s, tNow);
    drawRot(s, tNow);
  }

  function tick(ts) {
    if (lastTs !== null && playing) {
      const dt = Math.min((ts - lastTs) / 1000, 0.1);
      frac = (frac + (dt * speed) / 16.0) % 1.0;
      renderFrame();
    }
    lastTs = ts;
    requestAnimationFrame(tick);
  }

  btnPlay.addEventListener("click", () => {
    playing = !playing;
    btnPlay.textContent = playing ? "Pause" : "Play";
    btnPlay.classList.toggle("active", playing);
  });
  btnTwin.addEventListener("click", () => {
    showTwin = !showTwin;
    btnTwin.classList.toggle("active", showTwin);
    renderFrame();
  });
  document.querySelectorAll("button.spd").forEach(btn => {
    btn.addEventListener("click", () => {
      speed = parseFloat(btn.dataset.speed);
      document.querySelectorAll("button.spd").forEach(b => b.classList.remove("active"));
      btn.classList.add("active");
    });
  });
  slider.addEventListener("input", () => {
    frac = parseInt(slider.value, 10) / 1000;
    renderFrame();
  });
  window.addEventListener("resize", () => {
    renderFrame();
    drawDvChart();
    drawCorrChart();
  });

  renderFrame();
  drawDvChart();
  drawCorrChart();
  requestAnimationFrame(tick);
})();
</script>
</body>
</html>
"""


def build():
    root = Path(__file__).resolve().parent.parent
    data_path = root / "data" / "mars_cycler_solution.json"
    compact = json.dumps(json.loads(data_path.read_text(encoding="utf-8")), separators=(",", ":"))
    html = MARS_HTML_TEMPLATE.replace("__MARS_DATA_JSON__", compact)
    out_html = root / "web" / "mars.html"
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html, encoding="utf-8")
    print(f"Built {out_html} ({len(html):,} bytes)")


if __name__ == "__main__":
    build()
