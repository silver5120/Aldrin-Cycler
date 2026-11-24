"""
Build a self-contained interactive HTML5 Earth-Moon Aldrin Cycler simulator (web/index.html)
embedding the exact SciPy CR3BP numerical solution from data/cycler_solution.json.
"""

from __future__ import annotations

import json
from pathlib import Path


HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8" />
<meta name="viewport" content="width=device-width, initial-scale=1" />
<title>Earth–Moon Aldrin–Genova Shamrock Cycler</title>
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
    transition: opacity 0.15s ease, border-color 0.15s ease;
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
    height: 310px;
    display: block;
    border-radius: calc(var(--radius, 10px) - 4px);
    background: var(--background, #0b1020);
  }
  canvas.chart-canvas {
    width: 100%;
    height: 220px;
    display: block;
  }
  .legend {
    display: flex;
    flex-wrap: wrap;
    gap: 10px;
    font-size: 11px;
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
    canvas.orbit-canvas { height: 270px; }
  }
</style>
</head>
<body>
<div class="wrap">
  <div class="toolbar">
    <div class="btn-group">
      <button id="btnPlay" class="active">Pause</button>
      <button id="btnModeShamrock" class="active">Shamrock C-2-R (37.4 m/s)</button>
      <button id="btnModeArenstorf">Ballistic Arenstorf (0 m/s)</button>
    </div>
    <div class="btn-group">
      <button class="spd" data-speed="0.5">0.5×</button>
      <button class="spd active" data-speed="1">1×</button>
      <button class="spd" data-speed="2.5">2.5×</button>
      <span id="phaseBadge" class="phase-pill"><span class="dot" id="phaseDot"></span><span id="phaseText">Outbound Fig-8</span></span>
    </div>
    <div class="scrubber-row">
      <span class="mono" style="font-size:11.5px; color: var(--muted-foreground, #94a3b8);" id="dayReadout">Day 0.00 / 25.96</span>
      <input id="timeSlider" type="range" min="0" max="1000" value="0" />
    </div>
  </div>

  <div class="telemetry">
    <div class="stat">
      <div class="stat-label">Earth Altitude</div>
      <div class="stat-val mono" id="valEarthAlt">3,000 km</div>
    </div>
    <div class="stat">
      <div class="stat-label">Moon Altitude</div>
      <div class="stat-val mono" id="valMoonAlt">379,100 km</div>
    </div>
    <div class="stat">
      <div class="stat-label">ECI Speed</div>
      <div class="stat-val mono" id="valSpeed">9.13 km/s</div>
    </div>
    <div class="stat">
      <div class="stat-label">Jacobi C_J</div>
      <div class="stat-val mono" id="valJacobi">2.1181</div>
    </div>
  </div>

  <div class="views">
    <div class="panel">
      <div class="panel-head">
        <span class="panel-title">Rotating (Synodic) Frame</span>
        <span class="panel-sub">Earth &amp; Moon fixed · 1 L* = 384,400 km</span>
      </div>
      <canvas id="rotCanvas" class="orbit-canvas"></canvas>
      <div class="legend" id="rotLegend"></div>
    </div>

    <div class="panel">
      <div class="panel-head">
        <span class="panel-title">Earth-Centered Inertial (ECI) Frame</span>
        <span class="panel-sub">Moon orbits Earth · 27.32 d sidereal</span>
      </div>
      <canvas id="eciCanvas" class="orbit-canvas"></canvas>
      <div class="legend">
        <span class="legend-item"><span class="swatch" style="background:var(--chart-1,#38bdf8)"></span>Active Cycle</span>
        <span class="legend-item"><span class="swatch" style="background:#64748b"></span>Lunar Orbit (384.4k km)</span>
        <span class="legend-item"><span class="swatch" style="background:var(--warning,#f59e0b)"></span>3-Cycle Precession</span>
      </div>
    </div>
  </div>

  <div class="panel">
    <div class="panel-head">
      <span class="panel-title">Tsiolkovsky Cumulative IMLEO Trade Study (1–25 Crewed Lunar Flights)</span>
      <span class="panel-sub" id="tradeSub">Break-even vs Lean 18 t Capsule: ~6.8 flights · vs 45 t Heavy Hab: ~1.0 flight</span>
    </div>
    <canvas id="tradeCanvas" class="chart-canvas"></canvas>
    <div class="legend">
      <span class="legend-item"><span class="swatch" style="background:var(--destructive,#fb7185)"></span>Direct w/ 35 t Shielded Hab (134 t/flight)</span>
      <span class="legend-item"><span class="swatch" style="background:var(--warning,#f59e0b)"></span>Lean Orion/Apollo 18 t Capsule (62 t/flight)</span>
      <span class="legend-item"><span class="swatch" style="background:var(--success,#34d399)"></span>Aldrin Cycler (87 t setup + 49 t/flight)</span>
    </div>
  </div>
</div>

<script>
const DATA = __CYCLER_DATA_JSON__;

(function () {
  const cs = getComputedStyle(document.documentElement);
  const cssVar = (name, fb) => (cs.getPropertyValue(name) || "").trim() || fb;

  const COLORS = {
    bg: cssVar("--background", "#0b1020"),
    card: cssVar("--card", "#11182c"),
    border: cssVar("--border", "#23304d"),
    fg: cssVar("--foreground", "#e6edf7"),
    muted: cssVar("--muted-foreground", "#94a3b8"),
    outbound: cssVar("--chart-1", "#38bdf8"),
    inbound: cssVar("--success", "#34d399"),
    phase1: cssVar("--warning", "#f59e0b"),
    phase2: cssVar("--chart-4", "#a78bfa"),
    arenstorf: cssVar("--destructive", "#fb7185"),
    earth: "#3b82f6",
    moon: "#cbd5e1",
    burn: "#ef4444",
  };

  const MU = DATA.constants.mu;
  let mode = "shamrock"; // "shamrock" | "arenstorf"
  let playing = true;
  let speed = 1.0;
  let frac = 0.0; // normalized time [0, 1)
  let lastTs = null;

  const rotCanvas = document.getElementById("rotCanvas");
  const eciCanvas = document.getElementById("eciCanvas");
  const tradeCanvas = document.getElementById("tradeCanvas");
  const slider = document.getElementById("timeSlider");
  const dayReadout = document.getElementById("dayReadout");
  const valEarthAlt = document.getElementById("valEarthAlt");
  const valMoonAlt = document.getElementById("valMoonAlt");
  const valSpeed = document.getElementById("valSpeed");
  const valJacobi = document.getElementById("valJacobi");
  const phaseDot = document.getElementById("phaseDot");
  const phaseText = document.getElementById("phaseText");
  const rotLegend = document.getElementById("rotLegend");
  const btnPlay = document.getElementById("btnPlay");
  const btnShamrock = document.getElementById("btnModeShamrock");
  const btnArenstorf = document.getElementById("btnModeArenstorf");

  function updateLegend() {
    if (mode === "shamrock") {
      rotLegend.innerHTML =
        `<span class="legend-item"><span class="swatch" style="background:${COLORS.outbound}"></span>Outbound (3.4 d)</span>` +
        `<span class="legend-item"><span class="swatch" style="background:${COLORS.inbound}"></span>Inbound (3.4 d)</span>` +
        `<span class="legend-item"><span class="swatch" style="background:${COLORS.phase1}"></span>Phasing 1 (9.6 d)</span>` +
        `<span class="legend-item"><span class="swatch" style="background:${COLORS.phase2}"></span>Phasing 2 (9.6 d)</span>` +
        `<span class="legend-item"><span class="swatch" style="background:${COLORS.burn}"></span>±18.7 m/s Burn</span>`;
    } else {
      rotLegend.innerHTML =
        `<span class="legend-item"><span class="swatch" style="background:${COLORS.arenstorf}"></span>Arenstorf 3-Lobe (Δv = 0 m/s, 74.1 d)</span>` +
        `<span class="legend-item" style="color:${COLORS.warning}">Min Earth Alt: 171,635 km · Unstable (ν=147)</span>`;
    }
  }

  function resizeCanvas(canvas) {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const rect = canvas.getBoundingClientRect();
    const w = Math.max(200, Math.round(rect.width));
    const h = Math.max(160, Math.round(rect.height));
    if (canvas.width !== w * dpr || canvas.height !== h * dpr) {
      canvas.width = w * dpr;
      canvas.height = h * dpr;
    }
    const ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx, w, h };
  }

  // Binary search linear interpolation by mission time t_days
  function sampleAtTime(orbit, tTarget) {
    const ts = orbit.t_days;
    const n = ts.length;
    if (tTarget <= ts[0]) return { i: 0, alpha: 0 };
    if (tTarget >= ts[n - 1]) return { i: n - 2, alpha: 1 };
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (ts[mid] <= tTarget) lo = mid;
      else hi = mid;
    }
    const denom = ts[hi] - ts[lo] || 1e-9;
    return { i: lo, alpha: (tTarget - ts[lo]) / denom };
  }

  function lerp(arr, i, a) {
    return arr[i] + (arr[i + 1] - arr[i]) * a;
  }

  function mapPt(x, y, w, h, lim) {
    const scale = Math.min(w, h) / (2.0 * lim);
    return [w * 0.5 + x * scale, h * 0.5 - y * scale, scale];
  }

  function drawGrid(ctx, w, h, lim) {
    ctx.fillStyle = COLORS.bg;
    ctx.fillRect(0, 0, w, h);
    const [cx, cy, scale] = mapPt(0, 0, w, h, lim);
    ctx.strokeStyle = "rgba(148,163,184,0.12)";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, cy); ctx.lineTo(w, cy);
    ctx.moveTo(cx, 0); ctx.lineTo(cx, h);
    ctx.stroke();

    ctx.setLineDash([3, 4]);
    ctx.strokeStyle = "rgba(148,163,184,0.14)";
    [0.5, 1.0].forEach(r => {
      ctx.beginPath();
      ctx.arc(cx, cy, r * scale, 0, Math.PI * 2);
      ctx.stroke();
    });
    ctx.setLineDash([]);
    return scale;
  }

  function drawRotating(orbit, idxInfo, tNow) {
    const { ctx, w, h } = resizeCanvas(rotCanvas);
    const lim = mode === "shamrock" ? 1.24 : 1.38;
    const scale = drawGrid(ctx, w, h, lim);
    const toXY = (x, y) => [w * 0.5 + x * scale, h * 0.5 - y * scale];

    // Full orbit path
    if (mode === "shamrock") {
      const segColors = [COLORS.outbound, COLORS.inbound, COLORS.phase1, COLORS.phase2];
      for (let c = 0; c < 4; c++) {
        ctx.strokeStyle = segColors[c];
        ctx.globalAlpha = 0.35;
        ctx.lineWidth = 1.8;
        ctx.beginPath();
        let started = false;
        for (let k = 0; k < orbit.x_rot.length; k++) {
          if (orbit.phase[k] === c || (started && k > 0 && orbit.phase[k - 1] === c)) {
            const [px, py] = toXY(orbit.x_rot[k], orbit.y_rot[k]);
            if (!started) { ctx.moveTo(px, py); started = true; }
            else ctx.lineTo(px, py);
            if (orbit.phase[k] !== c) break;
          }
        }
        ctx.stroke();
      }
      ctx.globalAlpha = 1.0;

      // Burn markers
      const sp = orbit.metrics.state_perigee_out;
      [[sp[0], sp[1]], [sp[0], -sp[1]]].forEach(([bx, by]) => {
        const [px, py] = toXY(bx, by);
        ctx.fillStyle = COLORS.burn;
        ctx.beginPath();
        ctx.arc(px, py, 3.8, 0, Math.PI * 2);
        ctx.fill();
      });
    } else {
      ctx.strokeStyle = COLORS.arenstorf;
      ctx.globalAlpha = 0.4;
      ctx.lineWidth = 1.8;
      ctx.beginPath();
      for (let k = 0; k < orbit.x_rot.length; k++) {
        const [px, py] = toXY(orbit.x_rot[k], orbit.y_rot[k]);
        if (k === 0) ctx.moveTo(px, py);
        else ctx.lineTo(px, py);
      }
      ctx.stroke();
      ctx.globalAlpha = 1.0;
    }

    // Lagrange points
    ctx.font = "10px ui-monospace, monospace";
    ctx.textAlign = "center";
    Object.entries(DATA.lagrange_points).forEach(([name, [lx, ly]]) => {
      const [px, py] = toXY(lx, ly);
      ctx.strokeStyle = "rgba(148,163,184,0.45)";
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(px - 3, py); ctx.lineTo(px + 3, py);
      ctx.moveTo(px, py - 3); ctx.lineTo(px, py + 3);
      ctx.stroke();
      ctx.fillStyle = COLORS.muted;
      ctx.fillText(name, px, py + (ly >= 0 ? -6 : 12));
    });

    // Earth & Moon
    const [ex, ey] = toXY(-MU, 0);
    const [mx, my] = toXY(1.0 - MU, 0);
    ctx.fillStyle = COLORS.earth;
    ctx.beginPath(); ctx.arc(ex, ey, 7.5, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = COLORS.moon;
    ctx.beginPath(); ctx.arc(mx, my, 4.8, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = "#93c5fd";
    ctx.font = "bold 10.5px system-ui, sans-serif";
    ctx.fillText("Earth", ex, ey + 18);
    ctx.fillStyle = "#e2e8f0";
    ctx.fillText("Moon", mx, my + 16);

    // Active trail & spacecraft
    const segColors = [COLORS.outbound, COLORS.inbound, COLORS.phase1, COLORS.phase2];
    const activeCol = mode === "shamrock" ? segColors[orbit.phase[idxInfo.i]] : COLORS.arenstorf;
    const scx = lerp(orbit.x_rot, idxInfo.i, idxInfo.alpha);
    const scy = lerp(orbit.y_rot, idxInfo.i, idxInfo.alpha);

    const trailWindowDays = mode === "shamrock" ? 3.2 : 8.0;
    const period = orbit.metrics.period_days;
    const trailSteps = 45;
    ctx.lineWidth = 2.6;
    ctx.lineCap = "round";
    for (let s = 1; s <= trailSteps; s++) {
      let t0 = tNow - trailWindowDays * (s / trailSteps);
      let t1 = tNow - trailWindowDays * ((s - 1) / trailSteps);
      if (t0 < 0) t0 += period;
      if (t1 < 0) t1 += period;
      const p0 = sampleAtTime(orbit, t0);
      const p1 = sampleAtTime(orbit, t1);
      const [x0, y0] = toXY(lerp(orbit.x_rot, p0.i, p0.alpha), lerp(orbit.y_rot, p0.i, p0.alpha));
      const [x1, y1] = toXY(lerp(orbit.x_rot, p1.i, p1.alpha), lerp(orbit.y_rot, p1.i, p1.alpha));
      ctx.strokeStyle = activeCol;
      ctx.globalAlpha = 0.9 * (1 - s / trailSteps);
      ctx.beginPath();
      ctx.moveTo(x0, y0);
      ctx.lineTo(x1, y1);
      ctx.stroke();
    }
    ctx.globalAlpha = 1.0;

    const [spx, spy] = toXY(scx, scy);
    ctx.fillStyle = "#ffffff";
    ctx.strokeStyle = activeCol;
    ctx.lineWidth = 2.4;
    ctx.beginPath(); ctx.arc(spx, spy, 5.2, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
  }

  function drawECI(orbit, idxInfo, tNow) {
    const { ctx, w, h } = resizeCanvas(eciCanvas);
    const lim = mode === "shamrock" ? 1.24 : 1.38;
    const scale = drawGrid(ctx, w, h, lim);
    const toXY = (x, y) => [w * 0.5 + x * scale, h * 0.5 - y * scale];

    // Lunar orbit circle
    const [cx, cy] = toXY(0, 0);
    ctx.strokeStyle = "rgba(148,163,184,0.35)";
    ctx.setLineDash([4, 4]);
    ctx.lineWidth = 1.2;
    ctx.beginPath(); ctx.arc(cx, cy, scale, 0, Math.PI * 2); ctx.stroke();
    ctx.setLineDash([]);

    // 3-cycle precession backdrop for Shamrock
    if (mode === "shamrock" && orbit.multi_x_eci) {
      ctx.strokeStyle = COLORS.phase1;
      ctx.globalAlpha = 0.16;
      ctx.lineWidth = 1.2;
      ctx.beginPath();
      for (let k = 0; k < orbit.multi_x_eci.length; k++) {
        const [px, py] = toXY(orbit.multi_x_eci[k], orbit.multi_y_eci[k]);
        if (k === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
      }
      ctx.stroke();
      ctx.globalAlpha = 1.0;
    }

    // Current cycle ECI faint track
    ctx.strokeStyle = mode === "shamrock" ? COLORS.outbound : COLORS.arenstorf;
    ctx.globalAlpha = 0.3;
    ctx.lineWidth = 1.6;
    ctx.beginPath();
    for (let k = 0; k < orbit.x_eci.length; k++) {
      const [px, py] = toXY(orbit.x_eci[k], orbit.y_eci[k]);
      if (k === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
    }
    ctx.stroke();
    ctx.globalAlpha = 1.0;

    // Active trail (clamped to current cycle start to avoid ECI seam jump)
    const segColors = [COLORS.outbound, COLORS.inbound, COLORS.phase1, COLORS.phase2];
    const activeCol = mode === "shamrock" ? segColors[orbit.phase[idxInfo.i]] : COLORS.arenstorf;
    const trailWindowDays = mode === "shamrock" ? 3.5 : 8.0;
    const trailSteps = 45;
    ctx.lineWidth = 2.6;
    ctx.lineCap = "round";
    for (let s = 1; s <= trailSteps; s++) {
      const t0 = tNow - trailWindowDays * (s / trailSteps);
      const t1 = tNow - trailWindowDays * ((s - 1) / trailSteps);
      if (t0 < 0) continue;
      const p0 = sampleAtTime(orbit, t0);
      const p1 = sampleAtTime(orbit, t1);
      const [x0, y0] = toXY(lerp(orbit.x_eci, p0.i, p0.alpha), lerp(orbit.y_eci, p0.i, p0.alpha));
      const [x1, y1] = toXY(lerp(orbit.x_eci, p1.i, p1.alpha), lerp(orbit.y_eci, p1.i, p1.alpha));
      ctx.strokeStyle = activeCol;
      ctx.globalAlpha = 0.9 * (1 - s / trailSteps);
      ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
    }
    ctx.globalAlpha = 1.0;

    // Earth
    ctx.fillStyle = COLORS.earth;
    ctx.beginPath(); ctx.arc(cx, cy, 7.5, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = "#93c5fd";
    ctx.font = "bold 10.5px system-ui, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("Earth", cx, cy + 18);

    // Moving Moon
    const mx = lerp(orbit.xm_eci, idxInfo.i, idxInfo.alpha);
    const my = lerp(orbit.ym_eci, idxInfo.i, idxInfo.alpha);
    const [mpx, mpy] = toXY(mx, my);
    ctx.fillStyle = COLORS.moon;
    ctx.beginPath(); ctx.arc(mpx, mpy, 5.5, 0, Math.PI * 2); ctx.fill();
    ctx.fillStyle = "#e2e8f0";
    ctx.fillText("Moon", mpx + (mx >= 0 ? 18 : -18), mpy + (my >= 0 ? -8 : 14));

    // Spacecraft
    const scx = lerp(orbit.x_eci, idxInfo.i, idxInfo.alpha);
    const scy = lerp(orbit.y_eci, idxInfo.i, idxInfo.alpha);
    const [spx, spy] = toXY(scx, scy);
    ctx.fillStyle = "#ffffff";
    ctx.strokeStyle = activeCol;
    ctx.lineWidth = 2.4;
    ctx.beginPath(); ctx.arc(spx, spy, 5.2, 0, Math.PI * 2); ctx.fill(); ctx.stroke();
  }

  function drawTradeChart() {
    const { ctx, w, h } = resizeCanvas(tradeCanvas);
    const tr = DATA.trade_study;
    const be = tr.breakeven_vs_lean_flights;
    const beHeavy = tr.breakeven_vs_heavy_flights;
    document.getElementById("tradeSub").textContent =
      `Break-even vs Lean 18 t Capsule: ~${be.toFixed(1)} flights · vs 45 t Heavy Hab: ~${beHeavy.toFixed(1)} flights`;
    const pad = { l: 54, r: 18, t: 14, b: 38 };
    const pw = w - pad.l - pad.r;
    const ph = h - pad.t - pad.b;

    ctx.clearRect(0, 0, w, h);
    const maxY = Math.ceil(Math.max(...tr.cum_direct_imleo_t) / 500) * 500;
    const xMap = n => pad.l + ((n - 1) / 24) * pw;
    const yMap = v => pad.t + ph - (v / maxY) * ph;

    // Grid & axes
    ctx.strokeStyle = "rgba(148,163,184,0.16)";
    ctx.fillStyle = COLORS.muted;
    ctx.font = "10.5px ui-monospace, monospace";
    ctx.lineWidth = 1;

    for (let g = 0; g <= 4; g++) {
      const val = (maxY * g) / 4;
      const y = yMap(val);
      ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(w - pad.r, y); ctx.stroke();
      ctx.textAlign = "right";
      ctx.fillText(val.toLocaleString() + " t", pad.l - 6, y + 3.5);
    }
    [1, 5, 10, 15, 20, 25].forEach(n => {
      const x = xMap(n);
      ctx.beginPath(); ctx.moveTo(x, pad.t); ctx.lineTo(x, h - pad.b); ctx.stroke();
      ctx.textAlign = "center";
      ctx.fillText(String(n), x, h - 18);
    });
    ctx.textAlign = "center";
    ctx.fillText("Crewed Lunar Round-Trips", pad.l + pw * 0.5, h - 3);

    function plotLine(vals, color, dashed) {
      ctx.strokeStyle = color;
      ctx.lineWidth = 2.4;
      if (dashed) ctx.setLineDash([6, 4]);
      ctx.beginPath();
      vals.forEach((v, idx) => {
        const x = xMap(tr.flights[idx]);
        const y = yMap(v);
        if (idx === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      });
      ctx.stroke();
      ctx.setLineDash([]);
    }

    plotLine(tr.cum_direct_imleo_t, COLORS.arenstorf, false);
    plotLine(tr.cum_apollo_lean_imleo_t, COLORS.phase1, true);
    plotLine(tr.cum_cycler_imleo_t, COLORS.inbound, false);

    // Annotate break-even vs lean Apollo/Orion capsule
    if (isFinite(be) && be <= 25) {
      const bx = xMap(be);
      const by = yMap(tr.imleo_apollo_lean_t * be);
      ctx.strokeStyle = COLORS.phase1;
      ctx.setLineDash([3, 3]);
      ctx.lineWidth = 1.2;
      ctx.beginPath(); ctx.moveTo(bx, pad.t); ctx.lineTo(bx, h - pad.b); ctx.stroke();
      ctx.setLineDash([]);

      ctx.fillStyle = COLORS.phase1;
      ctx.beginPath(); ctx.arc(bx, by, 4.5, 0, Math.PI * 2); ctx.fill();
      ctx.fillStyle = COLORS.fg;
      ctx.font = "bold 11px system-ui, sans-serif";
      ctx.textAlign = "left";
      ctx.fillText(`Break-even vs Lean Capsule: ~${be.toFixed(1)} flights`, bx + 8, by + 16);
    }
  }

  function renderFrame() {
    const orbit = DATA[mode];
    const period = orbit.metrics.period_days;
    const tNow = frac * period;
    const idxInfo = sampleAtTime(orbit, tNow);

    const altE = lerp(orbit.alt_earth_km, idxInfo.i, idxInfo.alpha);
    const altM = lerp(orbit.alt_moon_km, idxInfo.i, idxInfo.alpha);
    const spd = lerp(orbit.v_eci_km_s, idxInfo.i, idxInfo.alpha);
    const cj = lerp(orbit.jacobi, idxInfo.i, idxInfo.alpha);

    dayReadout.textContent = `Day ${tNow.toFixed(2).padStart(5, "0")} / ${period.toFixed(2)}`;
    valEarthAlt.textContent = `${Math.round(altE).toLocaleString()} km`;
    valMoonAlt.textContent = `${Math.round(altM).toLocaleString()} km`;
    valSpeed.textContent = `${spd.toFixed(2)} km/s`;
    valJacobi.textContent = cj.toFixed(4);
    slider.value = Math.round(frac * 1000);

    if (mode === "shamrock") {
      const ph = orbit.phase[idxInfo.i];
      const meta = [
        ["Outbound Fig-8 (Earth → Moon)", COLORS.outbound],
        ["Inbound Fig-8 (Moon → Earth)", COLORS.inbound],
        ["Phasing Loop 1 (Sub-Lunar)", COLORS.phase1],
        ["Phasing Loop 2 (Sub-Lunar)", COLORS.phase2],
      ][ph];
      let label = meta[0];
      let col = meta[1];
      if (Math.abs(tNow - orbit.metrics.burn1_time_days) < 0.45) {
        label = "★ Burn 1: −18.7 m/s (Oberth)";
        col = COLORS.burn;
      } else if (tNow < 0.35 || Math.abs(tNow - period) < 0.45) {
        label = "★ Burn 2: +18.7 m/s (Re-inject)";
        col = COLORS.burn;
      } else if (Math.abs(tNow - orbit.metrics.perilune_time_days) < 0.35) {
        label = "★ Lunar Perilune (3,000 km)";
        col = COLORS.inbound;
      }
      phaseDot.style.background = col;
      phaseText.textContent = label;
    } else {
      phaseDot.style.background = COLORS.arenstorf;
      phaseText.textContent = "Ballistic Arenstorf (Δv = 0)";
    }

    drawRotating(orbit, idxInfo, tNow);
    drawECI(orbit, idxInfo, tNow);
  }

  function tick(ts) {
    if (lastTs !== null && playing) {
      const dt = Math.min((ts - lastTs) / 1000, 0.1);
      const cycleSeconds = mode === "shamrock" ? 11.0 : 16.0;
      frac = (frac + (dt * speed) / cycleSeconds) % 1.0;
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

  btnShamrock.addEventListener("click", () => {
    mode = "shamrock"; frac = 0;
    btnShamrock.classList.add("active");
    btnArenstorf.classList.remove("active");
    updateLegend(); renderFrame();
  });

  btnArenstorf.addEventListener("click", () => {
    mode = "arenstorf"; frac = 0;
    btnArenstorf.classList.add("active");
    btnShamrock.classList.remove("active");
    updateLegend(); renderFrame();
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
    drawTradeChart();
  });

  updateLegend();
  renderFrame();
  drawTradeChart();
  requestAnimationFrame(tick);
})();
</script>
</body>
</html>
"""


def build():
    root = Path(__file__).resolve().parent.parent
    data_path = root / "data" / "cycler_solution.json"
    data_text = data_path.read_text(encoding="utf-8")
    # Compact JSON for fast inline load
    compact = json.dumps(json.loads(data_text), separators=(",", ":"))
    html = HTML_TEMPLATE.replace("__CYCLER_DATA_JSON__", compact)
    out_html = root / "web" / "index.html"
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html, encoding="utf-8")
    print(f"Built {out_html} ({len(html):,} bytes)")


if __name__ == "__main__":
    build()
