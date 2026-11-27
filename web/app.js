"use strict";

// The page runs the real bioseq Python package in the browser with Pyodide.
// web/bridge.py turns each analysis into JSON; this file only draws it.

const PYODIDE_URL = "https://cdn.jsdelivr.net/pyodide/v0.26.4/full/";
const BIOSEQ_FILES = ["__init__.py", "fasta.py", "sequence.py", "motifs.py", "properties.py", "composition.py"];
const SEQ_DISPLAY_LIMIT = 12000;   // longer sequences are shown truncated (copy/download stay full)
const TABLE_ROW_LIMIT = 500;       // same idea for long tables

const TOOL_FOR_TAB = {
  overview: "basics", translation: "translation", orfs: "orfs", codons: "codons",
  motif: "motif", kmers: "kmers", composition: "composition", properties: "properties",
  quality: "quality", report: "report",
};
const DATASET_TABS = ["quality", "report"];  // tabs about all records, not the selected one
const RECORD_LIST_LIMIT = 200;

const state = {
  records: [],
  selected: null,
  tab: "overview",
  minOrf: "30",
  frame: "1",
  toStop: false,
  motif: "GAATTC",
  k: "3",
  top: "10",
  window: "100",
  na: "50",
  oligo: "250",
  trimQ: "20",
  trimWindow: "4",
  minLength: "36",
  format: "fasta",
};

let bridge = null;
let copies = [];
let downloads = {};

const $ = (sel, root = document) => root.querySelector(sel);
const $$ = (sel, root = document) => [...root.querySelectorAll(sel)];

// ----------------------------------------------------------------- helpers

function esc(value) {
  return String(value).replace(/[&<>"']/g, c => (
    { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
}

function fmt(n, digits = 0) {
  return Number(n).toLocaleString("en-US", { minimumFractionDigits: digits, maximumFractionDigits: digits });
}

function plural(n, word, many = word + "s") {
  return `${fmt(n)} ${n === 1 ? word : many}`;
}

function safeName(id) {
  return String(id).replace(/[^\w.-]+/g, "_") || "sequence";
}

function debounce(fn, ms) {
  let timer;
  return (...args) => { clearTimeout(timer); timer = setTimeout(() => fn(...args), ms); };
}

function call(fn, ...args) {
  return JSON.parse(bridge[fn](...args));
}

function toCsv(rows) {
  return rows.map(row => row.map(value => {
    const s = value == null ? "" : String(value);
    return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
  }).join(",")).join("\n") + "\n";
}

function toFasta(entries) {
  return entries.map(([id, seq]) =>
    `>${id}\n${(seq.match(/.{1,60}/g) || [""]).join("\n")}`).join("\n") + "\n";
}

function saveFile(name, text, type) {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const a = Object.assign(document.createElement("a"), { href: url, download: name });
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function copyButton(text, label = "Copy") {
  copies.push(text);
  return `<button class="small" data-copy="${copies.length - 1}">${label}</button>`;
}

function downloadButton(label, makeFile) {
  const key = String(Object.keys(downloads).length);
  downloads[key] = makeFile;
  return `<button class="small" data-download="${key}">${esc(label)}</button>`;
}

function flash(button, text) {
  const original = button.textContent;
  button.textContent = text;
  setTimeout(() => { button.textContent = original; }, 1200);
}

// ------------------------------------------------------------ status & msgs

function setStatus(kind, label) {
  const el = $("#status");
  el.className = `status ${kind}`;
  $(".label", el).textContent = label;
}

function showInputMsg(kind, text) {
  const el = $("#input-msg");
  el.hidden = !text;
  el.className = `msg ${kind}`;
  el.textContent = text || "";
}

function setEngineReady(ready) {
  $$("#btn-analyze, #btn-example").forEach(b => { b.disabled = !ready; });
}

// ------------------------------------------------------------------ start-up

async function fetchText(path) {
  const response = await fetch(path, { cache: "no-cache" });
  if (!response.ok) throw new Error(`could not fetch ${path} (${response.status})`);
  return response.text();
}

async function start() {
  bindEvents();
  const fromHash = location.hash.slice(1);
  if (fromHash in TOOL_FOR_TAB) state.tab = fromHash;
  updateTabs();
  if (location.protocol === "file:") {
    setStatus("error", "Needs a web server");
    $("#output").innerHTML = `<div class="empty-state"><strong>Open this page through a web server</strong>
      Browsers block it from loading the Python files when opened directly from disk.
      In the project folder run <code>python3 main.py --web</code>.</div>`;
    return;
  }
  setStatus("loading", "Loading engine…");
  try {
    if (typeof loadPyodide !== "function") throw new Error("the Pyodide script did not load");
    const [pyodide, sources, bridgeSource] = await Promise.all([
      loadPyodide({ indexURL: PYODIDE_URL }),
      Promise.all(BIOSEQ_FILES.map(f => fetchText(`../bioseq/${f}`))),
      fetchText("bridge.py"),
    ]);
    pyodide.FS.mkdirTree("/app/bioseq");
    BIOSEQ_FILES.forEach((f, i) => pyodide.FS.writeFile(`/app/bioseq/${f}`, sources[i]));
    pyodide.FS.writeFile("/app/bridge.py", bridgeSource);
    pyodide.runPython("import sys; sys.path.insert(0, '/app')");
    bridge = pyodide.pyimport("bridge");

    renderMotifChips(call("info").motifs);
    setStatus("ready", "Ready");
    setEngineReady(true);
    if ($("#seq-input").value.trim()) analyzeInput();
    else await loadExample();
  } catch (err) {
    console.error(err);
    setStatus("error", "Engine failed to load");
    $("#output").innerHTML = `<div class="empty-state"><strong>The analysis engine could not start</strong>
      ${esc(err.message || err)}.<br>Check your internet connection and reload the page.</div>`;
  }
}

// --------------------------------------------------------------- input side

async function loadExample() {
  try {
    $("#seq-input").value = await fetchText("../data/example.fasta");
    analyzeInput();
  } catch (err) {
    showInputMsg("error", `Could not load the example: ${err.message}`);
  }
}

function analyzeInput() {
  if (!bridge) return;
  const result = call("load_text", $("#seq-input").value);
  if (!result.ok) {
    showInputMsg("error", result.error);
    return;
  }
  const wasFastq = state.format === "fastq";
  state.records = result.records;
  state.format = result.format;
  if (!state.records.some(r => r.id === state.selected)) state.selected = state.records[0].id;
  const rna = state.records.filter(r => r.rna).length;
  const loaded = result.format === "fastq"
    ? `Loaded ${plural(state.records.length, "read")} from FASTQ.` +
      (result.truncated ? ` Only the first ${fmt(result.max_reads)} reads were loaded.` : "")
    : `Loaded ${plural(state.records.length, "sequence")}.`;
  showInputMsg(result.truncated ? "warn" : "ok", loaded +
    (rna ? ` ${rna === state.records.length ? "RNA was" : `${plural(rna, "RNA sequence")} were`} read as DNA (U → T).` : ""));
  if (result.format === "fastq" && !wasFastq) state.tab = "quality";
  renderRecords();
  renderTool();
}

async function readFile(file) {
  if (!file) return;
  try {
    $("#seq-input").value = await file.text();
    analyzeInput();
  } catch (err) {
    showInputMsg("error", `Could not read ${file.name}: ${err.message}`);
  }
}

function renderRecords() {
  const list = $("#record-list");
  $("#record-count").textContent = state.records.length ? plural(state.records.length, "sequence") : "";
  if (!state.records.length) {
    list.innerHTML = `<li class="muted empty">Nothing loaded yet.</li>`;
    return;
  }
  const shown = state.records.slice(0, RECORD_LIST_LIMIT);
  if (!shown.some(r => r.id === state.selected)) {
    const sel = state.records.find(r => r.id === state.selected);
    if (sel) shown.push(sel);
  }
  list.innerHTML = shown.map(r => `
    <li><button data-id="${esc(r.id)}" class="${r.id === state.selected ? "active" : ""}">
      <span class="rec-id">${esc(r.id)}</span>
      <span class="rec-meta">${fmt(r.length)} ${r.rna ? "nt · RNA" : "bp"} · ${fmt(r.gc, 1)}% GC</span>
    </button></li>`).join("") +
    (state.records.length > shown.length
      ? `<li class="muted empty">…and ${fmt(state.records.length - shown.length)} more. All of them are included in Read quality and All sequences.</li>`
      : "");
}

function renderMotifChips(motifs) {
  $("#motif-chips").innerHTML = Object.entries(motifs).map(([name, motif]) =>
    `<button type="button" data-motif="${esc(name)}">${esc(name)} <span class="chip-seq">${esc(motif)}</span></button>`
  ).join("");
}

// ------------------------------------------------------------- tool output

function params() {
  return {
    min_orf: state.minOrf, frame: state.frame, to_stop: state.toStop,
    motif: state.motif, k: state.k, top: state.top,
    window: state.window, na: state.na, oligo: state.oligo,
    trim_q: state.trimQ, trim_window: state.trimWindow, min_length: state.minLength,
  };
}

function updateTabs() {
  if (location.hash.slice(1) !== state.tab) history.replaceState(null, "", `#${state.tab}`);
  $$(".tabs [data-tab]").forEach(t => {
    const active = t.dataset.tab === state.tab;
    t.setAttribute("aria-selected", String(active));
    if (active && t.scrollIntoView) t.scrollIntoView({ block: "nearest", inline: "nearest" });
  });
  $$(".controls").forEach(c => { c.hidden = !c.dataset.for.split(" ").includes(state.tab); });
}

function renderTool() {
  updateTabs();
  const out = $("#output");
  if (!bridge) return;
  copies = [];
  downloads = {};
  if (!state.records.length) {
    out.innerHTML = `<div class="empty-state"><strong>No sequences loaded</strong>
      Paste a sequence on the left and press Analyze, or try the example.</div>`;
    return;
  }
  const result = call("run", TOOL_FOR_TAB[state.tab], state.selected || "", JSON.stringify(params()));
  if (!result.ok) {
    out.innerHTML = `<div class="msg error">${esc(result.error)}</div>`;
    return;
  }
  out.innerHTML = RENDER[state.tab](result);
}

function head(title, actions = "") {
  const rec = state.records.find(r => r.id === state.selected);
  const sub = DATASET_TABS.includes(state.tab) ? "" :
    `<span class="sub">${esc(rec.id)} · ${fmt(rec.length)} bp</span>`;
  return `<div class="out-head"><h2>${esc(title)}${sub}</h2><div class="actions">${actions}</div></div>`;
}

function formatSeq(seq, { stops = false } = {}) {
  if (!seq) return `<pre class="seq"><span class="pos">(empty)</span></pre>`;
  const shown = seq.slice(0, SEQ_DISPLAY_LIMIT);
  const width = String(shown.length).length;
  const lines = [];
  for (let i = 0; i < shown.length; i += 60) {
    let chunk = esc(shown.slice(i, i + 60).match(/.{1,10}/g).join(" "));
    if (stops) chunk = chunk.replace(/\*/g, '<span class="stop">*</span>');
    lines.push(`<span class="pos">${String(i + 1).padStart(width)}</span>  ${chunk}`);
  }
  const more = seq.length > shown.length
    ? `<div class="truncated">Showing the first ${fmt(shown.length)} of ${fmt(seq.length)} characters. Copy or download for the full sequence.</div>`
    : "";
  return `<pre class="seq">${lines.join("\n")}</pre>${more}`;
}

function seqBlock(title, sub, seq, fastaId, opts) {
  const actions = copyButton(seq) +
    downloadButton("FASTA", () => [`${safeName(fastaId)}.fasta`, toFasta([[fastaId, seq]]), "text/plain"]);
  return `<div class="seq-block">
    <div class="seq-head"><h3>${esc(title)} <span class="sub">${esc(sub)}</span></h3><div class="actions">${actions}</div></div>
    ${formatSeq(seq, opts)}
  </div>`;
}

function tile(label, value, unit = "", note = "") {
  return `<div class="tile"><div class="t-label">${label}</div>
    <div class="t-value">${value}<span class="t-unit">${unit}</span></div>
    ${note ? `<div class="t-note">${note}</div>` : ""}</div>`;
}

function strandBadge(strand) {
  return `<span class="strand ${strand === "+" ? "plus" : "minus"}">${strand === "+" ? "+" : "−"}</span>`;
}

function table(headers, rows, { rowAttrs } = {}) {
  const shown = rows.slice(0, TABLE_ROW_LIMIT);
  const ths = headers.map(h => `<th class="${h.num ? "num" : ""}">${h.label}</th>`).join("");
  const trs = shown.map((cells, i) =>
    `<tr ${rowAttrs ? rowAttrs(i) : ""}>${cells.map((c, j) =>
      `<td class="${headers[j].cls || (headers[j].num ? "num" : "")}">${c}</td>`).join("")}</tr>`).join("");
  const more = rows.length > shown.length
    ? `<p class="note">Showing the first ${fmt(shown.length)} of ${fmt(rows.length)} rows. Download the CSV for all of them.</p>`
    : "";
  return `<div class="table-wrap"><table><thead><tr>${ths}</tr></thead><tbody>${trs}</tbody></table></div>${more}`;
}

function barCell(fraction, label) {
  return `<div class="bar-cell"><span class="bar"><span style="width:${(fraction * 100).toFixed(1)}%"></span></span>${label}</div>`;
}

function niceStep(raw) {
  const power = 10 ** Math.floor(Math.log10(Math.max(raw, 1)));
  return [1, 2, 5, 10].map(m => m * power).find(step => step >= raw);
}

// Linear map of features along the sequence: + strand above the axis, - below.
function featureMap(length, features, { ticks = false } = {}) {
  if (!length) return "";
  const W = Math.max(300, Math.floor($("#output").clientWidth - 40));
  const pad = 14, laneH = 16, featH = 11, maxLanes = 6;
  const scale = (W - 2 * pad) / length;
  const x = pos => pad + pos * scale;

  const pack = strand => {
    const laneEnds = [];
    const placed = features.filter(f => f.strand === strand).map(f => {
      const x0 = x(f.start - 1);
      const x1 = Math.max(x(f.end), x0 + (ticks ? 2 : 3));
      let lane = ticks ? 0 : laneEnds.findIndex(end => end + 3 < x0);
      if (lane < 0) lane = laneEnds.length < maxLanes ? laneEnds.length : maxLanes - 1;
      laneEnds[lane] = Math.max(laneEnds[lane] ?? 0, x1);
      return { f, x0, x1, lane };
    });
    return { placed, lanes: Math.max(1, laneEnds.length) };
  };
  const plus = pack("+");
  const minus = pack("-");

  const axisY = 6 + plus.lanes * laneH + 6;
  const minusTop = axisY + 24;
  const height = minusTop + minus.lanes * laneH + 4;

  const shape = ({ f, x0, x1, lane }, strand) => {
    const y = strand === "+"
      ? axisY - 6 - (lane + 1) * laneH + (laneH - featH) / 2
      : minusTop + lane * laneH + (laneH - featH) / 2;
    const cls = `feat ${strand === "+" ? "plus" : "minus"}`;
    const title = `<title>${esc(f.label)}</title>`;
    if (ticks) {
      return `<g>${title}<rect class="${cls}" x="${x0.toFixed(1)}" y="${y}" width="${(x1 - x0).toFixed(1)}" height="${featH}" rx="1"/></g>`;
    }
    const w = x1 - x0, tip = Math.min(7, w / 2), b = y + featH, m = y + featH / 2;
    const pts = strand === "+"
      ? [[x0, y], [x1 - tip, y], [x1, m], [x1 - tip, b], [x0, b]]
      : [[x1, y], [x0 + tip, y], [x0, m], [x0 + tip, b], [x1, b]];
    return `<g>${title}<polygon class="${cls}" points="${pts.map(p => p.map(v => v.toFixed(1)).join(",")).join(" ")}"/></g>`;
  };

  const step = niceStep(length / Math.max(2, Math.floor(W / 110)));
  const tickMarks = [1];
  for (let p = step; p <= length; p += step) tickMarks.push(p);
  const axis = tickMarks.map(p => {
    const tx = x(p === 1 ? 0 : p);
    const anchor = tx > W - 30 ? "end" : tx < 30 ? "start" : "middle";
    return `<line class="tick" x1="${tx}" x2="${tx}" y1="${axisY}" y2="${axisY + 4}"/>
      <text x="${tx}" y="${axisY + 16}" text-anchor="${anchor}">${fmt(p)}</text>`;
  }).join("");

  return `<div class="map">
    <div class="legend"><span><span class="swatch plus"></span>+ strand</span><span><span class="swatch minus"></span>− strand</span><span>Hover for details</span></div>
    <svg viewBox="0 0 ${W} ${height}" width="${W}" height="${height}" role="img" aria-label="Feature map">
      <line class="axis" x1="${pad}" x2="${W - pad}" y1="${axisY}" y2="${axisY}"/>
      ${axis}
      ${plus.placed.map(p => shape(p, "+")).join("")}
      ${minus.placed.map(p => shape(p, "-")).join("")}
    </svg></div>`;
}

// Line chart of one value along the sequence, with optional shaded regions.
function lineChart(length, points, { label, min, max, digits = 0, unit = "", zero = false, bands = [] }) {
  if (!points.length) return "";
  const W = Math.max(300, Math.floor($("#output").clientWidth - 40));
  const H = 120, left = 44, right = 10, top = 8, bottom = 22;
  const x = pos => left + (pos / length) * (W - left - right);
  const y = v => top + (1 - (v - min) / (max - min)) * (H - top - bottom);
  const path = points.map((p, i) => `${i ? "L" : "M"}${x(p.pos).toFixed(1)},${y(p.value).toFixed(1)}`).join("");
  const shaded = bands.map(b =>
    `<rect class="band" x="${x(b.start - 1).toFixed(1)}" y="${top}" width="${Math.max(2, x(b.end) - x(b.start - 1)).toFixed(1)}" height="${H - top - bottom}"><title>${esc(b.label)}</title></rect>`).join("");
  const gridValues = zero ? [min, 0, max] : [min, (min + max) / 2, max];
  const grid = gridValues.map(v => `<line class="${zero && v === 0 ? "zero" : "grid"}" x1="${left}" x2="${W - right}" y1="${y(v)}" y2="${y(v)}"/>
    <text x="${left - 6}" y="${y(v) + 4}" text-anchor="end">${fmt(v, digits)}${unit}</text>`).join("");
  const step = niceStep(length / Math.max(2, Math.floor(W / 110)));
  const ticks = [];
  for (let p = step; p < length; p += step) ticks.push(p);
  const axis = ticks.map(p => `<text x="${x(p)}" y="${H - 6}" text-anchor="middle">${fmt(p)}</text>`).join("");
  return `<div class="chart"><div class="chart-label">${label}</div>
    <svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" role="img" aria-label="${esc(label)} along the sequence">
      ${shaded}${grid}${axis}<path class="line" d="${path}"/>
    </svg></div>`;
}

// Per-position quality box plot in the style of FastQC: whiskers 10th-90th
// percentile, box 25th-75th, median tick, and a line through the means.
function qualityPlot(positions, maxLength) {
  const W = Math.max(300, Math.floor($("#output").clientWidth - 40));
  const H = 230, left = 40, right = 10, top = 8, bottom = 26;
  const qMax = Math.max(41, ...positions.map(p => p.p90));
  const plotW = W - left - right;
  const x = pos => left + (pos / maxLength) * plotW;
  const y = q => top + (1 - q / qMax) * (H - top - bottom);
  const zones = [["bad", 0, 20], ["ok", 20, 28], ["good", 28, qMax]].map(([cls, lo, hi]) =>
    `<rect class="zone ${cls}" x="${left}" y="${y(hi)}" width="${plotW}" height="${y(lo) - y(hi)}"/>`).join("");
  const grid = [0, 10, 20, 30, 40].filter(q => q <= qMax).map(q =>
    `<line class="grid" x1="${left}" x2="${W - right}" y1="${y(q)}" y2="${y(q)}"/>
     <text x="${left - 6}" y="${y(q) + 4}" text-anchor="end">${q}</text>`).join("");
  const boxes = positions.map(p => {
    const x0 = x(p.start - 1), x1 = x(p.end), w = x1 - x0;
    const bx = x0 + w * 0.15, bw = Math.max(1, w * 0.7), cx = x0 + w / 2;
    const tip = `<title>${p.start === p.end ? `Base ${fmt(p.start)}` : `Bases ${fmt(p.start)}–${fmt(p.end)}`}: mean Q${fmt(p.mean, 1)}, median Q${p.median}, middle half Q${p.q1}–Q${p.q3} · ${plural(p.reads, "read")}</title>`;
    return `<g>${tip}
      <line class="whisker" x1="${cx}" x2="${cx}" y1="${y(p.p90)}" y2="${y(p.p10)}"/>
      <rect class="box" x="${bx}" y="${y(p.q3)}" width="${bw}" height="${Math.max(1, y(p.q1) - y(p.q3))}"/>
      <line class="median" x1="${bx}" x2="${bx + bw}" y1="${y(p.median)}" y2="${y(p.median)}"/></g>`;
  }).join("");
  const mean = positions.map((p, i) => `${i ? "L" : "M"}${(x(p.start - 1) + x(p.end)) / 2},${y(p.mean)}`).join("");
  const step = niceStep(maxLength / Math.max(2, Math.floor(W / 90)));
  const ticks = [1];
  for (let p = step; p <= maxLength; p += step) ticks.push(p);
  const axis = ticks.map(p => `<text x="${x(p === 1 ? 0.5 : p - 0.5)}" y="${H - 8}" text-anchor="middle">${fmt(p)}</text>`).join("");
  return `<div class="chart quality-plot"><div class="chart-label">Quality along the reads (Phred score by base position)</div>
    <div class="legend"><span><span class="swatch box"></span>Middle 50% of reads</span><span><span class="swatch whisker"></span>10th–90th percentile</span><span><span class="swatch mean"></span>Mean</span><span>Hover for details</span></div>
    <svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" role="img" aria-label="Quality by position in the read">
      ${zones}${grid}${boxes}<path class="mean" d="${mean}"/>${axis}
    </svg></div>`;
}

function histogram(bars, { label, xLabel }) {
  const W = Math.max(300, Math.floor($("#output").clientWidth - 40));
  const H = 140, left = 40, right = 10, top = 8, bottom = 26;
  const max = Math.max(1, ...bars.map(b => b.value));
  const bw = (W - left - right) / bars.length;
  const y = v => top + (1 - v / max) * (H - top - bottom);
  const rects = bars.map((b, i) => `<g><title>${esc(b.title)}</title>
    <rect class="hbar ${b.cls || ""}" x="${left + i * bw + 1}" y="${y(b.value)}" width="${Math.max(1, bw - 2)}" height="${H - bottom - y(b.value)}"/></g>`).join("");
  const every = Math.max(1, Math.ceil(bars.length / Math.floor((W - left) / 32)));
  const labels = bars.map((b, i) => i % every ? "" :
    `<text x="${left + (i + 0.5) * bw}" y="${H - 8}" text-anchor="middle">${esc(b.label)}</text>`).join("");
  return `<div class="chart"><div class="chart-label">${label}</div>
    <svg viewBox="0 0 ${W} ${H}" width="${W}" height="${H}" role="img" aria-label="${esc(label)}">
      <line class="grid" x1="${left}" x2="${W - right}" y1="${y(max)}" y2="${y(max)}"/>
      <text x="${left - 6}" y="${y(max) + 4}" text-anchor="end">${fmt(max)}</text>
      <text x="${left - 6}" y="${H - bottom + 4}" text-anchor="end">0</text>
      ${rects}${labels}
    </svg><div class="axis-label">${xLabel}</div></div>`;
}

const RENDER = {
  overview(r) {
    const bases = ["A", "C", "G", "T"];
    const counts = { ...r.counts, N: r.n };
    const shownBases = r.n ? [...bases, "N"] : bases;
    const pct = b => (r.length ? (100 * counts[b] / r.length) : 0);
    const tiles = [
      tile("Length", fmt(r.length), "bp"),
      tile("GC content", fmt(r.gc, 1), "%"),
      ...shownBases.map(b => tile(`<span class="swatch b-${b}"></span>${b}`, fmt(counts[b]), "", `${fmt(pct(b), 1)}%`)),
    ].join("");
    const bar = shownBases.map(b =>
      `<span class="b-${b}" style="flex:${counts[b]}" title="${b}: ${fmt(pct(b), 1)}%"></span>`).join("");
    const rec = state.selected;
    return head("Overview") +
      `<div class="tiles">${tiles}</div>
      <div class="composition" aria-hidden="true">${bar}</div>
      <div class="section">
        ${seqBlock("Reverse complement", "5′→3′, opposite strand", r.revcomp, `${rec}_revcomp`)}
        ${seqBlock("mRNA", "transcribed, T → U", r.mrna, `${rec}_mRNA`)}
      </div>`;
  },

  translation(r) {
    const rec = state.selected;
    if (r.frames) {
      const blocks = Object.entries(r.frames).map(([frame, protein]) =>
        seqBlock(`Frame ${frame.replace("-", "−")}`, `${fmt(protein.length)} aa`, protein, `${rec}_frame${frame}`, { stops: true })
      ).join("");
      const all = downloadButton("All frames (FASTA)", () =>
        [`${safeName(rec)}_six_frames.fasta`, toFasta(Object.entries(r.frames).map(([f, p]) => [`${rec}_frame${f}`, p])), "text/plain"]);
      return head("Six-frame translation", all) +
        `<p class="summary">Frames +1 to +3 read the sequence as given, starting at bases 1, 2 and 3. Frames −1 to −3 do the same on the reverse complement. <span class="stop-key">Stop codons are shown as <strong style="color:var(--stop)">*</strong>.</span></p>${blocks}`;
    }
    return head("Translation") +
      `<p class="summary">Standard genetic code, reading frame +${r.frame}${state.toStop ? ", stopping at the first stop codon" : ""}. Stop codons are shown as <strong style="color:var(--stop)">*</strong>, and codons containing N as X.</p>` +
      seqBlock(`Protein, frame +${r.frame}`, `${fmt(r.protein.length)} aa`, r.protein, `${rec}_frame${r.frame}`, { stops: true });
  },

  orfs(r) {
    const rec = state.selected;
    if (!r.orfs.length) {
      return head("Open reading frames") +
        `<div class="empty-state"><strong>No ORFs of at least ${fmt(r.min_orf)} nt</strong>Try a lower minimum length above.</div>`;
    }
    const csv = downloadButton("CSV", () => [`${safeName(rec)}_orfs.csv`, toCsv([
      ["start", "end", "strand", "length_nt", "length_aa", "mol_weight_da", "protein"],
      ...r.orfs.map(o => [o.start, o.end, o.strand, o.length, o.aa, o.mw, o.protein]),
    ]), "text/csv"]);
    const fasta = downloadButton("Proteins (FASTA)", () => [`${safeName(rec)}_orf_proteins.fasta`,
      toFasta(r.orfs.map(o => [`${rec}_orf_${o.start}-${o.end}(${o.strand}) length=${o.aa}aa`, o.protein])), "text/plain"]);
    const features = r.orfs.map(o => ({
      start: o.start, end: o.end, strand: o.strand,
      label: `ORF ${o.start}–${o.end} (${o.strand}) · ${o.length} nt · ${o.aa} aa`,
    }));
    const rows = r.orfs.map((o, i) => [
      i + 1,
      `${fmt(o.start)}–${fmt(o.end)}`,
      strandBadge(o.strand),
      fmt(o.length),
      fmt(o.aa),
      o.mw == null ? "–" : fmt(o.mw / 1000, 2),
      o.protein.length > 60
        ? `<details><summary>${esc(o.protein.slice(0, 40))}…</summary>${esc(o.protein)}</details>`
        : esc(o.protein),
    ]);
    return head("Open reading frames", csv + fasta) +
      `<p class="summary"><strong>${plural(r.orfs.length, "ORF")}</strong> of at least ${fmt(r.min_orf)} nt, searched on both strands. Positions refer to the forward strand.</p>` +
      featureMap(r.length, features) +
      table([
        { label: "#", num: true }, { label: "Position" }, { label: "Strand" },
        { label: "Length (nt)", num: true }, { label: "Protein (aa)", num: true },
        { label: "Mol. weight (kDa)", num: true }, { label: "Protein sequence", cls: "protein" },
      ], rows);
  },

  codons(r) {
    const rec = state.selected;
    if (!r.n_orfs) {
      return head("Codon usage") +
        `<div class="empty-state"><strong>No ORFs of at least ${fmt(r.min_orf)} nt</strong>Codon usage is counted inside ORFs. Try a lower minimum length above.</div>`;
    }
    const csv = downloadButton("CSV", () => [`${safeName(rec)}_codon_usage.csv`, toCsv([
      ["amino_acid", "codon", "count", "fraction_of_amino_acid"],
      ...r.groups.flatMap(g => g.codons.map(c => [g.name, c.codon, c.count, c.fraction == null ? "" : c.fraction.toFixed(4)])),
    ]), "text/csv"]);
    const rows = r.groups.map(g => `<tr>
      <th><span class="aa3">${g.name}</span><span class="aa1">${g.aa}</span><div class="aa-total">${plural(g.total, "codon")}</div></th>
      <td><div class="codon-cells">${g.codons.map(c => `
        <div class="codon-cell ${c.count ? "" : "zero"}" title="${c.codon}: ${c.count} of ${g.total} ${g.name} codons">
          <span class="codon">${c.codon}</span><span class="count">${c.count}</span>
          <span class="bar"><span style="width:${((c.fraction || 0) * 100).toFixed(1)}%"></span></span>
          <span class="frac">${c.fraction == null ? "–" : Math.round(c.fraction * 100) + "%"}</span>
        </div>`).join("")}</div></td></tr>`).join("");
    return head("Codon usage", csv) +
      `<p class="summary"><strong>${plural(r.total, "codon")}</strong> counted in ${plural(r.n_orfs, "ORF")} of at least ${fmt(r.min_orf)} nt, stop codons included. The bars show each codon's share among the codons for the same amino acid.</p>
      <div class="table-wrap"><table class="codon-table"><tbody>${rows}</tbody></table></div>`;
  },

  motif(r) {
    const rec = state.selected;
    const label = r.query.toUpperCase() === r.motif ? r.motif : `${r.query} (${r.motif})`;
    if (!r.total) {
      return head("Motif search") +
        `<div class="empty-state"><strong>No matches for ${esc(label)}</strong>Searched both strands of ${esc(rec)}. Try another motif, or pick one of the sites above.</div>`;
    }
    const plusHits = r.hits.filter(h => h.strand === "+").length;
    const csv = downloadButton("CSV", () => [`${safeName(rec)}_${safeName(r.motif)}_hits.csv`, toCsv([
      ["start", "end", "strand", "site", ...(r.enzyme ? ["cut_after"] : [])],
      ...r.hits.map(h => [h.start, h.end, h.strand, h.site, ...(r.enzyme ? [h.cut] : [])]),
    ]), "text/csv"]);
    const enz = r.enzyme;
    const cutText = h => `${fmt(h.cut)}^${fmt(h.cut + 1)}`;
    const features = r.hits.map(h => ({ ...h,
      label: `${h.site} at ${h.start}–${h.end} (${h.strand})${enz ? ` · cuts between ${cutText(h)}` : ""}` }));
    const rows = r.hits.map((h, i) => [i + 1, `${fmt(h.start)}–${fmt(h.end)}`, strandBadge(h.strand), esc(h.site),
      ...(enz ? [cutText(h)] : [])]);
    const capped = r.total > r.hits.length
      ? ` Only the first ${fmt(r.hits.length)} are listed.` : "";
    const enzNote = enz
      ? `<p class="note"><strong>${esc(enz.name)}</strong> cuts <span class="mono">${esc(enz.site.replace("^", " ^ "))}</span> and leaves ${esc(enz.ends === "blunt" ? "blunt ends" : `a ${enz.ends}`)}. <strong>Cut</strong> gives the two forward-strand bases the cut falls between: + rows show the top-strand cut, − rows the bottom-strand cut.</p>`
      : "";
    return head("Motif search", csv) +
      `<p class="summary"><strong>${plural(r.total, "match", "matches")}</strong> for <strong>${esc(label)}</strong>: ${fmt(plusHits)} on the + strand and ${fmt(r.hits.length - plusHits)} on the − strand.${capped} Sites that read the same on both strands, like GAATTC, appear once per strand.</p>` +
      enzNote +
      featureMap(r.length, features, { ticks: true }) +
      table([{ label: "#", num: true }, { label: "Position" }, { label: "Strand" }, { label: "Site (5′→3′)", cls: "mono" },
        ...(enz ? [{ label: "Cut", cls: "mono" }] : [])], rows);
  },

  kmers(r) {
    const rec = state.selected;
    if (!r.total) {
      return head("K-mer counts") +
        `<div class="empty-state"><strong>The sequence is shorter than k = ${r.k}</strong>Choose a smaller k above.</div>`;
    }
    const csv = downloadButton("CSV", () => [`${safeName(rec)}_${r.k}-mers.csv`, toCsv([
      ["kmer", "count", "percent"], ...r.top.map(t => [t.kmer, t.count, t.percent.toFixed(3)]),
    ]), "text/csv"]);
    const possible = r.k <= 12 ? ` out of ${fmt(4 ** r.k)} possible` : "";
    const max = r.top.length ? r.top[0].count : 1;
    const rows = r.top.map((t, i) => [i + 1, esc(t.kmer), fmt(t.count), barCell(t.count / max, `${fmt(t.percent, 1)}%`)]);
    return head(`${r.k}-mer counts`, csv) +
      `<p class="summary"><strong>${fmt(r.total)}</strong> overlapping ${r.k}-mers, <strong>${fmt(r.distinct)}</strong> distinct${possible}. The ${plural(r.top.length, "most common one", "most common ones")} are listed below.</p>` +
      table([{ label: "Rank", num: true }, { label: "K-mer", cls: "mono" }, { label: "Count", num: true }, { label: "Share" }], rows);
  },

  composition(r) {
    const rec = state.selected;
    const csv = downloadButton("Windows (CSV)", () => [`${safeName(rec)}_composition.csv`, toCsv([
      ["start", "end", "gc_percent", "gc_skew", "entropy_bits"],
      ...r.windows.map(w => [w.start, w.end, w.gc, w.gc_skew, w.entropy]),
    ]), "text/csv"]);
    const center = w => (w.start - 1 + w.end) / 2;
    const series = key => r.windows.map(w => ({ pos: center(w), value: w[key] }));
    const bands = r.islands.map(i => ({ ...i, label: `CpG island ${fmt(i.start)}–${fmt(i.end)} · ${fmt(i.gc, 1)}% GC · obs/exp ${fmt(i.obs_exp, 2)}` }));
    const skewMax = Math.max(0.1, ...r.windows.map(w => Math.abs(w.gc_skew)));
    const skewLim = Math.ceil(skewMax * 10) / 10;
    const islandRows = r.islands.map((i, n) => [n + 1, `${fmt(i.start)}–${fmt(i.end)}`, fmt(i.length), fmt(i.gc, 1), fmt(i.obs_exp, 2)]);
    return head("Composition", csv) +
      `<div class="tiles">
        ${tile("GC skew", fmt(r.gc_skew, 3), "", "(G − C) / (G + C)")}
        ${tile("AT skew", fmt(r.at_skew, 3), "", "(A − T) / (A + T)")}
        ${tile("Complexity", fmt(r.entropy, 2), "bits", "Shannon entropy, max 2")}
        ${tile("CpG obs/exp", fmt(r.cpg_obs_exp, 2), "", "about 0.2–0.3 in bulk human DNA")}
        ${tile("CpG islands", fmt(r.islands.length), "", "≥ 200 bp, ≥ 50% GC, obs/exp ≥ 0.6")}
      </div>
      <p class="summary">Values in sliding windows of ${fmt(r.window)} bp, every ${plural(r.step, "base")}. ${r.islands.length ? "Shaded regions are CpG islands." : ""}</p>
      ${lineChart(r.length, series("gc"), { label: "GC content (%)", min: 0, max: 100, unit: "%", bands })}
      ${lineChart(r.length, series("gc_skew"), { label: "GC skew", min: -skewLim, max: skewLim, digits: 1, zero: true })}
      ${lineChart(r.length, series("entropy"), { label: "Complexity (bits): dips mark repeats and low-complexity stretches", min: 0, max: 2, digits: 1 })}` +
      (r.islands.length
        ? `<h3 class="section-title">CpG islands</h3>` + table([{ label: "#", num: true }, { label: "Position" },
            { label: "Length (bp)", num: true }, { label: "GC %", num: true }, { label: "Obs/exp CpG", num: true }], islandRows)
        : "") +
      `<p class="note"><strong>GC skew</strong> changes sign near the origin and terminus of replication in many bacteria and phages.
        <strong>CpG islands</strong> use the Gardiner-Garden &amp; Frommer criteria: qualifying 200 bp windows are merged, then trimmed so each island meets the criteria as a whole.
        They are meaningful in vertebrate DNA, where CpG is otherwise rare. Bacteria and phages don't deplete CpG, so long stretches of their genomes qualify.</p>`;
  },

  properties(r) {
    const o = r.longest_orf;
    const orfTile = o
      ? tile("Longest ORF protein", o.mw == null ? "–" : fmt(o.mw / 1000, 2), o.mw == null ? "" : "kDa",
        `${fmt(o.aa)} aa · ORF ${fmt(o.start)}–${fmt(o.end)} (${o.strand})${o.mw == null ? " · contains unknown residues" : ""}`)
      : tile("Longest ORF protein", "–", "", `No ORFs of at least ${fmt(r.min_orf)} nt`);
    const nnTile = r.tm_nn == null
      ? tile("Tm, nearest-neighbor", "–", "", esc(r.tm_nn_error))
      : tile("Tm, nearest-neighbor", fmt(r.tm_nn, 1), "°C", `SantaLucia 1998 · ${fmt(r.na)} mM Na⁺ · ${fmt(r.oligo)} nM primer`);
    return head("Melting temperature & molecular weight") +
      `<div class="tiles wide">
        ${nnTile}
        ${tile("Tm, basic formula", fmt(r.tm, 1), "°C", r.tm_method)}
        ${tile("DNA weight, single strand", fmt(r.mw_single / 1000, 2), "kDa", `${fmt(r.mw_single, 2)} g/mol`)}
        ${tile("DNA weight, double strand", fmt(r.mw_double / 1000, 2), "kDa", `${fmt(r.mw_double, 2)} g/mol`)}
        ${orfTile}
      </div>
      <p class="note">
        <strong>Nearest-neighbor Tm</strong> is the method primer suppliers use. It sums the stacking energy of each
        neighboring base pair (SantaLucia 1998 unified parameters) and corrects for the Na⁺ and primer concentrations set above,
        assuming the primer is in excess over its target. Use it to design and match primers.<br>
        <strong>Basic formula Tm</strong> is a quick estimate that ignores salt and concentration:
        the Wallace rule 2·(A+T) + 4·(G+C) below 14 nt, and 64.9 + 41·(G+C − 16.4)/length from 14 nt up.
        ${r.length > 60 ? "Neither value is meaningful for a sequence this long. Both are meant for oligos up to about 60 nt." : ""}<br>
        <strong>DNA weight</strong> uses anhydrous nucleotide weights (OligoCalc), without a 5′ phosphate.
        <strong>Protein weight</strong> uses average amino acid residue weights (ExPASy).
      </p>`;
  },

  quality(r) {
    if (!r.fastq) {
      return head("Read quality") +
        `<div class="empty-state"><strong>Load a FASTQ file to see read quality</strong>
          FASTQ files from sequencing machines hold a quality score for every base. Try <code>data/samples/lambda_reads.fastq</code>.</div>`;
    }
    const s = r.summary, t = r.trim, a = t.after;
    const lengthText = s.min_length === s.max_length ? `${fmt(s.max_length)}` : `${fmt(s.min_length)}–${fmt(s.max_length)}`;
    const trimmed = downloadButton("Trimmed reads (FASTQ)", () => {
      const res = call("trimmed_fastq", JSON.stringify(params()));
      return ["trimmed_reads.fastq", res.ok ? res.text : "", "text/plain"];
    });
    const csv = downloadButton("Per-position CSV", () => ["quality_by_position.csv", toCsv([
      ["start", "end", "reads", "mean", "p10", "q1", "median", "q3", "p90"],
      ...r.positions.map(p => [p.start, p.end, p.reads, p.mean, p.p10, p.q1, p.median, p.q3, p.p90]),
    ]), "text/csv"]);
    const bars = r.read_quality.map(b => ({
      value: b.reads, label: `${b.q}`, cls: b.q < 20 ? "bad" : b.q < 28 ? "ok" : "good",
      title: `${plural(b.reads, "read")} with mean quality Q${b.q}`,
    }));
    const pctKept = s.bases ? 100 * a.bases / s.bases : 0;
    return head("Read quality", csv) +
      `<div class="tiles">
        ${tile("Reads", fmt(s.reads))}
        ${tile("Read length", lengthText, "bp", s.min_length === s.max_length ? "all reads" : `mean ${fmt(s.mean_length, 1)} bp`)}
        ${tile("Mean quality", `Q${fmt(s.mean_quality, 1)}`, "", "over all bases")}
        ${tile("≥ Q30", fmt(s.q30, 1), "%", "of bases (1 error in 1,000)")}
        ${tile("≥ Q20", fmt(s.q20, 1), "%", "of bases (1 error in 100)")}
        ${tile("GC · N", `${fmt(s.gc, 1)}%`, "", `${fmt(s.n, 2)}% N bases`)}
      </div>
      <p class="summary">A Phred score Q means a 1 in 10<sup>Q/10</sup> chance that the base is wrong.
        Green is good (Q28 and up), amber is acceptable, red is poor (below Q20). Illumina reads usually lose quality toward their 3′ end.${r.bin_size > 1 ? ` Positions are grouped in bins of ${fmt(r.bin_size)} bases.` : ""}</p>
      ${qualityPlot(r.positions, s.max_length)}
      ${histogram(bars, { label: "Reads by mean quality", xLabel: "Mean Phred score of the read" })}
      <h3 class="section-title">Quality trimming</h3>
      <p class="summary">Each read is cut at the first ${t.window}-base window whose mean quality falls below Q${t.threshold}, like Trimmomatic's SLIDINGWINDOW. Reads shorter than ${fmt(t.min_length)} bp after trimming are dropped. Change the settings above.</p>
      <div class="tiles">
        ${tile("Reads kept", fmt(a.reads), "", `${fmt(s.reads ? 100 * a.reads / s.reads : 0, 1)}% of ${fmt(s.reads)}`)}
        ${tile("Bases kept", fmt(pctKept, 1), "%", `${fmt(a.bases)} of ${fmt(s.bases)}`)}
        ${tile("Mean length", fmt(a.mean_length, 1), "bp", "after trimming")}
        ${tile("≥ Q30 after", fmt(a.q30, 1), "%", `was ${fmt(s.q30, 1)}%`)}
        ${tile("Mean quality after", `Q${fmt(a.mean_quality, 1)}`, "", `was Q${fmt(s.mean_quality, 1)}`)}
      </div>
      <div class="actions trim-actions">${trimmed}</div>
      <p class="note">Each read is also loaded as a sequence on the left, so the other tabs work on single reads. The whole set is included in <strong>All sequences</strong>.</p>`;
  },

  report(r) {
    const csv = downloadButton("CSV", () => ["bioseq_report.csv", toCsv([
      ["id", "length_bp", "gc_percent", "A", "C", "G", "T", "orfs", "longest_orf_start", "longest_orf_end", "longest_orf_strand", "longest_orf_aa"],
      ...r.rows.map(x => [x.id, x.length, x.gc, x.counts.A, x.counts.C, x.counts.G, x.counts.T, x.n_orfs,
        x.longest_orf?.start, x.longest_orf?.end, x.longest_orf?.strand, x.longest_orf?.aa]),
    ]), "text/csv"]);
    const rows = r.rows.map(x => [
      `<strong>${esc(x.id)}</strong>`, fmt(x.length), fmt(x.gc, 1),
      fmt(x.counts.A), fmt(x.counts.C), fmt(x.counts.G), fmt(x.counts.T),
      fmt(x.n_orfs),
      x.longest_orf ? `${fmt(x.longest_orf.start)}–${fmt(x.longest_orf.end)} ${strandBadge(x.longest_orf.strand)} · ${fmt(x.longest_orf.aa)} aa` : "–",
    ]);
    return head("All sequences", csv) +
      `<p class="summary">${plural(r.rows.length, "sequence")}, with ORFs of at least ${fmt(r.min_orf)} nt. Click a row to open that sequence.</p>` +
      table([
        { label: "Sequence" }, { label: "Length (bp)", num: true }, { label: "GC %", num: true },
        { label: "A", num: true }, { label: "C", num: true }, { label: "G", num: true }, { label: "T", num: true },
        { label: "ORFs", num: true }, { label: "Longest ORF" },
      ], rows, { rowAttrs: i => `class="clickable" data-select="${esc(r.rows[i].id)}"` });
  },
};

// ------------------------------------------------------------------- events

function bindEvents() {
  $("#btn-analyze").addEventListener("click", analyzeInput);
  $("#btn-example").addEventListener("click", loadExample);
  $("#btn-clear").addEventListener("click", () => {
    $("#seq-input").value = "";
    showInputMsg("", "");
    $("#seq-input").focus();
  });
  $("#seq-input").addEventListener("keydown", e => {
    if (e.key === "Enter" && (e.ctrlKey || e.metaKey)) analyzeInput();
  });
  $("#file-input").addEventListener("change", e => {
    readFile(e.target.files[0]);
    e.target.value = "";
  });

  const panel = $(".input-panel");
  panel.addEventListener("dragover", e => { e.preventDefault(); panel.classList.add("dragging"); });
  panel.addEventListener("dragleave", e => { if (!panel.contains(e.relatedTarget)) panel.classList.remove("dragging"); });
  panel.addEventListener("drop", e => {
    e.preventDefault();
    panel.classList.remove("dragging");
    readFile(e.dataTransfer.files[0]);
  });

  $("#record-list").addEventListener("click", e => {
    const b = e.target.closest("button[data-id]");
    if (!b) return;
    state.selected = b.dataset.id;
    if (DATASET_TABS.includes(state.tab)) state.tab = "overview";
    renderRecords();
    renderTool();
  });

  $(".tabs").addEventListener("click", e => {
    const t = e.target.closest("[data-tab]");
    if (!t) return;
    state.tab = t.dataset.tab;
    renderTool();
  });

  $("#frame-picker").addEventListener("click", e => {
    const b = e.target.closest("[data-frame]");
    if (!b) return;
    state.frame = b.dataset.frame;
    $$("#frame-picker button").forEach(x => x.classList.toggle("active", x === b));
    const all = state.frame === "all";
    $("#to-stop").disabled = all;
    $("#to-stop-wrap").classList.toggle("disabled", all);
    renderTool();
  });
  $("#to-stop").addEventListener("change", e => { state.toStop = e.target.checked; renderTool(); });

  const rerender = debounce(renderTool, 300);
  const bindNumber = (sel, key) => $(sel).addEventListener("input", e => { state[key] = e.target.value; rerender(); });
  bindNumber("#min-orf", "minOrf");
  bindNumber("#k-input", "k");
  bindNumber("#top-input", "top");
  bindNumber("#window-input", "window");
  bindNumber("#na-input", "na");
  bindNumber("#oligo-input", "oligo");
  bindNumber("#trim-q-input", "trimQ");
  bindNumber("#trim-window-input", "trimWindow");
  bindNumber("#min-length-input", "minLength");

  $("#motif-form").addEventListener("submit", e => {
    e.preventDefault();
    state.motif = $("#motif-input").value;
    renderTool();
  });
  $("#motif-chips").addEventListener("click", e => {
    const b = e.target.closest("[data-motif]");
    if (!b) return;
    $("#motif-input").value = b.dataset.motif;
    state.motif = b.dataset.motif;
    renderTool();
  });

  $("#output").addEventListener("click", e => {
    const el = e.target.closest("[data-copy], [data-download], [data-select]");
    if (!el) return;
    if ("copy" in el.dataset) {
      navigator.clipboard.writeText(copies[Number(el.dataset.copy)])
        .then(() => flash(el, "Copied"), () => flash(el, "Copy failed"));
    } else if ("download" in el.dataset) {
      saveFile(...downloads[el.dataset.download]());
    } else if ("select" in el.dataset) {
      state.selected = el.dataset.select;
      state.tab = "overview";
      renderRecords();
      renderTool();
    }
  });

  let lastWidth = window.innerWidth;
  window.addEventListener("resize", debounce(() => {
    if (window.innerWidth !== lastWidth && ["orfs", "motif", "composition", "quality"].includes(state.tab)) renderTool();
    lastWidth = window.innerWidth;
  }, 200));
}

start();
