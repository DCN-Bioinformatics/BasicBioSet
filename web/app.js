"use strict";

// The page runs the real bioseq Python package in the browser with Pyodide.
// web/bridge.py turns each analysis into JSON; this file only draws it.

const PYODIDE_URL = "https://cdn.jsdelivr.net/pyodide/v0.26.4/full/";
const BIOSEQ_FILES = ["__init__.py", "fasta.py", "sequence.py", "motifs.py", "properties.py"];
const SEQ_DISPLAY_LIMIT = 12000;   // longer sequences are shown truncated (copy/download stay full)
const TABLE_ROW_LIMIT = 500;       // same idea for long tables

const TOOL_FOR_TAB = {
  overview: "basics", translation: "translation", orfs: "orfs", codons: "codons",
  motif: "motif", kmers: "kmers", properties: "properties", report: "report",
};

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
  state.records = result.records;
  if (!state.records.some(r => r.id === state.selected)) state.selected = state.records[0].id;
  showInputMsg("ok", `Loaded ${plural(state.records.length, "sequence")}.`);
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
  list.innerHTML = state.records.map(r => `
    <li><button data-id="${esc(r.id)}" class="${r.id === state.selected ? "active" : ""}">
      <span class="rec-id">${esc(r.id)}</span>
      <span class="rec-meta">${fmt(r.length)} bp · ${fmt(r.gc, 1)}% GC</span>
    </button></li>`).join("");
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
  };
}

function updateTabs() {
  if (location.hash.slice(1) !== state.tab) history.replaceState(null, "", `#${state.tab}`);
  $$(".tabs [data-tab]").forEach(t => t.setAttribute("aria-selected", String(t.dataset.tab === state.tab)));
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
  const sub = state.tab === "report" ? "" :
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
      ["start", "end", "strand", "site"], ...r.hits.map(h => [h.start, h.end, h.strand, h.site]),
    ]), "text/csv"]);
    const features = r.hits.map(h => ({ ...h, label: `${h.site} at ${h.start}–${h.end} (${h.strand})` }));
    const rows = r.hits.map((h, i) => [i + 1, `${fmt(h.start)}–${fmt(h.end)}`, strandBadge(h.strand), esc(h.site)]);
    const capped = r.total > r.hits.length
      ? ` Only the first ${fmt(r.hits.length)} are listed.` : "";
    return head("Motif search", csv) +
      `<p class="summary"><strong>${plural(r.total, "match", "matches")}</strong> for <strong>${esc(label)}</strong>: ${fmt(plusHits)} on the + strand and ${fmt(r.hits.length - plusHits)} on the − strand.${capped} Sites that read the same on both strands, like GAATTC, appear once per strand.</p>` +
      featureMap(r.length, features, { ticks: true }) +
      table([{ label: "#", num: true }, { label: "Position" }, { label: "Strand" }, { label: "Site (5′→3′)", cls: "mono" }], rows);
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

  properties(r) {
    const o = r.longest_orf;
    const orfTile = o
      ? tile("Longest ORF protein", o.mw == null ? "–" : fmt(o.mw / 1000, 2), o.mw == null ? "" : "kDa",
        `${fmt(o.aa)} aa · ORF ${fmt(o.start)}–${fmt(o.end)} (${o.strand})${o.mw == null ? " · contains unknown residues" : ""}`)
      : tile("Longest ORF protein", "–", "", `No ORFs of at least ${fmt(r.min_orf)} nt`);
    return head("Melting temperature & molecular weight") +
      `<div class="tiles wide">
        ${tile("Melting temperature", fmt(r.tm, 1), "°C", r.tm_method)}
        ${tile("DNA weight, single strand", fmt(r.mw_single / 1000, 2), "kDa", `${fmt(r.mw_single, 2)} g/mol`)}
        ${tile("DNA weight, double strand", fmt(r.mw_double / 1000, 2), "kDa", `${fmt(r.mw_double, 2)} g/mol`)}
        ${orfTile}
      </div>
      <p class="note">
        <strong>Melting temperature</strong> is a quick estimate meant for short oligos such as PCR primers.
        Below 14 nt it uses the Wallace rule, 2·(A+T) + 4·(G+C). From 14 nt up it uses 64.9 + 41·(G+C − 16.4)/length.
        Neither formula accounts for salt or primer concentration${r.length > 60 ? ", and the value isn't meaningful for a sequence this long" : ""}.<br>
        <strong>DNA weight</strong> uses anhydrous nucleotide weights (OligoCalc), without a 5′ phosphate.
        <strong>Protein weight</strong> uses average amino acid residue weights (ExPASy).
      </p>`;
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
    if (state.tab === "report") state.tab = "overview";
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
    if (window.innerWidth !== lastWidth && (state.tab === "orfs" || state.tab === "motif")) renderTool();
    lastWidth = window.innerWidth;
  }, 200));
}

start();
