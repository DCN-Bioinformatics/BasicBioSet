"""JSON bridge between the web page and the bioseq package.

The page runs this file inside the browser with Pyodide. Every public
function takes plain strings and returns a JSON string, so the JavaScript
side never handles Python objects. It also runs in normal Python, which is
how tests/test_web.py checks it.
"""

import json
import os
import tempfile

from bioseq import (
    AMINO_ACID_NAMES,
    CODON_TABLE,
    COMMON_MOTIFS,
    RESTRICTION_ENZYMES,
    at_skew,
    base_counts,
    cpg_obs_exp,
    cut_position,
    dna_molecular_weight,
    fastq_summary,
    find_cpg_islands,
    find_enzyme,
    find_motif,
    find_orfs,
    format_fastq,
    gc_content,
    gc_skew,
    is_fastq,
    is_rna,
    kmer_counts,
    melting_temperature,
    melting_temperature_nn,
    orf_codon_usage,
    overhang,
    parse_fastq,
    protein_molecular_weight,
    quality_by_position,
    read_fasta,
    resolve_motif,
    reverse_complement,
    shannon_entropy,
    six_frame_translation,
    sliding_window,
    transcribe,
    translate,
    trim_read,
    validate,
)

MAX_HITS = 2000  # cap on motif hits sent to the page
MAX_WINDOWS = 800  # cap on sliding-window points sent to the page
MAX_READS = 20_000  # FASTQ reads loaded (keeps the browser responsive)
MAX_QUALITY_BINS = 150  # cap on per-position quality points

_records: dict[str, str] = {}
_qualities: dict[str, list[int]] = {}  # Phred scores per read, for FASTQ input


def _ok(**data) -> str:
    return json.dumps({"ok": True, **data})


def _error(message: str) -> str:
    return json.dumps({"ok": False, "error": message})


def _positive_number(params: dict, key: str, default: float, label: str) -> float:
    try:
        value = float(params.get(key, default))
    except (TypeError, ValueError):
        value = 0.0
    if not value > 0:
        raise ValueError(f"{label} must be a number above 0.")
    return value


def _positive_int(params: dict, key: str, default: int, label: str) -> int:
    try:
        value = int(params.get(key, default))
    except (TypeError, ValueError):
        value = 0
    if value < 1:
        raise ValueError(f"{label} must be a whole number of at least 1.")
    return value


def _protein_mw(protein: str) -> float | None:
    """Protein weight, or None if it has residues we can't weigh (X from N codons)."""
    try:
        return round(protein_molecular_weight(protein), 2)
    except ValueError:
        return None


def _orf_list(seq: str, min_orf: int) -> list[dict]:
    return [
        {"start": start + 1, "end": end, "strand": strand, "length": end - start,
         "protein": protein, "aa": len(protein), "mw": _protein_mw(protein)}
        for start, end, strand, protein in find_orfs(seq, min_length=min_orf)
    ]


# --------------------------------------------------------------------- input

def info() -> str:
    """Static data the page needs once at start-up."""
    return _ok(motifs=COMMON_MOTIFS, enzymes=RESTRICTION_ENZYMES)


def load_text(text: str) -> str:
    """Parse pasted or uploaded text.

    Text containing '>' is read as FASTA (with bioseq's own reader). Anything
    else is treated as one bare sequence; spaces, line breaks and position
    numbers (as in GenBank listings) are ignored. RNA is accepted and read as
    DNA (U -> T); such records are flagged with "rna": true. Text whose first
    line starts with '@' is read as FASTQ (at most MAX_READS reads); every read
    becomes a record and its quality scores are kept for the 'quality' tool.
    """
    global _records, _qualities
    text = text.strip()
    if not text:
        return _error("Paste a DNA sequence or a FASTA file first.")

    if is_fastq(text):
        return _load_fastq(text)
    if ">" in text:
        with tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False) as handle:
            handle.write(text)
            path = handle.name
        try:
            records = read_fasta(path)
        except ValueError:
            return _error("Found sequence text before the first '>' header line. "
                          "In FASTA format every sequence starts with a '>name' line.")
        finally:
            os.remove(path)
    else:
        records = {"sequence": "".join(c for c in text if not (c.isspace() or c.isdigit()))}

    if not records:
        return _error("No sequences found.")
    checked = {}
    rna = set()
    for seq_id, seq in records.items():
        try:
            checked[seq_id] = validate(seq)
        except ValueError as e:
            return _error(f"{seq_id}: {e}. Only A, C, G, T (or U) and N are allowed.")
        if is_rna(seq):
            rna.add(seq_id)
    _records, _qualities = checked, {}
    return _ok(format="fasta", records=_record_list(rna))


def _record_list(rna: set) -> list[dict]:
    return [{"id": seq_id, "length": len(seq), "gc": round(gc_content(seq), 1), "rna": seq_id in rna}
            for seq_id, seq in _records.items()]


def _load_fastq(text: str) -> str:
    global _records, _qualities
    try:
        reads = parse_fastq(text.splitlines(), max_reads=MAX_READS + 1)
    except ValueError as e:
        return _error(f"This looks like FASTQ, but it could not be read. {e}.")
    truncated = len(reads) > MAX_READS
    reads = reads[:MAX_READS]
    if not reads:
        return _error("No reads found.")
    records, qualities, rna = {}, {}, set()
    for read_id, seq, scores in reads:
        unique, n = read_id, 1
        while unique in records:  # keep reads that share an id
            n += 1
            unique = f"{read_id}_{n}"
        try:
            records[unique] = validate(seq)
        except ValueError as e:
            return _error(f"{read_id}: {e}. Only A, C, G, T (or U) and N are allowed.")
        if is_rna(seq):
            rna.add(unique)
        qualities[unique] = scores
    _records, _qualities = records, qualities
    return _ok(format="fastq", truncated=truncated, max_reads=MAX_READS, records=_record_list(rna))


# --------------------------------------------------------------------- tools

def _basics(seq: str, params: dict) -> dict:
    return {
        "length": len(seq),
        "counts": base_counts(seq),
        "n": seq.count("N"),
        "gc": round(gc_content(seq), 2),
        "revcomp": reverse_complement(seq),
        "mrna": transcribe(seq),
    }


def _translation(seq: str, params: dict) -> dict:
    frame = str(params.get("frame", "1"))
    if frame == "all":
        return {"frames": six_frame_translation(seq)}
    if frame not in ("1", "2", "3"):
        raise ValueError("Reading frame must be 1, 2, 3 or all.")
    return {"frame": frame,
            "protein": translate(seq[int(frame) - 1:], to_stop=bool(params.get("to_stop")))}


def _orfs(seq: str, params: dict) -> dict:
    min_orf = _positive_int(params, "min_orf", 30, "Minimum ORF length")
    return {"min_orf": min_orf, "length": len(seq), "orfs": _orf_list(seq, min_orf)}


def _codons(seq: str, params: dict) -> dict:
    min_orf = _positive_int(params, "min_orf", 30, "Minimum ORF length")
    usage = orf_codon_usage(seq, min_length=min_orf)
    synonyms: dict[str, list[str]] = {}
    for codon, aa in CODON_TABLE.items():
        synonyms.setdefault(aa, []).append(codon)
    groups = []
    for aa in sorted(synonyms, key=lambda a: AMINO_ACID_NAMES[a]):
        total = sum(usage[c] for c in synonyms[aa])
        groups.append({
            "aa": aa, "name": AMINO_ACID_NAMES[aa], "total": total,
            "codons": [{"codon": c, "count": usage[c],
                        "fraction": usage[c] / total if total else None}
                       for c in synonyms[aa]],
        })
    return {"min_orf": min_orf, "n_orfs": len(find_orfs(seq, min_length=min_orf)),
            "total": sum(usage.values()), "groups": groups}


def _motif(seq: str, params: dict) -> dict:
    query = str(params.get("motif", "")).strip()
    motif = resolve_motif(query)
    enzyme = find_enzyme(query)
    hits = find_motif(seq, motif)
    rows = []
    for start, end, strand in hits[:MAX_HITS]:
        row = {"start": start + 1, "end": end, "strand": strand,
               "site": seq[start:end] if strand == "+" else reverse_complement(seq[start:end])}
        if enzyme:
            row["cut"] = cut_position(start, end, strand, enzyme[1])
        rows.append(row)
    result = {"query": query, "motif": motif, "length": len(seq),
              "total": len(hits), "hits": rows}
    if enzyme:
        result["enzyme"] = {"name": enzyme[0], "site": enzyme[1], "ends": overhang(enzyme[1])}
    return result


def _kmers(seq: str, params: dict) -> dict:
    k = _positive_int(params, "k", 3, "k")
    top = _positive_int(params, "top", 10, "Number of k-mers to show")
    counts = kmer_counts(seq, k)
    total = sum(counts.values())
    return {"k": k, "total": total, "distinct": len(counts),
            "top": [{"kmer": kmer, "count": n, "percent": 100 * n / total}
                    for kmer, n in list(counts.items())[:top]]}


def _properties(seq: str, params: dict) -> dict:
    min_orf = _positive_int(params, "min_orf", 30, "Minimum ORF length")
    na = _positive_number(params, "na", 50, "Na+ concentration")
    oligo = _positive_number(params, "oligo", 250, "Primer concentration")
    orfs = _orf_list(seq, min_orf)
    longest = max(orfs, key=lambda o: o["length"]) if orfs else None
    try:
        tm_nn, tm_nn_error = round(melting_temperature_nn(seq, na, oligo), 1), None
    except ValueError as e:
        tm_nn, tm_nn_error = None, str(e)
    return {
        "length": len(seq),
        "tm_nn": tm_nn,
        "tm_nn_error": tm_nn_error,
        "na": na,
        "oligo": oligo,
        "tm": round(melting_temperature(seq), 1),
        "tm_method": "Wallace rule" if len(seq) < 14 else "GC formula",
        "mw_single": round(dna_molecular_weight(seq), 2),
        "mw_double": round(dna_molecular_weight(seq, double_stranded=True), 2),
        "min_orf": min_orf,
        "longest_orf": longest,
    }


def _composition(seq: str, params: dict) -> dict:
    window = _positive_int(params, "window", 100, "Window size")
    window = min(window, len(seq))
    step = max(1, window // 4, (len(seq) - window) // MAX_WINDOWS + 1)
    rows = sliding_window(seq, window, step) if seq else []
    islands = [{"start": start + 1, "end": end, "length": end - start,
                "gc": round(gc_content(seq[start:end]), 1),
                "obs_exp": round(cpg_obs_exp(seq[start:end]), 2)}
               for start, end in find_cpg_islands(seq)]
    return {
        "length": len(seq),
        "window": window,
        "step": step,
        "gc": round(gc_content(seq), 2),
        "gc_skew": round(gc_skew(seq), 4),
        "at_skew": round(at_skew(seq), 4),
        "entropy": round(shannon_entropy(seq), 3),
        "cpg_obs_exp": round(cpg_obs_exp(seq), 2),
        "windows": [{"start": r["start"] + 1, "end": r["end"], "gc": round(r["gc"], 2),
                     "gc_skew": round(r["gc_skew"], 4), "entropy": round(r["entropy"], 3)}
                    for r in rows],
        "islands": islands,
    }


def _trim_params(params: dict) -> tuple[int, int, int]:
    return (_positive_int(params, "trim_q", 20, "Quality threshold"),
            _positive_int(params, "trim_window", 4, "Trimming window"),
            _positive_int(params, "min_length", 36, "Minimum read length"))


def _trimmed_reads(params: dict) -> list[tuple[str, str, list[int]]]:
    threshold, window, min_length = _trim_params(params)
    kept = []
    for read_id, seq in _records.items():
        trimmed, scores = trim_read(seq, _qualities[read_id], threshold, window)
        if len(trimmed) >= min_length:
            kept.append((read_id, trimmed, scores))
    return kept


def _reads() -> list[tuple[str, str, list[int]]]:
    return [(read_id, seq, _qualities[read_id]) for read_id, seq in _records.items()]


def _rounded(summary: dict) -> dict:
    return {k: round(v, 2) if isinstance(v, float) else v for k, v in summary.items()}


def _quality(params: dict) -> dict:
    if not _qualities:
        return {"fastq": False}
    threshold, window, min_length = _trim_params(params)
    reads = _reads()
    before = fastq_summary(reads)
    bin_size = -(-before["max_length"] // MAX_QUALITY_BINS)  # ceiling division
    positions = [{**row, "mean": round(row["mean"], 2)}
                 for row in quality_by_position([q for _, _, q in reads], bin_size)]
    read_means = [round(sum(q) / len(q)) for _, _, q in reads if q]
    histogram = [{"q": q, "reads": read_means.count(q)}
                 for q in range(min(read_means, default=0), max(read_means, default=0) + 1)]
    lengths = sorted({len(seq) for _, seq, _ in reads})
    after = fastq_summary(_trimmed_reads(params))
    return {
        "fastq": True,
        "summary": _rounded(before),
        "bin_size": bin_size,
        "positions": positions,
        "read_quality": histogram,
        "variable_length": len(lengths) > 1,
        "trim": {"threshold": threshold, "window": window, "min_length": min_length,
                 "after": _rounded(after)},
    }


def trimmed_fastq(params_json: str = "{}") -> str:
    """The loaded FASTQ reads after quality trimming, as FASTQ text."""
    try:
        if not _qualities:
            return _error("Load a FASTQ file first.")
        kept = _trimmed_reads(json.loads(params_json or "{}"))
        return _ok(text=format_fastq(kept), reads=len(kept))
    except ValueError as e:
        return _error(str(e))


def _report(params: dict) -> dict:
    min_orf = _positive_int(params, "min_orf", 30, "Minimum ORF length")
    rows = []
    for seq_id, seq in _records.items():
        orfs = _orf_list(seq, min_orf)
        longest = max(orfs, key=lambda o: o["length"]) if orfs else None
        rows.append({"id": seq_id, "length": len(seq), "gc": round(gc_content(seq), 1),
                     "counts": base_counts(seq), "n_orfs": len(orfs), "longest_orf": longest})
    return {"min_orf": min_orf, "rows": rows}


TOOLS = {
    "basics": _basics,
    "translation": _translation,
    "orfs": _orfs,
    "codons": _codons,
    "motif": _motif,
    "kmers": _kmers,
    "properties": _properties,
    "composition": _composition,
}


def run(tool: str, seq_id: str, params_json: str = "{}") -> str:
    """Run one tool on a loaded sequence (or on all of them for 'report')."""
    try:
        params = json.loads(params_json or "{}")
        if tool == "report":
            return _ok(**_report(params))
        if tool == "quality":
            return _ok(**_quality(params))
        if tool not in TOOLS:
            return _error(f"Unknown tool: {tool}")
        if seq_id not in _records:
            return _error("Load a sequence first.")
        return _ok(**TOOLS[tool](_records[seq_id], params))
    except ValueError as e:
        return _error(str(e))
