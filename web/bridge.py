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
    base_counts,
    dna_molecular_weight,
    find_motif,
    find_orfs,
    gc_content,
    kmer_counts,
    melting_temperature,
    orf_codon_usage,
    protein_molecular_weight,
    read_fasta,
    resolve_motif,
    reverse_complement,
    six_frame_translation,
    transcribe,
    translate,
    validate,
)

MAX_HITS = 2000  # cap on motif hits sent to the page

_records: dict[str, str] = {}


def _ok(**data) -> str:
    return json.dumps({"ok": True, **data})


def _error(message: str) -> str:
    return json.dumps({"ok": False, "error": message})


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
    return _ok(motifs=COMMON_MOTIFS)


def load_text(text: str) -> str:
    """Parse pasted or uploaded text.

    Text containing '>' is read as FASTA (with bioseq's own reader). Anything
    else is treated as one bare sequence; spaces, line breaks and position
    numbers (as in GenBank listings) are ignored.
    """
    global _records
    text = text.strip()
    if not text:
        return _error("Paste a DNA sequence or a FASTA file first.")

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
    for seq_id, seq in records.items():
        try:
            checked[seq_id] = validate(seq)
        except ValueError as e:
            hint = " (RNA? replace U with T)" if "U" in seq.upper() else ""
            return _error(f"{seq_id}: {e}{hint}. Only A, C, G, T and N are allowed.")
    _records = checked
    return _ok(records=[
        {"id": seq_id, "length": len(seq), "gc": round(gc_content(seq), 1)}
        for seq_id, seq in _records.items()
    ])


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
    hits = find_motif(seq, motif)
    rows = [{"start": start + 1, "end": end, "strand": strand,
             "site": seq[start:end] if strand == "+" else reverse_complement(seq[start:end])}
            for start, end, strand in hits[:MAX_HITS]]
    return {"query": query, "motif": motif, "length": len(seq),
            "total": len(hits), "hits": rows}


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
    orfs = _orf_list(seq, min_orf)
    longest = max(orfs, key=lambda o: o["length"]) if orfs else None
    return {
        "length": len(seq),
        "tm": round(melting_temperature(seq), 1),
        "tm_method": "Wallace rule" if len(seq) < 14 else "GC formula",
        "mw_single": round(dna_molecular_weight(seq), 2),
        "mw_double": round(dna_molecular_weight(seq, double_stranded=True), 2),
        "min_orf": min_orf,
        "longest_orf": longest,
    }


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
}


def run(tool: str, seq_id: str, params_json: str = "{}") -> str:
    """Run one tool on a loaded sequence (or on all of them for 'report')."""
    try:
        params = json.loads(params_json or "{}")
        if tool == "report":
            return _ok(**_report(params))
        if tool not in TOOLS:
            return _error(f"Unknown tool: {tool}")
        if seq_id not in _records:
            return _error("Load a sequence first.")
        return _ok(**TOOLS[tool](_records[seq_id], params))
    except ValueError as e:
        return _error(str(e))
