"""Composition along a sequence: skew, sliding windows, complexity, CpG islands."""

import math
from collections import Counter

from .sequence import validate

# Gardiner-Garden & Frommer (1987) CpG island criteria.
CPG_MIN_LENGTH = 200
CPG_MIN_GC = 50.0
CPG_MIN_OBS_EXP = 0.6


def _skew(a: int, b: int) -> float:
    return (a - b) / (a + b) if a + b else 0.0


def gc_skew(seq: str) -> float:
    """(G - C) / (G + C), from -1 to 1. 0 when the sequence has no G or C."""
    seq = validate(seq)
    return _skew(seq.count("G"), seq.count("C"))


def at_skew(seq: str) -> float:
    """(A - T) / (A + T), from -1 to 1. 0 when the sequence has no A or T."""
    seq = validate(seq)
    return _skew(seq.count("A"), seq.count("T"))


def shannon_entropy(seq: str) -> float:
    """Shannon entropy of the base composition in bits (N ignored).

    2.0 means all four bases are equally common; low values flag
    low-complexity stretches such as poly-A tails or CACACA repeats.
    """
    counts = Counter(b for b in validate(seq) if b != "N")
    total = sum(counts.values())
    return -sum(n / total * math.log2(n / total) for n in counts.values()) if total else 0.0


def cpg_obs_exp(seq: str) -> float:
    """Observed / expected CpG ratio: CG count * length / (C count * G count)."""
    seq = validate(seq)
    c, g = seq.count("C"), seq.count("G")
    return seq.count("CG") * len(seq) / (c * g) if c and g else 0.0


def sliding_window(seq: str, window: int = 100, step: int | None = None) -> list[dict]:
    """GC content, GC skew and entropy in windows along the sequence.

    Returns one dict per window: start (0-based), end (exclusive), gc (%),
    gc_skew, at_skew and entropy. A window larger than the sequence shrinks
    to the whole sequence. step defaults to window // 2.
    """
    seq = validate(seq)
    if window < 1 or (step is not None and step < 1):
        raise ValueError("Window and step must be at least 1")
    window = min(window, len(seq))
    if not window:
        return []
    step = step or max(1, window // 2)
    rows = []
    for start in range(0, len(seq) - window + 1, step):
        chunk = seq[start:start + window]
        rows.append({
            "start": start,
            "end": start + window,
            "gc": 100 * (chunk.count("G") + chunk.count("C")) / window,
            "gc_skew": _skew(chunk.count("G"), chunk.count("C")),
            "at_skew": _skew(chunk.count("A"), chunk.count("T")),
            "entropy": shannon_entropy(chunk),
        })
    return rows


def _is_cpg_island(length: int, c: int, g: int, cg: int) -> bool:
    return (c > 0 and g > 0 and 100 * (c + g) / length >= CPG_MIN_GC
            and cg * length / (c * g) >= CPG_MIN_OBS_EXP)


def _trim_island(seq: str, start: int, end: int, min_length: int) -> tuple[int, int] | None:
    """Shrink a merged region from its ends until it meets the island criteria again.

    Merging overlapping windows can leave a region whose overall ratio dips
    just below a threshold. Each step drops the end base that helps least
    (A/T before G/C, and never half of a CpG before a lone base).
    """
    region = seq[start:end]
    c, g, cg = region.count("C"), region.count("G"), region.count("CG")
    while end - start >= min_length:
        if _is_cpg_island(end - start, c, g, cg):
            return start, end
        left_cost = (seq[start] in "GC") + (seq[start:start + 2] == "CG")
        right_cost = (seq[end - 1] in "GC") + (seq[end - 2:end] == "CG")
        if left_cost <= right_cost:
            base, pair = seq[start], seq[start:start + 2]
            start += 1
        else:
            base, pair = seq[end - 1], seq[end - 2:end]
            end -= 1
        c -= base == "C"
        g -= base == "G"
        cg -= pair == "CG"
    return None


def find_cpg_islands(seq: str, min_length: int = CPG_MIN_LENGTH) -> list[tuple[int, int]]:
    """Find CpG islands with the Gardiner-Garden & Frommer criteria.

    Every window of min_length bp with GC >= 50% and observed/expected
    CpG >= 0.6 is marked, overlapping marked windows are merged, and each
    merged region is trimmed at its ends until it meets the criteria as a
    whole. Returns (start, end) tuples, 0-based with end exclusive.
    """
    seq = validate(seq)
    w = min_length
    if w < 2 or len(seq) < w:
        return []
    # Running counts for the window seq[i:i + w], updated one base at a time.
    c, g = seq[:w].count("C"), seq[:w].count("G")
    cg = seq[:w].count("CG")
    merged: list[tuple[int, int]] = []
    for i in range(len(seq) - w + 1):
        if i:
            out, new = seq[i - 1], seq[i + w - 1]
            c += (new == "C") - (out == "C")
            g += (new == "G") - (out == "G")
            cg += (seq[i + w - 2:i + w] == "CG") - (seq[i - 1:i + 1] == "CG")
        if _is_cpg_island(w, c, g, cg):
            if merged and i <= merged[-1][1]:
                merged[-1] = (merged[-1][0], i + w)
            else:
                merged.append((i, i + w))
    trimmed = (_trim_island(seq, start, end, w) for start, end in merged)
    return [island for island in trimmed if island]
