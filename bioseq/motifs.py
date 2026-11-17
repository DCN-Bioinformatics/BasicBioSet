"""Motif search and k-mer counting."""

import re
from collections import Counter

from .sequence import reverse_complement, validate

# IUPAC nucleotide codes -> the bases each one matches.
# A motif N also matches an N in the sequence; concrete bases never do.
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
    "N": "ACGTN",
}

# A few well-known sites, usable by name in the menu.
COMMON_MOTIFS = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "NotI": "GCGGCCGC",
    "TATA-box": "TATAWAW",
}


def resolve_motif(name_or_motif: str) -> str:
    """Return the motif for a COMMON_MOTIFS name (any case), else the input uppercased."""
    by_name = {name.lower(): motif for name, motif in COMMON_MOTIFS.items()}
    text = name_or_motif.strip()
    return by_name.get(text.lower(), text.upper())


def _motif_regex(motif: str) -> re.Pattern:
    """Compile an IUPAC motif into a regex that also finds overlapping hits."""
    motif = motif.upper()
    if not motif:
        raise ValueError("Motif is empty")
    invalid = set(motif) - set(IUPAC)
    if invalid:
        raise ValueError(f"Invalid motif characters: {', '.join(sorted(invalid))}")
    pattern = "".join(c if len(IUPAC[c]) == 1 else f"[{IUPAC[c]}]" for c in motif)
    return re.compile(f"(?=({pattern}))")  # lookahead so matches can overlap


def find_motif(
    seq: str, motif: str, both_strands: bool = True
) -> list[tuple[int, int, str]]:
    """Find every occurrence of a motif (IUPAC codes allowed).

    Returns (start, end, strand) tuples sorted by position. As with find_orfs,
    coordinates are 0-based on the forward strand with end exclusive; for "-"
    hits, reverse_complement(seq[start:end]) is the matching sequence.
    Palindromic sites such as GAATTC are reported once per strand.
    """
    seq = validate(seq)
    regex = _motif_regex(motif)
    size = len(motif)
    hits = [(m.start(), m.start() + size, "+") for m in regex.finditer(seq)]
    if both_strands:
        n = len(seq)
        hits += [(n - m.start() - size, n - m.start(), "-")
                 for m in regex.finditer(reverse_complement(seq))]
    return sorted(hits)


def kmer_counts(seq: str, k: int) -> dict[str, int]:
    """Count every overlapping substring of length k.

    Returns {kmer: count}, most frequent first (ties in alphabetical order).
    A sequence shorter than k gives an empty dict.
    """
    if k < 1:
        raise ValueError("k must be at least 1")
    seq = validate(seq)
    counts = Counter(seq[i:i + k] for i in range(len(seq) - k + 1))
    return dict(sorted(counts.items(), key=lambda item: (-item[1], item[0])))
