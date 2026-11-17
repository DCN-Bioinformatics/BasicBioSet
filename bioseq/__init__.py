"""bioseq: a small toolkit for basic DNA sequence analysis."""

from .fasta import read_fasta
from .motifs import COMMON_MOTIFS, find_motif, kmer_counts, resolve_motif
from .properties import (
    dna_molecular_weight,
    melting_temperature,
    protein_molecular_weight,
)
from .sequence import (
    AMINO_ACID_NAMES,
    CODON_TABLE,
    base_counts,
    codon_usage,
    find_orfs,
    gc_content,
    orf_codon_usage,
    reverse_complement,
    six_frame_translation,
    transcribe,
    translate,
    validate,
)

__all__ = [
    "read_fasta",
    "COMMON_MOTIFS",
    "find_motif",
    "kmer_counts",
    "resolve_motif",
    "dna_molecular_weight",
    "melting_temperature",
    "protein_molecular_weight",
    "AMINO_ACID_NAMES",
    "CODON_TABLE",
    "base_counts",
    "codon_usage",
    "find_orfs",
    "gc_content",
    "orf_codon_usage",
    "reverse_complement",
    "six_frame_translation",
    "transcribe",
    "translate",
    "validate",
]
