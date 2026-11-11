"""Basic DNA sequence operations."""

from collections import Counter

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

VALID_BASES = set("ACGTN")


def validate(seq: str) -> str:
    """Return the sequence in uppercase, or raise ValueError if it has invalid bases."""
    seq = seq.upper()
    invalid = set(seq) - VALID_BASES
    if invalid:
        raise ValueError(f"Invalid DNA bases: {', '.join(sorted(invalid))}")
    return seq


def base_counts(seq: str) -> dict[str, int]:
    """Count each nucleotide (A, C, G, T)."""
    counts = Counter(validate(seq))
    return {base: counts.get(base, 0) for base in "ACGT"}


def gc_content(seq: str) -> float:
    """Percentage of G and C bases in the sequence."""
    seq = validate(seq)
    if not seq:
        return 0.0
    return 100 * (seq.count("G") + seq.count("C")) / len(seq)


def reverse_complement(seq: str) -> str:
    """Reverse complement of a DNA strand."""
    return validate(seq).translate(COMPLEMENT)[::-1]


def transcribe(seq: str) -> str:
    """DNA -> mRNA (T becomes U)."""
    return validate(seq).replace("T", "U")


def translate(seq: str, to_stop: bool = False) -> str:
    """DNA -> protein using the standard genetic code.

    Unknown codons (e.g. containing N) become 'X'. Stop codons are '*'.
    If to_stop is True, translation ends at the first stop codon.
    """
    seq = validate(seq)
    protein = []
    for i in range(0, len(seq) - 2, 3):
        aa = CODON_TABLE.get(seq[i:i + 3], "X")
        if aa == "*" and to_stop:
            break
        protein.append(aa)
    return "".join(protein)


def find_orfs(seq: str, min_length: int = 30) -> list[tuple[int, int, str]]:
    """Find open reading frames (ATG ... stop) on the forward strand.

    Returns (start, end, protein) tuples, 0-based with end exclusive.
    min_length is the minimum ORF length in nucleotides, including the stop codon.
    """
    seq = validate(seq)
    orfs = []
    for frame in range(3):
        i = frame
        while i < len(seq) - 2:
            if seq[i:i + 3] == "ATG":
                for j in range(i, len(seq) - 2, 3):
                    if CODON_TABLE.get(seq[j:j + 3]) == "*":
                        if j + 3 - i >= min_length:
                            orfs.append((i, j + 3, translate(seq[i:j])))
                        i = j  # continue scanning after this ORF
                        break
            i += 3
    return sorted(orfs)
