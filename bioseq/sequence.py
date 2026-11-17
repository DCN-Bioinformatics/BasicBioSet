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

AMINO_ACID_NAMES = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
    "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val", "*": "Stop",
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


def six_frame_translation(seq: str) -> dict[str, str]:
    """Translate all six reading frames.

    Keys "+1", "+2", "+3" start at offsets 0, 1, 2 of the forward strand;
    "-1", "-2", "-3" do the same on the reverse complement.
    """
    seq = validate(seq)
    rc = reverse_complement(seq)
    frames = {f"+{f + 1}": translate(seq[f:]) for f in range(3)}
    frames.update({f"-{f + 1}": translate(rc[f:]) for f in range(3)})
    return frames


def codon_usage(seq: str) -> dict[str, int]:
    """Count the codons of a coding sequence, reading in frame from its first base.

    Returns all 64 codons (unused ones as 0). Codons containing N and a
    trailing partial codon are skipped.
    """
    seq = validate(seq)
    counts = Counter(seq[i:i + 3] for i in range(0, len(seq) - 2, 3))
    return {codon: counts.get(codon, 0) for codon in CODON_TABLE}


def find_orfs(
    seq: str, min_length: int = 30, both_strands: bool = True
) -> list[tuple[int, int, str, str]]:
    """Find open reading frames (ATG ... stop) in all three reading frames.

    Searches the forward strand and, if both_strands is True, the reverse
    complement too. Returns (start, end, strand, protein) tuples sorted by
    position, where strand is "+" or "-". Coordinates are always 0-based on
    the forward strand with end exclusive, so seq[start:end] covers the ORF
    (for "-" ORFs, reverse_complement(seq[start:end]) gives the coding sequence).
    min_length is the minimum ORF length in nucleotides, including the stop codon.
    """
    seq = validate(seq)
    orfs = [(start, end, "+", protein)
            for start, end, protein in _scan_orfs(seq, min_length)]
    if both_strands:
        n = len(seq)
        orfs += [(n - end, n - start, "-", protein)
                 for start, end, protein in _scan_orfs(reverse_complement(seq), min_length)]
    return sorted(orfs)


def orf_codon_usage(
    seq: str, min_length: int = 30, both_strands: bool = True
) -> dict[str, int]:
    """Codon usage summed over every ORF find_orfs reports (stop codons included).

    "-" strand ORFs are read from their reverse complement, so every codon is
    counted in its coding orientation.
    """
    seq = validate(seq)
    usage = dict.fromkeys(CODON_TABLE, 0)
    for start, end, strand, _ in find_orfs(seq, min_length, both_strands):
        coding = seq[start:end] if strand == "+" else reverse_complement(seq[start:end])
        for codon, n in codon_usage(coding).items():
            usage[codon] += n
    return usage


def _scan_orfs(seq: str, min_length: int) -> list[tuple[int, int, str]]:
    """Find ORFs on one strand of an already-validated sequence.

    Returns (start, end, protein) tuples relative to the given strand.
    """
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
    return orfs
