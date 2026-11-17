"""Physical properties: melting temperature and molecular weight."""

from .sequence import reverse_complement, validate

# Anhydrous nucleotide weights in g/mol (as used by OligoCalc).
# N uses the average of the four bases.
DNA_WEIGHTS = {"A": 313.21, "C": 289.18, "G": 329.21, "T": 304.2, "N": 308.95}

# Average amino acid residue weights in g/mol (as used by ExPASy ProtParam).
PROTEIN_WEIGHTS = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388,
    "E": 129.1155, "Q": 128.1307, "G": 57.0519, "H": 137.1411, "I": 113.1594,
    "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766, "P": 97.1167,
    "S": 87.0782, "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}
WATER = 18.01524


def melting_temperature(seq: str) -> float:
    """Rough melting temperature (°C) of a short oligo such as a PCR primer.

    Fewer than 14 nt: Wallace rule, 2*(A+T) + 4*(G+C).
    14 nt or more:    64.9 + 41*(G+C - 16.4) / length.
    Ns count toward the length but not toward A/T or G/C. These are quick
    estimates that ignore salt and primer concentration.
    """
    seq = validate(seq)
    if not seq:
        raise ValueError("Cannot compute the melting temperature of an empty sequence")
    at = seq.count("A") + seq.count("T")
    gc = seq.count("G") + seq.count("C")
    if len(seq) < 14:
        return float(2 * at + 4 * gc)
    return 64.9 + 41 * (gc - 16.4) / len(seq)


def dna_molecular_weight(seq: str, double_stranded: bool = False) -> float:
    """Molecular weight of DNA in g/mol.

    Single strand: sum of nucleotide weights - 61.96 (no 5' phosphate).
    Double strand: the strand plus its reverse complement.
    """
    seq = validate(seq)
    if not seq:
        return 0.0
    weight = sum(DNA_WEIGHTS[base] for base in seq) - 61.96
    if double_stranded:
        weight += dna_molecular_weight(reverse_complement(seq))
    return weight


def protein_molecular_weight(protein: str) -> float:
    """Molecular weight of a protein in g/mol: residue weights + one water.

    A trailing stop '*' (as produced by translate) is ignored; any other
    character outside the 20 standard amino acids raises ValueError.
    """
    protein = protein.upper().rstrip("*")
    if not protein:
        return 0.0
    unknown = set(protein) - set(PROTEIN_WEIGHTS)
    if unknown:
        raise ValueError(f"Cannot weigh residues: {', '.join(sorted(unknown))}")
    return sum(PROTEIN_WEIGHTS[aa] for aa in protein) + WATER
