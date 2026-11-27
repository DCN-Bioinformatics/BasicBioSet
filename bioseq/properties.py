"""Physical properties: melting temperature and molecular weight."""

import math

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


# SantaLucia (1998) unified nearest-neighbor parameters for DNA/DNA duplexes,
# in 1 M NaCl: 5'->3' dinucleotide -> (dH kcal/mol, dS cal/(K*mol)).
# The other 6 dinucleotides use their reverse complement (TT = AA, etc.).
NN_PARAMS = {
    "AA": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "GT": (-8.4, -22.4), "CT": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "CG": (-10.6, -27.2), "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9),
}
NN_INIT = {"GC": (0.1, -2.8), "AT": (2.3, 4.1)}  # per terminal base pair
NN_SYMMETRY_DS = -1.4  # self-complementary duplexes
GAS_CONSTANT = 1.9872  # cal/(K*mol)


def melting_temperature_nn(seq: str, na_mm: float = 50.0, oligo_nm: float = 250.0) -> float:
    """Melting temperature (°C) of a primer by the nearest-neighbor method.

    Uses SantaLucia (1998) unified parameters with his entropy salt
    correction, for a primer annealing to its perfect complement:
      Tm = dH / (dS + R*ln(Ct)) - 273.15
    where Ct is the primer concentration. As in PCR, the primer is assumed to
    be in large excess over its target strand. na_mm is the Na+ concentration
    in mM, oligo_nm the primer concentration in nM. Raises ValueError for
    sequences with N or under 2 nt.
    """
    seq = validate(seq)
    if len(seq) < 2:
        raise ValueError("The nearest-neighbor method needs at least 2 nt")
    if "N" in seq:
        raise ValueError("The nearest-neighbor method can't handle N bases")
    if na_mm <= 0 or oligo_nm <= 0:
        raise ValueError("Salt and primer concentrations must be positive")

    dh = ds = 0.0
    for end in (seq[0], seq[-1]):
        h, s = NN_INIT["GC" if end in "GC" else "AT"]
        dh, ds = dh + h, ds + s
    for i in range(len(seq) - 1):
        pair = seq[i:i + 2]
        h, s = NN_PARAMS.get(pair) or NN_PARAMS[reverse_complement(pair)]
        dh, ds = dh + h, ds + s

    if seq == reverse_complement(seq):
        ds += NN_SYMMETRY_DS
    ds += 0.368 * (len(seq) - 1) * math.log(na_mm / 1000)
    return dh * 1000 / (ds + GAS_CONSTANT * math.log(oligo_nm * 1e-9)) - 273.15


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
