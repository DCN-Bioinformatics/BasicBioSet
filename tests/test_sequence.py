import tempfile
import unittest
from pathlib import Path

from bioseq import (
    base_counts,
    codon_usage,
    dna_molecular_weight,
    find_motif,
    find_orfs,
    gc_content,
    kmer_counts,
    melting_temperature,
    protein_molecular_weight,
    read_fasta,
    reverse_complement,
    six_frame_translation,
    transcribe,
    translate,
)


class TestSequence(unittest.TestCase):
    def test_base_counts(self):
        self.assertEqual(base_counts("AACGTT"), {"A": 2, "C": 1, "G": 1, "T": 2})

    def test_gc_content(self):
        self.assertAlmostEqual(gc_content("GCGC"), 100.0)
        self.assertAlmostEqual(gc_content("ATGC"), 50.0)
        self.assertEqual(gc_content(""), 0.0)

    def test_reverse_complement(self):
        self.assertEqual(reverse_complement("ATGC"), "GCAT")
        self.assertEqual(reverse_complement("aacg"), "CGTT")

    def test_transcribe(self):
        self.assertEqual(transcribe("ATGT"), "AUGU")

    def test_translate(self):
        self.assertEqual(translate("ATGGCTTAA"), "MA*")
        self.assertEqual(translate("ATGGCTTAAGGG", to_stop=True), "MA")
        self.assertEqual(translate("ATGNNN"), "MX")

    def test_invalid_bases(self):
        with self.assertRaises(ValueError):
            gc_content("ATGZ")

    def test_find_orfs(self):
        seq = "CC" + "ATG" + "GCT" * 3 + "TAA" + "GG"
        self.assertEqual(find_orfs(seq, min_length=6), [(2, 17, "+", "MAAA")])
        self.assertEqual(find_orfs(seq, min_length=100), [])

    def test_find_orfs_reverse_strand(self):
        forward = "CC" + "ATG" + "GCT" * 3 + "TAA" + "GG"
        seq = reverse_complement(forward)  # same gene, now on the "-" strand
        self.assertEqual(find_orfs(seq, min_length=6), [(2, 17, "-", "MAAA")])
        self.assertEqual(find_orfs(seq, min_length=6, both_strands=False), [])
        start, end, _, protein = find_orfs(seq, min_length=6)[0]
        self.assertEqual(translate(reverse_complement(seq[start:end]), to_stop=True), protein)

    def test_find_orfs_both_strands(self):
        forward = "ATG" + "GCT" * 3 + "TAA"
        reverse = reverse_complement("ATG" + "AAA" * 2 + "TGA")
        seq = forward + "CC" + reverse
        self.assertEqual(find_orfs(seq, min_length=6), [
            (0, 15, "+", "MAAA"),
            (17, 29, "-", "MKK"),
        ])

    def test_six_frame_translation(self):
        self.assertEqual(six_frame_translation("ATGGCTTAA"), {
            "+1": "MA*", "+2": "WL", "+3": "GL",
            "-1": "LSH", "-2": "*A", "-3": "KP",
        })

    def test_codon_usage(self):
        usage = codon_usage("ATGGCTGCTTAA")
        self.assertEqual(len(usage), 64)
        self.assertEqual((usage["ATG"], usage["GCT"], usage["TAA"]), (1, 2, 1))
        self.assertEqual(sum(usage.values()), 4)
        self.assertEqual(sum(codon_usage("ATGNNNGC").values()), 1)  # N codon and partial codon skipped


class TestMotifs(unittest.TestCase):
    def test_palindromic_site_found_on_both_strands(self):
        self.assertEqual(find_motif("GAATTCAAGAATTC", "GAATTC"), [
            (0, 6, "+"), (0, 6, "-"), (8, 14, "+"), (8, 14, "-"),
        ])

    def test_reverse_strand_hit(self):
        self.assertEqual(find_motif("AAACCC", "AAA"), [(0, 3, "+")])
        self.assertEqual(find_motif("AAACCC", "TTT"), [(0, 3, "-")])
        self.assertEqual(find_motif("AAACCC", "TTT", both_strands=False), [])

    def test_iupac_and_overlaps(self):
        self.assertEqual(find_motif("TATAAAT", "tatawaw", both_strands=False), [(0, 7, "+")])
        self.assertEqual(find_motif("AAAA", "AA", both_strands=False),
                         [(0, 2, "+"), (1, 3, "+"), (2, 4, "+")])

    def test_invalid_motif(self):
        for motif in ("", "GAZ"):
            with self.assertRaises(ValueError):
                find_motif("ACGT", motif)

    def test_kmer_counts(self):
        self.assertEqual(list(kmer_counts("ATATG", 2).items()), [("AT", 2), ("TA", 1), ("TG", 1)])
        self.assertEqual(kmer_counts("AC", 3), {})
        with self.assertRaises(ValueError):
            kmer_counts("ACGT", 0)


class TestProperties(unittest.TestCase):
    def test_melting_temperature(self):
        self.assertEqual(melting_temperature("ATGC"), 12.0)  # Wallace: 2*2 + 4*2
        self.assertAlmostEqual(melting_temperature("ATGC" * 5), 51.78)  # 64.9 + 41*(10-16.4)/20
        with self.assertRaises(ValueError):
            melting_temperature("")

    def test_dna_molecular_weight(self):
        self.assertAlmostEqual(dna_molecular_weight("A"), 251.25)
        self.assertAlmostEqual(dna_molecular_weight("A", double_stranded=True), 251.25 + 242.24)
        self.assertEqual(dna_molecular_weight(""), 0.0)

    def test_protein_molecular_weight(self):
        self.assertAlmostEqual(protein_molecular_weight("MA*"), 131.1926 + 71.0788 + 18.01524)
        self.assertEqual(protein_molecular_weight(""), 0.0)
        with self.assertRaises(ValueError):
            protein_molecular_weight("MXA")


class TestFasta(unittest.TestCase):
    def test_read_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.fasta"
            path.write_text(">a first record\nACGT\nacgt\n\n>b\nTTTT\n")
            self.assertEqual(read_fasta(path), {"a": "ACGTACGT", "b": "TTTT"})


if __name__ == "__main__":
    unittest.main()
