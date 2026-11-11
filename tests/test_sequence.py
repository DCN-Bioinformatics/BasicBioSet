import tempfile
import unittest
from pathlib import Path

from bioseq import (
    base_counts,
    find_orfs,
    gc_content,
    read_fasta,
    reverse_complement,
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
        self.assertEqual(find_orfs(seq, min_length=6), [(2, 17, "MAAA")])
        self.assertEqual(find_orfs(seq, min_length=100), [])


class TestFasta(unittest.TestCase):
    def test_read_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.fasta"
            path.write_text(">a first record\nACGT\nacgt\n\n>b\nTTTT\n")
            self.assertEqual(read_fasta(path), {"a": "ACGTACGT", "b": "TTTT"})


if __name__ == "__main__":
    unittest.main()
