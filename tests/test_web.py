"""Tests for web/bridge.py, the layer the website uses to call bioseq."""

import json
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "web"))

import bridge  # noqa: E402

EXAMPLE = (Path(__file__).resolve().parent.parent / "data" / "example.fasta").read_text()


def call(fn, *args):
    return json.loads(getattr(bridge, fn)(*args))


class TestLoadText(unittest.TestCase):
    def test_fasta(self):
        result = call("load_text", EXAMPLE)
        self.assertTrue(result["ok"])
        self.assertEqual([r["id"] for r in result["records"]], ["seq1", "seq2", "seq3"])
        self.assertEqual(result["records"][1]["gc"], 0.0)

    def test_bare_sequence_ignores_spaces_and_numbers(self):
        result = call("load_text", "  1 atgcatgcat gcatg\n 16 caTGC  ")
        self.assertTrue(result["ok"])
        self.assertEqual(result["records"], [{"id": "sequence", "length": 20, "gc": 50.0}])

    def test_errors(self):
        self.assertFalse(call("load_text", "   ")["ok"])
        self.assertIn("before the first", call("load_text", "ACGT\n>a\nACGT")["error"])
        rna = call("load_text", ">r\nAUGC")
        self.assertFalse(rna["ok"])
        self.assertIn("replace U with T", rna["error"])

    def test_failed_load_keeps_previous_records(self):
        call("load_text", ">keep\nACGT")
        call("load_text", ">bad\nACGZ")
        self.assertTrue(call("run", "basics", "keep")["ok"])


class TestRun(unittest.TestCase):
    def setUp(self):
        call("load_text", EXAMPLE)

    def run_tool(self, tool, seq_id="seq1", **params):
        return call("run", tool, seq_id, json.dumps(params))

    def test_basics(self):
        r = self.run_tool("basics")
        self.assertEqual(r["length"], 124)
        self.assertEqual(r["counts"], {"A": 33, "C": 18, "G": 37, "T": 36})

    def test_translation(self):
        self.assertEqual(self.run_tool("translation", frame="1", to_stop=True)["protein"][:3], "PMA")
        self.assertEqual(len(self.run_tool("translation", frame="all")["frames"]), 6)
        self.assertFalse(self.run_tool("translation", frame="7")["ok"])

    def test_orfs(self):
        orfs = self.run_tool("orfs", min_orf="30")["orfs"]
        self.assertEqual(len(orfs), 1)
        self.assertEqual((orfs[0]["start"], orfs[0]["end"], orfs[0]["strand"]), (4, 120, "+"))
        self.assertGreater(orfs[0]["mw"], 0)

    def test_bad_number_is_a_friendly_error(self):
        r = self.run_tool("orfs", min_orf="")
        self.assertFalse(r["ok"])
        self.assertIn("Minimum ORF length", r["error"])

    def test_codons(self):
        r = self.run_tool("codons", min_orf="30")
        self.assertEqual((r["n_orfs"], r["total"]), (1, 39))
        self.assertEqual(len(r["groups"]), 21)

    def test_motif(self):
        r = self.run_tool("motif", "seq3", motif="notI")
        self.assertEqual((r["motif"], r["total"]), ("GCGGCCGC", 4))
        self.assertFalse(self.run_tool("motif", motif="GAZ")["ok"])

    def test_kmers(self):
        r = self.run_tool("kmers", "seq2", k="2", top="2")
        self.assertEqual((r["total"], r["distinct"]), (43, 4))
        self.assertEqual([t["kmer"] for t in r["top"]], ["AT", "TA"])

    def test_properties(self):
        r = self.run_tool("properties")
        self.assertEqual(r["tm_method"], "GC formula")
        self.assertEqual(r["longest_orf"]["aa"], 38)

    def test_protein_with_unknown_residue_has_no_weight(self):
        call("load_text", ">n\nATGNNNGCTGCTGCTGCTGCTGCTGCTTAA")
        orf = self.run_tool("orfs", "n", min_orf="9")["orfs"][0]
        self.assertIn("X", orf["protein"])
        self.assertIsNone(orf["mw"])

    def test_report(self):
        rows = self.run_tool("report", "", min_orf="30")["rows"]
        self.assertEqual([r["n_orfs"] for r in rows], [1, 0, 0])

    def test_unknown_sequence_or_tool(self):
        self.assertFalse(self.run_tool("basics", "nope")["ok"])
        self.assertFalse(self.run_tool("nope")["ok"])


if __name__ == "__main__":
    unittest.main()
