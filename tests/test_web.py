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
        self.assertEqual(result["records"], [{"id": "sequence", "length": 20, "gc": 50.0, "rna": False}])
        self.assertEqual(result["format"], "fasta")

    def test_errors(self):
        self.assertFalse(call("load_text", "   ")["ok"])
        self.assertIn("before the first", call("load_text", "ACGT\n>a\nACGT")["error"])
        bad = call("load_text", ">r\nAUGZ")
        self.assertFalse(bad["ok"])
        self.assertIn("(or U)", bad["error"])

    def test_rna_is_read_as_dna(self):
        result = call("load_text", ">r\nAUGGCUUAA\n>d\nATG")
        self.assertTrue(result["ok"])
        self.assertEqual([(r["id"], r["rna"]) for r in result["records"]], [("r", True), ("d", False)])
        protein = call("run", "translation", "r", json.dumps({"frame": "1"}))["protein"]
        self.assertEqual(protein, "MA*")

    def test_failed_load_keeps_previous_records(self):
        call("load_text", ">keep\nACGT")
        call("load_text", ">bad\nACGZ")
        self.assertTrue(call("run", "basics", "keep")["ok"])


class TestFastq(unittest.TestCase):
    READS = "@a x\nACGTACGTAC\n+\nIIIIII####\n@a\nGGGGCCCC\n+\nIIIIIIII\n"

    def test_load_fastq(self):
        r = call("load_text", self.READS)
        self.assertTrue(r["ok"])
        self.assertEqual(r["format"], "fastq")
        self.assertFalse(r["truncated"])
        self.assertEqual([x["id"] for x in r["records"]], ["a", "a_2"])  # duplicate id kept
        self.assertTrue(call("run", "basics", "a")["ok"])  # reads work in the other tools

    def test_fasta_resets_qualities(self):
        call("load_text", self.READS)
        call("load_text", ">s\nACGT")
        self.assertEqual(call("run", "quality", "s"), {"ok": True, "fastq": False})

    def test_quality_and_trimming(self):
        call("load_text", self.READS)
        r = call("run", "quality", "", json.dumps({"trim_q": "20", "trim_window": "2", "min_length": "5"}))
        self.assertEqual((r["summary"]["reads"], r["summary"]["bases"]), (2, 18))
        self.assertEqual(len(r["positions"]), 10)
        self.assertEqual(r["trim"]["min_length"], 5)
        self.assertEqual((r["trim"]["after"]["reads"], r["trim"]["after"]["bases"]), (2, 14))
        out = call("trimmed_fastq", json.dumps({"trim_q": "20", "trim_window": "2", "min_length": "7"}))
        self.assertEqual(out["reads"], 1)
        self.assertEqual(out["text"], "@a_2\nGGGGCCCC\n+\nIIIIIIII\n")
        self.assertFalse(call("run", "quality", "", json.dumps({"trim_q": "0"}))["ok"])

    def test_bad_fastq(self):
        r = call("load_text", "@a\nACGT\n+\nII")
        self.assertFalse(r["ok"])
        self.assertIn("looks like FASTQ", r["error"])

    def test_read_cap(self):
        many = "".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(bridge.MAX_READS + 5))
        r = call("load_text", many)
        self.assertTrue(r["truncated"])
        self.assertEqual(len(r["records"]), bridge.MAX_READS)

    def test_sample_file(self):
        text = (Path(__file__).resolve().parent.parent / "data" / "samples" / "lambda_reads.fastq").read_text()
        r = call("load_text", text)
        self.assertEqual(len(r["records"]), 300)
        q = call("run", "quality", "", "{}")
        self.assertEqual(q["summary"]["max_length"], 150)
        self.assertLess(q["positions"][-1]["mean"], q["positions"][0]["mean"])  # 3' quality drop


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

    def test_enzyme_hits_have_cut_positions(self):
        r = self.run_tool("motif", "seq3", motif="notI")
        self.assertEqual(r["enzyme"], {"name": "NotI", "site": "GC^GGCCGC", "ends": "4-nt 5' overhang"})
        for hit in r["hits"]:
            expected = hit["start"] + 1 if hit["strand"] == "+" else hit["end"] - 2
            self.assertEqual(hit["cut"], expected)
        self.assertNotIn("enzyme", self.run_tool("motif", "seq3", motif="GCGGCCGC"))

    def test_composition(self):
        r = self.run_tool("composition", window="20")
        self.assertEqual(r["window"], 20)
        self.assertEqual(r["windows"][0]["start"], 1)
        self.assertTrue(all(w["end"] - w["start"] + 1 == 20 for w in r["windows"]))
        self.assertEqual(r["islands"], [])
        self.assertFalse(self.run_tool("composition", window="0")["ok"])

    def test_composition_window_capped_for_long_sequences(self):
        call("load_text", ">long\n" + "ACGT" * 5000)
        r = self.run_tool("composition", "long", window="10")
        self.assertLessEqual(len(r["windows"]), bridge.MAX_WINDOWS)

    def test_kmers(self):
        r = self.run_tool("kmers", "seq2", k="2", top="2")
        self.assertEqual((r["total"], r["distinct"]), (43, 4))
        self.assertEqual([t["kmer"] for t in r["top"]], ["AT", "TA"])

    def test_properties(self):
        r = self.run_tool("properties")
        self.assertEqual(r["tm_method"], "GC formula")
        call("load_text", ">n\nACGTNACGT")
        n = self.run_tool("properties", "n")
        self.assertIsNone(n["tm_nn"])
        self.assertIn("N bases", n["tm_nn_error"])
        call("load_text", ">p\nGTAAAACGACGGCCAGT")
        p = self.run_tool("properties", "p", na="50", oligo="250")
        self.assertEqual(p["tm_nn"], 53.9)
        self.assertGreater(self.run_tool("properties", "p", na="200")["tm_nn"], p["tm_nn"])
        self.assertFalse(self.run_tool("properties", "p", oligo="0")["ok"])
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
