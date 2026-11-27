import tempfile
import unittest
from pathlib import Path

from bioseq import (
    RESTRICTION_ENZYMES,
    at_skew,
    base_counts,
    codon_usage,
    cpg_obs_exp,
    cut_position,
    dna_molecular_weight,
    find_cpg_islands,
    find_enzyme,
    find_motif,
    find_orfs,
    gc_content,
    gc_skew,
    is_rna,
    kmer_counts,
    melting_temperature,
    melting_temperature_nn,
    overhang,
    protein_molecular_weight,
    read_fasta,
    reverse_complement,
    shannon_entropy,
    six_frame_translation,
    sliding_window,
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

    def test_rna_is_read_as_dna(self):
        self.assertTrue(is_rna("augc"))
        self.assertFalse(is_rna("ATGC"))
        self.assertEqual(translate("AUGGCUUAA"), "MA*")
        self.assertEqual(reverse_complement("AUGC"), "GCAT")
        self.assertEqual(gc_content("GCAU"), 50.0)

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

    def test_enzymes_are_palindromes_with_a_cut(self):
        for name, site in RESTRICTION_ENZYMES.items():
            self.assertEqual(site.count("^"), 1, name)
            plain = site.replace("^", "")
            self.assertEqual(plain, reverse_complement(plain), name)

    def test_find_enzyme(self):
        self.assertEqual(find_enzyme(" ecori "), ("EcoRI", "G^AATTC"))
        self.assertIsNone(find_enzyme("GAATTC"))

    def test_cut_position_and_overhang(self):
        seq = "AAGAATTCAA"
        hits = find_motif(seq, "GAATTC")
        self.assertEqual(hits, [(2, 8, "+"), (2, 8, "-")])
        # Top strand G^AATTC cuts after base 3, bottom strand after base 7.
        self.assertEqual([cut_position(*h, "G^AATTC") for h in hits], [3, 7])
        self.assertEqual(overhang("G^AATTC"), "4-nt 5' overhang")
        self.assertEqual(overhang("CTGCA^G"), "4-nt 3' overhang")
        self.assertEqual(overhang("GAT^ATC"), "blunt")

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

    def test_melting_temperature_nn(self):
        # Reference values from Biopython Tm_NN(DNA_NN3, Na=50, dnac1=250, dnac2=0, saltcorr=5).
        self.assertAlmostEqual(melting_temperature_nn("GTAAAACGACGGCCAGT"), 53.86, places=2)
        self.assertAlmostEqual(melting_temperature_nn("TAATACGACTCACTATAGGG"), 48.68, places=2)
        self.assertAlmostEqual(melting_temperature_nn("GCGGCCGC"), 39.58, places=2)  # self-complementary
        self.assertGreater(melting_temperature_nn("GTAAAACGACGGCCAGT", na_mm=200),
                           melting_temperature_nn("GTAAAACGACGGCCAGT"))
        for bad in ("A", "ACGTN"):
            with self.assertRaises(ValueError):
                melting_temperature_nn(bad)
        with self.assertRaises(ValueError):
            melting_temperature_nn("ACGT", oligo_nm=0)

    def test_dna_molecular_weight(self):
        self.assertAlmostEqual(dna_molecular_weight("A"), 251.25)
        self.assertAlmostEqual(dna_molecular_weight("A", double_stranded=True), 251.25 + 242.24)
        self.assertEqual(dna_molecular_weight(""), 0.0)

    def test_protein_molecular_weight(self):
        self.assertAlmostEqual(protein_molecular_weight("MA*"), 131.1926 + 71.0788 + 18.01524)
        self.assertEqual(protein_molecular_weight(""), 0.0)
        with self.assertRaises(ValueError):
            protein_molecular_weight("MXA")


class TestComposition(unittest.TestCase):
    def test_skews(self):
        self.assertEqual(gc_skew("GGGC"), 0.5)
        self.assertEqual(at_skew("ATTT"), -0.5)
        self.assertEqual(gc_skew("ATAT"), 0.0)

    def test_shannon_entropy(self):
        self.assertEqual(shannon_entropy("ACGT"), 2.0)
        self.assertEqual(shannon_entropy("CACACA"), 1.0)
        self.assertEqual(shannon_entropy("AAAA"), 0.0)
        self.assertEqual(shannon_entropy("NNN"), 0.0)

    def test_cpg_obs_exp(self):
        self.assertEqual(cpg_obs_exp("CGCG"), 2.0)  # 2 CpG * 4 / (2 C * 2 G)
        self.assertEqual(cpg_obs_exp("CCGG"), 1.0)
        self.assertEqual(cpg_obs_exp("AAAA"), 0.0)

    def test_sliding_window(self):
        rows = sliding_window("GGGGAAAA", window=4, step=4)
        self.assertEqual([(r["start"], r["end"], r["gc"]) for r in rows], [(0, 4, 100.0), (4, 8, 0.0)])
        self.assertEqual(len(sliding_window("ACGT", window=100)), 1)  # shrinks to the sequence
        self.assertEqual(sliding_window("", window=10), [])
        with self.assertRaises(ValueError):
            sliding_window("ACGT", window=0)

    def test_cpg_islands(self):
        island = "CG" * 150  # 300 bp, 100% GC, plenty of CpG
        seq = "AT" * 200 + island + "AT" * 200
        # A 200 bp window half in the flank is still 50% GC, so the island
        # reaches 100 bp into each AT flank: that is how the window rule works.
        self.assertEqual(find_cpg_islands(seq), [(300, 800)])
        self.assertEqual(find_cpg_islands("GGCCA" * 80), [])  # GC-rich but no CpG
        self.assertEqual(find_cpg_islands("CG" * 50), [])     # shorter than 200 bp

    def test_cpg_islands_meet_criteria_after_merging(self):
        seq = ("CGAT" * 40 + "GGCCA" * 30) * 4  # alternating CpG-rich and CpG-poor
        for start, end in find_cpg_islands(seq):
            region = seq[start:end]
            self.assertGreaterEqual(end - start, 200)
            self.assertGreaterEqual(gc_content(region), 50)
            self.assertGreaterEqual(cpg_obs_exp(region), 0.6)


class TestFasta(unittest.TestCase):
    def test_read_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.fasta"
            path.write_text(">a first record\nACGT\nacgt\n\n>b\nTTTT\n")
            self.assertEqual(read_fasta(path), {"a": "ACGTACGT", "b": "TTTT"})


if __name__ == "__main__":
    unittest.main()
