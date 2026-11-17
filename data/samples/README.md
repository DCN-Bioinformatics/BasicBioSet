# Sample sequences for trying bioseq

To load a file on the website, use **Upload file** or drag it onto the left panel.
Every expected result below was checked with bioseq.

## Real sequences (downloaded from NCBI)

### `human_genes.fasta`: human insulin and p53 mRNAs

|                       | INS (NM_000207.3)           | TP53 (NM_000546.6)           |
|-----------------------|-----------------------------|------------------------------|
| Length                | 465 bp                      | 2,512 bp                     |
| GC content            | 63.9 %                      | 53.4 %                       |
| Longest ORF           | 60–392 (+), 110 aa          | 143–1324 (+), 393 aa         |
| Protein starts with   | `MALWMRLLPLLALLALWGPDPAAA…` | `MEEPQSDPSVEPPLSQETFSDLWK…`  |
| Protein weight        | 11.98 kDa                   | 43.65 kDa                    |

- **All sequences:** compares both genes side by side.
- **ORFs:** the longest ORF is the real protein: preproinsulin (110 aa) and p53 (393 aa).
  Set the minimum length to **300**. Insulin then keeps only its real coding sequence.
  p53 keeps 3 ORFs: the real one plus two shorter ORFs that appear by chance.
- **Tm & weight:** p53 is calculated at 43.7 kDa, yet it runs at about 53 kDa on SDS-PAGE gels,
  which is where its name comes from. The calculated weight is the true one.

### `lambda_phage.fasta`: complete genome of phage λ (48,502 bp)

- **Motif search**, using the chips (each site appears once per strand, so matches = 2 × sites):
  - **EcoRI:** 5 sites, at 21,226 · 26,104 · 31,747 · 39,168 · 44,972
  - **BamHI:** 5 sites, at 5,505 · 22,346 · 27,972 · 34,499 · 41,732
  - **HindIII:** 6 sites. Textbook λ HindIII maps show 7. The 7th site (~37,584) is
    missing in this GenBank sequence, which reads `AAGCTC` there instead of `AAGCTT`.
    The λ strain used in labs carries point mutations in that region (the cI gene).
    This is a real example of a reference sequence differing from a lab strain.
- **ORFs:** 697 at the default 30 nt. The table shows the first 500, and the CSV download has all of them.
  At a minimum of 300 nt, 69 remain.
- **Overview:** long sequences are shown truncated after 12,000 bp. Copy and FASTA still give the full sequence.

## Designed sequences

### `feature_tour.fasta`: one record per feature

| Record              | Open this tab      | You should see |
|---------------------|--------------------|----------------|
| `two_strand_genes`  | ORFs               | 2 ORFs: 21–68 on the + strand (`MAEKLVREALGKELA`) and 94–138 on the − strand (`MSGRPEVLQKAFGE`), drawn above and below the axis |
| `restriction_sites` | Motif search       | One site per chip: EcoRI at 19, BamHI at 47, HindIII at 75, NotI at 103 (each listed on both strands). TATA-box: 1 match at 133, + strand only, because it isn't palindromic. Try the IUPAC motif `RGATCY`: it finds the BamHI site. |
| `codon_bias`        | Codon usage        | Leucine is CTG 7 times and CTT once (88 % / 13 %) |
| `unknown_bases`     | Overview, ORFs     | An extra **N** tile (3). The protein reads `MAKEXGLFEARK`, where X comes from the `NNN` codon, and its weight shows "–" |
| `microsatellites`   | K-mers             | k = 2: `CA` and `AC` dominate. k = 4: `CACA` and `ACAC`, then `GATA` from the (GATA)₈ repeat |

### `primers.fasta`: open the **Tm & weight** tab

| Primer             | Length | Tm      | Formula      |
|--------------------|--------|---------|--------------|
| M13_forward_-20    | 17 nt  | 47.1 °C | GC formula   |
| M13_reverse        | 17 nt  | 44.6 °C | GC formula   |
| T7_promoter        | 20 nt  | 47.7 °C | GC formula   |
| short_12mer        | 12 nt  | 36.0 °C | Wallace rule |

Primer suppliers usually quote higher Tm values for these primers. Their calculators use more
detailed (nearest-neighbour) models that also account for salt and primer concentration.
Bonus: `short_12mer` ends in a HindIII site, which you can find with Motif search.

### `genbank_style_paste.txt`: pasting a GenBank listing

Open the file, copy everything, paste it into the box and press **Analyze**. The position numbers
and spaces are ignored, and it loads as one 240 bp sequence (the start of the insulin mRNA, 63.3 % GC).

### `rna_example.fasta`: what an error looks like

Upload it to see the error message. It points out the U bases and suggests replacing U with T.
