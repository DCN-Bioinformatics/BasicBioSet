# bioseq website: step-by-step guide

*What each part does, what it is for, and the biology behind it.*

**[▶ Open the website](https://dcn-bioinformatics.github.io/BasicBioSet/web/)** · [Lee esta guía en español](Guia_web_bioseq.md)

## Contents

- [1. What bioseq is](#1-what-bioseq-is)
- [2. Biology basics](#2-biology-basics)
  - [2.1 DNA](#21-dna)
  - [2.2 RNA](#22-rna)
  - [2.3 Genes, codons and proteins](#23-genes-codons-and-proteins)
  - [2.4 File formats: FASTA and FASTQ](#24-file-formats-fasta-and-fastq)
- [3. Getting started](#3-getting-started)
  - [3.1 Opening the website](#31-opening-the-website)
  - [3.2 Layout](#32-layout)
  - [3.3 The Your sequences panel](#33-the-your-sequences-panel)
  - [3.4 Common buttons in the results](#34-common-buttons-in-the-results)
- [4. The analysis tabs, one by one](#4-the-analysis-tabs-one-by-one)
  - [The Overview tab](#the-overview-tab)
  - [The Translation tab](#the-translation-tab)
  - [The ORFs tab](#the-orfs-tab)
  - [The Codon usage tab](#the-codon-usage-tab)
  - [The Motif search tab](#the-motif-search-tab)
  - [The K-mers tab](#the-k-mers-tab)
  - [The Composition tab](#the-composition-tab)
  - [The Tm & weight tab](#the-tm--weight-tab)
  - [The Read quality tab](#the-read-quality-tab)
  - [The All sequences tab](#the-all-sequences-tab)
- [5. End-to-end examples](#5-end-to-end-examples)
  - [5.1 Analyzing an unknown gene](#51-analyzing-an-unknown-gene)
  - [5.2 Preparing a cloning](#52-preparing-a-cloning)
  - [5.3 Checking sequencing data](#53-checking-sequencing-data)
  - [5.4 Included example files](#54-included-example-files)
- [6. Glossary](#6-glossary)

## 1. What bioseq is

**bioseq** is a toolkit for analyzing DNA and RNA sequences. On the website you paste or upload a sequence and instantly get its composition, its translation into protein, the possible genes it contains, where restriction enzymes cut it, its melting temperature and, for reads from a sequencing machine, how good those reads are.

All the analysis happens **inside your browser**: the website downloads a version of Python (Pyodide) and runs, on your own computer, the same code the terminal version uses. This means:

- **Your sequences never leave your computer.** Nothing is uploaded to a server (the page footer says so too).
- The website and the terminal menu (`python3 main.py`) give **exactly the same results**.
- The first visit downloads about 10 MB and takes a moment; later visits are faster.

## 2. Biology basics

This section covers the minimum you need to understand what the website shows. If you already know the topic, skip to section 3.

### 2.1 DNA

**DNA** (deoxyribonucleic acid) is the molecule that stores the genetic information of almost every living thing. It is a long chain made of four building blocks called **nucleotides** or **bases**:

| Letter | Base | Pairs with | Type |
|---|---|---|---|
| A | Adenine | T (2 hydrogen bonds) | Purine |
| T | Thymine | A | Pyrimidine |
| G | Guanine | C (3 hydrogen bonds) | Purine |
| C | Cytosine | G | Pyrimidine |

DNA has **two strands** wound into a double helix. Each base on one strand pairs with its complement on the other (A with T, G with C), so knowing one strand tells you the other. Strands have a direction: they are read from **5′ to 3′** ("five prime to three prime"), and the opposite strand runs the other way.

DNA length is measured in **base pairs (bp)** or, for a single strand, in **nucleotides (nt)**. Besides A, C, G and T, the website accepts the letter **N**, which means "unknown base" (for example, a position the sequencer could not read).

### 2.2 RNA

**RNA** is a working copy of a stretch of DNA. It is similar, but uses **U (uracil)** instead of T. Copying DNA into RNA is called **transcription**. Messenger RNA (**mRNA**) carries the instructions for building a protein.

The website accepts RNA: if you paste a sequence containing U, it is read as DNA automatically (U becomes T) and is tagged "RNA" in the sequence list.

### 2.3 Genes, codons and proteins

Proteins are chains of **amino acids** (there are 20 kinds). The mRNA is read three bases at a time; each group of three is a **codon** and stands for one amino acid. This mapping is the **genetic code**:

- There are 4 × 4 × 4 = **64 possible codons**. 61 code for amino acids and 3 are **stop codons** (TAA, TAG, TGA), which mean "end of the protein".
- The codon **ATG** codes for methionine (M) and also marks the **start** of most genes.
- Several codons can give the same amino acid (leucine has 6, for example). The code is said to be "degenerate".
- Going from RNA to protein is called **translation**. The website writes amino acids with their one-letter code (M, A, L…) and stops as an asterisk `*`.

Because codons are read in threes, there are three ways to split a sequence into codons, depending on whether you start at base 1, 2 or 3. Each is a **reading frame**. Reading the opposite strand as well gives **six frames** in total.

### 2.4 File formats: FASTA and FASTQ

**FASTA** is the most common format for storing sequences. Each sequence starts with a header line beginning with `>` followed by its name, and the sequence goes underneath:

```
>my_gene optional description
ATGGCTAGCAAAGGAGAAGAACTTTTC...
```

**FASTQ** is the format sequencing machines produce. Each **read** (a DNA fragment read by the machine) takes 4 lines: the name (starting with `@`), the sequence, a `+` line, and a **quality** line with one character per base. That character encodes how sure the machine is that it read the base correctly (see the Read quality tab in section 4).

## 3. Getting started

### 3.1 Opening the website

1. Open a terminal in the project folder and run `python3 main.py --web`.
2. Open `http://127.0.0.1:8000/web/` in your browser.
3. Wait until the status indicator at the top right no longer says **Starting…**. When the engine is ready, the buttons become active.

If port 8000 is busy, use `python3 main.py --web --port 8080`. To let other computers on your network open it, add `--host 0.0.0.0`. If the project is published on GitHub Pages, just open its `.../web/` address.

> **Important:** Double-clicking `web/index.html` to open it from disk does not work: the browser blocks loading the Python files. Always open it through a server (`--web` or GitHub Pages).

### 3.2 Layout

The screen has three areas:

- **Top bar:** the bioseq name and the analysis engine's status indicator.
- **Left panel, Your sequences:** where you enter data and choose which sequence to analyze.
- **Right panel:** the analysis tabs, each tab's controls, and the results.

### 3.3 The Your sequences panel

This is where data is loaded. You can enter it in three ways:

| Element | What it does |
|---|---|
| Text box | Paste a bare sequence (letters only), FASTA text with one or more sequences, or FASTQ text. You can also **drag and drop** a file onto it. |
| **Analyze** | Reads the text in the box and loads the sequences. Numbers and spaces in a bare sequence are ignored, so you can paste sequences copied from GenBank with their numbering. |
| **Upload file** | Opens a file from your computer. Accepts .fasta, .fa, .fna, .ffn, .fas, .seq, .fastq, .fq and .txt. |
| **Example** | Loads an example sequence so you can try the website without your own data. |
| **Clear** | Empties the text box. |
| **Loaded** list | Shows each loaded sequence with its name, length and GC %. Click one to make the tabs analyze it. RNA sequences carry an "RNA" tag. |

If the text has an error (letters that are not bases, or text before the first `>` line), a message under the buttons explains what to fix. Only A, C, G, T (or U) and N are allowed.

With a **FASTQ** file, the website loads up to 20,000 reads (to keep the browser responsive) and opens the Read quality tab directly. Each read also appears in the list, so you can analyze reads one by one with the other tabs.

### 3.4 Common buttons in the results

- **Copy:** copies a sequence to the clipboard.
- **FASTA / CSV:** download the result. CSV files open in Excel, LibreOffice or Google Sheets. Downloads always hold the complete data even when the screen shows a shortened version (sequences longer than 12,000 characters are shown truncated).
- **Maps and charts:** hover over a bar, point or marker to see a label with the details.

## 4. The analysis tabs, one by one

The first eight tabs analyze **the sequence selected** in the list on the left. The last two (Read quality and All sequences) work on **all loaded sequences** at once.

| Tab | In one sentence |
|---|---|
| Overview | Length, how many A, C, G, T it has, GC %, reverse complement and mRNA. |
| Translation | Turns DNA into protein in one or all six reading frames. |
| ORFs | Finds possible genes: stretches from ATG to a stop codon. |
| Codon usage | Counts which codons the genes found use for each amino acid. |
| Motif search | Finds a short sequence or a restriction enzyme's cut site. |
| K-mers | Counts the most frequent "words" of k letters. |
| Composition | Charts of GC, GC skew and complexity along the sequence; CpG islands. |
| Tm & weight | Melting temperature of a primer and molecular weight. |
| Read quality | Quality of reads from a sequencer (FASTQ) and quality trimming. |
| All sequences | Summary table of every loaded sequence. |

### The Overview tab

#### What does it do?

The basic profile of the sequence: how many bases it has and of which kind, its GC content, its reverse complement and its messenger RNA version.

#### What is it for?

- Check at a glance that the loaded sequence is the right one (expected length, no unexpected Ns).
- Get the **reverse complement**, which you need, for example, to design the reverse primer of a PCR or to look for a gene on the other strand.
- Get the matching **mRNA**.

#### A little biology

**GC content** is the percentage of bases that are G or C. It matters because a G–C pair is held by 3 hydrogen bonds and an A–T pair by only 2: the more GC, the more stable the DNA and the higher the temperature needed to separate its strands. Each species has a typical GC (humans are around 41%, some bacteria exceed 70%), so it also hints at where a fragment comes from.

The **reverse complement** is the other DNA strand, also written 5′ to 3′: reverse the sequence and swap each base for its partner (A↔T, G↔C).

#### How to use it, step by step

1. Load a sequence and select it in the list on the left.
2. Open the **Overview** tab (it is shown by default).

#### How to read the results

- **Length:** length in base pairs (bp).
- **GC content:** percentage of G + C.
- **A, C, G, T (and N if present):** count of each base and its percentage. The colored bar shows each one's share.
- **Reverse complement:** the reverse complement strand, with buttons to copy it or download it as FASTA.
- **mRNA:** the transcribed sequence (T becomes U).

### The Translation tab

#### What does it do?

Translates the DNA sequence into protein with the standard genetic code, in the reading frame you choose or in all six at once.

#### What is it for?

- Find out which protein a gene encodes.
- Find which reading frame the gene is in: the right frame is usually the one that gives a long protein with no asterisks (stops) in the middle.
- Check whether a mutation changes an amino acid or introduces a premature stop.

#### A little biology

Recall (section 2.3) that the sequence is read three bases at a time. Frames **+1, +2 and +3** start reading at base 1, 2 and 3 of the sequence as written; frames **−1, −2 and −3** do the same on the reverse complement. A real gene lies in only one of them.

#### How to use it, step by step

1. Select the sequence and open **Translation**.
2. Under **Reading frame** choose **+1**, **+2**, **+3** or **All six**.
3. Optional: tick **Stop at first stop codon** so the protein ends at the first stop, as it does in the cell. This option is not available when viewing all six frames.

#### How to read the results

- The protein is shown with one-letter amino acid codes.
- Stop codons appear as a **red asterisk** (`*`).
- Codons containing an N, which cannot be translated, appear as **X**.
- In six-frame mode you can download everything as a single FASTA file.

### The ORFs tab

#### What does it do?

Automatically finds **ORFs** (open reading frames): stretches that start with an ATG codon and end at a stop codon, in all three frames on both strands.

#### What is it for?

- Locate possible **genes** in an unknown sequence (for example, a virus genome).
- See where each gene starts and ends, which strand it is on and what protein it would make.
- Get the proteins as FASTA to search them later in databases such as BLAST or UniProt.

#### A little biology

A protein-coding gene is, at its core, an ORF: it starts at ATG and runs without interruption to a stop. Many very short ORFs appear by chance and are not real genes, which is why a minimum length is set. The longer an ORF, the more likely it is a real gene (in eukaryotic genomes, introns split genes into pieces, so this simple search does not find them whole).

#### How to use it, step by step

1. Select the sequence and open **ORFs**.
2. Set **Minimum ORF length** in nucleotides; the default is 30 nt (10 amino acids). For real genomes, raise it to 300 nt (100 amino acids) or more.

#### How to read the results

- A **map** of the sequence with each ORF drawn as an arrow on the + or − strand (hover for details).
- A **table** with: start and end position (counting from 1 on the forward strand), strand (+ or −), length in nucleotides, protein length in amino acids, its **molecular weight in kDa** and the protein sequence (click to expand long ones).
- Buttons to download the table as **CSV** and the proteins as **FASTA**.
- With many ORFs, the table shows the first 500; the CSV holds all of them.
- If none are found, the website says so and suggests lowering the minimum length.

### The Codon usage tab

#### What does it do?

Counts how often each of the 64 codons appears inside the ORFs found (with the same minimum length as the ORFs tab) and groups them by amino acid.

#### What is it for?

- See which codons an organism "prefers" for each amino acid.
- **Optimize genes** for expression in another organism (for example, producing a human protein in the bacterium E. coli): it helps to use the codons the target organism uses most.
- Spot genes acquired from another species (horizontal gene transfer), which often have codon usage that differs from the rest of the genome.

#### A little biology

Since several codons encode the same amino acid (synonymous codons), each species uses them at different frequencies, partly because it has more transfer RNA for some than for others. This preference is called **codon usage bias**.

#### How to use it, step by step

1. Select the sequence and open **Codon usage**.
2. Adjust the **minimum ORF length** if needed (it is the same value as in the ORFs tab).

#### How to read the results

- One row per amino acid (three-letter name and one-letter code) with the total number of codons counted. Stops are included too.
- For each codon: how many times it appears and what **percentage** it represents among that amino acid's codons, with a bar. Codons that never appear are dimmed.
- Download as **CSV**.

### The Motif search tab

#### What does it do?

Finds every occurrence of a short sequence (a **motif**) on both strands. It also recognizes the names of **restriction enzymes** and some common motifs, and for enzymes it shows where they cut.

#### What is it for?

- Plan a **cloning**: find out whether an enzyme cuts your sequence, how many times and where.
- Find known regulatory signals, such as the **TATA box** in promoters.
- Search for any sequence of interest (a primer, a mutation, a tag).

#### A little biology

**Restriction enzymes** are bacterial proteins that cut DNA only where they find a specific sequence, their **recognition site** (for example, EcoRI recognizes GAATTC and cuts between the G and the A). In the lab they are used as "scissors" to cut and paste DNA fragments. Depending on where they cut, they leave **sticky ends** (a single-stranded overhang that makes joining fragments easier) or **blunt ends**.

**IUPAC codes** let you search for motifs with variable positions: R = A/G, Y = C/T, W = A/T, S = C/G, K = G/T, M = A/C and N = any base. For example, the TATA box is searched as `TATAWAW`.

#### How to use it, step by step

1. Select the sequence and open **Motif search**.
2. Type into **Motif or site name** a sequence (`GAATTC`), a motif with IUPAC codes (`TATAWAW`) or an enzyme name (`EcoRI`).
3. Press **Search**, or click one of the quick buttons below it for the most common enzymes and motifs.

#### How to read the results

- A summary with the number of matches on the + and − strands.
- A **map** with a marker at each match.
- A **table** with the position, strand and the exact sequence found (useful with IUPAC codes). If you searched for an enzyme, a **Cut** column shows which two bases it cuts between, and a note explains the kind of ends it leaves.
- Enzymes included: EcoRI, BamHI, HindIII, NotI, EcoRV, XhoI, PstI, SmaI, KpnI, SacI, XbaI, SalI, NcoI, NdeI, SpeI and BglII, plus the TATA-box motif.
- Most restriction sites are **palindromes** (they read the same on both strands, like GAATTC), so each site is found once on the + strand and once on the −: the number of matches is twice the number of sites.

> **Tip:** To clone a gene into a plasmid, choose enzymes that do **not** cut inside the gene. Search each candidate and discard those with matches in it.

### The K-mers tab

#### What does it do?

Counts every "word" of **k** consecutive letters in the sequence (overlapping) and shows the most frequent ones. For example, with k = 3 the words in `ATGCA` are ATG, TGC and GCA.

#### What is it for?

- Detect **repeats** and sequences that occur more often than expected.
- Compare sequences or identify organisms by their k-mer "signature" (many modern genome assembly and classification tools are built on this).
- See how diverse a sequence is: how many distinct words it has out of the possible ones.

#### How to use it, step by step

1. Select the sequence and open **K-mers**.
2. Choose the **k-mer length** (1 to 20). With k = 1 you get the frequency of each base; with k = 2, of each pair (dinucleotides), and so on.
3. Choose under **Show top** how many words you want to see (the most frequent).

#### How to read the results

- The total number of k-mers, how many are distinct and how many are possible (4 to the power of k).
- A table sorted by frequency with the k-mer, its count and its percentage, with a bar.
- Download as **CSV**.

### The Composition tab

#### What does it do?

Analyzes how the sequence's composition changes **along its length**, by splitting it into sliding windows, and looks for **CpG islands**.

#### What is it for?

- See GC-rich or GC-poor regions (they may correspond to genes, to DNA of foreign origin, or to regions that are hard to sequence or amplify by PCR).
- Locate the **origin of replication** in bacterial and viral genomes using GC skew.
- Detect **repeats and low-complexity stretches**.
- Find **CpG islands**, which in vertebrates often mark the start of genes.

#### A little biology

**Sliding window:** instead of one value for the whole sequence, the value is computed over a stretch (the window, for example 100 bp), then the window moves along a little and the value is computed again. This produces a chart along the sequence.

**GC skew** = (G − C) / (G + C). It measures whether one strand has more G than C, ranging from −1 to 1. In many bacteria and phages it changes sign near the origin and terminus of replication, because the two strands are copied differently. **AT skew** is the same with A and T.

**Complexity (Shannon entropy):** measures how varied the composition is, in bits. It is 2 when all four bases appear equally and drops toward 0 in repetitive stretches (for example, AAAAAAA or ATATATAT).

**CpG islands:** "CpG" is a C followed by a G on the same strand. In vertebrate DNA this dinucleotide is rare (its cytosine is often methylated and mutates easily), except in regions called CpG islands, which are frequently found in gene **promoters**. The website uses the classic Gardiner-Garden & Frommer criteria: at least **200 bp**, **≥ 50% GC** and an **observed/expected CpG ratio ≥ 0.6**. Bacteria and phages do not deplete CpG, so many "islands" show up there with no biological meaning.

#### How to use it, step by step

1. Select the sequence and open **Composition**.
2. Set the **Window size** in bp. Small windows show more detail; large windows smooth out noise in long sequences.

#### How to read the results

- Tiles with the overall values: GC skew, AT skew, complexity (maximum 2 bits), observed/expected CpG ratio (about 0.2–0.3 in typical human DNA) and number of CpG islands.
- Three charts along the sequence: **GC %** (with CpG islands shaded), **GC skew** and **complexity** (dips mark repeats).
- If there are CpG islands, a table with their position, length, GC % and observed/expected ratio.
- Download the per-window values as **CSV**.

### The Tm & weight tab

#### What does it do?

Calculates the sequence's **melting temperature (Tm)** and its **molecular weight**, as single-stranded DNA, double-stranded DNA, and for the protein of the longest ORF.

#### What is it for?

- **Design PCR primers:** the two primers of a pair should have similar Tm values, usually between 52 and 65 °C, and the PCR annealing temperature is set a few degrees below.
- Convert between mass and number of molecules (for example, to prepare DNA dilutions or work out moles of an oligonucleotide).
- Estimate the size of the expected protein (for example, to know which band it will appear at on an electrophoresis gel).

#### A little biology

**Tm** is the temperature at which half of the double-stranded DNA molecules have separated into single strands. It depends on the sequence (G–C pairs add more stability), on the salt concentration (Na⁺ stabilizes the double helix) and on the primer concentration.

**Nearest-neighbor method:** the one primer suppliers use. It adds up the energy of each pair of neighboring bases (SantaLucia 1998 parameters) and corrects for salt and primer concentration. This is the recommended value.

**Basic formula:** a quick estimate that ignores salt. Below 14 nt it uses the Wallace rule, 2·(A+T) + 4·(G+C); from 14 nt up, 64.9 + 41·(G+C − 16.4)/length.

Both calculations are meant for **oligonucleotides up to about 60 nt**. For longer sequences the website warns that the value is not meaningful.

#### How to use it, step by step

1. Load the primer sequence (for example, paste it on its own into the box) and open **Tm & weight**.
2. Set **Na⁺** (sodium concentration, in mM) and **Primer** (primer concentration, in nM) to match your protocol. The defaults (50 mM and 250 nM) are typical for PCR.
3. The **minimum ORF length** field only affects the protein weight tile.

#### How to read the results

- **Tm, nearest-neighbor:** the recommended Tm, in °C (not computed if the sequence has Ns or fewer than 2 bases; the website explains why).
- **Tm, basic formula:** the quick estimate, stating which formula was used.
- **DNA weight, single / double strand:** weight of single- and double-stranded DNA, in kDa and g/mol.
- **Longest ORF protein:** weight of the protein from the longest ORF, in kDa, with its length and position.

> **Tip:** The file `data/samples/primers.fasta` contains example primers for practicing with this tab.

### The Read quality tab

#### What does it do?

Analyzes the **quality** of the reads in a **FASTQ** file from a sequencer, lets you **trim** the low-quality parts and download the clean reads. This tab only works with FASTQ files; with FASTA it shows a notice.

This is the tab that answers the question "can I trust my sequencing data?".

#### What is it for?

- Check whether a sequencing run went well **before** using the data in other analyses (assembly, variant calling, etc.). Poor-quality data gives wrong results.
- See **where** quality drops within the reads (usually at the end).
- **Clean** the data: cut off bad ends and discard reads that end up too short.

#### A little biology

**What is sequencing?** Sequencing machines (such as Illumina's) do not read a genome in one go: they break it into millions of fragments and read each fragment base by base, producing **reads** of, say, 100–300 bases. At each position the machine measures a signal (usually light) and decides which base it is. Sometimes the signal is ambiguous and the machine can get it wrong.

**How is quality measured?** That is why, next to each base, the machine records how sure it is, as a **Phred quality score (Q)**. It is defined as Q = −10 · log₁₀(P), where P is the probability that the base was read wrongly. In practical terms: **a score of Q means a 1 in 10 to the power of Q/10 chance of error**. Every 10 extra points make an error 10 times less likely:

| Phred quality | Chance of error | Accuracy | Interpretation |
|---|---|---|---|
| Q10 | 1 in 10 | 90% | Poor |
| Q20 | 1 in 100 | 99% | Acceptable (common minimum threshold) |
| Q30 | 1 in 1,000 | 99.9% | Good (the Illumina benchmark) |
| Q40 | 1 in 10,000 | 99.99% | Excellent |

In the FASTQ file, each score is stored as one character (Phred+33 encoding: the character `I` means Q40, `5` means Q20, `+` means Q10). The website does this conversion for you.

**Why does quality drop at the end of reads?** In sequencing by synthesis, with every cycle some molecules in a cluster get ahead or fall behind, so the signal gets mixed and weaker. That is why Illumina reads normally lose quality toward their 3′ end.

**What is trimming?** It means cutting off the low-quality end of each read. The website uses the **sliding-window** method of the Trimmomatic program (SLIDINGWINDOW): it scans the read from the start with a 4-base window and, as soon as the window's mean quality falls below the threshold (Q20), cuts the read there (keeping the leading bases of that window that are themselves above the threshold). It then discards reads left shorter than the minimum (36 bp), since they are too short to be useful.

#### How to use it, step by step

1. Load a **FASTQ** file (with **Upload file** or by dragging it in). The website opens this tab automatically. To try it, use `data/samples/lambda_reads.fastq`.
2. Review the summary tiles and the two charts.
3. Adjust the trimming controls: **Trim below** (default Q20), **Window** (default 4 bases) and **Keep reads of at least** (default 36 bp).
4. Check the "after trimming" tiles to see how much data is kept and how much quality improves.
5. Press **Trimmed reads (FASTQ)** to download the trimmed reads, ready for other programs.

#### How to read the results

- **Reads:** number of reads loaded.
- **Read length:** length of the reads (or the range, if they vary, with the mean).
- **Mean quality:** mean quality over all bases. Above Q30 is very good.
- **≥ Q30** and **≥ Q20:** percentage of bases with quality of at least 30 (1 error in 1,000) and 20 (1 error in 100). This is the most widely used indicator: a good Illumina run usually has more than 75–80% of bases at ≥ Q30.
- **GC · N:** GC percentage across all reads and percentage of N (unread) bases. A high N percentage signals problems.
- **Quality along the reads chart:** in the style of the FastQC program. For each position within the read it shows a box with the distribution of qualities (the box spans the first to third quartile, the center line is the median, the whiskers run from the 10th to the 90th percentile). The background is colored: **green** = good (Q28 or more), **amber** = acceptable, **red** = poor (below Q20). If the boxes sink into the red zone toward the right, trimming is a good idea. In long reads, positions are grouped into bins.
- **Reads by mean quality histogram:** how many reads have each mean quality. Ideally there is a single peak in the green zone; a tail in the red means a group of bad reads.
- **Quality trimming:** tiles with reads kept, % of bases kept, mean length after trimming, and % ≥ Q30 and mean quality **after**, compared with the values before.
- Downloads: **Per-position CSV** (the data behind the box plot) and **Trimmed reads (FASTQ)** (the clean reads).

> **Tip:** If trimming throws away a very high percentage of reads, the problem is the sequencing run, not the trimming. If little is lost and % ≥ Q30 goes up, trimming has done its job.

> **Also:** Each read in the FASTQ is also loaded as a sequence in the list on the left, so you can select a specific read and analyze it with the other tabs.

### The All sequences tab

#### What does it do?

Shows a summary table of **all** loaded sequences at once (not just the selected one).

#### What is it for?

- Compare several sequences at a glance (length, GC, number of genes).
- Get a complete report to save or share. It is the same report the terminal version produces with `python3 main.py file.fasta`.

#### How to use it, step by step

1. Load a file with several sequences (for example, `data/samples/human_genes.fasta`, with the human insulin and p53 mRNAs).
2. Open **All sequences**.
3. Adjust the **minimum ORF length** to change which ORFs are counted.

#### How to read the results

- One row per sequence with: name, length, GC %, counts of A, C, G and T, number of ORFs found and the **longest ORF** (position, strand and protein length).
- Download the whole table as **CSV**.

## 5. End-to-end examples

### 5.1 Analyzing an unknown gene

1. Paste the sequence and press **Analyze**.
2. In **Overview**, check the length and GC %.
3. In **ORFs**, find the longest ORF: it is probably the gene. Note its position and strand.
4. Download the protein as FASTA and search it in BLAST to find out what it is.
5. In **Translation**, check the protein in the matching frame.

### 5.2 Preparing a cloning

1. Load the gene sequence.
2. In **Motif search**, try the candidate enzymes and discard those that cut inside the gene.
3. Design the primers with the chosen enzymes and check their temperature in **Tm & weight**: both should have a similar Tm.

### 5.3 Checking sequencing data

1. Upload the FASTQ file (for example, `data/samples/lambda_reads.fastq`).
2. In **Read quality**, look at % ≥ Q30 and the per-position quality chart.
3. If quality drops at the end, adjust the trimming and download the trimmed reads.
4. Use that clean file in the next steps of your analysis.

### 5.4 Included example files

| File (in data/samples/) | Use it to try… |
|---|---|
| `lambda_phage.fasta` | Genome of phage lambda: ORFs, composition (GC skew), motifs. |
| `lambda_reads.fastq` | Sequencing reads: the Read quality tab and trimming. |
| `human_genes.fasta` | Human insulin and p53 mRNAs: All sequences, ORFs, protein weights. |
| `primers.fasta` | Short primers: Tm & weight. |
| `rna_example.fasta` | RNA sequences (with U). |
| `feature_tour.fasta` | A tour of several of the website's features. |
| `genbank_style_paste.txt` | A sequence with numbers and spaces, as copied from GenBank. |

## 6. Glossary

| Term | Meaning |
|---|---|
| Amino acid | One of the 20 building blocks of proteins. |
| bp / nt | Base pairs / nucleotides: units of DNA length. |
| Codon | A group of three bases that stands for an amino acid or a stop. |
| CpG / CpG island | A C followed by a G; a CpG-rich region, often in the promoters of vertebrate genes. |
| FASTA / FASTQ | Sequence file formats; FASTQ includes the quality of each base. |
| GC content | Percentage of G and C bases in a sequence. |
| k-mer | A subsequence of k consecutive bases. |
| Motif | A short, biologically meaningful sequence searched for within another. |
| ORF | Open reading frame: from an ATG to a stop codon; a possible gene. |
| PCR | A technique for making millions of copies of a DNA fragment. |
| Phred (Q) | Quality score of a base: Q = −10·log₁₀(probability of error). |
| Primer | A short DNA fragment (about 18–30 nt) that marks where copying starts in a PCR. |
| Read | A DNA fragment read by a sequencer. |
| Reading frame | One of the three ways (six with the other strand) to split a sequence into codons. |
| Restriction enzyme | A protein that cuts DNA at a specific sequence. |
| Reverse complement | The other DNA strand, read 5′ to 3′. |
| Tm | Melting temperature: the temperature at which half of the double-stranded DNA has separated. |
| Transcription / translation | DNA to RNA / RNA to protein. |
| Trimming | Removing the low-quality parts of reads. |
