# bioseq

A small, beginner-friendly bioinformatics project in pure Python (no external dependencies).
It reads DNA sequences from a FASTA file and reports basic statistics.

## Features

- FASTA file parsing
- Nucleotide counts and GC content
- Reverse complement
- Transcription (DNA → mRNA) and translation (DNA → protein, standard genetic code)
- Six-frame translation
- Open reading frame (ORF) detection on both strands
- Codon usage table (from the ORFs found)
- Motif search on both strands, with IUPAC codes (e.g. `TATAWAW`) and named sites (EcoRI, BamHI, HindIII, NotI)
- K-mer counting
- Melting temperature estimate (Wallace rule / GC formula) for primers
- Molecular weight of DNA (single or double strand) and of proteins
- Interactive terminal menu to explore each feature
- A website version for people who prefer not to use the terminal

## Requirements

Python 3.10+

## Usage

```bash
python3 main.py                          # interactive main menu
python3 main.py my_sequences.fasta       # full report of your own file (no menu)
python3 main.py --min-orf 60             # full report, only ORFs >= 60 nt
python3 main.py my.fasta --menu          # open the menu with your own file loaded
```

## Website

The `web/` folder holds a website version of the whole toolkit. It runs the same
`bioseq` Python code inside the browser using [Pyodide](https://pyodide.org), so:

- users only need a link, with nothing to install and no terminal;
- sequences never leave the user's computer (there is no server doing the analysis);
- the website and the terminal menu always give the same results.

Users can paste sequences, upload or drag in FASTA files, view feature maps of
ORFs and motif hits, and download results as CSV or FASTA. The first visit
downloads about 10 MB for the Python engine. Later visits are cached and faster.

**Run it locally:**

```bash
python3 main.py --web              # opens http://127.0.0.1:8000/web/
python3 main.py --web --port 8080  # if port 8000 is busy
```

Opening `web/index.html` directly from disk won't work, because browsers block it
from loading the Python files. Use the command above.

**Publish it for free with GitHub Pages:** in the repository go to
*Settings → Pages*, choose *Deploy from a branch*, pick `main` and `/ (root)`,
then share `https://<your-user>.github.io/<repo>/web/`. The empty `.nojekyll` file
in the repository root is required. Without it, GitHub skips `bioseq/__init__.py`.

## Tests

```bash
python3 -m unittest discover tests
```

## Project layout

```
bioseq/
  fasta.py       FASTA reader
  sequence.py    sequence operations (GC, reverse complement, translation, ORFs, codon usage)
  motifs.py      motif search and k-mer counting
  properties.py  melting temperature and molecular weight
data/
  example.fasta  sample sequences
web/
  index.html     website (open through python3 main.py --web)
  app.js         page logic and rendering
  style.css      styles (light and dark mode)
  bridge.py      JSON layer between the page and bioseq
tests/
  test_sequence.py
  test_web.py    tests for web/bridge.py
main.py          command-line entry point (report, menu, --web)
```

## Ideas for next steps

- Sliding-window GC content
- Restriction digest (cut positions and fragment sizes)
- Pairwise alignment (Needleman-Wunsch / Smith-Waterman)
- Compare against [Biopython](https://biopython.org/) results
- Plot GC content along a sequence with matplotlib
