# bioseq

## Start here · Empieza aquí

1. **Open the website · Abre la web:** **[▶ dcn-bioinformatics.github.io/BasicBioSet/web](https://dcn-bioinformatics.github.io/BasicBioSet/web/)**
2. **Read the guide · Lee la guía:** [📘 English guide](guias/Bioseq_web_guide.docx) · [📗 Guía en español](guias/Guia_web_bioseq.docx)

- **English:** the website runs in your browser, with nothing to install, and your sequences never leave your computer. The step-by-step guide (a Word file in the [`guias/`](guias/) folder) explains every tab, what it is for, and the biology behind it. On GitHub, click the file, then **Download** (⬇) to open it in Word.
- **Español:** la web funciona en tu navegador, sin instalar nada, y tus secuencias nunca salen de tu ordenador. La guía paso a paso (un archivo Word en la carpeta [`guias/`](guias/)) explica cada pestaña, para qué sirve y la biología que hay detrás. En GitHub, haz clic en el archivo y luego en **Download** (⬇) para abrirlo con Word.

---

A small, beginner-friendly bioinformatics project in pure Python (no external dependencies).
It reads DNA sequences from a FASTA file and reports basic statistics.

## Features

- FASTA file parsing (DNA or RNA: U is read as T)
- FASTQ reads: quality summary (mean quality, % of bases ≥ Q20/Q30), per-position quality plot, sliding-window quality trimming and export of the trimmed reads
- Nucleotide counts and GC content
- Reverse complement
- Transcription (DNA → mRNA) and translation (DNA → protein, standard genetic code)
- Six-frame translation
- Open reading frame (ORF) detection on both strands
- Codon usage table (from the ORFs found)
- Motif search on both strands, with IUPAC codes (e.g. `TATAWAW`) and named sites
- 16 common restriction enzymes (EcoRI, BamHI, HindIII, NotI, EcoRV, XhoI, PstI, SmaI, …) with cut positions and the ends they leave
- K-mer counting
- Primer melting temperature by the nearest-neighbor method (SantaLucia 1998, with salt and primer concentration), plus the quick Wallace rule / GC formula
- Composition along the sequence: GC content, GC/AT skew and complexity (Shannon entropy) in sliding windows, and CpG islands
- Molecular weight of DNA (single or double strand) and of proteins
- Interactive terminal menu to explore each feature
- A website version for people who prefer not to use the terminal

## Requirements

Python 3.10+

## Usage

```bash
python3 main.py                          # interactive main menu
python3 main.py my_sequences.fasta       # full report of your own file (no menu)
python3 main.py my_reads.fastq           # read-quality summary of a FASTQ file
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
python3 main.py --web --host 0.0.0.0  # let other computers on your network open it
```

Opening `web/index.html` directly from disk won't work, because browsers block it
from loading the Python files. Use the command above.

**Publish it for free with GitHub Pages:** in the repository go to
*Settings → Pages*, choose *Deploy from a branch*, pick `main` and `/ (root)`,
then share `https://<your-user>.github.io/<repo>/web/`. The empty `.nojekyll` file
in the repository root is required. Without it, GitHub skips `bioseq/__init__.py`.

**Host it on your own server:** see [docs/SELF_HOSTING.md](docs/SELF_HOSTING.md)
for nginx, Apache and quick-test setups, HTTPS, updates and troubleshooting.

## Tests

```bash
python3 -m unittest discover tests
```

## Project layout

```
bioseq/
  fasta.py       FASTA and FASTQ readers, read quality and trimming
  sequence.py    sequence operations (GC, reverse complement, translation, ORFs, codon usage)
  motifs.py      motif search, restriction enzymes and k-mer counting
  properties.py  melting temperature and molecular weight
  composition.py skew, sliding windows, complexity and CpG islands
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

- Restriction digest (fragment sizes and a virtual gel)
- Pairwise alignment (Needleman-Wunsch / Smith-Waterman)
- Compare more results against [Biopython](https://biopython.org/) (the nearest-neighbor Tm already matches it)
