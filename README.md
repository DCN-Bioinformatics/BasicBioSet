# bioseq

A small, beginner-friendly bioinformatics project in pure Python (no external dependencies).
It reads DNA sequences from a FASTA file and reports basic statistics.

## Features

- FASTA file parsing
- Nucleotide counts and GC content
- Reverse complement
- Transcription (DNA → mRNA) and translation (DNA → protein, standard genetic code)
- Open reading frame (ORF) detection on the forward strand

## Requirements

Python 3.10+

## Usage

```bash
python3 main.py                          # analyzes data/example.fasta
python3 main.py my_sequences.fasta       # your own file
python3 main.py --min-orf 60             # only report ORFs >= 60 nt
```

## Tests

```bash
python3 -m unittest discover tests
```

## Project layout

```
bioseq/
  fasta.py       FASTA reader
  sequence.py    sequence operations (GC, reverse complement, translation, ORFs)
data/
  example.fasta  sample sequences
tests/
  test_sequence.py
main.py          command-line entry point
```

## Ideas for next steps

- Search ORFs on the reverse strand too
- Add k-mer counting or motif search
- Compare against [Biopython](https://biopython.org/) results
- Plot GC content along a sequence with matplotlib
