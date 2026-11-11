"""Analyze the DNA sequences in a FASTA file.

Usage:
    python3 main.py [path/to/file.fasta] [--min-orf N]
"""

import argparse

from bioseq import (
    base_counts,
    find_orfs,
    gc_content,
    read_fasta,
    reverse_complement,
)


def preview(seq: str, width: int = 50) -> str:
    return seq if len(seq) <= width else seq[:width] + "..."


def main():
    parser = argparse.ArgumentParser(description="Basic DNA sequence analysis")
    parser.add_argument("fasta", nargs="?", default="data/example.fasta",
                        help="FASTA file to analyze (default: data/example.fasta)")
    parser.add_argument("--min-orf", type=int, default=30,
                        help="minimum ORF length in nucleotides (default: 30)")
    args = parser.parse_args()

    records = read_fasta(args.fasta)
    print(f"Read {len(records)} sequence(s) from {args.fasta}\n")

    for seq_id, seq in records.items():
        counts = base_counts(seq)
        print(f"== {seq_id} ==")
        print(f"  Length:              {len(seq)} bp")
        print(f"  Base counts:         " + "  ".join(f"{b}={n}" for b, n in counts.items()))
        print(f"  GC content:          {gc_content(seq):.1f}%")
        print(f"  Reverse complement:  {preview(reverse_complement(seq))}")

        orfs = find_orfs(seq, min_length=args.min_orf)
        if orfs:
            print(f"  ORFs (>= {args.min_orf} nt):")
            for start, end, protein in orfs:
                print(f"    {start + 1}-{end}  ({end - start} nt)  {preview(protein)}")
        else:
            print(f"  ORFs (>= {args.min_orf} nt): none")
        print()


if __name__ == "__main__":
    main()
