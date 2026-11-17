"""Analyze DNA sequences, either interactively or as a one-shot report.

Usage:
    python3 main.py                          # interactive menu
    python3 main.py path/to/file.fasta       # full report, no menu
    python3 main.py --min-orf 60             # full report of the default file
    python3 main.py --web                    # open the website version locally
"""

import argparse
import functools
import http.server
import sys
import webbrowser
from pathlib import Path

from bioseq import (
    AMINO_ACID_NAMES,
    CODON_TABLE,
    COMMON_MOTIFS,
    base_counts,
    dna_molecular_weight,
    find_motif,
    find_orfs,
    gc_content,
    kmer_counts,
    melting_temperature,
    orf_codon_usage,
    protein_molecular_weight,
    read_fasta,
    resolve_motif,
    reverse_complement,
    six_frame_translation,
    transcribe,
    translate,
    validate,
)

DEFAULT_FASTA = "data/example.fasta"
DEFAULT_MIN_ORF = 30
MAX_LISTED = 50  # cap on how many motif hits get printed
DEFAULT_PORT = 8000


def preview(seq: str, width: int = 50) -> str:
    return seq if len(seq) <= width else seq[:width] + "..."


def wrap(seq: str, width: int = 60, indent: str = "  ") -> str:
    return "\n".join(indent + seq[i:i + width] for i in range(0, len(seq), width)) or indent + "(empty)"


# ---------------------------------------------------------------- report mode

def print_stats(seq_id: str, seq: str) -> None:
    counts = base_counts(seq)
    print(f"== {seq_id} ==")
    print(f"  Length:              {len(seq)} bp")
    print(f"  Base counts:         " + "  ".join(f"{b}={n}" for b, n in counts.items()))
    print(f"  GC content:          {gc_content(seq):.1f}%")


def print_orfs(seq: str, min_orf: int) -> None:
    orfs = find_orfs(seq, min_length=min_orf)
    if orfs:
        print(f"  ORFs (>= {min_orf} nt):")
        for start, end, strand, protein in orfs:
            print(f"    {start + 1}-{end} ({strand})  ({end - start} nt)  {preview(protein)}")
    else:
        print(f"  ORFs (>= {min_orf} nt): none")


def report(records: dict[str, str], min_orf: int) -> None:
    for seq_id, seq in records.items():
        print_stats(seq_id, seq)
        print(f"  Reverse complement:  {preview(reverse_complement(seq))}")
        print_orfs(seq, min_orf)
        print()


# ------------------------------------------------------------ interactive mode

def ask(prompt: str) -> str:
    return input(prompt).strip()


def ask_int(prompt: str, default: int) -> int:
    value = ask(f"{prompt} [{default}]: ")
    if not value:
        return default
    if value.isdigit() and int(value) > 0:
        return int(value)
    print(f"Not a positive number, using {default}.")
    return default


def load(path: str) -> dict[str, str] | None:
    try:
        records = read_fasta(path)
    except (OSError, ValueError) as e:
        print(f"Could not read {path}: {e}")
        return None
    print(f"Loaded {len(records)} sequence(s) from {path}")
    return records


def choose_sequence(records: dict[str, str]) -> tuple[str, str] | None:
    """Let the user pick a loaded record or type a sequence. Returns (id, seq)."""
    ids = list(records)
    print("\nWhich sequence?")
    for n, seq_id in enumerate(ids, 1):
        print(f"  {n}) {seq_id} ({len(records[seq_id])} bp)")
    print("  t) Type / paste a sequence")
    choice = ask("> ").lower()
    if choice == "t":
        seq = "".join(ask("Sequence: ").split())
        if not seq:
            return None
        try:
            return "typed", validate(seq)
        except ValueError as e:
            print(f"Error: {e}")
            return None
    if choice.isdigit() and 1 <= int(choice) <= len(ids):
        seq_id = ids[int(choice) - 1]
        return seq_id, records[seq_id]
    print("Invalid choice.")
    return None


# Each mode gets the chosen sequence plus the shared settings dict.

def mode_stats(seq_id, seq, settings):
    print_stats(seq_id, seq)


def mode_revcomp(seq_id, seq, settings):
    print(f"Reverse complement of {seq_id}:")
    print(wrap(reverse_complement(seq)))


def mode_transcribe(seq_id, seq, settings):
    print(f"mRNA of {seq_id}:")
    print(wrap(transcribe(seq)))


def mode_translate(seq_id, seq, settings):
    frame = ask("Reading frame (1, 2 or 3) [1]: ") or "1"
    if frame not in ("1", "2", "3"):
        print("Invalid frame.")
        return
    to_stop = ask("Stop at first stop codon? (y/N): ").lower() == "y"
    protein = translate(seq[int(frame) - 1:], to_stop=to_stop)
    print(f"Protein of {seq_id} (frame {frame}):")
    print(wrap(protein))


def mode_six_frames(seq_id, seq, settings):
    print(f"Six-frame translation of {seq_id}  (* = stop codon)")
    for frame, protein in six_frame_translation(seq).items():
        print(f"  Frame {frame}:")
        print(wrap(protein, indent="    "))


def mode_orfs(seq_id, seq, settings):
    settings["min_orf"] = ask_int("Minimum ORF length in nt", settings["min_orf"])
    print(f"== {seq_id} ==")
    print_orfs(seq, settings["min_orf"])


def mode_codon_usage(seq_id, seq, settings):
    min_orf = settings["min_orf"] = ask_int("Minimum ORF length in nt", settings["min_orf"])
    orfs = find_orfs(seq, min_length=min_orf)
    if not orfs:
        print(f"No ORFs >= {min_orf} nt in {seq_id}, so no coding sequence to count.")
        return

    usage = orf_codon_usage(seq, min_length=min_orf)
    synonyms: dict[str, list[str]] = {}
    for codon, aa in CODON_TABLE.items():
        synonyms.setdefault(aa, []).append(codon)

    print(f"Codon usage of {seq_id}: {sum(usage.values())} codons "
          f"from {len(orfs)} ORF(s) >= {min_orf} nt")
    print("  (count, and share among the codons for the same amino acid)")
    for aa in sorted(synonyms, key=lambda a: AMINO_ACID_NAMES[a]):
        total = sum(usage[c] for c in synonyms[aa])
        cells = "  ".join(
            f"{c} {usage[c]:>3} ({usage[c] / total:.2f})" if total else f"{c} {0:>3}   (-)"
            for c in synonyms[aa]
        )
        print(f"  {AMINO_ACID_NAMES[aa]:<4} {aa}  {cells}")


def mode_motif(seq_id, seq, settings):
    print("Motif to search (IUPAC codes allowed, e.g. GAATTC or TATAWAW),")
    print("or one of: " + ", ".join(f"{name} ({m})" for name, m in COMMON_MOTIFS.items()))
    raw = ask("> ")
    if not raw:
        return
    motif = resolve_motif(raw)
    label = motif if motif == raw.upper() else f"{raw} ({motif})"

    hits = find_motif(seq, motif)
    print(f"{label} in {seq_id}: {len(hits)} hit(s)")
    for start, end, strand in hits[:MAX_LISTED]:
        site = seq[start:end] if strand == "+" else reverse_complement(seq[start:end])
        print(f"  {start + 1}-{end} ({strand})  {site}")
    if len(hits) > MAX_LISTED:
        print(f"  ... and {len(hits) - MAX_LISTED} more")


def mode_kmers(seq_id, seq, settings):
    k = ask_int("k-mer length", 3)
    top = ask_int("How many to show", 10)
    counts = kmer_counts(seq, k)
    total = sum(counts.values())
    if not total:
        print(f"{seq_id} is shorter than k={k}.")
        return
    print(f"{k}-mers in {seq_id}: {total} total, {len(counts)} distinct")
    for kmer, n in list(counts.items())[:top]:
        print(f"  {kmer}  {n:>5}  ({100 * n / total:.1f}%)")


def mode_properties(seq_id, seq, settings):
    method = "Wallace rule" if len(seq) < 14 else "GC formula"
    print(f"== {seq_id} ==")
    print(f"  Melting temperature: {melting_temperature(seq):.1f} °C  "
          f"({method}; a rough estimate meant for short primers)")
    print(f"  Molecular weight:    {dna_molecular_weight(seq):,.2f} g/mol (single strand)")
    print(f"                       {dna_molecular_weight(seq, double_stranded=True):,.2f} g/mol (double strand)")
    orfs = find_orfs(seq, min_length=settings["min_orf"])
    if orfs:
        start, end, strand, protein = max(orfs, key=lambda o: o[1] - o[0])
        try:
            weight = f"{protein_molecular_weight(protein):,.2f} g/mol"
        except ValueError:  # protein contains X (codon with an N)
            weight = "n/a (unknown residues)"
        print(f"  Longest ORF protein: {weight} ({len(protein)} aa, ORF {start + 1}-{end} {strand})")
    else:
        print(f"  Longest ORF protein: no ORFs >= {settings['min_orf']} nt")


SECTIONS = [
    ("Basics", [
        ("2", "Base counts & GC content", mode_stats),
        ("3", "Reverse complement", mode_revcomp),
        ("4", "Transcribe (DNA -> mRNA)", mode_transcribe),
    ]),
    ("Translation & genes", [
        ("5", "Translate (one frame)", mode_translate),
        ("6", "Six-frame translation", mode_six_frames),
        ("7", "Find ORFs", mode_orfs),
        ("8", "Codon usage (from ORFs)", mode_codon_usage),
    ]),
    ("Search & properties", [
        ("9", "Motif search", mode_motif),
        ("10", "K-mer counting", mode_kmers),
        ("11", "Melting temperature & molecular weight", mode_properties),
    ]),
]
MODES = {key: func for _, items in SECTIONS for key, _, func in items}


def print_menu(path: str, min_orf: int) -> None:
    print(f"\n===== bioseq main menu =====   file: {path}  |  min ORF: {min_orf} nt")
    print("   1) Full report of all sequences")
    for title, items in SECTIONS:
        print(f"  {title}")
        for key, label, _ in items:
            print(f"  {key:>2}) {label}")
    print("   L) Load a different FASTA file")
    print("   0) Exit")


def menu(path: str, min_orf: int) -> None:
    records = load(path) or {}
    settings = {"min_orf": min_orf}
    while True:
        print_menu(path, settings["min_orf"])
        choice = ask("> ").lower()

        try:
            if choice == "0":
                return
            elif choice == "1":
                print()
                if records:
                    report(records, settings["min_orf"])
                else:
                    print("No sequences loaded.")
            elif choice in MODES:
                picked = choose_sequence(records)
                if picked is not None:
                    print()
                    MODES[choice](*picked, settings)
            elif choice == "l":
                new_path = ask("Path to FASTA file: ")
                new_records = load(new_path) if new_path else None
                if new_records is not None:
                    path, records = new_path, new_records
            else:
                print("Invalid choice, pick an option from the menu.")
        except ValueError as e:  # e.g. an invalid motif or an empty sequence
            print(f"Error: {e}")


# ------------------------------------------------------------------- web mode

class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format, *args):  # keep the terminal clean
        pass


def serve_web(port: int) -> None:
    """Serve the project folder locally and open the web app in a browser."""
    root = Path(__file__).resolve().parent
    handler = functools.partial(QuietHandler, directory=str(root))
    try:
        server = http.server.ThreadingHTTPServer(("127.0.0.1", port), handler)
    except OSError as e:
        sys.exit(f"Could not start the web server on port {port}: {e}\n"
                 f"Try another port, e.g. python3 main.py --web --port {port + 1}")
    url = f"http://127.0.0.1:{port}/web/"
    print(f"bioseq web app running at {url}", flush=True)
    print("Press Ctrl+C to stop.", flush=True)
    webbrowser.open(url)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.")
    finally:
        server.server_close()


def main():
    parser = argparse.ArgumentParser(description="Basic DNA sequence analysis")
    parser.add_argument("fasta", nargs="?",
                        help=f"FASTA file to analyze (default: {DEFAULT_FASTA}); "
                             "giving a file skips the menu and prints a full report")
    parser.add_argument("--min-orf", type=int,
                        help=f"minimum ORF length in nucleotides (default: {DEFAULT_MIN_ORF})")
    parser.add_argument("--menu", action="store_true",
                        help="open the interactive menu even when a file is given")
    parser.add_argument("--web", action="store_true",
                        help="open the website version in your browser")
    parser.add_argument("--port", type=int, default=DEFAULT_PORT,
                        help=f"port for --web (default: {DEFAULT_PORT})")
    args = parser.parse_args()

    if args.web:
        serve_web(args.port)
        return

    path = args.fasta or DEFAULT_FASTA
    min_orf = args.min_orf or DEFAULT_MIN_ORF

    if args.menu or (args.fasta is None and args.min_orf is None):
        try:
            menu(path, min_orf)
        except (EOFError, KeyboardInterrupt):
            print()
        print("Bye!")
        return

    records = read_fasta(path)
    print(f"Read {len(records)} sequence(s) from {path}\n")
    report(records, min_orf)


if __name__ == "__main__":
    sys.exit(main())
