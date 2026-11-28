"""Analyze DNA sequences, either interactively or as a one-shot report.

Usage:
    python3 main.py                          # interactive menu
    python3 main.py path/to/file.fasta       # full report, no menu
    python3 main.py --min-orf 60             # full report of the default file
    python3 main.py --web                    # open the website version locally
    python3 main.py --web --host 0.0.0.0     # let other computers reach it
"""

import argparse
import errno
import functools
import http.server
import socket
import sys
import webbrowser
from pathlib import Path

from bioseq import (
    AMINO_ACID_NAMES,
    CODON_TABLE,
    COMMON_MOTIFS,
    at_skew,
    base_counts,
    cpg_obs_exp,
    cut_position,
    dna_molecular_weight,
    fastq_summary,
    find_cpg_islands,
    find_enzyme,
    find_motif,
    find_orfs,
    format_fastq,
    gc_content,
    gc_skew,
    is_fastq,
    is_rna,
    kmer_counts,
    melting_temperature,
    melting_temperature_nn,
    orf_codon_usage,
    overhang,
    protein_molecular_weight,
    quality_by_position,
    read_fasta,
    read_fastq,
    resolve_motif,
    reverse_complement,
    shannon_entropy,
    six_frame_translation,
    sliding_window,
    transcribe,
    translate,
    trim_read,
    validate,
)

DEFAULT_FASTA = "data/example.fasta"
DEFAULT_MIN_ORF = 30
MAX_LISTED = 50  # cap on how many motif hits get printed
DEFAULT_PORT = 8000
DEFAULT_HOST = "127.0.0.1"
LOCAL_HOSTS = {"127.0.0.1", "localhost", "::1"}
SPARK = "▁▂▃▄▅▆▇█"


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


def read_sequences(path: str) -> tuple[dict[str, str], dict[str, list[int]]]:
    """Read a FASTA or FASTQ file (detected from its content).

    Returns ({id: sequence}, {id: phred scores}); the second dict is empty for
    FASTA. Reads that share an id get a _2, _3, ... suffix.
    """
    with open(path) as handle:
        fastq = is_fastq(handle.read(10_000))
    if not fastq:
        return read_fasta(path), {}
    records, qualities = {}, {}
    for read_id, seq, scores in read_fastq(path):
        unique, n = read_id, 1
        while unique in records:
            n += 1
            unique = f"{read_id}_{n}"
        records[unique], qualities[unique] = seq, scores
    return records, qualities


def load(path: str) -> tuple[dict[str, str], dict[str, list[int]]] | None:
    try:
        records, qualities = read_sequences(path)
    except (OSError, ValueError, UnicodeDecodeError) as e:
        print(f"Could not read {path}: {e}")
        return None
    kind = "read(s) from FASTQ file" if qualities else "sequence(s) from"
    print(f"Loaded {len(records)} {kind} {path}")
    rna = [seq_id for seq_id, seq in records.items() if is_rna(seq)]
    if rna:
        shown = ", ".join(rna[:5]) + (f" and {len(rna) - 5} more" if len(rna) > 5 else "")
        print(f"RNA read as DNA (U -> T): {shown}")
    return records, qualities


def choose_sequence(records: dict[str, str]) -> tuple[str, str] | None:
    """Let the user pick a loaded record or type a sequence. Returns (id, seq)."""
    ids = list(records)
    print("\nWhich sequence?")
    for n, seq_id in enumerate(ids[:MAX_LISTED], 1):
        print(f"  {n}) {seq_id} ({len(records[seq_id])} bp)")
    if len(ids) > MAX_LISTED:
        print(f"  ... {len(ids) - MAX_LISTED} more: type a number up to {len(ids)}")
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
    print("or a name. Restriction enzymes also show where they cut:")
    names = [f"{name} ({m})" for name, m in COMMON_MOTIFS.items()]
    for i in range(0, len(names), 4):
        print("  " + "  ".join(f"{n:<20}" for n in names[i:i + 4]).rstrip())
    raw = ask("> ")
    if not raw:
        return
    motif = resolve_motif(raw)
    label = motif if motif == raw.upper() else f"{raw} ({motif})"
    enzyme = find_enzyme(raw)

    hits = find_motif(seq, motif)
    print(f"{label} in {seq_id}: {len(hits)} hit(s)")
    if enzyme and hits:
        print(f"  {enzyme[0]} cuts {enzyme[1]} ({overhang(enzyme[1])}); "
              "'cut a^b' = between bases a and b (+ top strand, - bottom strand)")
    for start, end, strand in hits[:MAX_LISTED]:
        site = seq[start:end] if strand == "+" else reverse_complement(seq[start:end])
        cut = ""
        if enzyme:
            pos = cut_position(start, end, strand, enzyme[1])
            cut = f"  cut {pos}^{pos + 1}"
        print(f"  {start + 1}-{end} ({strand})  {site}{cut}")
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


def mode_composition(seq_id, seq, settings):
    window = ask_int("Window size in bp", min(100, max(1, len(seq))))
    print(f"== {seq_id} ==")
    print(f"  GC content:   {gc_content(seq):.1f}%")
    print(f"  GC skew:      {gc_skew(seq):+.3f}   (G-C)/(G+C)")
    print(f"  AT skew:      {at_skew(seq):+.3f}   (A-T)/(A+T)")
    print(f"  Complexity:   {shannon_entropy(seq):.2f} bits (Shannon entropy, max 2)")
    print(f"  CpG obs/exp:  {cpg_obs_exp(seq):.2f}")
    rows = sliding_window(seq, window)
    if len(rows) > 1:
        # Average neighboring windows so each line fits in about 60 characters.
        group = -(-len(rows) // 60)
        for key, label, lo, hi in (("gc", "GC %", 0, 100), ("gc_skew", "GC skew", -1, 1),
                                   ("entropy", "Complexity", 0, 2)):
            values = [sum(r[key] for r in rows[i:i + group]) / len(rows[i:i + group])
                      for i in range(0, len(rows), group)]
            line = "".join(SPARK[min(7, int((v - lo) / (hi - lo) * 8))] for v in values)
            print(f"  {label:<11} {line}  ({lo} to {hi})")
        print(f"  ({len(rows)} windows of {rows[0]['end']} bp, 1 to {len(seq):,})")
    islands = find_cpg_islands(seq)
    print(f"  CpG islands (>= 200 bp, >= 50% GC, obs/exp >= 0.6): {len(islands) or 'none'}")
    for start, end in islands[:MAX_LISTED]:
        region = seq[start:end]
        print(f"    {start + 1}-{end}  ({end - start} bp, {gc_content(region):.1f}% GC, "
              f"obs/exp {cpg_obs_exp(region):.2f})")


def mode_properties(seq_id, seq, settings):
    method = "Wallace rule" if len(seq) < 14 else "GC formula"
    print(f"== {seq_id} ==")
    try:
        print(f"  Tm, nearest-neighbor: {melting_temperature_nn(seq):.1f} °C  "
              "(SantaLucia 1998; 50 mM Na+, 250 nM primer)")
    except ValueError as e:
        print(f"  Tm, nearest-neighbor: n/a ({e})")
    print(f"  Tm, basic formula:    {melting_temperature(seq):.1f} °C  "
          f"({method}; a rough estimate)")
    if len(seq) > 60:
        print("  (Both Tm values are meant for primers up to about 60 nt.)")
    print(f"  Molecular weight:     {dna_molecular_weight(seq):,.2f} g/mol (single strand)")
    print(f"                        {dna_molecular_weight(seq, double_stranded=True):,.2f} g/mol (double strand)")
    orfs = find_orfs(seq, min_length=settings["min_orf"])
    if orfs:
        start, end, strand, protein = max(orfs, key=lambda o: o[1] - o[0])
        try:
            weight = f"{protein_molecular_weight(protein):,.2f} g/mol"
        except ValueError:  # protein contains X (codon with an N)
            weight = "n/a (unknown residues)"
        print(f"  Longest ORF protein:  {weight} ({len(protein)} aa, ORF {start + 1}-{end} {strand})")
    else:
        print(f"  Longest ORF protein:  no ORFs >= {settings['min_orf']} nt")


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
        ("12", "Composition: GC skew, complexity, CpG islands", mode_composition),
    ]),
]
MODES = {key: func for _, items in SECTIONS for key, _, func in items}


# ------------------------------------------------------------ FASTQ quality

def print_quality(records: dict[str, str], qualities: dict[str, list[int]]) -> None:
    reads = [(read_id, records[read_id], qualities[read_id]) for read_id in records]
    s = fastq_summary(reads)
    lengths = (f"{s['max_length']}" if s["min_length"] == s["max_length"]
               else f"{s['min_length']}-{s['max_length']} (mean {s['mean_length']:.1f})")
    print(f"  Reads:          {s['reads']:,}  ({s['bases']:,} bases)")
    print(f"  Read length:    {lengths} bp")
    print(f"  Mean quality:   Q{s['mean_quality']:.1f}")
    print(f"  Bases >= Q30:   {s['q30']:.1f}%   (>= Q20: {s['q20']:.1f}%)")
    print(f"  GC / N:         {s['gc']:.1f}% / {s['n']:.2f}%")
    bin_size = -(-s["max_length"] // 60)  # ceiling division: at most 60 columns
    rows = quality_by_position([q for _, _, q in reads], bin_size)
    line = "".join(SPARK[min(7, max(0, (r["median"] - 2) * 8 // 40))] for r in rows)
    print(f"  Median quality along the reads (Q2 to Q41, {bin_size} base(s) per column):")
    print(f"    {line}")
    worst = min(rows, key=lambda r: r["median"])
    print(f"    lowest median Q{worst['median']} at bases {worst['start']}-{worst['end']}")


def mode_quality(records: dict[str, str], qualities: dict[str, list[int]], path: str) -> None:
    if not qualities:
        print("Read quality needs a FASTQ file. Load one with L (try data/samples/lambda_reads.fastq).")
        return
    print_quality(records, qualities)
    if ask("\nTrim low-quality ends? (y/N): ").lower() != "y":
        return
    threshold = ask_int("Cut where a window's mean quality drops below Q", 20)
    window = ask_int("Window size in bases", 4)
    min_length = ask_int("Drop reads shorter than (bp)", 36)
    kept = []
    for read_id, seq in records.items():
        trimmed, scores = trim_read(seq, qualities[read_id], threshold, window)
        if len(trimmed) >= min_length:
            kept.append((read_id, trimmed, scores))
    before, after = sum(map(len, records.values())), sum(len(seq) for _, seq, _ in kept)
    print(f"Kept {len(kept):,} of {len(records):,} reads and "
          f"{100 * after / before if before else 0:.1f}% of bases.")
    if kept:
        summary = fastq_summary(kept)
        print(f"Quality after trimming: Q{summary['mean_quality']:.1f} mean, "
              f"{summary['q30']:.1f}% of bases >= Q30")
    stem = Path(path).name.rsplit(".", 1)[0]
    out = ask(f"Save trimmed reads to [{stem}.trimmed.fastq] (n to skip): ") or f"{stem}.trimmed.fastq"
    if out.lower() == "n":
        return
    try:
        Path(out).write_text(format_fastq(kept))
        print(f"Saved {len(kept):,} reads to {out}")
    except OSError as e:
        print(f"Could not save {out}: {e}")


def print_menu(path: str, min_orf: int) -> None:
    print(f"\n===== bioseq main menu =====   file: {path}  |  min ORF: {min_orf} nt")
    print("   1) Full report of all sequences")
    for title, items in SECTIONS:
        print(f"  {title}")
        for key, label, _ in items:
            print(f"  {key:>2}) {label}")
    print("  Reads")
    print("  13) Read quality & trimming (FASTQ files)")
    print("   L) Load a different FASTA or FASTQ file")
    print("   0) Exit")


def menu(path: str, min_orf: int) -> None:
    records, qualities = load(path) or ({}, {})
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
            elif choice == "13":
                print()
                mode_quality(records, qualities, path)
            elif choice in MODES:
                picked = choose_sequence(records)
                if picked is not None:
                    print()
                    MODES[choice](*picked, settings)
            elif choice == "l":
                new_path = ask("Path to FASTA or FASTQ file: ")
                loaded = load(new_path) if new_path else None
                if loaded is not None:
                    path, (records, qualities) = new_path, loaded
            else:
                print("Invalid choice, pick an option from the menu.")
        except ValueError as e:  # e.g. an invalid motif or an empty sequence
            print(f"Error: {e}")


# ------------------------------------------------------------------- web mode

class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format, *args):  # keep the terminal clean
        pass


def web_url(host: str, port: int) -> str:
    """URL to show for the web app. A wildcard address can't be opened in a
    browser, so it is shown as localhost (other computers use this machine's IP)."""
    if host in ("0.0.0.0", ""):
        host = "127.0.0.1"
    elif host == "::":
        host = "::1"
    if ":" in host:  # IPv6 literals need brackets in URLs
        host = f"[{host}]"
    return f"http://{host}:{port}/web/"


def serve_web(port: int, host: str = DEFAULT_HOST) -> None:
    """Serve the project folder and open the web app in a browser.

    By default it only listens on this computer. Pass another host (e.g.
    0.0.0.0) to let other computers on the network reach it; the browser is
    then not opened, since that is usually a server with no screen.
    """
    root = Path(__file__).resolve().parent
    handler = functools.partial(QuietHandler, directory=str(root))
    server_class = http.server.ThreadingHTTPServer
    if ":" in host:
        server_class = type("ThreadingHTTPServerV6", (server_class,),
                            {"address_family": socket.AF_INET6})
    try:
        server = server_class((host, port), handler)
    except OSError as e:
        hint = (f"Check the address given to --host ({host})."
                if isinstance(e, socket.gaierror) or e.errno == errno.EADDRNOTAVAIL
                else f"Try another port, e.g. python3 main.py --web --port {port + 1}")
        sys.exit(f"Could not start the web server on {host}:{port}: {e}\n{hint}")
    url = web_url(host, port)
    local = host in LOCAL_HOSTS
    print(f"bioseq web app running at {url}", flush=True)
    if not local:
        print(f"Listening on {host}:{port}: other computers can open "
              f"http://<this-machine's-address>:{port}/web/", flush=True)
        print("This simple server is meant for testing and local networks; "
              "for a public site see docs/SELF_HOSTING.md.", flush=True)
    print("Press Ctrl+C to stop.", flush=True)
    if local:
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
                        help=f"FASTA or FASTQ file to analyze (default: {DEFAULT_FASTA}); "
                             "giving a file skips the menu and prints a full report")
    parser.add_argument("--min-orf", type=int,
                        help=f"minimum ORF length in nucleotides (default: {DEFAULT_MIN_ORF})")
    parser.add_argument("--menu", action="store_true",
                        help="open the interactive menu even when a file is given")
    parser.add_argument("--web", action="store_true",
                        help="open the website version in your browser")
    parser.add_argument("--port", type=int, default=DEFAULT_PORT,
                        help=f"port for --web (default: {DEFAULT_PORT})")
    parser.add_argument("--host", default=DEFAULT_HOST,
                        help=f"address for --web to listen on (default: {DEFAULT_HOST}, "
                             "this computer only; use 0.0.0.0 to allow other computers)")
    args = parser.parse_args()

    if args.web:
        serve_web(args.port, args.host)
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

    records, qualities = read_sequences(path)
    if qualities:
        print(f"Read {len(records):,} reads from FASTQ file {path}\n")
        print_quality(records, qualities)
        return
    print(f"Read {len(records)} sequence(s) from {path}\n")
    report(records, min_orf)


if __name__ == "__main__":
    sys.exit(main())
