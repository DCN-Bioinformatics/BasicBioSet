"""Sequence file readers (FASTA and FASTQ) and read-quality tools for FASTQ."""

from collections import Counter
from collections.abc import Iterable
from pathlib import Path

PHRED_OFFSET = 33  # Sanger / Illumina 1.8+ quality encoding (Phred+33)
MAX_PHRED = 93


def read_fasta(path: str | Path) -> dict[str, str]:
    """Read a FASTA file and return {record_id: sequence}.

    The record id is the first word of the header line (after '>').
    """
    records: dict[str, str] = {}
    current_id = None
    chunks: list[str] = []

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    records[current_id] = "".join(chunks)
                current_id = line[1:].split()[0]
                chunks = []
            else:
                if current_id is None:
                    raise ValueError(f"{path}: sequence data before first header")
                chunks.append(line.upper())

    if current_id is not None:
        records[current_id] = "".join(chunks)
    return records


# --------------------------------------------------------------------- FASTQ

def is_fastq(text: str) -> bool:
    """True if the text looks like FASTQ: its first non-blank line starts with '@'."""
    for line in text.splitlines():
        if line.strip():
            return line.lstrip().startswith("@")
    return False


def phred_scores(quality: str, offset: int = PHRED_OFFSET) -> list[int]:
    """Decode a FASTQ quality string into Phred scores (Phred+33 by default)."""
    scores = [ord(c) - offset for c in quality]
    if any(not 0 <= q <= MAX_PHRED for q in scores):
        raise ValueError("Quality string has characters outside the Phred+33 range")
    return scores


def phred_string(scores: Iterable[int], offset: int = PHRED_OFFSET) -> str:
    """Encode Phred scores back into a FASTQ quality string."""
    return "".join(chr(q + offset) for q in scores)


def parse_fastq(lines: Iterable[str], max_reads: int | None = None) -> list[tuple[str, str, list[int]]]:
    """Parse FASTQ text (4 lines per read) into (read_id, sequence, phred_scores).

    The read id is the first word of the '@' header. Sequences are uppercased.
    Stops after max_reads reads if given. Raises ValueError, naming the line,
    for malformed records.
    """
    reads = []
    rows = [line.strip() for line in lines]
    while rows and not rows[-1]:
        rows.pop()
    i = 0
    while i < len(rows):
        if not rows[i]:
            i += 1
            continue
        if max_reads is not None and len(reads) >= max_reads:
            break
        header = rows[i]
        if not header.startswith("@") or len(header) < 2:
            raise ValueError(f"Line {i + 1}: expected a '@read_id' header")
        if i + 3 >= len(rows):
            raise ValueError(f"Line {i + 1}: incomplete record (FASTQ needs 4 lines per read)")
        seq, plus, quality = rows[i + 1].upper(), rows[i + 2], rows[i + 3]
        if not plus.startswith("+"):
            raise ValueError(f"Line {i + 3}: expected the '+' separator line")
        if len(seq) != len(quality):
            raise ValueError(f"Line {i + 4}: quality has {len(quality)} characters "
                             f"but the sequence has {len(seq)}")
        try:
            scores = phred_scores(quality)
        except ValueError as e:
            raise ValueError(f"Line {i + 4}: {e}") from None
        reads.append((header[1:].split()[0], seq, scores))
        i += 4
    return reads


def read_fastq(path: str | Path, max_reads: int | None = None) -> list[tuple[str, str, list[int]]]:
    """Read a FASTQ file into (read_id, sequence, phred_scores) tuples."""
    with open(path) as handle:
        return parse_fastq(handle, max_reads)


def format_fastq(reads: Iterable[tuple[str, str, list[int]]]) -> str:
    """Write (read_id, sequence, phred_scores) tuples as FASTQ text."""
    return "".join(f"@{read_id}\n{seq}\n+\n{phred_string(scores)}\n"
                   for read_id, seq, scores in reads)


def trim_read(seq: str, scores: list[int], threshold: int = 20,
              window: int = 4) -> tuple[str, list[int]]:
    """Sliding-window quality trimming, like Trimmomatic's SLIDINGWINDOW.

    Scans from the 5' end and cuts the read at the first window whose mean
    quality is below threshold, keeping the bases before that window plus any
    leading bases of the window that are themselves >= threshold.
    """
    if window < 1:
        raise ValueError("Window must be at least 1")
    window = min(window, len(scores))
    total = sum(scores[:window])
    for start in range(len(scores) - window + 1):
        if start:
            total += scores[start + window - 1] - scores[start - 1]
        if total < threshold * window:
            end = start
            while end < start + window and scores[end] >= threshold:
                end += 1
            return seq[:end], scores[:end]
    return seq, scores


def _quantile(counts: Counter, total: int, fraction: float) -> int:
    """The value at a given fraction of a histogram {value: count}."""
    target = fraction * (total - 1)
    seen = 0
    for value in sorted(counts):
        seen += counts[value]
        if seen > target:
            return value
    return max(counts)


def quality_by_position(score_lists: list[list[int]], bin_size: int = 1) -> list[dict]:
    """Quality distribution along the reads, like FastQC's per-base quality plot.

    Positions are grouped into bins of bin_size bases. Returns one dict per
    bin: start and end (1-based, inclusive), reads (how many reads reach it),
    mean, and the 10th, 25th, 50th, 75th and 90th percentiles.
    """
    if bin_size < 1:
        raise ValueError("Bin size must be at least 1")
    bins: list[Counter] = []
    for scores in score_lists:
        for pos, q in enumerate(scores):
            b = pos // bin_size
            if b == len(bins):
                bins.append(Counter())
            bins[b][q] += 1
    lengths = Counter(len(s) for s in score_lists)
    max_len = max(lengths, default=0)
    rows = []
    for b, counts in enumerate(bins):
        total = sum(counts.values())
        start = b * bin_size
        rows.append({
            "start": start + 1,
            "end": min(start + bin_size, max_len),
            "reads": sum(n for length, n in lengths.items() if length > start),
            "mean": sum(q * n for q, n in counts.items()) / total,
            **{key: _quantile(counts, total, f) for key, f in
               (("p10", .1), ("q1", .25), ("median", .5), ("q3", .75), ("p90", .9))},
        })
    return rows


def fastq_summary(reads: list[tuple[str, str, list[int]]]) -> dict:
    """Headline numbers for a set of reads.

    Returns reads, bases, min/mean/max length, mean quality over all bases,
    the percentage of bases at Q20 and Q30 or above, GC % and N %.
    """
    lengths = [len(seq) for _, seq, _ in reads]
    bases = sum(lengths)
    all_q = Counter(q for _, _, scores in reads for q in scores)
    letters = Counter(b for _, seq, _ in reads for b in seq)
    acgt = sum(letters[b] for b in "ACGT")
    return {
        "reads": len(reads),
        "bases": bases,
        "min_length": min(lengths, default=0),
        "max_length": max(lengths, default=0),
        "mean_length": bases / len(reads) if reads else 0.0,
        "mean_quality": sum(q * n for q, n in all_q.items()) / bases if bases else 0.0,
        "q20": 100 * sum(n for q, n in all_q.items() if q >= 20) / bases if bases else 0.0,
        "q30": 100 * sum(n for q, n in all_q.items() if q >= 30) / bases if bases else 0.0,
        "gc": 100 * (letters["G"] + letters["C"]) / acgt if acgt else 0.0,
        "n": 100 * letters["N"] / bases if bases else 0.0,
    }
