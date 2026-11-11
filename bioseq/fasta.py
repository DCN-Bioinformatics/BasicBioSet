"""Minimal FASTA file reader."""

from pathlib import Path


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
