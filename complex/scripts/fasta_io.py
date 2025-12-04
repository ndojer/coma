"""
A script for reading and writing files in FASTA format.
"""
from typing import List, Tuple

def read_fasta(path: str) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    with open(path, "r", encoding="utf-8") as fh:
        header, buf = None, []
        for ln in fh:
            ln = ln.rstrip()
            if ln.startswith(">"):
                if header is not None:
                    records.append((header, "".join(buf).upper()))
                header, buf = ln[1:], []
            else:
                buf.append(ln)
        if header is not None:
            records.append((header, "".join(buf).upper()))
    return records

def write_fasta_record(handle, header: str, sequence: str, width: int = 60) -> None:
    handle.write(f">{header}\n")
    for i in range(0, len(sequence), width):
        handle.write(sequence[i:i+width] + "\n")
