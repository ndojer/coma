"""
The main script responsible for coordinating the pipeline execution.
"""
from __future__ import annotations

from typing import List

from cli import parse_args
from fasta_io import read_fasta, write_fasta_record
from generator_inversion import generate_inversions
from generator_duplication import generate_duplications
from generator_translocation import generate_translocations
from event_builder import build_events
from mutator import apply_mutations
from csv_io import write_csv
from sv_models import SVSpec, SVType
from generator_indel import generate_indels

def main() -> None:
    args = parse_args()

    seqs = read_fasta(args.fasta)
    if not seqs:
        raise SystemExit(f"No sequences found in {args.fasta}")

    sv_specs: List[SVSpec] = []

    sv_specs += generate_inversions(args, seqs)

    sv_specs += generate_duplications(
        args,
        seqs,
        sv_id_start=len(sv_specs) + 1,
        existing_specs=sv_specs,
    )

    sv_specs += generate_translocations(
        args,
        seqs,
        sv_id_start=len(sv_specs) + 1,
        existing_specs=sv_specs,
    )

    sv_specs += generate_indels(
        args, seqs, sv_id_start=len(sv_specs)+1, existing_specs=sv_specs
    )

    if not sv_specs:
        raise SystemExit("No structural variants generated – check parameters.")

    events = build_events(sv_specs, seqs)
    mutated = apply_mutations(seqs, events, sv_specs)

    with open(args.output, "w", encoding="utf-8") as fout:
        for hdr, seq in mutated:
            write_fasta_record(fout, hdr, seq)

    write_csv(args.csv, sv_specs)

if __name__ == "__main__":
    main()
