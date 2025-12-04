#!/usr/bin/env python3
"""
A script modifying FASTA sequences based on a list of deletions and insertions.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Tuple

from sv_models import MutationEvent, SVSpec, SVType

__all__ = ["apply_mutations"]

def apply_mutations(seqs: List[Tuple[str, str]], events: List[MutationEvent], sv_specs: List[SVSpec]) -> List[Tuple[str, str]]:
    genome: Dict[str, List[str]] = {hdr: list(seq) for hdr, seq in seqs}
    shift_by_chr: Dict[str, int] = defaultdict(int)
    spec_by_id: Dict[int, SVSpec] = {spec.sv_id: spec for spec in sv_specs}

    for ev in events:
        spec = spec_by_id[ev.sv_id]
        chr_shift = shift_by_chr[ev.chr]

        if ev.is_deletion:
            start_mut = ev.start_ref + chr_shift
            end_mut = start_mut + ev.len_ref
            del genome[ev.chr][start_mut:end_mut]
            shift_by_chr[ev.chr] = chr_shift - ev.len_ref
            if spec.sv_type != SVType.INV:
                spec.extra_pos_mut = start_mut

            if spec.sv_type == SVType.DEL:
                spec.start_mut = start_mut
                spec.end_mut   = start_mut

            continue

        spec = spec_by_id[ev.sv_id]
        if spec.sv_type == SVType.INV:
            orig_len = spec.end_ref - spec.start_ref
            start_mut = ev.start_ref + chr_shift + orig_len
        else:
            start_mut = ev.start_ref + chr_shift

        ins_seq = ev.seq_ins
        genome[ev.chr][start_mut:start_mut] = list(ins_seq)
        shift_by_chr[ev.chr] = chr_shift + ev.len_ins

        end_mut_ins = start_mut + ev.len_ins
        old_start, old_end = spec.start_mut, spec.end_mut
        assert spec.start_mut is None
        if spec.start_mut is None or start_mut < spec.start_mut:
            spec.start_mut = start_mut
        assert spec.end_mut is None
        if spec.end_mut is None or end_mut_ins > spec.end_mut:
            spec.end_mut = end_mut_ins

        if spec.sv_type == SVType.DUP_TANDEM:
            orig_len = spec.end_ref - spec.start_ref
            new_start = end_mut_ins - orig_len - ev.len_ins
            spec.extra_pos_mut = spec.start_mut
            spec.start_mut = new_start

    return [(hdr, "".join(genome[hdr])) for hdr, _ in seqs]
