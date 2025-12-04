"""
A script generating indexes for indels.
"""

from __future__ import annotations
import random, sys
from typing import Dict, List, Tuple, Optional
from sv_models import SVSpec, SVType

__all__ = ["generate_indels"]

def _interval_free(intervals, start, end, min_distance):
    return all(end + min_distance <= s or start - min_distance >= e for s, e in intervals)

def _rand_dna(rng: random.Random, length: int) -> str:
    alphabet = ("A","C","G","T")
    return "".join(rng.choice(alphabet) for _ in range(length))

def generate_indels(args, seqs: List[Tuple[str,str]], *, sv_id_start: int = 1, existing_specs: Optional[List[SVSpec]] = None) -> List[SVSpec]:
    ins_n = args.ins_count or 0
    del_n = args.del_count or 0
    if ins_n + del_n == 0:
        return []

    rng = random.Random()
    chroms = [hdr for hdr, _ in seqs]
    chrom_len: Dict[str,int] = {hdr: len(seq) for hdr, seq in seqs}
    occupied: Dict[str, List[Tuple[int,int]]] = {c: [] for c in chroms}

    if existing_specs:
        for spec in existing_specs:
            occupied[spec.chr_ref].append((spec.start_ref, spec.end_ref))

            if spec.sv_type == SVType.DUP_TANDEM:
                L = spec.end_ref - spec.start_ref
                occupied[spec.chr_mut].append((spec.end_ref, spec.end_ref + L))

            if spec.sv_type in {
                SVType.TRANS_INTRACHR, SVType.TRANS_INTERCHR,
                SVType.TRANS_INTRA_INV, SVType.TRANS_INTER_INV
            } and spec.insert_pos_ref is not None:
                seg_len = spec.end_ref - spec.start_ref
                occupied[spec.chr_mut].append((spec.insert_pos_ref, spec.insert_pos_ref + seg_len))

            if spec.sv_type == SVType.INS and spec.insert_pos_ref is not None and spec.insert_seq:
                ins_len = len(spec.insert_seq)
                occupied[spec.chr_mut].append((spec.insert_pos_ref, spec.insert_pos_ref + ins_len))

    specs: List[SVSpec] = []
    sv_id = sv_id_start

    if del_n > 0:
        min_len = args.del_min_length or 50
        max_len = args.del_max_length or 1_000
        tries, max_tries = 0, del_n * 40
        while len([s for s in specs if s.sv_type == SVType.DEL]) < del_n and tries < max_tries:
            tries += 1
            chr = rng.choice(chroms)
            L = rng.randint(min_len, min(max_len, chrom_len[chr] - 1))
            start = rng.randint(0, chrom_len[chr] - L)
            end = start + L
            if not _interval_free(occupied[chr], start, end, args.min_distance):
                continue
            specs.append(SVSpec(
                sv_id=sv_id, sv_type=SVType.DEL,
                chr_ref=chr, chr_mut=chr,
                start_ref=start, end_ref=end,
                orientation="+"
            ))
            occupied[chr].append((start, end))
            sv_id += 1
        if len([s for s in specs if s.sv_type == SVType.DEL]) < del_n:
            sys.stderr.write(f"[WARN] placed {len([s for s in specs if s.sv_type == SVType.DEL])}/{del_n} deletions\n")

    if ins_n > 0:
        min_len = args.ins_min_length or 50
        max_len = args.ins_max_length or 1_000
        tries, max_tries = 0, ins_n * 40
        while len([s for s in specs if s.sv_type == SVType.INS]) < ins_n and tries < max_tries:
            tries += 1
            chr = rng.choice(chroms)
            L = rng.randint(min_len, max_len)
            pos = rng.randint(0, chrom_len[chr])
            if not _interval_free(occupied[chr], pos, pos + L, args.min_distance):
                continue
            seq_ins = _rand_dna(rng, L)
            specs.append(SVSpec(
                sv_id=sv_id, sv_type=SVType.INS,
                chr_ref=chr, chr_mut=chr,
                start_ref=pos, end_ref=pos,
                insert_pos_ref=pos, insert_seq=seq_ins,
                orientation="+"
            ))
            occupied[chr].append((pos, pos + L))
            sv_id += 1
        if len([s for s in specs if s.sv_type == SVType.INS]) < ins_n:
            sys.stderr.write(f"[WARN] placed {len([s for s in specs if s.sv_type == SVType.INS])}/{ins_n} insertions\n")

    return specs
