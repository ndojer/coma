"""
A script generating indexes of translocations.
"""
from __future__ import annotations

import random
import sys
from typing import Dict, List, Optional, Tuple

from sv_models import SVSpec, SVType

__all__ = ["generate_translocations"]

def _choose_segment(rng: random.Random, chrom_len: Dict[str, int], chrom: str, *, min_len: int, max_len: int) -> Tuple[int, int]:
    length = rng.randint(min_len, min(max_len, chrom_len[chrom] - 1))
    start = rng.randint(0, chrom_len[chrom] - length)
    end = start + length
    return start, end


def _interval_free(intervals: List[Tuple[int, int]], start: int, end: int, min_distance: int) -> bool:
    return all(end + min_distance <= s or start - min_distance >= e for s, e in intervals)

def generate_translocations(args, seqs: List[Tuple[str, str]], *, sv_id_start: int = 1, existing_specs: Optional[List[SVSpec]] = None) -> List[SVSpec]:
    total = (args.trl_count or 0) + (args.trl_interchr or 0) + (args.invtrl_count or 0) + (args.trl_interchr_inv or 0)
    if total == 0:
        return []

    rng = random.Random()

    chrom_len: Dict[str, int] = {hdr: len(seq) for hdr, seq in seqs}
    chroms = [hdr for hdr, _ in seqs]

    occupied: Dict[str, List[Tuple[int, int]]] = {c: [] for c in chroms}

    if existing_specs:
        for spec in existing_specs:
            occupied[spec.chr_ref].append((spec.start_ref, spec.end_ref))
            if spec.sv_type == SVType.DUP_TANDEM:
                seg_len = spec.end_ref - spec.start_ref
                occupied[spec.chr_mut].append((spec.end_ref, spec.end_ref + seg_len))

    specs: List[SVSpec] = []
    sv_id = sv_id_start

    def _place(n: int, sv_type: SVType, interchr: bool, inverted: bool, min_len: int, max_len: int) -> None:
        nonlocal sv_id
        tries, max_tries = 0, n * 60
        while n > 0 and tries < max_tries:
            tries += 1
            chr_src = rng.choice(chroms)
            start_src, end_src = _choose_segment(rng, chrom_len, chr_src, min_len=min_len, max_len=max_len)
            seg_len = end_src - start_src

            if interchr:
                chr_dst = rng.choice([c for c in chroms if c != chr_src])
            else:
                chr_dst = chr_src

            max_pos_dst = chrom_len[chr_dst] - seg_len
            if max_pos_dst <= 0:
                continue
            pos_dst = rng.randint(0, max_pos_dst)

            if not _interval_free(occupied[chr_src], start_src, end_src, args.min_distance):
                continue
            if not _interval_free(occupied[chr_dst], pos_dst, pos_dst + seg_len, args.min_distance):
                continue
            if chr_src == chr_dst:
                if not (end_src + args.min_distance <= pos_dst or pos_dst + seg_len + args.min_distance <= start_src):
                    continue

            orientation = "-" if inverted else "+"
            spec = SVSpec(
                sv_id=sv_id,
                sv_type=sv_type,
                chr_ref=chr_src,
                chr_mut=chr_dst,
                start_ref=start_src,
                end_ref=end_src,
                insert_pos_ref=pos_dst,
                orientation=orientation,
            )
            specs.append(spec)

            occupied[chr_src].append((start_src, end_src))
            occupied[chr_dst].append((pos_dst, pos_dst + seg_len))
            sv_id += 1
            n -= 1

        if n > 0:
            sys.stderr.write(f"[WARN] placed only {sv_id - sv_id_start}/{total} translocations (constraints too tight?)\n")

    _place(
        n=args.trl_count or 0,
        sv_type=SVType.TRANS_INTRACHR,
        interchr=False,
        inverted=False,
        min_len=args.trl_min_length or 100,
        max_len=args.trl_max_length or 5_000,
    )

    _place(
        n=args.trl_interchr or 0,
        sv_type=SVType.TRANS_INTERCHR,
        interchr=True,
        inverted=False,
        min_len=args.trl_min_length or 100,
        max_len=args.trl_max_length or 5_000,
    )

    _place(
        n=args.invtrl_count or 0,
        sv_type=SVType.TRANS_INTRA_INV,
        interchr=False,
        inverted=True,
        min_len=args.invtrl_min_length or args.trl_min_length or 100,
        max_len=args.invtrl_max_length or args.trl_max_length or 5_000,
    )

    _place(
        n=args.trl_interchr_inv or 0,
        sv_type=SVType.TRANS_INTER_INV,
        interchr=True,
        inverted=True,
        min_len=args.invtrl_min_length or args.trl_min_length or 100,
        max_len=args.invtrl_max_length or args.trl_max_length or 5_000,
    )

    return specs