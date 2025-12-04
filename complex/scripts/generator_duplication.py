"""
A script generating indexes of tandem duplications.
"""

from __future__ import annotations
import random
import sys
from typing import Dict, List, Tuple
from sv_models import SVSpec, SVType

__all__ = ["generate_duplications"]

def generate_duplications(args, seqs: List[Tuple[str, str]], *, sv_id_start: int = 1, existing_specs: List[SVSpec] | None = None) -> List[SVSpec]:

    if not args.dup_count:
        return []

    rng = random.Random()

    chrom_len: Dict[str, int] = {hdr: len(seq) for hdr, seq in seqs}
    occupied: Dict[str, List[Tuple[int, int]]] = {hdr: [] for hdr, _ in seqs}

    if existing_specs:
        for spec in existing_specs:
            occupied[spec.chr_ref].append((spec.start_ref, spec.end_ref))

    min_len = args.dup_min_length or 100
    max_len = args.dup_max_length or 2_000

    specs: List[SVSpec] = []
    sv_id, tries, max_tries = sv_id_start, 0, args.dup_count * 40

    while len(specs) < args.dup_count and tries < max_tries:
        tries += 1
        chrom = rng.choice(seqs)[0]
        length = rng.randint(min_len, min(max_len, chrom_len[chrom] // 2))
        start = rng.randint(0, chrom_len[chrom] - length)
        end = start + length

        total_end = end + length

        if any(
            not (total_end + args.min_distance <= s or start - args.min_distance >= e)
            for s, e in occupied[chrom]
        ):
            continue

        specs.append(
            SVSpec(
                sv_id=sv_id,
                sv_type=SVType.DUP_TANDEM,
                chr_ref=chrom,
                chr_mut=chrom,
                start_ref=start,
                end_ref=end,
                orientation="+",
            )
        )
        occupied[chrom].append((start, total_end))
        sv_id += 1

    if len(specs) < args.dup_count:
        sys.stderr.write(f"[WARN] placed {len(specs)}/{args.dup_count} duplications (constraints too tight?)\n")

    return specs
