"""
A script generating indexes of inversions.
"""
import random
import sys
from typing import List, Tuple, Dict

from sv_models import SVSpec, SVType

__all__ = ["generate_inversions"]

def generate_inversions(args, seqs: List[Tuple[str, str]]) -> List[SVSpec]:
    if not args.inv_count:
        return []

    rng = random.Random()
    chrom_len = {hdr: len(seq) for hdr, seq in seqs}
    occupied: Dict[str, List[Tuple[int, int]]] = {hdr: [] for hdr, _ in seqs}

    min_len = args.inv_min_length or 50
    max_len = args.inv_max_length or 1_000
    specs: List[SVSpec] = []
    sv_id, tries, max_tries = 1, 0, args.inv_count * 20

    while len(specs) < args.inv_count and tries < max_tries:
        tries += 1
        chrom = rng.choice(seqs)[0]
        length = rng.randint(min_len, min(max_len, chrom_len[chrom] // 2))
        start = rng.randint(0, chrom_len[chrom] - length)
        end = start + length

        if any(not (end + args.min_distance <= s or start - args.min_distance >= e)
               for s, e in occupied[chrom]):
            continue

        specs.append(
            SVSpec(sv_id, SVType.INV, chrom, chrom, start, end, orientation="-")
        )
        occupied[chrom].append((start, end))
        sv_id += 1

    if len(specs) < args.inv_count:
        sys.stderr.write(
            f"[WARN] placed {len(specs)}/{args.inv_count} inversions\n"
        )
    return specs
