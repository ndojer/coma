"""
A script accepting command line arguments.
"""
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        prog="generate_sv.py",
        description=(
            "Generate a mutated FASTA containing structural variants: "
            "inversions, tandem duplications, intrachromosomal translocations, inverted intrachromosomal translocations, "
            "interchromosomal translocations, inverted interchromosomal translocations."
        ),
    )

    parser.add_argument("fasta", help="path to the FASTA file with the reference sequence")
    parser.add_argument("--output", required=True, help="output FASTA file (modified sequence)")
    parser.add_argument("--csv", required=True, help="output CSV file (list of SVs)")
    parser.add_argument("--min_distance", type=int, required=True, help="minimum distance between SV breakpoints (nt)")

    parser.add_argument("--inv_count", type=int, help="number of inversions")
    parser.add_argument("--dup_count", type=int, help="number of tandem duplications")
    parser.add_argument("--trl_count", type=int, help="number of intrachromosomal translocations")
    parser.add_argument("--trl_interchr", type=int, default=0, help="number of interchromosomal translocations")
    parser.add_argument("--invtrl_count", type=int, help="number of intrachromosomal translocations with inversion")
    parser.add_argument("--trl_interchr_inv", type=int, default=0, help="number of interchromosomal translocations with inversion")

    parser.add_argument("--inv_min_length", type=int, help="minimum length of inversion")
    parser.add_argument("--inv_max_length", type=int, help="maximum length of inversion")

    parser.add_argument("--dup_min_length", type=int, help="minimum length of tandem duplication")
    parser.add_argument("--dup_max_length", type=int, help="maximum length of tandem duplication")

    parser.add_argument("--trl_min_length", type=int, help="minimum length of translocation segment")
    parser.add_argument("--trl_max_length", type=int, help="maximum length of translocation segment")

    parser.add_argument("--invtrl_min_length", type=int, help="minimum length of translocation with inversion")
    parser.add_argument("--invtrl_max_length", type=int, help="maximum length of translocation with inversion")

    parser.add_argument("--ins_count", type=int, help="number of insertions")
    parser.add_argument("--ins_min_length", type=int, help="minimum length of insertion")
    parser.add_argument("--ins_max_length", type=int, help="maximum length of insertion")

    parser.add_argument("--del_count", type=int, help="number of deletions")
    parser.add_argument("--del_min_length", type=int, help="minimum length of deletion")
    parser.add_argument("--del_max_length", type=int, help="maximum length of deletion")

    return parser.parse_args()
