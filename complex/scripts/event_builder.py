"""
A script converting the generated SVs into insertion/deletion notation and preparing the appropriate lists responsible for the sequence modification process.
"""
from __future__ import annotations

from typing import Dict, List, Tuple

from sv_models import MutationEvent, SVSpec, SVType

__all__ = ["build_events"]

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(s: str) -> str:
    """
    Return reverse-complement of a DNA sequence.
    """
    return s.translate(_COMP)[::-1]


def build_events(sv_specs: List[SVSpec], seqs: List[Tuple[str, str]]) -> List[MutationEvent]:
    """
    Convert SV list into a sorted list of insertions and deletions.
    """
    seq_dict: Dict[str, str] = {hdr: seq for hdr, seq in seqs}

    events: List[MutationEvent] = []

    for sv in sv_specs:
        if sv.sv_type == SVType.INV:
            _add_inversion_events(sv, seq_dict, events)

        elif sv.sv_type == SVType.DUP_TANDEM:
            _add_duplication_events(sv, seq_dict, events)

        elif sv.sv_type in {
            SVType.TRANS_INTRACHR,
            SVType.TRANS_INTERCHR,
            SVType.TRANS_INTRA_INV,
            SVType.TRANS_INTER_INV,
        }:
            _add_translocation_events(sv, seq_dict, events)

        elif sv.sv_type == SVType.DEL:
            events.append(
                MutationEvent(
                    sv_id=sv.sv_id,
                    chr=sv.chr_ref,
                    start_ref=sv.start_ref,
                    end_ref=sv.end_ref,
                    seq_ins="",
                )
            )

        elif sv.sv_type == SVType.INS:
            if sv.insert_pos_ref is None or not sv.insert_seq:
                raise ValueError("INS requires insert_pos_ref and insert_seq")
            events.append(
                MutationEvent(
                    sv_id=sv.sv_id,
                    chr=sv.chr_mut,
                    start_ref=sv.insert_pos_ref,
                    end_ref=sv.insert_pos_ref,
                    seq_ins=sv.insert_seq,
                )
            )

        else:
            raise NotImplementedError(
                f"build_events: unsupported SVType {sv.sv_type} – generator created an unknown variant"
            )

    events.sort()
    return events


def _add_inversion_events(sv: SVSpec, genome: Dict[str, str], out: List[MutationEvent]):
    """
    Decompose an inversion into deletion and insertion with a reverse-complemented fragment.
    """
    segment_rc = _revcomp(genome[sv.chr_ref][sv.start_ref: sv.end_ref])

    out.append(
        MutationEvent(
            sv_id=sv.sv_id,
            chr=sv.chr_ref,
            start_ref=sv.start_ref,
            end_ref=sv.end_ref,
            seq_ins="",
        )
    )

    out.append(
        MutationEvent(
            sv_id=sv.sv_id,
            chr=sv.chr_ref,
            start_ref=sv.start_ref,
            end_ref=sv.start_ref,
            seq_ins=segment_rc,
        )
    )


def _add_duplication_events(sv: SVSpec, genome: Dict[str, str], out: List[MutationEvent]):
    """
    Decompose a tandem duplication into an insertion at the duplication boundary.
    """
    segment = genome[sv.chr_ref][sv.start_ref: sv.end_ref]

    out.append(
        MutationEvent(
            sv_id=sv.sv_id,
            chr=sv.chr_ref,
            start_ref=sv.end_ref,
            end_ref=sv.end_ref,
            seq_ins=segment,
        )
    )


def _add_translocation_events(sv: SVSpec, genome: Dict[str, str], out: List[MutationEvent]):
    """
    Decompose a translocation into a deletion at the source and an insertion at the target.
    """
    if sv.insert_pos_ref is None:
        raise ValueError("SVSpec.insert_pos_ref is required for translocations")

    segment = genome[sv.chr_ref][sv.start_ref: sv.end_ref]
    if sv.orientation == "-":
        segment = _revcomp(segment)

    out.append(
        MutationEvent(
            sv_id=sv.sv_id,
            chr=sv.chr_ref,
            start_ref=sv.start_ref,
            end_ref=sv.end_ref,
            seq_ins="",
        )
    )

    out.append(
        MutationEvent(
            sv_id=sv.sv_id,
            chr=sv.chr_mut,
            start_ref=sv.insert_pos_ref,
            end_ref=sv.insert_pos_ref,
            seq_ins=segment,
        )
    )
