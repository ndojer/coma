"""
A script containing the classes used in the pipeline.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Optional

__all__ = ["SVType","SVSpec","MutationEvent"]

class SVType(str, Enum):

    INV = "INVERSION"
    DUP_TANDEM = "DUPLICATION"

    TRANS_INTRACHR = "TRANSLOCATION"
    TRANS_INTERCHR = "TRANSLOCATION_INTERCHR"

    TRANS_INTRA_INV = "INV_TRANSLOCATION"
    TRANS_INTER_INV = "TRANSLOCATION_INTERCHR_INV"

    DEL = "DELETION"
    INS = "INSERTION"

    def __str__(self) -> str:
        return self.value

@dataclass(slots=True)
class SVSpec:
    sv_id: int
    sv_type: SVType

    chr_ref: str
    chr_mut: str

    start_ref: int
    end_ref: int

    insert_pos_ref: Optional[int] = None
    insert_seq: Optional[str] = None

    start_mut: Optional[int] = None
    end_mut: Optional[int] = None

    extra_pos_mut: Optional[int] = None

    orientation: str = "+"

    def __repr__(self) -> str:
        return (
            f"SVSpec(id={self.sv_id}, type={self.sv_type}, "
            f"src={self.chr_ref}:{self.start_ref}-{self.end_ref}, "
            f"dst={self.chr_mut}:{self.insert_pos_ref}, "
            f"mut={self.start_mut}-{self.end_mut}, orient={self.orientation})"
        )

@dataclass(slots=True)
class MutationEvent:
    sv_id: int
    chr: str
    start_ref: int
    end_ref: int
    seq_ins: str = ""

    @property
    def is_deletion(self) -> bool:
        return self.seq_ins == ""

    @property
    def len_ref(self) -> int:
        return self.end_ref - self.start_ref

    @property
    def len_ins(self) -> int:
        return len(self.seq_ins)

    def __lt__(self, other: "MutationEvent") -> bool:
        if self.chr != other.chr:
            return self.chr < other.chr
        return self.start_ref < other.start_ref
