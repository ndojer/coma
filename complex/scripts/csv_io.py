"""
A script saving the generated SV in CSV format.
"""
import csv
from typing import List

from sv_models import SVSpec

def none2dash(v):
    return '-' if v is None else v

def write_csv(path: str, specs: List[SVSpec]) -> None:
    fields = [
        "SV_id", "chr_ref", "chr_mut", "type",
        "start_ref", "end_ref", "insert_ref", 
        "start_mut", "end_mut", "del/mid_mut",
        "orientation",
    ]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for s in specs:
            w.writerow({
                "SV_id": s.sv_id,
                "chr_ref": s.chr_ref,
                "chr_mut": s.chr_mut,
                "type": s.sv_type.value,
                "start_ref": s.start_ref,
                "end_ref": s.end_ref,
                "insert_ref": none2dash(s.insert_pos_ref),
                "start_mut": s.start_mut,
                "end_mut": s.end_mut,
                "del/mid_mut": none2dash(s.extra_pos_mut),
                "orientation": s.orientation,
            })
