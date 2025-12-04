import argparse
import pandas as pd
import csv

def load_sv_regions(csv_path):
    df = pd.read_csv(csv_path)
    sv_regions = []
    for _, row in df.iterrows():
        sv_regions.append({
            "chr": str(row["chr_ref"]).lower().replace("chr", ""),
            "start": int(row["start_ref"]),
            "end": int(row["end_ref"])
        })
    return sv_regions

def split_sdata_by_overlap(sdata_path, csv_path, overlap_out, non_overlap_out, connection_out):
    with open(sdata_path, "r", encoding="utf-8") as f:
        sdata_lines = f.readlines()

    header_lines = []
    fragment_lines = []
    for line in sdata_lines:
        if line.startswith("#") or "Fragment ID" in line:
            header_lines.append(line)
        else:
            fragment_lines.append(line)

    sv_regions = load_sv_regions(csv_path)

    overlapping = []
    non_overlapping = []
    connections = []

    for line in fragment_lines:
        parts = line.strip().split("\t")
        if len(parts) < 5:
            continue
        frag_id = parts[0]
        frag_chr = parts[1].lower().replace("chr", "")
        frag_start = int(parts[3])
        frag_end = int(parts[4])

        matched_any = False
        for sv in sv_regions:
            if sv["chr"] == frag_chr and not (frag_end < sv["start"] or frag_start > sv["end"]):
                if not matched_any:
                    overlapping.append(line)
                    matched_any = True
                connections.append({
                    "FragmentID": frag_id,
                    "Chr": f"chr{frag_chr}",
                    "Frag_Start": frag_start,
                    "Frag_End": frag_end,
                    "SV_Chr": f"chr{sv['chr']}",
                    "SV_Start": sv["start"],
                    "SV_End": sv["end"]
                })

        if not matched_any:
            non_overlapping.append(line)

    with open(overlap_out, "w", encoding="utf-8") as f:
        f.writelines(header_lines + overlapping)

    with open(non_overlap_out, "w", encoding="utf-8") as f:
        f.writelines(header_lines + non_overlapping)

    with open(connection_out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["FragmentID", "Chr", "Frag_Start", "Frag_End", "SV_Chr", "SV_Start", "SV_End"])
        writer.writeheader()
        writer.writerows(connections)

    print(f"Saved {len(overlapping)} overlapping fragments: {overlap_out}")
    print(f"Saved {len(non_overlapping)} non-overlapping fragments: {non_overlap_out}")
    print(f"Saved {len(connections)} SV connections: {connection_out}")

def main():
    parser = argparse.ArgumentParser(description="Split .sdata file by SV overlap and generate SV connection report")
    parser.add_argument("sdata", help="Path to .sdata file")
    parser.add_argument("csv", help="Path to .csv file with structural variants")
    parser.add_argument("--out_overlap", default="query_overlap.sdata", help="Output file for overlapping (default: query_overlap.sdata)")
    parser.add_argument("--out_nonoverlap", default="query_non_overlap.sdata", help="Output file for non-overlapping (default: query_non_overlap.sdata)")
    parser.add_argument("--connection_out", default="query_connection.csv", help="Output CSV file for SV-fragment connections")

    args = parser.parse_args()

    split_sdata_by_overlap(args.sdata, args.csv, args.out_overlap, args.out_nonoverlap, args.connection_out)

if __name__ == "__main__":
    main()
