import pysam
import sys
import argparse
from collections import defaultdict

# -------------------------
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam")
    parser.add_argument("classification")
    parser.add_argument("output")
    parser.add_argument("--min-flank", type=int, default=200)
    return parser.parse_args()

# -------------------------
def load_classifications(tsv):
    d = {}

    with open(tsv) as f:
        header = f.readline().strip().split('\t')
        col = {name: i for i, name in enumerate(header)}

        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < len(header):
                continue

            read_id = fields[col["read_id"]]
            decision = fields[col["decision"]]

            if decision != "KEEP":
                continue

            variant = fields[col["variant"]]
            d[read_id] = variant

    return d

# -------------------------
def extract_flanks(bam_path, classifications, min_flank, output):
    bam = pysam.AlignmentFile(bam_path, "rb")

    flanks_per_read = defaultdict(list)

    for read in bam.fetch():
        if read.query_name not in classifications:
            continue
        if read.cigartuples is None:
            continue

        read_id = read.query_name
        variant = classifications[read_id]
        seq     = read.query_sequence
        cigar   = read.cigartuples

        # 5' flank
        if not variant.startswith("5p_Truncated"):
            if cigar[0][0] == 4:
                l = cigar[0][1]
                if l >= min_flank:
                    flanks_per_read[read_id].append(
                        ("flank5", l, seq[:l])
                    )

        # 3' flank
        if not variant.startswith("3p_Truncated"):
            if cigar[-1][0] == 4:
                l = cigar[-1][1]
                if l >= min_flank:
                    flanks_per_read[read_id].append(
                        ("flank3", l, seq[-l:])
                    )

    bam.close()

    # -------------------------
    # WRITE OUTPUT
    # -------------------------
    n_written = 0

    with open(output, "w") as out:
        for read_id, flanks in flanks_per_read.items():
            if not flanks:
                continue
                
            variant = classifications[read_id]

            sides = {f[0] for f in flanks}
            
            # Categorization according to the flanks
            if variant.startswith("Full_Length_Like"):
                category = "complete" if {"flank5", "flank3"} <= sides else "partial"
            else:
                category = "partial" 

            for side, l, seq in flanks:
                out.write(
                    f">{read_id}|{side}_{l}bp|{category}|{variant}\n{seq}\n"
                )
            n_written += 1 

    print(f"[done] {n_written} reads avec flancs écrits")

# -------------------------
def main():
    args = parse_args()

    classifications = load_classifications(args.classification)
    print(f"[info] {len(classifications)} reads retained")

    extract_flanks(
        args.bam,
        classifications,
        args.min_flank,
        args.output,
    )

# -------------------------
if __name__ == "__main__":
    main()
