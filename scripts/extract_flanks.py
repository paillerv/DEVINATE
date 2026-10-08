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
    parser.add_argument("--min-flank", type=int, default=100)
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
        
        ref_len = bam.get_reference_length(read.reference_name)

        # -------------------------
        # 5' FLANK
        # -------------------------
        if not variant.startswith("5p_Truncated"):
            if cigar[0][0] == 4:  # Soft-clipping au début
                l = cigar[0][1]
                ref_start = read.reference_start  # Position de début sur la référence
                
                # On retire la portion qui empiète à l'intérieur du TE
                if ref_start < l:
                    true_flank_len = l - ref_start
                    if true_flank_len >= min_flank:
                        seq_flank = seq[ref_start:l]
                        flanks_per_read[read_id].append(
                            ("flank5", true_flank_len, seq_flank)
                        )

        # -------------------------
        # 3' FLANK
        # -------------------------
        if not variant.startswith("3p_Truncated"):
            if cigar[-1][0] == 4:  # Soft-clipping à la fin
                l = cigar[-1][1]
                ref_end = read.reference_end  # Fin de l'alignement sur la référence
                
                # Distance entre la fin de l'alignement et la fin du TE de référence
                remaining_ref = ref_len - ref_end
                
                # On retire la portion qui empiète sur la fin du TE
                if l > remaining_ref:
                    true_flank_len = l - remaining_ref
                    if true_flank_len >= min_flank:
                        seq_flank = seq[-true_flank_len:]
                        flanks_per_read[read_id].append(
                            ("flank3", true_flank_len, seq_flank)
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
