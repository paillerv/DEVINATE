#!/usr/bin/env python3
import sys
import os
import pysam

# Discovery parameters
CLUSTER_THRESHOLD = 50
MIN_VAR_LENGTH    = 150     # DEL or INS
MIN_AF            = 0.05
END_MARGIN        = 200

if len(sys.argv) < 3:
    print("Usage: python3 discover_variants_only.py <input.bam> <output_variants.txt>")
    sys.exit(1)

input_bam_path = sys.argv[1]
out_file = sys.argv[2]

bam_in = pysam.AlignmentFile(input_bam_path, "rb")

if not bam_in.lengths:
    print("[ERROR] Cannot read the BAM header.")
    sys.exit(1)

te_length = bam_in.lengths[0]
te_name = bam_in.references[0]

clusters = []
total_full_length_reads = 0

for read in bam_in.fetch():
    if read.is_unmapped or not read.cigartuples:
        continue
        
    # Extract only full spanning reads on the TE (5' to 3')
    if (read.reference_start <= END_MARGIN) and (read.reference_end >= (te_length - END_MARGIN)):
        total_full_length_reads += 1
        current_pos = read.reference_start
        
        for operation, length in read.cigartuples:
            var_type = None
            
            # DEL (CIGAR op 2)
            if operation == 2 and length >= MIN_VAR_LENGTH:
                var_type = "DEL"
                v_start = current_pos
                v_end = current_pos + length
                
            # INS (CIGAR op 1)
            elif operation == 1 and length >= MIN_VAR_LENGTH:
                var_type = "INS"
                v_start = current_pos
                v_end = current_pos  # An insertion doesnt increment the current position

            if var_type:
                matched_cluster = False
                for c in clusters:
                    if c["type"] != var_type:
                        continue
                    
                    mean_c_start = c["starts_sum"] / c["count"]
                    mean_c_end = c["ends_sum"] / c["count"]

                    # For DEL : check the overlap start/end
                    # For INS : v_start == v_end, check the proximity for an insertion site
                    # Check is the variant is in the tolerance area (CLUSTER_THRESHOLD = 50)
                    if abs(v_start - mean_c_start) <= CLUSTER_THRESHOLD and abs(v_end - mean_c_end) <= CLUSTER_THRESHOLD:
                        c["starts_sum"] += v_start
                        c["ends_sum"] += v_end
                        c["lengths_sum"] += length
                        c["count"] += 1
                        matched_cluster = True
                        break
                
                # If no cluster matches, a new one is created
                if not matched_cluster:
                    clusters.append({
                        "type": var_type,
                        "starts_sum": v_start,
                        "ends_sum": v_end,
                        "lengths_sum": length,
                        "count": 1
                    })

            # Advance of the position relative to the reference  (DEL, MATCH, EQUAL, DIFF)
            if operation in (0, 2, 7, 8):
                current_pos += length

bam_in.close()

if total_full_length_reads == 0:
    print("[WARNING] No overlapping read found.")
    min_reads_required = 0
else:
    min_reads_required = max(1, int(total_full_length_reads * MIN_AF))

print(f"[INFO] Overlapping reads found : {total_full_length_reads}")
print(f"[INFO] Required threshold ({MIN_AF*100}%): {min_reads_required} read(s) required to consider a structural variant.")

# Export the variants with the number of supporting reads
exported_count = 0

with open(out_file, "w") as f:

    for c in sorted(
        clusters,
        key=lambda x: x["count"],
        reverse=True
    ):

        if c["count"] >= min_reads_required:

            mean_start = int(
                c["starts_sum"] / c["count"]
            )

            mean_end = int(
                c["ends_sum"] / c["count"]
            )

            mean_len = int(
                c["lengths_sum"] / c["count"]
            )


            if c["type"] == "DEL":

                var_id = (
                    f"{te_name}:DEL:"
                    f"{mean_start}-{mean_end}"
                )

            else:

                var_id = (
                    f"{te_name}:INS:"
                    f"{mean_start}+{mean_len}bp"
                )


            # number of reads supporting the variant
            support = c["count"]


            f.write(
                f"{var_id}:{support}\n"
            )


            exported_count += 1

print(f"[INFO] Variant detection finished. {exported_count} variant(s) exported in {out_file}")
