import sys
import subprocess
import re
from collections import defaultdict

# -------------------------
# PARAMETERS
# -------------------------
BED  = sys.argv[1]
BAM  = sys.argv[2]
CLSF = sys.argv[3]

TOLERANCE = 10000

# -------------------------
# Utils
# -------------------------

def get_ref_length(cigar):
    length = 0
    for l, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        if op in ('M', 'D', 'N', '=', 'X'):
            length += int(l)
    return length

def get_aln_type(flag):
    if flag & 2048:
        return "SUPPLEMENTARY"
    elif flag & 256:
        return "SECONDARY"
    else:
        return "PRIMARY"

# -------------------------
# 1. Classification
# -------------------------

classification = {}

with open(CLSF) as f:
    header = f.readline().strip().split('\t')
    col_idx = {name: i for i, name in enumerate(header)}
    
    for line in f:
        if not line.strip():
            continue
        fields = line.strip().split('\t')
        if len(fields) <= max(col_idx.values()):
            continue
        
        read_id = fields[col_idx['read_id']]
        variant = fields[col_idx['variant']]
        
        classification[read_id] = variant

# -------------------------
# 2. Flanks (MULTI)
# -------------------------

flanks = defaultdict(list)

with open(BED) as f:
    for line in f:
        if line.startswith("chr"):
            continue

        fields = line.strip().split("\t")

        if len(fields) == 8:
            chrom, start, end, read_id, variant, category, mode, dist = fields
            flank_status = "NA"
        elif len(fields) == 9:
            (
                chrom,
                start,
                end,
                read_id,
                variant,
                category,
                mode,
                flank_status,
                dist,
            ) = fields
        else:
            continue

        if "," in chrom:
            chroms = chrom.split(",")
            starts = start.split(",")
            ends = end.split(",")
            for c, s, e in zip(chroms, starts, ends):
                flank_pos = (int(s) + int(e)) // 2
                flanks[read_id].append(
                    {"chrom": c, "pos": flank_pos, "status": flank_status}
                )
        else:
            flank_pos = (int(start) + int(end)) // 2
            flanks[read_id].append(
                {
                    "chrom": chrom,
                    "pos": flank_pos,
                    "status": flank_status,
                }
            )

# -------------------------
# 3. Alignements (MULTI)
# -------------------------

alignments = defaultdict(list)

result = subprocess.run(
    ["samtools", "view", "-F", "4", "-q", "30", BAM],
    capture_output=True, text=True
)

for line in result.stdout.splitlines():
    f = line.split('\t')

    read_id = f[0]
    flag    = int(f[1])
    chrom   = f[2]
    pos     = int(f[3]) - 1
    mapq    = int(f[4])
    cigar   = f[5]

    if cigar == "*":
        continue

    ref_len = get_ref_length(cigar)
    end     = pos + ref_len

    alignments[read_id].append({
        "chrom": chrom,
        "start": pos,
        "end": end,
        "mapq": mapq,
        "type": get_aln_type(flag)
    })

# -------------------------
# 4. Comparaison ALL vs ALL
# -------------------------

print("read_id\tvariant\tflank_chr\tflank_pos\tread_chr\tread_start\tread_end\tmapq\taln_type\tmatch\tdistance")

for read_id in flanks:

    if read_id not in alignments:
        continue

    best = None
    best_dist = float("inf")
    best_score = (float("inf"), float("inf"))
    best_same_chr = False
    best_flank = None

    # ALL flanks vs ALL alignments
    for f in flanks[read_id]:
        for aln in alignments[read_id]:

            same_chr = (f["chrom"] == aln["chrom"])

            if same_chr:
                if f["pos"] < aln["start"]:
                    dist = aln["start"] - f["pos"]
                elif f["pos"] > aln["end"]:
                    dist = f["pos"] - aln["end"]
                else:
                    dist = 0
            else:
                dist = float("inf")

            # priority : 1. same chrom, 2. minimal distance
            score = (0 if same_chr else 1, dist)

            if best is None or score < best_score:
                best = aln
                best_score = score
                best_dist = dist if same_chr else -1
                best_same_chr = same_chr
                best_flank = f

    if best is None:
        continue

    # -------------------------
    # Classification
    # -------------------------

    if not best_same_chr:
        match = "DIFF_CHR"
    elif best_dist == 0:
        match = "MATCH"
    elif best_dist <= TOLERANCE:
        match = "NEAR"
    else:
        match = "FAR"

    variant_val = classification.get(read_id, "NA")

    print(
        f"{read_id}\t{variant_val}\t"
        f"{best_flank['chrom']}\t{best_flank['pos']}\t"
        f"{best['chrom']}\t{best['start']}\t{best['end']}\t"
        f"{best['mapq']}\t{best['type']}\t{match}\t{best_dist}"
    )
