import sys
from collections import defaultdict

# -------------------------
# PARAMETERS
# -------------------------
MAX_CLUSTER_DIST = 10000

insertions = []

# -------------------------
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith("chr"):
            continue

        fields = line.strip().split("\t")
        if len(fields) < 8:
            continue

        chrom, start, end, read_id, variant, support, flank_status, dist = fields

        if "," in chrom:
            chrom = chrom.split(",")[0]
        if "," in start:
            start = start.split(",")[0]
        if "," in end:
            end = end.split(",")[0]

        insertions.append(
            {
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "read_id": read_id,
                "variant": variant,
                "support": support,
                "flank_status": flank_status,
                "dist": int(dist) if dist != "NA" else 0,
            }
        )

# -------------------------
# Regroup
# -------------------------

by_chrom = defaultdict(list)
for ins in insertions:
    by_chrom[ins['chrom']].append(ins)

clusters = []

for chrom, ins_list in by_chrom.items():

    ins_list.sort(key=lambda x: x['start'])

    current_cluster = []
    cluster_start = None
    cluster_end = None

    for ins in ins_list:

        if not current_cluster:
            current_cluster = [ins]
            cluster_start = ins['start']
            cluster_end = ins['end']
            continue

        if ins['start'] <= cluster_end + MAX_CLUSTER_DIST:
            current_cluster.append(ins)
            cluster_end = max(cluster_end, ins['end'])
        else:
            clusters.append({
                'chrom': chrom,
                'start': cluster_start,
                'end': cluster_end,
                'reads': current_cluster
            })

            current_cluster = [ins]
            cluster_start = ins['start']
            cluster_end = ins['end']

    if current_cluster:
        clusters.append({
            'chrom': chrom,
            'start': cluster_start,
            'end': cluster_end,
            'reads': current_cluster
        })

# -------------------------
# Summary
# -------------------------

final_clusters = []

for cl in clusters:

    reads = cl['reads']

    variants     = sorted(set(x['variant'] for x in reads))
    supports     = sorted(set(x['support'] for x in reads))
    flank_stats = sorted(set(x['flank_status'] for x in reads))

    count = len(reads)
    paired_count = sum(1 for x in reads if x['support'] == 'paired')

    if paired_count >= 1:
        confidence = "HIGH"
    elif 'unpaired' in supports:
        confidence = "AMB"
    else:
        confidence = "LOW"

    final_clusters.append({
        'chrom': cl['chrom'],
        'start': cl['start'],
        'end': cl['end'],
        'count': count,
        'paired_count': paired_count,
        'confidence': confidence,
        'variants': variants,
        'supports': supports,
        'flank_status': flank_stats
    })


final_clusters.sort(key=lambda x: (x['chrom'], x['start']))

# -------------------------
# Output
# -------------------------

print("chr\tstart\tend\tcount\tpaired_count\tconfidence\tvariants\tsupports\tflank_status")

for cl in final_clusters:
    print(
        f"{cl['chrom']}\t{cl['start']}\t{cl['end']}\t"
        f"{cl['count']}\t{cl['paired_count']}\t{cl['confidence']}\t"
        f"{';'.join(cl['variants'])}\t"
        f"{';'.join(cl['supports'])}\t"
        f"{';'.join(cl['flank_status'])}"
    )
