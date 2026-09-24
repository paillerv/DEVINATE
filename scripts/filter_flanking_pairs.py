import sys
import subprocess
import re
from collections import defaultdict

BAM = sys.argv[1]

# -------------------------
# PARAMETERS
# -------------------------
MAPQ_MIN        = 40
MAX_DIST        = 20000

MIN_ALIGNED_LEN = 150
MIN_COVERAGE    = 0.7
MAX_MISMATCH    = 0.10

CANONICAL_CHROMS = {'2L', '2R', '3L', '3R', 'X', '4'}

# -------------------------
def get_ref_length(cigar):
    length = 0
    for l, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):
        if op in ('M', 'D', 'N', '=', 'X'):
            length += int(l)
    return length

def get_nm(fields):
    for field in fields:
        if field.startswith("NM:i:"):
            return int(field.split(":")[2])
    return None

# -------------------------
reads = defaultdict(list)

result = subprocess.run(
    ["samtools", "view", "-F", "4", "-F", "256", "-F", "2048", "-q", str(MAPQ_MIN), BAM],
    capture_output=True, text=True
)

for line in result.stdout.splitlines():
    f = line.strip().split('\t')

    name   = f[0]
    flag   = int(f[1])
    chrom  = f[2]
    pos    = int(f[3]) - 1
    mapq   = int(f[4])
    cigar  = f[5]

    if cigar == '*':
        continue
    if chrom not in CANONICAL_CHROMS:
        continue

    parts = name.split('|')

    # -------------------------
    # PARSE HEADER
    # -------------------------
    read_id = parts[0]
    side    = parts[1] if len(parts) > 1 else None

    # if variant present 
    if len(parts) >= 4:
        variant = parts[3]
    else:
        variant = "TE"

    match = re.search(r'flank[35]_(\d+)bp', side)
    if not match:
        continue

    flank_len = int(match.group(1))
    ref_len   = get_ref_length(cigar)

    # -------------------------
    # FILTERS
    # -------------------------
    if ref_len < MIN_ALIGNED_LEN:
        continue

    coverage = ref_len / flank_len if flank_len > 0 else 0
    if coverage < MIN_COVERAGE:
        continue

    nm = get_nm(f[11:])
    if nm is not None and ref_len > 0:
        if (nm / ref_len) > MAX_MISMATCH:
            continue

    # -------------------------
    # breakpoint
    # -------------------------
    end        = pos + ref_len
    is_reverse = (flag & 16) != 0

    if side.startswith('flank5'):
        bp = end if not is_reverse else pos
    else:
        bp = pos if not is_reverse else end

    reads[read_id].append({
        'chrom'    : chrom,
        'bp'       : bp,
        'side'     : side,
        'mapq'     : mapq,
        'variant'  : variant
    })

# -------------------------
# OUTPUT
# -------------------------
print("chr\tstart\tend\tread_id\tvariant\tsupport\tflank_status\tdist")

for read_id, alns in reads.items():

    f5 = [a for a in alns if a["side"].startswith("flank5")]
    f3 = [a for a in alns if a["side"].startswith("flank3")]

    variant = alns[0]["variant"]

    flank_status = "2flanks" if (f5 and f3) else "1flank"

    # ── paired ─────────────────────────────
    if f5 and f3:
        best_pair = None
        best_score = -1

        for a5 in f5:
            for a3 in f3:

                if a5["chrom"] != a3["chrom"]:
                    continue

                dist = abs(a5["bp"] - a3["bp"])

                if dist > MAX_DIST:
                    continue

                score = a5["mapq"] + a3["mapq"]

                if score > best_score:
                    best_score = score
                    best_pair = (a5, a3)

        if best_pair:
            a5, a3 = best_pair
            start = min(a5["bp"], a3["bp"])
            end_bp = max(a5["bp"], a3["bp"])
            dist = abs(a5["bp"] - a3["bp"])

            print(
                f"{a5['chrom']}\t{start}\t{end_bp}\t{read_id}\t"
                f"{variant}\tpaired\t{flank_status}\t{dist}"
            )
            continue

        # fallback: 2 flanks found but unpaired
        best_f5 = max(f5, key=lambda x: x["mapq"])
        best_f3 = max(f3, key=lambda x: x["mapq"])

        chroms = f"{best_f5['chrom']},{best_f3['chrom']}"
        starts = f"{best_f5['bp']},{best_f3['bp']}"
        ends = f"{best_f5['bp']+1},{best_f3['bp']+1}"

        print(
            f"{chroms}\t{starts}\t{ends}\t{read_id}\t"
            f"{variant}\tunpaired\t{flank_status}\tNA"
        )
        continue

    # a single flank
    best = max(alns, key=lambda x: x["mapq"])
    print(
        f"{best['chrom']}\t{best['bp']}\t{best['bp']+1}\t{read_id}\t"
        f"{variant}\tsingle\t{flank_status}\tNA"
    )
