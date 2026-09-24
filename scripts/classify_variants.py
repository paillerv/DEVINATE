#!/usr/bin/env python3

import os
import sys
from collections import defaultdict
import pysam

# -------------------------
# PARAMETERS
# -------------------------

MIN_ID = 0.80
MIN_COV = 0.50
MAX_INDEL_RATE = 0.10
MAX_UNEXPECTED_DEL = 200
MAX_UNEXPECTED_INS = 200

# Continuous scoring & ambiguity parameters
WINDOW = 100
SCORE_THRESH = 0.50
DELTA_AMBIGUOUS = 0.05

# Insertions
INSERTION_TOLERANCE = 50
CLUSTER_THRESHOLD=50

# -------------------------
# INPUTS
# -------------------------

if len(sys.argv) < 3:
    print("Usage: python classify_variants.py <bam> <variants.txt> [outdir]")
    sys.exit(1)

bam_path = sys.argv[1]
variants_file = sys.argv[2]
OUTDIR = sys.argv[3] if len(sys.argv) > 3 else "results"

# -------------------------
# BAM & REF LENGTH
# -------------------------

bam = pysam.AlignmentFile(bam_path, "rb")

REF_LEN = bam.lengths[0] if bam.lengths else None
if REF_LEN is None:
    print("[ERROR] Cannot detect TE length from BAM header")
    sys.exit(1)

print(f"[INFO] TE reference length : {REF_LEN} bp")

# -------------------------
# OUTPUT DIRECTORIES & FILES
# -------------------------

BAMDIR = os.path.join(OUTDIR, "bam")
os.makedirs(BAMDIR, exist_ok=True)

tsv = open(os.path.join(OUTDIR, "classification.tsv"), "w")
log = open(os.path.join(OUTDIR, "summary.log"), "w")

out_bams = {}
counts = defaultdict(int)

stats_acc = defaultdict(
    lambda: {
        "identity": [],
        "coverage": [],
        "indel_rate": [],
        "read_len": [],
        "ref_span": [],
    }
)

# -------------------------
# LOAD VARIANTS
# -------------------------

VARIANTS = {}

with open(variants_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        # Split the variants file
        parts = line.rsplit(":", 3)
        if len(parts) < 3:
            print(f"[WARNING] Wrong line : {line}")
            continue

        # Extract 
        if len(parts) == 4:
            te, vtype, coords, support_str = parts
            try:
                support = int(support_str)
            except ValueError:
                support = 0
        else:
            te, vtype, coords = parts[0], parts[1], parts[2]
            support = 0

        if vtype == "DEL":
            try:
                start, end = map(int, coords.split("-"))
                VARIANTS[line] = {
                    "type": "DEL",
                    "start": start,
                    "end": end,
                    "length": abs(end - start),
                    "support": support,
                }
            except Exception:
                print(f"[WARNING] Skipping malformed DEL: {line}")

        elif vtype == "INS":
            try:
                pos, length = coords.split("+")
                pos = int(pos)
                length = int(length.replace("bp", ""))
                VARIANTS[line] = {
                    "type": "INS",
                    "pos": pos,
                    "length": length,
                    "support": support,
                }
            except Exception:
                print(f"[WARNING] Skipping malformed INS: {line}")

# -------------------------
# HELPER FUNCTIONS: BAM WRITERS
# -------------------------

def sanitize_filename(name):
    return name.replace(":", "_").replace(",", "_")


def get_bam(name):
    if name not in out_bams:
        fname = f"{sanitize_filename(name)}.bam"
        path = os.path.join(BAMDIR, fname)
        out_bams[name] = pysam.AlignmentFile(path, "wb", template=bam)
    return out_bams[name]

# -------------------------
# SCORING & VARIANT EVALUATION
# -------------------------
def score_variant(ds, de, vs, ve, window=WINDOW):
    """Compute the score :overlap * start_match * end_match"""
    dlen = de - ds
    vlen = ve - vs
    ov = max(0, min(de, ve) - max(ds, vs))
    if ov == 0:
        return 0.0

    overlap_score = ov / max(dlen, vlen)
    start_score = max(0.0, 1.0 - (abs(ds - vs) / window))
    end_score = max(0.0, 1.0 - (abs(de - ve) / window))

    return overlap_score * start_score * end_score

def evaluate_read_dels(read_dels):
    total_matched_del_len = 0
    assigned_labels = []

    for ds, de in read_dels:
        candidates = []
        for var_id, info in VARIANTS.items():
            if info["type"] == "DEL":
                sc = score_variant(ds, de, info["start"], info["end"])
                if sc >= SCORE_THRESH:
                    candidates.append((var_id, sc, info["length"]))

        if not candidates:
            continue

        candidates.sort(key=lambda x: x[1], reverse=True)
        best_var, best_score, best_len = candidates[0]

        # Resolve ambiguities with the DELTA_AMBIGUOUS score 
        valid = [c for c in candidates if c[1] >= SCORE_THRESH]
        if len(valid) > 1 and (best_score - valid[1][1]) < DELTA_AMBIGUOUS:
            amb_ids = [
                c[0] for c in valid if (best_score - c[1]) < DELTA_AMBIGUOUS
            ]
            assigned_labels.append(f"ambiguous:{','.join(sorted(amb_ids))}")
        else:
            assigned_labels.append(best_var)

        total_matched_del_len += best_len

    return assigned_labels, total_matched_del_len


def detect_insertions(insertions):
    found = []
    total_matched_ins_len = 0
    for var_id, info in VARIANTS.items():
        if info["type"] != "INS":
            continue
        pos = info["pos"]
        size = info["length"]

        for ipos, ilen in insertions:
            if (
                abs(ipos - pos) <= INSERTION_TOLERANCE
                and abs(ilen - size) <= INSERTION_TOLERANCE
            ):
                found.append(var_id)
                total_matched_ins_len += size
                break
    return found, total_matched_ins_len

def classify_and_filter_read(
    pos_start: int,
    pos_end: int,
    te_len: int,
    softclip_left: int,
    softclip_right: int,
    mapq: int,
    identity: float,
    min_cov: float = 0.50,
    min_clip: int = 300,
    min_mapq: int = 30,
    min_id: float = 0.85,
    tolerance: int = 300
) -> tuple[str, str, str]:

    """Return (STATUT_DECISION, CATEGORIE_BIOLOGIQUE, RAISON_DROP)"""

    aligned_len = pos_end - pos_start + 1
    coverage = aligned_len / te_len
    max_clip = max(softclip_left, softclip_right)

    # --- 1. QUALITY FILTERS ---
    if coverage < min_cov:
        return "DROP", "N/A", "coverage_too_low"
    if max_clip < min_clip:
        return "DROP", "N/A", "softclip_too_short"
    if mapq < min_mapq:
        return "DROP", "N/A", "low_MAPQ"
    if identity < min_id:
        return "DROP", "N/A", "low_identity"

    # --- 2. CLASSIFICATION ---
    is_5p_complete = pos_start <= tolerance
    is_3p_complete = pos_end >= (te_len - tolerance)

    if is_5p_complete and is_3p_complete:
        category = "Full_Length"
    elif not is_5p_complete and is_3p_complete:
        category = "5p_Truncated"
    elif is_5p_complete and not is_3p_complete:
        category = "3p_Truncated"
    else:
        return "DROP", "Internal_Fragment", "internal_fragment"

    return "KEEP", category, "."

# -------------------------
# METRICS COMPUTATION
# -------------------------

def compute_identity(aln):
    """Identity based on 'de' minimap2's tag or 'NM'"""
    if aln.has_tag("de"):
        return max(0.0, 1.0 - aln.get_tag("de"))

    if aln.has_tag("NM"):
        nm = aln.get_tag("NM")
        aln_len = aln.query_alignment_length
        return max(0.0, 1.0 - (nm / aln_len)) if aln_len else 0.0

    return 0.0

# -------------------------
#  RENAME VARIANTS
# -------------------------

def simplify_variant_name(variant):
    """From Gypsy:ZAM-fl:DEL:5964-6591:72 to DEL:5964-6591 ; TIRANT:INS:5000+300bp:42 -> INS:5000+300b"""
    parts = variant.rsplit(":", 3)

    if len(parts) == 4:
        _, vtype, coords, support = parts

        if vtype in ("DEL", "INS"):
            return f"{vtype}:{coords}"

    return variant

# -------------------------
# HEADER TSV
# -------------------------

tsv.write(
    "read_id\tvariant\tcategory\tidentity\tcoverage\t"
    "indel_rate\tread_len\tref_start\tref_end\tref_span\t"
    "drop_reason\tdecision\n"
)

# -------------------------
# MAIN LOOP
# -------------------------

for aln in bam:
    if aln.is_unmapped or aln.cigarstring is None:
        continue
    if aln.is_secondary or aln.is_supplementary:
        continue

    read_id = aln.query_name

    # Extraction CIGAR & Soft-clips (op == 4)
    curr_pos = aln.reference_start
    read_dels = []
    insertions = []
    total_indel_bp = 0
    
    softclip_left = aln.cigartuples[0][1] if aln.cigartuples[0][0] == 4 else 0
    softclip_right = aln.cigartuples[-1][1] if aln.cigartuples[-1][0] == 4 else 0

    for op, length in aln.cigartuples:
        if op == 2:  # Deletion (D)
            read_dels.append((curr_pos, curr_pos + length))
            total_indel_bp += length
            curr_pos += length
        elif op == 1:  # Insertion (I)
            insertions.append((curr_pos, length))
            total_indel_bp += length
        elif op in (0, 7, 8):  # Match
            curr_pos += length
        elif op == 3:  # Skip
            curr_pos += length

    ref_span = aln.reference_end - aln.reference_start
    pos_start = aln.reference_start
    pos_end = aln.reference_end

    # Extraction of the coordinates 1-based
    ref_start = aln.reference_start + 1 if aln.reference_start is not None else "."
    ref_end = aln.reference_end if aln.reference_end is not None else "."

    # Evaluate of the variants
    del_labels, sum_matched_dels = evaluate_read_dels(read_dels)
    ins_labels, sum_matched_ins = detect_insertions(insertions)
    all_detected = del_labels + ins_labels

    # Global metrics (identity, coverage....)
    identity = compute_identity(aln)
    read_len = aln.infer_read_length() or 0
    coverage = ref_span / REF_LEN if REF_LEN else 0.0

    sum_matched_indels = sum_matched_dels + sum_matched_ins
    unexplained_indels = max(0, total_indel_bp - sum_matched_indels)
    aln_len = aln.query_alignment_length or 1
    indel_rate = unexplained_indels / aln_len

    # --- CLASSIFICATION + FILTERING ---
    status, category, initial_drop = classify_and_filter_read(
        pos_start=pos_start,
        pos_end=pos_end,
        te_len=REF_LEN,
        softclip_left=softclip_left,
        softclip_right=softclip_right,
        mapq=aln.mapping_quality,
        identity=identity,
        min_cov=MIN_COV,
        min_clip=15
    )

    variant_category = "Full_Length_Like" if category == "Full_Length" else category

    if not all_detected:
        base_variant = variant_category
    else:
        unique_vars = list(dict.fromkeys(
            simplify_variant_name(v) for v in all_detected
        ))
        vars_str = "+".join(unique_vars)
        base_variant = f"{variant_category}:{vars_str}"

    output_category = "FLL" if category == "Full_Length" else category

    drop = []
    if initial_drop != ".":
        drop.append(initial_drop)

    if indel_rate > MAX_INDEL_RATE:
        drop.append("high_indel_rate")

    # Check unexpected deletions
    for ds, de in read_dels:
        dlen = de - ds
        if dlen >= MAX_UNEXPECTED_DEL:
            is_known = any(
                info["type"] == "DEL" and score_variant(ds, de, info["start"], info["end"]) >= SCORE_THRESH
                for info in VARIANTS.values()
            )
            if not is_known:
                drop.append("unexpected_deletion")
                break

    # Check unexpected insertions
    for ins_pos, ins_len in insertions:
        if ins_len >= MAX_UNEXPECTED_INS:
            is_known_ins = any(
                info["type"] == "INS" and abs(info["pos"] - ins_pos) <= CLUSTER_THRESHOLD 
                for info in VARIANTS.values()
            )
            if not is_known_ins:
                drop.append("unexpected_insertion")
                break

    decision = "KEEP" if not drop else "DROP"

    # Writing TSV output
    tsv.write(
        f"{read_id}\t"
        f"{base_variant}\t"
        f"{output_category}\t"
        f"{identity:.3f}\t"
        f"{coverage:.3f}\t"
        f"{indel_rate:.3f}\t"
        f"{read_len}\t"
        f"{ref_start}\t"
        f"{ref_end}\t"
        f"{ref_span}\t"
        f"{','.join(drop) if drop else '.'}\t"
        f"{decision}\n"
    )

    counts[(base_variant, decision)] += 1

    if decision == "KEEP":
        get_bam(base_variant).write(aln)

# -------------------------
# FERMETURE ET INDEXATION
# -------------------------

bam.close()
tsv.close()

print("[INFO] BAM indexing...")

for name, f in out_bams.items():
    path = os.path.join(BAMDIR, f"{sanitize_filename(name)}.bam")
    f.close()
    try:
        pysam.index(path)
    except Exception as e:
        print(f"[WARNING] Indexation fail for {path}: {e}")

# -------------------------
# GÉNÉRATION DU SUMMARY LOG
# -------------------------

log.write("=== SUMMARY COUNTS ===\n\n")
log.write(f"{'variant':<55} {'decision':<10} {'n':>8}\n")
log.write("-" * 75 + "\n")

total = 0
for (variant, decision), n in sorted(counts.items()):
    log.write(f"{variant:<55} {decision:<10} {n:>8}\n")
    total += n

log.write("-" * 75 + "\n")
log.write(f"{'TOTAL':<55} {'':10} {total:>8}\n\n")

log.write("=== STATS KEEP READS ONLY ===\n")


def summary(values):
    if not values:
        return "n=0"
    values = sorted(values)
    med = values[len(values) // 2]
    return f"n={len(values)} min={min(values):.3f} median={med:.3f} max={max(values):.3f}"


for variant in sorted(stats_acc):
    acc = stats_acc[variant]
    log.write(f"\n--- {variant} ---\n")
    log.write(f"identity   : {summary(acc['identity'])}\n")
    log.write(f"coverage   : {summary(acc['coverage'])}\n")
    log.write(f"indel_rate : {summary(acc['indel_rate'])}\n")
    log.write(f"read_len   : {summary(acc['read_len'])}\n")
    log.write(f"ref_span   : {summary(acc['ref_span'])}\n")

log.write("\n\n=== DROP REASONS ===\n\n")
drop_counts = defaultdict(int)

with open(os.path.join(OUTDIR, "classification.tsv")) as f:
    f.readline()  # Skip header
    for line in f:
        fields = line.rstrip().split("\t")
        if fields[-1] == "DROP":
            for reason in fields[-2].split(","):
                drop_counts[reason] += 1

for reason, n in sorted(drop_counts.items(), key=lambda x: -x[1]):
    log.write(f"{reason:<30} {n:>8}\n")

log.close()

# -------------------------
# RAPPORT FINAL TERMINAL
# -------------------------

print("\n==============================")
print(" classify_variants terminé ")
print("==============================")
print(f"Variants generated (KEEP) : {len(out_bams)}")
print(f"Output repertory             : {OUTDIR}")
print("==============================")
