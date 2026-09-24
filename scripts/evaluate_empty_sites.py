import os
import sys
import pysam
from collections import defaultdict, Counter

# ------------------------------------------------
# PARAMETERS
# ------------------------------------------------

MIN_MAPQ = 20
SITE_FLANK = 50
SOFTCLIP_PROXIMITY = 15
INDEL_PROXIMITY = 50
MAX_LOCAL_INDEL = 100
LOCAL_WINDOW = 100

# ------------------------------------------------ 
# FUNCTIONS 
# ------------------------------------------------

def load_classification(classification_tsv):
    valid_te_reads = {}

    valid_categories = {
        "Full_Length",
        "Full_Length_Like",
        "FLL",
        "5p_Truncated",
        "3p_Truncated"
    }

    with open(classification_tsv, "r") as f:
        header = f.readline().strip().split("\t")
        col_idx = {name: i for i, name in enumerate(header)}

        required = {"read_id", "category", "decision"}
        missing = required - set(col_idx)

        if missing:
            raise ValueError(
                f"Colonness classification.tsv : "
                f"{', '.join(sorted(missing))}"
            )

        for line in f:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(col_idx.values()):
                continue

            read_id = fields[col_idx["read_id"]]
            category = fields[col_idx["category"]]
            decision = fields[col_idx["decision"]]

            if decision == "KEEP":
                is_valid = any(
                    category == cat or category.startswith(cat + ":")
                    for cat in valid_categories
                )
                if is_valid:
                    valid_te_reads[read_id] = category

    return valid_te_reads

def load_clusters(clusters_tsv):
    clusters = []
    with open(clusters_tsv, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if (
                parts[0].lower() == "chr"
                or parts[1].lower() == "start"
            ):
                continue

            try:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue

            clusters.append({
                "chrom": chrom,
                "start": start,
                "end": end
            })

    return clusters


def load_insertions(insertions_bed, valid_te_reads):
    insertions = []
    with open(insertions_bed, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")

            if len(parts) < 4:
                continue
            if (
                parts[0].lower() in ("chr", "chrom")
                or parts[1].lower() == "start"
            ):
                continue

            chrom = parts[0]
            read_id = parts[3]

            # Only TE+ validated
            if read_id not in valid_te_reads:
                continue

            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue

            variant = parts[4] if len(parts) > 4 else "NA"
            support = parts[6] if len(parts) > 6 else "NA"
            flank_status = parts[7] if len(parts) > 7 else "NA"

            insertions.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "read_id": read_id,
                "variant": variant,
                "support": support,
                "flank_status": flank_status
            })

    return insertions


def assign_insertions_to_clusters(insertions, clusters):
    cluster_insertions = defaultdict(list)
    for ins in insertions:
        candidates = []
        for idx, cl in enumerate(clusters):
            if cl["chrom"] != ins["chrom"]:
                continue
            if (
                cl["start"] <= ins["start"] <= cl["end"]
            ):
                candidates.append(idx)

        if not candidates:
            continue

        if len(candidates) == 1:
            best_idx = candidates[0]

        else:
            best_idx = min(
                candidates,
                key=lambda i: abs(
                    clusters[i]["start"] - ins["start"]
                )
            )
        cluster_insertions[best_idx].append(ins)
    return cluster_insertions


def determine_representative_site(insertions):
    if not insertions:
        return None, None, Counter()
    position_counts = Counter(
        ins["start"]
        for ins in insertions
    )
    max_count = max(position_counts.values())

    candidate_positions = [
        pos
        for pos, count in position_counts.items()
        if count == max_count
    ]

    if len(candidate_positions) == 1:
        representative_start = candidate_positions[0]
    else:
        support_scores = defaultdict(int)
        for ins in insertions:
            if ins["start"] not in candidate_positions:
                continue
            score = 0
            if ins["support"] == "paired":
                score += 2
            if ins["flank_status"] == "2flanks":
                score += 1
            support_scores[ins["start"]] += score
        representative_start = max(
            candidate_positions,
            key=lambda p: support_scores[p]
        )

    matching = [
        ins for ins in insertions
        if ins["start"] == representative_start
    ]
    end_counts = Counter(
        ins["end"]
        for ins in matching
    )
    representative_end = end_counts.most_common(1)[0][0]
    return (
        representative_start,
        representative_end,
        position_counts
    )

# ------------------------------------------------
# OUTILS CIGAR
# ------------------------------------------------

def get_local_cigar_events(
    read,
    site_start,
    softclip_proximity=SOFTCLIP_PROXIMITY,
    indel_proximity=INDEL_PROXIMITY,
    local_window=LOCAL_WINDOW
):

    softclip_near_site = False
    large_indel_near_site = False
    local_events = []

    if not read.cigartuples:
        return softclip_near_site, large_indel_near_site, local_events

    ref_pos = read.reference_start

    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            # M / = / X : Match / Mismatch 
            block_start = ref_pos
            block_end = ref_pos + length
            if (
                block_start <= site_start + local_window
                and block_end >= site_start - local_window
            ):
                local_events.append(("M", block_start, block_end, length))

            ref_pos = block_end

        elif op == 1:
            distance = abs(ref_pos - site_start)
            if distance <= local_window:
                local_events.append(("I", ref_pos, ref_pos, length))
            if distance <= indel_proximity and length > MAX_LOCAL_INDEL:
                large_indel_near_site = True

        elif op in (2, 3):
            block_start = ref_pos
            block_end = ref_pos + length
            if (
                block_start <= site_start + local_window
                and block_end >= site_start - local_window
            ):
                local_events.append(("D" if op == 2 else "N", block_start, block_end, length))
            near_breakpoint = (
                block_start - indel_proximity <= site_start <= block_end + indel_proximity
            )

            if near_breakpoint and (length > MAX_LOCAL_INDEL or op == 3):
                large_indel_near_site = True

            ref_pos = block_end

        elif op == 4:
            distance = abs(ref_pos - site_start)
            if distance <= local_window:
                local_events.append(("S", ref_pos, ref_pos, length))
            if distance <= softclip_proximity and length > softclip_proximity:
                softclip_near_site = True

        elif op == 5:
            distance = abs(ref_pos - site_start)
            if distance <= local_window:
                local_events.append(("H", ref_pos, ref_pos, length))

        elif op == 6:
            pass

    return softclip_near_site, large_indel_near_site, local_events


def is_clean_empty_read(
    read,
    site_start,
    site_end
):
    required_start = max(
        0,
        site_start - SITE_FLANK
    )

    required_end = site_end + SITE_FLANK

    # --------------------------------------------------------
    # 1. Coverage
    # --------------------------------------------------------
    if read.reference_start > required_start:
        return False, "insufficient_left_flank"
    if read.reference_end < required_end:
        return False, "insufficient_right_flank"

    # --------------------------------------------------------
    # 2. CIGAR analysis
    # --------------------------------------------------------
    (
        softclip_near_site,
        large_indel_near_site,
        local_events
    ) = get_local_cigar_events(
        read,
        site_start
    )

    # --------------------------------------------------------
    # 3. Breakpoints softclip
    # --------------------------------------------------------
    if softclip_near_site:
        return False, "softclip_near_site"

    # --------------------------------------------------------
    # 4. Indel close to the breakpoin
    # --------------------------------------------------------
    if large_indel_near_site:
        return False, "large_indel_near_site"

    # --------------------------------------------------------
    # 5. Read clean
    # --------------------------------------------------------
    return True, "clean"

# ------------------------------------------------
# EMPTY READS
# ------------------------------------------------
def get_empty_reads(
    bam,
    chrom,
    site_start,
    site_end,
    valid_te_reads
):

    empty_ids = set()
    ambiguous_ids = set()
    ambiguous_reasons = Counter()

    fetch_start = max(
        0,
        site_start - SITE_FLANK
    )

    fetch_end = site_end + SITE_FLANK

    try:
        fetched_reads = bam.fetch(
            chrom,
            fetch_start,
            fetch_end
        )

    except ValueError:
        return (
            empty_ids,
            ambiguous_ids,
            ambiguous_reasons
        )
    for read in fetched_reads:

        # ----------------------------------------------------
        #  BAM filters
        # ----------------------------------------------------
        if read.is_unmapped:
            continue
        if read.is_secondary:
            continue
        if read.is_supplementary:
            continue
        if read.mapping_quality < MIN_MAPQ:
            continue
        read_id = read.query_name

        # ----------------------------------------------------
        # Known TE+ read
        # ----------------------------------------------------
        if read_id in valid_te_reads:
            continue

        # ----------------------------------------------------
        # EMPTY reads
        # ----------------------------------------------------
        is_empty, reason = is_clean_empty_read(
            read,
            site_start,
            site_end
        )
        if is_empty:
            empty_ids.add(read_id)
        else:
            if reason in {
                "softclip_near_site",
                "large_indel_near_site"
            }:
                ambiguous_ids.add(read_id)
                ambiguous_reasons[reason] += 1

    return (
        empty_ids,
        ambiguous_ids,
        ambiguous_reasons
    )

# ------------------------------------------------
# MAIN
# ------------------------------------------------
def main():
    if len(sys.argv) < 6:
        print(
            "Usage:\n"
            "python evaluate_insertion_sites_v4.py "
            "<clusters.tsv> "
            "<insertions.bed> "
            "<classification.tsv> "
            "<all_reads_bam> "
            "<output.tsv>"
        )
        sys.exit(1)

    clusters_tsv = sys.argv[1]
    insertions_bed = sys.argv[2]
    classification_tsv = sys.argv[3]
    bam_path = sys.argv[4]
    output_file_path = sys.argv[5]

    # --------------------------------------------------------
    # Checks
    # --------------------------------------------------------
    for path, label in [
        (clusters_tsv, "clusters.tsv"),
        (insertions_bed, "insertions.bed"),
        (classification_tsv, "classification.tsv"),
        (bam_path, "BAM")
    ]:

        if not os.path.exists(path):
            print(
                f"[ERREUR] {label} introuvable : {path}",
                file=sys.stderr
            )
            sys.exit(1)

    sample_name = (
        os.path.basename(clusters_tsv)
        .split(".")[0]
    )

    # ------------------------------------------------
    # 1. CLASSIFICATION TE+
    # ------------------------------------------------
    print(
        "[INFO] Load TE classification"
    )

    valid_te_reads = load_classification(
        classification_tsv
    )

    print(
        f"[INFO] {len(valid_te_reads)} reads TE+ found "
        f"(KEEP + Full_Length/5p_Truncated/3p_Truncated)."
    )

    # ------------------------------------------------
    # 2. CLUSTERS
    # ------------------------------------------------
    print(
        "[INFO] Load insertions..."
    )

    clusters = load_clusters(
        clusters_tsv
    )

    print(
        f"[INFO] {len(clusters)} clusters found."
    )

    # ------------------------------------------------
    # 3. INSERTIONS.BED
    # ------------------------------------------------
    print(
        "[INFO] Load insertions TE+..."
    )

    insertions = load_insertions(
        insertions_bed,
        valid_te_reads
    )

    print(
        f"[INFO] {len(insertions)} insertions "
        f"from  TE+ reads valides."
    )

    # ------------------------------------------------
    # 4. ASSOCIATION INSERTIONS -> CLUSTERS
    # ------------------------------------------------
    cluster_insertions = assign_insertions_to_clusters(
        insertions,
        clusters
    )

    print(
        f"[INFO] {len(cluster_insertions)} clusters "
        f"contain at least on TE+ insertions."
    )

    # ------------------------------------------------
    # 5. BAM
    # ------------------------------------------------
    bam = pysam.AlignmentFile(
        bam_path,
        "rb"
    )

    # ------------------------------------------------
    # 6. OUTPUT
    # ------------------------------------------------
    with open(output_file_path, "w") as out_f:

        out_f.write(
            "Sample\t"
            "Chrom\t"
            "Insertion_Start\t"
            "Insertion_End\t"
            "Site_Support\t"
            "N_TE_Full_Length\t"
            "N_TE_5p_Truncated\t"
            "N_TE_3p_Truncated\t"
            "Alt_Reads_TE+\t"
            "Ref_Reads_Empty\t"
            "Excluded_Ambiguous\t"
            "Ambiguous_SoftClip\t"
            "Ambiguous_LargeIndel\t"
            "Total_Informative\t"
            "VAF_Insertion\n"
        )

        # ------------------------------------------------
        # 7. Cluster analysis
        # ------------------------------------------------

        for idx, cl in enumerate(clusters):

            chrom = cl["chrom"]
            cluster_start = cl["start"]
            cluster_end = cl["end"]

            cluster_ins = cluster_insertions.get(
                idx,
                []
            )

            # ------------------------------------------------
            # No TE+ read associated
            # ------------------------------------------------

            if not cluster_ins:

                out_f.write(
                    f"{sample_name}\t"
                    f"{chrom}\t"
                    f"{cluster_start}\t"
                    f"{cluster_end}\t"
                    f"NA\t"
                    f"NA\t"
                    f"0\t"
                    f"0\t"
                    f"0\t"
                    f"0\t"
                    f"0\t"
                    f"0\t"
                    f"0\t"
                    f"0\t"
                    f"0\t"
                    f"0\t"
                    f"0.0%\n"
                )

                continue

            # ------------------------------------------------
            # Representative site
            # ------------------------------------------------

            (
                site_start,
                site_end,
                position_counts
            ) = determine_representative_site(
                cluster_ins
            )

            # ------------------------------------------------
            # IDs TE+ of the insertions
            # ------------------------------------------------

            alt_ids = {
                ins["read_id"]
                for ins in cluster_ins
            }

            # ------------------------------------------------
            # TE categories
            # ------------------------------------------------

            full_length_ids = {
                ins["read_id"]
                for ins in cluster_ins
                if ins["variant"].startswith("Full_Length")
                or ins["variant"].startswith("Full_Length_Like")
                or ins["variant"] == "FLL"
            }

            truncated_5p_ids = {
                ins["read_id"]
                for ins in cluster_ins
                if ins["variant"].startswith("5p_Truncated")
            }

            truncated_3p_ids = {
                ins["read_id"]
                for ins in cluster_ins
                if ins["variant"].startswith("3p_Truncated")
            }

            # ------------------------------------------------
            # EMPTY
            # ------------------------------------------------

            (
                empty_ids,
                ambiguous_ids,
                ambiguous_reasons
            ) = get_empty_reads(
                bam,
                chrom,
                site_start,
                site_end,
                valid_te_reads
            )

            # ------------------------------------------------
            # Computations
            # ------------------------------------------------

            n_alt = len(alt_ids)
            n_empty = len(empty_ids)
            n_ambiguous = len(ambiguous_ids)
            total = n_alt + n_empty
            vaf = (
                n_alt / total * 100
                if total > 0
                else 0
            )

            # ------------------------------------------------
            # Observed positions
            # ------------------------------------------------

            observed_sites = ";".join(
                f"{pos}:{count}"
                for pos, count
                in sorted(position_counts.items())
            )

            site_support = position_counts[site_start]

            # ------------------------------------------------
            # Exclusion explanation
            # ------------------------------------------------

            n_softclip = ambiguous_reasons.get(
                "softclip_near_site",
                0
            )

            n_large_indel = ambiguous_reasons.get(
                "large_indel_near_site",
                0
            )

            # ------------------------------------------------
            # OUTPUT
            # ------------------------------------------------

            out_f.write(
                f"{sample_name}\t"
                f"{chrom}\t"
                f"{site_start}\t"
                f"{site_end}\t"
                f"{site_support}\t"
                f"{len(full_length_ids)}\t"
                f"{len(truncated_5p_ids)}\t"
                f"{len(truncated_3p_ids)}\t"
                f"{n_alt}\t"
                f"{n_empty}\t"
                f"{n_ambiguous}\t"
                f"{n_softclip}\t"
                f"{n_large_indel}\t"
                f"{total}\t"
                f"{vaf:.1f}%\n"
            )

    bam.close()

    print(
        f"[SUCCÈS] Final report generated : "
        f"{output_file_path}"
    )


if __name__ == "__main__":
    main()
