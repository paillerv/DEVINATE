#!/usr/bin/env python3

import sys
import csv


# ============================================================
# USAGE
# ============================================================

if len(sys.argv) != 4:
    print(
        f"Usage: {sys.argv[0]} <clusters.tsv> <vaf.tsv> <output.tsv>",
        file=sys.stderr
    )
    sys.exit(1)

CLUSTERS_FILE = sys.argv[1]
VAF_FILE = sys.argv[2]
OUTPUT_FILE = sys.argv[3]


# ============================================================
# 1. SIMPLIFY VARIANTS
# ============================================================

def simplify_variants(variants_string):

    if not variants_string or variants_string == ".":
        return variants_string

    variants = variants_string.split(";")

    # --------------------------------------------------------
    # Check whether at least one Full_Length_Like is present
    # --------------------------------------------------------

    has_fll = any(
        variant == "Full_Length_Like"
        or variant.startswith("Full_Length_Like:")
        for variant in variants
    )

    # --------------------------------------------------------
    # No FLL:
    # keep truncated variants unchanged
    # --------------------------------------------------------

    if not has_fll:
        return ";".join(variants)

    # --------------------------------------------------------
    # FLL is present:
    #
    # Keep:
    #   Full_Length_Like
    #   Full_Length_Like:DEL:...
    #   Full_Length_Like:INS:...
    #
    # Remove:
    #   3p_Truncated
    #   3p_Truncated:...
    #   5p_Truncated
    #   5p_Truncated:...
    # --------------------------------------------------------

    simplified = []

    for variant in variants:

        # Full-Length without variant
        if variant == "Full_Length_Like":
            simplified.append("FLL")

        # Full-Length with DEL / INS variant
        elif variant.startswith("Full_Length_Like:"):
            simplified.append(
                variant.replace(
                    "Full_Length_Like:",
                    "FLL:",
                    1
                )
            )

        # Remove 3p truncated variants
        elif variant.startswith("3p_Truncated"):
            continue

        # Remove 5p truncated variants
        elif variant.startswith("5p_Truncated"):
            continue

        # Keep anything unexpected
        else:
            simplified.append(variant)

    # --------------------------------------------------------
    # Remove exact duplicates while preserving order
    # --------------------------------------------------------

    seen = set()
    result = []

    for variant in simplified:

        if variant not in seen:
            result.append(variant)
            seen.add(variant)

    return ";".join(result)


# ============================================================
# 2. READ VAF FILE
# ============================================================

print(
    "[INFO] Reading VAF file...",
    file=sys.stderr
)

vaf_data = {}

with open(VAF_FILE, "r", newline="") as f:

    reader = csv.DictReader(
        f,
        delimiter="\t"
    )

    if not reader.fieldnames:
        raise ValueError(
            "VAF file has no header."
        )

    required_columns = {
        "chr",
        "start",
        "end",
        "VAF_with_drops"
    }

    missing = (
        required_columns
        - set(reader.fieldnames)
    )

    if missing:
        raise ValueError(
            "Missing columns in VAF file: "
            + ", ".join(sorted(missing))
        )

    for row in reader:

        key = (
            row["chr"],
            row["start"],
            row["end"]
        )

        if key in vaf_data:

            print(
                f"[WARNING] Duplicate VAF key: {key}. "
                "Keeping the last occurrence.",
                file=sys.stderr
            )

        vaf_data[key] = row


print(
    f"[INFO] VAF entries loaded: {len(vaf_data)}",
    file=sys.stderr
)


# ============================================================
# 3. OUTPUT COLUMNS
# ============================================================

output_columns = [
    "chr",
    "start",
    "end",
    "count",
    "paired_count",
    "confidence",
    "variants",
    "VAF_with_drops"
]


# ============================================================
# 4. MERGE CLUSTERS + VAF
# ============================================================

print(
    f"[INFO] Writing merged output: {OUTPUT_FILE}",
    file=sys.stderr
)

n_clusters = 0


with open(CLUSTERS_FILE, "r", newline="") as infile, \
     open(OUTPUT_FILE, "w", newline="") as outfile:

    reader = csv.DictReader(
        infile,
        delimiter="\t"
    )

    if not reader.fieldnames:
        raise ValueError(
            "Clusters file has no header."
        )

    required_columns = {
        "chr",
        "start",
        "end",
        "count",
        "paired_count",
        "confidence",
        "variants"
    }

    missing = (
        required_columns
        - set(reader.fieldnames)
    )

    if missing:
        raise ValueError(
            "Missing columns in clusters file: "
            + ", ".join(sorted(missing))
        )

    writer = csv.DictWriter(
        outfile,
        fieldnames=output_columns,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()

    for row in reader:

        n_clusters += 1

        # ----------------------------------------------------
        # Simplify variants
        # ----------------------------------------------------

        row["variants"] = simplify_variants(
            row["variants"]
        )

        # ----------------------------------------------------
        # Build genomic key
        # ----------------------------------------------------

        key = (
            row["chr"],
            row["start"],
            row["end"]
        )

        # ----------------------------------------------------
        # Find corresponding VAF entry
        # ----------------------------------------------------

        vaf_row = vaf_data.get(key)

        # ----------------------------------------------------
        # Prepare output row
        # ----------------------------------------------------

        output_row = {
            "chr": row["chr"],
            "start": row["start"],
            "end": row["end"],
            "count": row["count"],
            "paired_count": row["paired_count"],
            "confidence": row["confidence"],
            "variants": row["variants"],
            "VAF_with_drops": ""
        }

        # ----------------------------------------------------
        # Add VAF_with_drops if a match exists
        # ----------------------------------------------------

        if vaf_row is not None:
            output_row["VAF_with_drops"] = vaf_row[
                "VAF_with_drops"
            ]

        writer.writerow(output_row)


# ============================================================
# 5. SUMMARY
# ============================================================

print(
    "",
    file=sys.stderr
)

print(
    "[DONE]",
    file=sys.stderr
)

print(
    f"  Clusters: {n_clusters}",
    file=sys.stderr
)

print(
    f"  Output:   {OUTPUT_FILE}",
    file=sys.stderr
)
