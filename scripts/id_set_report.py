#!/usr/bin/env python3
"""
ID-set-level abundance reporter for MiDAS / Danish WWTP datasets.

Instead of grouping by genus, this script summarizes a user-selected set
of species IDs (e.g. comammox MAGs).

Usage (basic):
    python id_set_report.py -i comammox.ids

Optional:
    python id_set_report.py -i comammox.ids \
        -L comammox-taxa \
        -f hq.spRep-abundance-in.midas_dk.datasets-gtdb.taxa.tsv \
        --global-end 84

Input assumptions for the abundance+taxonomy matrix:
- Tab-separated file.
- Column 1      : Species ID (e.g. 001-SG-02.bin.1.500)
- Columns 2..N-1: Sample abundances (first global WWTPs, then Danish WWTPs)
- Column N      : GTDB taxonomy string.
"""

import sys
import argparse
import numpy as np
from collections import defaultdict


# --------------------------
# Argument parser
# --------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize species-level abundance patterns for a selected set of species IDs."
    )
    parser.add_argument(
        "-i", "--id-file",
        required=True,
        help="File containing species IDs of interest (first column = ID)."
    )
    parser.add_argument(
        "-L", "--label",
        default="selected IDs",
        help="Label/name for this ID set (used in the summary header)."
    )
    parser.add_argument(
        "-f", "--file",
        default="hq.spRep-abundance-in.midas_dk.datasets-gtdb.taxa.tsv",
        help="Input TSV file containing abundances and GTDB taxonomy."
    )
    parser.add_argument(
        "--global-end",
        type=int,
        default=84,
        help="Last global-sample column index (1-based). Columns 2..global-end are global WWTPs; the rest are Danish."
    )
    return parser.parse_args()


# --------------------------
# Helpers
# --------------------------

def load_id_set(path):
    ids = set()
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                # Take the first column as ID
                id_ = line.split()[0]
                ids.add(id_)
    except OSError as e:
        sys.stderr.write(f"Error opening ID file '{path}': {e}\n")
        sys.exit(1)
    return ids


# --------------------------
# Main
# --------------------------

def main():
    args = parse_args()
    id_file = args.id_file
    label = args.label
    fname = args.file
    GLOBAL_END = args.global_end

    selected_ids = load_id_set(id_file)
    if not selected_ids:
        sys.stderr.write(f"No IDs loaded from {id_file}.\n")
        sys.exit(1)

    # species_pos: species -> list of (abundance>0, wwtp_name, type)
    species_pos = defaultdict(list)

    # species_all: species -> list of abundances across all samples (including zeros)
    species_all = defaultdict(list)

    # species_tax: species -> GTDB taxonomy string
    species_tax = {}

    # wwtp_types: wwtp_name -> "global" or "dk" (only for samples with abundance > 0)
    wwtp_types = {}

    # Track global maximum abundance across this ID set
    max_val = -1.0
    max_species = None
    max_wwtp = None

    # --------------------------
    # Read abundance+taxonomy file
    # --------------------------
    try:
        f = open(fname)
    except OSError as e:
        sys.stderr.write(f"Error opening file '{fname}': {e}\n")
        sys.exit(1)

    with f:
        header = f.readline().rstrip("\n")
        if not header:
            sys.stderr.write("Input file appears to be empty or missing a header line.\n")
            sys.exit(1)

        cols = header.split("\t")
        sample_names = cols[1:-1]          # sample columns are 2..N-1
        n_samples_total = len(sample_names)

        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) != len(cols):
                # Skip malformed lines
                continue

            sp = parts[0]
            if sp not in selected_ids:
                continue

            tax = parts[-1]
            abund_strs = parts[1:-1]

            species_tax[sp] = tax

            # Iterate over all sample columns
            for j, a in enumerate(abund_strs):
                try:
                    val = float(a)
                except ValueError:
                    val = 0.0   # treat non-numeric as zero

                species_all[sp].append(val)

                wwtp = sample_names[j]
                # Original column index (1-based): 1=Species, 2=first sample, ...
                col_index = j + 2
                sample_type = "global" if col_index <= GLOBAL_END else "dk"

                if val > 0:
                    species_pos[sp].append((val, wwtp, sample_type))
                    wwtp_types[wwtp] = sample_type

                    if val > max_val:
                        max_val = val
                        max_species = sp
                        max_wwtp = wwtp

    # --------------------------
    # ID-set-level stats
    # --------------------------

    # Only consider species that have at least one non-zero abundance
    species_list = [sp for sp in species_pos if len(species_pos[sp]) > 0]
    nspecies = len(species_list)

    wwtps_detected = set(wwtps for wwtps in wwtp_types.keys())
    nsamples_detected = len(wwtps_detected)

    global_wwtps = {w for w, t in wwtp_types.items() if t == "global"}
    dk_wwtps = {w for w, t in wwtp_types.items() if t == "dk"}

    # --------------------------
    # Species-level stats table
    # --------------------------

    # Each record in species_stats:
    # (sp, tax,
    #  n_pos, n_pos_global, n_pos_dk,
    #  mn_pos, q1_pos, median_pos, q3_pos, mx_pos,
    #  mn_all, q1_all, median_all, q3_all, mx_all,
    #  sort_key)
    species_stats = []

    for sp in species_list:
        tax = species_tax[sp]

        # Positive abundances (>0)
        pos_vals = np.array([v[0] for v in species_pos[sp]], dtype=float)
        # All abundances (including zeros)
        all_vals = np.array(species_all[sp], dtype=float)

        n_pos = len(pos_vals)
        n_pos_global = sum(1 for v in species_pos[sp] if v[2] == "global")
        n_pos_dk = sum(1 for v in species_pos[sp] if v[2] == "dk")

        # Stats based on positive abundances only
        mn_pos = float(pos_vals.min())
        mx_pos = float(pos_vals.max())
        q1_pos, q2_pos, q3_pos = np.percentile(pos_vals, [25, 50, 75])

        # Stats based on all samples (including zeros)
        mn_all = float(all_vals.min())
        mx_all = float(all_vals.max())
        q1_all, q2_all, q3_all = np.percentile(all_vals, [25, 50, 75])

        species_stats.append(
            (
                sp, tax,
                n_pos, n_pos_global, n_pos_dk,
                mn_pos, q1_pos, q2_pos, q3_pos, mx_pos,
                mn_all, q1_all, q2_all, q3_all, mx_all,
                mx_pos  # sort key
            )
        )

    # Sort species by maximum positive abundance (descending)
    species_stats.sort(key=lambda x: x[-1], reverse=True)

    # If we have any species, identify the most prevalent (max n_pos)
    if nspecies > 0:
        most_prev_entry = max(species_stats, key=lambda x: x[2])  # x[2] = n_pos
        most_prev_sp = most_prev_entry[0]
        most_prev_tax = most_prev_entry[1]
        most_prev_total = most_prev_entry[2]
        most_prev_global = most_prev_entry[3]
        most_prev_dk = most_prev_entry[4]
    else:
        most_prev_sp = None
        most_prev_tax = None
        most_prev_total = 0
        most_prev_global = 0
        most_prev_dk = 0

    # --------------------------
    # Formatting helpers
    # --------------------------

    LINE_WIDTH = 80

    def center_title(text: str) -> str:
        """Return a centered title line with '*' padding."""
        return text.center(LINE_WIDTH, "*")

    # --------------------------
    # OUTPUT
    # --------------------------

    print(center_title(" overall summary "))

    if nspecies == 0:
        print(f"###   Group {label}:   ###")
        print(f"0 HQ species from the selected ID set detected with abundance > 0 "
              f"in 0/{n_samples_total} datasets.")
        print("*note: 1 dataset for 1 global WWTP and 3 datasets for 1 Danish WWTP.")
        # Still print empty tables for consistency
        print(center_title(" species-level statistics "))
        print("SpeciesID\tGTDB_taxonomy\tN_pos\tN_pos_global\tN_pos_dk\t"
              "Min_pos\tQ1_pos\tMedian_pos\tQ3_pos\tMax_pos\t"
              "Min_all\tQ1_all\tMedian_all\tQ3_all\tMax_all")
        print(center_title(" species abundance distribution "))
        print("SpeciesID\tGTDB_taxonomy\tAbundance\tWWTP")
        return

    # Taxonomy of the highest-abundance species
    best_tax = species_tax.get(max_species, "NA")

    # ---- description block in your style ----
    print(f"###   Group {label}:   ###")
    print(
        f"{nspecies} HQ species detected in "
        f"{nsamples_detected}/{n_samples_total} datasets (abundance > 0), "
        f"({len(global_wwtps)} global + {len(dk_wwtps)} Danish datasets)."
    )
    # Here we assume max_val is already in percent units (e.g. 8.0732 means 8.0732%)
    print(
        f"The highest abundance is {max_val:.4f}% in {max_wwtp} "
        f"(species {max_species}, taxa: {best_tax})."
    )
    print(
        f"The most prevalent HQ species is {most_prev_sp} "
        f"(taxa: {most_prev_tax}), detected in "
        f"{most_prev_total}/{n_samples_total} datasets "
        f"({most_prev_global} global + {most_prev_dk} Danish datasets)."
    )
    print("*note: 1 dataset for 1 global WWTP and 3 datasets for 1 Danish WWTP.")

    # Section: species-level statistics
    print(center_title(" species-level statistics "))

    print("SpeciesID\tGTDB_taxonomy\tN_pos\tN_pos_global\tN_pos_dk\t"
          "Min_pos\tQ1_pos\tMedian_pos\tQ3_pos\tMax_pos\t"
          "Min_all\tQ1_all\tMedian_all\tQ3_all\tMax_all")

    for (
        sp, tax,
        n_pos, n_pos_global, n_pos_dk,
        mn_pos, q1_pos, q2_pos, q3_pos, mx_pos,
        mn_all, q1_all, q2_all, q3_all, mx_all,
        _sort_key
    ) in species_stats:
        print(
            f"{sp}\t{tax}\t"
            f"{n_pos}\t{n_pos_global}\t{n_pos_dk}\t"
            f"{mn_pos}\t{q1_pos}\t{q2_pos}\t{q3_pos}\t{mx_pos}\t"
            f"{mn_all}\t{q1_all}\t{q2_all}\t{q3_all}\t{mx_all}"
        )

    # Section: species abundance distribution (long table)
    print(center_title(" species abundance distribution "))

    print("SpeciesID\tGTDB_taxonomy\tAbundance\tWWTP")

    for sp, tax, *_rest in species_stats:
        vals = species_pos[sp]  # list of (abundance, wwtp, type) with abundance > 0
        vals_sorted = sorted(vals, key=lambda v: v[0], reverse=True)
        for val, wwtp, sample_type in vals_sorted:
            print(f"{sp}\t{tax}\t{val}\t{wwtp}")


if __name__ == "__main__":
    main()

