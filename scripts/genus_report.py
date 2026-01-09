#!/usr/bin/env python3
"""
Genus-level abundance reporter for MiDAS / Danish WWTP datasets.

Usage (basic):
    python genus_report.py -g Accumulibacter

Optional:
    python genus_report.py -g Accumulibacter \
        -f hq.spRep-abundance-in.midas_dk.datasets-gtdb.taxa.tsv \
        --global-end 84

Input assumptions:
- Tab-separated file.
- Column 1      : Species ID
- Columns 2..N-1: Sample abundances (first global WWTPs, then Danish WWTPs)
- Column N      : GTDB taxonomy with a "g__Genus" field.
"""

import sys
import argparse
import numpy as np
from collections import defaultdict


# --------------------------
# Argument Parser
# --------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize species-level abundance patterns for a given genus."
    )
    parser.add_argument(
        "-g", "--genus",
        required=True,
        help="Genus name to match within taxonomy (e.g. 'Accumulibacter')."
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
# Main program
# --------------------------

def main():
    args = parse_args()
    genus = args.genus
    fname = args.file
    GLOBAL_END = args.global_end

    # species_pos: species -> list of (abundance>0, wwtp_name, type)
    species_pos = defaultdict(list)

    # species_all: species -> list of abundances across all samples (including zeros)
    species_all = defaultdict(list)

    # species_tax: species -> GTDB taxonomy string
    species_tax = {}

    # wwtp_types: wwtp_name -> "global" or "dk" (only for samples with abundance > 0)
    wwtp_types = {}

    # Track global maximum abundance across this genus
    max_val = -1.0
    max_species = None
    max_wwtp = None

    # --------------------------
    # Read file
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
            tax = parts[-1]
            abund_strs = parts[1:-1]

            # Extract genus from taxonomy: locate "g__"
            gpos = tax.find("g__")
            if gpos == -1:
                continue
            g_part = tax[gpos + 3:].split(";", 1)[0]  # substring after g__ up to next ';'

            # Only keep species whose genus field contains the requested genus string
            if genus not in g_part:
                continue

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
    # Genus-level stats
    # --------------------------

    # Only consider species that have at least one non-zero abundance
    species_list = [sp for sp in species_pos if len(species_pos[sp]) > 0]
    nspecies = len(species_list)

    wwtps_detected = set(wwtp_types.keys())
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

    # Section: overall summary
    print(center_title(" overall summary "))

    if nspecies == 0:
        # No species with abundance > 0 for this genus
        print(
            f"Genus {genus}: 0 HQ species detected with abundance > 0 "
            f"in 0/{n_samples_total} samples."
        )

        print(center_title(" species-level statistics "))
        print("SpeciesID\tGTDB_taxonomy\tN_pos\tN_pos_global\tN_pos_dk\t"
              "Min_pos\tQ1_pos\tMedian_pos\tQ3_pos\tMax_pos\t"
              "Min_all\tQ1_all\tMedian_all\tQ3_all\tMax_all")

        print(center_title(" species abundance distribution "))
        print("SpeciesID\tGTDB_taxonomy\tAbundance\tWWTP")
        return

    # taxonomy of the highest-abundance species
    best_tax = species_tax.get(max_species, "NA")

    # Summary line with highest abundance + most prevalent species info
    print(
        f"###   Genus {genus}:   ###\n{nspecies} HQ species detected in "
        f"{nsamples_detected}/{n_samples_total} datasets (abundance > 0), "
        f"({len(global_wwtps)} global + {len(dk_wwtps)} Danish datasets).\n" 
        f"The highest abundance is {max_val}% in {max_wwtp} (species {max_species}, taxa: {best_tax}).\n"
        f"The most prevalent HQ species is {most_prev_sp}, taxa: {most_prev_tax}, detected in {most_prev_total}/{n_samples_total} datasets "
        f"({most_prev_global} global + {most_prev_dk} Danish datasets).\n"
        f"*note: 1 dataset for 1 global wwtp and 3 datasets for 1 Danish wwtp."
    )

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

