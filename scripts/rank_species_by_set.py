#!/usr/bin/env python3
import argparse
import sys
import re

def sample_to_wwtp_id(sample_name: str) -> str:
    parts = sample_name.split("-")
    if len(parts) >= 2:
        return "-".join(parts[-2:])
    return sample_name

def load_wwtp_metadata(meta_path: str):
    """
    Return dict: WWTP_id -> (Continent, Country/Region, Process_type)
    """
    meta = {}
    with open(meta_path) as f:
        header = f.readline().rstrip("\n")
        cols = header.split("\t")

        idx_id = cols.index("WWTP_id")
        idx_cont = cols.index("Continent")
        idx_country = cols.index("Country/Region")
        idx_proc = cols.index("Process_type")

        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            meta[parts[idx_id]] = (
                parts[idx_cont],
                parts[idx_country],
                parts[idx_proc]
            )
    return meta

def extract_species_token_from_tax(tax: str) -> str:
    spos = tax.find("s__")
    if spos == -1:
        return ""
    s_part = tax[spos + 3:]
    if ";" in s_part:
        s_part = s_part.split(";", 1)[0]
    return s_part.strip()

def is_seqcode_candidate(tax: str) -> str:
    """
    Y if:
      - s__ (empty)
      - s__Genus sp123456
    else N
    """
    s_tok = extract_species_token_from_tax(tax)

    if s_tok == "":
        return "Y"

    if re.search(r"\bsp\d+\b", s_tok):
        return "Y"

    return "N"

def load_species_set(path: str):
    """
    Load species IDs from file, one per line.
    """
    species = set()
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            species.add(line.split()[0])
    return species

def main():
    ap = argparse.ArgumentParser(
        description="Rank a given set of species by highest abundance across global WWTPs."
    )
    ap.add_argument("-f", "--file",
                    default="hq.spRep-abundance-in.midas_dk.datasets-gtdb.taxa.tsv")
    ap.add_argument("-m", "--metadata",
                    default="wwtp.metadata.tsv")
    ap.add_argument("--global-end",
                    type=int, default=84)
    ap.add_argument("-S", "--species-set",
                    required=True,
                    help="File containing species IDs (one per line).")
    ap.add_argument("--process-filter",
                    default="none",
                    help='Process_type to keep, e.g. "C,N,DN,P", or "none".')
    ap.add_argument("--detect-cutoff",
                    type=float, default=0.0,
                    help="Detection cutoff: abundance > cutoff. Default 0.")
    args = ap.parse_args()

    species_set = load_species_set(args.species_set)
    if not species_set:
        sys.stderr.write("ERROR: species set is empty.\n")
        sys.exit(1)

    do_filter = (args.process_filter.lower() != "none")
    cutoff = args.detect_cutoff

    meta = load_wwtp_metadata(args.metadata)

    with open(args.file) as f:
        header = f.readline().rstrip("\n")
        cols = header.split("\t")
        sample_names_all = cols[1:-1]

        global_indices = [
            j for j in range(len(sample_names_all))
            if (j + 2) <= args.global_end
        ]
        global_sample_names = [sample_names_all[j] for j in global_indices]

        keep = [True] * len(global_indices)
        if do_filter:
            for k, sname in enumerate(global_sample_names):
                w_id = sample_to_wwtp_id(sname)
                keep[k] = (meta.get(w_id, ("", "", ""))[2] == args.process_filter)

        results = []

        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(cols):
                continue

            sp = parts[0]
            if sp not in species_set:
                continue

            tax = parts[-1]
            abund_strs = parts[1:-1]

            max_ab = -1.0
            max_sample = "NA"
            max_cont = "NA"
            max_country = "NA"
            detected_wwtps = []

            for k, j in enumerate(global_indices):
                if not keep[k]:
                    continue

                sname = sample_names_all[j]
                w_id = sample_to_wwtp_id(sname)
                cont, country, proc = meta.get(w_id, ("Unknown", "Unknown", "Unknown"))

                try:
                    val = float(abund_strs[j])
                except ValueError:
                    val = 0.0

                if val > cutoff:
                    detected_wwtps.append(sname)

                if val > max_ab:
                    max_ab = val
                    max_sample = sname
                    max_cont = cont
                    max_country = country

            if max_ab < 0:
                max_ab = 0.0
                detected_wwtps = []

            results.append((
                sp,
                max_ab,
                max_sample,
                max_cont,
                max_country,
                len(detected_wwtps),
                ",".join(detected_wwtps),
                tax,
                is_seqcode_candidate(tax)
            ))

    results.sort(key=lambda x: (-x[1], -x[5], x[0]))

    print(
        "SpeciesID\tHighest_abundance\tHighest_abundance_WWTP\t"
        "Highest_abundance_Continent\tHighest_abundance_Country\t"
        "Number_of_WWTP_detected\tWWTP_detected\tGTDB_taxonomy\tseqcode_candidates"
    )

    for r in results:
        print(
            f"{r[0]}\t{r[1]:.6f}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]}\t{r[7]}\t{r[8]}"
        )

if __name__ == "__main__":
    main()

