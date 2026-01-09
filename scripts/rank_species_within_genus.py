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
            if len(parts) <= max(idx_id, idx_cont, idx_country, idx_proc):
                continue
            meta[parts[idx_id]] = (
                parts[idx_cont],
                parts[idx_country],
                parts[idx_proc]
            )
    return meta

def extract_genus_from_tax(tax: str) -> str:
    gpos = tax.find("g__")
    if gpos == -1:
        return ""
    return tax[gpos + 3:].split(";", 1)[0]

def extract_species_token_from_tax(tax: str) -> str:
    """
    Return the species token after 's__' (could be empty).
    Example:
      ...;s__Accumulibacter propinquus -> 'Accumulibacter propinquus'
      ...;s__Accumulibacter sp017302455 -> 'Accumulibacter sp017302455'
      ...;s__ -> ''
    """
    spos = tax.find("s__")
    if spos == -1:
        return ""
    s_part = tax[spos + 3:]
    # truncate at next ';' if present
    if ";" in s_part:
        s_part = s_part.split(";", 1)[0]
    return s_part.strip()

def is_seqcode_candidate(tax: str) -> str:
    """
    Y if species name is missing or placeholder (SeqCode candidate / not a proper Latin binomial).
    Criteria (based on user's definition):
      - s__ (empty)
      - s__<Genus> sp<digits> (e.g. 'Accumulibacter sp017302455')
    Otherwise N (e.g. 'Accumulibacter propinquus')
    """
    s_tok = extract_species_token_from_tax(tax)

    # Case 1: empty species label
    if s_tok == "":
        return "Y"

    # Case 2: placeholder 'Genus sp123...' or 'Genus sp012345'
    # Allow "sp" + digits, optionally separated by space
    if re.search(r"\bsp\d+\b", s_tok):
        return "Y"

    # Otherwise treat as having a real species epithet
    return "N"

def main():
    ap = argparse.ArgumentParser(
        description="Rank species within a genus by highest abundance across global WWTPs."
    )
    ap.add_argument("-f", "--file",
                    default="hq.spRep-abundance-in.midas_dk.datasets-gtdb.taxa.tsv")
    ap.add_argument("-m", "--metadata",
                    default="wwtp.metadata.tsv")
    ap.add_argument("--global-end",
                    type=int, default=84)
    ap.add_argument("-i", "--input-genus",
                    required=True)
    ap.add_argument("--process-filter",
                    default="none",
                    help='Process_type to keep, e.g. "C,N,DN,P", or "none"')
    ap.add_argument("--detect-cutoff",
                    type=float, default=0.0,
                    help="Detection cutoff (abundance > cutoff). Default: 0")
    args = ap.parse_args()

    genus = args.input_genus
    if genus.startswith("g__"):
        genus = genus[3:]

    do_filter = (args.process_filter.lower() != "none")
    cutoff = args.detect_cutoff

    meta = load_wwtp_metadata(args.metadata)

    with open(args.file) as f:
        header = f.readline().rstrip("\n")
        if not header:
            sys.stderr.write("ERROR: abundance file missing header.\n")
            sys.exit(1)

        cols = header.split("\t")
        sample_names_all = cols[1:-1]

        global_indices = [
            j for j in range(len(sample_names_all))
            if (j + 2) <= args.global_end
        ]
        global_sample_names = [sample_names_all[j] for j in global_indices]

        # process-type mask
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
            tax = parts[-1]

            if extract_genus_from_tax(tax) != genus:
                continue

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

            seqcode = is_seqcode_candidate(tax)

            results.append((
                sp,
                max_ab,
                max_sample,
                max_cont,
                max_country,
                len(detected_wwtps),
                ",".join(detected_wwtps),
                tax,
                seqcode
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

