#!/usr/bin/env python3
"""Validation and TSV summaries for run_overlap.sh. Python standard library only."""
import argparse
import csv
import gzip
import math
import os
import re
import sys
from collections import Counter, defaultdict
from contextlib import ExitStack
from pathlib import Path

RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
PREFIXES = ["d__", "p__", "c__", "o__", "f__", "g__", "s__"]
META_COLUMNS = ["GenomeID", "path", "taxonomy"] + RANKS
RAW_COLUMNS = ["Ref_file", "Query_file", "ANI", "Align_fraction_ref", "Align_fraction_query"]
UNASSIGNED_TAXA = {"", "na", "n/a", "none", "null", "-", "unclassified", "unassigned", "unknown", "not classified"}


def fail(message):
    raise ValueError(message)


def tsv_text(columns, rows):
    import io
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(columns)
    writer.writerows(rows)
    return stream.getvalue()


def write_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.read_text(encoding="utf-8") == text:
        return
    temp = path.with_name(path.name + ".tmp")
    temp.write_text(text, encoding="utf-8")
    temp.replace(path)


def write_tsv(path, columns, records):
    write_text(path, tsv_text(columns, [[record.get(c, "") for c in columns] for record in records]))


def parse_taxonomy(text, context):
    text = text.strip()
    parsed = dict.fromkeys(RANKS, "")
    if text.casefold() in UNASSIGNED_TAXA:
        return text, parsed
    # GTDB-Tk may report a domain-level, unclassified result without rank prefixes.
    # Preserve this source classification and use only its explicitly stated domain.
    unclassified_domain = re.fullmatch(r"unclassified\s+(bacteria|archaea)", text, re.I)
    if unclassified_domain:
        parsed["domain"] = "d__" + {"bacteria": "Bacteria", "archaea": "Archaea"}[unclassified_domain[1].casefold()]
        return text, parsed
    tokens = [token.strip() for token in text.split(";") if token.strip()]
    seen = set()
    for token in tokens:
        prefix = token[:3]
        if prefix not in PREFIXES or prefix in seen:
            fail(f"{context}: invalid/duplicate GTDB rank: {token!r}")
        seen.add(prefix)
        rank = RANKS[PREFIXES.index(prefix)]
        value = token[3:].strip()
        if value.casefold() not in UNASSIGNED_TAXA:
            parsed[rank] = prefix + value
    # Keep the original taxonomy column; rank columns contain parsed assignments.
    return text, parsed


def fasta_probe(path):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = False
        for _ in range(10000):
            line = handle.readline()
            if not line:
                break
            line = line.strip()
            if not line:
                continue
            if not header:
                if not line.startswith(">"):
                    fail(f"Not a FASTA header: {path}")
                header = True
            elif not line.startswith(">"):
                if re.fullmatch(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+", line):
                    return
                fail(f"Not a nucleotide FASTA sequence: {path}")
    fail(f"No nucleotide sequence found near the start of FASTA: {path}")


def load_metadata(filename):
    with open(filename, encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != META_COLUMNS:
            fail(f"Not a normalized metadata table: {filename}")
        rows = list(reader)
    if not rows:
        fail(f"Empty normalized metadata: {filename}")
    return rows, {row["path"]: row for row in rows}


def check_raw(filename):
    with open(filename, encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"), [])
    if any(column not in header for column in RAW_COLUMNS):
        fail(f"Invalid skani header in {filename}; expected {', '.join(RAW_COLUMNS)}")


def audit_log(filename):
    suspicious = []
    with open(filename, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if re.search(r"warn|error", line, re.I) and re.search(
                r"fail|could not|unable to|skipp|no (?:valid )?(?:sequence|record|k.?mer)|"
                r"(?:empty|invalid) (?:file|fasta|genome)|too (?:small|short)", line, re.I
            ):
                suspicious.append(line.strip())
    if suspicious:
        fail(f"Potential skipped/invalid genome(s) in {filename}. Inspect this log before interpreting no-hit genomes.\n" + "\n".join(suspicious[:10]))


def agreement(query, reference, rank):
    if not query[rank] or not reference[rank]:
        return "unassigned"
    return "same" if query[rank] == reference[rank] else "different"


def pct(number, total):
    return f"{100 * number / total:.4f}" if total else "NA"


def number(value):
    return f"{value:.2f}"


def lineage(row, rank):
    return ";".join(row[r] or p for r, p in zip(RANKS, PREFIXES) if RANKS.index(r) <= RANKS.index(rank))


def validate_labels(label_a, label_b):
    for label in (label_a, label_b):
        if not label.strip() or any(c in label for c in "\t\r\n"):
            fail("Catalog labels must be nonempty single-line values without tabs")
    if label_a == label_b:
        fail("Use distinct display labels for catalogs A and B")


def prepare(args):
    validate_labels(args.label_a, args.label_b)
    root = Path(args.root).resolve()
    prepared, probed = [], set()
    for side, filename, label in [("a", args.catalog_a, args.label_a), ("b", args.catalog_b, args.label_b)]:
        source = Path(filename)
        if not source.is_absolute():
            source = root / source
        rows, seen_ids, seen_paths, stats = [], set(), set(), []
        with source.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = [x.strip() for x in next(reader, [])]
            if header not in (["GenomeID", "path"], ["GenomeID", "path", "taxonomy"]):
                fail(f"{source}: expected TAB-delimited columns GenomeID, path and optional taxonomy")
            for line_no, values in enumerate(reader, 2):
                if not values or all(not value.strip() for value in values):
                    continue
                if len(values) != len(header):
                    fail(f"{source}:{line_no}: expected {len(header)} columns, found {len(values)}")
                genome_id, raw_path = [value.strip() for value in values[:2]]
                if not genome_id or not raw_path or any(c in genome_id + raw_path for c in "\t\r\n"):
                    fail(f"{source}:{line_no}: empty ID/path or embedded tab/newline")
                path = Path(raw_path)
                if not path.is_absolute():
                    path = root / path
                path = str(path.resolve())
                if any(c in path for c in "\t\r\n"):
                    fail(f"{source}:{line_no}: resolved FASTA path contains a tab/newline")
                if genome_id in seen_ids or path in seen_paths:
                    fail(f"{source}:{line_no}: duplicate GenomeID or resolved FASTA path within catalog {side.upper()}")
                seen_ids.add(genome_id)
                seen_paths.add(path)
                if not Path(path).is_file() or not os.access(path, os.R_OK):
                    fail(f"{source}:{line_no}: missing/unreadable FASTA: {path}")
                stat = os.stat(path)
                if stat.st_size == 0:
                    fail(f"Empty FASTA: {path}")
                if path not in probed:
                    fasta_probe(path)
                    probed.add(path)
                text = values[2] if len(values) == 3 else ""
                if any(c in text for c in "\t\r\n"):
                    fail(f"{source}:{line_no}: taxonomy contains a tab/newline")
                taxonomy, parsed = parse_taxonomy(text, f"{source}:{line_no}")
                rows.append(dict(GenomeID=genome_id, path=path, taxonomy=taxonomy, **parsed))
                stats.append((path, stat.st_size, stat.st_mtime_ns))
        if not rows:
            fail(f"Empty metadata: {source}")
        prepared.append((side, source, label, rows, stats))
    # Validate both catalogs before replacing any normalized inputs.
    out, counts = Path(args.out), []
    for side, source, label, rows, stats in prepared:
        name = f"catalog_{side}"
        write_tsv(out / f"{name}.metadata.tsv", META_COLUMNS, rows)
        write_text(out / f"{name}.paths.txt", "".join(row["path"] + "\n" for row in rows))
        write_text(out / f"{name}.files.tsv", tsv_text(["path", "size_bytes", "mtime_ns"], sorted(stats)))
        species_labels = Counter(row["species"] for row in rows if row["species"])
        counts.append(dict(catalog_side=side, catalog_label=label, n_representatives=len(rows),
                           unassigned_species_labels=sum(not row["species"] for row in rows),
                           repeated_species_labels=sum(count > 1 for count in species_labels.values()),
                           source_metadata=str(source.resolve())))
    write_tsv(out / "catalog_counts.tsv", list(counts[0]), counts)


PAIR_COLUMNS = [
    "catalog_a_label", "catalog_b_label", "catalog_a_id", "catalog_b_id", "ANI",
    "AF_catalog_a", "AF_catalog_b", "AF_max", "AF_min", "passes_primary",
    "catalog_a_taxonomy", "catalog_b_taxonomy", "species_taxonomy_agreement", "genus_taxonomy_agreement"
]


def summarize(args):
    validate_labels(args.label_a, args.label_b)
    rows, paths = {}, {}
    rows["a"], paths["a"] = load_metadata(args.catalog_a)
    rows["b"], paths["b"] = load_metadata(args.catalog_b)
    by_id = {side: {row["GenomeID"]: row for row in rows[side]} for side in ("a", "b")}
    indices = {side: {row["GenomeID"]: i for i, row in enumerate(rows[side])} for side in ("a", "b")}
    labels = {"a": args.label_a, "b": args.label_b}
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    check_raw(args.raw)
    reported = {side: Counter() for side in ("a", "b")}
    matches = {side: Counter() for side in ("a", "b")}
    low_af = {side: Counter() for side in ("a", "b")}
    best_any = {side: {} for side in ("a", "b")}
    best_match = {side: {} for side in ("a", "b")}
    near_ids = {side: set() for side in ("a", "b")}
    rules = [("primary", args.ani, args.af, args.af_mode),
             ("ANI95_maxAF50", 95.0, 50.0, "max"),
             ("ANI95_maxAF65", 95.0, 65.0, "max"),
             ("ANI95_bothAF65", 95.0, 65.0, "both")]
    sensitivity = {name: [set(), set(), 0] for name, *_ in rules}
    seen_pairs = set()

    def selected_af(pair):
        return pair["af_max"] if args.af_mode == "max" else pair["af_min"]

    def best_key(pair, other):
        return (-pair["ani"], -selected_af(pair), -pair["af_min"], pair[other])

    def update_best(table, genome, pair, other):
        if genome not in table or best_key(pair, other) < best_key(table[genome], other):
            table[genome] = pair

    def annotate(pair, passes=True):
        a, b = by_id["a"][pair["a"]], by_id["b"][pair["b"]]
        return dict(catalog_a_label=args.label_a, catalog_b_label=args.label_b,
                    catalog_a_id=pair["a"], catalog_b_id=pair["b"], ANI=number(pair["ani"]),
                    AF_catalog_a=number(pair["af_a"]), AF_catalog_b=number(pair["af_b"]),
                    AF_max=number(pair["af_max"]), AF_min=number(pair["af_min"]), passes_primary=int(passes),
                    catalog_a_taxonomy=a["taxonomy"], catalog_b_taxonomy=b["taxonomy"],
                    species_taxonomy_agreement=agreement(a, b, "species"), genus_taxonomy_agreement=agreement(a, b, "genus"))

    with ExitStack() as stack:
        writers = {}
        for name in ("matched_pairs", "boundary_review"):
            writers[name] = csv.DictWriter(stack.enter_context((out / f"{name}.tsv").open("w", encoding="utf-8", newline="")),
                                           fieldnames=PAIR_COLUMNS, delimiter="\t", lineterminator="\n")
            writers[name].writeheader()
        raw = stack.enter_context(open(args.raw, encoding="utf-8", newline=""))
        for line_no, record in enumerate(csv.DictReader(raw, delimiter="\t"), 2):
            if not any(record.values()):
                continue
            apath, bpath = record.get("Query_file"), record.get("Ref_file")
            if apath not in paths["a"] or bpath not in paths["b"]:
                fail(f"{args.raw}:{line_no}: FASTA path missing from catalog metadata")
            aid, bid = paths["a"][apath]["GenomeID"], paths["b"][bpath]["GenomeID"]
            code = indices["a"][aid] * len(rows["b"]) + indices["b"][bid]
            if code in seen_pairs:
                fail(f"Repeated pair {aid} / {bid}; use whole-genome inputs without --qi/--ri")
            seen_pairs.add(code)
            values = [float(record[column]) for column in ("ANI", "Align_fraction_query", "Align_fraction_ref")]
            if any(not math.isfinite(value) or not 0 <= value <= 100 for value in values):
                fail(f"{args.raw}:{line_no}: ANI/AF must be finite percentages from 0 to 100")
            ani, af_a, af_b = values
            pair = dict(a=aid, b=bid, ani=ani, af_a=af_a, af_b=af_b, af_max=max(af_a, af_b), af_min=min(af_a, af_b))
            passes = ani >= args.ani and selected_af(pair) >= args.af
            high_ani_low_af = ani >= args.ani and selected_af(pair) < args.af
            near = ((abs(ani - args.ani) <= 0.25 and selected_af(pair) >= min(50, args.af)) or
                    (ani >= args.ani and abs(selected_af(pair) - args.af) <= 2))
            for side, other in (("a", "b"), ("b", "a")):
                genome = pair[side]
                reported[side][genome] += 1
                update_best(best_any[side], genome, pair, other)
                if high_ani_low_af:
                    low_af[side][genome] += 1
                if near:
                    near_ids[side].add(genome)
                if passes:
                    matches[side][genome] += 1
                    update_best(best_match[side], genome, pair, other)
            for name, ani_cut, af_cut, mode in rules:
                af = pair["af_max"] if mode == "max" else pair["af_min"]
                if ani >= ani_cut and af >= af_cut:
                    sensitivity[name][0].add(aid)
                    sensitivity[name][1].add(bid)
                    sensitivity[name][2] += 1
            if near:
                writers["boundary_review"].writerow(annotate(pair, passes))
            if passes:
                writers["matched_pairs"].writerow(annotate(pair))

    statuses = {}
    for side, other in (("a", "b"), ("b", "a")):
        statuses[side] = []
        for row in rows[side]:
            genome = row["GenomeID"]
            matched = bool(matches[side][genome])
            reason = "" if matched else ("high_ANI_insufficient_AF" if low_af[side][genome] else
                                        "below_ANI_threshold" if reported[side][genome] else "no_reported_hit")
            status = dict(catalog_side=side, catalog_label=labels[side], other_catalog_label=labels[other], **row,
                          status="matched" if matched else "unmatched", unmatched_reason=reason,
                          n_reported_hits=reported[side][genome], n_matching_genomes=matches[side][genome],
                          multiple_matches=int(matches[side][genome] > 1), near_threshold=int(genome in near_ids[side]))
            for prefix, table in (("best_match", best_match[side]), ("best_reported", best_any[side])):
                pair = table.get(genome)
                status[f"{prefix}_other_id"] = pair[other] if pair else ""
                for field, key in (("ANI", "ani"), ("AF_catalog_a", "af_a"), ("AF_catalog_b", "af_b"), ("AF_max", "af_max"), ("AF_min", "af_min")):
                    status[f"{prefix}_{field}"] = number(pair[key]) if pair else ""
                status[f"{prefix}_other_taxonomy"] = by_id[other][pair[other]]["taxonomy"] if pair else ""
            statuses[side].append(status)
        columns = list(statuses[side][0])
        write_tsv(out / f"catalog_{side}_status.tsv", columns, statuses[side])
        for name in ("matched", "unmatched"):
            selected = [row for row in statuses[side] if row["status"] == name]
            write_tsv(out / f"catalog_{side}_{name}.tsv", columns, selected)
            write_text(out / f"catalog_{side}_{name}_ids.txt", "".join(row["GenomeID"] + "\n" for row in selected))
        write_tsv(out / f"best_matches_{side}.tsv", PAIR_COLUMNS,
                  [annotate(best_match[side][row["GenomeID"]]) for row in rows[side] if row["GenomeID"] in best_match[side]])
    reciprocal = [pair for pair in best_match["a"].values() if best_match["b"][pair["b"]]["a"] == pair["a"]]
    write_tsv(out / "reciprocal_best_matches.tsv", PAIR_COLUMNS, [annotate(pair) for pair in reciprocal])

    conflicts = 0
    with ExitStack() as stack:
        writers = {}
        for name in ("ambiguous_matches", "taxonomy_conflicts"):
            writers[name] = csv.DictWriter(stack.enter_context((out / f"{name}.tsv").open("w", encoding="utf-8", newline="")),
                                           fieldnames=PAIR_COLUMNS, delimiter="\t", lineterminator="\n")
            writers[name].writeheader()
        handle = stack.enter_context((out / "matched_pairs.tsv").open(encoding="utf-8", newline=""))
        for row in csv.DictReader(handle, delimiter="\t"):
            if matches["a"][row["catalog_a_id"]] > 1 or matches["b"][row["catalog_b_id"]] > 1:
                writers["ambiguous_matches"].writerow(row)
            if row["species_taxonomy_agreement"] == "different" or row["genus_taxonomy_agreement"] == "different":
                conflicts += 1
                writers["taxonomy_conflicts"].writerow(row)

    taxon_rows, label_rows, only_taxa = [], [], []
    for side in ("a", "b"):
        groups = defaultdict(lambda: [0, 0])
        for row in statuses[side]:
            for rank in RANKS:
                group = groups[(rank, row[rank] or "Unclassified", lineage(row, rank))]
                group[0] += 1
                group[1] += row["status"] == "matched"
        for (rank, taxon, lin), (total, matched) in sorted(groups.items(), key=lambda item: (RANKS.index(item[0][0]), item[0][2])):
            taxon_rows.append(dict(catalog_side=side, catalog_label=labels[side], rank=rank, taxon=taxon, lineage=lin,
                                   representatives=total, matched_representatives=matched, unmatched_representatives=total - matched,
                                   overlap_pct=pct(matched, total)))
    for rank in RANKS:
        counts = {side: Counter(lineage(row, rank) for row in rows[side] if row[rank]) for side in ("a", "b")}
        a, b = set(counts["a"]), set(counts["b"])
        label_rows.append(dict(catalog_a_label=args.label_a, catalog_b_label=args.label_b, rank=rank,
                               catalog_a_assigned_taxa=len(a), catalog_b_assigned_taxa=len(b), shared_taxa=len(a & b),
                               catalog_a_only_taxa=len(a - b), catalog_b_only_taxa=len(b - a),
                               catalog_a_unassigned=sum(not row[rank] for row in rows["a"]),
                               catalog_b_unassigned=sum(not row[rank] for row in rows["b"])))
        for side, unique in (("a", a - b), ("b", b - a)):
            for lin in sorted(unique):
                only_taxa.append(dict(catalog_side=side, catalog_label=labels[side], rank=rank, taxon=lin.split(";")[-1],
                                      lineage=lin, representatives=counts[side][lin]))
    write_tsv(out / "overlap_by_taxon.tsv", list(taxon_rows[0]), taxon_rows)
    write_tsv(out / "taxonomy_label_overlap.tsv", list(label_rows[0]), label_rows)
    write_tsv(out / "catalog_only_taxa.tsv", ["catalog_side", "catalog_label", "rank", "taxon", "lineage", "representatives"], only_taxa)

    def counts_summary(a_matched, b_matched, pairs):
        return dict(catalog_a_total=len(rows["a"]), catalog_a_matched=a_matched,
                    catalog_a_unmatched=len(rows["a"]) - a_matched, catalog_a_overlap_pct=pct(a_matched, len(rows["a"])),
                    catalog_b_total=len(rows["b"]), catalog_b_matched=b_matched,
                    catalog_b_unmatched=len(rows["b"]) - b_matched, catalog_b_overlap_pct=pct(b_matched, len(rows["b"])),
                    qualifying_pairs=pairs)

    summary = dict(catalog_a_label=args.label_a, catalog_b_label=args.label_b, taxonomy_release=args.release,
                   ANI_threshold_pct=args.ani, AF_threshold_pct=args.af, AF_mode=args.af_mode,
                   **counts_summary(len(best_match["a"]), len(best_match["b"]), sum(matches["a"].values())),
                   catalog_a_multiple_matches=sum(count > 1 for count in matches["a"].values()),
                   catalog_b_multiple_matches=sum(count > 1 for count in matches["b"].values()),
                   reciprocal_best_pairs=len(reciprocal), taxonomy_conflict_pairs=conflicts)
    write_tsv(out / "summary.tsv", list(summary), [summary])
    sensitivity_rows = []
    for rule, ani_cut, af_cut, mode in rules:
        a, b, pairs = sensitivity[rule]
        sensitivity_rows.append(dict(catalog_a_label=args.label_a, catalog_b_label=args.label_b, rule=rule,
                                     ANI_threshold_pct=ani_cut, AF_threshold_pct=af_cut, AF_mode=mode,
                                     **counts_summary(len(a), len(b), pairs)))
    write_tsv(out / "threshold_sensitivity.tsv", list(sensitivity_rows[0]), sensitivity_rows)
    print(f"{args.label_a}: {len(best_match['a'])}/{len(rows['a'])} matched; "
          f"{args.label_b}: {len(best_match['b'])}/{len(rows['b'])} matched.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="command", required=True)
    p = subs.add_parser("prepare")
    for flag in ("root", "out", "catalog-a", "catalog-b"):
        p.add_argument("--" + flag, required=True)
    p.add_argument("--label-a", default="Catalog_A")
    p.add_argument("--label-b", default="Catalog_B")
    p.set_defaults(func=prepare)
    p = subs.add_parser("check-parameters")
    for flag in ("ani", "af", "screen", "search-af"):
        p.add_argument("--" + flag, type=float, required=True)
    p.add_argument("--label-a", default="Catalog_A")
    p.add_argument("--label-b", default="Catalog_B")

    def check_parameters(args):
        validate_labels(args.label_a, args.label_b)
        for value in (args.ani, args.af, args.screen, args.search_af):
            if not math.isfinite(value) or not 0 < value <= 100:
                fail("ANI/AF parameters must be percentages in (0,100]")
        if args.screen >= min(args.ani, 95):
            fail("SCREEN_ANI must be lower than the primary ANI cutoff and 95")
        if args.search_af > min(args.af, 50):
            fail("SEARCH_MIN_AF must not exceed the primary AF cutoff or 50")
    p.set_defaults(func=check_parameters)
    p = subs.add_parser("check-version")
    p.add_argument("version")

    def check_version(args):
        match = re.search(r"(\d+)\.(\d+)\.(\d+)", args.version)
        if not match or tuple(map(int, match.groups())) < (0, 3, 1):
            fail(f"Use skani >=0.3.1; found {args.version!r}")
    p.set_defaults(func=check_version)
    for command, function in (("check-raw", check_raw), ("audit-log", audit_log)):
        p = subs.add_parser(command)
        p.add_argument("filename")
        p.set_defaults(func=lambda args, f=function: f(args.filename))
    p = subs.add_parser("summarize")
    for flag in ("catalog-a", "catalog-b", "raw", "out"):
        p.add_argument("--" + flag, required=True)
    p.add_argument("--label-a", default="Catalog_A")
    p.add_argument("--label-b", default="Catalog_B")
    p.add_argument("--release", default="unspecified")
    p.add_argument("--ani", type=float, default=95.0)
    p.add_argument("--af", type=float, default=65.0)
    p.add_argument("--af-mode", choices=("max", "both"), default="max")
    p.set_defaults(func=summarize)
    args = parser.parse_args()
    try:
        if args.command == "summarize":
            if any(not math.isfinite(value) or not 0 < value <= 100 for value in (args.ani, args.af)):
                fail("ANI/AF thresholds must be percentages in (0,100]")
        args.func(args)
    except (ValueError, OSError, csv.Error, KeyError, TypeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())


