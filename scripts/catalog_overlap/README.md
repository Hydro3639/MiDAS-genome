# Genome catalog overlap identification workflow

Compare two genome catalogs using skani ANI and aligned fractions. The workflow accepts any two catalogs through a Bash main script and editable configuration. Catalog A is searched against sketches of catalog B. Both catalogs receive matched/unmatched reports, and all qualifying pairs are retained.

## Files and dependencies

The main files are `run_overlap.sh`, `config.sh` and `scripts/overlap_tables.py`. The Python helper uses its standard library. `examples/` contains metadata templates for users to edit.

Run on Linux with Bash 4 or newer, Python 3.7 or newer, skani 0.3.1 or newer, GNU, `awk`, `cmp` and `flock`. Set `SKANI` and `PYTHON` in `config.sh` if the executables use different names or locations. The default thread count is eight.

To select a Conda environment from `config.sh`, set its name separately from the executable.

```bash
SKANI_ENV="catalog-overlap"
SKANI="skani"
```

The main script uses `conda run` for every skani invocation, including version checks, help, sketching and searching. Manual activation is unnecessary, and `conda` must be available in the launching shell. `SKANI_ENV="PATH"` is the default and uses the active shell's executable. An environment prefix can also be supplied; relative prefixes containing `/` use `PROJECT_ROOT`. `PYTHON` continues to select the Python executable for the table-processing helper. See the [Conda run documentation](https://docs.conda.io/projects/conda/en/stable/commands/run.html).

## Quick start

Supply two metadata tables and point the main script to them.

```bash
bash run_overlap.sh run \
  --catalog-a metadata/catalog_a.tsv \
  --catalog-b metadata/catalog_b.tsv \
  --root . \
  --output overlap_results \
  --threads 8
```

Optional display labels appear in the reports. Output filenames and column names stay stable when labels change.

```bash
bash run_overlap.sh run \
  --catalog-a metadata/catalog_a.tsv \
  --catalog-b metadata/catalog_b.tsv \
  --root . \
  --label-a "Collection A" \
  --label-b "Collection B" \
  --output overlap_results
```

Alternatively, fill in `config.sh`, or copy it to `config.local.sh`.

```bash
CATALOG_A_METADATA="metadata/catalog_a.tsv"
CATALOG_B_METADATA="metadata/catalog_b.tsv"
LABEL_A="Catalog_A"
LABEL_B="Catalog_B"
OUT_DIR="results"
```

Then run `bash run_overlap.sh run`, or `bash run_overlap.sh run --config config.local.sh`. Run `bash run_overlap.sh --help` for all main-script options. Command-line options override the configuration values.

## Input format and paths

Each metadata table is tab-separated with the exact header `GenomeID` and `path`. A third column named `taxonomy` is optional.

| Column | Content |
| --- | --- |
| `GenomeID` | Unique genome identifier within that catalog |
| `path` | Readable nucleotide FASTA path, including gzip FASTAs |
| `taxonomy` | Optional rank-prefixed lineage, for example `d__Bacteria;p__Example;g__Example;s__Example species` |

IDs and resolved FASTA paths must be unique within each catalog. IDs can occur in both catalogs. Both inputs must contain at least one genome. Taxonomy can be omitted from either table or left blank. When provided, taxonomy uses the `d__`, `p__`, `c__`, `o__`, `f__`, `g__` and `s__` rank prefixes for the optional taxonomy reports.

The example metadata files are templates. Replace their paths with the locations of your own FASTAs before running the analysis.

## Matching rules

The default rule requires ANI ≥95% and aligned fraction ≥65% in at least one direction. The directional AFs describe the fractions aligned in catalog A and catalog B for the same comparison.

| Setting | Meaning |
| --- | --- |
| `ANI_THRESHOLD=95` or `--ani 95` | Minimum ANI percentage, inclusive |
| `AF_THRESHOLD=65` or `--af 65` | Minimum aligned fraction percentage, inclusive |
| `AF_MODE=max` or `--af-mode max` | At least one directional AF reaches the threshold |
| `AF_MODE=both` or `--af-mode both` | Both directional AFs reach the threshold |

Search settings are broader than the final matching rule. `SCREEN_ANI` defaults to 80 and `SEARCH_MIN_AF` to 15. The workflow checks that screening also covers the included sensitivity summaries at ANI95/maxAF50, ANI95/maxAF65 and ANI95/bothAF65. Adjust screening in the configuration when using lower primary thresholds.

Each input row is a catalog unit. To compare species representation, supply one representative per species from each catalog. The workflow performs cross-catalog comparison and does not dereplicate the supplied catalogs internally. It reports candidates returned by skani under the configured screening settings. An unmatched genome has no qualifying reported counterpart in the supplied other catalog under the selected rule.

All qualifying pairs determine matched counts. A genome can match multiple genomes in the other catalog, so matched counts for A and B can differ. Best matches are selected by highest ANI, then highest AF under the selected mode, then highest minimum directional AF, then the other genome ID. Best-match tables can reuse the same counterpart. Reciprocal best matches select pairs that choose each other under this ranking.

Taxonomy labels do not determine sequence overlap. Optional taxonomy reports describe label agreement, conflicts and taxonomic groups represented among the matched/unmatched genomes.

## Outputs

Reports are written in `OUT_DIR/04_summary/`.

| File | Content |
| --- | --- |
| `summary.tsv` | Total, matched and unmatched counts and overlap percentages for each catalog |
| `catalog_a_status.tsv`, `catalog_b_status.tsv` | Status, candidate counts, reasons and best counterparts for every genome |
| `catalog_a_matched.tsv`, `catalog_b_matched.tsv` | Genomes with at least one qualifying counterpart |
| `catalog_a_unmatched.tsv`, `catalog_b_unmatched.tsv` | Genomes with no qualifying counterpart |
| `catalog_*_matched_ids.txt`, `catalog_*_unmatched_ids.txt` | Corresponding genome ID lists |
| `matched_pairs.tsv` | Every qualifying pair with ANI, both AFs and optional taxonomy annotations |
| `best_matches_a.tsv`, `best_matches_b.tsv` | One best qualifying pair for each matched genome on the indicated side |
| `reciprocal_best_matches.tsv` | Pairs selected as best from both sides |
| `ambiguous_matches.tsv` | Qualifying pairs where either genome has multiple partners |
| `boundary_review.tsv` | Reported pairs near the primary ANI/AF thresholds, with their pass/fail flags |
| `threshold_sensitivity.tsv` | Counts under the primary rule and three additional ANI95/AF rules |
| `overlap_by_taxon.tsv` | Matched/unmatched genome counts by rank and lineage for both catalogs |
| `taxonomy_label_overlap.tsv` | Shared and catalog-specific assigned taxonomy labels by rank |
| `catalog_only_taxa.tsv` | Assigned taxonomy labels occurring in only one input catalog |
| `taxonomy_conflicts.tsv` | Qualifying pairs with different assigned genus or species labels |

In pair tables, `catalog_a_id` and `catalog_b_id` identify the two genomes. `AF_catalog_a` is skani's query AF, and `AF_catalog_b` is its reference AF. Status tables use `GenomeID` for the genome on the indicated side and `best_match_other_id` for its best qualifying counterpart. Unmatched reasons distinguish no reported hit, ANI below threshold and high ANI with insufficient AF. Boundary review includes ANI within 0.25 percentage points of the cutoff with AF ≥min(50, primary AF), or passing ANI with AF within two percentage points of its cutoff.

Other output folders contain normalized inputs (`01_input/`), reference sketches (`02_sketch/`), the raw search (`03_search/`), logs (`logs/`) and input/settings records (`record/`).

## Stages and reuse

Available stages are `prepare`, `sketch`, `search`, `summarize` and `run`. Supply the same configuration or input arguments for stage-by-stage execution. `prepare` checks metadata, readable/nonempty files and the beginning of each FASTA. It can run before skani is installed. `run` executes all stages.

Repeated runs reuse compatible sketches and raw searches. Reuse compares plain records of FASTA paths, sizes, modification times and relevant settings. Use a new output directory when changing catalog B genomes or sketch settings. Changes to catalog A genomes or search settings cause a new search. Final thresholds, taxonomy and display labels can be changed while reusing a compatible raw search.
