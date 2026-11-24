This modular genome quality classification workflow starts from genome FASTAs and combines CheckM2 quality estimates with RNA annotations. GTDB-Tk domain assignment determines the Barrnap and tRNAscan-SE models throughout the workflow. The Bash main script provides the stage overview and argument handling. Bash modules run the individual tools, and Python helpers merge reports, parse GTDB domains, summarize RNA calls and assign quality classes.

The main files are `genome_qc_pipeline.sh`, `config.sh`, `lib/common.sh`, the stage scripts in `modules/`, and four helpers in `scripts/`.

| Python helper | Processing |
| --- | --- |
| `merge_checkm2.py` | Merge CheckM2 chunk reports using the manifest genome IDs |
| `parse_gtdb_domains.py` | Build a bacterial/archaeal domain map from GTDB-Tk summaries |
| `summarize_rna.py` | Summarize Bakta, Barrnap and tRNAscan-SE outputs and export RNA coordinates |
| `aggregate.py` | Combine quality/RNA measurements and write one genome quality classification |

Dependencies are Bash 4 or newer, Python 3.7 or newer with its standard library, GNU coreutils, `cmp`, `flock`, CheckM2, GTDB-Tk, Bakta, Barrnap and tRNAscan-SE. CheckM2, GTDB-Tk and Bakta require their reference databases. Each software environment can be an active `PATH`, a Conda environment name or a Conda prefix.

Fill in these settings in `config.sh`.

```bash
GENOMES_DIR="genomes"
OUTDIR="results"
CHECKM2_DB="databases/checkm2/uniref100.KO.1.dmnd"
GTDBTK_DB="databases/gtdbtk"
BAKTA_DB="databases/bakta"
```

CheckM2 expects the database file, and GTDB-Tk and Bakta expect database directories. `GTDBTK_DATA_PATH` can also supply the GTDB-Tk database default.

The genome folder accepts `.fa`, `.fna`, `.fasta` and gzip equivalents, including file symlinks. Each FASTA file represents one genome. IDs are derived from filenames by removing the terminal FASTA suffix and normalizing unsupported characters to underscores. Resulting IDs must be unique. CheckM2 report IDs are matched to this manifest, preserving IDs that themselves end in a FASTA-like suffix.

Tools are located on the active `PATH` by default. If separate Conda environments are used, configure their names or prefixes.

```bash
CHECKM2_ENV="checkm2"
GTDBTK_ENV="gtdbtk"
BAKTA_ENV="bakta"
BARRNAP_ENV="barrnap"
TRNASCAN_ENV="trnascan"
```

Adjust threads, concurrent genome jobs, databases, optional GTDB-Tk arguments and annotation settings in `config.sh` for the execution environment. The defaults allocate eight threads to CheckM2/GTDB-Tk and four concurrent RNA annotation jobs with one thread per task.

Run the full workflow after configuration.

```bash
bash genome_qc_pipeline.sh full
```

A separate configuration file can be selected.

```bash
bash genome_qc_pipeline.sh full --config config.local.sh
```

The main input and resource settings can also be supplied through arguments.

```bash
bash genome_qc_pipeline.sh full \
  --genomes genomes \
  --output results \
  --checkm2-db databases/checkm2/uniref100.KO.1.dmnd \
  --gtdbtk-db databases/gtdbtk \
  --bakta-db databases/bakta \
  --threads 8 --jobs 4
```

Command-line paths use the launch directory. Command-line values override the corresponding configuration defaults. Run `bash genome_qc_pipeline.sh --help` for the full main-script interface.

Stages can be run separately in dependency order.

```bash
bash genome_qc_pipeline.sh prepare
bash genome_qc_pipeline.sh checkm2
bash genome_qc_pipeline.sh gtdbtk
bash genome_qc_pipeline.sh domains
bash genome_qc_pipeline.sh bakta
bash genome_qc_pipeline.sh barrnap
bash genome_qc_pipeline.sh trnascan
bash genome_qc_pipeline.sh summarize
bash genome_qc_pipeline.sh aggregate
bash genome_qc_pipeline.sh status
```

The `rna` action runs Bakta, Barrnap, tRNAscan-SE and RNA summarization. Run `gtdbtk` and `domains` before executing the RNA stages separately. A GTDB domain map is required for Barrnap, tRNAscan-SE and RNA summarization.

RNA model selection always follows the GTDB-Tk domain assignment. Missing or conflicting domain assignments are recorded as `Unknown`.

| GTDB domain | Barrnap model | tRNAscan-SE model |
| --- | --- | --- |
| Bacteria | `bac` | `-B` |
| Archaea | `arc` | `-A` |
| Unknown or unassigned | Both `bac` and `arc` | Both `-B` and `-A` |

For unknown domains, rRNA calls from both model sets contribute to presence detection, and the tRNAscan-SE result covering the most standard amino-acid types is selected. Bakta is run on every genome. `BAKTA_RNA_ONLY=1` enables RNA-focused Bakta annotation, and `BAKTA_RNA_ONLY=0` runs full Bakta.

The main results are written below `OUTDIR`.

| Output | Content |
| --- | --- |
| `work/genomes.tsv` | Input genome IDs and FASTA paths |
| `work/gtdb_domain.tsv` | Domain map used for RNA model selection |
| `results/checkm2/checkm2_summary.tsv` | CheckM2 completeness and contamination |
| `results/rna/rna_summary.tsv` | Per-genome RNA measurements |
| `results/rna/rrna_hits.tsv` | Raw RNA coordinates and diagnostics |
| `results/final/all_genomes_quality.tsv` | GTDB domain, selected RNA models, measurements and one `bin_status` per genome |
| `results/final/status_counts.tsv` | One row per status with `bin_status` and `n_genomes` |

Default HQ candidates have CheckM2 completeness ≥90% and contamination ≤5%. HQ also requires 5S, 16S and 23S rRNA, with each gene type detected by Bakta or Barrnap, and tRNAs covering at least 18 of the 20 standard amino acids according to either Bakta or the selected tRNAscan-SE model. Each predictor is evaluated separately. A single detected tRNA does not satisfy this requirement.

MQ requires completeness ≥50% and contamination ≤10% after HQ assignment. LQ has completeness below the MQ threshold and contamination within the MQ limit. Higher contamination is labeled `Contaminated`. Missing or invalid CheckM2 measurements are labeled `UNASSESSED`. Thresholds and `MIN_UNIQUE_TRNA` can be configured, with the tRNA threshold restricted to 1–20 amino-acid types.

The RNA summary and final table use the same tRNA field names and one classification rule.

| Field | Definition |
| --- | --- |
| `gtdb_domain` | GTDB-Tk domain assignment used for model selection |
| `barrnap_model_used`, `trnascan_model_used` | RNA models contributing to the summary |
| `bakta_trna_unique_standard20` | Standard amino-acid types represented by Bakta tRNA calls |
| `trnascan_trna_unique_standard20` | Standard amino-acid types represented by the selected tRNAscan-SE result |
| `trna_complete` | Either predictor reaches `MIN_UNIQUE_TRNA` |
| `rrna_complete` | 5S, 16S and 23S types are all detected |
| `bin_status` | Single HQ, MQ, LQ, Contaminated or UNASSESSED classification in the final table |
| `hq_fail_reasons` | Reasons HQ requirements were not met, using the configured thresholds |

Bakta gene/product labels and tRNAscan-SE isotypes are normalized to the same 20-amino-acid set. Initiator methionine is counted as methionine, and Ile2 as isoleucine. Unknown, suppressor and nonstandard amino-acid types do not increase this count. Raw tRNA totals and label counts remain available as diagnostics. Bakta label handling follows its [tRNA annotation implementation](https://github.com/oschwengers/bakta/blob/main/bakta/features/t_rna.py).

RNA hit diagnostics include partial/truncated predictions and contig-edge information. rRNA presence includes reported partial predictions, so `rrna_complete` records detection of all three types rather than proving full-length genes or an intact operon.

Successful steps are reused when the same configured analysis is repeated. Plain input records compare genome file paths, sizes, modification times and annotation settings. Sequence checksums are not generated. Failed RNA tasks are recorded in `work/*_failed.txt` and can be retried through their stage command. Use the same configuration file and arguments for a resumed run. Assign a separate `OUTDIR` when changing genome inputs, database paths or annotation settings, and rerun `prepare` before continuing a changed analysis. Final threshold changes can be applied by rerunning `summarize` and `aggregate`.

For results from an earlier release, rerun `summarize` and then `aggregate` with the same configuration and arguments to rebuild the current RNA and final-table schemas. The aggregation helper checks its required input columns and rejects outdated RNA tables. Existing raw annotation files can be reused.
