## Intro
This repository contains scripts developed for the [**MiDAS**](https://midasfieldguide.org/guide) genome project. The project aims to recover [MIMAG-standard](https://www.nature.com/articles/nbt.3893/tables/1) high-quality (HQ) metagenome-assembled genomes (MAGs) for the core taxa from activated sludge (AS) samples collected from wastewater treatment plants (WWTPs) worldwide, building on the global survey conducted in the [**MiDAS 4**](https://www.nature.com/articles/s41467-022-29438-7) study.

In total, 83 AS samples were subjected to Nanopore-based ultra-deep long-read sequencing and included in the analysis. All samples were obtained in compliance with the [**Nagoya Protocol**](https://www.cbd.int/abs/default.shtml) or with approvals from relevant local authorities in regions/countries that are not parties to the Protocol. We gratefully acknowledge our collaborators for their efforts in securing the necessary permissions and enabling the generation of this global dataset.



## Structure

| Folder | Content |
| --- | --- |
| `scripts/` | General analysis scripts and reusable computational workflows that are **not figure-specific**. |
| [<code>└──&nbsp;genome_quality_classification/</code>](scripts/genome_quality_classification/) | Modular workflow for genome quality assessment and classification based on CheckM2 quality estimates, GTDB-Tk domain assignments, and rRNA and tRNA annotations from Bakta, Barrnap, and tRNAscan-SE. |
| `analysis/` | Figure-oriented analysis layer used to reproduce manuscript results. |
| <code>├──&nbsp;scripts/</code> | R scripts used to generate individual manuscript figures. |
| <code>├──&nbsp;sources/</code> | Processed input data used directly for figure generation. |
| <code>└──&nbsp;figures/</code> | Figures generated for the manuscript. |
| `data/` | Lightweight metadata and reference tables used as inputs across analyses. |
| `envs/` | Documentation of the software environments used in this project, including package versions and environment setup instructions to support reproducibility. |
| `README.md` | Overview of the project, repository structure, and analysis workflow. |
| `LICENSE` | The [license](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-a-repository/licensing-a-repository) for this repository. |
