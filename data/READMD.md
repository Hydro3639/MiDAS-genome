## Description of datasets

The complete datasets generated in this study are provided in the Supplementary files.  
Input files used directly for figure generation are available in the `analysis/sources` directory.

This directory contains curated, general-purpose datasets produced as part of the MiDAS genome project, intended to support downstream analyses and data reuse.

### `midas.hq.spRep.abun.profile.tsv`

**Description**  

Abundance profiles of MiDAS high-quality (HQ) species representatives across activated sludge samples.

- **Data type**: abundance profile based on sylph
- **Resolution**: HQ species representative (spRep) of the MiDAS genome catalog
- **Sample coverage**:  
  - 83 globally collected WWTP samples (one sample per WWTP)  
  - [69 Denmark-collected samples](https://www.nature.com/articles/s41467-021-22203-2) (three samples per WWTP)
- **Value definition**: relative abundance per sample (%)

### `wwtp_sample_sim.metadata.tsv`

**Description**  

Metadata table describing wastewater treatment plant (WWTP) samples included in the MiDAS project.

-   **Data type**: sample metadata (simplified version)
-   **Sample coverage**: global WWTP samples included in MiDAS analyses
-   **Value definition**:
    -   `WWTP_id`: unique identifier of the wastewater treatment plant
    -   `Continent`: continent of sample origin
    -   `Country/Region`: country or region of the WWTP
    -   `Process_type`: treatment process configuration (e.g. `C`, `C,N`, `C,N,DN`, `C,N,DN,P`)

----------

### `midas.hq.mags.gtdb.r226.tsv`

**Description** 

GTDB taxonomy assignments for MiDAS **high-quality (HQ) MAGs**, classified against GTDB release r226.

-   **Data type**: taxonomy annotation table
-   **Resolution**: metagenome-assembled genome (MAG)
-   **Sample coverage**: all MiDAS HQ MAGs
-   **Value definition**:
    -   `user_genome`: MiDAS MAG identifier
    -   `classification`: GTDB taxonomic string from domain to species
----------

### `midas.mq.mags.gtdb.r226.tsv`

**Description**  

GTDB taxonomy assignments for MiDAS **medium-quality (MQ) MAGs**, classified against GTDB release r226.
-   **Data type**: taxonomy annotation table
-   **Resolution**: metagenome-assembled genome (MAG)
-   **Sample coverage**: all MiDAS MQ MAGs
-   **Value definition**:
    -   `user_genome`: MiDAS MAG identifier
    -   `classification`: GTDB taxonomic string from domain to species
