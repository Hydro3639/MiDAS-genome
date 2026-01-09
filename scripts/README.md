## Scripts

This directory contains general analysis scripts used in the MiDAS genome project.  
The figure-specific scripts located in `analysis/scripts`, the scripts in this directory are **not tied to individual figures** and are intended for exploratory analyses, data processing, and intermediate result generation.

## Example scripts

This section provides brief usage examples for commonly used analysis scripts in this directory.  
For full command-line options, run each script with the `-h` flag.

----------

### `genus_report.py`

Summarize **species-level abundance patterns** for a given genus across global and Danish WWTP samples.

**Use case**  
Use this script to explore how species within a specific genus (e.g. _Accumulibacter_) are distributed across samples and regions.

**Basic usage**
```
./genus_report.py \
  -g Accumulibacter \
  -f ../data/midas.hq.spRep.abun.profile.tsv \
  --global-end 84
```

**Key arguments**

-   `-g / --genus`  
    Genus name to match within the GTDB taxonomy (case-sensitive).
    
-   `-f / --file`  
    Input TSV file containing species-level abundances and taxonomy.
    
-   `--global-end`  
    Last column index (1-based) corresponding to globally collected WWTP samples.  
    Columns after this index are treated as Denmark-collected samples.
    
----------

### `id_set_report.py`

Summarize **species-level abundance patterns** for a user-defined set of species IDs.

**Use case**

Use this script to analyze abundance patterns for:

-   Candidate species for naming
-   Functional guilds
-   Manually curated species sets

**Basic usage**

```
./id_set_report.py \
  -i pao_species_ids.txt \
  -L PAOs \
  -f ../data/midas.hq.spRep.abun.profile.tsv \
  --global-end 84
``` 

**Key arguments**

-   `-i / --id-file`  
    Text file containing species IDs of interest (one ID per line).
    
-   `-L / --label`  
    Label for the ID set, used in the output summary header.
    
-   `-f / --file`  
    Input TSV file containing species-level abundances and taxonomy.
    
-   `--global-end`  
    Last column index (1-based) corresponding to global WWTP samples.

## `rank_species_within_genus.py`

Rank **species within a given genus** by their highest observed abundance across global WWTP samples.

**Use case**  

Use this script to identify the most abundant or most widely detected species within a genus, for example, when:

-   Prioritizing species for naming or curation
-   Comparing dominant species across treatment configurations
-   Selecting representative species for downstream analyses
   
**Basic usage**

```
./rank_species_within_genus.py \
  -i Accumulibacter \
  -f ../data/midas.hq.spRep.abun.profile.tsv \
  -m ../data/wwtp_sample_sim.metadata.tsv \
  --global-end 84
```
**Optional filters**

```
./rank_species_within_genus.py \
  -i Accumulibacter \
  -f ../data/midas.hq.spRep.abun.profile.tsv \
  -m ../data/wwtp_sample_sim.metadata.tsv \
  --process-filter "C,N,DN,P" \
  --detect-cutoff 0.1 \
  --global-end 84
``` 

**Key arguments**

-   `-i / --input-genus`  
    Genus name to query (e.g. `Accumulibacter`).
    
-   `-f / --file`  
    Input TSV file containing species-level abundances and taxonomy.
    
-   `-m / --metadata`  
    Sample metadata table including treatment process information.
    
-   `--global-end`  
    Last column index (1-based) corresponding to global WWTP samples.
    
-   `--process-filter`  
    Restrict analysis to samples with a specific `Process_type`  
    (use `"none"` to disable filtering).
    
-   `--detect-cutoff`  
    Minimum abundance threshold (%) for species detection  
    (default: `0`, i.e. any detected abundance).
    

## `rank_species_by_set.py`

Rank a **user-defined set of species** by their highest observed abundance across global WWTP samples.

**Use case**  

Use this script to rank species belonging to:

-   Functional guilds (e.g. PAOs, GAOs)
-   Candidate species sets for naming
-   Curated taxonomic or ecological groups
   
**Basic usage**

```
./rank_species_by_set.py \
  -S comammox_species_ids.txt \
  -f ../data/wwtp_sample_sim.metadata.tsv \
  -m sample_metadata.tsv \
  --global-end 84
``` 

**Optional filters**

```
./rank_species_by_set.py \
  -S comammox_species_ids.txt \
  -f ../data/midas.hq.spRep.abun.profile.tsv \
  -m ../data/wwtp_sample_sim.metadata.tsv \
  --process-filter "C,N,DN,P" \
  --detect-cutoff 0.1 \
  --global-end 84
``` 

**Key arguments**

-   `-S / --species-set`  
    File containing species IDs of interest (one ID per line).
    
-   `-f / --file`  
    Input TSV file containing species-level abundances and taxonomy.
    
-   `-m / --metadata`  
    Sample metadata table including treatment process information.
    
-   `--global-end`  
    Last column index (1-based) corresponding to global WWTP samples.
    
-   `--process-filter`  
    Restrict analysis to samples with a specific `Process_type`  
    (use `"none"` to disable filtering).
    
-   `--detect-cutoff`  
    Minimum abundance threshold for species detection  
    (default: `0`).
