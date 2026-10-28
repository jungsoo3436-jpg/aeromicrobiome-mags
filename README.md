# aeromicrobiome-mags

Code and derived data for:

**Park J, Fowler SJ.** Genome-resolved comparison shows that apparent selection in
near-surface bioaerosols reflects lineage, not survival traits. *(manuscript in preparation)*

This repository reproduces every figure, table and statistic in the paper from the
derived tables in `data/`. Raw sequencing reads are deposited separately in the
European Nucleotide Archive under project accession **PRJEB109287**.

## What the study did

Aerosols were collected 1.5 m above ground at five sites in Metro Vancouver, Canada
(three wetlands, a dairy farm and a landfill), alongside two soil samples per site, or
two water samples at Burns Bog. Each aerosol metagenome was co-assembled with its two
paired source metagenomes, yielding 525 dereplicated bacterial metagenome-assembled
genomes (MAGs). Each MAG was then compared with its own source environment at three
levels: depth-normalised abundance, survival-associated gene content, and population
structure.

## Repository layout

```
analysis/
  MAG_air_trait_analysis.R   single script; reproduces all figures and tables
data/
  README.md                  provenance and contents of every input file
  ...                        derived input tables (see data/README.md)
hpc/
  README.md                  upstream pipeline steps, tools and parameters
figures/
  Fig1-4, FigS1-4            published figures as produced by the script
```

## Running the analysis

Requires **R 4.4.0** or later and the following packages:

```r
install.packages(c("data.table", "dplyr", "tidyr", "stringr", "ggplot2", "scales",
                   "patchwork", "lme4", "lmerTest", "ape", "phytools"))
```

Then, from the repository root:

```bash
AERO_BASE=$(pwd) Rscript analysis/MAG_air_trait_analysis.R
```

or set `BASE` in Section 0 of the script and source it in RStudio.

Outputs are written to `mag_based/corrected_analysis/` relative to `AERO_BASE`,
and comprise the eight figure PDFs plus a `.tsv` for every statistic reported in
the paper. Runtime is a few minutes on a laptop; the phylogenetic and population
genomics sections are the slowest steps.

### Note on paths

The script was written against the original working tree, in which the derived
tables sat one directory above a `mag_based/` subfolder. Inputs in this repository
have been collected into `data/`, so the file paths in Section 0 of the script
(`F_COV`, `F_META`, `F_KO`, `F_READS`, `DIR_IS` and the tree and alignment files)
must be pointed at `data/` before the script will run unmodified.

## Software used

| Step | Tool | Version |
|---|---|---|
| Read quality filtering | KneadData (Trimmomatic) | 0.12.0 |
| Co-assembly | MEGAHIT | 1.2.9 |
| Mapping for binning | Bowtie2 | 2.5.4 |
| Binning | MetaBAT 2 | 2.17 |
| Genome quality | CheckM | 1.1.3 |
| Dereplication | dRep | 3.6.1 |
| Coverage mapping | minimap2 | 2.28 |
| Coverage calculation | CoverM | as cited |
| Taxonomy | GTDB-Tk (release 226) | 2.4.1 |
| Gene prediction | Prodigal | 2.6.3 |
| Annotation | eggNOG-mapper (eggNOG 5.0) | 2.1.12 |
| Population genomics | inStrain | 1.12.0 |
| Phylogeny | FastTree | 2.1.11 |
| Statistics and figures | R | 4.4.0 |

## Licence

Code in `analysis/` is released under the MIT Licence (`LICENSE`).
Derived data in `data/` are released under CC BY 4.0 (`LICENSE-data`).

## Contact

Jungsoo Park, jungsoop@sfu.ca
S. Jane Fowler, sjfowler@sfu.ca
Department of Biological Sciences, Simon Fraser University, Burnaby, Canada
