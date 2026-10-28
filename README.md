# aeromicrobiome-mags

Code and derived data for:

**Park J, Fowler SJ.** Genome-resolved comparison with paired source environments
shows that lineage, not survival traits, distinguishes airborne bacteria.
*(manuscript in preparation)*

This repository reproduces every figure, table and statistic in the paper from the
derived tables in `data/`. The raw sequencing reads and the MAGs themselves are
deposited together in the European Nucleotide Archive under project accession
**PRJEB109287**.

## Companion repository

This is the genome-resolved half of a two-paper study. The read-based community
and gene-centric analyses of the same aerosol and paired source metagenomes are in
[aeromicrobiome-readbased](https://github.com/jungsoo3436-jpg/aeromicrobiome-readbased),
which covers Kraken2/Bracken taxonomic profiling, diversity and ordination,
site-blocked DESeq2 differential abundance, and curated marker gene profiling
against a non-redundant gene catalog. Both repositories draw on the same
ENA deposit, PRJEB109287, which holds both the raw reads and the MAGs.

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
  kraken_reports/            gzipped Kraken2 reports, one per library
  ...                        derived input tables (see data/README.md)
hpc/
  README.md                  upstream pipeline steps, tools and parameters
figures/
  Fig1-4, FigS1-4            published figures as produced by the script
results/
  Fig*, FigS*                the same figures, as written by a re-run
  *.tsv                      every statistic reported in the paper
  sessionInfo.txt            R and package versions of that run
```

## Running the analysis

Requires **R 4.4.0** or later and the following packages:

```r
install.packages(c("data.table", "dplyr", "tidyr", "stringr", "ggplot2", "scales",
                   "patchwork", "lme4", "lmerTest", "jsonlite", "ape", "phytools"))
```

Then, from the repository root:

```bash
Rscript analysis/MAG_air_trait_analysis.R
```

or set the working directory to the repository root and source the script in
RStudio. Every path in Section 0 is relative to the repository root, so nothing
needs editing. Set the `AERO_BASE` environment variable only if you need to run
against a copy of the inputs held elsewhere.

Outputs are written to `results/`: the eight figure PDFs, a `.tsv` for every
statistic reported in the paper, and `sessionInfo.txt`. Runtime is a few minutes
on a laptop; the phylogenetic and population genomics sections are the slowest
steps. The script runs end to end and exits cleanly; the committed `results/`
is what it produced, on macOS 27.0.1 under R 4.5.3 (see
`results/sessionInfo.txt` for every package version).

The same run on Linux x86_64 under R 4.3.3 reproduced every statistic to the
precision reported in the paper, with the bootstrap tables byte-identical.
Residual differences are confined to the tenth significant digit and beyond,
and come from the iterative mixed-model optimiser and platform BLAS.

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
