# Upstream pipeline

The steps below produced the derived tables in `../data/` from the raw reads
deposited in the European Nucleotide Archive under accession PRJEB109287. Tool
versions and all non-default parameters are given as reported in the manuscript
Methods. The original job scripts are not yet included in this repository; they
are available from the corresponding author on request and will be added here.

## 1. Read quality filtering

KneadData 0.12.0, default settings, which trim adapters and low-quality bases with
Trimmomatic. No reference database was used for read removal.

## 2. Per-site co-assembly

MEGAHIT 1.2.9. For each site, the quality-filtered reads of the aerosol metagenome
and its two paired source metagenomes were co-assembled together, with k-mer sizes
from 21 to 121 and a minimum contig length of 1,500 bp.

## 3. Mapping and binning

Reads from each sample were mapped back to the co-assembly of their own site with
Bowtie2 2.5.4. Contigs of at least 1,500 bp were binned with MetaBAT 2 v2.17 using
default parameters.

## 4. Genome quality and dereplication

Completeness and contamination were estimated with CheckM 1.1.3. Bins with
completeness of at least 50% and contamination below 10% were retained as MAGs;
those above 90% completeness and below 5% contamination are referred to as high
quality. MAGs were dereplicated across sites with dRep 3.6.1 at 95% average
nucleotide identity, with the default 75% completeness filter.

## 5. Coverage

Quality-filtered reads from all 17 metagenomes, including the two field controls,
were mapped to the 532 dereplicated MAGs with minimap2 2.28 in short-read mode.
Mean coverage per contig was calculated with CoverM using default settings and
averaged per MAG, weighting contigs by length.

## 6. Taxonomy and closest reference

GTDB-Tk 2.4.1 against GTDB release 226. Reference genome metadata were retrieved
from NCBI BioSample records using NCBI Datasets.

## 7. Gene prediction and annotation

Prodigal 2.6.3 in single-genome mode, then eggNOG-mapper 2.1.12 against the
eggNOG 5.0 database.

## 8. Phylogeny

FastTree 2.1.11 on the GTDB-Tk alignment of 120 bacterial marker genes, under the
WAG model with gamma-distributed rate variation, rooted at its midpoint.

## 9. Population genomics

inStrain 1.12.0 profile and compare, run on the Dairy Farm aerosol and the two Dairy
Farm soil metagenomes against the 233 Dairy Farm MAGs, using the minimap2 alignments
from step 5. Because those alignments did not record read pairing, all aligned reads
were used rather than only properly paired reads. Profiles were also generated with
reads rarefied to 5x and 10x per position.
