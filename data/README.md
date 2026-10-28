# Derived input data

Every file here is a derived table produced by the upstream pipeline described in
`../hpc/README.md`. Raw reads are in the European Nucleotide Archive under
accession PRJEB109287.

Sample codes: `A_` aerosol, `S_` soil, `W_` water, `AC_` field control.
Site codes: `BB` Burns Bog, `DF` Dairy Farm, `MK` Minnekhada, `RM` Richmond Marsh,
`VL` Vancouver Landfill.

| File | Produced by | Contents |
|---|---|---|
| `site_derep_all894_mean_cov.REBUILT.tsv` | CoverM, contig mode | Length-weighted mean coverage of each dereplicated MAG in every metagenome. One row per MAG, one column per sample. |
| `site_coassemble_HQ_MQ_checkm_gtdb.v2_summery.csv` | CheckM 1.1.3 and GTDB-Tk 2.4.1 | Per-bin completeness, contamination, genome size, N50, GC content, GTDB classification, closest reference accession and its ANI. |
| `read_counts_per_sample.tsv` | samtools flagstat | Total and primary mapped reads per sample, used to express coverage per 100 million quality-filtered reads. |
| `eggnog_allMAGs.KO_long.tsv` | eggNOG-mapper 2.1.12 | Long-format MAG to KEGG Orthology assignments, one row per MAG and KO. |
| `kegg_ko_list.tsv` | KEGG | KO identifiers and definitions, used to verify that the 47 a priori marker KOs correspond to the intended gene symbols. |
| `closest_ref_summary.jsonl` | NCBI Datasets | One JSON record per closest reference genome, carrying the BioSample isolation source, host, environmental context and metagenome source fields used for habitat assignment. |
| `bac120_525_fasttree.nwk` | FastTree 2.1.11 | Maximum-likelihood tree of the 525 bacterial MAGs from the GTDB-Tk alignment of 120 bacterial marker genes, WAG model with gamma rate variation. |
| `gtdbtk_bac120_msa_HQ.fasta.gz`, `gtdbtk_bac120_msa_MQ.fasta.gz` | GTDB-Tk 2.4.1 | Concatenated marker gene alignments for the high and medium quality MAG sets. |
| `controls/genome_mapped_reads_17samples.tsv` | minimap2 and CoverM | Reads mapped per MAG in all 17 metagenomes, including the two field controls, used for the control comparison in Supplementary Figure S1. |
| `controls/AC_DF2.flagstat.txt`, `controls/AC_RM2.flagstat.txt` | samtools flagstat | Mapping statistics for the two field control libraries. |
| `instrain_DF/*_genome_info.tsv` | inStrain 1.12.0 profile | Per-MAG coverage, breadth and nucleotide diversity for the Dairy Farm aerosol and two soil metagenomes. |
| `instrain_DF/rar5/`, `instrain_DF/rar10/` | inStrain 1.12.0 profile | As above, with reads rarefied to 5x and 10x per position. |
| `instrain_DF/compare_genomeWide_compare.tsv` | inStrain 1.12.0 compare | Pairwise popANI and the fraction of each MAG compared, for aerosol versus soil and soil versus soil. |

## Note on one filename

The analysis script reads `kegg_ko_list.tsv`. In the original working tree this file
was named `kegg_ko_list.txt`, and because the read is guarded by `file.exists()` the
marker gene verification step was skipped silently. The file is supplied here with the
extension the script expects, so that step now runs.
