## ============================================================================
## Air association of metagenome-assembled genomes (MAGs): source or survival?
## Park and Fowler, analysis script for the ISME Journal manuscript
## ----------------------------------------------------------------------------
## One script reproduces every number, figure and table in the paper from the
## input files listed in SECTION 0.
##
## Air association
##   Coverage (CoverM, minimap2 alignments) is divided by the quality-filtered
##   reads in each metagenome and expressed per 100 million reads. Each MAG is
##   scored at the site whose co-assembly produced it as
##   log2[(aerosol + 1e-4) / (mean of the two source metagenomes + 1e-4)].
##   >= +1 is air-enriched, <= -1 air-depleted, otherwise intermediate.
##   Sites whose aerosol metagenome has < 1% of reads mapped to the MAGs (MK and
##   RM) are excluded from all air association analyses.
##
## Outputs (written to results/)                       Manuscript item
##   Fig1_air_vs_source.pdf                           Figure 1
##   Fig2_phylum_composition.pdf                      Figure 2
##   Fig3_trait_odds_ratios.pdf                       Figure 3
##   Fig4_lineage_and_reference_habitat.pdf           Figure 4
##   FigS1_field_controls.pdf                         Supplementary Figure S1
##   FigS2_genome_novelty.pdf                         Supplementary Figure S2
##   FigS3_cutoff_sensitivity.pdf                     Supplementary Figure S3
##   FigS4_population_genomics.pdf                    Supplementary Figure S4
##   RQ3_closest_reference_habitat.tsv                Supplementary Table S1
##   SuppTable_marker_KO_verification.tsv             Supplementary Table S2
##   SuppTable_read_recruitment_and_composition.tsv   Table 1 values
##   00_QC_*, *.tsv                                   QC checks and every statistic
##                                                    behind every figure
##   sessionInfo.txt                                  R and package versions
##
## How to run
##   Run from the repository root; every path below is relative to it and
##   needs no editing. From a shell:  Rscript analysis/MAG_air_trait_analysis.R
##   In RStudio: set the working directory to the repository root and Source.
##   Required packages: data.table, dplyr, tidyr, stringr, ggplot2, scales,
##   patchwork, lme4, lmerTest, jsonlite, ape, phytools.
## ============================================================================
## SECTION 0. SETTINGS  (the only section you normally need to edit)
## ============================================================================

## Repository root. Defaults to the current directory, so the script runs as-is
## from the repository root; AERO_BASE overrides it if you need to run elsewhere.
BASE <- Sys.getenv("AERO_BASE", ".")

## Input files (relative to BASE)
F_COV   <- "data/site_derep_all894_mean_cov.REBUILT.tsv"           # corrected coverage
F_META  <- "data/site_coassemble_HQ_MQ_checkm_gtdb.v2_summery.csv" # CheckM + GTDB-Tk, 894 bins
F_KO    <- "data/eggnog_allMAGs.KO_long.tsv"                       # MAG_id <tab> KO
F_READS <- "data/read_counts_per_sample.tsv"  # samtools flagstat on the CoverM BAMs:
                                              # sample, primary, primary_mapped

## Output folder (created if missing). The eight published figure PDFs are also
## kept in figures/; a re-run writes the complete set, every statistic table and
## sessionInfo.txt here.
OUTDIR  <- "results"

## Air association = within-site log2 ratio of depth-normalised coverage:
##   own-site air sample / mean of own-site soil or water samples.
## Primary cutoff: at least 2-fold higher (enriched) or lower (depleted).
LFC_CUT     <- 1
## Cutoffs for the sensitivity analysis (Supplementary Figure S3): 1.5-, 2-, 3-, 4-, 9-fold.
## As an air share of (air + paired source) these are 60, 67, 75, 80 and 90%;
## 9-fold matches the old ">=90% air" rule.
LFC_CUTOFFS <- c(log2(1.5), 1, log2(3), 2, log2(9))
## Pseudocount added to normalised coverage before the log ratio
## (about the 5th percentile of non-zero normalised coverage)
PSEUDO      <- 1e-4
## Sites whose air sample has fewer than this % of reads mapping to the MAGs
## are reported but excluded from the enriched-vs-depleted tests, because
## every genome there sits at the detection floor (MK and RM: <0.1%).
MIN_AIR_MAPPED_PCT <- 1


## Bootstrap / permutation settings
B_BOOT <- 5000
set.seed(1)

## Figure text size. Axis text, axis titles, legends and facet labels all
## scale from this, so raising it enlarges every figure consistently.
## ggplot derives axis/legend text as 0.8 x base size.
BASE_FONT <- 15

## Colours (unchanged from the manuscript figures)
COL_ENR <- "#009E73"
COL_DEP <- "#7A7A7A"
COL_INT <- "#BDBDBD"
ENV_COLS <- c(Air = "#009E73", Soil = "#E69F00", Water = "#0072B2")
## Phylum colours shared by Figures 2 and 4; all other phyla are grey "Other"
PHY_COLS <- c(Bacteroidota = "#E69F00", Pseudomonadota = "#56B4E9",
              Actinomycetota = "#CC79A7", Bacillota = "#0072B2",
              Myxococcota = "#D55E00", Chlamydiota = "#F0E442",
              Planctomycetota = "#999933", Verrucomicrobiota = "#882255",
              Other = "grey80")
## Habitat colours for Figure 4 (ColorBrewer Set1; saturated, so they read
## differently from the pastel phylum colours)
HAB_COLS <- c("Animal-associated" = "#E41A1C", "Engineered/waste" = "#984EA3",
              "Food" = "#FF7F00", "Environmental" = "#377EB8", "Unknown" = "grey85")
PHY_PALETTE <- c("#E58606", "#5D69B1", "#52BCA3", "#99C945", "#CC61B0",
                 "#24796C", "#DAA51B", "#2F8AC4", "#764E9F", "#ED645A",
                 "#CC3A8E", "#A5AA99", "#9630EE", "#B79B6E", "#5B1A18")


## ============================================================================
## SECTION 1. PACKAGES AND HELPERS
## ============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(lme4)        # mixed models with phylum as a random effect
  library(lmerTest)    # P values for linear mixed models (Satterthwaite)
})

## Global figure defaults: minimum readable sizes for every plot, applied
## before any figure is built so later figures inherit them automatically.
theme_set(theme_classic(base_size = BASE_FONT))
theme_update(
  axis.title   = element_text(size = BASE_FONT),
  axis.text    = element_text(size = BASE_FONT * 0.85, colour = "black"),
  legend.title = element_text(size = BASE_FONT),
  legend.text  = element_text(size = BASE_FONT * 0.85),
  strip.text   = element_text(size = BASE_FONT * 0.9),
  plot.title   = element_text(size = BASE_FONT * 1.15, face = "bold"),
  plot.subtitle= element_text(size = BASE_FONT * 0.85)
)

setwd(path.expand(BASE))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(OUTDIR, f)

## Collapse GTDB suffixes (Bacteroidota_A -> Bacteroidota)
collapse_phylum <- function(x) {
  x <- ifelse(is.na(x) | x == "" | x == "NA", "Unclassified", x)
  x <- str_replace(x, "_[A-Z]+$", "")
  str_replace(x, "_[0-9]+$", "")
}

hhi <- function(v) { p <- table(v) / length(v); sum(p^2) }

## ============================================================================
## SECTION 2. LOAD INPUTS
## ============================================================================
cov <- fread(F_COV)
setnames(cov, 1, "MAG_id")
SAMPLES   <- setdiff(names(cov), "MAG_id")
AIR_COLS  <- grep("^A_", SAMPLES, value = TRUE)
SOIL_COLS <- grep("^S_", SAMPLES, value = TRUE)
WAT_COLS  <- grep("^W_", SAMPLES, value = TRUE)

## Metadata: CheckM2 quality and GTDB-Tk taxonomy for all 894 bins; keep the
## 532 dereplicated genomes present in the coverage table.
meta <- fread(F_META, na.strings = c("N/A", "NA", "")) %>%
  transmute(MAG_id = as.character(MAG_id),
            quality, completeness, contamination,
            genome_size = as.numeric(genome_size),
            n50_contigs = as.numeric(n50_contigs),
            gc_percent  = as.numeric(gc_percent),
            gtdb_classification,
            gtdb_closest_ani = as.numeric(gtdb_closest_ani)) %>%
  filter(MAG_id %in% cov$MAG_id)
stopifnot(nrow(meta) == nrow(cov), !any(is.na(meta$genome_size)))

ko_long <- fread(F_KO, header = FALSE, col.names = c("MAG_id", "KO")) %>%
  distinct(MAG_id, KO)

## Sequencing depth: primary reads per sample (each read counted once)
reads <- fread(F_READS) %>% mutate(pct_mapped = 100 * primary_mapped / primary)
stopifnot(setequal(reads$sample, SAMPLES))
reads <- reads[match(SAMPLES, reads$sample), ]
cat("\n=== Reads per sample and % mapped to the 532 MAGs ===\n"); print(reads)

## Depth-normalised coverage: mean coverage per 100 million primary reads
cov_norm <- as.data.frame(cov)
for (s in SAMPLES) cov_norm[[s]] <- cov_norm[[s]] / (reads$primary[reads$sample == s] / 1e8)

cat("\nLoaded:", nrow(cov), "genomes x", length(SAMPLES), "samples;",
    n_distinct(ko_long$MAG_id), "genomes with KO annotation\n")


## ============================================================================
## SECTION 3. BUILD THE MAG TABLE
## ============================================================================
## All air-association metrics use depth-normalised coverage (cov_norm).
##   lfc_own      : PRIMARY. log2[(own-site air + PSEUDO) / (mean own-site soil
##                  or water + PSEUDO)]. Each site has one air sample and two
##                  soil or water samples, so this compares a genome with itself
##                  at the place it was assembled; differences in airborne
##                  biomass between sites cancel out.
##   air_fraction : descriptive only. Mean normalised coverage over the five
##                  air samples / (that + mean over soil + mean over water).
##   majority_env : environment type holding >50% of that pooled signal (Fig. 1A).
cov_df <- as.data.frame(cov)                     # raw coverage (overview plots, QC)
site_of_genome <- sub("\\..*$", "", cov_df$MAG_id)
site_of_sample <- substr(sub("^[A-Z]_", "", SAMPLES), 1, 2)

own_air <- own_env <- numeric(nrow(cov_norm))
for (i in seq_len(nrow(cov_norm))) {
  s <- site_of_genome[i]
  a <- SAMPLES[site_of_sample == s & SAMPLES %in% AIR_COLS]
  e <- SAMPLES[site_of_sample == s & !SAMPLES %in% AIR_COLS]
  stopifnot(length(a) == 1, length(e) >= 1)
  own_air[i] <- cov_norm[i, a]
  own_env[i] <- mean(unlist(cov_norm[i, e]))
}

## Sites whose air sample barely recruits MAG reads
air_map <- reads %>% filter(sample %in% AIR_COLS) %>%
  mutate(site = substr(sub("^A_", "", sample), 1, 2))
ASSESSABLE_SITES <- air_map$site[air_map$pct_mapped >= MIN_AIR_MAPPED_PCT]
cat("\nSites with an assessable air sample (>=", MIN_AIR_MAPPED_PCT, "% reads mapped):",
    paste(ASSESSABLE_SITES, collapse = ", "), "\n")

ko_totals <- ko_long %>% count(MAG_id, name = "total_KO")

mags <- cov_norm %>%
  transmute(
    MAG_id,
    site       = site_of_genome,
    air_norm   = rowMeans(across(all_of(AIR_COLS))),
    soil_norm  = rowMeans(across(all_of(SOIL_COLS))),
    water_norm = rowMeans(across(all_of(WAT_COLS)))
  ) %>%
  mutate(
    total_cov      = rowSums(cov_df[, SAMPLES]),           # raw, for overview plots
    pooled         = air_norm + soil_norm + water_norm,
    air_fraction   = ifelse(pooled > 0, air_norm   / pooled, NA_real_),
    soil_fraction  = ifelse(pooled > 0, soil_norm  / pooled, NA_real_),
    water_fraction = ifelse(pooled > 0, water_norm / pooled, NA_real_),
    own_air_norm   = own_air,
    own_env_norm   = own_env,
    lfc_own        = log2((own_air + PSEUDO) / (own_env + PSEUDO)),
    assessable     = site %in% ASSESSABLE_SITES,
    majority_env   = case_when(air_fraction   > 0.5 ~ "Air",
                               soil_fraction  > 0.5 ~ "Soil",
                               water_fraction > 0.5 ~ "Water",
                               TRUE ~ "No majority")
  ) %>%
  left_join(meta, by = "MAG_id") %>%
  left_join(ko_totals, by = "MAG_id") %>%
  mutate(
    domain = str_remove(str_extract(gtdb_classification, "d__[^;]+"), "^d__"),
    phylum = collapse_phylum(str_remove(str_extract(gtdb_classification, "p__[^;]+"), "^p__"))
  )

bact_all <- mags %>% filter(domain == "Bacteria")      # all sites (descriptive)
bact     <- bact_all %>% filter(assessable)            # sites used for the tests

classify <- function(df, cut) {
  df %>% mutate(air_group = case_when(
    lfc_own >=  cut ~ "Air-enriched",
    lfc_own <= -cut ~ "Air-depleted",
    TRUE ~ "Intermediate"))
}
bact     <- classify(bact, LFC_CUT)
bact_all <- classify(bact_all, LFC_CUT) %>%
  mutate(air_group = ifelse(assessable, air_group, "Not assessable"))

fwrite(mags, out("MAG_table_corrected.tsv"), sep = "\t")


## ============================================================================
## SECTION 4. QC  (stop if the corrected table does not behave correctly)
## ============================================================================
peak_sample <- SAMPLES[max.col(as.matrix(cov_df[, SAMPLES]), ties.method = "first")]
peak_site   <- substr(sub("^[A-Z]_", "", peak_sample), 1, 2)
own_peak    <- mean(peak_site == site_of_genome)

qc_site <- tibble(site = site_of_genome, own = peak_site == site_of_genome) %>%
  group_by(site) %>% summarise(n = n(), own_peak = sum(own), pct = 100 * mean(own))

cat("\n=== QC 1: genomes peaking in a sample from their own site ===\n")
print(qc_site)
cat(sprintf("Overall: %.1f%% (was 32.9%% with the scrambled table)\n", 100 * own_peak))
stopifnot(own_peak > 0.95)

chk <- cov_df[cov_df$MAG_id == "MK.bin.161", c("S_MK1", "S_MK2")]
cat("\n=== QC 2: MK.bin.161 coverage in MK soils (expected ~7.26 and ~10.83) ===\n")
print(chk)
stopifnot(abs(chk$S_MK1 - 7.263895) < 1e-3, abs(chk$S_MK2 - 10.834081) < 1e-3)

site_tab <- function(x) paste(names(table(x)), table(x), sep = ":", collapse = ", ")
counts_by_cutoff <- bind_rows(lapply(LFC_CUTOFFS, function(cut) {
  b <- classify(bact, cut)
  tibble(fold = round(2^cut, 2), log2_cut = round(cut, 3),
         n_enriched = sum(b$air_group == "Air-enriched"),
         n_intermediate = sum(b$air_group == "Intermediate"),
         n_depleted = sum(b$air_group == "Air-depleted"),
         enriched_by_site = site_tab(b$site[b$air_group == "Air-enriched"]),
         depleted_by_site = site_tab(b$site[b$air_group == "Air-depleted"]))
}))
cat("\n=== QC 3: group sizes at each cutoff ===\n")
print(counts_by_cutoff)

fwrite(qc_site,          out("00_QC_own_site_peak.tsv"), sep = "\t")
fwrite(counts_by_cutoff, out("00_QC_group_sizes_by_cutoff.tsv"), sep = "\t")


## ============================================================================
## SECTION 5. FIGURE 1  Air vs paired source
## ============================================================================
## Each MAG is compared with itself: abundance in the air sample vs the mean of
## the soil or water samples collected with it (its paired source). Sites are
## only the pairing unit and are not compared with each other.
##   A  Air vs paired-source abundance for every assessable bacterial MAG
##   B  Number of MAGs in each class
## Per-site breakdown is written to a table (Supp.), not plotted.
f1_df <- bact %>%
  mutate(air_group = factor(air_group, levels = c("Air-enriched", "Intermediate", "Air-depleted")),
         x = own_env_norm + PSEUDO, y = own_air_norm + PSEUDO)
grp_cols <- c("Air-enriched" = COL_ENR, "Intermediate" = COL_INT, "Air-depleted" = COL_DEP)

p1a <- ggplot(f1_df, aes(x = x, y = y, colour = air_group)) +
  geom_abline(slope = 1, intercept = 0, colour = "black", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = c(-1, 1) * LFC_CUT * log10(2),
              linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_point(size = 1.4, alpha = 0.8) +
  scale_x_log10(labels = label_log()) + scale_y_log10(labels = label_log()) +
  scale_colour_manual(values = grp_cols, name = NULL) +
  coord_equal() +
  labs(x = "Coverage in paired soil or water (per 100 M reads)",
       y = "Coverage in air (per 100 M reads)", title = "A") +
  theme_classic(base_size = BASE_FONT) + theme(legend.position = "bottom")

f1b_df <- f1_df %>% count(air_group)
p1b <- ggplot(f1b_df, aes(x = air_group, y = n, fill = air_group)) +
  geom_col(colour = "black", linewidth = 0.25, width = 0.65) +
  geom_text(aes(label = n), vjust = -0.4, size = BASE_FONT * 0.28) +
  scale_fill_manual(values = grp_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = NULL, y = "Number of bacterial MAGs", title = "B") +
  theme_classic(base_size = BASE_FONT)

pdf(out("Fig1_air_vs_source.pdf"), width = 12, height = 6)
print(p1a + p1b + plot_layout(widths = c(1.4, 1)))
dev.off()

## Supplementary table: recovery and classes per sampling pair
supp_pairs <- bact_all %>%
  group_by(site) %>%
  summarise(MAGs = n(), HQ = sum(quality == "HQ"), MQ = sum(quality == "MQ"),
            assessable = first(assessable),
            air_enriched = sum(air_group == "Air-enriched"),
            intermediate = sum(air_group == "Intermediate"),
            air_depleted = sum(air_group == "Air-depleted"), .groups = "drop") %>%
  left_join(air_map %>% select(site, air_reads = primary, air_pct_mapped = pct_mapped), by = "site")
fwrite(supp_pairs, out("SuppTable_pairs_recovery_and_classes.tsv"), sep = "\t")
print(supp_pairs)


## ============================================================================
## SECTION 6. FIGURE 2  Phylum composition and dominance (HHI)
## ============================================================================
## Two tests are reported:
##   (1) bootstrap of each group separately (manuscript method), and
##   (2) a size-matched null: random subsamples of the air-depleted group of the
##       SAME size as the air-enriched group. HHI rises automatically in small
##       samples, so (2) is the fair comparison when group sizes differ.
two <- bact %>% filter(air_group != "Intermediate")
ph_enr <- two$phylum[two$air_group == "Air-enriched"]
ph_dep <- two$phylum[two$air_group == "Air-depleted"]
n_enr  <- length(ph_enr)

hhi_enr <- hhi(ph_enr); hhi_dep <- hhi(ph_dep)
boot_enr <- replicate(B_BOOT, hhi(sample(ph_enr, replace = TRUE)))
boot_dep <- replicate(B_BOOT, hhi(sample(ph_dep, replace = TRUE)))
delta_boot <- boot_enr - boot_dep
p_boot <- 2 * min(mean(delta_boot <= 0), mean(delta_boot >= 0))

null_size <- replicate(B_BOOT, hhi(sample(ph_dep, n_enr, replace = FALSE)))
p_size_matched <- (sum(null_size >= hhi_enr) + 1) / (B_BOOT + 1)

hhi_tbl <- tibble(
  n_enriched = n_enr, n_depleted = length(ph_dep),
  HHI_enriched = hhi_enr, HHI_depleted = hhi_dep, delta_HHI = hhi_enr - hhi_dep,
  boot_CI_low = quantile(delta_boot, 0.025), boot_CI_high = quantile(delta_boot, 0.975),
  p_bootstrap = p_boot,
  size_matched_null_mean = mean(null_size),
  size_matched_null_95 = quantile(null_size, 0.95),
  p_size_matched = p_size_matched
)
cat("\n=== HHI ===\n"); print(t(hhi_tbl))
fwrite(hhi_tbl, out("Fig2_HHI_stats.tsv"), sep = "\t")

phylum_tbl <- two %>%
  count(air_group, phylum) %>%
  group_by(air_group) %>% mutate(n_total = sum(n), pct = 100 * n / n_total) %>%
  ungroup() %>% arrange(air_group, desc(pct))
fwrite(phylum_tbl, out("Fig2_phylum_composition.tsv"), sep = "\t")

phy_order <- names(PHY_COLS)
p2_df <- phylum_tbl %>%
  mutate(phylum = ifelse(phylum %in% phy_order, phylum, "Other")) %>%
  group_by(air_group, n_total, phylum) %>% summarise(pct = sum(pct), .groups = "drop") %>%
  mutate(air_group = factor(air_group, levels = c("Air-enriched", "Air-depleted")))

p2a <- ggplot(p2_df %>% mutate(phylum = factor(phylum, levels = phy_order),
                               air_group = factor(paste0(air_group, "\n(n = ", n_total, ")"),
                                                  levels = paste0(levels(air_group), "\n(n = ",
                                                                  n_total[match(levels(air_group), air_group)], ")"))),
              aes(x = air_group, y = pct / 100, fill = phylum)) +
  geom_col(color = "black", linewidth = 0.25, width = 0.7) +
  scale_fill_manual(values = PHY_COLS, name = "Phylum", drop = TRUE) +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Fraction of MAGs", title = "A") +
  theme_classic(base_size = BASE_FONT)

null_df <- tibble(HHI = null_size)
p2b <- ggplot(null_df, aes(x = HHI)) +
  geom_histogram(bins = 40, fill = COL_DEP, colour = "white") +
  geom_vline(xintercept = hhi_enr, colour = COL_ENR, linewidth = 1) +
  annotate("text", x = hhi_enr, y = Inf, vjust = 1.5, hjust = 1.05, colour = COL_ENR, size = BASE_FONT * 0.25,
           label = sprintf("Air-enriched\nHHI = %.2f", hhi_enr)) +
  annotate("text", x = hhi_enr, y = Inf, vjust = 5.2, hjust = 1.05, colour = COL_ENR, size = BASE_FONT * 0.25,
           label = sprintf("italic(p) == %.3f", p_size_matched), parse = TRUE) +
  labs(x = sprintf("HHI of %d randomly drawn air-depleted MAGs", n_enr),
       y = "Number of random sets", title = "B") +
  theme_classic(base_size = BASE_FONT)

pdf(out("Fig2_phylum_composition.pdf"), width = 12, height = 6)
print(p2a + p2b + plot_layout(widths = c(1, 1.1)))
dev.off()

## Within-site check: most air-enriched MAGs come from DF, so repeat the
## size-matched test using only DF genomes (removes between-site differences
## in phylum composition from the null).
df_two <- two %>% filter(site == "DF")
e_df <- df_two$phylum[df_two$air_group == "Air-enriched"]
d_df <- df_two$phylum[df_two$air_group == "Air-depleted"]
if (length(e_df) >= 3 && length(d_df) > length(e_df)) {
  null_df_site <- replicate(B_BOOT, hhi(sample(d_df, length(e_df), replace = FALSE)))
  hhi_df_tbl <- tibble(site = "DF", n_enriched = length(e_df), n_depleted = length(d_df),
                       HHI_enriched = hhi(e_df), HHI_depleted = hhi(d_df),
                       size_matched_null_mean = mean(null_df_site),
                       p_size_matched = (sum(null_df_site >= hhi(e_df)) + 1) / (B_BOOT + 1))
  cat("\n=== HHI within DF only ===\n"); print(t(hhi_df_tbl))
  fwrite(hhi_df_tbl, out("Fig2_HHI_withinDF.tsv"), sep = "\t")
}


## ============================================================================
## SECTION 7. CURATED MARKER GENE MAP  (47 KOs, 10 categories; Supplementary Table S2)
## ============================================================================
ko_trait_map <- tribble(
  ~Trait_Category, ~Gene, ~KO,
  "DNA_Repair_Genome_Stability", "recA", "K03553",
  "DNA_Repair_Genome_Stability", "uvrA", "K03701",
  "DNA_Repair_Genome_Stability", "uvrB", "K03702",
  "DNA_Repair_Genome_Stability", "uvrC", "K03703",
  "DNA_Repair_Genome_Stability", "uvrD", "K03657",
  "DNA_Repair_Genome_Stability", "mutS", "K03555",
  "DNA_Repair_Genome_Stability", "mutL", "K03572",
  "DNA_Repair_Genome_Stability", "lexA", "K01356",
  "Oxidative_Stress_Defense", "sodA/sodB", "K04564",
  "Oxidative_Stress_Defense", "katG", "K03782",
  "Oxidative_Stress_Defense", "ahpC", "K03386",
  "Oxidative_Stress_Defense", "trxB", "K00384",
  "Oxidative_Stress_Defense", "gor",  "K00383",
  "Resource_Storage_Compounds", "ppk1", "K00937",
  "Resource_Storage_Compounds", "ppk2", "K22468",
  "Resource_Storage_Compounds", "ppx",  "K01524",
  "Resource_Storage_Compounds", "phaA", "K00626",
  "Resource_Storage_Compounds", "phaB", "K00023",
  "Resource_Storage_Compounds", "phaC", "K03821",
  "Resource_Storage_Compounds", "glgA", "K00703",
  "Resource_Storage_Compounds", "glgC", "K00975",
  "Spore_Dormancy", "spo0A",   "K07699",
  "Spore_Dormancy", "spoIIAA", "K06378",
  "Spore_Dormancy", "spoIIAB", "K06379",
  "Spore_Dormancy", "spoIIE",  "K06382",
  "Pigmentation_UV_Protection", "crtB", "K02291",
  "Pigmentation_UV_Protection", "phr",  "K01669",
  "Pigmentation_UV_Protection", "crtI", "K10027",
  "EPS_Biofilm_Matrix", "bcsA", "K00694",
  "EPS_Biofilm_Matrix", "pgaA", "K11935",
  "EPS_Biofilm_Matrix", "pgaC", "K11936",
  "Desiccation_Osmoprotection", "otsA", "K00697",
  "Desiccation_Osmoprotection", "otsB", "K01087",
  "Desiccation_Osmoprotection", "betA", "K00108",
  "Desiccation_Osmoprotection", "ectA", "K06718",
  "Desiccation_Osmoprotection", "ectC", "K06720",
  "Acid_pH_Stress_Tolerance", "gadA/gadB", "K01580",
  "Acid_pH_Stress_Tolerance", "gadC",      "K20265",
  "Acid_pH_Stress_Tolerance", "ureC",      "K01428",
  "Membrane_Homeostasis_Mechanosensing", "mscL", "K03282",
  "Membrane_Homeostasis_Mechanosensing", "mscS", "K03442",
  "Membrane_Homeostasis_Mechanosensing", "cfa",  "K00574",
  "C1_and_Trace_Gas_Metabolism", "coxL", "K03520",
  "C1_and_Trace_Gas_Metabolism", "coxM", "K03519",
  "C1_and_Trace_Gas_Metabolism", "coxS", "K03518",
  "C1_and_Trace_Gas_Metabolism", "xoxF", "K23995",
  "C1_and_Trace_Gas_Metabolism", "mxaF", "K14028"
)
stopifnot(nrow(ko_trait_map) == 47, !any(duplicated(ko_trait_map$KO)))

## KO numbers verified against the KEGG KO list (rest.kegg.jp/list/ko,
## downloaded 22 Oct 2026). Corrections made on that date:
##   uvrD K03555 -> K03657 (K03555 is mutS); mutS K02335 -> K03555 (K02335 is polA);
##   katG K03781 -> K03782 (K03781 is katE); trxB K03671 -> K00384 (K03671 is trxA);
##   sodB K04565 removed (K04565 is Cu-Zn SOD; KEGG has one Fe-Mn SOD KO, K04564,
##   now labelled sodA/sodB); gadC K01581 -> K20265 (K01581 is ornithine
##   decarboxylase); mscL K03823 -> K03282; mscS K03824 -> K03442;
##   coxL and coxS were swapped (coxL = K03520, coxS = K03518);
##   hhyL K22216 and hhyS K22217 were not hydrogenase KOs (ilvD1, lactonase).
##   KEGG has no KO specific to the high-affinity group 1h [NiFe]-hydrogenase,
##   so hhyL/hhyS are dropped rather than replaced, keeping the list a priori
##   (47 genes).
F_KEGG <- "data/kegg_ko_list.tsv"
if (file.exists(F_KEGG)) {
  kegg <- fread(F_KEGG, header = FALSE, sep = "\t", col.names = c("KO", "kegg_definition"))
  ## KEGG uses a different symbol for three markers (same enzyme)
  kegg_synonym <- c(phaA = "atoB", phaB = "phbB", "sodA/sodB" = "SOD2")
  ko_check <- ko_trait_map %>% left_join(kegg, by = "KO") %>%
    mutate(symbol_found = mapply(function(g, d) any(sapply(c(strsplit(g, "/")[[1]],
                                  if (g %in% names(kegg_synonym)) kegg_synonym[[g]]),
                                  function(x) grepl(paste0("\\b", x, "\\b"), d, ignore.case = TRUE))),
                                 Gene, kegg_definition))
  fwrite(ko_check, out("SuppTable_marker_KO_verification.tsv"), sep = "\t")
  if (!all(ko_check$symbol_found)) {
    print(ko_check %>% filter(!symbol_found))
    stop("A marker KO does not match its gene symbol in KEGG; check the table above.")
  }
  cat("All", nrow(ko_check), "marker KOs match their gene symbols in KEGG.\n")
}

trait_labels <- c(
  DNA_Repair_Genome_Stability         = "DNA repair",
  Oxidative_Stress_Defense            = "Oxidative stress defence",
  Resource_Storage_Compounds          = "Storage compounds",
  Spore_Dormancy                      = "Sporulation",
  Pigmentation_UV_Protection          = "Pigmentation and UV protection",
  EPS_Biofilm_Matrix                  = "Extracellular polymeric substances",
  Desiccation_Osmoprotection          = "Desiccation and osmoprotection",
  Acid_pH_Stress_Tolerance            = "Acid tolerance",
  Membrane_Homeostasis_Mechanosensing = "Membrane homeostasis and mechanosensing",
  C1_and_Trace_Gas_Metabolism         = "C1 and trace gas metabolism"
)

## Presence matrix: one row per bacterial MAG, one column per marker gene (0/1)
pres_long <- ko_long %>%
  inner_join(ko_trait_map, by = "KO") %>%
  distinct(MAG_id, Gene) %>% mutate(present = 1L)
pres <- expand_grid(MAG_id = bact_all$MAG_id, Gene = ko_trait_map$Gene) %>%
  left_join(pres_long, by = c("MAG_id", "Gene")) %>%
  mutate(present = coalesce(present, 0L)) %>%
  left_join(ko_trait_map %>% select(Gene, Trait_Category), by = "Gene")

detected_genes <- pres %>% group_by(Gene) %>% summarise(n = sum(present)) %>%
  filter(n > 0) %>% pull(Gene)
cat("\nMarker genes detected in at least one bacterial MAG:", length(detected_genes), "of", nrow(ko_trait_map), "\n")
cat("Not detected:", paste(setdiff(ko_trait_map$Gene, detected_genes), collapse = ", "), "\n")


## ============================================================================
## SECTION 11B. RQ2  FILTERING HYPOTHESIS, PHYLOGENY-AWARE TESTS
## ============================================================================
## H (a priori): genomes enriched in air carry more atmospheric stress-tolerance
## traits than expected from their phylogeny and genome size.
##   1. Genome characteristics by group (annotation depth must be reported)
##   2. Mixed models: phylum as a random effect; annotation depth (log total
##      KOs), completeness and site as fixed effects. Two predictors:
##        group : air-enriched vs air-depleted (odds ratio, OR)
##        lfc   : continuous log2 air / paired source, all assessable MAGs
##   3. Within-phylum replication (Bacillota, Actinomycetota)
##   4. Smallest detectable effect, so a null result can be interpreted
MIN_PHY_N  <- 10                         # phyla with >= this many MAGs per group
WITHIN_PHY <- c("Bacillota", "Actinomycetota")

## ---- 1. Genome characteristics -------------------------------------------
two_g <- bact %>% filter(air_group != "Intermediate") %>%
  mutate(grp = factor(air_group, levels = c("Air-depleted", "Air-enriched")))
gchar <- two_g %>%
  select(grp, completeness, contamination, genome_size, gc_percent, total_KO) %>%
  pivot_longer(-grp, names_to = "variable") %>%
  group_by(variable, grp) %>%
  summarise(n = n(), median = median(value), q1 = quantile(value, 0.25),
            q3 = quantile(value, 0.75), .groups = "drop")
gtest <- bind_rows(lapply(c("completeness", "contamination", "genome_size", "gc_percent", "total_KO"),
  function(v) {
    y <- two_g[[v]]; if (v %in% c("genome_size", "total_KO")) y <- log(y)
    fit <- suppressMessages(lmerTest::lmer(y ~ grp + site + (1 | phylum), data = two_g %>% mutate(y = y)))
    co  <- summary(fit)$coefficients["grpAir-enriched", ]
    fit0 <- lm(y ~ grp + site, data = two_g %>% mutate(y = y))
    tibble(variable = v, scale = ifelse(v %in% c("genome_size", "total_KO"), "log", "raw"),
           effect_site_adj = coef(fit0)["grpAir-enriched"],
           p_site_adj = summary(fit0)$coefficients["grpAir-enriched", 4],
           effect_site_phylum_adj = co["Estimate"],
           t_site_phylum_adj = co["t value"],
           p_site_phylum_adj = co["Pr(>|t|)"])  # Satterthwaite df (lmerTest)
  }))
cat("\n=== Genome characteristics by group ===\n"); print(gchar, n = Inf); print(gtest)
fwrite(gchar, out("RQ2_genome_characteristics.tsv"), sep = "\t")
fwrite(gtest, out("RQ2_genome_characteristics_tests.tsv"), sep = "\t")

## ---- 2. Mixed models -------------------------------------------------------
mm_dat <- bact %>%
  mutate(grp = factor(ifelse(air_group == "Air-enriched", "Air-enriched",
                             ifelse(air_group == "Air-depleted", "Air-depleted", NA)),
                      levels = c("Air-depleted", "Air-enriched")),
         z_logKO = as.numeric(scale(log(total_KO))),
         z_compl = as.numeric(scale(completeness)),
         phylum  = factor(phylum))

fit_one <- function(df, formula, term, family = binomial) {
  ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  w <- NULL
  fit <- withCallingHandlers(
    tryCatch({
      f0 <- suppressMessages(glmer(formula, data = df, family = family, control = ctrl))
      list(fit = f0, co = summary(f0)$coefficients, singular = isSingular(f0))
    }, error = function(e) NULL),
    warning = function(x) { w <<- conditionMessage(x); invokeRestart("muffleWarning") })
  if (is.null(fit)) return(tibble(est = NA_real_, se = NA_real_, p = NA_real_, warn = "failed"))
  if (fit$singular && is.null(w)) w <- "phylum variance ~0 (singular)"
  co <- fit$co
  if (!term %in% rownames(co)) return(tibble(est = NA_real_, se = NA_real_, p = NA_real_, warn = "no term"))
  tibble(est = co[term, "Estimate"], se = co[term, "Std. Error"], p = co[term, "Pr(>|z|)"],
         warn = ifelse(is.null(w), "", w))
}

## Category level: genes present out of genes in the category (binomial)
cat_counts <- pres %>% filter(Gene %in% detected_genes) %>%
  group_by(MAG_id, Trait_Category) %>%
  summarise(k = sum(present), m = n(), .groups = "drop")

mm_category <- function(df, pred, counts = cat_counts) {
  term <- if (pred == "group") "grpAir-enriched" else "lfc_own"
  rhs  <- if (pred == "group") "grp" else "lfc_own"
  d0 <- counts %>% inner_join(df, by = "MAG_id")
  if (pred == "group") d0 <- d0 %>% filter(!is.na(grp))
  d0 %>% group_by(Trait_Category) %>%
    group_modify(function(x, k) {
      f <- as.formula(paste("cbind(k, m - k) ~", rhs, "+ z_logKO + z_compl + site + (1 | phylum)"))
      fit_one(x, f, term)
    }) %>% ungroup() %>%
    mutate(predictor = pred, padj = p.adjust(p, "BH"))
}

mm_gene <- function(df, pred, ptab = pres, genes = detected_genes) {
  term <- if (pred == "group") "grpAir-enriched" else "lfc_own"
  rhs  <- if (pred == "group") "grp" else "lfc_own"
  d0 <- ptab %>% filter(Gene %in% genes) %>% inner_join(df, by = "MAG_id")
  if (pred == "group") d0 <- d0 %>% filter(!is.na(grp))
  d0 %>% group_by(Trait_Category, Gene) %>%
    group_modify(function(x, k) {
      if (sum(x$present) < 5 || sum(1 - x$present) < 5) return(tibble())
      f <- as.formula(paste("present ~", rhs, "+ z_logKO + z_compl + site + (1 | phylum)"))
      fit_one(x, f, term) %>%
        mutate(prev_enr = if (pred == "group") mean(x$present[x$grp == "Air-enriched"]) else NA_real_,
               prev_dep = if (pred == "group") mean(x$present[x$grp == "Air-depleted"]) else NA_real_)
    }) %>% ungroup() %>%
    mutate(predictor = pred, padj = p.adjust(p, "BH"))
}

add_or <- function(x) x %>% mutate(OR = exp(est), OR_lo = exp(est - 1.96 * se), OR_hi = exp(est + 1.96 * se))
mm_cat_grp  <- add_or(mm_category(mm_dat, "group"))
mm_cat_lfc  <- add_or(mm_category(mm_dat, "lfc"))
mm_gene_grp <- add_or(mm_gene(mm_dat, "group"))
mm_gene_lfc <- add_or(mm_gene(mm_dat, "lfc"))

cat("\n=== RQ2 mixed model, category level (air-enriched vs depleted) ===\n")
print(mm_cat_grp %>% select(Trait_Category, OR, OR_lo, OR_hi, p, padj, warn) %>% arrange(p), n = Inf)
cat("\n=== RQ2 mixed model, category level (continuous log2 ratio) ===\n")
print(mm_cat_lfc %>% select(Trait_Category, OR, OR_lo, OR_hi, p, padj, warn) %>% arrange(p), n = Inf)
cat("\n=== RQ2 mixed model, gene level, BH < 0.1 ===\n")
print(bind_rows(mm_gene_grp, mm_gene_lfc) %>% filter(padj < 0.1) %>%
        select(predictor, Gene, OR, OR_lo, OR_hi, padj))
fwrite(bind_rows(mm_cat_grp, mm_cat_lfc),   out("RQ2_mixed_category_stats.tsv"), sep = "\t")
fwrite(bind_rows(mm_gene_grp, mm_gene_lfc), out("RQ2_mixed_gene_stats.tsv"), sep = "\t")

## ---- 3. Within-phylum replication (no random effect needed) ----------------
within_phy <- bind_rows(lapply(WITHIN_PHY, function(ph) {
  d <- mm_dat %>% filter(phylum == ph, !is.na(grp))
  ne <- sum(d$grp == "Air-enriched"); nd <- sum(d$grp == "Air-depleted")
  cat_counts %>% inner_join(d, by = "MAG_id") %>%
    group_by(Trait_Category) %>%
    group_modify(function(x, k) {
      fit <- tryCatch(suppressWarnings(glm(cbind(k, m - k) ~ grp + z_logKO + z_compl,
                                           family = binomial, data = x)), error = function(e) NULL)
      if (is.null(fit)) return(tibble(est = NA_real_, se = NA_real_, p = NA_real_))
      co <- summary(fit)$coefficients["grpAir-enriched", ]
      tibble(est = co["Estimate"], se = co["Std. Error"], p = co["Pr(>|z|)"])
    }) %>% ungroup() %>%
    mutate(phylum = ph, n_enr = ne, n_dep = nd, padj = p.adjust(p, "BH"))
})) %>% add_or()
cat("\n=== RQ2 within-phylum, category level ===\n")
print(within_phy %>% select(phylum, Trait_Category, n_enr, n_dep, OR, OR_lo, OR_hi, padj), n = Inf)
fwrite(within_phy, out("RQ2_within_phylum_category_stats.tsv"), sep = "\t")

## ---- 4. Smallest detectable effect (gene level) ----------------------------
## Two-proportion arcsine approximation, two-sided alpha 0.05, power 0.80,
## unadjusted. Gives the smallest prevalence in air-enriched MAGs (and OR)
## that this design would detect, given the prevalence in air-depleted MAGs.
n1 <- sum(mm_dat$grp == "Air-enriched", na.rm = TRUE)
n2 <- sum(mm_dat$grp == "Air-depleted", na.rm = TRUE)
h  <- (qnorm(0.975) + qnorm(0.80)) * sqrt(1 / n1 + 1 / n2)
mde <- mm_gene_grp %>%
  transmute(Trait_Category, Gene, prev_dep, prev_enr,
            detectable_up   = sin(pmin(pi / 2, asin(sqrt(prev_dep)) + h / 2))^2,
            detectable_down = sin(pmax(0,      asin(sqrt(prev_dep)) - h / 2))^2,
            OR_detectable_up   = (detectable_up / (1 - detectable_up)) / (prev_dep / (1 - prev_dep)),
            OR_detectable_down = (detectable_down / (1 - detectable_down)) / (prev_dep / (1 - prev_dep)))
cat(sprintf("\n=== Smallest detectable difference, n = %d vs %d (80%% power) ===\n", n1, n2))
print(summary(mde$detectable_up - mde$prev_dep))
fwrite(mde, out("RQ2_detectable_effects.tsv"), sep = "\t")

## ---- Fig. 3: category-level odds ratios ------------------------------------
f3_or <- bind_rows(
  mm_cat_grp %>% mutate(model = sprintf("All phyla, phylum random effect (%d vs %d)", n1, n2)),
  within_phy %>% mutate(model = sprintf("Within %s (%d vs %d)", phylum, n_enr, n_dep))
) %>%
  filter(is.finite(OR_lo), is.finite(OR_hi), OR_hi / OR_lo < 1e3) %>%
  mutate(Trait = recode(Trait_Category, !!!trait_labels),
         sig = !is.na(padj) & padj < 0.05)
ord <- mm_cat_grp %>% arrange(est) %>% pull(Trait_Category) %>% recode(!!!trait_labels)
f3_or <- f3_or %>% mutate(Trait = factor(Trait, levels = ord))
## One panel per model, sharing the category axis. Model and group sizes are
## in the panel titles, so no colour legend is needed. Filled points: BH < 0.05.
## A category with no estimate in a panel (complete separation) is marked "n.e.".
f3_or <- f3_or %>%
  mutate(panel = ifelse(is.na(phylum), sprintf("All phyla\n(%d vs %d)", n1, n2),
                        sprintf("%s\n(%d vs %d)", phylum, n_enr, n_dep)))
panel_lv <- c(sprintf("All phyla\n(%d vs %d)", n1, n2),
              setdiff(unique(f3_or$panel), sprintf("All phyla\n(%d vs %d)", n1, n2)))
f3_or <- f3_or %>% mutate(panel = factor(panel, levels = panel_lv))
f3_ne <- expand_grid(panel = factor(panel_lv, levels = panel_lv), Trait = factor(ord, levels = ord)) %>%
  anti_join(f3_or %>% distinct(panel, Trait), by = c("panel", "Trait"))
X_LIM <- c(0.02, 50)
p_f3or <- ggplot(f3_or, aes(x = OR, y = Trait)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = pmax(OR_lo, X_LIM[1]), xmax = pmin(OR_hi, X_LIM[2])), height = 0.25) +
  geom_point(aes(fill = sig), shape = 21, size = 3, colour = "black") +
  geom_text(data = f3_ne, aes(x = 1, y = Trait, label = "n.e."), size = BASE_FONT * 0.25,
            colour = "grey40", inherit.aes = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "black", `FALSE` = "white"), guide = "none") +
  scale_x_log10(limits = X_LIM, breaks = c(0.1, 1, 10), labels = c("0.1", "1", "10")) +
  facet_grid(. ~ panel) +
  labs(x = "Odds ratio of marker gene presence (log scale)",
       y = NULL) +
  theme_classic(base_size = BASE_FONT) +
  theme(axis.text = element_text(size = BASE_FONT, colour = "black"),
        axis.title.x = element_text(size = BASE_FONT),
        strip.text = element_text(size = BASE_FONT, face = "bold"),
        strip.background = element_blank(),
        panel.spacing = unit(1.2, "lines"),
        panel.grid.major.y = element_line(colour = "grey92"))
pdf(out("Fig3_trait_odds_ratios.pdf"), width = 13, height = 6.5)
print(p_f3or)
dev.off()

## ---- 5. Sensitivity to the enrichment cut-off (Supplementary Figure S3) ----
## The category mixed model is refitted with air-enriched and air-depleted
## MAGs defined at each cut-off in LFC_CUTOFFS (1.5- to 9-fold).
sens_mm <- bind_rows(lapply(LFC_CUTOFFS, function(cut) {
  d <- mm_dat %>%
    mutate(grp = factor(case_when(lfc_own >=  cut ~ "Air-enriched",
                                  lfc_own <= -cut ~ "Air-depleted",
                                  TRUE ~ NA_character_),
                        levels = c("Air-depleted", "Air-enriched")))
  add_or(mm_category(d, "group")) %>%
    mutate(fold = 2^cut,
           n_enr = sum(d$grp == "Air-enriched", na.rm = TRUE),
           n_dep = sum(d$grp == "Air-depleted", na.rm = TRUE))
}))
cat("\n=== Cut-off sensitivity, category mixed model ===\n")
print(sens_mm %>% select(fold, n_enr, n_dep, Trait_Category, OR, padj) %>%
        arrange(Trait_Category, fold), n = Inf)
fwrite(sens_mm, out("FigS3_cutoff_sensitivity_stats.tsv"), sep = "\t")

fs2 <- sens_mm %>%
  filter(is.finite(OR_lo), is.finite(OR_hi), OR_hi / OR_lo < 1e3) %>%
  mutate(Trait = factor(recode(Trait_Category, !!!trait_labels), levels = rev(ord)),
         cut_lab = factor(sprintf("%g", round(fold, 1)), levels = sprintf("%g", round(2^LFC_CUTOFFS, 1))),
         sig = !is.na(padj) & padj < 0.05)
n_lab <- sens_mm %>% distinct(fold, n_enr, n_dep) %>%
  mutate(lab = sprintf("%g-fold: %d vs %d", round(fold, 1), n_enr, n_dep)) %>% pull(lab)
p_fs2 <- ggplot(fs2, aes(x = cut_lab, y = OR, shape = sig)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = OR_lo, ymax = OR_hi), width = 0.25, linewidth = 0.8) +
  geom_point(size = 3.8, stroke = 1.2, fill = "white") +
  scale_y_log10() +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c(`TRUE` = expression("BH" ~ italic(p) < 0.05), `FALSE` = "n.s."), name = NULL) +
  facet_wrap(~ Trait, ncol = 5, labeller = label_wrap_gen(width = 22)) +
  labs(x = "Fold difference used to define air-enriched and air-depleted MAGs",
       y = "Odds ratio of marker gene presence (95% CI)") +
  theme_classic(base_size = BASE_FONT) +
  theme(axis.text = element_text(size = BASE_FONT, colour = "black"),
        axis.title = element_text(size = BASE_FONT * 1.1),
        strip.text = element_text(size = BASE_FONT * 1.0, face = "bold"),
        strip.background = element_blank(),
        legend.position = "top", legend.text = element_text(size = BASE_FONT * 1.2),
        legend.key.size = unit(1.6, "lines"),
        panel.spacing = unit(1, "lines")) +
  guides(shape = guide_legend(override.aes = list(size = 5)))
pdf(out("FigS3_cutoff_sensitivity.pdf"), width = 15, height = 8)
print(p_fs2)
dev.off()


## ============================================================================
## SECTION 11D. RQ3  SOURCE HYPOTHESIS: HABITAT OF THE CLOSEST REFERENCE GENOME
## ============================================================================
## For MAGs with a species-level GTDB match (closest-reference ANI >= 95%),
## the habitat of that reference is taken from its NCBI BioSample record
## (datasets summary genome accession ... --as-json-lines; NCBI datasets
## v18.39.0, run 22 Oct 2026). Rules are keyword-based, applied in this order,
## and written out in full so they can be inspected:
##   1 Animal-associated : host is an animal or human, or source mentions an
##                         animal body site (gut, faeces, rumen, oral, skin ...)
##   2 Engineered/waste  : digester, sludge, wastewater, landfill, bioreactor,
##                         compost, bedding, food waste
##   3 Food              : cheese, dairy product, food
##   4 Environmental     : soil, sediment, water, marine, lake, spring, plant, alga
##   5 Unknown           : no informative field
F_REFSUM <- "data/closest_ref_summary.jsonl"
if (file.exists(F_REFSUM)) {
  suppressPackageStartupMessages(library(jsonlite))
  null_na <- function(x) if (is.null(x) || length(x) == 0) NA_character_ else as.character(x[1])
  refs <- bind_rows(lapply(readLines(F_REFSUM), function(l) {
    d  <- fromJSON(l, simplifyVector = FALSE)
    at <- d$assembly_info$biosample$attributes
    av <- setNames(vapply(at, function(z) null_na(z$value), ""),
                   tolower(vapply(at, function(z) null_na(z$name), "")))
    g  <- function(k) if (k %in% names(av)) av[[k]] else NA_character_
    tibble(gtdb_closest = d$accession, ref_organism = null_na(d$organism$organism_name),
           isolation_source = g("isolation_source"), host = g("host"),
           env_broad = g("env_broad_scale"), env_local = g("env_local_scale"),
           env_medium = g("env_medium"), metagenome_source = g("metagenome_source"))
  }))
  uninformative <- "^(missing|not applicable|not collected|na|n/a|none|unknown|mixed|metagenome|from type material|env[0-9]+|envo[:_][0-9]+)$"
  clean <- function(x) ifelse(is.na(x) | grepl(uninformative, trimws(x), ignore.case = TRUE), NA, x)
  refs <- refs %>% mutate(across(isolation_source:metagenome_source, clean))

  ## Keywords match at the start of a word (\\b), and short words that occur inside
  ## other words (colon/colony, oral/littoral, ear/near) must match as whole words.
  kw_animal  <- "\\b(gut|f(a)?eces|faecal|fecal|stool|rumen|ruminant|gastrointestinal|intestin|caec|colon\\b|oral\\b|tooth|teeth|dental|saliva|skin|ears?\\b|nasal|urine|vagin|milk|udder|joint|blood|wound|clinical|insect|grub|crustacean|mucus)"
  kw_plant_host <- "paris polyphylla|pseudopediastrum|plant|alga"
  kw_engineered <- "digest|sludge|wastewater|landfill|leachate|bioreactor|reactor|compost|bedding|fermentation bed|food waste|manure"
  kw_food    <- "cheese|dairy product|food|fermented"
  kw_env     <- "soil|sediment|water|marine|\\bsea|lake|river|spring|subsurface|mangrove|root|rhizosphere|oligotrophic|plant|alga|drainage"

  classify_habitat <- function(r) {
    txt <- tolower(paste(na.omit(c(r$isolation_source, r$env_broad, r$env_local,
                                   r$env_medium, r$metagenome_source)), collapse = " | "))
    h <- tolower(ifelse(is.na(r$host), "", r$host))
    if (nzchar(h) && !grepl(kw_plant_host, h)) return("Animal-associated")
    if (grepl(kw_animal, txt, perl = TRUE) && !grepl("manure", txt)) return("Animal-associated")
    if (grepl(kw_engineered, txt)) return("Engineered/waste")
    if (grepl(kw_food, txt)) return("Food")
    if (grepl(kw_env, txt, perl = TRUE) || grepl(kw_plant_host, h)) return("Environmental")
    "Unknown"
  }
  refs$habitat <- vapply(seq_len(nrow(refs)), function(i) classify_habitat(refs[i, ]), "")

  meta_acc <- fread(F_META) %>% transmute(MAG_id, gtdb_closest)
  hab <- bact %>% filter(air_group != "Intermediate", gtdb_closest_ani >= 95) %>%
    left_join(meta_acc, by = "MAG_id") %>% left_join(refs, by = "gtdb_closest") %>%
    mutate(habitat = coalesce(habitat, "Unknown"))
  fwrite(hab %>% select(MAG_id, site, air_group, lfc_own, phylum, gtdb_classification,
                        gtdb_closest, gtdb_closest_ani, ref_organism, isolation_source, host,
                        env_broad, env_local, env_medium, metagenome_source, habitat),
         out("RQ3_closest_reference_habitat.tsv"), sep = "\t")

  hab_tab <- hab %>% count(air_group, habitat) %>% group_by(air_group) %>%
    mutate(n_group = sum(n), pct = 100 * n / n_group) %>% ungroup()
  cat("\n=== RQ3 habitat of closest reference (species-level matches) ===\n")
  print(hab_tab, n = Inf)
  fwrite(hab_tab, out("RQ3_habitat_by_group.tsv"), sep = "\t")

  ## Test: animal-associated vs other informative habitats (Unknown excluded)
  hab_inf <- hab %>% filter(habitat != "Unknown") %>%
    mutate(animal = habitat == "Animal-associated")
  rq3_test <- bind_rows(
    tibble(subset = "BB, DF, VL", data = list(hab_inf)),
    tibble(subset = "DF only",    data = list(hab_inf %>% filter(site == "DF"))),
    tibble(subset = "BB, DF, VL; Unknown counted as not animal",
           data = list(hab %>% mutate(animal = habitat == "Animal-associated")))) %>%
    rowwise() %>%
    mutate(n_enr = sum(data$air_group == "Air-enriched"),
           n_dep = sum(data$air_group == "Air-depleted"),
           pct_animal_enr = 100 * mean(data$animal[data$air_group == "Air-enriched"]),
           pct_animal_dep = 100 * mean(data$animal[data$air_group == "Air-depleted"]),
           fisher_OR = unname(fisher.test(table(data$air_group, data$animal))$estimate),
           fisher_p  = fisher.test(table(data$air_group, data$animal))$p.value) %>%
    ungroup() %>% select(-data)
  cat("\n=== RQ3 test: animal-associated reference, air-enriched vs depleted ===\n")
  print(rq3_test)
  fwrite(rq3_test, out("RQ3_habitat_test.tsv"), sep = "\t")

  hab_cols <- c("Animal-associated" = "#D55E00", "Engineered/waste" = "#CC79A7",
                "Food" = "#F0E442", "Environmental" = "#0072B2", "Unknown" = "grey85")
}


## ============================================================================
## SECTION 11E. RQ3  PHYLOGENETIC SIGNAL OF AIR ASSOCIATION
## ============================================================================
## Tree: the 525 bacterial MAGs, from the GTDB-Tk bac120 marker alignment
## (HQ and MQ runs combined; identical sequences for MAGs in both), built with
## FastTree 2.1.11:  FastTree -wag -gamma bac_msa_525.faa > bac120_525_fasttree.nwk
## The alignment is written below if the tree file is missing.
## Tests (assessable sites only):
##   (a) Pagel's lambda of the within-site log2 ratio, raw and after removing
##       site means (lambda ~1 = trait follows the tree; ~0 = no signal),
##       also for DF and for VL + BB separately
##   (b) clustering of air-enriched MAGs: mean pairwise patristic distance
##       among them vs random sets of the same size (permutation)
F_TREE <- "data/bac120_525_fasttree.nwk"
F_MSA  <- c("data/gtdbtk_bac120_msa_HQ.fasta.gz", "data/gtdbtk_bac120_msa_MQ.fasta.gz")
suppressPackageStartupMessages({ library(ape); library(phytools) })
if (!file.exists(F_TREE) && all(file.exists(F_MSA))) {
  read_msa <- function(f) { x <- readLines(gzfile(f)); h <- grep("^>", x)
    setNames(vapply(seq_along(h), function(i) paste(x[(h[i] + 1):(c(h[-1], length(x) + 1)[i] - 1)], collapse = ""), ""),
             sub("^>(\\S+).*", "\\1", x[h])) }
  msa <- c(read_msa(F_MSA[2]), read_msa(F_MSA[1])); msa <- msa[!duplicated(names(msa), fromLast = TRUE)]
  ids <- bact_all$MAG_id; stopifnot(all(ids %in% names(msa)))
  writeLines(paste0(">", ids, "\n", msa[ids]), "bac_msa_525.faa")
  cat("\nTree file missing. Alignment written to bac_msa_525.faa; build the tree with:\n",
      "  FastTree -wag -gamma bac_msa_525.faa >", F_TREE, "\n")
}
if (file.exists(F_TREE)) {
  tr <- midpoint.root(read.tree(F_TREE))
  stopifnot(all(bact_all$MAG_id %in% tr$tip.label))
  tr_a <- keep.tip(tr, bact$MAG_id)

  lfc_raw <- setNames(bact$lfc_own, bact$MAG_id)
  lfc_res <- setNames(resid(lm(lfc_own ~ site, data = bact)), bact$MAG_id)
  lam_raw <- phylosig(tr_a, lfc_raw[tr_a$tip.label], method = "lambda", test = TRUE)
  lam_res <- phylosig(tr_a, lfc_res[tr_a$tip.label], method = "lambda", test = TRUE)

  D <- cophenetic(tr_a)
  enr <- bact$MAG_id[bact$air_group == "Air-enriched"]
  mpd <- function(x) mean(D[x, x][upper.tri(D[x, x])])
  obs_mpd  <- mpd(enr)
  null_mpd <- replicate(999, mpd(sample(tr_a$tip.label, length(enr))))
  ## same test within DF only (removes site from the null)
  df_ids <- bact$MAG_id[bact$site == "DF"]; enr_df <- intersect(enr, df_ids)
  obs_mpd_df  <- mpd(enr_df)
  null_mpd_df <- replicate(999, mpd(sample(df_ids, length(enr_df))))

  phylo_tbl <- tibble(
    test = c("Pagel lambda, log2 ratio", "Pagel lambda, log2 ratio minus site mean",
             "MPD air-enriched vs random (BB, DF, VL)", "MPD air-enriched vs random (DF only)"),
    n = c(length(lfc_raw), length(lfc_res), length(enr), length(enr_df)),
    statistic = c(lam_raw$lambda, lam_res$lambda, obs_mpd, obs_mpd_df),
    null_mean = c(NA, NA, mean(null_mpd), mean(null_mpd_df)),
    p = c(lam_raw$P, lam_res$P, (sum(null_mpd <= obs_mpd) + 1) / 1000,
          (sum(null_mpd_df <= obs_mpd_df) + 1) / 1000))
  ## Lambda by site subset (site means removed), so the signal does not rest on DF
  lam_sub <- bind_rows(lapply(list("DF", c("VL", "BB")), function(ss) {
    d  <- bact %>% filter(site %in% ss)
    x  <- setNames(d$lfc_own - ave(d$lfc_own, d$site), d$MAG_id)
    tt <- keep.tip(tr, d$MAG_id)
    l  <- phylosig(tt, x[tt$tip.label], method = "lambda", test = TRUE)
    tibble(test = paste0("Pagel lambda, log2 ratio minus site mean (", paste(ss, collapse = " + "), ")"),
           n = nrow(d), statistic = l$lambda, null_mean = NA_real_, p = l$P)
  }))
  phylo_tbl <- bind_rows(phylo_tbl, lam_sub)
  cat("\n=== RQ3 phylogenetic signal ===\n"); print(phylo_tbl)
  fwrite(phylo_tbl, out("RQ3_phylogenetic_signal.tsv"), sep = "\t")

  ## Fig. 4 (left): circular tree of the 525 bacterial MAGs in the style of
  ## iTOL colour ranges: black branches over pastel phylum sectors, an outer
  ## ring with each MAG's log2 air ratio (MK and RM MAGs: no segment) and
  ## points marking air-enriched MAGs. Legends sit in the corners.
  ## Fig. 4 (right): habitat of the closest species-level reference genome.
  ntip <- length(tr$tip.label)
  tip_phy <- bact_all$phylum[match(tr$tip.label, bact_all$MAG_id)]
  tip_phy <- ifelse(tip_phy %in% names(PHY_COLS), tip_phy, "Other")
  tip_grp <- bact_all$air_group[match(tr$tip.label, bact_all$MAG_id)]
  tip_lfc <- bact_all$lfc_own[match(tr$tip.label, bact_all$MAG_id)]
  tip_ok  <- bact_all$assessable[match(tr$tip.label, bact_all$MAG_id)]
  LFC_MAX  <- 4
  ring_pal <- colorRampPalette(c(COL_DEP, "white", COL_ENR))(101)
  ring_col <- ifelse(tip_ok, ring_pal[round((pmax(-LFC_MAX, pmin(LFC_MAX, tip_lfc)) + LFC_MAX) /
                                              (2 * LFC_MAX) * 100) + 1], NA)

  pdf(out("Fig4_lineage_and_reference_habitat.pdf"), width = 15, height = 8.5)
  layout(matrix(1:2, nrow = 1), widths = c(1.75, 1.2))
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  R <- max(node.depth.edgelength(tr))
  plot(tr, type = "fan", show.tip.label = FALSE, edge.width = 0.6, edge.color = "grey15",
       x.lim = c(-2.05, 1.25) * R, y.lim = c(-1.5, 1.85) * R)
  lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  a_tip <- atan2(lp$yy[1:ntip], lp$xx[1:ntip])
  r0  <- max(sqrt(lp$xx[1:ntip]^2 + lp$yy[1:ntip]^2)) * 1.01
  dth <- pi / ntip
  seg <- function(r1, r2, th, col) polygon(c(r1, r2, r2, r1) * cos(c(th - dth, th - dth, th + dth, th + dth)),
                                          c(r1, r2, r2, r1) * sin(c(th - dth, th - dth, th + dth, th + dth)),
                                          col = col, border = NA)
  for (k in seq_len(ntip)) {
    if (tip_phy[k] != "Other") seg(0.30 * r0, r0, a_tip[k], adjustcolor(PHY_COLS[[tip_phy[k]]], alpha.f = 0.55))
    if (!is.na(ring_col[k])) seg(r0 * 1.02, r0 * 1.12, a_tip[k], ring_col[k])
  }
  e <- which(tip_grp == "Air-enriched")
  points(r0 * 1.17 * cos(a_tip[e]), r0 * 1.17 * sin(a_tip[e]), pch = 21, bg = COL_ENR, col = "black", cex = 1.1)
  ## phylum legend (top left)
  ph_show <- setdiff(names(PHY_COLS), "Other")
  legend(x = -0.4 * R, y = 1.9 * R, xjust = 0.5, ncol = 2, legend = ph_show, fill = adjustcolor(PHY_COLS[ph_show], alpha.f = 0.55),
         border = "grey40", bty = "n", cex = BASE_FONT / 12, title = "Phylum", title.adj = 0, xpd = NA)
  ## ring key and enriched marker (bottom left)
  kx <- -2.05 * R; ky <- -1.30 * R; kw <- 0.6 * R; kh <- 0.06 * R
  text(kx, ky + kh + 0.09 * R, "log2 ratio, aerosol to source", adj = 0, cex = BASE_FONT / 12)
  rect(kx + (0:100) / 101 * kw, ky, kx + (1:101) / 101 * kw, ky + kh, col = ring_pal, border = NA)
  rect(kx, ky, kx + kw, ky + kh, border = "black")
  text(kx + c(0, 0.5, 1) * kw, ky - 0.08 * R, c(-LFC_MAX, 0, LFC_MAX), cex = BASE_FONT / 13)
  points(kx + 0.03 * R, ky - 0.25 * R, pch = 21, bg = COL_ENR, col = "black", cex = 1.5)
  text(kx + 0.09 * R, ky - 0.25 * R, "Air-enriched MAG", adj = 0, cex = BASE_FONT / 12)

  if (exists("hab_tab")) {
    hab_lv <- names(HAB_COLS)
    m <- hab_tab %>% mutate(habitat = factor(habitat, levels = hab_lv)) %>%
      select(air_group, habitat, pct) %>%
      pivot_wider(names_from = air_group, values_from = pct, values_fill = 0) %>%
      arrange(habitat)
    mm4 <- as.matrix(m[, c("Air-enriched", "Air-depleted")]); rownames(mm4) <- as.character(m$habitat)
    ng <- hab_tab %>% distinct(air_group, n_group)
    lab4 <- sprintf("%s\n(n = %d)", colnames(mm4), ng$n_group[match(colnames(mm4), ng$air_group)])
    par(mar = c(6, 6, 3, 14), las = 1)
    bx <- barplot(mm4, col = HAB_COLS[rownames(mm4)], border = "white", names.arg = rep("", ncol(mm4)),
                  ylab = "", cex.axis = BASE_FONT / 13, width = 1, space = c(0.3, 1.3),
                  xlim = c(0, 3.4), ylim = c(0, 100))
    mtext("MAGs with a species-level reference (%)", side = 2, line = 3.8, las = 0, cex = BASE_FONT / 13)
    mtext(lab4, side = 1, at = bx, line = 3, cex = BASE_FONT / 15)
    legend(x = par("usr")[2] + 0.35, y = 100, xpd = NA, ncol = 1, bty = "n", y.intersp = 1.3,
           legend = sub("Engineered/waste", "Engineered or waste", rev(rownames(mm4))),
           fill = rev(HAB_COLS[rownames(mm4)]), border = "white", cex = BASE_FONT / 12,
           title = "Reference habitat", title.adj = 0)
  }
  dev.off()
}


## ============================================================================
## SECTION 11F. RQ4  POPULATION GENOMICS OF AIRBORNE VS SOURCE POPULATIONS
## ============================================================================
## Dairy Farm only: the one pair with enough air coverage for population-level
## analysis. inStrain v1.12.0 profiled each DF sample against the 233 DF
## genomes (--pairing_filter all_reads, because minimap2 flagged almost all
## reads unpaired), then compared samples pairwise.
##
## Two questions, each with its own confound control:
##   (a) Are airborne populations less genetically diverse than their source
##       populations (an aerosolisation bottleneck)?
##       Confound: nucleotide diversity rises with coverage, and air coverage
##       is lower. Controlled by rarefying both samples to the same depth
##       (--rarefied_coverage) and by restricting to coverage-matched genomes.
##   (b) Are airborne populations genetically distinct from source populations?
##       Measured as popANI over positions covered in both samples.
##       Control: the soil-vs-soil comparison gives the normal difference
##       between two samples of the SAME habitat.
##       Confound: popANI tracks the fraction of genome compared; handled by
##       restricting to pairs with a matched compared fraction.
DIR_IS      <- "data/instrain_DF"             # folder holding the inStrain tables
IS_SAMPLES  <- c(air = "A_DF1", soil1 = "S_DF1", soil2 = "S_DF2")
IS_RAR      <- c(5, 10)                       # rarefaction depths rerun
MIN_COV_IS  <- 5                              # min coverage per sample
MIN_BREADTH <- 0.5                            # min breadth_minCov per sample
MIN_COMPARED<- 0.5                            # min percent_compared for popANI

read_gi <- function(f) if (file.exists(f)) fread(f) else NULL

## ---- (a) nucleotide diversity, air vs soil --------------------------------
nd_compare <- function(subdir, col, min_cov = MIN_COV_IS) {
  gi <- lapply(IS_SAMPLES, function(s)
    read_gi(file.path(DIR_IS, subdir, sprintf("%s_genome_info.tsv", s))))
  if (any(vapply(gi, is.null, TRUE))) return(NULL)
  keep <- lapply(gi, function(x)
    x$genome[x$coverage >= min_cov & x$breadth_minCov >= MIN_BREADTH])
  ok <- Reduce(intersect, keep[c("air", "soil1")])        # air vs soil1
  if (length(ok) < 5) return(NULL)
  d <- tibble(genome = ok,
              air  = gi$air[[col]][match(ok, gi$air$genome)],
              soil = gi$soil1[[col]][match(ok, gi$soil1$genome)],
              air_cov  = gi$air$coverage[match(ok, gi$air$genome)],
              soil_cov = gi$soil1$coverage[match(ok, gi$soil1$genome)]) %>%
    mutate(log2_cov_ratio = log2(air_cov / soil_cov), delta = air - soil)
  d
}

nd_test <- function(d, label) {
  if (is.null(d) || nrow(d) < 5) return(tibble())
  m <- d %>% filter(abs(log2_cov_ratio) < 0.5)            # coverage-matched
  tibble(analysis = label, n = nrow(d),
         median_air = median(d$air), median_soil = median(d$soil),
         n_air_lower = sum(d$air < d$soil),
         p_all = suppressWarnings(wilcox.test(d$air, d$soil, paired = TRUE)$p.value),
         rho_cov_confound = suppressWarnings(cor(d$log2_cov_ratio, d$delta, method = "spearman")),
         n_cov_matched = nrow(m),
         p_cov_matched = if (nrow(m) >= 5)
           suppressWarnings(wilcox.test(m$air, m$soil, paired = TRUE)$p.value) else NA_real_)
}

nd_raw <- nd_compare(".", "nucl_diversity")
nd_tbl <- bind_rows(
  nd_test(nd_raw, "raw nucleotide diversity (default profile)"),
  nd_test(nd_compare(".", "nucl_diversity_rarefied"), "rarefied to 50x (default)"),
  bind_rows(lapply(IS_RAR, function(r)
    bind_rows(
      nd_test(nd_compare(sprintf("rar%d", r), "nucl_diversity_rarefied"),
              sprintf("rarefied to %dx, genomes >=%dx", r, MIN_COV_IS)),
      nd_test(nd_compare(sprintf("rar%d", r), "nucl_diversity_rarefied", min_cov = 2 * r),
              sprintf("rarefied to %dx, genomes >=%dx in both", r, 2 * r))))))
if (nrow(nd_tbl)) {
  cat("\n=== RQ4a: nucleotide diversity, airborne vs source populations (DF) ===\n")
  print(nd_tbl, width = Inf)
  fwrite(nd_tbl, out("RQ4_nucleotide_diversity_tests.tsv"), sep = "\t")
}

## ---- (b) popANI, with the soil-vs-soil control ----------------------------
F_CMP <- file.path(DIR_IS, "compare_genomeWide_compare.tsv")
popani <- NULL
if (file.exists(F_CMP)) {
  cmp <- fread(F_CMP) %>%
    mutate(s1 = str_extract(name1, "^[A-Z]_[A-Z]{2}[0-9]"),
           s2 = str_extract(name2, "^[A-Z]_[A-Z]{2}[0-9]"),
           pair = paste(s1, s2, sep = "_"),
           type = ifelse(str_detect(s1, "^A_") | str_detect(s2, "^A_"),
                         "Air vs soil", "Soil vs soil"))
  key_air  <- paste(IS_SAMPLES[["air"]],   IS_SAMPLES[["soil1"]], sep = "_")
  key_ctrl <- paste(IS_SAMPLES[["soil1"]], IS_SAMPLES[["soil2"]], sep = "_")
  wide <- cmp %>% filter(pair %in% c(key_air, key_ctrl)) %>%
    select(genome, pair, popANI, percent_compared) %>%
    pivot_wider(names_from = pair, values_from = c(popANI, percent_compared))
  pc_a <- paste0("percent_compared_", key_air); pc_c <- paste0("percent_compared_", key_ctrl)
  pa_a <- paste0("popANI_", key_air);           pa_c <- paste0("popANI_", key_ctrl)
  popani <- wide %>%
    filter(.data[[pc_a]] >= MIN_COMPARED, .data[[pc_c]] >= MIN_COMPARED) %>%
    transmute(genome, air_soil = .data[[pa_a]], soil_soil = .data[[pa_c]],
              pc_air_soil = .data[[pc_a]], pc_soil_soil = .data[[pc_c]],
              delta = air_soil - soil_soil, delta_pc = pc_air_soil - pc_soil_soil)
  matched <- popani %>% filter(abs(delta_pc) < 0.1)
  popani_tbl <- tibble(
    n = nrow(popani),
    median_air_soil = median(popani$air_soil), median_soil_soil = median(popani$soil_soil),
    n_air_soil_lower = sum(popani$air_soil < popani$soil_soil),
    p_all = suppressWarnings(wilcox.test(popani$air_soil, popani$soil_soil, paired = TRUE)$p.value),
    rho_compared_confound = cor(popani$delta_pc, popani$delta, method = "spearman"),
    n_compared_matched = nrow(matched),
    p_compared_matched = if (nrow(matched) >= 5)
      suppressWarnings(wilcox.test(matched$air_soil, matched$soil_soil, paired = TRUE)$p.value) else NA_real_,
    same_pop_air_soil = sum(popani$air_soil >= 0.99999),
    same_pop_soil_soil = sum(popani$soil_soil >= 0.99999))
  cat("\n=== RQ4b: popANI, airborne vs source, with soil-vs-soil control (DF) ===\n")
  print(popani_tbl, width = Inf)
  fwrite(popani_tbl, out("RQ4_popANI_tests.tsv"), sep = "\t")
  fwrite(popani,     out("RQ4_popANI_per_genome.tsv"), sep = "\t")
}

## ---- Supplementary Figure S4 ---------------------------------------------------------------
## Both panels are paired, one point per genome, compared against the 1:1 line:
##   left  : nucleotide diversity in air vs soil (strictest depth match)
##   right : divergence between air and soil populations vs between the two
##           soil populations (within-habitat baseline)
## Points are coloured by air association class, as in Figure 1.
grp_cols5 <- c("Air-enriched" = COL_ENR, "Intermediate" = COL_INT, "Air-depleted" = COL_DEP)
add_grp <- function(d) d %>%
  mutate(MAG_id = sub("\\.fa(sta)?$", "", genome)) %>%
  left_join(bact %>% select(MAG_id, air_group), by = "MAG_id") %>%
  mutate(air_group = factor(air_group, levels = names(grp_cols5)))
theme5 <- theme_classic(base_size = BASE_FONT) +
  theme(axis.text = element_text(size = BASE_FONT, colour = "black"),
        axis.title = element_text(size = BASE_FONT * 1.05),
        legend.position = "top", legend.text = element_text(size = BASE_FONT),
        plot.margin = margin(5, 15, 5, 5))
if (!is.null(nd_raw) || !is.null(popani)) {
  pl <- list()
  if (!is.null(nd_raw)) {
    ## strictest depth-matched comparison: rarefied to the largest depth in
    ## IS_RAR, genomes with at least twice that coverage in both samples.
    ## All other comparisons are in RQ4_nucleotide_diversity_tests.tsv.
    r_max <- max(IS_RAR)
    best <- nd_compare(sprintf("rar%d", r_max), "nucl_diversity_rarefied", min_cov = 2 * r_max)
    use  <- add_grp(if (!is.null(best)) best else nd_raw)
    lim  <- range(c(use$soil, use$air)) * c(0.9, 1.05)
    pl$a <- ggplot(use, aes(x = soil, y = air, fill = air_group)) +
      geom_abline(slope = 1, intercept = 0, colour = "black", linewidth = 0.4) +
      geom_point(shape = 21, size = 4.5, alpha = 0.8, colour = "black", stroke = 0.3) +
      scale_fill_manual(values = grp_cols5, drop = FALSE, name = NULL) +
      scale_x_continuous(limits = lim, breaks = seq(0, 0.02, 0.004)) +
      scale_y_continuous(limits = lim, breaks = seq(0, 0.02, 0.004)) +
      coord_equal() +
      labs(x = "Nucleotide diversity in soil", y = "Nucleotide diversity in air") +
      theme5
  }
  if (!is.null(popani)) {
    ## divergence = 1 - popANI per 10 kb (popANI substitutions per 10 kb);
    ## pseudo-log axes so identical populations (0) stay on the plot
    pd <- add_grp(popani) %>%
      mutate(div_air = (1 - air_soil) * 1e4, div_soil = (1 - soil_soil) * 1e4)
    plim <- c(0, max(c(pd$div_air, pd$div_soil)) * 1.1)
    pl$b <- ggplot(pd, aes(x = div_soil, y = div_air, fill = air_group)) +
      geom_abline(slope = 1, intercept = 0, colour = "black", linewidth = 0.4) +
      geom_point(shape = 21, size = 4.5, alpha = 0.8, colour = "black", stroke = 0.3) +
      scale_fill_manual(values = grp_cols5, drop = FALSE, name = NULL) +
      scale_x_continuous(trans = pseudo_log_trans(sigma = 0.01, base = 10), limits = plim,
                         breaks = c(0, 0.01, 0.1, 1, 10), labels = c("0", "0.01", "0.1", "1", "10")) +
      scale_y_continuous(trans = pseudo_log_trans(sigma = 0.01, base = 10), limits = plim,
                         breaks = c(0, 0.01, 0.1, 1, 10), labels = c("0", "0.01", "0.1", "1", "10")) +
      coord_equal() +
      labs(x = "Divergence, soil vs soil\n(popANI substitutions per 10 kb)",
           y = "Divergence, air vs soil\n(popANI substitutions per 10 kb)") +
      theme5
  }
  pdf(out("FigS4_population_genomics.pdf"), width = 13, height = 7)
  print(if (length(pl) == 2) (pl$a + pl$b + plot_layout(guides = "collect") &
                                theme(legend.position = "top")) else pl[[1]])
  dev.off()
}


## ============================================================================
## SECTION 11G. READ RECRUITMENT, SAMPLE COMPOSITION AND FIELD CONTROLS
## ============================================================================
## (a) Per-sample table: reads, % recruited to the 532 MAGs, and kraken2
##     composition (from the companion read-based analysis).
## (b) Field controls: collection fluid placed in the sampler and left open at
##     the site for the same duration without running the pump (passive
##     deposition plus handling). Libraries were obtained at DF and RM only.
##     Controls were mapped with the same minimap2 settings; reads per contig
##     come from samtools idxstats, summed per genome.
## Sample naming: the VL air library is labelled AC_VL2 in the read-based
## pipeline and A_VL2 here; it is the VL air sample, not a control.
F_IDX   <- "data/controls/genome_mapped_reads_17samples.tsv"
DIR_KR  <- "data/kraken_reports"            # gzipped Kraken2 reports, one per library
CONTROLS <- c("AC_DF2", "AC_RM2")
KR_NAME  <- c(A_VL2 = "AC_VL2")                 # kraken file name differences

read_flag_primary <- function(f) {
  x <- readLines(f)
  c(primary = as.numeric(sub(" .*", "", grep(" primary$", x, value = TRUE)[1])),
    mapped  = as.numeric(sub(" .*", "", grep(" primary mapped", x, value = TRUE)[1])))
}

if (file.exists(F_IDX)) {
  idx <- fread(F_IDX, header = FALSE,
               col.names = c("genome", "sample", "mapped_reads", "genome_length")) %>%
    filter(genome != "genome") %>%
    mutate(across(c(mapped_reads, genome_length), as.numeric))

  ctrl_depth <- bind_rows(lapply(CONTROLS, function(cs) {
    f <- file.path("data/controls", paste0(cs, ".flagstat.txt"))
    v <- read_flag_primary(f); tibble(sample = cs, primary = v[["primary"]], primary_mapped = v[["mapped"]])
  }))
  depth <- bind_rows(reads %>% select(sample, primary, primary_mapped), ctrl_depth)

  ## ---- (a) per-sample recruitment and composition --------------------------
  kr_one <- function(smp) {
    fn <- if (smp %in% names(KR_NAME)) KR_NAME[[smp]] else smp
    f <- file.path(DIR_KR, paste0(fn, ".kreport.gz"))
    if (!file.exists(f)) return(NULL)
    ## Read the gzipped report through a base-R connection. fread() can take a
    ## .gz path directly only if the optional R.utils package is installed, and
    ## this keeps the dependency list unchanged.
    k <- fread(text = readLines(gzfile(f, "rt")),
               header = FALSE, sep = "\t", quote = "",
               col.names = c("pct", "clade", "taxon", "rank", "taxid", "name")) %>%
      mutate(name = trimws(name))
    tot <- sum(k$clade[k$name %in% c("unclassified", "root")])
    g <- function(n) { v <- k$clade[k$name == n]; if (length(v)) 100 * v[1] / tot else 0 }
    tibble(sample = smp, kraken_pairs = tot, pct_unclassified = g("unclassified"),
           pct_bacteria = g("Bacteria"), pct_archaea = g("Archaea"), pct_fungi = g("Fungi"),
           pct_human = g("Homo sapiens"), pct_viruses = g("Viruses"))
  }
  comp <- bind_rows(lapply(depth$sample, kr_one))
  if (nrow(comp) == 0) {
    warning("No Kraken2 reports found in ", DIR_KR,
            "; the composition columns of SuppTable_read_recruitment_and_composition ",
            "will be NA.", call. = FALSE)
    comp <- tibble(sample = character(), kraken_pairs = numeric(),
                   pct_unclassified = numeric(), pct_bacteria = numeric(),
                   pct_archaea = numeric(), pct_fungi = numeric(),
                   pct_human = numeric(), pct_viruses = numeric())
  }
  recruit <- depth %>%
    mutate(type = case_when(sample %in% CONTROLS ~ "Field control",
                            str_detect(sample, "^A_") ~ "Air",
                            str_detect(sample, "^S_") ~ "Soil", TRUE ~ "Water"),
           site = substr(sub("^[A-Z]+_", "", sample), 1, 2),
           pct_mapped_to_MAGs = 100 * primary_mapped / primary,
           assessable_air = ifelse(type == "Air", pct_mapped_to_MAGs >= MIN_AIR_MAPPED_PCT, NA)) %>%
    left_join(comp, by = "sample") %>%
    mutate(pct_bacterial_reads_recruited = ifelse(pct_bacteria > 0,
                                                  pmin(100, 100 * pct_mapped_to_MAGs / pct_bacteria), NA)) %>%
    arrange(factor(type, levels = c("Air", "Field control", "Soil", "Water")), site)
  cat("\n=== Read recruitment and composition per sample ===\n")
  print(recruit %>% select(sample, type, primary, pct_mapped_to_MAGs, pct_bacteria,
                           pct_fungi, pct_human, pct_unclassified) %>%
          mutate(across(where(is.double), ~ round(.x, 2))), n = Inf)
  fwrite(recruit, out("SuppTable_read_recruitment_and_composition.tsv"), sep = "\t")

  ## ---- (b) field controls vs air ------------------------------------------
  rpkm <- idx %>% left_join(depth %>% select(sample, primary), by = "sample") %>%
    mutate(rpkm = mapped_reads / (genome_length / 1e3) / (primary / 1e6)) %>%
    select(genome, sample, mapped_reads, rpkm)
  ctrl_pairs <- c(AC_DF2 = "A_DF1", AC_RM2 = "A_RM1")

  ctrl_top <- bind_rows(lapply(CONTROLS, function(cs) {
    rpkm %>% filter(sample == cs) %>% mutate(share_pct = 100 * mapped_reads / sum(mapped_reads)) %>%
      slice_max(mapped_reads, n = 15) %>%
      left_join(mags %>% select(genome = MAG_id, site, domain, gtdb_classification), by = "genome")
  }))
  fwrite(ctrl_top, out("Controls_top_genomes.tsv"), sep = "\t")

  ctrl_sim <- bind_rows(lapply(CONTROLS, function(cs) {
    w <- rpkm %>% select(genome, sample, rpkm) %>% pivot_wider(names_from = sample, values_from = rpkm)
    others <- setdiff(names(w), c("genome", CONTROLS))
    tibble(control = cs, sample = others,
           spearman_rho = vapply(others, function(o) cor(w[[cs]], w[[o]], method = "spearman"), 0))
  }))
  fwrite(ctrl_sim, out("Controls_profile_similarity.tsv"), sep = "\t")

  enr_df <- bact %>% filter(air_group == "Air-enriched", site == "DF") %>% pull(MAG_id)
  ctrl_cmp <- rpkm %>% filter(sample %in% c("A_DF1", "AC_DF2")) %>%
    select(genome, sample, rpkm) %>% pivot_wider(names_from = sample, values_from = rpkm) %>%
    mutate(air_enriched = genome %in% enr_df,
           in_bact = genome %in% bact_all$MAG_id,
           air_over_control = A_DF1 / AC_DF2)
  ctrl_sum <- ctrl_cmp %>% filter(air_enriched) %>%
    summarise(n = n(), median_ratio = median(air_over_control, na.rm = TRUE),
              min_ratio = min(air_over_control, na.rm = TRUE),
              n_zero_in_control = sum(AC_DF2 == 0))
  cat("\n=== DF air-enriched bacterial MAGs: active air sample vs passive field control ===\n")
  print(ctrl_sum)
  fwrite(ctrl_cmp, out("Controls_DF_air_vs_control.tsv"), sep = "\t")
  fwrite(ctrl_sum, out("Controls_DF_air_vs_control_summary.tsv"), sep = "\t")

  s1_cols <- c("Air-enriched" = COL_ENR, "Intermediate" = COL_INT, "Air-depleted" = COL_DEP)
  s1_df <- ctrl_cmp %>% filter(in_bact, str_detect(genome, "^DF\\.")) %>%
    left_join(bact %>% select(MAG_id, air_group), by = c("genome" = "MAG_id")) %>%
    mutate(air_group = factor(air_group, levels = names(s1_cols)),
           x = AC_DF2 + 1e-3, y = A_DF1 + 1e-3) %>%
    arrange(air_group == "Air-enriched")
  p_ctrl <- ggplot(s1_df, aes(x = x, y = y, fill = air_group)) +
    geom_abline(slope = 1, intercept = 0, colour = "black", linewidth = 0.4) +
    geom_point(shape = 21, size = 3.5, alpha = 0.8, colour = "black", stroke = 0.3) +
    scale_x_log10(labels = label_number(drop0trailing = TRUE)) +
    scale_y_log10(labels = label_number(drop0trailing = TRUE)) + coord_equal() +
    scale_fill_manual(values = s1_cols, name = NULL) +
    labs(x = "Field control (RPKM)", y = "Aerosol metagenome (RPKM)") +
    theme_classic(base_size = BASE_FONT) +
    theme(axis.text = element_text(size = BASE_FONT, colour = "black"),
          axis.title = element_text(size = BASE_FONT * 1.05),
          legend.position = "top", legend.text = element_text(size = BASE_FONT))
  pdf(out("FigS1_field_controls.pdf"), width = 7.5, height = 7.5)
  print(p_ctrl)
  dev.off()
}


## ============================================================================
## SECTION 12. SUPPLEMENTARY FIGURE S2  Genome novelty (closest-reference ANI, GTDB-Tk)
## ============================================================================
## Missing ANI = no reference genome close enough for GTDB-Tk to report one.
## The CDF uses the full group size as denominator, so MAGs without an ANI
## value are not silently dropped.
ANI_SPECIES <- 95
nov <- bact %>%
  filter(air_group != "Intermediate") %>%
  mutate(air_group = factor(air_group, levels = c("Air-enriched", "Air-depleted")))

nov_summary <- nov %>% group_by(air_group) %>%
  summarise(n = n(),
            n_missing_ANI = sum(is.na(gtdb_closest_ani)),
            pct_missing   = 100 * n_missing_ANI / n,
            n_with_ANI    = sum(!is.na(gtdb_closest_ani)),
            n_below_95    = sum(gtdb_closest_ani < ANI_SPECIES, na.rm = TRUE),
            pct_below_95_of_assigned = ifelse(n_with_ANI > 0, 100 * n_below_95 / n_with_ANI, NA_real_),
            .groups = "drop")
cat("\n=== Genome novelty ===\n"); print(nov_summary)
fwrite(nov_summary, out("FigS2_genome_novelty_summary.tsv"), sep = "\t")

## One stacked bar per group: share of MAGs whose closest GTDB reference is
## the same species (ANI >= 95%), a related species (ANI < 95%), or not
## reported (no reference close enough for GTDB-Tk to report an ANI).
nov_cls <- c("Same species (ANI of 95% or more)" = "#08519C",
             "Related species (ANI < 95%)" = "#9ECAE1",
             "No reported ANI" = "white")
nov_bar <- nov %>%
  mutate(cls = case_when(is.na(gtdb_closest_ani) ~ "No reported ANI",
                         gtdb_closest_ani >= ANI_SPECIES ~ "Same species (ANI of 95% or more)",
                         TRUE ~ "Related species (ANI < 95%)"),
         cls = factor(cls, levels = rev(names(nov_cls)))) %>%
  count(air_group, cls) %>% group_by(air_group) %>%
  mutate(N = sum(n), pct = 100 * n / N) %>% ungroup() %>%
  mutate(xlab = factor(sprintf("%s\n(n = %d)", air_group, N),
                       levels = unique(sprintf("%s\n(n = %d)", air_group, N)[order(air_group)])))
fwrite(nov_bar, out("FigS2_genome_novelty_classes.tsv"), sep = "\t")
p_nov <- ggplot(nov_bar, aes(x = xlab, y = pct, fill = cls)) +
  geom_col(colour = "black", linewidth = 0.3, width = 0.6) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = BASE_FONT * 0.3,
            colour = ifelse(nov_bar$cls == names(nov_cls)[1], "white", "black")) +
  scale_fill_manual(values = nov_cls, breaks = names(nov_cls), name = "Closest GTDB reference") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "MAGs (%)") +
  theme_classic(base_size = BASE_FONT) +
  theme(axis.text = element_text(size = BASE_FONT, colour = "black"),
        axis.title = element_text(size = BASE_FONT * 1.05),
        legend.text = element_text(size = BASE_FONT), legend.title = element_text(size = BASE_FONT))
pdf(out("FigS2_genome_novelty.pdf"), width = 8, height = 6)
print(p_nov)
dev.off()


## ============================================================================
## SECTION 14. RECORD
## ============================================================================
writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
cat("\n=== Done. All outputs are in:", file.path(BASE, OUTDIR), "===\n")
print(list.files(OUTDIR))
