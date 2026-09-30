# Mitochondrial and nuclear genotype comparisons
# Compare Gen5 communities with genotype- and stage-specific Gen0 profiles.
# Fit a mixed model with vial intercepts and a vial-mean sensitivity model.
# Research code example; processed research data are not included.
# Run from the folder containing README.md. See CHANGES.md before using outputs.

# Input and output settings: edit these paths when data become available.
INPUT_RDS <- file.path("data", "06_phyloseq_clean_CHAP1_NOHOST.rds")
OUTPUT_ROOT <- "outputs"
PDF_DEVICE <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

if (!file.exists(INPUT_RDS)) {
  stop("Processed research data are not included. Add the phyloseq RDS at ",
       INPUT_RDS, " or edit INPUT_RDS at the top of this script.", call. = FALSE)
}

required_packages <- c("phyloseq", "dplyr", "tidyr", "stringr", "vegan", "ggplot2", "emmeans", "nlme", "readr", "tibble", "scales", "openxlsx", "car")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install the required packages: ", paste(missing_packages, collapse = ", "),
       ". See README.md for installation instructions.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(vegan)
  library(ggplot2)
  library(emmeans)
  library(nlme)
  library(readr)
  library(tibble)
  library(scales)
  library(openxlsx)
})

set.seed(1234)

# Use car through its namespace to avoid masking dplyr functions.
# car is required; car:: calls avoid attaching its recode function.

# ------------------------------
# 2) Create organized output folders
# ------------------------------

main_outdir <- file.path(OUTPUT_ROOT, "Fig4_mitonuclear_divergence")

dir.create(main_outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(main_outdir, "input_processed"), showWarnings = FALSE)
dir.create(file.path(main_outdir, "figure_files"), showWarnings = FALSE)
dir.create(file.path(main_outdir, "figure_source_data"), showWarnings = FALSE)
dir.create(file.path(main_outdir, "statistics"), showWarnings = FALSE)

fig_dir   <- file.path(main_outdir, "figure_files")
data_dir  <- file.path(main_outdir, "figure_source_data")
stat_dir  <- file.path(main_outdir, "statistics")
input_dir <- file.path(main_outdir, "input_processed")

# ------------------------------
# 3) User options
# ------------------------------

SAVE_TABLES <- TRUE
MIN_GEN0_VIALS <- 1
DODGE_W <- 0.70

# ------------------------------
# 4) Helper functions
# ------------------------------

to_rel <- function(psx) {
  transform_sample_counts(psx, function(x) if (sum(x) > 0) x / sum(x) else x)
}

p_to_star <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ "ns"
  )
}

safe_tab_to_csv <- function(x, filename) {
  write_csv(
    as.data.frame(x) %>% rownames_to_column("term"),
    filename
  )
}

pretty_genotype <- function(x) {
  dplyr::recode(
    as.character(x),
    "oreore" = "(ore);OreR",
    "simore" = "(simw501);OreR",
    "oreaut" = "(ore);Aut",
    "simaut" = "(simw501);Aut",
    .default = as.character(x)
  )
}

# ------------------------------
# 5) Load processed phyloseq object
# ------------------------------

ps.clean <- readRDS(INPUT_RDS)

if (!inherits(ps.clean, "phyloseq")) {
  stop("INPUT_RDS must contain a phyloseq object.", call. = FALSE)
}
required_metadata <- c("genotype", "generation", "life_stage", "vial")
missing_metadata <- setdiff(required_metadata, names(data.frame(sample_data(ps.clean))))
if (length(missing_metadata) > 0) {
  stop("Missing sample metadata: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
}

saveRDS(
  ps.clean,
  file = file.path(input_dir, "06_phyloseq_clean_CHAP1_NOHOST.rds")
)

# ============================================================
# PART A: SUPPLEMENTARY FIGURE S2
# Generation 0 host-associated microbiome composition
# Bray-Curtis PCoA by mito-nuclear genotype
#
# ============================================================

genotype_labels_parsed <- c(
  "(ore);OreR"     = "'(ore);OreR'",
  "(simw501);OreR" = "'(simw'^501*');OreR'",
  "(ore);Aut"      = "'(ore);Aut'",
  "(simw501);Aut"  = "'(simw'^501*');Aut'"
)

ps.gen0 <- subset_samples(
  ps.clean,
  generation == "Generation_0" &
    life_stage %in% c("Larvae", "Female") &
    genotype %in% c("oreore", "simore", "oreaut", "simaut")
)

ps.gen0 <- prune_samples(sample_sums(ps.gen0) > 0, ps.gen0)
ps.gen0 <- prune_taxa(taxa_sums(ps.gen0) > 0, ps.gen0)

cat("\nSamples included in Supplementary Fig. S2:\n")
print(
  data.frame(sample_data(ps.gen0)) %>%
    count(genotype, life_stage, vial)
)

cat("\nNumber of Gen0 host samples:", nsamples(ps.gen0), "\n")
cat("Number of taxa:", ntaxa(ps.gen0), "\n")

ps.gen0.ra <- transform_sample_counts(ps.gen0, function(x) x / sum(x))
ps.gen0.ra <- prune_taxa(taxa_sums(ps.gen0.ra) > 0, ps.gen0.ra)

otu.gen0 <- as(otu_table(ps.gen0.ra), "matrix")
if (taxa_are_rows(ps.gen0.ra)) otu.gen0 <- t(otu.gen0)

meta.gen0 <- data.frame(sample_data(ps.gen0.ra)) %>%
  rownames_to_column("SampleID") %>%
  mutate(
    genotype = as.character(genotype),
    life_stage = as.character(life_stage),
    vial = as.character(vial),

    Genotype = pretty_genotype(genotype),
    Genotype = factor(
      Genotype,
      levels = c(
        "(ore);OreR",
        "(simw501);OreR",
        "(ore);Aut",
        "(simw501);Aut"
      )
    ),

    Stage = dplyr::recode(
      life_stage,
      "Larvae" = "Larvae",
      "Female" = "Adult females"
    ),
    Stage = factor(Stage, levels = c("Larvae", "Adult females")),

    Mito = case_when(
      genotype %in% c("oreore", "oreaut") ~ "(ore)",
      genotype %in% c("simore", "simaut") ~ "(simw501)",
      TRUE ~ NA_character_
    ),
    Nuclear = case_when(
      genotype %in% c("oreore", "simore") ~ "OreR",
      genotype %in% c("oreaut", "simaut") ~ "Aut",
      TRUE ~ NA_character_
    ),
    Mito = factor(Mito, levels = c("(ore)", "(simw501)")),
    Nuclear = factor(Nuclear, levels = c("OreR", "Aut"))
  )

otu.gen0 <- otu.gen0[meta.gen0$SampleID, , drop = FALSE]

bray.gen0 <- vegan::vegdist(otu.gen0, method = "bray")

pcoa.gen0 <- cmdscale(bray.gen0, eig = TRUE, k = 2)

ord_df <- as.data.frame(pcoa.gen0$points) %>%
  rownames_to_column("SampleID") %>%
  rename(PCoA1 = V1, PCoA2 = V2) %>%
  left_join(meta.gen0, by = "SampleID")

eig <- pcoa.gen0$eig
eig_pos <- eig[eig > 0]

pc1_var <- round(100 * eig_pos[1] / sum(eig_pos), 1)
pc2_var <- round(100 * eig_pos[2] / sum(eig_pos), 1)

cat("\nPCoA variance explained:\n")
cat("PCoA1:", pc1_var, "%\n")
cat("PCoA2:", pc2_var, "%\n")

cat("\nOverall Gen0 PERMANOVA: Bray-Curtis ~ Mito * Nuclear * Stage\n")
# Pools share vials: the original free-permutation design needs review.
# by = "margin" respects term hierarchy; inspect which terms are tested.
perm_all_s2 <- vegan::adonis2(
  bray.gen0 ~ Mito * Nuclear * Stage,
  data = meta.gen0,
  permutations = 999,
  by = "margin"
)
print(perm_all_s2)

cat("\nOverall Gen0 betadisper by Genotype x Stage\n")
bd_all_s2 <- vegan::betadisper(
  bray.gen0,
  group = interaction(meta.gen0$Genotype, meta.gen0$Stage, drop = TRUE)
)
bd_all_perm_s2 <- vegan::permutest(bd_all_s2, permutations = 999)
print(bd_all_perm_s2)

stage_stats_s2 <- list()

for (stg in levels(meta.gen0$Stage)) {

  keep <- meta.gen0$Stage == stg

  otu_sub <- otu.gen0[keep, , drop = FALSE]
  meta_sub <- meta.gen0[keep, , drop = FALSE] %>% droplevels()

  bray_sub <- vegan::vegdist(otu_sub, method = "bray")

  cat("\nStage-specific Gen0 PERMANOVA:", stg, "\n")
  perm_sub <- vegan::adonis2(
    bray_sub ~ Mito * Nuclear,
    data = meta_sub,
    permutations = 999,
    by = "margin"
  )
  print(perm_sub)

  cat("\nStage-specific Gen0 betadisper:", stg, "\n")
  bd_sub <- vegan::betadisper(
    bray_sub,
    group = interaction(meta_sub$Mito, meta_sub$Nuclear, drop = TRUE)
  )
  bd_sub_perm <- vegan::permutest(bd_sub, permutations = 999)
  print(bd_sub_perm)

  stage_stats_s2[[stg]] <- list(
    permanova = perm_sub,
    betadisper = bd_sub_perm
  )
}

genotype_cols <- c(
  "(ore);OreR" = "#0072B2",
  "(simw501);OreR" = "#D55E00",
  "(ore);Aut" = "#009E73",
  "(simw501);Aut" = "#CC79A7"
)

genotype_shapes <- c(
  "(ore);OreR" = 21,
  "(simw501);OreR" = 24,
  "(ore);Aut" = 22,
  "(simw501);Aut" = 23
)

S2_gen0_pcoa <- ggplot(
  ord_df,
  aes(x = PCoA1, y = PCoA2, fill = Genotype, shape = Genotype)
) +
  geom_hline(yintercept = 0, color = "grey90", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "grey90", linewidth = 0.35) +
  geom_point(size = 3.8, color = "black", stroke = 0.45, alpha = 0.92) +
  facet_wrap(~ Stage, nrow = 1) +
  scale_fill_manual(
    values = genotype_cols, breaks = names(genotype_cols),
    labels = parse(text = genotype_labels_parsed[names(genotype_cols)]),
    name = "Mito-nuclear genotype"
  ) +
  scale_shape_manual(
    values = genotype_shapes, breaks = names(genotype_shapes),
    labels = parse(text = genotype_labels_parsed[names(genotype_shapes)]),
    name = "Mito-nuclear genotype"
  ) +
  guides(
    shape = "none",
    fill = guide_legend(nrow = 1, byrow = TRUE,
                        override.aes = list(shape = unname(genotype_shapes), color = "black", size = 3.8, alpha = 0.95))
  ) +
  labs(x = paste0("PCoA1 (", pc1_var, "%)"), y = paste0("PCoA2 (", pc2_var, "%)")) +
  coord_equal() +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
    strip.text = element_text(face = "bold", size = 13.5, color = "black"),
    axis.title = element_text(face = "bold", color = "black", size = 14),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.text = element_text(color = "black", size = 11),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 11.5),
    legend.text = element_text(size = 10.5, color = "black"),
    legend.key.size = unit(0.55, "cm"),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.65),
    axis.line = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 10, r = 12, b = 10, l = 12)
  )

S2_gen0_pcoa

ggsave(file.path(fig_dir, "Supplementary_Figure_S2_Gen0_host_PCoA.pdf"), S2_gen0_pcoa,
       width = 9.8, height = 5.4, units = "in", device = PDF_DEVICE, bg = "white")
ggsave(file.path(fig_dir, "Supplementary_Figure_S2_Gen0_host_PCoA.png"), S2_gen0_pcoa,
       width = 9.8, height = 5.4, units = "in", dpi = 600, bg = "white")
ggsave(file.path(fig_dir, "Supplementary_Figure_S2_Gen0_host_PCoA.jpeg"), S2_gen0_pcoa,
       width = 9.8, height = 5.4, units = "in", dpi = 600, bg = "white", device = "jpeg")

write.csv(ord_df, file.path(data_dir, "Supplementary_Figure_S2_Gen0_host_PCoA_coordinates.csv"), row.names = FALSE)
write.csv(as.data.frame(perm_all_s2) %>% rownames_to_column("term"),
          file.path(stat_dir, "Supplementary_Figure_S2_Gen0_PERMANOVA_overall.csv"), row.names = FALSE)
write.csv(as.data.frame(bd_all_perm_s2$tab) %>% rownames_to_column("term"),
          file.path(stat_dir, "Supplementary_Figure_S2_Gen0_betadisper_overall.csv"), row.names = FALSE)

for (stg in names(stage_stats_s2)) {
  safe_stg <- gsub(" ", "_", stg)
  write.csv(as.data.frame(stage_stats_s2[[stg]]$permanova) %>% rownames_to_column("term"),
            file.path(stat_dir, paste0("Supplementary_Figure_S2_Gen0_PERMANOVA_", safe_stg, ".csv")), row.names = FALSE)
  write.csv(as.data.frame(stage_stats_s2[[stg]]$betadisper$tab) %>% rownames_to_column("term"),
            file.path(stat_dir, paste0("Supplementary_Figure_S2_Gen0_betadisper_", safe_stg, ".csv")), row.names = FALSE)
}

# ============================================================
# PART B: FIGURE 4 MAIN ANALYSIS
# Gen0 host composition + Gen5 divergence from Gen0
#
# ============================================================

ps.bc <- subset_samples(
  ps.clean,
  !is.na(generation) &
    generation %in% c("Generation_0", "Generation_5") &
    life_stage %in% c("Larvae", "Female") &
    genotype %in% c("oreore", "oreaut", "simore", "simaut")
)

ps.bc <- prune_samples(sample_sums(ps.bc) > 0, ps.bc)
ps.bc.ra <- to_rel(ps.bc)

cat("\n================ 1) BASIC SUBSET SUMMARY ================\n")
cat("Number of samples:", nsamples(ps.bc.ra), "\n")
cat("Number of taxa:", ntaxa(ps.bc.ra), "\n\n")

meta0 <- data.frame(sample_data(ps.bc.ra))

cat("Samples by genotype:\n")
print(table(meta0$genotype), quote = FALSE)
cat("\nSamples by generation:\n")
print(table(meta0$generation), quote = FALSE)
cat("\nSamples by life stage:\n")
print(table(meta0$life_stage), quote = FALSE)

otu.bc <- as(otu_table(ps.bc.ra), "matrix")
if (taxa_are_rows(ps.bc.ra)) otu.bc <- t(otu.bc)

meta.bc <- data.frame(sample_data(ps.bc.ra)) %>%
  mutate(
    SampleID = sample_names(ps.bc.ra),
    Mito_raw = case_when(
      genotype %in% c("oreore", "oreaut") ~ "ore",
      genotype %in% c("simore", "simaut") ~ "simw501",
      TRUE ~ NA_character_
    ),
    Mito = factor(ifelse(Mito_raw == "ore", "(ore)", "(simw501)"), levels = c("(ore)", "(simw501)")),
    Nuclear = case_when(
      genotype %in% c("oreore", "simore") ~ "OreR",
      genotype %in% c("oreaut", "simaut") ~ "Aut",
      TRUE ~ NA_character_
    ),
    Nuclear = factor(Nuclear, levels = c("OreR", "Aut")),
    Stage = factor(life_stage, levels = c("Larvae", "Female")),
    vial = as.character(vial),
    Group = paste(genotype, life_stage, sep = "_"),
    vial_id = paste(genotype, life_stage, vial, sep = "__")
  )

otu.bc <- otu.bc[meta.bc$SampleID, , drop = FALSE]

group_structure <- meta.bc %>%
  count(genotype, generation, life_stage, vial) %>%
  arrange(genotype, life_stage, generation, vial)

write_csv(group_structure, file.path(data_dir, "Figure4_sample_grouping_structure.csv"))

gen0_ids <- meta.bc %>% filter(generation == "Generation_0") %>% pull(SampleID)
otu_gen0 <- otu.bc[gen0_ids, , drop = FALSE]

meta_gen0 <- meta.bc %>%
  filter(SampleID %in% gen0_ids) %>%
  mutate(
    Mito = factor(Mito, levels = c("(ore)", "(simw501)")),
    Nuclear = factor(Nuclear, levels = c("OreR", "Aut")),
    Stage = factor(Stage, levels = c("Larvae", "Female"))
  )

bray_gen0 <- vegan::vegdist(otu_gen0, method = "bray")

perm_gen0_all <- vegan::adonis2(
  bray_gen0 ~ Mito * Nuclear * Stage,
  data = meta_gen0,
  permutations = 999,
  by = "margin"
)

cat("\nGen0 host composition PERMANOVA: Bray-Curtis ~ Mito * Nuclear * Stage\n")
print(perm_gen0_all)

meta_gen0 <- meta_gen0 %>%
  mutate(GenotypeStage = interaction(Mito, Nuclear, Stage, drop = TRUE))

bd_gen0_all <- vegan::permutest(
  vegan::betadisper(bray_gen0, group = meta_gen0$GenotypeStage),
  permutations = 999
)

cat("\nGen0 dispersion check by Mito x Nuclear x Stage group\n")
print(bd_gen0_all)

gen0_stage_results <- list()

for (stg in c("Larvae", "Female")) {

  keep_ids <- meta_gen0 %>% filter(Stage == stg) %>% pull(SampleID)
  otu_sub <- otu_gen0[keep_ids, , drop = FALSE]
  meta_sub <- meta_gen0 %>% filter(SampleID %in% keep_ids) %>% droplevels()
  bray_sub <- vegan::vegdist(otu_sub, method = "bray")

  perm_sub <- vegan::adonis2(bray_sub ~ Mito * Nuclear, data = meta_sub, permutations = 999, by = "margin")

  bd_sub <- vegan::permutest(
    vegan::betadisper(bray_sub, interaction(meta_sub$Mito, meta_sub$Nuclear, drop = TRUE)),
    permutations = 999
  )

  gen0_stage_results[[stg]] <- list(permanova = perm_sub, betadisper = bd_sub)

  cat("\nGen0 host composition PERMANOVA within", stg, ": Bray-Curtis ~ Mito * Nuclear\n")
  print(perm_sub)
  cat("\nGen0 dispersion check within", stg, "by Mito x Nuclear group\n")
  print(bd_sub)
}

safe_tab_to_csv(perm_gen0_all, file.path(stat_dir, "Figure4_final_Gen0_hostComposition_PERMANOVA_all.csv"))
safe_tab_to_csv(bd_gen0_all$tab, file.path(stat_dir, "Figure4_final_Gen0_hostComposition_betadisper_all.csv"))

for (stg in names(gen0_stage_results)) {
  safe_tab_to_csv(gen0_stage_results[[stg]]$permanova, file.path(stat_dir, paste0("Figure4_final_Gen0_hostComposition_PERMANOVA_", stg, ".csv")))
  safe_tab_to_csv(gen0_stage_results[[stg]]$betadisper$tab, file.path(stat_dir, paste0("Figure4_final_Gen0_hostComposition_betadisper_", stg, ".csv")))
}

gen0_groups <- meta.bc %>%
  filter(generation == "Generation_0") %>%
  group_by(Group) %>%
  summarise(n_gen0_samples = n(), n_gen0_vials = n_distinct(vial), .groups = "drop") %>%
  filter(n_gen0_vials >= MIN_GEN0_VIALS)

write_csv(gen0_groups, file.path(data_dir, "Figure4_final_Gen0_centroid_groups.csv"))

centroids <- lapply(gen0_groups$Group, function(g) {
  ids <- meta.bc$SampleID[meta.bc$Group == g & meta.bc$generation == "Generation_0"]
  cvec <- colMeans(otu.bc[ids, , drop = FALSE])
  list(Group = g, centroid = cvec)
})

centroids_df <- data.frame(Group = gen0_groups$Group, stringsAsFactors = FALSE)
centroids_df$centroid <- lapply(centroids, function(x) x$centroid)

gen5_ids <- meta.bc %>% filter(generation == "Generation_5") %>% pull(SampleID)

results_list <- list()

for (samp in gen5_ids) {
  rowm <- meta.bc[meta.bc$SampleID == samp, , drop = FALSE]
  g <- as.character(rowm$Group)
  if (!(g %in% centroids_df$Group)) next
  cen <- centroids_df$centroid[[which(centroids_df$Group == g)]]
  d <- as.numeric(vegdist(rbind(cen, otu.bc[samp, , drop = FALSE]), method = "bray")[1])

  results_list[[length(results_list) + 1]] <- data.frame(
    SampleID = samp, genotype = as.character(rowm$genotype),
    Mito = as.character(rowm$Mito), Nuclear = as.character(rowm$Nuclear),
    Stage = as.character(rowm$Stage), vial = as.character(rowm$vial),
    vial_id = as.character(rowm$vial_id), Group = g, Distance = d,
    stringsAsFactors = FALSE
  )
}

results <- bind_rows(results_list) %>%
  mutate(
    Mito = factor(Mito, levels = c("(ore)", "(simw501)")),
    Nuclear = factor(Nuclear, levels = c("OreR", "Aut")),
    Stage = factor(Stage, levels = c("Larvae", "Female")),
    vial = factor(vial),
    vial_id = factor(vial_id)
  )

cell_sizes <- results %>% count(Stage, Nuclear, Mito) %>% arrange(Stage, Nuclear, Mito)
cell_medians <- results %>%
  group_by(Stage, Nuclear, Mito) %>%
  summarise(median_distance = median(Distance, na.rm = TRUE), mean_distance = mean(Distance, na.rm = TRUE),
            sd_distance = sd(Distance, na.rm = TRUE), .groups = "drop") %>%
  arrange(Stage, Nuclear, Mito)

write_csv(results, file.path(data_dir, "Figure4_final_per_sample_distances.csv"))
write_csv(cell_sizes, file.path(data_dir, "Figure4_final_cell_sizes.csv"))
write_csv(cell_medians, file.path(data_dir, "Figure4_final_cell_medians_means.csv"))

mod_lme <- nlme::lme(
  Distance ~ Mito * Nuclear * Stage,
  random = ~ 1 | vial_id,
  data = results,
  method = "REML",
  na.action = na.omit,
  control = nlme::lmeControl(msMaxIter = 200, opt = "optim")
)

cat("\n================ 6) MIXED MODEL SUMMARY ================\n")
print(summary(mod_lme))

lme_anova <- anova(mod_lme)
print(lme_anova)

lme_anova_df <- as.data.frame(lme_anova) %>% rownames_to_column("term")
write_csv(lme_anova_df, file.path(stat_dir, "Figure4_final_LME_anova.csv"))

# BH correction below treats the four stage-by-nuclear contrasts as one family.
# Confidence intervals remain unadjusted, pointwise intervals.
emm_lme <- emmeans(mod_lme, ~ Mito | Nuclear * Stage)

contr_obj <- contrast(emm_lme, method = list("simw501 - ore" = c(-1, 1)))

contr_df <- as.data.frame(
  summary(contr_obj, infer = c(TRUE, TRUE), adjust = "none")
) %>%
  arrange(Stage, Nuclear) %>%
  mutate(
    p.raw = p.value,
    p.value = p.adjust(p.raw, method = "BH"),
    adjustment_family = "Four mitochondrial contrasts from the full mixed model",
    p.signif = p_to_star(p.value),
    interpretation = case_when(
      estimate > 0 ~ "simw501 greater divergence than ore",
      estimate < 0 ~ "simw501 lower divergence than ore",
      TRUE ~ "no difference"
    )
  )

write_csv(contr_df, file.path(stat_dir, "Figure4_final_emmeans_simw501_minus_ore.csv"))

run_focus_test <- function(df, stage_val, nuclear_val) {
  sub_df <- df %>% filter(Stage == stage_val, Nuclear == nuclear_val)

  mod <- nlme::lme(
    Distance ~ Mito,
    random = ~ 1 | vial_id,
    data = sub_df,
    method = "REML",
    control = nlme::lmeControl(msMaxIter = 200, opt = "optim")
  )

  emm <- emmeans(mod, ~ Mito)

  cont <- as.data.frame(
    summary(contrast(emm, method = list("simw501 - ore" = c(-1, 1))), infer = c(TRUE, TRUE), adjust = "none")
  ) %>%
    mutate(
      Stage = stage_val, Nuclear = nuclear_val,
      p.signif = p_to_star(p.value),
      interpretation = case_when(
        estimate > 0 ~ "simw501 greater divergence than ore",
        estimate < 0 ~ "simw501 lower divergence than ore",
        TRUE ~ "no difference"
      )
    )

  return(cont)
}

res_larvae_orer <- run_focus_test(results, "Larvae", "OreR")
res_larvae_aut  <- run_focus_test(results, "Larvae", "Aut")
res_female_orer <- run_focus_test(results, "Female", "OreR")
res_female_aut  <- run_focus_test(results, "Female", "Aut")

focused_results <- bind_rows(res_larvae_orer, res_larvae_aut, res_female_orer, res_female_aut) %>%
  arrange(Stage, Nuclear) %>%
  mutate(
    p.raw = p.value,
    p.value = p.adjust(p.raw, method = "BH"),
    adjustment_family = "Four mitochondrial contrasts from separately fitted models",
    p.signif = p_to_star(p.value)
  )

write_csv(focused_results, file.path(stat_dir, "Figure4_final_focused_tests_simw501_minus_ore.csv"))

# -----------------------
# Median-centered variance comparison for the derived distance response.
# Pools share vials; this is an exploratory check, not a residual diagnostic.
# -----------------------

cat("\nLevene's test - Distance variance by Mito x Nuclear x Stage\n")
levene_distance <- car::leveneTest(
  Distance ~ Mito * Nuclear * Stage, data = results, center = median
)
print(levene_distance)

levene_distance_df <- as.data.frame(levene_distance) %>%
  rownames_to_column("term") %>%
  mutate(
    Analysis = "Levene's test",
    Response = "Gen5-to-Gen0 divergence Distance",
    Model = "Distance ~ Mito * Nuclear * Stage"
  )

write_csv(levene_distance_df, file.path(stat_dir, "Figure4_final_Levene_Distance_dispersion.csv"))

# ------------------------------
# B8) Main Figure 4 plot
# ------------------------------

plot_results <- results %>%
  mutate(
    Stage_label = dplyr::recode(as.character(Stage), "Larvae" = "Larvae", "Female" = "Adult females"),
    Stage_label = factor(Stage_label, levels = c("Larvae", "Adult females")),
    Nuclear = factor(Nuclear, levels = c("OreR", "Aut")),
    Mito = factor(Mito, levels = c("(ore)", "(simw501)"))
  )

mito_fill_cols <- c("(ore)" = "#B8C6A6", "(simw501)" = "#D9B382")
mito_line_cols <- c("(ore)" = "#3F5F3A", "(simw501)" = "#8A5A00")

p_final <- ggplot(
  plot_results,
  aes(x = Nuclear, y = Distance, fill = Mito, color = Mito)
) +
  geom_boxplot(position = position_dodge(width = DODGE_W), width = 0.56, outlier.shape = NA, linewidth = 0.65, alpha = 0.85) +
  geom_point(aes(shape = Mito), position = position_jitterdodge(jitter.width = 0.10, dodge.width = DODGE_W),
             size = 2.35, alpha = 0.88, stroke = 0.55) +
  stat_summary(fun = median, geom = "point", position = position_dodge(width = DODGE_W),
               shape = 23, size = 3.0, fill = "white", color = "black", stroke = 0.65, show.legend = FALSE) +
  facet_wrap(~ Stage_label, nrow = 1) +
  scale_fill_manual(values = mito_fill_cols, labels = c("(ore)" = "(ore)", "(simw501)" = expression("(simw"^501*")")), name = "Mitochondrial genotype") +
  scale_color_manual(values = mito_line_cols, labels = c("(ore)" = "(ore)", "(simw501)" = expression("(simw"^501*")")), name = "Mitochondrial genotype") +
  scale_shape_manual(values = c("(ore)" = 21, "(simw501)" = 24), labels = c("(ore)" = "(ore)", "(simw501)" = expression("(simw"^501*")")), name = "Mitochondrial genotype") +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1.0, by = 0.2), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Nuclear genotype", y = "Bray-Curtis distance to Generation 0 centroid") +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 24), color = mito_line_cols, fill = mito_fill_cols, size = 3)), color = "none", shape = "none") +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold", color = "black", size = 12),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_text(face = "bold"),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.45),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.65),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.65),
    strip.text = element_text(face = "bold", size = 12.5, color = "black"),
    panel.spacing = unit(1.2, "lines"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11.5),
    legend.text = element_text(size = 10.5),
    legend.key.size = unit(0.65, "cm"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 10, r = 18, b = 10, l = 10)
  )

p_final

ggsave(file.path(fig_dir, "Figure4_final_Gen5_to_Gen0_MitoNuclear_Divergence_ISME.pdf"), plot = p_final,
       width = 8.2, height = 4.6, units = "in", device = PDF_DEVICE, bg = "white")
ggsave(file.path(fig_dir, "Figure4_final_Gen5_to_Gen0_MitoNuclear_Divergence_ISME.png"), plot = p_final,
       width = 8.2, height = 4.6, units = "in", dpi = 600, bg = "white")
ggsave(file.path(fig_dir, "Figure4_final_Gen5_to_Gen0_MitoNuclear_Divergence_ISME.jpeg"), plot = p_final,
       width = 8.2, height = 4.6, units = "in", dpi = 600, bg = "white", device = "jpeg")

# ============================================================
# PART C: SUPPLEMENTARY FIGURE S5
# Vial-level sensitivity analysis for Figure 4
# ============================================================

ps.s5 <- subset_samples(
  ps.clean,
  !is.na(generation) &
    generation %in% c("Generation_0", "Generation_5") &
    life_stage %in% c("Larvae", "Female") &
    genotype %in% c("oreore", "simore", "oreaut", "simaut")
)

ps.s5 <- prune_samples(sample_sums(ps.s5) > 0, ps.s5)
ps.s5 <- prune_taxa(taxa_sums(ps.s5) > 0, ps.s5)

ps.s5.ra <- transform_sample_counts(ps.s5, function(x) if (sum(x) > 0) x / sum(x) else x)
ps.s5.ra <- prune_taxa(taxa_sums(ps.s5.ra) > 0, ps.s5.ra)

otu.s5 <- as(otu_table(ps.s5.ra), "matrix")
if (taxa_are_rows(ps.s5.ra)) otu.s5 <- t(otu.s5)

meta.s5 <- data.frame(sample_data(ps.s5.ra)) %>%
  rownames_to_column("SampleID") %>%
  mutate(
    genotype = as.character(genotype), generation = as.character(generation),
    life_stage = as.character(life_stage), vial = as.character(vial),
    Mito = case_when(genotype %in% c("oreore", "oreaut") ~ "(ore)", genotype %in% c("simore", "simaut") ~ "(simw501)", TRUE ~ NA_character_),
    Mito = factor(Mito, levels = c("(ore)", "(simw501)")),
    Nuclear = case_when(genotype %in% c("oreore", "simore") ~ "OreR", genotype %in% c("oreaut", "simaut") ~ "Aut", TRUE ~ NA_character_),
    Nuclear = factor(Nuclear, levels = c("OreR", "Aut")),
    Stage = dplyr::recode(life_stage, "Larvae" = "Larvae", "Female" = "Adult females"),
    Stage = factor(Stage, levels = c("Larvae", "Adult females")),
    Generation = dplyr::recode(generation, "Generation_0" = "Gen0", "Generation_5" = "Gen5"),
    Generation = factor(Generation, levels = c("Gen0", "Gen5")),
    Group = paste(genotype, life_stage, sep = "_"),
    VialID = paste(genotype, generation, life_stage, vial, sep = "__")
  )

otu.s5 <- otu.s5[meta.s5$SampleID, , drop = FALSE]

vial_meta <- meta.s5 %>%
  group_by(VialID, genotype, generation, Mito, Nuclear, Stage, Generation, Group, vial) %>%
  summarise(n_pools_collapsed = n(), .groups = "drop") %>%
  arrange(genotype, Stage, Generation, vial)

otu_vial <- rowsum(otu.s5, group = meta.s5$VialID, reorder = FALSE)
n_per_vial <- as.vector(table(factor(meta.s5$VialID, levels = rownames(otu_vial))))
otu_vial <- sweep(otu_vial, 1, n_per_vial, "/")
otu_vial <- sweep(otu_vial, 1, rowSums(otu_vial), "/")
otu_vial[is.na(otu_vial)] <- 0
otu_vial <- otu_vial[vial_meta$VialID, , drop = FALSE]

write_csv(vial_meta, file.path(data_dir, "Supplementary_Figure_S5_vial_level_metadata.csv"))

gen0_groups_s5 <- vial_meta %>%
  filter(Generation == "Gen0") %>%
  group_by(Group) %>%
  summarise(n_gen0_vials = n_distinct(vial), .groups = "drop")

centroids_s5 <- lapply(gen0_groups_s5$Group, function(g) {
  ids <- vial_meta$VialID[vial_meta$Group == g & vial_meta$Generation == "Gen0"]
  cvec <- colMeans(otu_vial[ids, , drop = FALSE])
  list(Group = g, centroid = cvec)
})

centroids_df_s5 <- data.frame(Group = gen0_groups_s5$Group, stringsAsFactors = FALSE)
centroids_df_s5$centroid <- lapply(centroids_s5, function(x) x$centroid)

gen5_vials <- vial_meta %>% filter(Generation == "Gen5") %>% pull(VialID)

results_list_s5 <- list()

for (v in gen5_vials) {
  rowm <- vial_meta[vial_meta$VialID == v, , drop = FALSE]
  g <- as.character(rowm$Group[1])
  if (!(g %in% centroids_df_s5$Group)) next
  cen <- centroids_df_s5$centroid[[which(centroids_df_s5$Group == g)]]
  d <- as.numeric(vegan::vegdist(rbind(cen, otu_vial[v, , drop = FALSE]), method = "bray")[1])

  results_list_s5[[length(results_list_s5) + 1]] <- data.frame(
    VialID = v, genotype = as.character(rowm$genotype[1]),
    Mito = as.character(rowm$Mito[1]), Nuclear = as.character(rowm$Nuclear[1]),
    Stage = as.character(rowm$Stage[1]), vial = as.character(rowm$vial[1]),
    Group = g, Distance = d, n_pools_collapsed = rowm$n_pools_collapsed[1],
    stringsAsFactors = FALSE
  )
}

s5_dist <- bind_rows(results_list_s5) %>%
  mutate(
    Mito = factor(Mito, levels = c("(ore)", "(simw501)")),
    Nuclear = factor(Nuclear, levels = c("OreR", "Aut")),
    Stage = factor(Stage, levels = c("Larvae", "Adult females")),
    vial = factor(vial)
  )

write_csv(s5_dist, file.path(data_dir, "Supplementary_Figure_S5_vial_level_distances.csv"))

cell_summary_s5 <- s5_dist %>%
  group_by(Stage, Nuclear, Mito) %>%
  summarise(n_vials = n(), mean_distance = mean(Distance, na.rm = TRUE), median_distance = median(Distance, na.rm = TRUE),
            sd_distance = sd(Distance, na.rm = TRUE), se_distance = sd_distance / sqrt(n_vials), .groups = "drop") %>%
  arrange(Stage, Nuclear, Mito)

write_csv(cell_summary_s5, file.path(data_dir, "Supplementary_Figure_S5_cell_summary.csv"))

mod_s5 <- lm(Distance ~ Mito * Nuclear * Stage, data = s5_dist)
anova_s5 <- anova(mod_s5)

write_csv(as.data.frame(anova_s5) %>% rownames_to_column("term"), file.path(stat_dir, "Supplementary_Figure_S5_model_anova.csv"))

emm_s5 <- emmeans(mod_s5, ~ Mito | Nuclear * Stage)
contr_s5 <- contrast(emm_s5, method = list("simw501 - ore" = c(-1, 1)))

contr_s5_df <- as.data.frame(summary(contr_s5, infer = c(TRUE, TRUE), adjust = "none")) %>%
  arrange(Stage, Nuclear) %>%
  mutate(
    p.raw = p.value,
    p.value = p.adjust(p.raw, method = "BH"),
    adjustment_family = "Four mitochondrial contrasts from the vial-mean sensitivity model",
    interpretation = case_when(estimate > 0 ~ "simw501 greater divergence than ore", estimate < 0 ~ "simw501 lower divergence than ore", TRUE ~ "no difference")
  )

write_csv(contr_s5_df, file.path(stat_dir, "Supplementary_Figure_S5_emmeans_contrasts.csv"))

plot_df_s5 <- s5_dist %>%
  mutate(Stage = factor(Stage, levels = c("Larvae", "Adult females")), Nuclear = factor(Nuclear, levels = c("OreR", "Aut")), Mito = factor(Mito, levels = c("(ore)", "(simw501)")))

mito_fill_cols_s5 <- c("(ore)" = "#B8C6A6", "(simw501)" = "#D9B382")
mito_line_cols_s5 <- c("(ore)" = "#3F5F3A", "(simw501)" = "#8A5A00")

S5_vial_sensitivity <- ggplot(plot_df_s5, aes(x = Nuclear, y = Distance, fill = Mito, color = Mito)) +
  geom_point(aes(shape = Mito), position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.65), size = 3.1, alpha = 0.95, stroke = 0.65) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.65), shape = 23, size = 3.4, fill = "white", color = "black", stroke = 0.70, show.legend = FALSE) +
  stat_summary(fun.data = mean_se, geom = "errorbar", position = position_dodge(width = 0.65), width = 0.15, linewidth = 0.65, color = "black") +
  facet_wrap(~ Stage, nrow = 1) +
  scale_fill_manual(values = mito_fill_cols_s5, breaks = c("(ore)", "(simw501)"), labels = c("(ore)" = "(ore)", "(simw501)" = expression(simw^501)), name = "Mitochondrial genotype") +
  scale_color_manual(values = mito_line_cols_s5, breaks = c("(ore)", "(simw501)"), labels = c("(ore)" = "(ore)", "(simw501)" = expression(simw^501)), name = "Mitochondrial genotype") +
  scale_shape_manual(values = c("(ore)" = 21, "(simw501)" = 24), breaks = c("(ore)", "(simw501)"), labels = c("(ore)" = "(ore)", "(simw501)" = expression(simw^501)), name = "Mitochondrial genotype") +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1.0, by = 0.2), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Nuclear genotype", y = "Bray-Curtis distance to Generation 0 centroid") +
  guides(fill = guide_legend(override.aes = list(shape = c(21, 24), color = mito_line_cols_s5, fill = mito_fill_cols_s5, size = 3)), color = "none", shape = "none") +
  theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(face = "bold", color = "black", size = 12),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_text(face = "bold"),
    axis.line = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
    strip.text = element_text(face = "bold", size = 12.5, color = "black"),
    panel.spacing = unit(1.0, "lines"),
    legend.position = "top",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10.5, color = "black"),
    legend.key.size = unit(0.55, "cm"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 10, r = 12, b = 10, l = 12)
  )

S5_vial_sensitivity

ggsave(file.path(fig_dir, "Supplementary_Figure_S5_vial_level_Fig4_sensitivity.pdf"), plot = S5_vial_sensitivity,
       width = 8.8, height = 4.9, units = "in", device = PDF_DEVICE, bg = "white")
ggsave(file.path(fig_dir, "Supplementary_Figure_S5_vial_level_Fig4_sensitivity.png"), plot = S5_vial_sensitivity,
       width = 8.8, height = 4.9, units = "in", dpi = 600, bg = "white")
ggsave(file.path(fig_dir, "Supplementary_Figure_S5_vial_level_Fig4_sensitivity.jpeg"), plot = S5_vial_sensitivity,
       width = 8.8, height = 4.9, units = "in", dpi = 600, bg = "white", device = "jpeg")

# ============================================================
# PART D + E: STATS SUMMARY + SUPPLEMENTARY TABLE S3
# ============================================================

capture.output(
  {
    cat("FIGURE 4 FINAL STATS SUMMARY\n")
    cat("============================\n\n")

    cat("Gen0 host composition PERMANOVA: Bray-Curtis ~ Mito * Nuclear * Stage\n")
    print(perm_gen0_all)
    cat("\nGen0 dispersion check by Mito x Nuclear x Stage group\n")
    print(bd_gen0_all)

    for (stg in names(gen0_stage_results)) {
      cat("\nGen0 host composition PERMANOVA within", stg, ": Bray-Curtis ~ Mito * Nuclear\n")
      print(gen0_stage_results[[stg]]$permanova)
      cat("\nGen0 dispersion check within", stg, "by Mito x Nuclear group\n")
      print(gen0_stage_results[[stg]]$betadisper)
    }

    cat("\nMixed model summary:\n")
    print(summary(mod_lme))
    cat("\nMixed model ANOVA:\n")
    print(lme_anova)
    cat("\nPrimary emmeans contrasts: simw501 - ore\n")
    print(contr_df)
    cat("\nFocused tests: simw501 - ore\n")
    print(focused_results)

    cat("\nLevene's test - Distance variance by Mito x Nuclear x Stage\n")
    print(levene_distance)
  },
  file = file.path(stat_dir, "Figure4_final_stats_summary.txt")
)

# Add new Levene sheet to Supplementary Table S3

supp_s3_levene <- levene_distance_df

wb <- createWorkbook()

addWorksheet(wb, "Gen0_PERMANOVA_all")
writeData(wb, "Gen0_PERMANOVA_all", as.data.frame(perm_gen0_all) %>% rownames_to_column("term"))

addWorksheet(wb, "Gen0_betadisper_all")
writeData(wb, "Gen0_betadisper_all", as.data.frame(bd_gen0_all$tab) %>% rownames_to_column("term"))

# Include stage-specific tests in the exported workbook.

addWorksheet(wb, "Gen0_PERMANOVA_Larvae")
writeData(wb, "Gen0_PERMANOVA_Larvae", as.data.frame(gen0_stage_results[["Larvae"]]$permanova) %>% rownames_to_column("term"))

addWorksheet(wb, "Gen0_betadisper_Larvae")
writeData(wb, "Gen0_betadisper_Larvae", as.data.frame(gen0_stage_results[["Larvae"]]$betadisper$tab) %>% rownames_to_column("term"))

addWorksheet(wb, "Gen0_PERMANOVA_Female")
writeData(wb, "Gen0_PERMANOVA_Female", as.data.frame(gen0_stage_results[["Female"]]$permanova) %>% rownames_to_column("term"))

addWorksheet(wb, "Gen0_betadisper_Female")
writeData(wb, "Gen0_betadisper_Female", as.data.frame(gen0_stage_results[["Female"]]$betadisper$tab) %>% rownames_to_column("term"))

addWorksheet(wb, "Fig4_LME_anova")
writeData(wb, "Fig4_LME_anova", lme_anova_df)

addWorksheet(wb, "Fig4_focused_contrasts")
writeData(wb, "Fig4_focused_contrasts", focused_results)

addWorksheet(wb, "Fig4_Levene_Distance")
writeData(wb, "Fig4_Levene_Distance", supp_s3_levene)

addWorksheet(wb, "S5_vial_sensitivity_ANOVA")
writeData(wb, "S5_vial_sensitivity_ANOVA", as.data.frame(anova_s5) %>% rownames_to_column("term"))

addWorksheet(wb, "S5_vial_sensitivity_contrasts")
writeData(wb, "S5_vial_sensitivity_contrasts", contr_s5_df)

header_style <- createStyle(textDecoration = "bold", halign = "center")
for (sheet in names(wb)) {
  addStyle(wb, sheet = sheet, style = header_style, rows = 1, cols = 1:50, gridExpand = TRUE)
  setColWidths(wb, sheet = sheet, cols = 1:50, widths = "auto")
}

saveWorkbook(wb, file = file.path(stat_dir, "Supplementary_Table_S3_Figure4_mitonuclear_divergence_statistics.xlsx"), overwrite = TRUE)

write_csv(supp_s3_levene, file.path(stat_dir, "Supplementary_Table_S3I_Fig4_Levene_Distance.csv"))

write_csv(as.data.frame(gen0_stage_results[["Larvae"]]$permanova) %>% rownames_to_column("term"),
          file.path(stat_dir, "Supplementary_Table_S3J_Gen0_PERMANOVA_Larvae.csv"))
write_csv(as.data.frame(gen0_stage_results[["Larvae"]]$betadisper$tab) %>% rownames_to_column("term"),
          file.path(stat_dir, "Supplementary_Table_S3K_Gen0_betadisper_Larvae.csv"))
write_csv(as.data.frame(gen0_stage_results[["Female"]]$permanova) %>% rownames_to_column("term"),
          file.path(stat_dir, "Supplementary_Table_S3L_Gen0_PERMANOVA_Female.csv"))
write_csv(as.data.frame(gen0_stage_results[["Female"]]$betadisper$tab) %>% rownames_to_column("term"),
          file.path(stat_dir, "Supplementary_Table_S3M_Gen0_betadisper_Female.csv"))

cat("\nDONE. Genotype divergence outputs saved:\n")
cat(file.path(stat_dir, "Figure4_final_Levene_Distance_dispersion.csv"), "\n")
cat(file.path(stat_dir, "Supplementary_Table_S3I_Fig4_Levene_Distance.csv"), "\n")

# Record software versions used for this run.
capture.output(sessionInfo(), file = file.path(main_outdir, "sessionInfo.txt"))
