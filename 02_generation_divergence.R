# Temporal community divergence
# Compare host microbiome composition across generations in oreore.
# Average pooled samples within each vial before community comparisons.
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

required_packages <- c("phyloseq", "dplyr", "tidyr", "tibble", "vegan", "ggplot2", "readr", "patchwork", "openxlsx")
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
  library(tibble)
  library(vegan)
  library(ggplot2)
  library(readr)
  library(patchwork)
  library(openxlsx)
})

# Seed for reproducible permutation tests.
set.seed(1234)

# -----------------------
# 2) Create output folders
# -----------------------

main_outdir <- file.path(OUTPUT_ROOT, "Fig2_oreOreR_generation_divergence")

dir.create(main_outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(main_outdir, "input_processed"), showWarnings = FALSE)
dir.create(file.path(main_outdir, "figure_files"), showWarnings = FALSE)
dir.create(file.path(main_outdir, "figure_source_data"), showWarnings = FALSE)
dir.create(file.path(main_outdir, "statistics"), showWarnings = FALSE)

fig_dir   <- file.path(main_outdir, "figure_files")
data_dir  <- file.path(main_outdir, "figure_source_data")
stat_dir  <- file.path(main_outdir, "statistics")
input_dir <- file.path(main_outdir, "input_processed")

# -----------------------
# 3) Load phyloseq object
# -----------------------

ps.clean <- readRDS(INPUT_RDS)

if (!inherits(ps.clean, "phyloseq")) {
  stop("INPUT_RDS must contain a phyloseq object.", call. = FALSE)
}
required_metadata <- c("genotype", "generation", "life_stage", "vial", "replicate")
missing_metadata <- setdiff(required_metadata, names(data.frame(sample_data(ps.clean))))
if (length(missing_metadata) > 0) {
  stop("Missing sample metadata: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
}

saveRDS(
  ps.clean,
  file = file.path(input_dir, "06_phyloseq_clean_CHAP1_NOHOST.rds")
)

# -----------------------
# 4) Subset to compatible control genotype only
# -----------------------

ps.oo <- subset_samples(
  ps.clean,
  genotype == "oreore" &
    !is.na(generation) &
    life_stage %in% c("Larvae", "Female")
)

ps.oo <- prune_samples(sample_sums(ps.oo) > 0, ps.oo)

# -----------------------
# 5) Transform to relative abundance
# -----------------------

ps.ra <- transform_sample_counts(ps.oo, function(x) x / sum(x))
ps.ra <- prune_samples(sample_sums(ps.ra) > 0, ps.ra)

otu <- as(otu_table(ps.ra), "matrix")
if (taxa_are_rows(ps.ra)) otu <- t(otu)

meta <- data.frame(sample_data(ps.ra)) %>%
  rownames_to_column("phyloseq_sample_id") %>%
  mutate(
    generation = factor(generation, levels = paste0("Generation_", 0:5)),
    life_stage = factor(life_stage, levels = c("Larvae", "Female")),
    vial = as.character(vial),
    replicate = as.character(replicate),
    genotype = as.character(genotype)
  )

otu <- otu[meta$phyloseq_sample_id, , drop = FALSE]

# -----------------------
# 6) Collapse to vial-level means
# -----------------------

meta <- meta %>%
  mutate(
    vial_group = paste(genotype, generation, life_stage, vial, sep = "__")
  )

otu_vial <- rowsum(otu, group = meta$vial_group, reorder = FALSE)

n_per_group <- as.vector(
  table(factor(meta$vial_group, levels = rownames(otu_vial)))
)

otu_vial <- sweep(otu_vial, 1, n_per_group, "/")
otu_vial <- sweep(otu_vial, 1, rowSums(otu_vial), "/")
otu_vial[is.na(otu_vial)] <- 0

meta_vial <- meta %>%
  group_by(vial_group) %>%
  summarise(
    genotype = first(genotype),
    generation = first(generation),
    life_stage = first(life_stage),
    vial = first(vial),
    n_pools_collapsed = n(),
    .groups = "drop"
  ) %>%
  as.data.frame()

rownames(meta_vial) <- meta_vial$vial_group

otu_use <- otu_vial[meta_vial$vial_group, , drop = FALSE]

meta_use <- meta_vial %>%
  rownames_to_column("vial_sample_id") %>%
  mutate(
    sample_id = vial_sample_id,
    generation = factor(generation, levels = paste0("Generation_", 0:5)),
    life_stage = factor(life_stage, levels = c("Larvae", "Female")),
    life_stage_label = recode(
      as.character(life_stage),
      "Larvae" = "Larvae",
      "Female" = "Adult females"
    ),
    life_stage_label = factor(
      life_stage_label,
      levels = c("Larvae", "Adult females")
    )
  )

otu_use <- otu_use[meta_use$sample_id, , drop = FALSE]

cat("\n============================\n")
cat("SAMPLE COUNTS USED\n")
cat("============================\n")
print(table(meta_use$generation, meta_use$life_stage))

# ============================================================
# PANEL A: DISTANCE TO GENERATION 0 CENTROID
#
# ============================================================

centroids_G0 <- meta_use %>%
  filter(generation == "Generation_0") %>%
  group_by(life_stage) %>%
  summarise(gen0_ids = list(sample_id), .groups = "drop") %>%
  rowwise() %>%
  mutate(
    centroid = list(colMeans(otu_use[unlist(gen0_ids), , drop = FALSE]))
  ) %>%
  ungroup()

dist_G0_df_all <- meta_use %>%
  select(sample_id, generation, life_stage, life_stage_label, vial, n_pools_collapsed) %>%
  rowwise() %>%
  mutate(
    DistanceToGen0 = {
      cen <- centroids_G0$centroid[[which(centroids_G0$life_stage == life_stage)]]
      as.numeric(
        vegdist(
          rbind(cen, otu_use[sample_id, , drop = FALSE]),
          method = "bray"
        )[1]
      )
    }
  ) %>%
  ungroup()

dist_G0_df <- dist_G0_df_all %>%
  filter(generation %in% paste0("Generation_", 1:5)) %>%
  mutate(
    generation_label = factor(
      generation,
      levels = paste0("Generation_", 1:5),
      labels = paste0("Gen ", 1:5)
    )
  )

write_csv(
  dist_G0_df,
  file.path(data_dir, "Figure2_final_A_distance_to_G0_data.csv")
)

summary_G0_df <- dist_G0_df %>%
  group_by(generation, generation_label, life_stage_label) %>%
  summarise(
    median_dist = median(DistanceToGen0, na.rm = TRUE),
    mean_dist = mean(DistanceToGen0, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# PANEL B: DISTANCE TO GENERATION 1 CENTROID
#
# ============================================================

centroids_G1 <- meta_use %>%
  filter(generation == "Generation_1") %>%
  group_by(life_stage) %>%
  summarise(gen1_ids = list(sample_id), .groups = "drop") %>%
  rowwise() %>%
  mutate(
    centroid = list(colMeans(otu_use[unlist(gen1_ids), , drop = FALSE]))
  ) %>%
  ungroup()

dist_G1_df <- meta_use %>%
  filter(generation %in% paste0("Generation_", 2:5)) %>%
  select(sample_id, generation, life_stage, life_stage_label, vial, n_pools_collapsed) %>%
  rowwise() %>%
  mutate(
    DistanceToGen1 = {
      cen <- centroids_G1$centroid[[which(centroids_G1$life_stage == life_stage)]]
      as.numeric(
        vegdist(
          rbind(cen, otu_use[sample_id, , drop = FALSE]),
          method = "bray"
        )[1]
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    generation_label = factor(
      generation,
      levels = paste0("Generation_", 2:5),
      labels = paste0("Gen ", 2:5)
    )
  )

write_csv(
  dist_G1_df,
  file.path(data_dir, "Figure2_final_B_distance_to_G1_data.csv")
)

summary_G1_df <- dist_G1_df %>%
  group_by(generation, generation_label, life_stage_label) %>%
  summarise(
    median_dist = median(DistanceToGen1, na.rm = TRUE),
    mean_dist = mean(DistanceToGen1, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# PLOTTING
# ============================================================

stage_cols <- c(
  "Larvae" = "#D55E00",
  "Adult females" = "#0072B2"
)

shared_y_max <- max(
  dist_G0_df$DistanceToGen0,
  dist_G1_df$DistanceToGen1,
  na.rm = TRUE
) * 1.10

theme_fig2_boxed <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    axis.title = element_text(face = "bold", size = 11.5),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.text = element_text(color = "black", size = 10),
    axis.text.x = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.25),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.65),
    axis.line = element_blank(),
    plot.tag = element_text(face = "bold", size = 20),
    plot.tag.position = c(-0.035, 1.03),
    plot.margin = margin(t = 22, r = 18, b = 16, l = 42)
  )

pA <- ggplot(
  dist_G0_df,
  aes(x = generation_label, y = DistanceToGen0, fill = life_stage_label, color = life_stage_label)
) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.58, alpha = 0.22,
               outlier.shape = NA, linewidth = 0.65) +
  geom_point(aes(shape = life_stage_label),
             position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
             size = 2.6, alpha = 0.95, stroke = 0.75) +
  geom_line(data = summary_G0_df,
            aes(x = generation_label, y = median_dist, group = life_stage_label, color = life_stage_label),
            position = position_dodge(width = 0.75), linewidth = 0.9, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = stage_cols, name = "Life stage") +
  scale_color_manual(values = stage_cols, name = "Life stage") +
  scale_shape_manual(values = c("Larvae" = 17, "Adult females" = 16), name = "Life stage") +
  labs(x = "Generation", y = "Bray-Curtis distance to Generation 0 centroid",
       title = "Divergence from Generation 0", tag = "A") +
  coord_cartesian(ylim = c(0, shared_y_max), clip = "off") +
  theme_fig2_boxed

pB <- ggplot(
  dist_G1_df,
  aes(x = generation_label, y = DistanceToGen1, fill = life_stage_label, color = life_stage_label)
) +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.58, alpha = 0.22,
               outlier.shape = NA, linewidth = 0.65) +
  geom_point(aes(shape = life_stage_label),
             position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
             size = 2.6, alpha = 0.95, stroke = 0.75) +
  geom_line(data = summary_G1_df,
            aes(x = generation_label, y = median_dist, group = life_stage_label, color = life_stage_label),
            position = position_dodge(width = 0.75), linewidth = 0.9, linetype = "dashed", inherit.aes = FALSE) +
  scale_fill_manual(values = stage_cols, name = "Life stage") +
  scale_color_manual(values = stage_cols, name = "Life stage") +
  scale_shape_manual(values = c("Larvae" = 17, "Adult females" = 16), name = "Life stage") +
  labs(x = "Generation", y = "Bray-Curtis distance to Generation 1 centroid",
       title = "Divergence from Generation 1", tag = "B") +
  coord_cartesian(ylim = c(0, shared_y_max), clip = "off") +
  theme_fig2_boxed +
  theme(legend.position = "none")

pA
pB

ggsave(file.path(fig_dir, "Figure2_final_A_distance_to_G0_Gen1_to_Gen5_boxed.pdf"), plot = pA,
       width = 8.8, height = 4.8, units = "in", device = PDF_DEVICE)
ggsave(file.path(fig_dir, "Figure2_final_A_distance_to_G0_Gen1_to_Gen5_boxed.png"), plot = pA,
       width = 8.8, height = 4.8, units = "in", dpi = 600, bg = "white")
ggsave(file.path(fig_dir, "Figure2_final_A_distance_to_G0_Gen1_to_Gen5_boxed.jpeg"), plot = pA,
       width = 8.8, height = 4.8, units = "in", dpi = 600, bg = "white", device = "jpeg")

ggsave(file.path(fig_dir, "Figure2_final_B_distance_to_G1_Gen2_to_Gen5_boxed.pdf"), plot = pB,
       width = 8.8, height = 4.8, units = "in", device = PDF_DEVICE)
ggsave(file.path(fig_dir, "Figure2_final_B_distance_to_G1_Gen2_to_Gen5_boxed.png"), plot = pB,
       width = 8.8, height = 4.8, units = "in", dpi = 600, bg = "white")
ggsave(file.path(fig_dir, "Figure2_final_B_distance_to_G1_Gen2_to_Gen5_boxed.jpeg"), plot = pB,
       width = 8.8, height = 4.8, units = "in", dpi = 600, bg = "white", device = "jpeg")

combined_fig2 <- pA / pB
combined_fig2

ggsave(file.path(fig_dir, "Figure2_final_combined_A_B_boxed.pdf"), plot = combined_fig2,
       width = 9.2, height = 9.4, units = "in", device = PDF_DEVICE)
ggsave(file.path(fig_dir, "Figure2_final_combined_A_B_boxed.png"), plot = combined_fig2,
       width = 9.2, height = 9.4, units = "in", dpi = 600, bg = "white")
ggsave(file.path(fig_dir, "Figure2_final_combined_A_B_boxed.jpeg"), plot = combined_fig2,
       width = 9.2, height = 9.4, units = "in", dpi = 600, bg = "white", device = "jpeg")

# ============================================================
# STATS
# ============================================================

# -----------------------
# PERMANOVA: generation * life_stage
# -----------------------

bray_all <- vegdist(otu_use, method = "bray")

meta_stats <- meta_use %>%
  mutate(
    generation = factor(generation, levels = paste0("Generation_", 0:5)),
    life_stage = factor(life_stage, levels = c("Larvae", "Female"))
  )

perm_time <- adonis2(
  bray_all ~ generation * life_stage,
  data = meta_stats,
  permutations = 999,
  by = "margin"
)

perm_time_df <- as.data.frame(perm_time) %>%
  rownames_to_column("term")

write_csv(
  perm_time_df,
  file.path(stat_dir, "Figure2_final_PERMANOVA_generation_lifestage.csv")
)

# -----------------------
# Betadisper by life stage
# -----------------------

bd_stage <- betadisper(bray_all, group = meta_stats$life_stage)
bd_stage_test <- permutest(bd_stage, permutations = 999)

bd_stage_df <- data.frame(
  sample_id = names(bd_stage$distances),
  distance_to_lifestage_centroid = bd_stage$distances
) %>%
  left_join(
    meta_stats %>% select(sample_id, generation, life_stage, vial),
    by = "sample_id"
  )

write_csv(
  bd_stage_df,
  file.path(stat_dir, "Figure2_final_betadisper_lifestage.csv")
)

# -----------------------
# Dispersion differences across generation groups.

cat("\nBetadisper by generation alone\n")
bd_generation <- betadisper(bray_all, group = meta_stats$generation)
bd_generation_test <- permutest(bd_generation, permutations = 999)
print(bd_generation_test)

bd_generation_df <- data.frame(
  sample_id = names(bd_generation$distances),
  distance_to_generation_centroid = bd_generation$distances
) %>%
  left_join(
    meta_stats %>% select(sample_id, generation, life_stage, vial),
    by = "sample_id"
  )

write_csv(
  bd_generation_df,
  file.path(stat_dir, "Figure2_final_betadisper_generation.csv")
)

# -----------------------
# Dispersion differences across combined generation and stage groups.
# This is an omnibus group comparison, not a factorial interaction test.

cat("\nBetadisper by combined generation x life_stage groups\n")
meta_stats <- meta_stats %>%
  mutate(gen_by_stage = interaction(generation, life_stage, drop = TRUE))

bd_gen_stage <- betadisper(bray_all, group = meta_stats$gen_by_stage)
bd_gen_stage_test <- permutest(bd_gen_stage, permutations = 999)
print(bd_gen_stage_test)

bd_gen_stage_df <- data.frame(
  sample_id = names(bd_gen_stage$distances),
  distance_to_gen_stage_centroid = bd_gen_stage$distances
) %>%
  left_join(
    meta_stats %>% select(sample_id, generation, life_stage, vial),
    by = "sample_id"
  )

write_csv(
  bd_gen_stage_df,
  file.path(stat_dir, "Figure2_final_betadisper_generation_x_lifestage.csv")
)

# -----------------------
# Linear model on Panel A response
# -----------------------

lm_G0 <- lm(DistanceToGen0 ~ generation * life_stage, data = dist_G0_df)
lm_G0_anova <- anova(lm_G0)
lm_G0_anova_df <- as.data.frame(lm_G0_anova) %>% rownames_to_column("term")
write_csv(lm_G0_anova_df, file.path(stat_dir, "Figure2_final_distance_to_G0_lm.csv"))

# -----------------------
# Linear model on Panel B response
# -----------------------

lm_G1 <- lm(DistanceToGen1 ~ generation * life_stage, data = dist_G1_df)
lm_G1_anova <- anova(lm_G1)
lm_G1_anova_df <- as.data.frame(lm_G1_anova) %>% rownames_to_column("term")
write_csv(lm_G1_anova_df, file.path(stat_dir, "Figure2_final_distance_to_G1_lm.csv"))

# -----------------------
# Direct Gen1 vs Gen5 PERMANOVA within each life stage
# -----------------------

extract_adonis_term <- function(adonis_obj) {
  df <- as.data.frame(adonis_obj)
  if ("Model" %in% rownames(df)) return(df["Model", , drop = FALSE])
  keep <- setdiff(rownames(df), c("Residual", "Total"))
  df[keep[1], , drop = FALSE]
}

extract_betadisper_group <- function(bd_obj) {
  df <- as.data.frame(bd_obj$tab)
  if ("Groups" %in% rownames(df)) return(df["Groups", , drop = FALSE])
  df[1, , drop = FALSE]
}

gen1_gen5_results <- list()

for (stg in c("Larvae", "Female")) {
  keep_ids <- meta_stats %>%
    filter(life_stage == stg, generation %in% c("Generation_1", "Generation_5")) %>%
    pull(sample_id)

  otu_sub <- otu_use[keep_ids, , drop = FALSE]
  meta_sub <- meta_stats %>% filter(sample_id %in% keep_ids) %>% droplevels()
  bray_sub <- vegdist(otu_sub, method = "bray")

  perm_g1_g5 <- adonis2(bray_sub ~ generation, data = meta_sub, permutations = 999)
  bd_g1_g5 <- betadisper(bray_sub, group = meta_sub$generation)
  bd_g1_g5_test <- permutest(bd_g1_g5, permutations = 999)

  gen1_gen5_results[[stg]] <- list(permanova = perm_g1_g5, betadisper = bd_g1_g5_test)
}

gen1_gen5_table <- bind_rows(
  lapply(names(gen1_gen5_results), function(stg) {
    perm_row <- extract_adonis_term(gen1_gen5_results[[stg]]$permanova)
    bd_row <- extract_betadisper_group(gen1_gen5_results[[stg]]$betadisper)
    tibble(
      life_stage = stg,
      permanova_F = perm_row[1, "F"],
      permanova_R2 = perm_row[1, "R2"],
      permanova_P = perm_row[1, "Pr(>F)"],
      betadisper_F = bd_row[1, "F"],
      betadisper_P = bd_row[1, "Pr(>F)"]
    )
  })
) %>%
  mutate(life_stage_label = recode(life_stage, "Larvae" = "Larvae", "Female" = "Adult females"))

write_csv(gen1_gen5_table, file.path(stat_dir, "Figure2_final_Gen1_vs_Gen5_direct_PERMANOVA.csv"))

# -----------------------
# Print everything to console
# -----------------------

cat("\n============================\n")
cat("FIGURE 2 FINAL STATS SUMMARY\n")
cat("============================\n\n")

cat("Sample counts:\n")
print(table(meta_use$generation, meta_use$life_stage))

cat("\nPERMANOVA: Bray-Curtis ~ generation * life_stage\n")
print(perm_time)

cat("\nBetadisper by life stage:\n")
print(bd_stage_test)

cat("\nBetadisper by generation alone:\n")
print(bd_generation_test)

cat("\nBetadisper by combined generation x life_stage groups:\n")
print(bd_gen_stage_test)

cat("\nLinear model for Panel A: DistanceToGen0 ~ generation * life_stage\n")
print(lm_G0_anova)
print(summary(lm_G0))

cat("\nLinear model for Panel B: DistanceToGen1 ~ generation * life_stage\n")
print(lm_G1_anova)
print(summary(lm_G1))

cat("\nDirect Gen1 vs Gen5 PERMANOVA within each life stage:\n")
print(gen1_gen5_table)

for (stg in names(gen1_gen5_results)) {
  cat("\nLife stage:", stg, "\n")
  cat("PERMANOVA Gen1 vs Gen5:\n")
  print(gen1_gen5_results[[stg]]$permanova)
  cat("\nBetadisper Gen1 vs Gen5:\n")
  print(gen1_gen5_results[[stg]]$betadisper)
}

capture.output(
  {
    cat("FIGURE 2 FINAL STATS SUMMARY\n")
    cat("============================\n\n")

    cat("Sample counts:\n")
    print(table(meta_use$generation, meta_use$life_stage))

    cat("\nPERMANOVA: Bray-Curtis ~ generation * life_stage\n")
    print(perm_time)

    cat("\nBetadisper by life stage:\n")
    print(bd_stage_test)

    cat("\nBetadisper by generation alone:\n")
    print(bd_generation_test)

    cat("\nBetadisper by combined generation x life_stage groups:\n")
    print(bd_gen_stage_test)

    cat("\nLinear model for Panel A: DistanceToGen0 ~ generation * life_stage\n")
    print(lm_G0_anova)
    print(summary(lm_G0))

    cat("\nLinear model for Panel B: DistanceToGen1 ~ generation * life_stage\n")
    print(lm_G1_anova)
    print(summary(lm_G1))

    cat("\nDirect Gen1 vs Gen5 PERMANOVA within each life stage:\n")
    print(gen1_gen5_table)

    for (stg in names(gen1_gen5_results)) {
      cat("\nLife stage:", stg, "\n")
      cat("PERMANOVA Gen1 vs Gen5:\n")
      print(gen1_gen5_results[[stg]]$permanova)
      cat("\nBetadisper Gen1 vs Gen5:\n")
      print(gen1_gen5_results[[stg]]$betadisper)
    }
  },
  file = file.path(stat_dir, "Figure2_final_stats_summary.txt")
)

# ============================================================
# SUPPLEMENTARY TABLE S2
# ============================================================

supp_s2_sample_counts <- as.data.frame.matrix(
  table(meta_use$generation, meta_use$life_stage)
) %>%
  rownames_to_column("generation")

supp_s2_G0_lm <- as.data.frame(lm_G0_anova) %>%
  rownames_to_column("term") %>%
  mutate(
    Analysis = "Linear model ANOVA",
    Response = "Bray-Curtis distance to Generation 0 stage-specific centroid",
    Model = "DistanceToGen0 ~ generation * life_stage",
    term = recode(term, "generation" = "generation", "life_stage" = "life stage",
                  "generation:life_stage" = "generation x life stage", "Residuals" = "Residuals")
  ) %>%
  select(Analysis, Response, Model, term, Df, `Sum Sq`, `Mean Sq`, `F value`, `Pr(>F)`)

supp_s2_G1_lm <- as.data.frame(lm_G1_anova) %>%
  rownames_to_column("term") %>%
  mutate(
    Analysis = "Linear model ANOVA",
    Response = "Bray-Curtis distance to Generation 1 stage-specific centroid",
    Model = "DistanceToGen1 ~ generation * life_stage",
    term = recode(term, "generation" = "generation", "life_stage" = "life stage",
                  "generation:life_stage" = "generation x life stage", "Residuals" = "Residuals")
  ) %>%
  select(Analysis, Response, Model, term, Df, `Sum Sq`, `Mean Sq`, `F value`, `Pr(>F)`)

supp_s2_permanova <- as.data.frame(perm_time) %>%
  rownames_to_column("term") %>%
  mutate(
    Analysis = "PERMANOVA",
    Response = "Bray-Curtis dissimilarity among ore;OreR communities",
    Model = "bray_all ~ generation * life_stage",
    term = recode(term, "generation:life_stage" = "generation x life stage",
                  "Residual" = "Residual", "Total" = "Total")
  ) %>%
  select(Analysis, Response, Model, term, Df, SumOfSqs, R2, F, `Pr(>F)`)

supp_s2_betadisper <- as.data.frame(bd_stage_test$tab) %>%
  rownames_to_column("term") %>%
  mutate(
    Analysis = "Betadisper",
    Response = "Distance to life-stage centroid",
    Model = "betadisper(bray_all, group = life_stage)",
    term = recode(term, "Groups" = "life stage", "Residuals" = "Residuals")
  ) %>%
  select(Analysis, Response, Model, term, Df, `Sum Sq`, `Mean Sq`, F, `Pr(>F)`)

supp_s2_betadisper_generation <- as.data.frame(bd_generation_test$tab) %>%
  rownames_to_column("term") %>%
  mutate(
    Analysis = "Betadisper ",
    Response = "Distance to generation centroid",
    Model = "betadisper(bray_all, group = generation)",
    term = recode(term, "Groups" = "generation", "Residuals" = "Residuals")
  ) %>%
  select(Analysis, Response, Model, term, Df, `Sum Sq`, `Mean Sq`, F, `Pr(>F)`)

supp_s2_betadisper_gen_stage <- as.data.frame(bd_gen_stage_test$tab) %>%
  rownames_to_column("term") %>%
  mutate(
    Analysis = "Betadisper ",
    Response = "Distance to generation x life-stage centroid",
    Model = "betadisper(bray_all, group = interaction(generation, life_stage))",
    term = recode(term, "Groups" = "generation x life stage", "Residuals" = "Residuals")
  ) %>%
  select(Analysis, Response, Model, term, Df, `Sum Sq`, `Mean Sq`, F, `Pr(>F)`)

supp_s2_gen1_gen5 <- gen1_gen5_table %>%
  transmute(
    Analysis = "Direct Gen1 vs Gen5 PERMANOVA and betadisper",
    life_stage = life_stage_label,
    permanova_model = "Bray-Curtis ~ generation",
    permanova_F = permanova_F,
    permanova_R2 = permanova_R2,
    permanova_P = permanova_P,
    betadisper_model = "betadisper by generation",
    betadisper_F = betadisper_F,
    betadisper_P = betadisper_P
  )

wb <- createWorkbook()

addWorksheet(wb, "Sample_counts")
writeData(wb, "Sample_counts", supp_s2_sample_counts)

addWorksheet(wb, "Distance_to_G0_LM")
writeData(wb, "Distance_to_G0_LM", supp_s2_G0_lm)

addWorksheet(wb, "Distance_to_G1_LM")
writeData(wb, "Distance_to_G1_LM", supp_s2_G1_lm)

addWorksheet(wb, "PERMANOVA")
writeData(wb, "PERMANOVA", supp_s2_permanova)

addWorksheet(wb, "Betadisper_lifestage")
writeData(wb, "Betadisper_lifestage", supp_s2_betadisper)

addWorksheet(wb, "Betadisper_generation")
writeData(wb, "Betadisper_generation", supp_s2_betadisper_generation)

addWorksheet(wb, "Betadisper_gen_x_stage")
writeData(wb, "Betadisper_gen_x_stage", supp_s2_betadisper_gen_stage)

addWorksheet(wb, "Gen1_vs_Gen5")
writeData(wb, "Gen1_vs_Gen5", supp_s2_gen1_gen5)

header_style <- createStyle(textDecoration = "bold", halign = "center")

for (sheet in names(wb)) {
  addStyle(wb, sheet = sheet, style = header_style, rows = 1, cols = 1:50, gridExpand = TRUE)
  setColWidths(wb, sheet = sheet, cols = 1:50, widths = "auto")
}

saveWorkbook(
  wb,
  file = file.path(stat_dir, "Supplementary_Table_S2_oreOreR_generation_lifestage_statistics.xlsx"),
  overwrite = TRUE
)

write_csv(supp_s2_sample_counts, file.path(stat_dir, "Supplementary_Table_S2A_sample_counts.csv"))
write_csv(supp_s2_G0_lm, file.path(stat_dir, "Supplementary_Table_S2B_distance_to_G0_LM.csv"))
write_csv(supp_s2_G1_lm, file.path(stat_dir, "Supplementary_Table_S2C_distance_to_G1_LM.csv"))
write_csv(supp_s2_permanova, file.path(stat_dir, "Supplementary_Table_S2D_PERMANOVA.csv"))
write_csv(supp_s2_betadisper, file.path(stat_dir, "Supplementary_Table_S2E_betadisper_lifestage.csv"))
write_csv(supp_s2_betadisper_generation, file.path(stat_dir, "Supplementary_Table_S2G_betadisper_generation.csv"))
write_csv(supp_s2_betadisper_gen_stage, file.path(stat_dir, "Supplementary_Table_S2H_betadisper_gen_x_stage.csv"))
write_csv(supp_s2_gen1_gen5, file.path(stat_dir, "Supplementary_Table_S2F_Gen1_vs_Gen5.csv"))

# ============================================================
# DONE
# ============================================================

cat("\nDONE. Generation divergence outputs saved:\n")
cat(file.path(stat_dir, "Supplementary_Table_S2G_betadisper_generation.csv"), "\n")
cat(file.path(stat_dir, "Supplementary_Table_S2H_betadisper_gen_x_stage.csv"), "\n")

# Record software versions used for this run.
capture.output(sessionInfo(), file = file.path(main_outdir, "sessionInfo.txt"))
