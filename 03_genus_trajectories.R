# Genus abundance trajectories
# Plot dominant and filtered bacterial genera across generations and vials.
# Sum genera within each sample before averaging the Other category.
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

required_packages <- c("phyloseq", "dplyr", "tidyr", "ggplot2", "stringr", "forcats", "scales")
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
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(scales)
})

set.seed(1234)

# -----------------------
# 2) Create output folders
# -----------------------

main_outdir <- file.path(OUTPUT_ROOT, "Fig3_oreOreR_genus_trajectories")

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
# 3) Load processed phyloseq object
# -----------------------

ps.clean <- readRDS(INPUT_RDS)

if (!inherits(ps.clean, "phyloseq")) {
  stop("INPUT_RDS must contain a phyloseq object.", call. = FALSE)
}
required_metadata <- c("genotype", "generation", "life_stage", "vial")
missing_metadata <- setdiff(required_metadata, names(data.frame(sample_data(ps.clean))))
if (length(missing_metadata) > 0) {
  stop("Missing sample metadata: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
}
missing_ranks <- setdiff(c("Genus", "Family"), rank_names(ps.clean))
if (length(missing_ranks) > 0) {
  stop("Missing taxonomy ranks: ", paste(missing_ranks, collapse = ", "), call. = FALSE)
}

saveRDS(
  ps.clean,
  file = file.path(input_dir, "06_phyloseq_clean_CHAP1_NOHOST.rds")
)

# -------------------------
# Helper functions
# -------------------------

fix_genus <- function(df) {
  df %>%
    mutate(
      Genus  = as.character(Genus),
      Family = as.character(Family),
      Genus  = ifelse(is.na(Genus) | Genus == "", Family, Genus)
    )
}

to_rel_abund <- function(psx) {
  transform_sample_counts(psx, function(x) if (sum(x) > 0) x / sum(x) else x)
}

# ============================================================
# PART A: FIGURE 3 FINAL
# Dominant genus-level trajectories across replicate vials
# ============================================================

# -------------------------
# A1) FOOD inoculum
# -------------------------

ps.food <- subset_samples(
  ps.clean,
  genotype == "oreore" &
    life_stage == "food" &
    vial %in% c("Vial_1", "Vial_2")
)

ps.food <- prune_samples(sample_sums(ps.food) > 0, ps.food)

ps.food.ra    <- to_rel_abund(ps.food)
# The original taxonomy and Gardnerella exclusions are retained; see CHANGES.md.
ps.food.genus <- tax_glom(ps.food.ra, taxrank = "Genus", NArm = TRUE)

df.food <- psmelt(ps.food.genus) %>%
  fix_genus() %>%
  filter(Genus != "Gardnerella") %>%
  mutate(
    GenerationPlot = "Inoculum",
    Vial  = str_replace(as.character(vial), "_", " "),
    Stage = "food"
  )

# -------------------------
# A2) OREORE flies
# -------------------------

ps.ore <- subset_samples(
  ps.clean,
  genotype == "oreore" &
    generation %in% paste0("Generation_", 0:5) &
    life_stage %in% c("Larvae", "Female")
)

ps.ore <- prune_samples(sample_sums(ps.ore) > 0, ps.ore)

ps.ore.ra    <- to_rel_abund(ps.ore)
ps.ore.genus <- tax_glom(ps.ore.ra, taxrank = "Genus", NArm = TRUE)

df.fly <- psmelt(ps.ore.genus) %>%
  fix_genus() %>%
  filter(Genus != "Gardnerella") %>%
  mutate(
    GenerationPlot = recode(
      generation,
      "Generation_0" = "Gen 0",
      "Generation_1" = "Gen 1",
      "Generation_2" = "Gen 2",
      "Generation_3" = "Gen 3",
      "Generation_4" = "Gen 4",
      "Generation_5" = "Gen 5"
    ),
    Vial  = str_replace(as.character(vial), "_", " "),
    Stage = as.character(life_stage)
  )

# -------------------------
# A3) Top genera
# -------------------------

topN_fig3 <- 5

top_genera_fig3 <- bind_rows(
  df.fly  %>% select(Genus, Abundance),
  df.food %>% select(Genus, Abundance)
) %>%
  group_by(Genus) %>%
  summarise(total_abund = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_abund)) %>%
  slice_head(n = topN_fig3) %>%
  pull(Genus)

cat("\nTop genera shown individually in Figure 3:\n")
print(top_genera_fig3)

df.fly2  <- df.fly  %>% mutate(Genus2 = ifelse(Genus %in% top_genera_fig3, Genus, "Other"))
df.food2 <- df.food %>% mutate(Genus2 = ifelse(Genus %in% top_genera_fig3, Genus, "Other"))

# -------------------------
# A4) Expand food inoculum into both larval and adult female rows
# -------------------------

food_expanded <- df.food2 %>%
  mutate(
    Vial = str_replace(as.character(vial), "_", " ")
  ) %>%
  tidyr::crossing(
    StagePretty = factor(
      c("A  Larvae", "B  Adult females"),
      levels = c("A  Larvae", "B  Adult females")
    )
  )

fly_expanded <- df.fly2 %>%
  mutate(
    Vial = str_replace(as.character(vial), "_", " "),
    StagePretty = case_when(
      Stage == "Larvae" ~ "A  Larvae",
      Stage == "Female" ~ "B  Adult females"
    ),
    StagePretty = factor(
      StagePretty,
      levels = c("A  Larvae", "B  Adult females")
    )
  )

plot_df_fig3 <- bind_rows(fly_expanded, food_expanded) %>%
  mutate(
    GenerationPlot = factor(
      GenerationPlot,
      levels = c("Inoculum", "Gen 0", "Gen 1", "Gen 2", "Gen 3", "Gen 4", "Gen 5")
    ),
    Vial = factor(Vial, levels = c("Vial 1", "Vial 2", "Vial 3"))
  ) %>%
  # Sum constituent genera within each pool, then average pools within a vial.
  group_by(StagePretty, GenerationPlot, Vial, Sample, Genus2) %>%
  summarise(
    pool_abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(StagePretty, GenerationPlot, Vial, Genus2) %>%
  summarise(
    mean_abundance = mean(pool_abundance, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------
# A5) Colors
# -------------------------

genus_palette_fig3 <- c(
  "Acetobacter" = "#C23B22",
  "Lactococcus" = "#6A4C93",
  "Levilactobacillus" = "#F4D03F",
  "Lactiplantibacillus" = "#3A86FF",
  "Lactobacillus" = "#00A878",
  "Other" = "grey85"
)

present_genera_fig3 <- unique(as.character(plot_df_fig3$Genus2))
present_genera_fig3 <- c(
  names(genus_palette_fig3)[names(genus_palette_fig3) %in% present_genera_fig3],
  setdiff(present_genera_fig3, names(genus_palette_fig3))
)

plot_df_fig3$Genus2 <- factor(plot_df_fig3$Genus2, levels = present_genera_fig3)

pal_fig3 <- genus_palette_fig3[present_genera_fig3]

missing_fig3 <- setdiff(present_genera_fig3, names(genus_palette_fig3))
if (length(missing_fig3) > 0) {
  pal_fig3[missing_fig3] <- scales::hue_pal()(length(missing_fig3))
}

# -------------------------
# A6) Final Figure 3
# -------------------------

p_final_fig3 <- ggplot(
  plot_df_fig3,
  aes(
    x = GenerationPlot,
    y = mean_abundance,
    group = Genus2,
    color = Genus2
  )
) +
  geom_line(linewidth = 1.2, alpha = 0.95) +
  geom_point(size = 2.2, alpha = 0.95) +
  facet_grid(
    rows = vars(StagePretty),
    cols = vars(Vial),
    scales = "fixed"
  ) +
  scale_color_manual(values = pal_fig3, name = "Genus") +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.02)),
    labels = percent_format(accuracy = 1)
  ) +
  labs(
    x = "Generation",
    y = "Mean relative abundance",
    title = "Dominant bacterial genera across generations in the compatible control genotype"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(face = "bold", color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black"),

    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
    strip.text = element_text(face = "bold", size = 13, color = "black"),

    panel.spacing.x = unit(1.0, "lines"),
    panel.spacing.y = unit(1.2, "lines"),

    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),

    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11),
    legend.key.width = unit(1.4, "lines"),

    panel.grid = element_blank(),
    plot.margin = margin(t = 14, r = 12, b = 12, l = 12)
  )

p_final_fig3

# -------------------------
# A7) Save Figure 3 and source data
# -------------------------

ggsave(
  file.path(fig_dir, "Figure3_oreore_genus_trajectories_combined.png"),
  p_final_fig3,
  width = 12,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(fig_dir, "Figure3_oreore_genus_trajectories_combined.jpeg"),
  p_final_fig3,
  width = 12,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(fig_dir, "Figure3_oreore_genus_trajectories_combined.pdf"),
  p_final_fig3,
  width = 12,
  height = 8.5,
  device = PDF_DEVICE,
  bg = "white"
)

write.csv(
  plot_df_fig3,
  file.path(data_dir, "Figure3_oreore_genus_trajectories_combined_data.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(top_genera_fig3 = top_genera_fig3),
  file.path(data_dir, "Figure3_top5_genera_list.csv"),
  row.names = FALSE
)

# ============================================================
# PART B: SUPPLEMENTARY FIGURE S3
# Full genus-level trajectories across replicate vials
# Compatible control genotype: (ore);OreR
# ============================================================

# -------------------------
# B1) Subset compatible control genotype
# -------------------------
# Keeps only (ore);OreR, larvae + adult females, Generation 0-5

ps.s3 <- subset_samples(
  ps.clean,
  genotype == "oreore" &
    life_stage %in% c("Larvae", "Female") &
    generation %in% paste0("Generation_", 0:5)
)

ps.s3 <- prune_samples(sample_sums(ps.s3) > 0, ps.s3)
ps.s3 <- prune_taxa(taxa_sums(ps.s3) > 0, ps.s3)

cat("\nSamples included in Supplementary Fig. S3:\n")
print(
  data.frame(sample_data(ps.s3)) %>%
    count(life_stage, generation, vial)
)

cat("\nNumber of samples:", nsamples(ps.s3), "\n")
cat("Number of taxa:", ntaxa(ps.s3), "\n")

# -------------------------
# B2) Transform to relative abundance and aggregate by genus
# -------------------------

ps.s3.ra <- transform_sample_counts(ps.s3, function(x) x / sum(x))
ps.s3.ra <- prune_taxa(taxa_sums(ps.s3.ra) > 0, ps.s3.ra)

ps.s3.genus <- tax_glom(ps.s3.ra, taxrank = "Genus", NArm = TRUE)

df.s3 <- psmelt(ps.s3.genus) %>%
  fix_genus() %>%
  filter(Genus != "Gardnerella") %>%
  mutate(
    Genus = ifelse(
      Genus == "Burkholderia-Caballeronia-Paraburkholderia",
      "BCP",
      Genus
    ),
    Generation_num = as.numeric(str_extract(as.character(generation), "\\d+")),
    Generation = factor(
      paste0("G", Generation_num),
      levels = paste0("G", 0:5)
    ),
    Stage = dplyr::recode(
      as.character(life_stage),
      "Larvae" = "Larvae",
      "Female" = "Adult females"
    ),
    Stage = factor(Stage, levels = c("Larvae", "Adult females")),
    Vial = str_replace(as.character(vial), "Vial_", "Vial "),
    Vial = factor(Vial, levels = c("Vial 1", "Vial 2", "Vial 3"))
  )

# -------------------------
# B3) Vial-level means
# -------------------------
# This prevents pooled samples from same vial being treated as independent
# for trajectory visualization.

df.vial <- df.s3 %>%
  group_by(Stage, Vial, Generation, Generation_num, Genus) %>%
  summarise(
    Relative_abundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

# -------------------------
# B4) Apply genus filter, then select top genera for plotting
# -------------------------
# Filter rule:
# retain genera that reach >=0.5% relative abundance in at least two replicate vials
# within any Stage × Generation group.
# Because this figure is only (ore);OreR, genotype is fixed.

GENUS_MIN_ABUND <- 0.005
MIN_VIALS <- 2
topN_s3 <- 10

retained_genera <- df.vial %>%
  group_by(Stage, Generation, Genus) %>%
  summarise(
    n_vials_ge_threshold = sum(Relative_abundance >= GENUS_MIN_ABUND, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Genus) %>%
  summarise(
    keep = any(n_vials_ge_threshold >= MIN_VIALS),
    .groups = "drop"
  ) %>%
  filter(keep) %>%
  pull(Genus)

cat("\nGenera retained by 0.5% in at least 2 vials rule:\n")
print(retained_genera)

top_genera_s3 <- df.vial %>%
  filter(Genus %in% retained_genera) %>%
  group_by(Genus) %>%
  summarise(
    mean_abundance = mean(Relative_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = topN_s3) %>%
  pull(Genus)

cat("\nTop genera shown individually in Supplementary Fig. S3:\n")
print(top_genera_s3)

df.plot.s3 <- df.vial %>%
  filter(Genus %in% retained_genera) %>%
  mutate(
    Genus2 = ifelse(Genus %in% top_genera_s3, Genus, "Other")
  ) %>%
  group_by(Stage, Vial, Generation, Generation_num, Genus2) %>%
  summarise(
    Relative_abundance = sum(Relative_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Stage, Vial, Generation, Generation_num) %>%
  mutate(
    Relative_abundance = Relative_abundance / sum(Relative_abundance)
  ) %>%
  ungroup()

# -------------------------
# B5) Genus order and palette
# -------------------------

preferred_order <- c(
  "Acetobacter",
  "Lactococcus",
  "Levilactobacillus",
  "Lactiplantibacillus",
  "Lactobacillus",
  "BCP",
  "Limosilactobacillus",
  "Liquorilactobacillus",
  "Enterobacter",
  "Pseudoclavibacter",
  "Prevotella_9",
  "Other"
)

extra_genera <- setdiff(unique(as.character(df.plot.s3$Genus2)), preferred_order)

genus_order_s3 <- c(
  preferred_order[preferred_order %in% unique(as.character(df.plot.s3$Genus2))],
  extra_genera
)

genus_order_s3 <- c(setdiff(genus_order_s3, "Other"), "Other")
genus_order_s3 <- genus_order_s3[genus_order_s3 %in% unique(as.character(df.plot.s3$Genus2))]

df.plot.s3$Genus2 <- factor(df.plot.s3$Genus2, levels = genus_order_s3)

genus_palette_s3 <- c(
  "Acetobacter" = "#C23B22",
  "Lactococcus" = "#6A4C93",
  "Levilactobacillus" = "#F4D03F",
  "Lactiplantibacillus" = "#3A86FF",
  "Lactobacillus" = "#E69F00",
  "BCP" = "black",
  "Limosilactobacillus" = "#AA4499",
  "Liquorilactobacillus" = "#577590",
  "Enterobacter" = "#BC6C25",
  "Pseudoclavibacter" = "#F28482",
  "Prevotella_9" = "#90BE6D",
  "Other" = "grey85"
)

missing_cols_s3 <- setdiff(genus_order_s3, names(genus_palette_s3))

if (length(missing_cols_s3) > 0) {
  extra_cols_s3 <- scales::hue_pal()(length(missing_cols_s3))
  names(extra_cols_s3) <- missing_cols_s3
  genus_palette_s3 <- c(genus_palette_s3, extra_cols_s3)
}

# -------------------------
# B6) Final Supplementary Figure S3
# -------------------------

S3_trajectories <- ggplot(
  df.plot.s3,
  aes(
    x = Generation_num,
    y = Relative_abundance,
    color = Genus2,
    group = Genus2
  )
) +
  geom_line(
    linewidth = 0.9,
    alpha = 0.95
  ) +
  geom_point(
    size = 2.2,
    alpha = 0.95
  ) +
  facet_grid(
    Stage ~ Vial
  ) +
  scale_color_manual(
    values = genus_palette_s3,
    breaks = genus_order_s3,
    name = "Genus",
    guide = guide_legend(
      nrow = 2,
      byrow = TRUE
    )
  ) +
  scale_x_continuous(
    breaks = 0:5,
    labels = paste0("G", 0:5),
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(
    x = "Generation",
    y = "Relative abundance"
  ) +
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.7
    ),
    strip.text = element_text(
      face = "bold",
      size = 12,
      color = "black"
    ),

    axis.title = element_text(
      face = "bold",
      color = "black",
      size = 13
    ),
    axis.title.x = element_text(
      margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      margin = margin(r = 8)
    ),

    axis.text = element_text(
      color = "black",
      size = 10.5
    ),

    legend.position = "bottom",
    legend.title = element_text(
      face = "bold",
      size = 11
    ),
    legend.text = element_text(
      size = 9.5,
      color = "black"
    ),
    legend.key.size = unit(0.42, "cm"),
    legend.spacing.x = unit(0.22, "cm"),

    panel.grid.major.y = element_line(
      color = "grey90",
      linewidth = 0.25
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),

    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.65
    ),
    axis.line = element_blank(),

    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),

    plot.margin = margin(
      t = 10,
      r = 12,
      b = 10,
      l = 12
    )
  )

S3_trajectories

# -------------------------
# B7) Save Supplementary Figure S3 and source data
# -------------------------

ggsave(
  file.path(fig_dir, "Supplementary_Figure_S3_full_genus_trajectories_oreOreR.pdf"),
  S3_trajectories,
  width = 12.5,
  height = 7.2,
  units = "in",
  device = PDF_DEVICE,
  bg = "white"
)

ggsave(
  file.path(fig_dir, "Supplementary_Figure_S3_full_genus_trajectories_oreOreR.png"),
  S3_trajectories,
  width = 12.5,
  height = 7.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  file.path(fig_dir, "Supplementary_Figure_S3_full_genus_trajectories_oreOreR.jpeg"),
  S3_trajectories,
  width = 12.5,
  height = 7.2,
  units = "in",
  dpi = 600,
  bg = "white",
  device = "jpeg"
)

write.csv(
  df.plot.s3,
  file.path(data_dir, "Supplementary_Figure_S3_full_genus_trajectories_data.csv"),
  row.names = FALSE
)

write.csv(
  df.plot.s3,
  file.path(data_dir, "Supplementary_Figure_S3_plotted_top10_plus_other.csv"),
  row.names = FALSE
)

write.csv(
  df.vial,
  file.path(data_dir, "Supplementary_Figure_S3_full_vial_level_genus_relative_abundance.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(retained_genera = retained_genera),
  file.path(data_dir, "Supplementary_Figure_S3_retained_genera_filter_list.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(top_genera_s3 = top_genera_s3),
  file.path(data_dir, "Supplementary_Figure_S3_top10_genera_list.csv"),
  row.names = FALSE
)

# ============================================================
# PART C: SAVE SUMMARY / PARAMETERS
# ============================================================

capture.output(
  {
    cat("FIGURE 3 + SUPPLEMENTARY FIGURE S3 SUMMARY\n")
    cat("==========================================\n\n")

    cat("Input object:\n")
    cat("06_phyloseq_clean_CHAP1_NOHOST.rds\n\n")

    cat("Figure 3 topN:\n")
    cat(topN_fig3, "\n\n")

    cat("Top genera shown individually in Figure 3:\n")
    print(top_genera_fig3)

    cat("\nSupplementary Fig. S3 filtering parameters:\n")
    cat("GENUS_MIN_ABUND =", GENUS_MIN_ABUND, "\n")
    cat("MIN_VIALS =", MIN_VIALS, "\n")
    cat("topN_s3 =", topN_s3, "\n\n")

    cat("Samples included in Supplementary Fig. S3:\n")
    print(
      data.frame(sample_data(ps.s3)) %>%
        count(life_stage, generation, vial)
    )

    cat("\nNumber of samples in ps.s3:", nsamples(ps.s3), "\n")
    cat("Number of taxa in ps.s3:", ntaxa(ps.s3), "\n")

    cat("\nGenera retained by 0.5% in at least 2 vials rule:\n")
    print(retained_genera)

    cat("\nTop genera shown individually in Supplementary Fig. S3:\n")
    print(top_genera_s3)
  },
  file = file.path(stat_dir, "Figure3_and_SuppFigS3_summary.txt")
)

# ============================================================
# DONE
# ============================================================

cat("\nDONE. Files saved inside folder:\n")
cat(main_outdir, "\n\n")

cat("Figure files:\n")
cat(file.path(fig_dir, "Figure3_oreore_genus_trajectories_combined.pdf"), "\n")
cat(file.path(fig_dir, "Figure3_oreore_genus_trajectories_combined.png"), "\n")
cat(file.path(fig_dir, "Figure3_oreore_genus_trajectories_combined.jpeg"), "\n")
cat(file.path(fig_dir, "Supplementary_Figure_S3_full_genus_trajectories_oreOreR.pdf"), "\n")
cat(file.path(fig_dir, "Supplementary_Figure_S3_full_genus_trajectories_oreOreR.png"), "\n")
cat(file.path(fig_dir, "Supplementary_Figure_S3_full_genus_trajectories_oreOreR.jpeg"), "\n\n")

cat("Figure source data:\n")
cat(file.path(data_dir, "Figure3_oreore_genus_trajectories_combined_data.csv"), "\n")
cat(file.path(data_dir, "Figure3_top5_genera_list.csv"), "\n")
cat(file.path(data_dir, "Supplementary_Figure_S3_full_genus_trajectories_data.csv"), "\n")
cat(file.path(data_dir, "Supplementary_Figure_S3_plotted_top10_plus_other.csv"), "\n")
cat(file.path(data_dir, "Supplementary_Figure_S3_full_vial_level_genus_relative_abundance.csv"), "\n")
cat(file.path(data_dir, "Supplementary_Figure_S3_retained_genera_filter_list.csv"), "\n")
cat(file.path(data_dir, "Supplementary_Figure_S3_top10_genera_list.csv"), "\n\n")

cat("Summary:\n")
cat(file.path(stat_dir, "Figure3_and_SuppFigS3_summary.txt"), "\n\n")

cat("Processed input copy:\n")
cat(file.path(input_dir, "06_phyloseq_clean_CHAP1_NOHOST.rds"), "\n")

# Record software versions used for this run.
capture.output(sessionInfo(), file = file.path(main_outdir, "sessionInfo.txt"))
