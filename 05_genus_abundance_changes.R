# Genus abundance changes between Gen0 and Gen5
# Calculate vial-level abundance changes and genotype-by-stage interactions.
# Keep all retained genera in the statistical screen; top ten are for plotting.
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

required_packages <- c("phyloseq", "dplyr", "tidyr", "ggplot2", "readr", "tibble", "stringr", "scales", "openxlsx")
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
  library(readr)
  library(tibble)
  library(stringr)
  library(scales)
  library(openxlsx)
})

set.seed(1234)

# -------------------------
# 2) Create output folders
# -------------------------

main_outdir <- file.path(OUTPUT_ROOT, "Fig5_genus_level_delta")

dir.create(main_outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(main_outdir, "input_processed"), showWarnings = FALSE)
dir.create(file.path(main_outdir, "figure_files"), showWarnings = FALSE)
dir.create(file.path(main_outdir, "figure_source_data"), showWarnings = FALSE)
dir.create(file.path(main_outdir, "statistics"), showWarnings = FALSE)

fig_dir   <- file.path(main_outdir, "figure_files")
data_dir  <- file.path(main_outdir, "figure_source_data")
stat_dir  <- file.path(main_outdir, "statistics")
input_dir <- file.path(main_outdir, "input_processed")

# -------------------------
# 3) Load phyloseq object
# -------------------------

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

# ============================================================
# PART A: FIGURE 5 FINAL — TWO FOCAL GENERA
# ============================================================

# -------------------------
# A1) Choose focal genera
# -------------------------

focal_genera <- c(
  "Lactiplantibacillus",
  "Acetobacter"
)

# -------------------------
# A2) Subset Gen0 and Gen5 host samples
# -------------------------

ps.sub.fig5 <- subset_samples(
  ps.clean,
  life_stage %in% c("Larvae", "Female") &
    generation %in% c("Generation_0", "Generation_5") &
    genotype %in% c("oreore", "simore", "oreaut", "simaut")
)

ps.sub.fig5 <- prune_samples(sample_sums(ps.sub.fig5) > 0, ps.sub.fig5)

# -------------------------
# A3) Transform to relative abundance and collapse to genus
# -------------------------

ps.sub.fig5.ra <- transform_sample_counts(
  ps.sub.fig5,
  function(x) if (sum(x) > 0) x / sum(x) else x
)

ps.genus.fig5 <- tax_glom(ps.sub.fig5.ra, taxrank = "Genus", NArm = TRUE)

# -------------------------
# A4) Melt and clean metadata
# -------------------------

df_long_fig5 <- psmelt(ps.genus.fig5) %>%
  mutate(
    Genus = as.character(Genus),
    Family = as.character(Family),
    Genus = ifelse(is.na(Genus) | Genus == "", Family, Genus),

    generation = recode(
      generation,
      "Generation_0" = "Gen0",
      "Generation_5" = "Gen5"
    ),
    generation = factor(generation, levels = c("Gen0", "Gen5")),

    genotype = factor(
      genotype,
      levels = c("oreore", "simore", "oreaut", "simaut")
    ),

    genotype_label = recode(
      as.character(genotype),
      "oreore" = "(ore);OreR",
      "simore" = "(simw501);OreR",
      "oreaut" = "(ore);Aut",
      "simaut" = "(simw501);Aut"
    ),

    life_stage = factor(
      life_stage,
      levels = c("Larvae", "Female")
    ),

    life_stage_label = recode(
      as.character(life_stage),
      "Larvae" = "Larvae",
      "Female" = "Adult females"
    ),
    life_stage_label = factor(
      life_stage_label,
      levels = c("Larvae", "Adult females")
    ),

    vial = factor(vial),

    Mito = case_when(
      genotype %in% c("oreore", "oreaut") ~ "ore",
      genotype %in% c("simore", "simaut") ~ "simw501",
      TRUE ~ NA_character_
    ),

    Nuclear = case_when(
      genotype %in% c("oreore", "simore") ~ "OreR",
      genotype %in% c("oreaut", "simaut") ~ "Aut",
      TRUE ~ NA_character_
    ),

    Mito = factor(Mito, levels = c("ore", "simw501")),
    Nuclear = factor(Nuclear, levels = c("OreR", "Aut"))
  )

# -------------------------
# A5) Check sample/vial structure
# -------------------------

sample_count_table_fig5 <- df_long_fig5 %>%
  distinct(genotype, life_stage_label, generation, vial) %>%
  count(genotype, life_stage_label, generation, name = "n_vials") %>%
  arrange(life_stage_label, genotype, generation)

cat("\n================ FIGURE 5 SAMPLE COUNTS ================\n")
print(as.data.frame(sample_count_table_fig5), row.names = FALSE)

write_csv(
  sample_count_table_fig5,
  file.path(data_dir, "Figure5_two_genera_sample_counts.csv")
)

# -------------------------
# A6) Collapse to vial-level genus abundance
# -------------------------

df_vial_fig5 <- df_long_fig5 %>%
  group_by(
    Genus,
    genotype,
    genotype_label,
    Mito,
    Nuclear,
    life_stage,
    life_stage_label,
    generation,
    vial
  ) %>%
  summarise(
    Abundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  df_vial_fig5,
  file.path(data_dir, "Figure5_two_genera_full_vial_level_genus_relative_abundance.csv")
)

# -------------------------
# A7) Check focal genera exist
# -------------------------

available_genera <- sort(unique(df_vial_fig5$Genus))
missing_genera <- setdiff(focal_genera, available_genera)

if (length(missing_genera) > 0) {
  stop(
    paste(
      "These focal genera were not found in the genus table:",
      paste(missing_genera, collapse = ", ")
    )
  )
}

# -------------------------
# A8) Compute delta: Gen5 - Gen0
# -------------------------

df_delta <- df_vial_fig5 %>%
  filter(Genus %in% focal_genera) %>%
  select(
    Genus,
    genotype,
    genotype_label,
    Mito,
    Nuclear,
    life_stage,
    life_stage_label,
    generation,
    vial,
    Abundance
  ) %>%
  pivot_wider(
    names_from = generation,
    values_from = Abundance
  ) %>%
  drop_na(Gen0, Gen5) %>%
  mutate(
    delta = Gen5 - Gen0,
    Genus = factor(Genus, levels = focal_genera),
    Mito_label = recode(
      as.character(Mito),
      "ore" = "(ore)",
      "simw501" = "(simw501)"
    ),
    Mito_label = factor(
      Mito_label,
      levels = c("(ore)", "(simw501)")
    )
  ) %>%
  arrange(Genus, life_stage_label, Nuclear, genotype, vial)

cat("\n================ TWO-GENERA DELTA TABLE ================\n")
print(as.data.frame(df_delta), row.names = FALSE)

write_csv(
  df_delta,
  file.path(data_dir, "Figure5_two_genera_delta_raw.csv")
)

# -------------------------
# A9) Per-genus statistics
# Model: delta ~ Mito * Nuclear * life_stage
# -------------------------

stats_two <- df_delta %>%
  group_by(Genus) %>%
  group_modify(~{

    fit <- lm(delta ~ Mito * Nuclear * life_stage, data = .x)
    dr <- drop1(fit, test = "F")

    data.frame(
      p_3way = if ("Mito:Nuclear:life_stage" %in% rownames(dr)) {
        dr["Mito:Nuclear:life_stage", "Pr(>F)"]
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }) %>%
  ungroup() %>%
  mutate(
    padj_3way = p.adjust(p_3way, method = "BH"),
    label = case_when(
      is.na(padj_3way) ~ "NA",
      padj_3way < 0.001 ~ "BH adj. P < 0.001",
      padj_3way < 0.01 ~ "BH adj. P < 0.01",
      padj_3way < 0.05 ~ "BH adj. P < 0.05",
      TRUE ~ paste0("BH adj. P = ", signif(padj_3way, 2))
    )
  )

cat("\n================ TWO-GENERA MODEL STATS ================\n")
print(as.data.frame(stats_two), row.names = FALSE)

write_csv(
  stats_two,
  file.path(stat_dir, "Figure5_two_genera_stats.csv")
)

# -------------------------
# A10) Direction and magnitude summary
# -------------------------

direction_summary <- df_delta %>%
  group_by(Genus, genotype, genotype_label, Mito, Nuclear, life_stage_label) %>%
  summarise(
    mean_delta = mean(delta, na.rm = TRUE),
    median_delta = median(delta, na.rm = TRUE),
    sd_delta = sd(delta, na.rm = TRUE),
    consistency = mean(sign(delta) == sign(mean(delta, na.rm = TRUE))),
    direction = case_when(
      mean_delta > 0 ~ "Increase",
      mean_delta < 0 ~ "Decrease",
      TRUE ~ "Neutral"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    abs_mean_delta = abs(mean_delta),
    effect_strength = abs_mean_delta * consistency
  ) %>%
  arrange(Genus, life_stage_label, Nuclear, Mito)

cat("\n================ TWO-GENERA DIRECTION SUMMARY ================\n")
print(as.data.frame(direction_summary), row.names = FALSE)

write_csv(
  direction_summary,
  file.path(data_dir, "Figure5_two_genera_direction_summary.csv")
)

# -------------------------
# A11) Figure 5 plot
# -------------------------

mito_fill_cols <- c(
  "(ore)" = "#B8C6A6",
  "(simw501)" = "#D9B382"
)

mito_line_cols <- c(
  "(ore)" = "#3F5F3A",
  "(simw501)" = "#8A5A00"
)

y_lim_fig5 <- max(abs(df_delta$delta), na.rm = TRUE) * 1.15

p_two <- ggplot(
  df_delta,
  aes(
    x = Mito_label,
    y = delta,
    fill = Mito_label,
    color = Mito_label
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "black"
  ) +
  geom_point(
    aes(shape = Mito_label),
    position = position_jitter(width = 0.08, height = 0),
    size = 2.65,
    alpha = 0.90,
    stroke = 0.60
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 3.35,
    fill = "white",
    color = "black",
    stroke = 0.70,
    show.legend = FALSE
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.15,
    linewidth = 0.70,
    color = "black"
  ) +
  facet_grid(
    life_stage_label + Nuclear ~ Genus,
    scales = "fixed"
  ) +
  scale_fill_manual(
    values = mito_fill_cols,
    labels = c(
      "(ore)" = "(ore)",
      "(simw501)" = expression("(simw"^501*")")
    ),
    name = "Mitochondrial genotype"
  ) +
  scale_color_manual(
    values = mito_line_cols,
    labels = c(
      "(ore)" = "(ore)",
      "(simw501)" = expression("(simw"^501*")")
    ),
    name = "Mitochondrial genotype"
  ) +
  scale_shape_manual(
    values = c(
      "(ore)" = 21,
      "(simw501)" = 24
    ),
    labels = c(
      "(ore)" = "(ore)",
      "(simw501)" = expression("(simw"^501*")")
    ),
    name = "Mitochondrial genotype"
  ) +
  scale_y_continuous(
    limits = c(-y_lim_fig5, y_lim_fig5),
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  labs(
    x = "Mitochondrial genotype",
    y = expression(Delta * " relative abundance (Generation 5 - Generation 0)")
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = c(21, 24),
        color = mito_line_cols,
        fill = mito_fill_cols,
        size = 3
      )
    ),
    color = "none",
    shape = "none"
  ) +
  theme_classic(base_size = 12.5) +
  theme(
    axis.title = element_text(
      face = "bold",
      color = "black",
      size = 11.8
    ),
    axis.title.x = element_text(
      margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      margin = margin(r = 8)
    ),
    axis.text = element_text(
      color = "black",
      size = 10
    ),
    axis.text.x = element_text(
      face = "bold",
      angle = 22,
      hjust = 1
    ),
    axis.line = element_blank(),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.4
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.55
    ),
    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.55
    ),
    strip.text = element_text(
      face = "bold",
      size = 10.8,
      color = "black"
    ),
    panel.spacing.x = unit(0.95, "lines"),
    panel.spacing.y = unit(0.85, "lines"),
    legend.position = "top",
    legend.title = element_text(
      face = "bold",
      size = 10.5
    ),
    legend.text = element_text(
      size = 9.8
    ),
    legend.key.size = unit(0.55, "cm"),
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
      r = 10,
      b = 10,
      l = 10
    )
  )

p_two

ggsave(
  filename = file.path(fig_dir, "Figure5_two_genera_delta_ISME.pdf"),
  plot = p_two,
  width = 8.6,
  height = 7.4,
  units = "in",
  device = PDF_DEVICE,
  bg = "white"
)

ggsave(
  filename = file.path(fig_dir, "Figure5_two_genera_delta_ISME.png"),
  plot = p_two,
  width = 8.6,
  height = 7.4,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(fig_dir, "Figure5_two_genera_delta_ISME.jpeg"),
  plot = p_two,
  width = 8.6,
  height = 7.4,
  units = "in",
  dpi = 600,
  bg = "white",
  device = "jpeg"
)

# -------------------------
# A12) Save Figure 5 summary
# -------------------------

capture.output(
  {
    cat("FIGURE 5 TWO-GENERA FINAL SUMMARY\n")
    cat("=================================\n\n")

    cat("Sample counts:\n")
    print(sample_count_table_fig5)

    cat("\nFocal genera:\n")
    print(focal_genera)

    cat("\nModel statistics:\n")
    print(stats_two)

    cat("\nDirection summary:\n")
    print(direction_summary)
  },
  file = file.path(stat_dir, "Figure5_two_genera_summary.txt")
)

# ============================================================
# PART B: SUPPLEMENTARY FIGURE S4
# Expanded genus-level Gen5 - Gen0 delta screen
# ============================================================

# -------------------------
# B1) Parameters
# -------------------------

GENUS_MIN_ABUND <- 0.005
MIN_VIALS <- 2
topN <- 10

# -------------------------
# B2) Subset Gen0 and Gen5 host samples
# -------------------------

ps.sub.s4 <- subset_samples(
  ps.clean,
  life_stage %in% c("Larvae", "Female") &
    generation %in% c("Generation_0", "Generation_5") &
    genotype %in% c("oreore", "simore", "oreaut", "simaut")
)

ps.sub.s4 <- prune_samples(sample_sums(ps.sub.s4) > 0, ps.sub.s4)
ps.sub.s4 <- prune_taxa(taxa_sums(ps.sub.s4) > 0, ps.sub.s4)

cat("\n================ S4 SAMPLE STRUCTURE ================\n")
print(
  data.frame(sample_data(ps.sub.s4)) %>%
    count(genotype, life_stage, generation, vial)
)

cat("\nNumber of samples:", nsamples(ps.sub.s4), "\n")
cat("Number of taxa:", ntaxa(ps.sub.s4), "\n")

# -------------------------
# B3) Transform to relative abundance and collapse to genus
# -------------------------

ps.sub.s4.ra <- transform_sample_counts(
  ps.sub.s4,
  function(x) if (sum(x) > 0) x / sum(x) else x
)

ps.sub.s4.ra <- prune_taxa(taxa_sums(ps.sub.s4.ra) > 0, ps.sub.s4.ra)

ps.genus.s4 <- tax_glom(ps.sub.s4.ra, taxrank = "Genus", NArm = TRUE)

# -------------------------
# B4) Melt and clean metadata
# -------------------------

df_long_s4 <- psmelt(ps.genus.s4) %>%
  mutate(
    Genus = as.character(Genus),
    Family = as.character(Family),
    Genus = ifelse(is.na(Genus) | Genus == "", Family, Genus),

    Genus = ifelse(
      Genus == "Burkholderia-Caballeronia-Paraburkholderia",
      "BCP",
      Genus
    ),

    generation = recode(
      as.character(generation),
      "Generation_0" = "Gen0",
      "Generation_5" = "Gen5"
    ),
    generation = factor(generation, levels = c("Gen0", "Gen5")),

    genotype = factor(
      genotype,
      levels = c("oreore", "simore", "oreaut", "simaut")
    ),

    genotype_label = recode(
      as.character(genotype),
      "oreore" = "(ore);OreR",
      "simore" = "(simw501);OreR",
      "oreaut" = "(ore);Aut",
      "simaut" = "(simw501);Aut"
    ),
    genotype_label = factor(
      genotype_label,
      levels = c(
        "(ore);OreR",
        "(simw501);OreR",
        "(ore);Aut",
        "(simw501);Aut"
      )
    ),

    life_stage = factor(
      life_stage,
      levels = c("Larvae", "Female")
    ),

    life_stage_label = recode(
      as.character(life_stage),
      "Larvae" = "Larvae",
      "Female" = "Adult females"
    ),
    life_stage_label = factor(
      life_stage_label,
      levels = c("Larvae", "Adult females")
    ),

    vial = factor(vial),

    Mito = case_when(
      genotype %in% c("oreore", "oreaut") ~ "ore",
      genotype %in% c("simore", "simaut") ~ "simw501",
      TRUE ~ NA_character_
    ),

    Nuclear = case_when(
      genotype %in% c("oreore", "simore") ~ "OreR",
      genotype %in% c("oreaut", "simaut") ~ "Aut",
      TRUE ~ NA_character_
    ),

    Mito = factor(Mito, levels = c("ore", "simw501")),
    Nuclear = factor(Nuclear, levels = c("OreR", "Aut"))
  )

# -------------------------
# B5) Check sample/vial counts
# -------------------------

sample_count_table_s4 <- df_long_s4 %>%
  distinct(genotype, genotype_label, life_stage_label, generation, vial) %>%
  count(genotype, genotype_label, life_stage_label, generation, name = "n_vials") %>%
  arrange(life_stage_label, genotype, generation)

cat("\n================ S4 SAMPLE COUNTS ================\n")
print(as.data.frame(sample_count_table_s4), row.names = FALSE)

write_csv(
  sample_count_table_s4,
  file.path(data_dir, "Supplementary_Figure_S4_sample_counts.csv")
)

# -------------------------
# B6) Collapse to vial-level genus relative abundance
# -------------------------

df_vial_s4 <- df_long_s4 %>%
  group_by(
    Genus,
    genotype,
    genotype_label,
    Mito,
    Nuclear,
    life_stage,
    life_stage_label,
    generation,
    vial
  ) %>%
  summarise(
    Abundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  df_vial_s4,
  file.path(data_dir, "Supplementary_Figure_S4_full_vial_level_genus_relative_abundance.csv")
)

# -------------------------
# B7) Apply genus abundance filter
# -------------------------

retained_genera <- df_vial_s4 %>%
  group_by(genotype, life_stage_label, generation, Genus) %>%
  summarise(
    n_vials_ge_threshold = sum(Abundance >= GENUS_MIN_ABUND, na.rm = TRUE),
    max_abundance = max(Abundance, na.rm = TRUE),
    mean_abundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Genus) %>%
  summarise(
    keep = any(n_vials_ge_threshold >= MIN_VIALS),
    max_abundance_any_group = max(max_abundance, na.rm = TRUE),
    mean_abundance_all_groups = mean(mean_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(keep) %>%
  arrange(desc(mean_abundance_all_groups)) %>%
  pull(Genus)

cat("\n================ RETAINED GENERA ================\n")
print(retained_genera)

write_csv(
  data.frame(retained_genera = retained_genera),
  file.path(data_dir, "Supplementary_Figure_S4_retained_genera.csv")
)

# -------------------------
# B8) Select top retained genera for plotting
# -------------------------

top_genera <- df_vial_s4 %>%
  filter(Genus %in% retained_genera) %>%
  group_by(Genus) %>%
  summarise(
    mean_abundance = mean(Abundance, na.rm = TRUE),
    max_abundance = max(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = topN) %>%
  pull(Genus)

cat("\n================ TOP GENERA PLOTTED ================\n")
print(top_genera)

write_csv(
  data.frame(top_genera = top_genera),
  file.path(data_dir, "Supplementary_Figure_S4_top10_genera.csv")
)

# -------------------------
# B9) Compute delta for all retained genera
# -------------------------

# Keep every retained genus distinct; restrict to top genera only for plotting.
df_delta_all <- df_vial_s4 %>%
  filter(Genus %in% retained_genera) %>%
  select(
    Genus,
    genotype,
    genotype_label,
    Mito,
    Nuclear,
    life_stage,
    life_stage_label,
    generation,
    vial,
    Abundance
  ) %>%
  pivot_wider(
    names_from = generation,
    values_from = Abundance
  ) %>%
  drop_na(Gen0, Gen5) %>%
  mutate(
    delta = Gen5 - Gen0,
    Mito_label = recode(
      as.character(Mito),
      "ore" = "(ore)",
      "simw501" = "(simw501)"
    ),
    Mito_label = factor(
      Mito_label,
      levels = c("(ore)", "(simw501)")
    ),
    Genus = factor(Genus, levels = retained_genera)
  ) %>%
  arrange(Genus, life_stage_label, Nuclear, genotype, vial)

cat("\n================ DELTA TABLE: ALL RETAINED GENERA ================\n")
print(as.data.frame(df_delta_all), row.names = FALSE)

write_csv(
  df_delta_all,
  file.path(data_dir, "Supplementary_Figure_S4_delta_raw_all_retained.csv")
)

# -------------------------
# B10) Delta table for plotted top genera
# -------------------------

df_delta_top <- df_delta_all %>%
  filter(as.character(Genus) %in% top_genera) %>%
  mutate(
    Genus = factor(as.character(Genus), levels = top_genera)
  )

write_csv(
  df_delta_top,
  file.path(data_dir, "Supplementary_Figure_S4_delta_top10_plotted.csv")
)

# -------------------------
# B11) Per-genus statistics for all retained genera
# -------------------------

stats_all <- df_delta_all %>%
  mutate(
    Genus = as.character(Genus)
  ) %>%
  group_by(Genus) %>%
  group_modify(~{

    fit <- tryCatch(
      lm(delta ~ Mito * Nuclear * life_stage, data = .x),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      return(data.frame(
        p_3way = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    dr <- tryCatch(
      drop1(fit, test = "F"),
      error = function(e) NULL
    )

    if (is.null(dr)) {
      return(data.frame(
        p_3way = NA_real_,
        stringsAsFactors = FALSE
      ))
    }

    data.frame(
      p_3way = if ("Mito:Nuclear:life_stage" %in% rownames(dr)) {
        dr["Mito:Nuclear:life_stage", "Pr(>F)"]
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }) %>%
  ungroup() %>%
  mutate(
    padj_3way = p.adjust(p_3way, method = "BH"),
    label = case_when(
      is.na(padj_3way) ~ "NA",
      padj_3way < 0.001 ~ "BH adj. P < 0.001",
      padj_3way < 0.01 ~ "BH adj. P < 0.01",
      padj_3way < 0.05 ~ "BH adj. P < 0.05",
      TRUE ~ paste0("BH adj. P = ", signif(padj_3way, 2))
    )
  ) %>%
  arrange(padj_3way)

cat("\n================ ALL RETAINED GENERA MODEL STATS ================\n")
print(as.data.frame(stats_all), row.names = FALSE)

write_csv(
  stats_all,
  file.path(stat_dir, "Supplementary_Figure_S4_stats_all_retained.csv")
)

# -------------------------
# B12) Direction and magnitude summary for plotted genera
# -------------------------

direction_summary_top <- df_delta_top %>%
  group_by(Genus, genotype, genotype_label, Mito, Nuclear, life_stage_label) %>%
  summarise(
    mean_delta = mean(delta, na.rm = TRUE),
    median_delta = median(delta, na.rm = TRUE),
    sd_delta = sd(delta, na.rm = TRUE),
    consistency = mean(sign(delta) == sign(mean(delta, na.rm = TRUE))),
    direction = case_when(
      mean_delta > 0 ~ "Increase",
      mean_delta < 0 ~ "Decrease",
      TRUE ~ "Neutral"
    ),
    .groups = "drop"
  ) %>%
  mutate(
    abs_mean_delta = abs(mean_delta),
    effect_strength = abs_mean_delta * consistency
  ) %>%
  arrange(Genus, life_stage_label, Nuclear, Mito)

cat("\n================ DIRECTION SUMMARY: TOP GENERA ================\n")
print(as.data.frame(direction_summary_top), row.names = FALSE)

write_csv(
  direction_summary_top,
  file.path(data_dir, "Supplementary_Figure_S4_direction_summary_top10.csv")
)

# -------------------------
# B13) Supplementary Figure S4 plot
# -------------------------

y_lim_s4 <- max(abs(df_delta_top$delta), na.rm = TRUE) * 1.15

S4_delta_screen <- ggplot(
  df_delta_top,
  aes(
    x = Mito_label,
    y = delta,
    fill = Mito_label,
    color = Mito_label
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.42,
    color = "black"
  ) +
  geom_point(
    aes(shape = Mito_label),
    position = position_jitter(width = 0.08, height = 0),
    size = 2.05,
    alpha = 0.88,
    stroke = 0.55
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 2.75,
    fill = "white",
    color = "black",
    stroke = 0.65,
    show.legend = FALSE
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.15,
    linewidth = 0.55,
    color = "black"
  ) +
  facet_grid(
    life_stage_label + Nuclear ~ Genus,
    scales = "fixed"
  ) +
  scale_fill_manual(
    values = mito_fill_cols,
    breaks = c("(ore)", "(simw501)"),
    labels = c(
      "(ore)" = "(ore)",
      "(simw501)" = expression(simw^501)
    ),
    name = "Mitochondrial genotype"
  ) +
  scale_color_manual(
    values = mito_line_cols,
    breaks = c("(ore)", "(simw501)"),
    labels = c(
      "(ore)" = "(ore)",
      "(simw501)" = expression(simw^501)
    ),
    name = "Mitochondrial genotype"
  ) +
  scale_shape_manual(
    values = c(
      "(ore)" = 21,
      "(simw501)" = 24
    ),
    breaks = c("(ore)", "(simw501)"),
    labels = c(
      "(ore)" = "(ore)",
      "(simw501)" = expression(simw^501)
    ),
    name = "Mitochondrial genotype"
  ) +
  scale_x_discrete(
    labels = c(
      "(ore)" = "(ore)",
      "(simw501)" = expression(simw^501)
    )
  ) +
  scale_y_continuous(
    limits = c(-y_lim_s4, y_lim_s4),
    breaks = scales::pretty_breaks(n = 5),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  labs(
    x = "Mitochondrial genotype",
    y = expression(Delta * " relative abundance (Generation 5 - Generation 0)")
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = c(21, 24),
        color = mito_line_cols,
        fill = mito_fill_cols,
        size = 3
      )
    ),
    color = "none",
    shape = "none"
  ) +
  theme_classic(base_size = 10.8) +
  theme(
    axis.title = element_text(
      face = "bold",
      color = "black",
      size = 11.2
    ),
    axis.title.x = element_text(
      margin = margin(t = 8)
    ),
    axis.title.y = element_text(
      margin = margin(r = 8)
    ),
    axis.text = element_text(
      color = "black",
      size = 8.6
    ),
    axis.text.x = element_text(
      face = "bold",
      angle = 28,
      hjust = 1
    ),
    axis.line = element_blank(),
    axis.ticks = element_line(
      color = "black",
      linewidth = 0.35
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.45
    ),
    strip.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.45
    ),
    strip.text = element_text(
      face = "bold",
      size = 8.8,
      color = "black"
    ),
    panel.spacing.x = unit(0.45, "lines"),
    panel.spacing.y = unit(0.45, "lines"),
    legend.position = "top",
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    legend.text = element_text(
      size = 9.5
    ),
    legend.key.size = unit(0.50, "cm"),
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
      r = 10,
      b = 10,
      l = 10
    )
  )

S4_delta_screen

ggsave(
  filename = file.path(fig_dir, "Supplementary_Figure_S4_full_genus_delta_screen.pdf"),
  plot = S4_delta_screen,
  width = 14.5,
  height = 8.2,
  units = "in",
  device = PDF_DEVICE,
  bg = "white"
)

ggsave(
  filename = file.path(fig_dir, "Supplementary_Figure_S4_full_genus_delta_screen.png"),
  plot = S4_delta_screen,
  width = 14.5,
  height = 8.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(fig_dir, "Supplementary_Figure_S4_full_genus_delta_screen.jpeg"),
  plot = S4_delta_screen,
  width = 14.5,
  height = 8.2,
  units = "in",
  dpi = 600,
  bg = "white",
  device = "jpeg"
)

# -------------------------
# B14) Save S4 summary
# -------------------------

capture.output(
  {
    cat("SUPPLEMENTARY FIGURE S4 FINAL SUMMARY\n")
    cat("=====================================\n\n")

    cat("Genus filter:\n")
    cat("GENUS_MIN_ABUND =", GENUS_MIN_ABUND, "\n")
    cat("MIN_VIALS =", MIN_VIALS, "\n")
    cat("topN =", topN, "\n\n")

    cat("Sample counts:\n")
    print(sample_count_table_s4)

    cat("\nRetained genera:\n")
    print(retained_genera)

    cat("\nTop genera plotted:\n")
    print(top_genera)

    cat("\nModel statistics for retained genera:\n")
    print(stats_all)

    cat("\nDirection summary for plotted genera:\n")
    print(direction_summary_top)
  },
  file = file.path(stat_dir, "Supplementary_Figure_S4_summary.txt")
)

# ============================================================
# PART C: SUPPLEMENTARY TABLE S4
# Figure 5 + Supplementary Figure S4 genus-level delta statistics
# ============================================================

# -----------------------
# C1) Figure 5 two-genera model stats
# -----------------------

supp_s4_fig5_stats <- stats_two %>%
  mutate(
    Analysis = "Figure 5 two-genera linear model",
    Response = "Delta relative abundance: Generation 5 - Generation 0",
    Model = "delta ~ Mito * Nuclear * life_stage",
    Test = "Mito × Nuclear × developmental stage"
  ) %>%
  select(
    Analysis, Response, Model, Test,
    Genus, p_3way, padj_3way, label
  )

# -----------------------
# C2) Figure 5 two-genera direction summary
# -----------------------

supp_s4_fig5_direction <- direction_summary %>%
  mutate(
    Analysis = "Figure 5 two-genera direction summary",
    Response = "Delta relative abundance: Generation 5 - Generation 0"
  ) %>%
  select(
    Analysis, Response,
    Genus, genotype, genotype_label, Mito, Nuclear, life_stage_label,
    mean_delta, median_delta, sd_delta,
    consistency, direction, abs_mean_delta, effect_strength
  )

# -----------------------
# C3) Supplementary Fig. S4 retained genera
# -----------------------

supp_s4_retained_genera <- data.frame(
  Analysis = "Supplementary Fig. S4 retained genera",
  Filter = "Genera retained if >=0.5% relative abundance in at least 2 replicate vials within any genotype × life stage × generation group",
  Genus = retained_genera
)

# -----------------------
# C4) Supplementary Fig. S4 top 10 plotted genera
# -----------------------

supp_s4_top10_genera <- data.frame(
  Analysis = "Supplementary Fig. S4 top 10 plotted genera",
  Selection = "Top retained genera by mean relative abundance",
  Genus = top_genera
)

# -----------------------
# C5) Supplementary Fig. S4 all retained genera stats
# -----------------------

supp_s4_all_retained_stats <- stats_all %>%
  mutate(
    Analysis = "Supplementary Fig. S4 all retained genera linear model",
    Response = "Delta relative abundance: Generation 5 - Generation 0",
    Model = "delta ~ Mito * Nuclear * life_stage",
    Test = "Mito × Nuclear × developmental stage"
  ) %>%
  select(
    Analysis, Response, Model, Test,
    Genus, p_3way, padj_3way, label
  )

# -----------------------
# C6) Supplementary Fig. S4 top 10 direction summary
# -----------------------

supp_s4_top10_direction <- direction_summary_top %>%
  mutate(
    Analysis = "Supplementary Fig. S4 top 10 direction summary",
    Response = "Delta relative abundance: Generation 5 - Generation 0"
  ) %>%
  select(
    Analysis, Response,
    Genus, genotype, genotype_label, Mito, Nuclear, life_stage_label,
    mean_delta, median_delta, sd_delta,
    consistency, direction, abs_mean_delta, effect_strength
  )

# -----------------------
# C7) Write Excel workbook
# -----------------------

wb <- createWorkbook()

addWorksheet(wb, "Fig5_two_genera_stats")
writeData(wb, "Fig5_two_genera_stats", supp_s4_fig5_stats)

addWorksheet(wb, "Fig5_direction_summary")
writeData(wb, "Fig5_direction_summary", supp_s4_fig5_direction)

addWorksheet(wb, "S4_retained_genera")
writeData(wb, "S4_retained_genera", supp_s4_retained_genera)

addWorksheet(wb, "S4_top10_genera")
writeData(wb, "S4_top10_genera", supp_s4_top10_genera)

addWorksheet(wb, "S4_all_retained_stats")
writeData(wb, "S4_all_retained_stats", supp_s4_all_retained_stats)

addWorksheet(wb, "S4_top10_direction")
writeData(wb, "S4_top10_direction", supp_s4_top10_direction)

header_style <- createStyle(textDecoration = "bold", halign = "center")

for (sheet in names(wb)) {
  addStyle(
    wb,
    sheet = sheet,
    style = header_style,
    rows = 1,
    cols = 1:50,
    gridExpand = TRUE
  )
  setColWidths(wb, sheet = sheet, cols = 1:50, widths = "auto")
}

saveWorkbook(
  wb,
  file = file.path(stat_dir, "Supplementary_Table_S4_genus_level_delta_statistics.xlsx"),
  overwrite = TRUE
)

# Also save CSV versions of Supplementary Table S4 sheets
write_csv(
  supp_s4_fig5_stats,
  file.path(stat_dir, "Supplementary_Table_S4A_Fig5_two_genera_stats.csv")
)

write_csv(
  supp_s4_fig5_direction,
  file.path(stat_dir, "Supplementary_Table_S4B_Fig5_direction_summary.csv")
)

write_csv(
  supp_s4_retained_genera,
  file.path(stat_dir, "Supplementary_Table_S4C_S4_retained_genera.csv")
)

write_csv(
  supp_s4_top10_genera,
  file.path(stat_dir, "Supplementary_Table_S4D_S4_top10_genera.csv")
)

write_csv(
  supp_s4_all_retained_stats,
  file.path(stat_dir, "Supplementary_Table_S4E_S4_all_retained_stats.csv")
)

write_csv(
  supp_s4_top10_direction,
  file.path(stat_dir, "Supplementary_Table_S4F_S4_top10_direction.csv")
)

# ============================================================
# DONE
# ============================================================

cat("\nDONE. Files saved inside folder:\n")
cat(main_outdir, "\n\n")

cat("Figure files:\n")
cat(file.path(fig_dir, "Figure5_two_genera_delta_ISME.pdf"), "\n")
cat(file.path(fig_dir, "Figure5_two_genera_delta_ISME.png"), "\n")
cat(file.path(fig_dir, "Figure5_two_genera_delta_ISME.jpeg"), "\n")
cat(file.path(fig_dir, "Supplementary_Figure_S4_full_genus_delta_screen.pdf"), "\n")
cat(file.path(fig_dir, "Supplementary_Figure_S4_full_genus_delta_screen.png"), "\n")
cat(file.path(fig_dir, "Supplementary_Figure_S4_full_genus_delta_screen.jpeg"), "\n\n")

cat("Figure source data:\n")
cat(file.path(data_dir, "Figure5_two_genera_delta_raw.csv"), "\n")
cat(file.path(data_dir, "Figure5_two_genera_direction_summary.csv"), "\n")
cat(file.path(data_dir, "Figure5_two_genera_sample_counts.csv"), "\n")
cat(file.path(data_dir, "Figure5_two_genera_full_vial_level_genus_relative_abundance.csv"), "\n")
cat(file.path(data_dir, "Supplementary_Figure_S4_sample_counts.csv"), "\n")
cat(file.path(data_dir, "Supplementary_Figure_S4_full_vial_level_genus_relative_abundance.csv"), "\n")
cat(file.path(data_dir, "Supplementary_Figure_S4_delta_raw_all_retained.csv"), "\n")
cat(file.path(data_dir, "Supplementary_Figure_S4_delta_top10_plotted.csv"), "\n")
cat(file.path(data_dir, "Supplementary_Figure_S4_retained_genera.csv"), "\n")
cat(file.path(data_dir, "Supplementary_Figure_S4_top10_genera.csv"), "\n")
cat(file.path(data_dir, "Supplementary_Figure_S4_direction_summary_top10.csv"), "\n\n")

cat("Statistics:\n")
cat(file.path(stat_dir, "Supplementary_Table_S4_genus_level_delta_statistics.xlsx"), "\n")
cat(file.path(stat_dir, "Figure5_two_genera_stats.csv"), "\n")
cat(file.path(stat_dir, "Figure5_two_genera_summary.txt"), "\n")
cat(file.path(stat_dir, "Supplementary_Figure_S4_stats_all_retained.csv"), "\n")
cat(file.path(stat_dir, "Supplementary_Figure_S4_summary.txt"), "\n\n")

cat("Processed input copy:\n")
cat(file.path(input_dir, "06_phyloseq_clean_CHAP1_NOHOST.rds"), "\n")

# Record software versions used for this run.
capture.output(sessionInfo(), file = file.path(main_outdir, "sessionInfo.txt"))
