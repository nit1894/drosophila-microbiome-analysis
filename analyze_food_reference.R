#!/usr/bin/env Rscript
# Microbiome food-reference analysis
# Adapted from Nitin Bansal's research analysis for a portable code example.
# Run from the repository root; see README.md for inputs and assumptions.

parse_options <- function(args) {
  opts <- list(input = NULL, output = NULL, demo = FALSE,
               overwrite = FALSE, exclude_genus = character(0))
  i <- 1L
  while (i <= length(args)) {
    arg <- args[i]
    if (arg %in% c('--demo', '--overwrite')) {
      opts[[substring(arg, 3L)]] <- TRUE
    } else if (arg %in% c('--input', '--output', '--exclude-genus')) {
      if (i == length(args) || startsWith(args[i + 1L], '--')) {
        stop('Missing value after ', arg, call. = FALSE)
      }
      i <- i + 1L
      key <- gsub('-', '_', substring(arg, 3L))
      opts[[key]] <- if (key == 'exclude_genus') {
        trimws(strsplit(args[i], ',', fixed = TRUE)[[1]])
      } else args[i]
    } else {
      stop('Unknown option: ', arg, call. = FALSE)
    }
    i <- i + 1L
  }
  if (opts$demo && !is.null(opts$input)) {
    stop('Choose either --demo or --input.', call. = FALSE)
  }
  if (!opts$demo && is.null(opts$input)) {
    stop('Use --demo, or --input PATH_TO_PHYLOSEQ.rds. See README.md.', call. = FALSE)
  }
  opts
}

opts <- parse_options(commandArgs(trailingOnly = TRUE))
if (getRversion() < '4.1.0') stop('R 4.1 or later is required.', call. = FALSE)
file_arg <- grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub('^--file=', '', file_arg[1]), mustWork = TRUE)
} else NULL
project_root <- if (!is.null(script_path)) dirname(dirname(script_path)) else getwd()
main_outdir <- if (is.null(opts$output)) {
  file.path(project_root, 'outputs', if (opts$demo) 'demo' else 'research')
} else opts$output
if (dir.exists(main_outdir) &&
    length(list.files(main_outdir, all.files = TRUE, no.. = TRUE)) &&
    !opts$overwrite) {
  stop('Output directory is not empty. Choose another --output directory, ',
       'or explicitly use --overwrite.', call. = FALSE)
}
excluded_genera <- opts$exclude_genus
required_packages <- c('phyloseq', 'dplyr', 'tidyr', 'tibble', 'vegan',
                       'ggplot2', 'stringr', 'scales', 'openxlsx', 'emmeans', 'car')
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop('Missing packages: ', paste(missing_packages, collapse = ', '),
       '. Run Rscript scripts/install_packages.R first.', call. = FALSE)
}
if (utils::packageVersion('ggplot2') < '3.4.0') {
  stop('ggplot2 >= 3.4.0 is required for linewidth support.', call. = FALSE)
}
pdf_device <- if (capabilities('cairo')) grDevices::cairo_pdf else grDevices::pdf

load_demo <- function() {
  demo_dir <- file.path(project_root, 'data', 'demo')
  counts <- as.matrix(read.csv(file.path(demo_dir, 'counts.csv'),
                               row.names = 1, check.names = FALSE))
  taxonomy <- as.matrix(read.csv(file.path(demo_dir, 'taxonomy.csv'),
                                 row.names = 1, check.names = FALSE))
  metadata <- read.csv(file.path(demo_dir, 'metadata.csv'),
                       row.names = 1, check.names = FALSE,
                       stringsAsFactors = FALSE)
  phyloseq::phyloseq(
    phyloseq::otu_table(counts, taxa_are_rows = TRUE),
    phyloseq::tax_table(taxonomy), phyloseq::sample_data(metadata)
  )
}

suppressPackageStartupMessages({
  library(phyloseq)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(vegan)
  library(ggplot2)
  library(stringr)
  library(scales)
  library(openxlsx)
})

# Statistical dependencies are required so package availability cannot
# silently change the method or omit an analysis.
library(emmeans)


set.seed(1234)

# -----------------------
# 2) Create output folders
# -----------------------


dir.create(main_outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(main_outdir, "code"), showWarnings = FALSE)
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

ps.clean <- if (opts$demo) load_demo() else readRDS(opts$input)
if (!methods::is(ps.clean, 'phyloseq')) {
  stop('Input must be a phyloseq object.', call. = FALSE)
}
methods::validObject(ps.clean)
if (is.null(tax_table(ps.clean, errorIfNULL = FALSE)) ||
    !all(c('Genus', 'Family') %in% rank_names(ps.clean))) {
  stop('Taxonomy must include Genus and Family columns.', call. = FALSE)
}
sample_meta <- data.frame(sample_data(ps.clean))
required_columns <- c('genotype', 'generation', 'life_stage', 'vial', 'replicate')
if (!all(required_columns %in% colnames(sample_meta))) {
  stop('Missing metadata columns: ',
       paste(setdiff(required_columns, colnames(sample_meta)), collapse = ', '),
       call. = FALSE)
}
raw_counts <- as(otu_table(ps.clean), 'matrix')
if (!is.numeric(raw_counts) || any(!is.finite(raw_counts)) || any(raw_counts < 0)) {
  stop('Counts must be numeric, finite, and nonnegative.', call. = FALSE)
}
if (any(abs(raw_counts - round(raw_counts)) > 1e-7)) {
  stop('Supply untransformed counts rather than relative abundances.', call. = FALSE)
}
depths <- sample_sums(ps.clean)
qc <- data.frame(sample_id = names(depths), sequencing_depth = unname(depths),
                 retained = depths > 0)
write.csv(qc, file.path(stat_dir, 'sample_depth_qc.csv'), row.names = FALSE)
if (!any(depths > 0)) stop('Every sample has zero counts.', call. = FALSE)
if (any(depths == 0)) message('Removing ', sum(depths == 0), ' zero-count samples.')
# Remove empty samples BEFORE division by sample totals.
ps.clean <- prune_samples(depths > 0, ps.clean)
ps.clean <- prune_taxa(taxa_sums(ps.clean) > 0, ps.clean)
sample_meta <- data.frame(sample_data(ps.clean))
selected <- with(sample_meta, genotype %in% c('oreore', 'simore', 'oreaut', 'simaut') &
                   (life_stage == 'food' |
                    (generation == 'Generation_0' & life_stage %in% c('Larvae', 'Female'))))
if (anyNA(selected)) stop('Selected-sample metadata contain missing values.', call. = FALSE)
if (!any(selected)) stop('No matching food or G0 host samples.', call. = FALSE)
for (field in required_columns) {
  values <- as.character(sample_meta[selected, field])
  if (anyNA(values) || any(trimws(values) == '')) {
    stop('Missing or blank ', field, ' in selected samples.', call. = FALSE)
  }
}
write.csv(sample_meta[selected, , drop = FALSE],
          file.path(input_dir, 'selected_sample_metadata.csv'), row.names = TRUE)
capture.output(str(opts), file = file.path(input_dir, 'run_options.txt'))
if (length(excluded_genera)) {
  message('Composition plot excludes: ', paste(excluded_genera, collapse = ', '),
          '. Bray-Curtis distances still use all retained ASVs.')
}

saveRDS(
  ps.clean,
  file = file.path(input_dir, "06_phyloseq_clean_CHAP1_NOHOST.rds")
)

# ============================================================
# PART A: SUPPLEMENTARY FIGURE S1
# Food/source community composition only
# ============================================================

fix_genus <- function(df) {
  df %>%
    mutate(
      Genus  = as.character(Genus),
      Family = as.character(Family),
      Genus  = ifelse(is.na(Genus) | Genus == "",
                      ifelse(is.na(Family) | Family == "", "Unclassified",
                             paste0("Unclassified: ", Family)), Genus)
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

genotype_facet_labels <- c(
  "(ore);OreR"       = "'(ore);OreR'",
  "(simw501);OreR"   = "'(simw'^501*');OreR'",
  "(ore);Aut"        = "'(ore);Aut'",
  "(simw501);Aut"    = "'(simw'^501*');Aut'"
)

ps.food <- subset_samples(
  ps.clean,
  life_stage == "food" &
    genotype %in% c("oreore", "simore", "oreaut", "simaut")
)

if (nsamples(ps.food) == 0) stop("No food samples found.", call. = FALSE)
ps.food <- prune_samples(sample_sums(ps.food) > 0, ps.food)
ps.food <- prune_taxa(taxa_sums(ps.food) > 0, ps.food)

cat("\nSamples included in Supplementary Fig. S1:\n")
print(
  data.frame(sample_data(ps.food)) %>%
    count(genotype, life_stage, vial)
)

cat("\nNumber of food/source samples:", nsamples(ps.food), "\n")
cat("Number of taxa:", ntaxa(ps.food), "\n")

ps.food.ra <- transform_sample_counts(ps.food, function(x) x / sum(x))
ps.food.ra <- prune_taxa(taxa_sums(ps.food.ra) > 0, ps.food.ra)

ps.genus <- tax_glom(ps.food.ra, taxrank = "Genus", NArm = FALSE)

df <- psmelt(ps.genus) %>%
  fix_genus() %>%
  filter(!Genus %in% excluded_genera) %>%
  mutate(
    Genus = ifelse(
      Genus == "Burkholderia-Caballeronia-Paraburkholderia",
      "BCP",
      Genus
    ),
    genotype_pretty = pretty_genotype(genotype),
    genotype_pretty = factor(
      genotype_pretty,
      levels = c(
        "(ore);OreR",
        "(simw501);OreR",
        "(ore);Aut",
        "(simw501);Aut"
      )
    ),
    sample_short = case_when(
      !is.na(vial) & vial != "" ~ str_replace(as.character(vial), "Vial_", "V"),
      TRUE ~ as.character(Sample)
    )
  )

if (nrow(df) == 0L) stop('No taxa remain for the composition plot.', call. = FALSE)

topN <- 10

top_genera <- df %>%
  group_by(Genus) %>%
  summarise(total_abundance = sum(Abundance), .groups = "drop") %>%
  arrange(desc(total_abundance)) %>%
  slice_head(n = topN) %>%
  pull(Genus)

cat("\nTop genera in food/source samples:\n")
print(top_genera)

df.sum <- df %>%
  mutate(
    Genus2 = ifelse(Genus %in% top_genera, Genus, "Other")
  ) %>%
  group_by(genotype_pretty, sample_short, Genus2) %>%
  summarise(
    Relative_abundance = sum(Abundance),
    .groups = "drop"
  )

df.sum <- df.sum %>%
  group_by(genotype_pretty, sample_short) %>%
  mutate(display_total = sum(Relative_abundance))
if (any(!is.finite(df.sum$display_total)) || any(df.sum$display_total <= 0)) {
  stop('A source sample has no retained abundance for the composition plot.', call. = FALSE)
}
df.sum <- df.sum %>%
  mutate(
    Relative_abundance = Relative_abundance / display_total
  ) %>%
  ungroup() %>% select(-display_total)

cat("\nCheck: each food/source sample should sum to 1:\n")
print(
  df.sum %>%
    group_by(genotype_pretty, sample_short) %>%
    summarise(total = sum(Relative_abundance), .groups = "drop")
)

preferred_order <- c(
  "Acetobacter", "Lactococcus", "Levilactobacillus", "Lactiplantibacillus",
  "BCP", "Limosilactobacillus", "Liquorilactobacillus", "Enterobacter",
  "Pseudoclavibacter", "Prevotella_9", "Other"
)

extra_genera <- setdiff(unique(as.character(df.sum$Genus2)), preferred_order)

legend_order <- c(
  preferred_order[preferred_order %in% unique(as.character(df.sum$Genus2))],
  extra_genera
)

legend_order <- c(setdiff(legend_order, "Other"), "Other")
legend_order <- legend_order[legend_order %in% unique(as.character(df.sum$Genus2))]

df.sum$Genus2 <- factor(df.sum$Genus2, levels = legend_order)

genus_palette <- c(
  "Acetobacter" = "#C23B22", "Lactococcus" = "#6A4C93",
  "Levilactobacillus" = "#F4D03F", "Lactiplantibacillus" = "#3A86FF",
  "BCP" = "black", "Limosilactobacillus" = "#AA4499",
  "Liquorilactobacillus" = "#577590", "Enterobacter" = "#BC6C25",
  "Pseudoclavibacter" = "#F28482", "Prevotella_9" = "#90BE6D",
  "Other" = "grey85"
)

missing_cols <- setdiff(legend_order, names(genus_palette))
if (length(missing_cols) > 0) {
  extra_cols <- scales::hue_pal()(length(missing_cols))
  names(extra_cols) <- missing_cols
  genus_palette <- c(genus_palette, extra_cols)
}

S1_food <- ggplot(
  df.sum,
  aes(x = sample_short, y = Relative_abundance, fill = Genus2)
) +
  geom_col(width = 0.82, color = "white", linewidth = 0.25) +
  facet_grid(
    ~ genotype_pretty, scales = "free_x", space = "free_x",
    labeller = labeller(
      genotype_pretty = as_labeller(genotype_facet_labels, label_parsed)
    )
  ) +
  scale_fill_manual(
    values = genus_palette, breaks = legend_order, name = "Genus",
    guide = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(linewidth = 0))
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "Food/source sample", y = "Relative abundance",
       caption = if (opts$demo) "Synthetic demonstration data" else NULL) +
  theme_classic(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
    strip.text = element_text(face = "bold", size = 12.5, color = "black"),
    axis.title = element_text(face = "bold", color = "black", size = 14),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.text = element_text(color = "black", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11.5),
    legend.text = element_text(size = 10, color = "black"),
    legend.key.size = unit(0.42, "cm"),
    legend.spacing.x = unit(0.22, "cm"),
    panel.grid = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.65),
    axis.line = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 10, r = 12, b = 10, l = 12)
  )

S1_food

ggsave(file.path(fig_dir, "Supplementary_Figure_S1_food_source_composition.pdf"),
       S1_food, width = 11.5, height = 5.8, units = "in", device = pdf_device, bg = "white")
ggsave(file.path(fig_dir, "Supplementary_Figure_S1_food_source_composition.png"),
       S1_food, width = 11.5, height = 5.8, units = "in", dpi = 600, bg = "white")
ggsave(file.path(fig_dir, "Supplementary_Figure_S1_food_source_composition.jpeg"),
       S1_food, width = 11.5, height = 5.8, units = "in", dpi = 600, bg = "white", device = "jpeg")

write.csv(df.sum, file.path(data_dir, "Supplementary_Figure_S1_food_source_composition_data.csv"), row.names = FALSE)

# ============================================================
# PART B: FIGURE 1B FINAL
# Generation 0 distance from food/source community
# ============================================================

ps.ra <- transform_sample_counts(ps.clean, function(x) x / sum(x))
ps.ra <- prune_samples(sample_sums(ps.ra) > 0, ps.ra)

otu <- as(otu_table(ps.ra), "matrix")
if (taxa_are_rows(ps.ra)) otu <- t(otu)

meta <- data.frame(sample_data(ps.ra)) %>%
  rownames_to_column("sample_id_phyloseq") %>%
  mutate(
    genotype   = as.character(genotype),
    generation = as.character(generation),
    life_stage = as.character(life_stage),
    vial       = as.character(vial),
    replicate  = as.character(replicate)
  )

otu <- otu[meta$sample_id_phyloseq, , drop = FALSE]

# -----------------------
# B2) Define food/source and G0 host samples
# -----------------------

meta_fig <- meta %>%
  mutate(
    sample_type = case_when(
      life_stage == "food" & genotype != "Homogenized" ~ "Food/source",
      generation == "Generation_0" & life_stage == "Larvae" ~ "G0 larvae",
      generation == "Generation_0" & life_stage == "Female" ~ "G0 adult females",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(sample_type), genotype %in% c("oreore", "simore", "oreaut", "simaut")) %>%
  mutate(
    genotype = factor(
      genotype,
      levels = c("oreore", "simore", "oreaut", "simaut")
    ),
    genotype_label = dplyr::recode(
      as.character(genotype),
      "oreore" = "(ore);OreR",
      "simore" = "(simw501);OreR",
      "oreaut" = "(ore);Aut",
      "simaut" = "(simw501);Aut"
    ),
    genotype_label = factor(
      genotype_label,
      levels = c("(ore);OreR", "(simw501);OreR", "(ore);Aut", "(simw501);Aut")
    ),
    stage_label = case_when(
      sample_type == "G0 larvae" ~ "Larvae",
      sample_type == "G0 adult females" ~ "Adult females",
      TRUE ~ "Food/source"
    ),
    stage_label = factor(
      stage_label,
      levels = c("Larvae", "Adult females", "Food/source")
    ),
    
    # Mito = which mitochondrial haplotype (ore vs simw501)
    # Nuclear = which nuclear background (OreR vs Aut)
    # This lets us test Mito, Nuclear, and Mito:Nuclear separately,
    # in addition to (not instead of) the original 4-level genotype test.
    Mito = case_when(
      genotype %in% c("oreore", "oreaut") ~ "(ore)",
      genotype %in% c("simore", "simaut") ~ "(simw501)",
      TRUE ~ NA_character_
    ),
    Mito = factor(Mito, levels = c("(ore)", "(simw501)")),
    Nuclear = case_when(
      genotype %in% c("oreore", "simore") ~ "OreR",
      genotype %in% c("oreaut", "simaut") ~ "Aut",
      TRUE ~ NA_character_
    ),
    Nuclear = factor(Nuclear, levels = c("OreR", "Aut"))
  )

otu_fig <- otu[meta_fig$sample_id_phyloseq, , drop = FALSE]

# -----------------------
# B3) Collapse host samples to vial-level means
# -----------------------

meta_fig <- meta_fig %>%
  mutate(
    vial_group = case_when(
      sample_type == "Food/source" ~ paste("Food", genotype, sample_id_phyloseq, sep = "__"),
      sample_type != "Food/source" ~ paste(genotype, generation, life_stage, vial, sep = "__")
    )
  )

otu_vial <- rowsum(otu_fig, group = meta_fig$vial_group, reorder = FALSE)

n_per_group <- as.vector(table(factor(meta_fig$vial_group, levels = rownames(otu_vial))))
otu_vial <- sweep(otu_vial, 1, n_per_group, "/")
otu_vial <- sweep(otu_vial, 1, rowSums(otu_vial), "/")
otu_vial[is.na(otu_vial)] <- 0

meta_vial <- meta_fig %>%
  group_by(vial_group) %>%
  summarise(
    sample_type = first(sample_type),
    genotype = first(genotype),
    genotype_label = first(genotype_label),
    Mito = first(Mito),
    Nuclear = first(Nuclear),
    generation = first(generation),
    life_stage = first(life_stage),
    stage_label = first(stage_label),
    vial = first(vial),
    n_pools_collapsed = n(),
    .groups = "drop"
  ) %>%
  as.data.frame()

rownames(meta_vial) <- meta_vial$vial_group
otu_vial <- otu_vial[meta_vial$vial_group, , drop = FALSE]

cat("\nSample counts after vial-level collapse:\n")
print(table(meta_vial$sample_type))

cat("\nSample counts by genotype and stage:\n")
print(table(meta_vial$genotype_label, meta_vial$stage_label))

# -----------------------
# B5) Food/source consistency check
# -----------------------

food_ids <- rownames(meta_vial)[meta_vial$sample_type == "Food/source"]
food_meta <- meta_vial[food_ids, , drop = FALSE]
food_otu <- otu_vial[food_ids, , drop = FALSE]

food_n <- table(droplevels(food_meta$genotype))
if (length(food_n) != 4L || any(food_n < 2L)) {
  stop('The source comparison requires four genotypes and at least two food samples each.',
       call. = FALSE)
}
if (any(duplicated(paste(food_meta$genotype, food_meta$vial, sep = '__')))) {
  warning('Food samples share genotype/vial labels. Confirm biological independence; ',
          'food samples are not collapsed in the source PERMANOVA.')
}
write.csv(as.data.frame(food_n), file.path(stat_dir, 'food_replicate_counts.csv'),
          row.names = FALSE)
food_bray <- vegdist(food_otu, method = "bray")

cat("\nFood/source PERMANOVA: food community ~ genotype\n")
food_perm <- adonis2(food_bray ~ genotype, data = food_meta, permutations = 999)
print(food_perm)



# With fewer than three observations in any group, within-group
# dispersion is not supported adequately; skip rather than export a misleading P.
food_bd_perm <- NULL
food_dispersion_status <- if (any(food_n < 3L)) {
  'Skipped: fewer than three food observations in at least one genotype.'
} else {
  food_bd <- betadisper(food_bray, group = droplevels(food_meta$genotype), type = 'median')
  food_bd_perm <- permutest(food_bd, permutations = 999)
  'Calculated: spatial-median PERMDISP, 999 permutations.'
}
cat('\nFood/source dispersion:', food_dispersion_status, '\n')
if (!is.null(food_bd_perm)) print(food_bd_perm)
writeLines(food_dispersion_status, file.path(stat_dir, 'food_dispersion_status.txt'))

# -----------------------
# B6) Calculate distance to shared food/source centroid
# -----------------------

food_centroid <- colMeans(food_otu)

host_ids <- rownames(meta_vial)[meta_vial$sample_type %in% c("G0 larvae", "G0 adult females")]

dist_to_food <- sapply(host_ids, function(id) {
  as.numeric(vegdist(rbind(otu_vial[id, ], food_centroid), method = "bray")[1])
})


dist_df <- meta_vial[host_ids, , drop = FALSE] %>%
  rownames_to_column("sample_id") %>%
  mutate(
    distance_to_food = dist_to_food[sample_id],
    stage_label = factor(stage_label, levels = c("Larvae", "Adult females")),
    genotype_label = factor(
      genotype_label,
      levels = c("(ore);OreR", "(simw501);OreR", "(ore);Aut", "(simw501);Aut")
    ),
    Mito = factor(Mito, levels = c("(ore)", "(simw501)")),
    Nuclear = factor(Nuclear, levels = c("OreR", "Aut"))
  )

# -----------------------
# B7) Statistical models
# -----------------------

cat("\nMain model (GENOTYPE): distance_to_food ~ stage_label * genotype_label\n")
host_cells <- table(dist_df$stage_label, dist_df$genotype_label)
if (any(host_cells < 2L)) {
  stop('Each stage/genotype cell requires at least two independent host vials.', call. = FALSE)
}
if (any(!is.finite(dist_df$distance_to_food))) {
  stop('Non-finite reference distances.', call. = FALSE)
}
lm_dist <- lm(distance_to_food ~ stage_label * genotype_label, data = dist_df)
lm_dist_anova <- anova(lm_dist)

print(lm_dist_anova)
print(summary(lm_dist))

# Two explicit contrast families. BH is applied across all contrasts
# within each exported table, including contrasts across the by groups.
emm_stage <- emmeans(lm_dist, ~ stage_label | genotype_label)
emm_genotype <- emmeans(lm_dist, ~ genotype_label | stage_label)
emm_stage_by_genotype <- as.data.frame(summary(pairs(emm_stage), adjust = 'none'))
emm_genotype_by_stage <- as.data.frame(summary(pairs(emm_genotype), adjust = 'none'))
emm_stage_by_genotype$p_BH <- p.adjust(emm_stage_by_genotype$p.value, method = 'BH')
emm_genotype_by_stage$p_BH <- p.adjust(emm_genotype_by_stage$p.value, method = 'BH')
print(emm_stage_by_genotype)
print(emm_genotype_by_stage)
write.csv(emm_stage_by_genotype, file.path(stat_dir, 'stage_contrasts_BH.csv'), row.names = FALSE)
write.csv(emm_genotype_by_stage, file.path(stat_dir, 'genotype_contrasts_BH.csv'), row.names = FALSE)

p_stage <- lm_dist_anova["stage_label", "Pr(>F)"]
p_genotype <- lm_dist_anova["genotype_label", "Pr(>F)"]
p_interaction <- lm_dist_anova["stage_label:genotype_label", "Pr(>F)"]

format_p <- function(p) {
  if (is.na(p)) return("P = NA")
  if (p < 0.001) return("P < 0.001")
  paste0("P = ", signif(p, 2))
}

cat("\nP-values for reporting (GENOTYPE genotype-level model):\n")
cat("Stage:", format_p(p_stage), "\n")
cat("Genotype:", format_p(p_genotype), "\n")
cat("Stage x Genotype:", format_p(p_interaction), "\n")

# -----------------------

# -----------------------

cat("\nDecomposed model: distance_to_food ~ stage_label * Mito * Nuclear\n")
lm_dist_decomp <- lm(distance_to_food ~ stage_label * Mito * Nuclear, data = dist_df)
lm_dist_decomp_anova <- anova(lm_dist_decomp)

print(lm_dist_decomp_anova)
print(summary(lm_dist_decomp))

cat("\nP-values for reporting (NEW decomposed Mito x Nuclear model):\n")
for (term_name in rownames(lm_dist_decomp_anova)) {
  if (term_name != "Residuals") {
    cat(term_name, ":", format_p(lm_dist_decomp_anova[term_name, "Pr(>F)"]), "\n")
  }
}

# -----------------------

# -----------------------

cat("\nLevene's test - distance_to_food variance by developmental stage\n")
levene_stage <- car::leveneTest(distance_to_food ~ stage_label, data = dist_df, center = median)
levene_genotype <- car::leveneTest(distance_to_food ~ genotype_label, data = dist_df, center = median)
dist_df$stage_genotype <- interaction(dist_df$stage_label, dist_df$genotype_label, drop = TRUE)
levene_cells <- car::leveneTest(distance_to_food ~ stage_genotype, data = dist_df, center = median)
print(levene_stage)
print(levene_genotype)
print(levene_cells)
grDevices::pdf(file.path(stat_dir, 'linear_model_diagnostics.pdf'), width = 8, height = 8)
par(mfrow = c(2, 2))
plot(lm_dist, which = c(1, 2, 3, 5))
grDevices::dev.off()

# -----------------------
# B8) Build final clean Panel B figure
# -----------------------

stage_cols <- c("Larvae" = "#D55E00", "Adult females" = "#0072B2")

dist_df <- dist_df %>%
  mutate(
    genotype_label_clean = dplyr::recode(
      as.character(genotype_label),
      "(ore);OreR" = "(ore);\nOreR",
      "(simw501);OreR" = "(simw501);\nOreR",
      "(ore);Aut" = "(ore);\nAut",
      "(simw501);Aut" = "(simw501);\nAut"
    ),
    genotype_label_clean = factor(
      genotype_label_clean,
      levels = c("(ore);\nOreR", "(simw501);\nOreR", "(ore);\nAut", "(simw501);\nAut")
    )
  )

panel_B_final <- ggplot(
  dist_df,
  aes(x = genotype_label_clean, y = distance_to_food, fill = stage_label, color = stage_label)
) +
  geom_boxplot(
    position = position_dodge(width = 0.70), width = 0.52, alpha = 0.22,
    outlier.shape = NA, linewidth = 0.70
  ) +
  geom_point(
    aes(shape = stage_label),
    position = position_jitterdodge(jitter.width = 0.055, dodge.width = 0.70),
    size = 3.1, alpha = 0.96, stroke = 0.80
  ) +
  scale_fill_manual(values = stage_cols, name = "Generation 0 stage") +
  scale_color_manual(values = stage_cols, name = "Generation 0 stage") +
  scale_shape_manual(values = c("Larvae" = 17, "Adult females" = 16), name = "Generation 0 stage") +
  scale_x_discrete(
    labels = expression(
      atop("(ore);", "OreR"),
      atop((simw^501)*";", "OreR"),
      atop("(ore);", "Aut"),
      atop((simw^501)*";", "Aut")
    )
  ) +
  annotate(
    "text", x = 2.95, y = max(dist_df$distance_to_food, na.rm = TRUE) * 1.045,
    label = paste("Stage effect:", format_p(p_stage)), hjust = 0, vjust = 1, size = 4.0, color = "black"
  ) +
  labs(x = "Mitochondrial-nuclear genotype", y = "Bray-Curtis distance to\nmean food/source profile", title = NULL,
       caption = if (opts$demo) "Synthetic demonstration data" else NULL) +
  coord_cartesian(
    ylim = c(min(dist_df$distance_to_food, na.rm = TRUE) * 0.92,
             max(dist_df$distance_to_food, na.rm = TRUE) * 1.12),
    clip = "off"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold", size = 14, color = "black"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8), vjust = 0.45),
    axis.text = element_text(color = "black", size = 12),
    axis.text.x = element_text(face = "bold", angle = 0, hjust = 0.5, vjust = 0.5),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(color = "black", size = 11.5),
    legend.key.size = unit(0.58, "cm"),
    legend.box.margin = margin(t = -4, r = 0, b = 0, l = 0),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.28),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.70),
    axis.line = element_blank(),
    plot.margin = margin(t = 12, r = 18, b = 10, l = 18)
  )

panel_B_final

ggsave(file.path(fig_dir, "Figure1B_for_combined_Figure1.pdf"), plot = panel_B_final,
       width = 9.8, height = 4.2, units = "in", device = pdf_device, bg = "white", limitsize = FALSE)
ggsave(file.path(fig_dir, "Figure1B_for_combined_Figure1.png"), plot = panel_B_final,
       width = 9.8, height = 4.2, units = "in", dpi = 600, bg = "white", limitsize = FALSE)
ggsave(file.path(fig_dir, "Figure1B_for_combined_Figure1.jpeg"), plot = panel_B_final,
       width = 9.8, height = 4.2, units = "in", dpi = 600, bg = "white", device = "jpeg", limitsize = FALSE)

write.csv(dist_df, file.path(data_dir, "Figure1B_distance_to_food_data.csv"), row.names = FALSE)

# ============================================================
# PART C: SUPPLEMENTARY TABLE S1
# Food/source PERMANOVA + G0 distance-to-food linear models
# ============================================================

supp_s1_food_perm <- data.frame(
  Analysis = c("PERMANOVA", "PERMANOVA", "PERMANOVA"),
  Response = rep("Bray-Curtis dissimilarity among food/source samples", 3),
  Model = rep("food_bray ~ genotype", 3),
  Term = c("genotype", "Residual", "Total"),
  Df = food_perm$Df,
  SumOfSqs = food_perm$SumOfSqs,
  R2 = food_perm$R2,
  F = food_perm$F,
  P = food_perm$`Pr(>F)`
)

supp_s1_g0_lm <- data.frame(
  Analysis = "Linear model ANOVA",
  Response = "Bray-Curtis distance to mean food/source profile",
  Model = "distance_to_food ~ stage_label * genotype_label",
  Term = rownames(lm_dist_anova),
  Df = lm_dist_anova$Df,
  SumSq = lm_dist_anova$`Sum Sq`,
  MeanSq = lm_dist_anova$`Mean Sq`,
  F = lm_dist_anova$`F value`,
  P = lm_dist_anova$`Pr(>F)`
)

supp_s1_g0_lm$Term <- dplyr::recode(
  supp_s1_g0_lm$Term,
  "stage_label" = "developmental stage",
  "genotype_label" = "genotype",
  "stage_label:genotype_label" = "developmental stage x genotype",
  "Residuals" = "Residuals"
)

supp_s1_g0_lm_decomp <- data.frame(
  Analysis = "Linear model ANOVA (decomposed)",
  Response = "Bray-Curtis distance to mean food/source profile",
  Model = "distance_to_food ~ stage_label * Mito * Nuclear",
  Term = rownames(lm_dist_decomp_anova),
  Df = lm_dist_decomp_anova$Df,
  SumSq = lm_dist_decomp_anova$`Sum Sq`,
  MeanSq = lm_dist_decomp_anova$`Mean Sq`,
  F = lm_dist_decomp_anova$`F value`,
  P = lm_dist_decomp_anova$`Pr(>F)`
)

levene_stage_df <- as.data.frame(levene_stage)
levene_stage_df$Test <- "Levene's test: distance_to_food ~ stage_label"

levene_genotype_df <- as.data.frame(levene_genotype)
levene_genotype_df$Test <- "Levene's test: distance_to_food ~ genotype_label"

levene_cells_df <- as.data.frame(levene_cells)
levene_cells_df$Test <- "Levene's test: distance_to_food ~ stage x genotype groups"
supp_s1_levene <- bind_rows(
  levene_stage_df %>% rownames_to_column("Term"),
  levene_genotype_df %>% rownames_to_column("Term"),
  levene_cells_df %>% rownames_to_column("Term")
)

# -----------------------
# C3) Write Excel file (now with 4 sheets instead of 2)
# -----------------------

wb <- createWorkbook()

addWorksheet(wb, "Food_source_PERMANOVA")
writeData(wb, "Food_source_PERMANOVA", supp_s1_food_perm)

addWorksheet(wb, "G0_distance_to_food_LM")
writeData(wb, "G0_distance_to_food_LM", supp_s1_g0_lm)

addWorksheet(wb, "G0_LM_MitoNuclear_decomp")
writeData(wb, "G0_LM_MitoNuclear_decomp", supp_s1_g0_lm_decomp)

addWorksheet(wb, "Levene_dispersion_tests")
writeData(wb, "Levene_dispersion_tests", supp_s1_levene)

header_style <- createStyle(textDecoration = "bold", halign = "center")

for (sheet in c("Food_source_PERMANOVA", "G0_distance_to_food_LM",
                "G0_LM_MitoNuclear_decomp", "Levene_dispersion_tests")) {
  addStyle(wb, sheet, header_style, rows = 1, cols = 1:15, gridExpand = TRUE)
  setColWidths(wb, sheet, cols = 1:15, widths = "auto")
}

saveWorkbook(
  wb,
  file = file.path(stat_dir, "Supplementary_Table_S1_FoodSource_G0Distance_Statistics.xlsx"),
  overwrite = TRUE
)

write.csv(supp_s1_food_perm, file.path(stat_dir, "Supplementary_Table_S1A_Food_source_PERMANOVA.csv"), row.names = FALSE)
write.csv(supp_s1_g0_lm, file.path(stat_dir, "Supplementary_Table_S1B_G0_distance_to_food_LM.csv"), row.names = FALSE)
write.csv(supp_s1_g0_lm_decomp, file.path(stat_dir, "Supplementary_Table_S1C_G0_LM_MitoNuclear_decomp.csv"), row.names = FALSE)
write.csv(supp_s1_levene, file.path(stat_dir, "Supplementary_Table_S1D_Levene_dispersion_tests.csv"), row.names = FALSE)

# ============================================================
# PART D: SAVE FULL STATS OUTPUT TXT
# ============================================================

capture.output(
  {
    cat("============================================================\n")
    cat("FIGURE 1B + SUPPLEMENTARY FIGURE S1 + SUPPLEMENTARY TABLE S1\n")
    cat("Generation 0 food/source community and host divergence analysis\n")
    cat("PORTABLE EXAMPLE - inspect assumptions before interpreting results\n")
    cat("============================================================\n\n")
    
    cat("Samples included in Supplementary Fig. S1:\n")
    print(data.frame(sample_data(ps.food)) %>% count(genotype, life_stage, vial))
    
    cat("\nNumber of food/source samples:", nsamples(ps.food), "\n")
    cat("Number of taxa in ps.food:", ntaxa(ps.food), "\n")
    
    cat("\nTop genera in food/source samples:\n")
    print(top_genera)
    
    cat("\nSample counts after vial-level collapse:\n")
    print(table(meta_vial$sample_type))
    
    cat("\nSample counts by genotype and stage:\n")
    print(table(meta_vial$genotype_label, meta_vial$stage_label))
    
    cat("\nFood/source PERMANOVA: food community ~ genotype\n")
    print(food_perm)
    
    cat("\nFood/source beta-dispersion by genotype\n")
    cat(food_dispersion_status, "\n")
    if (!is.null(food_bd_perm)) print(food_bd_perm)
    
    cat("\n=== GENOTYPE MODEL ===\n")
    cat("Main model: distance_to_food ~ stage_label * genotype_label\n")
    print(lm_dist_anova)
    print(summary(lm_dist))
    
    cat("\nP-values for reporting (GENOTYPE genotype-level model):\n")
    cat("Stage:", format_p(p_stage), "\n")
    cat("Genotype:", format_p(p_genotype), "\n")
    cat("Stage x Genotype:", format_p(p_interaction), "\n")
    
    cat("\n=== DECOMPOSED MODEL ===\n")
    cat("Model: distance_to_food ~ stage_label * Mito * Nuclear\n")
    print(lm_dist_decomp_anova)
    print(summary(lm_dist_decomp))
    
    cat("\n=== LEVENE'S TESTS FOR DISPERSION ===\n")
    cat("Levene's test - distance_to_food variance by stage:\n")
    print(levene_stage)
    cat("\nLevene's test - distance_to_food variance by genotype:\n")
    print(levene_genotype)
    
    cat("\nStage contrasts; raw p.value and BH-adjusted p_BH:\n")
    print(emm_stage_by_genotype)
    cat("\nGenotype contrasts; raw p.value and BH-adjusted p_BH:\n")
    print(emm_genotype_by_stage)
    cat("\nLevene's test across all stage/genotype cells:\n")
    print(levene_cells)
    cat("\nExploratory model assumptions: independent host vials; ",
        "fixed shared mean food profile; Type I sequential ANOVA.\n")

  },
  file = file.path(stat_dir, "Figure1B_and_SuppFigS1_full_stats_output.txt")
)

# ============================================================
# DONE
# ============================================================

cat("\nDONE. Files saved inside folder:\n")
cat(main_outdir, "\n\n")
cat("Additional statistics files:\n")
cat(file.path(stat_dir, "Supplementary_Table_S1C_G0_LM_MitoNuclear_decomp.csv"), "\n")
cat(file.path(stat_dir, "Supplementary_Table_S1D_Levene_dispersion_tests.csv"), "\n")


# Preserve computational provenance without copying private inputs into the repository.
capture.output(sessionInfo(), file = file.path(stat_dir, 'sessionInfo.txt'))
if (!is.null(script_path)) file.copy(script_path, file.path(main_outdir, 'code'), overwrite = TRUE)
message(if (opts$demo) 'SYNTHETIC DEMO ONLY: ' else 'Research analysis: ',
        'outputs saved to ', normalizePath(main_outdir))
