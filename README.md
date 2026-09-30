# Microbiome analysis code — Nitin Bansal

Four R scripts from doctoral research on how host mitochondrial and nuclear
genotypes and developmental stage shape the Drosophila microbiome across generations.
The scripts show data subsetting, biological replicate aggregation, community
comparisons, statistical models, and figure and table generation.

**Data status:** The processed research dataset is not included at present.
These files are shared as code examples for review. They require the original
phyloseq input to run; there is no demo dataset or demo command in this package.

## Scripts

| File | Question and analysis | Main outputs |
|---|---|---|
| `02_generation_divergence.R` | How does oreore community composition change across generations? Vial means, Bray–Curtis distance to stage-specific Gen0/Gen1 mean profiles, PERMANOVA, dispersion checks, and linear models. | Figure 2, supplementary Table S2, source CSVs |
| `03_genus_trajectories.R` | Which bacterial genera change across generations and replicate vials? Dominant genera and an abundance-filtered genus view. | Figure 3, supplementary Figure S3, source CSVs |
| `04_mitonuclear_divergence.R` | Does Gen5 divergence from Gen0 differ by mitochondrial genotype, nuclear background, and stage? PCoA, PERMANOVA, vial-intercept mixed models, contrasts, and a vial-mean sensitivity analysis. | Figure 4, supplementary Figures S2/S5, Table S3 |
| `05_genus_abundance_changes.R` | Which genera contribute to Gen0-to-Gen5 changes? Vial-level abundance differences and per-genus interaction models with BH correction. | Figure 5, supplementary Figure S4, Table S4 |

Each script is independent and reads the same input. Numbers correspond to the
research figures, not a required execution order. These scripts start from a
processed phyloseq object; sequencing preprocessing and ASV inference are upstream.

## Required input

The default path is `data/06_phyloseq_clean_CHAP1_NOHOST.rds`. Add that file when
available, or edit `INPUT_RDS` at the top of each script. The RDS must contain a
phyloseq object with an ASV abundance table and sample metadata. The genus scripts
also require `Genus` and `Family` taxonomy ranks.

| Metadata field | Values used by these scripts |
|---|---|
| `genotype` | `oreore`, `simore`, `oreaut`, `simaut` |
| `generation` | `Generation_0` through `Generation_5` |
| `life_stage` | `Larvae`, `Female`; `food` is additionally used in the trajectory script |
| `vial` | Replicate lineage/vial identifier; trajectory panels expect `Vial_1`, `Vial_2`, `Vial_3` |
| `replicate` | Pool identifier, required by `02_generation_divergence.R` |

The study must include the relevant baseline and follow-up groups. Missing metadata,
unbalanced cells, and repeated lineage identifiers need scientific review before reuse.

## Packages

Use R with `phyloseq`, `dplyr`, `tidyr`, `tibble`, `vegan`, `ggplot2`, `readr`,
`patchwork`, `openxlsx`, `stringr`, `forcats`, `scales`, `emmeans`, `nlme`, and `car`.
The plots use the `linewidth` argument supported by ggplot2 3.4.0 and later.
Each script checks its own dependencies. `car` is used through `car::` to avoid
masking `dplyr::recode()`.

For a fresh R environment:

```r
install.packages(c("dplyr", "tidyr", "tibble", "vegan", "ggplot2", "readr",
                   "patchwork", "openxlsx", "stringr", "forcats", "scales",
                   "emmeans", "nlme", "car"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("phyloseq")
```

## Running after the input is available

Set the working directory to this folder, which contains `README.md` and the four
R scripts. Then run a script from RStudio with Source, or from a terminal:

```sh
Rscript 02_generation_divergence.R
Rscript 03_genus_trajectories.R
Rscript 04_mitonuclear_divergence.R
Rscript 05_genus_abundance_changes.R
```

Outputs are written under `outputs/` into separate figure-specific folders. They
include figures, source tables, statistics, a copy of the processed input, and
`sessionInfo.txt`. Existing outputs with the same names can be overwritten. The
PDF device uses Cairo when available and the standard PDF device otherwise.

## Interpretation and validation status

These are research scripts for a specific study design. The models and filtering
choices are retained from the original analyses, with the corrections listed in
`CHANGES.md`. Bray–Curtis references are arithmetic mean relative-abundance
profiles; the original figure labels use the word “centroid.” They are not PCoA
centroids. Genus abundance changes describe relative proportions, not absolute
bacterial quantities or a dedicated compositional differential-abundance test.

The source was reviewed and checked for balanced delimiters and portable paths.
The revised scripts have **not been executed in R** during this cleanup because
R and the processed research input were unavailable. Figures and statistics must
be regenerated and compared before the revised outputs are used in a manuscript.
The original free-permutation design and model assumptions also need review, as
explained in `CHANGES.md`.

Author: Nitin Bansal. Research code edited for sharing with assistance; substantive
corrections and remaining limitations are documented in `CHANGES.md`.
