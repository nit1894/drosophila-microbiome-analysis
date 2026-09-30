# Microbiome food-reference analysis in R

An analysis example adapted from my doctoral work on host-associated microbial communities. The workflow compares Generation 0 larval and adult-female communities with a shared food/source reference, accounting for the experimental sampling structure.

**Author:** Nitin Bansal, PhD  
**Main tools:** phyloseq, vegan, ggplot2, dplyr, car, emmeans, and openxlsx.

## What the code demonstrates

- Validating sequencing counts and sample metadata before analysis.
- Averaging relative-abundance profiles from multiple host pools within each experimental vial.
- Calculating Bray-Curtis dissimilarity to a shared mean food/source profile.
- Comparing developmental stage and mitochondrial-nuclear genotype with linear models.
- Inspecting model diagnostics, testing variance differences, and documenting limitations.
- Exporting figures, source data, model tables, and reproducibility information.

This example starts with **processed ASV counts and taxonomy**. It does not perform raw-read processing, DADA2 inference, shotgun metagenomics, or machine-learning training.

## Quick start

Install R 4.1 or later. Download this repository and open its folder as the working directory. In a terminal, run:

```bash
Rscript scripts/install_packages.R
Rscript scripts/analyze_food_reference.R --demo
Rscript scripts/check_demo_outputs.R
```

The first command installs the required R packages, including phyloseq through Bioconductor. Package installation requires internet access and, on some systems, development libraries for package compilation.

Alternatively, in RStudio with the repository folder as the working directory:

```r
source("scripts/install_packages.R")
system2(file.path(R.home("bin"), "Rscript"),
        c("scripts/analyze_food_reference.R", "--demo"))
source("scripts/check_demo_outputs.R")
```

The bundled demonstration has **80 simulated samples and 12 illustrative taxa**: eight food samples and 72 host pools that collapse to 24 host vials. Counts and sample IDs are synthetic; the genus and genotype labels illustrate the categories used by the analysis. Demonstration statistics are not research findings.

**Validation status:** the source and synthetic input structure have been checked statically, and the reference-distance calculations have been independently checked in Python. R and the original research RDS were unavailable in the preparation environment, so the R workflow has not yet been executed there. Run the demonstration and output checks in R before relying on this example or publishing it as a tested workflow.

## Analyze a research phyloseq object

```bash
Rscript scripts/analyze_food_reference.R \
  --input data/private/your_processed_phyloseq.rds \
  --output outputs/research
```

Use a fresh output directory. Reusing a nonempty directory requires the explicit `--overwrite` option. Existing output files with matching names will then be replaced; unrelated files are retained.

The phyloseq object must contain:

| Component | Required content |
|---|---|
| OTU table | Finite, nonnegative, untransformed integer ASV counts; either table orientation is supported |
| Taxonomy | `Genus` and `Family` columns, with taxon IDs matching the OTU table |
| Sample metadata | `genotype`, `generation`, `life_stage`, `vial`, and `replicate` |

Metadata values are case-sensitive:

- `genotype`: `oreore`, `simore`, `oreaut`, or `simaut`.
- `generation`: `Generation_0` for the host samples included here.
- `life_stage`: `food`, `Larvae`, or `Female`.
- `vial`: a biological-vial identifier within genotype, generation, and stage.
- `replicate`: a pool/sample identifier. Multiple host pools in a vial are averaged before analysis.

Food samples are selected by life stage and genotype, regardless of their generation label. Confirm that they are the intended source/reference samples. Each food sample is treated as an independent source observation in PERMANOVA; shared genotype/vial labels generate a warning and require review of the sampling design.

## Analysis choices and limitations

### Reference and experimental unit

Each retained sample is converted to relative abundance after zero-count samples are removed. Host-pool profiles are averaged equally within each vial. The food reference is the arithmetic mean of the retained food-sample relative-abundance profiles, with equal weight per source sample. It is a **mean composition in ASV space**, rather than a centroid from a PCoA embedding.

Bray-Curtis distances to that reference describe compositional divergence. They do not establish transmission from food, absolute bacterial growth, or a causal mechanism. Models condition on the estimated reference; uncertainty in the food reference is not propagated.

The linear models assume independent host vials. If larval and adult observations share the same biological vial, lineage, or another repeated-measure unit, that dependence must be represented in an appropriate model. The current metadata labels alone do not establish independence.

### Statistical models

The genotype model is `distance_to_food ~ stage_label * genotype_label`. A second parameterization is `distance_to_food ~ stage_label * Mito * Nuclear`. With a complete two-by-two genetic design and all interactions, these describe the same stage-by-genotype cell means; they are not independent confirmations of a result.

`anova(lm(...))` supplies sequential Type I sums of squares. Term order matters for an unbalanced design. The exported omnibus model p-values are unadjusted and the analyses are exploratory. Inspect the interactions and the cell sizes before interpreting main effects.

Two post-hoc contrast families are exported with raw and Benjamini-Hochberg-adjusted p-values:

1. Stage comparisons across all four genotypes: four contrasts.
2. Genotype comparisons across both stages: twelve contrasts.

Median-centered Levene tests are calculated for stage, genotype, and the eight stage-by-genotype cells. These assess variance in the scalar reference-distance response, whereas PERMDISP assesses multivariate dispersion. The source PERMDISP test is skipped if any genotype has fewer than three source observations. Source PERMANOVA is retained, but a nonsignificant result with small groups does not demonstrate equivalent source communities.

The response is bounded between zero and one. Linear-model residual, Q-Q, variance, and influence plots are exported for inspection. Variance-test results do not automatically select or validate a model; interpretation requires reviewing diagnostics and the sampling design.

### Composition plot

Taxa lacking a genus annotation are retained and labeled `Unclassified: FAMILY` or `Unclassified`. The ten most abundant displayed genera are shown separately; remaining labels are combined as `Other`.

No genus is removed by default. If a justified display exclusion is needed:

```bash
Rscript scripts/analyze_food_reference.R --demo \
  --exclude-genus Gardnerella --output outputs/demo_filtered
```

This option affects the food-composition plot and its source data, which are renormalized after exclusion. It does not remove ASVs from the reference-distance analysis. An exclusion needs a documented scientific reason; a taxon name alone does not establish contamination.

## Outputs

Outputs are saved under `outputs/demo/` or the requested directory:

- `figure_files/`: composition and reference-distance plots as PDF, PNG, and JPEG.
- `figure_source_data/`: plotted compositions and vial-level host distances.
- `statistics/`: ANOVA/variance-test tables, post-hoc contrasts, source PERMANOVA, model-diagnostic PDF, QC counts, full text results, and `sessionInfo.txt`.
- `input_processed/`: a processed input snapshot, selected metadata, and command options.
- `code/`: the analysis-script snapshot used for the run.

The distance-plot p-value is generated from the fitted model. A fixed random seed makes permutation tests and plotted jitter reproducible within a compatible package environment; installed versions are recorded in `sessionInfo.txt`. This repository does not yet include a tested dependency lockfile.

Research `.rds` files and generated outputs are excluded by `.gitignore`. The synthetic CSV inputs are included so the public example does not require unpublished research data. Review the actual files selected for upload when using GitHub's browser uploader.

## Changes from the research script

- Replaced the lab-specific `setwd()` and fixed input filename with command-line input/output options.
- Added a synthetic CSV demonstration and explicit input checks.
- Removed zero-count samples before normalization.
- Calculated the figure annotation from the fitted model.
- Required statistical packages instead of substituting Bartlett's test for Levene's test or omitting contrasts.
- Added explicit BH correction across each stated contrast family.
- Added residual diagnostics and a variance check across all stage-by-genotype cells.
- Preserved unclassified taxa in the composition plot and made genus exclusion optional.
- Skipped poorly supported source-dispersion tests and recorded the reason.

These changes improve portability and make analysis decisions explicit. They do not verify the original empirical results; run and review the adapted workflow on the research input separately before using it for a manuscript.

## Method documentation

- [phyloseq data import](https://joey711.github.io/phyloseq/import-data.html)
- [vegan dispersion analysis](https://vegandevs.github.io/vegan/reference/betadisper.html)
- [R linear-model ANOVA](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/anova.lm.html)
- [car Levene test](https://rdrr.io/cran/car/man/leveneTest.html)
- [emmeans multiplicity and grouped comparisons](https://rvlenth.github.io/emmeans/articles/confidence-intervals.html)
