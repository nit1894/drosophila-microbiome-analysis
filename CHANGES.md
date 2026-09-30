# Changes and points to review

Original uploaded text files were left unchanged. These R copies retain the study
subsets, abundance thresholds, model formulas, and output filenames. Changes that
can affect numerical results are identified below.

## Sharing and portability

- Replaced the laboratory-specific `setwd()` path with editable `INPUT_RDS` and
  `OUTPUT_ROOT` settings. Run from the folder containing the README.
- Put generated files under `outputs/`, added missing-input and dependency messages,
  checked required metadata and taxonomy ranks, and recorded `sessionInfo()`.
- Added a standard PDF-device fallback and removed empty `code/` folders.
- Replaced editing-history comments with analysis descriptions. A combined
  generation-by-stage dispersion comparison is described as an omnibus group test,
  rather than a factorial interaction test. Its computation is unchanged.

## Corrections that require a rerun

1. **Genus trajectories — Figure 3:** The original `Other` category averaged
   abundance over constituent genera as well as pools. It now sums the genera in
   each pool, then averages pools within a vial. Individually plotted genus means
   retain their intended definition; `Other` values will change.
2. **Expanded genus screen — supplementary Figure S4:** The original delta table
   set factor levels to the top ten plotted genera before testing all retained
   genera. Non-top-ten names became `NA`, so the statistical screen could lose or
   combine them. The revised table retains all genus identities; only the plotted
   table is restricted to the top ten. Screen P values and BH corrections may change.
3. **Mito-nuclear contrasts — Figure 4/S5:** `emmeans` adjusts within each `by`
   group by default. With one mitochondrial contrast per nuclear-by-stage group,
   the original `adjust="BH"` did not correct across all four contrasts. The copies
   export `p.raw` and apply BH explicitly across each set of four: full mixed-model
   contrasts, separately fitted focused contrasts, and the vial-mean sensitivity
   contrasts. Those are three separate analysis families. Confidence intervals
   remain pointwise and unadjusted. The model fits and contrast estimates are unchanged.
4. **Variance check — Figure 4:** `car` is now required, and median-centered
   `car::leveneTest()` is used explicitly. The original Bartlett fallback both changed
   the test and could fail during tabular export. With `car` installed this retains
   the intended variance test; without it the script now stops before analysis.

## Scientific choices retained for review

- **Independence and permutations:** Figure 4's Gen0 community tests use pooled
  samples with free permutations even when several pools come from one vial. The
  mixed-model random intercept does not repair the separate PERMANOVA design.
  Determine the independent experimental unit and appropriate whole-vial/lineage
  permutation scheme from the study design before interpreting those tests.
- **Longitudinal dependence:** Figure 2 averages pools within a vial, but uses
  ordinary linear models and free permutations across generations. If vial labels
  track persistent lineages, the longitudinal dependence still needs to be handled.
- **PERMANOVA term hierarchy:** `by="margin"` is preserved. Inspect the returned
  terms: a term involved in a retained higher-order interaction is not necessarily
  tested as a separate main effect. Do not report unreturned tests as if computed.
- **Shared baseline:** Distances use an estimated group mean baseline profile.
  The existing models treat that reference as fixed and do not propagate its
  uncertainty. In Figure 4 the pool-level baseline and vial-mean sensitivity baseline
  may weight vials differently if pool counts are unequal.
- **Scalar variance check:** Levene's test on Figure 4 pool distances is exploratory:
  pools share vials. It does not replace mixed-model residual diagnostics or evaluation
  of group-specific residual variances.
- **Genus annotations and filtering:** The original `tax_glom(..., NArm=TRUE)` is
  retained, so unassigned genus annotations are omitted before the family-label
  fallback. The trajectory script also excludes `Gardnerella`; its scientific
  rationale is not established by this code review. Supplementary Figure S3 then
  renormalizes over retained genera. Its proportions differ in denominator from
  the main trajectory plot. Review and explain these choices.
- **Across-generation matching:** Genus deltas pair Gen0 and Gen5 by vial label.
  Confirm that labels identify matched replicate lineages; a shared label alone
  does not establish a paired design.
- **Multiple analyses:** The two-focal-genus screen and expanded screen retain
  their separate BH families. A focal genus can have different adjusted P values
  in those analyses. Focused, full-model, and sensitivity contrasts are alternative
  analyses of the same observations, not independent confirmation.
- **Bounded/compositional outcomes:** Review model residuals and assumptions for
  Bray–Curtis distances and changes in relative abundance. Relative-abundance
  results alone cannot establish absolute abundance changes.

## Validation performed

Static delimiter/string checks, input/output-path inspection, dependency inventory,
source-to-copy diffs, and ZIP integrity checks. No R execution or research-data
validation was possible in this environment. This is a reviewed source package,
not an assertion that regenerated scientific outputs have been validated.

## Reference documentation

- [emmeans: grouping and multiplicity adjustments](https://rvlenth.github.io/emmeans/articles/confidence-intervals.html)
- [vegan: PERMANOVA and permutation controls](https://vegandevs.github.io/vegan/reference/adonis.html)
- [vegan: marginal tests and model hierarchy](https://vegandevs.github.io/vegan/reference/anova.cca.html)
