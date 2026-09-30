# Run once: Rscript scripts/install_packages.R
if (getRversion() < '4.1.0') stop('R 4.1 or later is required.')
options(repos = c(CRAN = 'https://cloud.r-project.org'))
cran_packages <- c('dplyr', 'tidyr', 'tibble', 'vegan', 'ggplot2',
                   'stringr', 'scales', 'openxlsx', 'emmeans', 'car')
missing <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing)) install.packages(missing)
if (requireNamespace('ggplot2', quietly = TRUE) &&
    utils::packageVersion('ggplot2') < '3.4.0') {
  install.packages('ggplot2')
}
if (!requireNamespace('phyloseq', quietly = TRUE)) {
  if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')
  BiocManager::install('phyloseq', ask = FALSE, update = FALSE)
}
required <- c(cran_packages, 'phyloseq')
still_missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(still_missing)) {
  stop('Installation incomplete: ', paste(still_missing, collapse = ', '))
}
message('Required packages are installed.')
