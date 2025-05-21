#' KIRC expression and survival data
#'
#' A subset of The Cancer Genome Atlas Kidney Renal Clear Cell Carcinoma
#' (TCGA-KIRC) cohort, pre-processed for the DegreeAFT vignette.
#'
#' @format A list with three components:
#' \describe{
#'   \item{expr}{numeric matrix, genes (rows) × samples (cols)}
#'   \item{survival}{data frame with columns of sample, OS, OS.time}
#'   \item{marker}{selected markers through univariate log rank test}
#' }
#' @source The Cancer Genome Atlas
"kirc"
