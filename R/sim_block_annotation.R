#' Simulated block data annotation.
#'
#' A dataset containing simulated block data annotation for
#'  use in the rolypoly vignette.
#'
#' @format A data frame with 1000 rows and 4 variables:
#' \describe{
#'   \item{chrom}{chromosome, we only use autosomes}
#'   \item{start}{base pair position of variant}
#'   \item{end}{rsid identifier of variant}
#'   \item{label}{effect size, univariate regression coefficient}
#' }
#' @source I generated these fields to link with SNP positions
"sim_block_annotation"
