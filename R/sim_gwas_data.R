#' Simulated GWAS summary statistics
#'
#' A dataset containing simulated genome-wide association summary statistics for
#'  use in the rolypoly vignette.
#'
#' @format A data frame with 14934 rows and 6 variables:
#' \describe{
#'   \item{chrom}{chromosome, we only use autosomes}
#'   \item{pos}{base pair position of variant}
#'   \item{rsid}{rsid identifier of variant}
#'   \item{beta}{effect size, univariate regression coefficient}
#'   \item{se}{standard error of effect size}
#'   \item{maf}{minor allele frequency}
#' }
#' @source rsids were from 1000g and I generated the other fields
"sim_gwas_data"
