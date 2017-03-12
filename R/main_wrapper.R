#' Main wrapper functions.
#'
#' The main endpoint user functions.
#' @name main_wrapper
NULL


#' Main rolypoly wrapper function.
#'
#' The entry point for rolypoly analysis. If no expression data, we assume
#' that we are running just the vegas score processing.
#'
#' @param rolypoly Previous rolypoly run to parts of pipeline.
#' @param gwas_data Gwas data for a trait, including snp annotations.
#' @param block_annotation Start and end points for blocks
#' @param block_data Information about blocks.
#' @param ld_folder Folder with LD information.
#' @param bootstrap_iters Number bootstrap iterations to perform for inference.
#' @param outlier_threshold Set to positive if we want to run robusted regression.
#' @param perform_cv If we want to interpret annotation effects do not set this to
#'  T. However, if our goal is prediction accuracy then set this to T.
#' @param gwas_z_filter Z-score filter for SNPs, helps prevent large effects
#'  biasing inference.
#' @param add_spline If we want to fit a spline to maf.
#' @param n_knots number of knots to add to the spline.
#' @param n_folds number of folds for cross validation
#' @param add_poly If we want to fit a polynomial to maf.
#' @param n_degree the degree of the polynomial.
#' @param run_light if we want to throw away bootstrap data, and save memory
#' @param gwas_link_parallel if user wants to run in gwas linking in parallel,
#'  registerDoParallel must have been run in advance.
#' @param bootstrap_parallel if user wants to run in bootstraps in parallel,
#'  registerDoParallel must have been run in advance.
#' @param keep_model if we should keep the regression model, can be large.
#' @param keep_gwas set to T if we want to include gwas in returned rolypoly object.
#' @param ... other arguments to pass to cv.glmnet
#' @return rolypoly object
#' @export
#' @examples
#' \dontrun{rolypoly_roll(rolypoly)}
rolypoly_roll <- function(rolypoly = NULL, gwas_data = NULL, block_annotation = NULL,
                          block_data = NULL, ld_folder = NULL, bootstrap_iters = 50,
                          outlier_threshold = -1, perform_cv = F, n_folds = 10,
                          gwas_z_filter = -1, add_spline = F,  n_knots = 1,
                          add_poly = F, n_degree = 2, run_light = T,
                          gwas_link_parallel = F, bootstrap_parallel = F, keep_model = F,
                          keep_gwas = F, ...) {
  if (is.null(rolypoly)) { rolypoly <- list(); class(rolypoly) <- 'rolypoly' }
  if (!is.null(gwas_data))
    rolypoly <- rolypoly_load_gwas(rolypoly, gwas_data, gwas_z_filter = gwas_z_filter,
                                   add_spline = add_spline, n_knots = n_knots,
                                   add_poly = add_poly, n_degree = n_degree)
  if (!is.null(block_annotation))
    rolypoly <- rolypoly_load_block_annotation(rolypoly, block_annotation)
  if (!is.null(ld_folder))
    rolypoly <- rolypoly_link_blocks_and_gwas(rolypoly, ld_folder, run_parallel = gwas_link_parallel)
  if (!keep_gwas) {rolypoly$gwas <- NULL}
  if (!is.null(block_data))
    rolypoly <- rolypoly_load_block_data(rolypoly, block_data)
  if (is.null(rolypoly$data)) {return(rolypoly)} # just building object

  if (is.null(rolypoly$raw_block_data)) {return(rolypoly)}

  rolypoly <- rolypoly_perform_inference(rolypoly, bootstrap_iters = bootstrap_iters,
                                         outlier_threshold = outlier_threshold,
                                         run_light = run_light, run_parallel = bootstrap_parallel)
  if (perform_cv) {
    rolypoly <- rolypoly_perform_regularized_inference(rolypoly, n_folds = n_folds, ...)
  }
  if (!keep_model) {
    rolypoly$full_results$model <- NULL
  }
  return(rolypoly)
}
