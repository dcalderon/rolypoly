#' Inference functions.
#'
#' Functions for inferring relevant annotations using the polyTest model.
#' @name inference
#' @importFrom foreach foreach %dopar% %do%
NULL


#' Run inference.
#'
#' Coordinates running inference.
#'
#' @param rolypoly rolypoly object
#' @param bootstrap_iters number of bootstrap iterations to perform
#' @param outlier_threshold threshold for performing robust regression, still experimental.
#' @param run_light if we throw out bootstrap data
#' @param run_parallel if we collect bootstraps in parallel
#' @return rolypoly object with inference information attached
#' @examples
#' \dontrun{run_inference(rolypoly)}
rolypoly_perform_inference <- function(rolypoly, bootstrap_iters = 50,
                                       outlier_threshold = -1, run_light = F,
                                       run_parallel = F) {
  if (is.null(rolypoly$data)) {
    warning('data has not been precomputed, returning without results')
    return(rolypoly)
  }

  message('starting inference procedure')
  # first identify outliers, might not even be needed
  if (outlier_threshold > 0) {
    stop('not yet implemented')
    rolypoly$robust_results <- robust_parameter_estimator(
      vectorize_rolypoly(rolypoly$data), outlier_threshold = outlier_threshold
    )
  }

  # then fit model
  rolypoly$full_results <- parameter_estimator(vectorize_rolypoly(rolypoly$data))

  # add on block values
  rolypoly$block_values <- calculate_block_values(
    rolypoly$raw_block_data, rolypoly$full_results$parameters)

  # add on heritability values
  rolypoly$block_heritability_contribution <-
    calculate_annotation_block_heritability(rolypoly$raw_block_data,
                                            rolypoly$full_results$parameters)

  # calculate expected block values (accounts for ld and snp error)
  rolypoly <- calculate_expected_block_values_given_ld(rolypoly, rolypoly$block_values)

  # Bootstrap error and 95% confidence interval estimates
  if (bootstrap_iters > 0) {
    rolypoly <- bootstrap_estimator(rolypoly, bootstrap_iters = bootstrap_iters,
                                    run_light = run_light, run_parallel = run_parallel)
  }
  return(rolypoly)
}

#' Find parameter estimates for the data.
#'
#' @param vectorized_rolypoly_data rolypoly data that has been vectorized
#' @return results of inference
#' @examples
#' \dontrun{parameter_estimator(rolypoly)}
parameter_estimator <- function(vectorized_rolypoly_data) {
  full_results <- list()
  m <- stats::lm(vectorized_rolypoly_data$y ~
            offset(vectorized_rolypoly_data$noise_per_snp) +
            vectorized_rolypoly_data$x);
  full_results$parameters <- stats::coef(m)
  annotation_names <- c('intercept', colnames(vectorized_rolypoly_data$x))
  names(full_results$parameters) <- annotation_names
  full_results$model <- m
  return(full_results)
}

#' Bootstrap parameter estimates for confidence intervals.
#'
#' @param rolypoly rolypoly object
#' @param bootstrap_iters number of bootstrap iterations to run
#' @param run_light if we throw away bootstrap data
#' @param run_parallel if we want to collect bootstrap data in parallel
#' @examples
#' \dontrun{bootstrap_estimator(rolypoly)}
bootstrap_estimator <- function(rolypoly, bootstrap_iters, run_light, run_parallel) {
  partitions_present <- unique(sapply(rolypoly$data, function(block) {block$partition}))

  # run things in parallel if user specified
  "%loop_function%" <- `%do%`
  if (run_parallel) { "%loop_function%" <- `%dopar%` }

  rolypoly$boot_data <-
    foreach(i = 1:bootstrap_iters, .combine = 'c', .inorder = F,
            .verbose = F, .packages = c('rolypoly')
    ) %loop_function% {
      message(paste('starting bootstrap iteration',  i, sep = ': '))
      boot_partitions <- sample(partitions_present, length(partitions_present), replace = T)
      boot_results <- parameter_estimator(
        vectorize_rolypoly(rolypoly$data[
          sapply(rolypoly$data, function(block) { block$partition %in% boot_partitions})
          ]))
      # return the important bits
      return(list(list(
        gammas = boot_results$parameters,
        block_values = calculate_block_values(rolypoly$raw_block_data,
                                              boot_results$parameters),
        block_heritability = calculate_annotation_block_heritability(
          rolypoly$raw_block_data, boot_results$parameters
        )
      )))
    }
  rolypoly$bootstrap_results <- make_results_df(
    sapply(rolypoly$boot_data, function(boot) {boot$gammas}),
    names(rolypoly$full_results$parameters),
    rolypoly$full_results$parameters
  )
  # # block value results
  rolypoly$bootstrap_block_values <- make_results_df(
    sapply(rolypoly$boot_data, function(boot) {boot$block_values}),
    rownames(rolypoly$raw_block_data),
    rolypoly$block_values
  )
  # # heritability contribution values
  rolypoly$bootstrap_block_heritability_contribution <- make_results_df(
    sapply(rolypoly$boot_data, function(boot) {boot$block_heritability}),
    colnames(rolypoly$raw_block_data),
    rolypoly$block_heritability_contribution
  )
  if (run_light) {rolypoly$boot_data <- NULL}
  return(rolypoly)
}

#' Find robust parameter estimates for the data.
#'
#' @param vectorized_rolypoly_data vectorized rolypoly data
#' @param outlier_threshold outlier threshold for robust inference
#' @examples
#' \dontrun{parameter_estimator(rolypoly)}
robust_parameter_estimator <- function(vectorized_rolypoly_data, outlier_threshold = 10) {
  full_results <- list()
  m <- MASS::rlm(vectorized_rolypoly_data$y ~ offset(vectorized_rolypoly_data$noise_per_snp) +
                   vectorized_rolypoly_data$x, psi = MASS::psi.huber, k = outlier_threshold)
  full_results$parameters <- stats::coef(m)
  annotation_names <- c('intercept', colnames(vectorized_rolypoly_data$x))
  names(full_results$parameters) <- annotation_names
  full_results$model <- m
  return(full_results)
}

#' Take a list of rolypoly data and vectorize it for inference.
#'
#' @param data the list of block information from rolypoly object
#' @return list of necessary information for inference
#' @examples
#' \dontrun{vectorize_rolypoly(data)}
vectorize_rolypoly <- function(data) {
  # use only blocks flagged for inference inclusion
  data <- data[sapply(data, function(block) {block$include_in_inference})]

  # unpack data
  y <- do.call('c', lapply(data, function(block) { block$y }))
  x <- do.call('rbind', lapply(data, function(block) { block$x }))
  rownames(x) <- do.call('c', lapply(data, function(block) { block$snps$rsid }))
  noise_per_snp <- do.call('c', lapply(data, function(block) { block$snps$se**2 }))

  # exclude na elements
  na_elements <- is.na(y) | apply(x, 1, function(x) { any(is.na(x))}) | is.na(noise_per_snp)
  return(list(
    y = y[!na_elements], x = x[!na_elements,],
    noise_per_snp = noise_per_snp[!na_elements])
  )
}

#' Caclulate predicted block values based on block information and model fit.
#'
#' From a model fit we can predict expected variance of a block based on information
#' we have about the block. In the example of gene expression this would equate
#' to predicting the importance of a gene based on its signature of expression.
#'
#' @param block_data block_data
#' @param params parameter fit
#' @return returns block values
#' @examples
#' \dontrun{calculate_gene_values(block_data, params)}
calculate_block_values <- function(block_data, params) {
  if (any(is.na(params))) {
    warning('NA pameters found, take results with grain of salt')
    params[is.na(params)] <- 0
  }
  # only include params for which we have block data
  block_values <- as.numeric(as.matrix(block_data) %*% params[colnames(block_data)])
  names(block_values) <- rownames(block_data)
  return(block_values)
}

#' Caclulate predicted block values based on block information and model fit.
#'
#' From a model fit we can predict expected variance of a block based on information
#' we have about the block. In the example of gene expression this would equate
#' to predicting the importance of a gene based on its signature of expression.
#'
#' @param rolypoly rolypoly object
#' @param block_values estimated block values.
#' @examples
#' \dontrun{calculate_expected_block_values_given_ld(rolypoly, block_values)}
calculate_expected_block_values_given_ld <- function(rolypoly, block_values) {
  expected_block_values <- lapply(rolypoly$data, function(block) {
    if (!block$include_in_inference) { return(NULL) }
    if (!(block$block_info$label %in% names(block_values))) { return(NULL)}
    return(list(
      # the expected value is equal to
      # the sum of the snp variance + sum of ld_squared * block value
      # The equivalent of the trace of ld %*% ld is sum(ld * ld) since ld symmetric
      g_hat = sum(block$snps$se**2) + block_values[block$block_info$label] *
        sum(block$ld_matrix_squared),
      g = block$gwas_block_score,
      label = block$block_info$label
    ))
  })
  rolypoly$expected_block_values <- data.table::rbindlist(
    expected_block_values[!sapply(expected_block_values, is.null)]
  )
  return(rolypoly)
}

#' Caclulate the contribution of block annotations to the heritability of a trait.
#'
#' A vector of indenpendent heritability contributions of block annotations is returned. Sum
#' the vector to get total explained heritability and divide by sum to get proportion.
#'
#' @param block_data functional information of blocks
#' @param params parameter fit
#' @examples
#' \dontrun{calculate_annotation_block_heritability(block_data, params)}
calculate_annotation_block_heritability <- function(block_data, params) {
  if (any(is.na(params))) {
    warning('NA pameters found, take results with grain of salt')
    params[is.na(params)] <- 0
  }
  contributions <- rowSums(params[colnames(block_data)] * t(block_data))
  names(contributions) <- colnames(block_data)
  return(contributions)
}

#' Helper function to make a summary table of results from bootstrap data.
#'
#' @param value_collection collection of bootstrapped value estimates
#' @param annotations vector of annotation names
#' @param model_estimates estimates for bias parameter estimates
#' @return data frame with results summary
#' @examples
#' \dontrun{make_results_df(value_collection)}
make_results_df <- function(value_collection, annotations, model_estimates) {

  # in the case we're calculating single parameter estimates
  if (is.null(dim(value_collection))) {value_collection <- matrix(value_collection, nrow = 1)}

  parameter_estimates <- data.frame(
    annotation = annotations,
    bootstrap_estimate = rowMeans(value_collection),
    bootstrap_error = apply(value_collection, 1, stats::sd)
  ) %>%
    dplyr::mutate(
      bt_value = bootstrap_estimate / bootstrap_error,
      bp_value = stats::pnorm(bt_value, lower.tail = F),
      bias_corrected_estimate = 2 * model_estimates - bootstrap_estimate
    )

  parameter_estimates$CI_lo <-
    apply(value_collection, 1, function(x) { stats::quantile(x, 0.025, na.rm = T) })
  parameter_estimates$CI_hi <-
    apply(value_collection, 1, function(x) { stats::quantile(x, 0.975, na.rm = T) })

  return(parameter_estimates)
}
