#' Inference functions that include regularization
#'
#' Functions for inferring relevant annotations using the polyTest model.
#' @name regularized_inference
NULL


#' Perform regularization inference.
#'
#' Use CV to find appropriate values of lambda for either feature selection
#' or for prediction.
#'
#' @param vectorized_rolypoly_data rolypoly data used for inference
#' @param n_folds number of folds for cross validation
#' @param ... other arguments to pass to cv.glmnet
#' @return results from cross validation
#' @examples
#' \dontrun{cv_regularized_parameter_estimator(vectorized_rolypoly_data)}
cv_regularized_parameter_estimator <- function(vectorized_rolypoly_data, n_folds = 10, ...) {
  full_results <- list()
  m <- glmnet::cv.glmnet(
    x = vectorized_rolypoly_data$x,
    y = vectorized_rolypoly_data$y,
    offset = vectorized_rolypoly_data$noise_per_snp,
    foldid = cut(1:length(vectorized_rolypoly_data$y), breaks = n_folds, labels = F),
    family = 'gaussian', ... = ...
  )

  # can choose coefficients with either: lambda.min, or lambda.1se
  full_results$parameters <- stats::coef(m, s = 'lambda.min') %>% as.numeric()
  annotation_names <- c('intercept', colnames(vectorized_rolypoly_data$x))
  names(full_results$parameters) <- annotation_names
  full_results$model <- m
  return(full_results)
}

#' Run inference with added regularization.
#'
#' If p-values are desired use the other inference function. This for
#' prediction purposes.
#'
#' @param rolypoly a rolypoly object
#' @param ... other arguments to pass to cv.glmnet
#' @return rolypoly object with regularization results
#' @examples
#' \dontrun{rolypoly_perform_regularized_inference(rolypoly)}
rolypoly_perform_regularized_inference <- function(rolypoly, ...) {

  message('performing cross validation')
  rolypoly$cv_regularized_full_results <- cv_regularized_parameter_estimator(
    vectorize_rolypoly(rolypoly$data), ...
  )

  # add on block values
  rolypoly$regularized_block_values <- calculate_block_values(
    rolypoly$raw_block_data, rolypoly$cv_regularized_full_results$parameters
  )
  # calculate expected block values (accounts for ld and snp error)
  rolypoly$regularized_expected_block_values <-
    calculate_expected_block_values_given_ld(
      rolypoly, rolypoly$regularized_block_values)$expected_block_values

  return(rolypoly)
}
