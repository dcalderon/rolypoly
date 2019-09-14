#' Plot functions.
#'
#' Functions for rolypoly results.
#' @name rolypoly_plots
#' @importFrom ggplot2 ggplot aes geom_point theme_bw xlab ylab geom_hline
#'  theme element_text geom_pointrange
NULL

#' Visualize parameter estimates after running inference
#'
#' @param rolypoly a rolypoly object
#' @return ggplot2 object
#' @export
#' @examples
#' \dontrun{plot_rolypoly_annotation_estimates(rolypoly)}
plot_rolypoly_annotation_estimates <- function(rolypoly) {
  df <- rolypoly$bootstrap_results %>%
    merge(data.frame(annotation = names(rolypoly$full_results$parameters),
                     effect = rolypoly$full_results$parameters)) %>%
    dplyr::mutate(fbt_value = effect / bootstrap_error)
  df$fdr <- stats::p.adjust(df$bp_value, method = 'BH')

  gene_annotations <- colnames(rolypoly$raw_block_data)
  p <- df %>% dplyr::filter(annotation %in% gene_annotations) %>% dplyr::arrange(-bt_value) %>%
    ggplot(aes(x = factor(annotation, levels = annotation),
                        y = effect, ymin = CI_lo, ymax = CI_hi)) +
    geom_pointrange(color = 'darkgray', fill = 'black', shape = 21) +
    theme_bw() + xlab('') + ylab('RolyPoly estimated effect') +
    theme(text = element_text(size = 14), #, family = 'GaramondNo8'),
          axis.text.x = element_text(angle = 55, hjust = 1, size = 14)) +
    geom_hline(yintercept = 0, linetype = 3)
  return(p)
}

#' Rank annotations by p-value after running inference
#'
#' @param rolypoly a rolypoly object
#' @return ggplot2 object
#' @export
#' @examples
#' \dontrun{plot_rolypoly_annotation_ranking(rolypoly)}
plot_rolypoly_annotation_ranking <- function(rolypoly) {
  df <- rolypoly$bootstrap_results %>%
    merge(data.frame(annotation = names(rolypoly$full_results$parameters),
                     effect = rolypoly$full_results$parameters)) %>%
    dplyr::mutate(fbt_value = effect / bootstrap_error)
  df$fdr <- stats::p.adjust(df$bp_value, method = 'BH')

  gene_annotations <- colnames(rolypoly$raw_block_data)
  p <- df %>% dplyr::filter(annotation %in% gene_annotations) %>% dplyr::arrange(-bt_value) %>%
    ggplot(aes(x = factor(annotation, levels = annotation), y = -log10(bp_value))) +
    geom_point() + theme_bw() + xlab('') + ylab('-log10(p value)') +
    theme(text = element_text(size = 14), #, family = 'GaramondNo8'),
                       axis.text.x = element_text(angle = 55, hjust = 1, size = 14)) +
    geom_hline(yintercept = -log10(0.05), linetype = 2)
  return(p)
}
