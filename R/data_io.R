#' Functions for opening and organizing data.
#'
#' We must open block annotation data, gwas data, snp annotations, gene annotations.
#' Here, you'll find functions that to this and organize these data into a rolypoly object.
#' @name data_io
#' @import data.table
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom dplyr row_number
NULL


#' Load annotations for blocks of LD, in some cases this is a gene annotation
#' with a window around a gene.
#'
#' @param rolypoly rolypoly data object
#' @param block_annotation annotation information for block
#' @param genes if these are genes
#' @return rolypoly data with block annotations attached
#' @examples \dontrun{rolypoly_load_block_annotation(rolypoly, block_annotation)}
rolypoly_load_block_annotation <- function(rolypoly, block_annotation, genes = T) {
  message('adding block annotations')
  if (class(block_annotation) != 'data.frame') { stop('require data frame of annotations') }
  if (all(!(c('chrom', 'start', 'end', 'label') %in% colnames(block_annotation)))) {
    stop('require chrom, start, and end, cols')
  }
  if (any(duplicated(block_annotation$label))) { stop('remove duplicate block labels and rerun') }
  # sort and add a random partition label for bootstrap
  rolypoly$blocks <- block_annotation %>% dplyr::arrange(chrom, start) %>%
    dplyr::mutate(partition = cut(row_number(), breaks = 100, labels = F))
  return(rolypoly)
}

#' Load gwas data
#'
#' @param rolypoly rolypoly data
#' @param gwas_data gwas data
#' @param snp_annotations if there are additional snp annotations included
#' @param gwas_z_filter if we want to remove large effect SNPs
#' @param add_spline for fitting a spline to maf
#' @param n_knots number of knots for spline
#' @param add_poly for fitting a polynomial to maf
#' @param n_degree degree of polynomial to fit
#' @return rolypoly object with gwas data loaded
#' @examples \dontrun{rolypoly_load_gwas(rolypoly, gwas_data)}
rolypoly_load_gwas <- function(rolypoly, gwas_data, snp_annotations = NULL,
                               gwas_z_filter = -1, add_spline = F, n_knots = 1,
                               add_poly = F, n_degree = 2) {
  message('adding gwas')
  if (class(gwas_data)[1] != 'data.frame') { stop('require data frame of gwas')}
  necessary_cols <- c('chrom', 'pos', 'rsid', 'beta', 'se', 'maf') # try to control for maf
  if (all(!(necessary_cols %in% colnames(gwas_data)))) { stop('require chrom, pos, beta, and se') }
  if (!is.null(snp_annotations)) {
    if (all(snp_annotations %in% colnames(gwas_data))) {
      message(paste('including the following snp annotations', snp_annotations, sep = ': '))
      necessary_cols <- c(necessary_cols, snp_annotations)
    } else { stop('could not find provided snp annotations') }
  }
  if (!is.numeric(gwas_z_filter) | !is.integer(gwas_z_filter)) { gwas_z_filter <- as.numeric(gwas_z_filter) }
  if (grepl('chr', gwas_data$chrom[1])) { stop('remove chr from chrom and rerun') }
  if ('maf' %in% colnames(gwas_data)) {
    message('filtering out SNPs with MAF < 1%')
    gwas_data <- gwas_data %>%
      dplyr::mutate(maf = ifelse(test = maf > 0.5, yes = 1 - maf, no = maf)) %>%
      dplyr::filter(maf > 0.01)
  }
  if (gwas_z_filter > 0) {
    window_size <- 1e5
    message(paste('removing SNPs with |z| > ', gwas_z_filter, sep = ''))
    gwas_data <- gwas_data %>% dplyr::mutate(abs_z = abs(beta/se)) %>%
      dplyr::filter(abs_z < gwas_z_filter)

    #     ### uncomment if we want to remove a window around a SNP
    #     # this should remove large SNPs and all SNPs within a window region around that SNP
    #     gwas_data <- lapply(1:22, function(chrom) {
    #       chrom_data <- gwas_data[gwas_data$chrom == chrom, ]
    #       large_snp_pos <- (chrom_data %>% dplyr::filter(abs_z > gwas_z_filter))$pos
    #       message(paste('removing ', length(large_snp_pos), ' snps from chromosome ',
    #                     chrom, ', and regions around', sep = ''))
    #       sapply(large_snp_pos, function(big_snp_pos) {
    #         chrom_data <<- chrom_data %>% dplyr::filter(pos > (big_snp_pos - (window_size/2)),
    #                                                     pos < (big_snp_pos + (window_size/2)))
    #       })
    #       return(chrom_data)
    #     }) %>% data.table::rbindlist() %>% dplyr::select(-abs_z) %>% data.frame
  }
  if (add_spline) {
    spline_df <- n_knots + 3
    # save spline function for later interpretation
    rolypoly$spline_function <- function(maf) {
      splines::bs(maf, knots = exp(seq(log(1e-2), log(0.45), length.out = n_knots)),
                  Boundary.knots = c(1e-3, 0.5), degree = 3)
    }
    maf_splines <- data.frame(rolypoly$spline_function(gwas_data$maf))
    colnames(maf_splines) <- paste('maf_spline', 1:spline_df, sep = '_')
    necessary_cols <- c(necessary_cols, colnames(maf_splines))
    rolypoly$snp_annotations <- colnames(maf_splines)
    gwas_data <- cbind(gwas_data, maf_splines)
  }
  if (add_poly) {
    rolypoly$poly_function <- function(maf) { stats::poly(maf, degree = n_degree, raw = T) }
    maf_poly <- data.frame(rolypoly$poly_function(gwas_data$maf))
    colnames(maf_poly) <- paste('maf_poly', 1:n_degree, sep = '_')
    necessary_cols <- c(necessary_cols, colnames(maf_poly))
    rolypoly$snp_annotations <- colnames(maf_poly)
    gwas_data <- cbind(gwas_data, maf_poly)
  }
  # convert to necessary data structure, ande only use
  gwas_data <- data.table::data.table(gwas_data[gwas_data$chrom %in% 1:22, necessary_cols]) %>%
    dplyr::mutate(chrom = as.numeric(chrom))

  # make sure we do not have factors for the rsid
  # using character decreases memory footprint enormously when copied to each gene
  gwas_data$rsid <- as.character(gwas_data$rsid)

  rolypoly$gwas <- gwas_data
  rolypoly$snp_annotations <- c(rolypoly$snp_annotations, snp_annotations)
  return(rolypoly)
}

#' Link blocks and gwas
#'
#' Takes block information, potentially independent LD blocks or gene blocks,
#' and gwas data and organizes the data for interenal processing
#'
#' @param rolypoly a rolypoly object
#' @param ld_folder path to a folder with ld data
#' @param r2_threshold LD threshold to look at data
#' @param run_parallel check if user wants to run in parallel
#' @return rolypoly object with data attached
#' @examples
#' \dontrun{rolypoly_link_blocks_and_gwas(rolypoly, ld_folder, r2_threshold)}
rolypoly_link_blocks_and_gwas <- function(rolypoly, ld_folder, r2_threshold = .2,
                                          run_parallel = F) {

  # nested lists of chrom and then individual blocks
  chrom_block_list <- lapply(split(rolypoly$blocks, f = rolypoly$blocks$chrom),
                             function(chrom) {split(chrom, f = chrom$label)})
  chrom_gwas_list <- lapply( # make sure we have this key set
    split(rolypoly$gwas, f = rolypoly$gwas$chrom), function(gwas) {
      gwas <- data.table::data.table(gwas)
      data.table::setkey(gwas, pos)
      return(gwas)
    })
  # prevent naming issues and indexing issues
  names(chrom_gwas_list) <- paste('chr', names(chrom_gwas_list), sep = '')

  use_TSS_info <- F
  # adding tss distance and maf as snp annotations
  if ('TSS' %in% rolypoly$blocks) {
    use_TSS_info <- T
    rolypoly$snp_annotations <- # appended if we have splines
      c(rolypoly$snp_annotations, 'upstream_tss_dist','downstream_tss_dist')
  }

  # run things in parallel if user specified
  "%loop_function%" <- `%do%`
  if (run_parallel) { "%loop_function%" <- `%dopar%` }

  message('beginning processs of linking gwas and block annotations')
  data <-
    foreach(
      chrom_block = chrom_block_list, .combine = 'c', .inorder = F,
      .verbose = F, .packages = c('rolypoly')
    ) %loop_function% {
      chrom <- chrom_block[[1]][1,]$chrom
      chrom_idx <- paste('chr', chrom, sep = '')
      chrom_ld_file_path <- paste(ld_folder, '/', chrom, '.Rds', sep = '')
      ld_data <- readRDS(chrom_ld_file_path)[(R**2 > r2_threshold), .(SNP_A, SNP_B, R)]
      setkey(ld_data, SNP_A) # just in case
      message(paste(utils::timestamp(quiet = T), ' - starting blocks on chrom: ', chrom, sep = ''))
      if (!(chrom_idx %in% names(chrom_gwas_list))) {
        warning(paste(chrom_idx, ' missing, could be problem', sep = ''))
        return(NULL)
      }
      if (is.null(chrom_gwas_list[[chrom_idx]])) {
        warning(paste(chrom_idx, ' data missing, could be problem', sep = ''))
        return(NULL)
      }
      chrom_data <- lapply(chrom_block, function(block) {
        # I think this is by data.table fast search convention
        snp_data <- chrom_gwas_list[[chrom_idx]][.(block$start:block$end), nomatch = 0L]
        # actually set this up properly, for now we're adding a 4th degree polynomial
        # might want to not include the absolute value here...
        if (use_TSS_info) {
          snp_data <- snp_data %>% dplyr::mutate(
            # eventually change this over to a spline as well
            tss_dist = (pos - block$TSS)/5e3,
            upstream_tss_dist = ifelse(test = tss_dist <= 0, yes = abs(tss_dist), no = 0),
            downstream_tss_dist = ifelse(test = tss_dist > 0, yes = tss_dist, no = 0)
          )
        }
        if (nrow(snp_data) == 0) { return(NULL) }
        beta_squared <- snp_data$beta^2
        # again, pretty sure this is binarize data table lookup
        sub_ld <- ld_data[.(snp_data$rsid), nomatch = 0L]
        if (nrow(sub_ld) == 0) {ld_matrix <- diag(1, nrow = nrow(snp_data))}
        else {ld_matrix <- make_ld_matrix(snp_data$rsid, sub_ld)}

        block_data <- list(
          block_info = block,
          snps = snp_data,
          ld_data = sub_ld,
          n_snps = nrow(snp_data),
          gwas_block_score = sum(beta_squared),
          # we probably do not need to save the plain ld matrix
          # ld_matrix = ld_matrix,
          # ld_matrix_squared = ld_matrix * ld_matrix,
          # using sparse matrices woop
          ld_matrix_squared = Matrix::Matrix(ld_matrix * ld_matrix),
          partition = block$partition,
          y = beta_squared
        )
        return(block_data)
      })
      return(chrom_data)
    }
  rolypoly$data <- data[!sapply(data, is.null)]
  names(data) <- sapply(data, function(block) {block$block_info$label})
  rolypoly$gwas <- rolypoly$gwas[, !grepl('maf_spline|maf_poly', colnames(rolypoly$gwas))]
  return(rolypoly)
}

#' Block annotations, usually gene model.
#'
#' @param rolypoly a rolypoly object
#' @param block_data a data frame of block information, usually gene expression.
#'  Requires rownames that are identitcal to block labels loaded previously.
#' @return a rolypoly object with block information loaded
#' @examples
#' \dontrun{rolypoly_load_block_data(rolypoly, block_data)}
rolypoly_load_block_data <- function(rolypoly, block_data) {
  if (is.null(rolypoly$data)) {
    message('no loaded block data, not linking functional data')
    return(rolypoly)
  }
  if (class(block_data) != 'data.frame') {block_data <- data.frame(block_data)}
  message('merging functional information about blocks')
  block_data <- data.frame(block_data) # important to get nrow, could be better
  data <- lapply(rolypoly$data, function(block) {
    # now we should add the block data, and LD transform the data
    if (!(block$block_info$label %in% rownames(block_data))) {
      block$include_in_inference <- F
      block$x <- NULL # to make sure we totally replace previous stuffs
      return(block)
    }
    x <- block_data[block$block_info$label,]
    if (nrow(x) != 1) { stop('remove duplicates from block data') }
    x <- as.matrix(x[rep(1,block$n_snps),])
    if (!is.null(rolypoly$snp_annotations)) {
      # snp information data table
      x <- cbind(as.matrix(block$snps[,rolypoly$snp_annotations, with = F]), x)
    }
    colnames(x) <- c(rolypoly$snp_annotations, colnames(block_data))

    # consider not using the as.matrix such that we can use sparse matrices
    # with open chromatin data
    block$x <- as.matrix(block$ld_matrix_squared %*% x)
    block$include_in_inference <- T
    return(block)
  })
  rolypoly$data <- data[!sapply(data, is.null)]
  # save later for block value prediction, especially bootstraps
  rolypoly$raw_block_data <- block_data
  return(rolypoly)
}

#' Add LD corrected block scores to rolypoly.
#'
#' @param rolypoly rolypoly data
#' @param fast_calculation if F then LD deconvolution else quadratic form.
#' @return rolypoly object with LD corrected gwas scores attached
#' @export
#' @examples \dontrun{rolypoly_add_ld_corrected_gwas_block_scores(rolypoly)}
rolypoly_add_ld_corrected_gwas_block_scores <- function(rolypoly, fast_calculation = T) {

  # warning('this function is still extremely slow')
  # stop('figure out ld_matrix from ld_matrix squared still')

  # no fallback for now
  if (!requireNamespace("CompQuadForm", quietly = TRUE)) {
    stop("CompQuadForm needed for this function to work. Please install it.",
         call. = FALSE)
  }

  g <- 0
  rolypoly$data <- lapply(rolypoly$data, function(block) {
    g <<- g + 1
    if ( (g %% 1e3) == 0) { message(paste('another k genes', g, sep = ': ')) }

    ######## THIS IS CODE IF WE USE THE COMPQUADFORM PACKAGE
    # # initialize values
    raw_score <- sum((block$snps$beta / block$snps$se)**2)
    # null matrix is just the ld, rounding because numerical
    # lambda <- eigen(block$ld_matrix, symmetric = T, only.values = T)$values
    lambda <- eigen(round(block$ld_matrix, 10), symmetric = T, only.values = T)$values
    # get corrected p value
    p_val <- CompQuadForm::imhof(raw_score, lambda = lambda, epsabs = 1e-15, epsrel = 1e-15)$Qq
    # give min p-val if too small, I wouldn't trust p-values below 1e-15
    if (p_val <= 1e-15) {
      block$gene_score_message <- 'p-value below 1e-15 threshold'
      p_val <- NA
    } else if (p_val > 1) {
      block$gene_score_message <- 'p-value above 1'
      p_val <- NA
    }

    #     ######## THIS IS CODE IF WE USE LD DECONVOLUTION METHOD, TOO SLOW for now
    #     # rounding because numerical
    #     ld_matrix <- round(block$ld_matrix, 10)
    #     # check if ld matrix is psd
    #     if (!matrixcalc::is.positive.definite(ld_matrix)) {
    #       # if it is not then we can calculate the nearest psd version
    #       ld_matrix <- (Matrix::nearPD(ld_matrix, corr = T, keepDiag = T,
    #                                    ensureSymmetry = F, doSym = T))$mat
    #     }
    #     # then we can calculate deconvolved score from z-scores
    #     zs <- (block$snps$beta / block$snps$se)**2
    #     # parens make faster
    #     block_score <- zs %*% (solve(ld_matrix) %*% zs)
    #     # then compute p_value
    #     block$corrected_gwas_block_score_pval <- stats::rchisq(block_score, df = block$n_snps)


    block$corrected_gwas_block_score_pval <- p_val
    return(block)
  })
  return(rolypoly)
}

#' Helper function to pull LD data from NCBI.
#'
#' Given the path of a gwas file open it into a data.table object
#'
#' @param all_snps The snps that were queried
#' @param ld_data A returned LD matrix with SNP, Proxy, and RSquared columns
#' @return an LD matrix where query snps will be the first columns in the correct order
#' @examples \dontrun{make_ld_matrix(all_snps, ld_data)}
make_ld_matrix <- function(all_snps, ld_data) {

  # the same snp shouldn't be have more than one entry
  mat_dim <- length(all_snps)
  ld_matrix <- diag(mat_dim)

  # handles singleton genes
  if (mat_dim == 1) { return(as.matrix(1)) } # no need to check psd
  if (mat_dim == 2) {
    r_value <- as.numeric(ld_data[(SNP_A == all_snps[1]) &
                                    (SNP_B == all_snps[2])]$R)

    if (length(r_value) == 1) { # spliting up ifs allowed to go through...
      if (!is.na(r_value)) {
        # keep symmetric
        ld_matrix[1, 2] <- r_value
        ld_matrix[2, 1] <- r_value
      }
    } # everything else stays zero
  } else {
    # figure out how to do this faster
    sapply(1:mat_dim, function(i){
      sapply((i+1):mat_dim, function(j) {
        r_value <- as.numeric(ld_data[(SNP_A == all_snps[i]) &
                                        (SNP_B == all_snps[j])]$R)
        if (length(r_value) == 1) {
          if (!is.na(r_value)) {
            # keep symmetric
            ld_matrix[i, j] <<- r_value
            ld_matrix[j, i] <<- r_value
          }
        } # everything else stays zero
      })
    })
  }

  return(ld_matrix)
}
