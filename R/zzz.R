# hack to avoid notes about data.table or ggplot not recognizing bindings
# http://stackoverflow.com/questions/9439256/
# how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
globalVariables(c('SNP_A', 'SNP_B', 'bootstrap_estimate', 'bootstrap_error',
                  'bt_value', 'effect', 'bootstrap_error', 'annotation',
                  'bp_value', 'pos', 'R', 'tss_dist', 'chrom', 'maf', 'se',
                  'abs_z', 'chrom_block', '.', 'i', 'CI_hi', 'CI_lo',
                  'geom_pointrange', 'start'))

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to rolypoly! Say hi @diegoisworking.")
}

.onLoad <- function(libname, pkgname) {
  op <- options()
  op.devtools <- list(
    devtools.path = "~/R-dev",
    devtools.install.args = "",
    devtools.name = "Diego",
    devtools.desc.author = '"Diego Calderon <dcal@stanford.edu> [aut, cre]"',
    devtools.desc.license = "GPL-3",
    devtools.desc.suggests = NULL,
    devtools.desc = list()
  )
  toset <- !(names(op.devtools) %in% names(op))
  if(any(toset)) options(op.devtools[toset])
}
