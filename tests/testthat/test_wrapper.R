context("Working wrapper")

test_that("wrapper works", {
  rp <- rolypoly_roll(
    gwas_data = sim_gwas_data,
    block_annotation = sim_block_annotation,
    block_data = sim_expression_data_normalized,
    ld_folder = system.file("extdata", "example_ld", package = "rolypoly"),
    bootstrap_iters = 2
  )
  expect_equal(class(rp), 'rolypoly')
})
