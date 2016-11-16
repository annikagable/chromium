library(tibble)

context("Triangle vizualization")

test_that("specifying_a_genomic_model_by_name_works", {

  data("chr2")
  chr2_binned <- bin_chrom(chr2, binSize = 10000)

  visualize_chrom(chr2_binned, chr = 2,
                  from = 300e4,
                  to = 600e4,
                  geneModels = "mm9")

  expect_true(file.exists("triangle_visualization.pdf"))

  if (file.exists("triangle_visualization.pdf")) {
    file.remove("triangle_visualization.pdf")
  }

  visualize_chrom(chr2_binned, chr = 2,
                  from = 300e4,
                  to = 600e4,
                  geneModels = "mmusculus")

  expect_true(file.exists("triangle_visualization.pdf"))

  if (file.exists("triangle_visualization.pdf")) {
    file.remove("triangle_visualization.pdf")
  }

  })

