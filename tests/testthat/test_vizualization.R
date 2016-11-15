
context("Triangle vizualization")

test_that("specifiying_a_genomic_model_but_no_genome_identifyer_works", {
  # chr2 <- import_chrom(bed = 'RFanno_HindIII.bed',
  #                      raflist = list('chr2.raf'),
  #                      workDir = file.path(system.file("extdata",
  #                                                      package = "chromium")))

  data("chr2")

  chr2_binned <- bin_chrom(chr2, binSize = 10000)
  data("ga_mm9_chr2")

  visualize_chrom(chr2_binned, chr = 2,
                  from = 330e4,
                  to = 810e4,
                  geneModels = ga_mm9_chr2)

  expect_true(file.exists("triangle_visualization.pdf"))
  if(file.exists("triangle_visualization.pdf")){
    file.remove("triangle_visualization.pdf")
  }
  })

