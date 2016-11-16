library("InteractionSet")

context("Ligation frequency conversion")

test_that("conversion_to_LFM_works_without_interaction_ID", {

  set.seed(1000)
  N <- 5
  all.starts <- round(runif(N, 1, 100))
  all.ends <- all.starts + round(runif(N, 5, 20))
  all.regions <- GRanges(rep(c("chrA", "chrB"), c(N - 3, 3)),
                         IRanges(all.starts, all.ends))

  Nr <- 2
  Nc <- 3
  all.anchor1 <- sample(N, Nr)
  all.anchor2 <- sample(N, Nc)
  counts <- matrix(rpois(Nr*Nc, lambda = 10), Nr, Nc)
  x <- ContactMatrix(counts, all.anchor1, all.anchor2, all.regions)

  # turn into a an InteractionSet

  x_i <- deflate(x)
  expect_true(all.equal(names(Iset_to_LFM(x_i)),c("LFM", "RFanno", "RFpairs")))

  })

