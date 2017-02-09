library("InteractionSet")

context("Conversion between matrix and InteractionSet")

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


test_that("conversion_from_Iset_to_LFM_works_on_example_data", {

  load(file.path('data', 'exampleLFM.rda'))
  load(file.path('data', 'exampleIset.rda'))

  RFanno <- read.table(file = system.file("extdata", "exRFanno.bed", package = "chromium"))
  names(RFanno) <- c('chr', 'start', 'end', 'RF_id')
  RFanno$chr <- as.character(RFanno$chr)

  RFpairs <- read.table(file = system.file("extdata", "exRFpairs.raf", package = "chromium"))
  names(RFpairs) <- c('RF1', 'RF2')

  converted <- Iset_to_LFM(exampleIset)
  expect_true(all.equal(converted[[1]], exampleLFM))
  expect_true(all.equal(converted[[2]], RFanno))
  expect_true(all.equal(converted[[3]], RFpairs))
})

test_that("conversion_from_LFM_to_Iset_works_on_example_data", {

  # Since saving the example Iset for some reason changes the assays attribute,
  # I'm creating an Iset from scratch to compare the converted Iset against.

  load(file.path('data', 'exampleLFM.rda'))
  converted <- LFM_to_Iset(LFM = exampleLFM, binned = TRUE)

  start <- rep.int(seq(1,30,by=3),3)
  allRegions <- GRanges(as.character(c(rep.int(1,10),rep.int(2,10),rep.int(3,10))),
                        IRanges(start, start + 2),
                        IDs = c(1:30))
  gi <- GInteractions(anchor1 = c(1,1:18),
                      anchor2 = c(2,5,7,11,12,11,14,10,15,16,15,19,20,14,17,19,17,20,20),
                      allRegions)
  iset <- InteractionSet(matrix(rep_len(1, 19)), gi, metadata = list(binSize = 3))

  expect_true(all.equal(converted, iset))
})
