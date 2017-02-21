library("InteractionSet")

context("Import of an InteractionSet")
# I am suppressing warnings because R 3.3.2 gives a warning if you test dataframes like "any(df < 0)" or
# "sum( df < 0 )" if the contents of the dataframe are numeric and not integer.

test_that("read_RFanno_and_RFpairs_do_not_change", {

  # import bed and raf files to make RFanno and RFpairs dataframes
  wd <- system.file("extdata", package = "chromium")
  bd <- "exRFanno.bed"
  rl <- c('exRFpairs.raf', 'largeRFpairs.raf')
  mylist <- suppressWarnings(read_bed_raf(bed = bd, raflist = rl, workDir = wd))

  expect_equal_to_reference(mylist, file = 'RFannoRFpairsList.rds')
})

test_that("string_RFpairs_give_error", {

  wd <- system.file("extdata", package = "chromium")
  bd <- "exRFanno.bed"

  rl <- c('characterRFpairs.raf')
  expect_error(suppressWarnings(read_bed_raf(bed = bd, raflist = rl, workDir = wd)), regexp = "(scan() expected 'a real', got)*")

  rl <- c('3colsRFpairs.raf')
  expect_error(suppressWarnings(read_bed_raf(bed = bd, raflist = rl, workDir = wd)), regexp = '*(more than two colums)*')

  rl <- c('missingRFpairs.raf')
  expect_error(suppressWarnings(read_bed_raf(bed = bd, raflist = rl, workDir = wd)), regexp = "*(contains missing restriction fragment IDs)*")

})

## NEED TO change something here. expect_error is not really checking the error message
