#' Make example LFM and Iset
#' my aims are
#' i) create simple examples to demonstrate all package functions on
#' ii) to use the examples to debug the input and conversion functions
#'


devtools::load_all()

#' Create example restriction fragment pairs

options(digits=22)
RFpairs <- data.frame(RF1 = c(1,1:18), RF2 = c(2,5,7,11,12,11,14,10,15,16,15,19,20,14,17,19,17,20,20))
RFpairs
RFpairsLarge <- data.frame(RF1 = c(9000001,1:18), RF2 = c(6000000000005,9e+09,7,11,12,11,14,10,15,16,15,19,20,14,17,19,17,20,20))
RFpairsLarge
RFpairs3Cols <- data.frame(RF1 = c(1,1:18), RF2 = c(2,5,7,11,12,11,14,10,15,16,15,19,20,14,17,19,17,20,20), extra = c(1:19))
RFpairs3Cols

#' Write the RFpairs to file

setwd(system.file('extdata', package = 'chromium'))

write.table(RFpairs, file = 'exRFpairs.raf', row.names = F, col.names = F, sep = "\t")
write.table(RFpairsLarge, file = 'largeRFpairs.raf', row.names = F, col.names = F, sep = "\t")
write.table(RFpairs, file = 'characterRFpairs.raf', row.names = F, col.names = T, sep = "\t")
write.table(RFpairs3Cols, file = '3colsRFpairs.raf', row.names = F, col.names = F, sep = "\t")
system("cut -f1 exRFpairs.raf > tmp1")
system("cut -f2 exRFpairs.raf | head -n5 > tmp2")
system("paste tmp1 tmp2 > missingRFpairs.raf")
system("rm tmp1 tmp2")

#' Create example genomic annotation of the restriction fragment pairs

chr <- as.character(c(rep.int(1,10),rep.int(2,10),rep.int(3,10)))
start <- rep.int(seq(1,30,by=3),3)
end <- start+2
RF_id <- c(1:30)
RFanno <- data.frame(chr,start,end,RF_id)
RFanno$chr <- as.character(RFanno$chr)

#' Write RFanno into a bed file

write.table(RFanno, file = 'exRFanno.bed', row.names = F, col.names = F, sep = "\t")

# #' Save RFanno and RFpairs in R object to use for testing purposes
#
# save(RFanno, RFpairs, file = "RFanno_RFpairs.rda")

#' From the annotation and the restriction fragment pairs, create an LFM (ligation frequency matrix)

exampleLFM <- Create_any_resolution_LFM(RFpairs,RFanno)
exampleLFM


#' Import RFpairs and RFanno from file into an InteractionSet

exampleIset <- import_chrom(bed = 'exRFanno.bed',
                                 workDir = system.file('extdata', package = 'chromium'),
                                 raflist = list("exRFpairs.raf"),
                                 binned = T)
interactions(exampleIset)

#' so far so good
#' Convert imported Iset into an LFM and check if it's the same as the LFM created directly from the pairs.

converted <- Iset_to_LFM(exampleIset)
all.equal(converted[[1]], exampleLFM)
all.equal(converted[[2]], RFanno)
all.equal(converted[[3]], RFpairs)

#' Now convert LFM into Iset in order to check that it works correctly

converted2 <- LFM_to_Iset(exampleLFM, binned = T)
all.equal(converted2, exampleIset)

#' Great, conversion works correctly on toy example.
#'
#' Save the data in their respective folders: The example data will be saved in data:

devtools::use_data(exampleLFM, exampleIset, overwrite = TRUE)

#' The example raw data (i.e. RFanno and RFpairs) is already saved in inst/extdata


###################
### old example ###
###################

# # import the example data
# Example_Iset <- chromium::import_chrom(bed = "RFanno_HindIII.bed",
#                                workDir = system.file("extdata", package = "chromium"),
#                                binned = F)
# InteractionSet::interactions(Example_Iset)
#
# # convert into a list of LFM, RFpairs and RFanno
# convertedList <- chromium::Iset_to_LFM(Example_Iset)
# str(convertedList)
#
# # converted list should equal list imported to LFM
# origList <- chromium:::read_bed_raf(bed = "RFanno_HindIII.bed",
#                                     workDir = system.file("extdata", package = "chromium"))
# str(origList)
#
# chr2_Iset <- chromium::subset_chrom(Example_Iset, chr = 2)
#
# small_chr2 <- chr2_Iset[1000000:1000010]
#
# InteractionSet::interactions(small_chr2)
#
# small_chr2_LFM <- chromium::Iset_to_LFM(small_chr2)[[1]]
#
# small_chr2_LFM
