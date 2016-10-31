#' Make example LFM and Iset

library(InteractionSet)

# # import the exmple data
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

###########################################
# Rather make toy example than a real one #
###########################################
# my aims are
# i) create simple examples to demonstrate all package functions on
# ii) to use the examples to debug the input and conversion functions

# create restriction fragment pairs
RFpairs <- data.frame(RF1 = c(1,1:18), RF2 = c(2,5,7,11,12,11,14,10,15,16,15,19,20,14,17,19,17,20,20))

# write the RFpairs to file
write.table(RFpairs, file = file.path(system.file('extdata', package = 'chromium'), 'exRFpairs.raf'), row.names = F, col.names = F, sep = "\t")

# genomic annotation of the restriction fragment pairs
chr <- as.character(c(rep.int(1,10),rep.int(2,10),rep.int(3,10)))
start <- rep.int(seq(1,30,by=3),3)
end <- start+2
RF_id <- c(1:30)
RFanno <- data.frame(chr,start,end,RF_id)
RFanno$chr <- as.character(exRFanno$chr)

# write RFanno into a bed file
write.table(RFanno, file = file.path(system.file('extdata', package = 'chromium'), 'exRFanno.bed'), row.names = F, col.names = F, sep = "\t")

# Make an LFM
toyLFM <- Create_any_resolution_LFM(RFpairs,RFanno)
toyLFM

# import RFpairs and RFanno from file

toyIset_imported <- import_chrom(bed = 'exRFanno.bed',
                                 workDir = system.file('extdata', package = 'chromium'),
                                 raflist = list("exRFpairs.raf"),
                                 binned = T)
interactions(toyIset_imported) # so far so good

