#' Make example LFM and Iset

Example_Iset <- chromium::import_chrom(bed = "RFanno_HindIII.bed",
                               workDir = system.file("extdata", package = "chromium"),
                               binned = F)
InteractionSet::interactions(Example_Iset)

chr2_Iset <- chromium::subset_chrom(Example_Iset, chr = 2)

small_chr2 <- chr2_Iset[1000000:1000010]

InteractionSet::interactions(small_chr2)

small_chr2_LFM <- chromium::Iset_to_LFM(small_chr2)[[1]]

small_chr2_LFM
