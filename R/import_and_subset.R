
#' Read preprocessed genomic interaction files into an InteractionSet object.
#'
#' @param bed File name of a bed file without header, with columns for chromosome (Ensembl format, no "chr"), start, end, and
#' restriction fragment ID of each restriction fragment.
#' @param workDir The directory where all input files are located. Defaults to the current working directory.
#' @param raflist A list of file names containing the interactions. The files are in raf format, meaning they contain two
#' columns (no header) with the interacting restriction fragment IDs. If raflist=NULL, all raf files in the directory will
#' be read.
#' @param binned Boolean, whether the matrix is already binned. Defaults to FALSE.
#' @return An InteractionSet object containing all interactions from all input raf files, and the bin size (if any) in the metadata.
#' @examples
#' # You can import several specific raf files from the same directory
#' # (the working directory if no workDir is specified).
#' # myInteractions <- import_chrom(bed = "annotation.bed",
#' #                                raflist = list("1.raf", "2.raf"))
#'
#' # Or just choose all raf files within a specific directory
#' myInteractions <- import_chrom(bed = "RFanno_HindIII.bed", workDir = system.file("extdata", package = "chromium"), binned = F)
#'
#' @export import_chrom
#'
import_chrom <- function(bed, raflist = NULL, workDir = getwd(), binned = FALSE){

  list <- read_bed_raf(bed, raflist, workDir)
  RFanno <- list[[1]]
  RFpairs <- list[[2]]

  Iset <- RFpairs_to_Iset(RFanno, RFpairs, binned)
  return(Iset)
}

#' Read preprocessed genomic interaction files to produce RFanno and RFpairs.
#'
#' @param bed File name of a bed file without header, with columns for chromosome (Ensembl format, no "chr"), start, end, and
#' restriction fragment ID of each restriction fragment.
#' @param workDir The directory where all input files are located. Defaults to the current working directory.
#' @param raflist A list of file names containing the interactions. The files are in raf format, meaning they contain two
#' columns (no header) with the interacting restriction fragment IDs. If raflist=NULL, all raf files in the directory will be
#' read.
#' @return A list consisting of an RFanno dataframe containing the restriction fragment annotation, and an RFpairs dataframe
#' containing the interacting restriction fragment pairs.
#' @examples
#' # RFanno <- read_bed_raf(bed = "annotation.bed", raflist = list("1.raf", "2.raf"))[[1]]
#' # RFpairs <- read_bed_raf(bed = "annotation.bed", raflist = list("1.raf", "2.raf"))[[2]]
#'
read_bed_raf = function(bed, raflist = NULL, workDir = getwd()){

  # create restriction fragment annotation
  message(paste0("Reading restriction fragment annotation from ", bed))
  RFanno <- utils::read.delim(file.path(workDir, bed), header = FALSE, sep = "\t")
  if(ncol(RFanno) != 4) stop("Your annotation file needs to have four columns: chr, start, end, and ID.")

  message("Sorting annotation by restriction fragment IDs (4th column).")
  RFanno <- RFanno[ SummarizedExperiment::order(RFanno[ ,4]), ]
  names(RFanno) <- c("chr", "start", "end", "RF_id")
  if(any(diff(RFanno[ ,4]) != 1)) stop("There are restriction fragment IDs missing in the annotation file.")
  if(any(RFanno[ ,2] > RFanno[ ,3])) stop("The are start values that are larger than their end values in the annotation file.")

  # create restriction fragment pairs
  if(is.null(raflist)){
    raflist <- list.files(workDir, pattern = "*.raf")
  }
  message("Reading interactions from your .raf file(s).")
  pairslist <- lapply(raflist, function(file){
    df <- utils::read.delim(file.path(workDir, file), header = FALSE, sep = "\t")
    stopifnot(ncol(df) == 2)
    return(df)
  })
  RFpairs <- do.call(rbind,pairslist)

  message("Sorting interactions by restriction fragment IDs.")
  # put higher ID in second column.
  RF_1 <- pmin(RFpairs[,1], RFpairs[,2])
  RF_2 <- pmax(RFpairs[,1], RFpairs[,2])
  RFpairs <- data.frame(RF_1 = RF_1, RF_2 = RF_2)

  # sort IDs first by first column, then by second column.
  RFpairs <- RFpairs[ SummarizedExperiment::order(RFpairs[,1], RFpairs[,2]), ]

  return(list(RFanno = RFanno, RFpairs = RFpairs))
}

#' Choose a chromosome and a region within which InteractionSet interactions and regions should be retained.
#'
#' Only those interactions and regions of an InteractionSet which lie within the specified range will be retained.
#' This is useful to reduce the object size of the InteractionSet.
#'
#' @param Iset The InteractionSet object to be subsetted.
#' @param chr The chromosome for which interactions and regions should be retained.
#' @param from Optional. The starting point of the interactions and regions to be retained. If NULL, the chromosome
#' beginning is selected.
#' @param to Optional. The ending point of the interactions and regions to be retained. If NULL, the chromosome end
#' is selected.
#' @return An InteractionSet containing only those regions and those interactions specified in the arguments.
#' @examples
#' # chr2_Iset <- subset_chrom(Iset, chr = 2)
#'
#' @export subset_chrom
#'
subset_chrom <- function(Iset, chr, from = NULL, to = NULL){
  if(is.null(S4Vectors::metadata(Iset)$binSize)) binned = FALSE else binned = TRUE
  subLFM <- Iset_region_to_LFM(Iset, chr, from, to)
  subIset <- LFM_to_Iset(subLFM, binned)
  S4Vectors::metadata(subIset)$binSize <- S4Vectors::metadata(Iset)$binSize
  return(subIset)
}



#' Choose a chromosome and a range within which InteractionSet interactions should be retained.
#'
#' The interactions of the InteractionSet will be subsetted, but all regions (i.e. the annotation) will be retained.
#' This is only useful if you want to add interactions in the other regions later. Otherwise, especially for direct
#' conversion to matrices, it is more useful to subset the whole InteractionSet, including the regions.
#'
#' @param Iset The InteractionSet object to be subsetted.
#' @param chr The chromosome for which interactions should be retained.
#' @param from Optional. The starting point of the interactions to be retained. If NULL, the chromosome beginning is selected.
#' @param to Optional. The ending point of the interactions to be retained. If NULL, the chromosome end is selected.
#' @return An InteractionSet object containing only those interactions specified, but all original genomic annotation.
#' @examples
#' # chr2_Iset <- subset_interactions(Iset, chr=2)
#'

subset_interactions <- function(Iset, chr, from = NULL, to = NULL){

  # maximal entry on the selected chromsome
  max <- max(SummarizedExperiment::end(InteractionSet::regions(Iset)[GenomicRanges::seqnames(InteractionSet::regions(Iset)) == chr]))

  if(is.null(from)) from <- 0
  if(is.null(to)) to <- max

  # input checking
  if(from < 0 || to > max) stop("Subset out of bounds. Provide a region within the range of the chromosome")
  if(from > to) stop ('"from" has to be smaller than "to"')

  # subsetting
  indlist <- lapply(InteractionSet::anchors(Iset), function(i){
    start <- methods::as(start(GenomicRanges::ranges(i)), "vector")
    chrom <- methods::as(GenomicRanges::seqnames(i), "vector")
    start >= from & start <= to & chrom == chr
  })
  indices <- which(indlist$first & indlist$second)
  subset <- Iset[indices]
  if(length(subset) < 1) warning("Your subset does not contain any interactions.", call. = FALSE)

  return(subset)
}

