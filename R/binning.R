#' Bin an existing InteractionSet
#'
#' For any InteractionSet object, bin the regions and the interactions to make them more feasible to visualize.
#'
#' @param Iset An InteractionSet object.
#' @param binSize The desired bin size.
#' @return The binned InteractionSet, with the bin size stored in metadata(myInteractionSet)$binSize
#'
#' @export bin_chrom

bin_chrom <- function(Iset, binSize){
  list <- Iset_to_LFM(Iset)
  RFanno <- list[[2]]
  RFpairs <- list[[3]]
  LFM <- Create_any_resolution_LFM(RFpairs, RFanno, binSize)
  ifelse(is.null(binSize), binned <-  F, binned <-  T)
  outIset <- LFM_to_Iset(LFM, binned = binned)
  S4Vectors::metadata(outIset)$binSize <- binSize
  return(outIset)
}

#' Create a ligation frequency matrix or InteractionSet binned at any resolution.
#'
#' Use this function to create a sparse ligation frequency matrix with genome annotation or an InteractionSet.
#' Binning is optional. By default, the restriction fragment (RF) resolution will be chosen.
#'
#' @param RFpairs A dataframe with two columns containing restriction fragment pairs. The restriction
#' fragment with the lower ID is always in the left column, and all restriction fragments are sorted
#' in ascending order first by the left column and then by the right colunn.
#' @param RFanno A dataframe giving the genomic annotation of the restriction fragments. Columns are:
#' chromosome, start of restriction fragment, end of restriction fragment, restriction fragment ID.
#' The table is sorted by chromosomal position and restriction fragments are numbered along their
#' genomic position.
#' @param binSize The size of a bin. Optional parameter; if no value is given, the function will default
#' to not binning (binSize=NULL) and instead giving the ligation frequency matrix at restriction fragment
#' resolution (calculation can be slow).
#' @return The ligation frequency matrix \code{LFM}, a sparse matrix annotated by genomic position,
#' or an InteractionSet.
#' @examples
#' library(chromium)
#'
#' # Restriction fragment pairs are given here (RFpairs)
#' toyRFpairs <- data.frame(RF1 = c(1:15), RF2 = c(5,7,11,12,11,14,10,15,16,15,19,20,14,17,19))
#'
#' # Genomic annotation of the restriction fragment pairs (RFanno)
#' chr <- c(rep.int(1,10), rep.int(2,10), rep.int(3,10))
#' start <- rep.int(seq(1,30,by = 3), 3)
#' end <- start + 2
#' RF_id <- c(1:30)
#' toyRFanno <- data.frame(chr, start, end, RF_id)
#'
#' # Create a toy ligation frequency matrix
#' toyMatrix <- Create_any_resolution_LFM(toyRFpairs, toyRFanno, binSize = 3)
#'
#' # In a more realistic example the bin size would range between kilobases and megabases
#'

Create_any_resolution_LFM <- function(RFpairs, RFanno, binSize = NULL){

  # Input checking, make sure that the input dataframes are ordered
  if(!exists("RFpairs")) stop("Please provide RFpairs.")
  if(!exists("RFanno")) stop("Please provide RFanno.")
  if(!ncol(RFpairs) == 2) stop("RFpairs has to have exactly two columns.")
  if(!all(RFpairs >= 0)) stop("All RFpairs entries have to be larger than or equal to zero.")
  if(!ncol(RFanno) >= 3) stop("RFanno requires at least three columns: chr, start, end.")
  if(!all(RFpairs[ ,2] <= nrow(RFanno))) stop("You have interactions not supported by RFanno.")
  if(!all(RFpairs[ ,1] <= RFpairs[, 2])) stop("All left hand mates of RFpairs must be 5' with respect to their right hand mates.")
  if(ncol(RFanno) >= 4) if(!all(diff(RFanno[,4]) == 1 )) stop("Your RFanno IDs are not incrementing by 1:
                                                            Gaps in RFanno detected.")


  if(!is.null(binSize)){
    if(!is.numeric(binSize)) stop("Bin size has to be numeric.")
    if(!binSize %% 1 == 0) stop("Please provide an integer bin size.")

    #Now link the non-numeric chromosomes to a number (e.g. mouse: chr X = 20, chr Y = 21)
    rep.mat <- data.frame(chr = unique(RFanno[ ,1]), repl = seq_along(unique(RFanno[ ,1])), stringsAsFactors=FALSE)

    #bin is a matrix with chr and bin, corresponding to the restriction fragments sorted by ID
    bin <- cbind(chr = rep.mat[match(RFanno[,1], rep.mat$chr),2], bin = 1 + RFanno[ ,2] %/% binSize)

    # bin_ref are the bins that are numbered continuously along all chromosomes
    # this ensures that the bin count doesn't start with 1 every time a new chromosome starts
    bin_ref <- NULL
    id <- 0
    for(i in seq_along(rep.mat[,2])){
      bin_ref <- c(bin_ref, id + bin[bin[, 1] == i, 2])
      id <- max(bin_ref)
    }

    # Reference the nth restriction fragment from the bin_ref vector: for every RF_id1 or RF_id2 entry (from table RFpairs),
    # writes out the referenced bins.
    # binPairs is a dataframe with two columns containing all the bins that interact with each other
    binPairs <- cbind(bin_ref[as.numeric(RFpairs[ ,1])], bin_ref[as.numeric(RFpairs[ ,2])])

    pairs <- binPairs
    dimension <- max(bin_ref)
  }else{
    pairs <- RFpairs
    dimension <- nrow(RFanno)
  }

  # Express the pairs of RF as unique identifiers
  # This may take a while (around minute for ~ 100M pairs)

  pairID <- RFpair_to_pairID(pairs, RFanno)

  # count how many times a given pair was present in our RFpairs.
  pairCount <- plyr::count(pairID)     #column "x": pairID, column "freq": count

  # convert IDs back to RF pairs
  frags <- pairID_to_RFpair(as.numeric(pairCount$x), RFanno)

  # We create a matrix m containing RF ids for both mates and the frequency of the corresponding ligation products
  # For a library with 100M pairs this step takes 5 minutes
  m <- cbind(frags, pairCount$freq) #

  # Now create a genome wide LFM for these connections.
  LFM <- Matrix::Matrix(0, nrow = dimension, ncol = dimension, sparse=TRUE)
  LFM[ m[, 1:2] ] <- m[,3]

  # Genomic annotation of the ligation frequency matrix.
  if(is.null(binSize)){
    # Here the rownames will hold the genomic annotation of each restriction fragment.
    rownames(LFM) <- paste(RFanno[,1], RFanno[,2], RFanno[,3], sep = "_")
    colnames(LFM) <- rownames(LFM)
  }else{
    matrixAnno <- do.call("c", lapply(unique(RFanno[, 1]), function(f){
            paste(f, seq(1, 1 + max(RFanno[RFanno[,1] == f,2]), by = binSize),
            seq(1, 1 + max(RFanno[RFanno[,1] == f,2]), by = binSize) + binSize - 1, sep = "_") }))
    rownames(LFM) <- colnames(LFM) <- matrixAnno
  }

  return(LFM)
}



#' Encode a restriction fragment pair or bin pair into one number, the pair ID.
#'
#' Using the restriction fragment pairs or bin pairs and their genomic annotation, a unique pair ID is generated.
#'
#' @param pairs A dataframe with two columns containing restriction fragment pairs or bin pairs. The mate with the lower ID is
#' always in the left column, and all pairs are sorted in ascending order first by the left column and then by the right colunn.
#' @param RFanno A dataframe giving the genomic annotation of the restriction fragments or bins. Columns are: chromosome, start of
#' restriction fragment or bin, end of restriction fragment or bin, restriction fragment or bin ID. The table is sorted by
#' chromosomal position and restriction fragments or bins are numbered along their genomic position.
#' @return A vector containing the pair IDs.
#' @examples
#' # restriction fragment pairs are given here
#' mypairs <- data.frame(c(1:15), c(16:30))
#'
#' # genomic annotation of the restriction fragment pairs
#' chr <- c(rep.int(1,10), rep.int(2,10), rep.int(3,10))
#' start <- rep.int(seq(1, 30, by = 3), 3)
#' end <- start + 2
#' RF_id <- c(1:30)
#' myanno <- data.frame(chr, start, end, RF_id)
#'
#' # generate the encoded pairs
#' code <- RFpair_to_pairID(mypairs, myanno)
#' # 16  47  78 109 140 171 202 233 264 295 326 357 388 419 450
#'
RFpair_to_pairID <- function(pairs, RFanno){
  a <- pairs[,1] - 1
  b <- pairs[,2] - 1
  base <- as.numeric(nrow(RFanno)) # number of restriction fragments
  return( a * base + b + 1)  # returns a vector of encoded pairs
}



#' Decode the pair ID into the bin pairs or restriction fragment pairs.
#'
#' Using the genomic annotation, the unique pair ID is decoded into the two restriction fragment IDs or bin IDs.
#'
#' @param pairID A vector with pair IDs that each encode a pair of restriction fragments or bins.
#' @param RFanno A dataframe giving the genomic annotation of the restriction fragments or bins. Columns are: chromosome, start of
#' restriction fragment or bin, end of restriction fragment or bin, restriction fragment or bin ID. The table is sorted by chromosomal
#' position and restriction fragments or bins are numbered along their genomic position.
#' @return A dataframe with two columns containing restriction fragment pairs or bin pairs.
#' @examples
#' # pair IDs are given here
#' code <- c(16,47,78,109,140,171,202,233,264,295,326,357,388,419,450)
#'
#' # genomic annotation of the restriction fragment pairs
#' chr <- c(rep.int(1,10), rep.int(2,10), rep.int(3,10))
#' start <- rep.int(seq(1, 30, by = 3),3)
#' end <- start + 2
#' RF_id <- c(1:30)
#' myanno <- data.frame(chr, start, end, RF_id)
#'
#' # generate the decoded pairs
#' mypairs <- pairID_to_RFpair(code, myanno)
#'
#' > mypairs
#'       [,1] [,2]
#' [1,]    1   16
#' [2,]    2   17
#' [3,]    3   18
#'  ..     ..  ..
#'
pairID_to_RFpair <- function(pairID, RFanno){
  stopifnot(is.null(dim(pairID)) || ncol(pairID) == 1, !any(pairID < 1))
  base <- as.numeric(nrow(RFanno))
  pairID <- pairID - 1
  return(cbind(pairID %/% base, pairID %% base) + 1)
}

