#' Convert a sparse ligation frequency matrix into an InteractionSet object
#'
#'
#' The data can be binned or not binned, normalized or not normalized.
#'
#' @param LFM A sparse ligation frequency matrix, with the genomic annotation as row and column names.
#' @param binned Boolean, whether the matrix is already binned. Defaults to TRUE.
#' @return An interaction set object containing the information from the LFM and the genomic annotation,
#' and with the binSize stored as metadata
#' @examples
#' LFM <- Matrix::Matrix(c(1,2,1,2,5,0,1,0,0), 3, sparse = TRUE)
#' Iset <- LFM_to_Iset(LFM, binned = FALSE)
#'
#' @export LFM_to_Iset

LFM_to_Iset <- function(LFM, binned = FALSE) {

    LFM <- sym_to_triangle(LFM)
    if (sum(LFM) == 0)  warning("You are converting an empty matrix.")
    if (is.null(rownames(LFM))) stop("No row names present. Row and column names of the
                                     need to be provided in chr_start_end format, e.g. 4_201000_202000.")

    # summarize normalized (or non-normalized) sparse matrix
    summLFM <- Matrix::summary(LFM)

    # order by RFpairs
    summLFM <- summLFM[SummarizedExperiment::order(summLFM[, 1], summLFM[, 2]), ]

    # get the annotation from the LFM rownames
    strings <- strsplit(rownames(LFM), split = "_")
    anno <- do.call(rbind, strings)
    anno <- as.data.frame(anno, stringsAsFactors = FALSE)
    names(anno) <- c("chr", "start", "end")
    anno[, 2:3] <- data.frame(lapply(anno[, 2:3], as.numeric))
    anno$ID <- seq_along(anno$start)
    ranges <- GenomicRanges::GRanges(anno$chr, IRanges::IRanges(anno$start, anno$end), IDs = anno$ID)

    if (binned == TRUE) {
        binSize <- (anno$end)[1] - (anno$start)[1] + 1
    } else {
        binSize <- NULL
    }

    # create GInteractions object for output
    outGI <- InteractionSet::GInteractions(summLFM[, 1], summLFM[, 2], ranges)

    # put the normalized matrix values into the InteractionSet object
    freq <- matrix(summLFM[, 3])
    Iset <- InteractionSet::InteractionSet(freq, outGI)

    S4Vectors::metadata(Iset)$binSize <- binSize

    return(Iset)
}

#' Convert an InteractionSet object into a ligation frequency matrix.
#'
#' The data can be binned or not binned, normalized or non--normalized.
#'
#' @param Iset An InteractionSet object
#' @param assay_number If the Iset contains multiple assays, this gives the
#' number of the assay to extract.
#' @param Iset_col The column of the assay to extract. This
#' column treated as corresponding to one ligation frequency
#' matrix. The function creates one ligation frequency matrix from one column of
#' the assay slot at a time.
#' @return A list consisting of: the sparse ligation frequency matrix with t
#' he genomic annotation as rownames and colnames,
#' the restriction fragment annotation (or bin annotation), and the restriction
#' fragment pairs (or bin pairs).
#' @examples
#' LFM <- Matrix::Matrix(c(1,2,1,2,5,0,1,0,0), 3, sparse = TRUE)
#' Iset <- LFM_to_Iset(LFM, binned = FALSE)
#' LFM_back <- Iset_to_LFM(Iset)[[1]]
#' LFM == LFM_back
#'
#' @export Iset_to_LFM

Iset_to_LFM <- function(Iset, iset_col = 1, assay_number = 1) {

    # Create RFanno object
    RFanno <- data.frame(InteractionSet::regions(Iset), stringsAsFactors = FALSE)
    # add IDs column if it  does not exist
    if( !("IDs" %in% colnames(RFanno)) ){
      RFanno$IDs <- seq_len(nrow(RFanno))
    }

    RFanno <- RFanno[,  c("seqnames", "start", "end", "IDs")]
    names(RFanno) <- c("chr", "start", "end", "RF_id")
    RFanno$chr <- as.character(RFanno$chr)
    RFanno$start <- as.numeric(RFanno$start)
    RFanno$end <- as.numeric(RFanno$end)
    RFanno <- RFanno[order(RFanno$RF_id), ]

    # Create RFpairs object
    Iset <- InteractionSet::swapAnchors(Iset)
    RFpairs <- as.data.frame(InteractionSet::anchors(Iset, id = TRUE))
    RFpairs <- RFpairs[order(RFpairs[, 1], RFpairs[, 2]), ]
    names(RFpairs) <- c("RF1", "RF2")

    # Store the interaction pairs and frequencies in a summary(LFM)
    summLFM <- as.data.frame(InteractionSet::anchors(Iset, id = TRUE))
    summLFM <- cbind(as.matrix(summLFM),
               SummarizedExperiment::assay(Iset, assay_number)[, iset_col])

    # Get matrix size
    matSize <- length(InteractionSet::regions(Iset))

    LFM <- Matrix::Matrix(0, nrow = matSize, ncol = matSize, sparse = TRUE)
    LFM[summLFM[, 1:2]] <- summLFM[, 3]

    rownames(LFM) <- paste(RFanno$chr, RFanno$start, RFanno$end, sep = "_")
    colnames(LFM) <- rownames(LFM)

    return(list(LFM = LFM, RFanno = RFanno, RFpairs = RFpairs))
}

#' Convert restriction fragment annotation and restriction fragment pair objects to an InteractionSet object.
#'
#' @param RFanno A dataframe giving the genomic annotation of the restriction fragments. Columns are: chromosome, start of
#' restriction fragment, end of restriction fragment, restriction fragment ID. The table is sorted by chromosomal position
#' and restriction fragments are numbered along their genomic position.
#' @param RFpairs A dataframe with two columns containing restriction fragment pairs. The restriction fragment with the
#' lower ID is always in the left column, and all restriction fragments are sorted in ascending order first by the left
#' column and then by the right colunn.
#' @param binned Boolean, whether the matrix is already binned. Defaults to FALSE.
#' @return An InteractionSet object.
#' @examples
#' # Iset <- RFpairs_to_Iset(RFanno, RFpairs)
#'

RFpairs_to_Iset <- function(RFanno, RFpairs, binned = FALSE) {

    ranges <- GenomicRanges::GRanges(RFanno$chr, IRanges::IRanges(RFanno$start, RFanno$end), IDs = RFanno$RF_id)

    pairID <- RFpair_to_pairID(RFpairs, RFanno)

    # count how many times a given pair was present in our RFpairs.
    pairCount <- plyr::count(pairID)  #column 'x': pairID, column 'freq': count

    # convert IDs back to RF pairs
    frags <- pairID_to_RFpair(as.numeric(pairCount$x), RFanno)

    # We create a matrix m containing RF ids for both mates and the frequency of the corresponding ligation products
    summLFM <- cbind(frags, pairCount$freq)  #

    # create GInteractions object for output
    outGI <- InteractionSet::GInteractions(summLFM[, 1], summLFM[, 2], ranges)

    # put the normalized matrix values into the InteractionSet object
    freq <- matrix(summLFM[, 3])
    Iset <- InteractionSet::InteractionSet(freq, outGI)

    if (binned) {
        S4Vectors::metadata(Iset)$binSize <- (RFanno$end)[1] - (RFanno$start)[1] + 1
    } else {
      S4Vectors::metadata(Iset)$binSize <- NULL
    }

    return(Iset)
}

#' Convert a sparse symmetric matrix to a sparse triangle matrix.
#'
#' @param mat A symmetric sparse matrix.
#' @return A triangular sparse matrix, containing only the upper triangle.
#' @examples
#' symmMatrix <- Matrix::Matrix(c(1,2,1,2,5,0,1,0,0), 3, sparse = TRUE)
#' triangleMatrix <- sym_to_triangle(symmMatrix)
#'
#' @export sym_to_triangle
#'
sym_to_triangle <- function(mat) {
    summ <- Matrix::summary(mat)
    summ <- summ[summ$i <= summ$j, ]
    trimat <- Matrix::sparseMatrix(i = summ$i, j = summ$j, x = summ$x, dims = dim(mat))
    rownames(trimat) <- rownames(mat)
    colnames(trimat) <- colnames(mat)
    return(trimat)
}

#' Convert a sparse triangle matrix to a sparse symmetric matrix.
#'
#' @param trimat A triangular sparse matrix (either upper or lower triangle can be filled).
#' @return A symmetric sparse matrix.
#' @examples
#' triangleMatrix <- Matrix::Matrix(c(1,0,0,2,5,0,1,0,0), 3, sparse = TRUE)
#' symmMatrix <- triangle_to_sym(triangleMatrix)
#'
#' @export triangle_to_sym
#'
triangle_to_sym <- function(trimat) {
    summ <- Matrix::summary(trimat)
    summ <- summ[summ$i <= summ$j, ]
    mat <- Matrix::sparseMatrix(i = summ$i, j = summ$j, x = summ$x, dims = dim(trimat), symmetric = TRUE)
    rownames(mat) <- rownames(trimat)
    colnames(mat) <- colnames(trimat)
    return(mat)
}

#' Convert an InteractionSet region into a ligation frequency matrix.
#'
#' Pick a genomic region from an InteractionSet object and convert it into a sparse triangular matrix.
#'
#' @param Iset An interactionSet object, with the bin size stored as metadata under the name 'binSize'.
#' @param chr The chromosome name to be plotted with or without 'chr', depending on your interactionSet object.
#' @param from Start of the region you want to select. If NULL, the start of the chromosome will be 'from'.
#' @param to End of the region you want to select. If NULL, the end of the chromosome will be 'to'.

#' @examples
#' # matrix <- Iset_region_to_matrix(Iset, chr = 11, from = 30000000, to = 30100000)
#'

Iset_region_to_LFM <- function(Iset, chr, from = NULL, to = NULL) {

    # maximal entry on the selected chromsome
    max <- max(SummarizedExperiment::end(InteractionSet::regions(Iset)[GenomicRanges::seqnames(InteractionSet::regions(Iset)) == chr]))

    if (is.null(from))
        from <- 1
    if (is.null(to))
        to <- max

    # input checking
    if (from < 0 || to > max)
        stop("Subset out of bounds. Provide a region within the range of the chromosome")
    if (from > to)
        stop("'from' has to be smaller than 'to'")

    ## convert interactionSet to matrix
    Iset <- InteractionSet::swapAnchors(Iset)
    region <- GenomicRanges::GRanges(chr, IRanges::IRanges(from, to))
    mat <- InteractionSet::inflate(Iset, region, region, sparse = TRUE,
                                   swap = FALSE, fill = SummarizedExperiment::assay(Iset))
    mat <- Matrix::as.matrix(mat)

    ## old version indices=which(start(ranges(regions(Iset))) >= from & start(ranges(regions(Iset))) <= to & as.vector(seqnames(regions(Iset))) == chr) mat =
    ## InteractionSet::inflate(Iset, indices, indices, sparse=TRUE, swap=FALSE, fill=SummarizedExperiment::assay(Iset))

    ## produce rownames and colnames for the matrix
    olap = IRanges::subsetByOverlaps(InteractionSet::regions(Iset), region)
    names <- as.data.frame(olap)
    rownames(mat) <- colnames(mat) <- paste(names$seqnames, names$start, names$end, sep = "_")

    return(mat)
}

#' Convert an InteractionSet into a sparse, upper triangle ligation frequency matrix.
#'
#' Wrapper for Iset_to_LFM (no subsetting) and Iset_region_to_LFM (subsetting by chromosome and, if required, by region).
#'
#' @param Iset An interactionSet object, with the bin size stored as metadata under the name 'binSize'.
#' @param chr Optional. The chromosome name to be plotted with or without 'chr', depending on your interactionSet object.
#' @param from Optional. Start of the region you want to select. If NULL, the start of the chromosome will be 'from'.
#' @param to Optional. End of the region you want to select. If NULL, the end of the chromosome will be 'to'.
#' @return A sparse, upper triangle ligation frequency matrix.
#' @examples
#' # LFM <- easy_Iset_to_LFM(Iset)
#' # chr2_1_300000 <- easy_Iset_to_LFM(Iset, chr = 2, from = 1, to = 300000)
#'
easy_Iset_to_LFM <- function(Iset, chr = NULL, from = NULL, to = NULL) {
    if (is.null(chr)) {
        if (!is.null(from) || !is.null(to))
            stop("'from' and 'to' cannot be defined unless 'chr' is defined.")
        LFM <- Iset_to_LFM(Iset)[[1]]
    } else LFM <- Iset_region_to_LFM(Iset, chr, from, to)
    return(LFM)
}
