#' Normalize an InteractionSet with iterative proportional fitting.
#'
#' A binned or unbinned InteractionSet is normalized using iterative proportional fitting (IPF), an implementation of the
#' iterative matrix balancing algorithm by Pukelsheim and Simeone (2009).
#'
#' @param Iset An InteractionSet with binned or unbinned data. If the data is binned, metadata(Iset)$binSize must contain
#' the bin size.
#' @param numberOfIterations The number of iterations for normalization. Usually, 20 should be more than enough.
#' @param plotLoss Optional boolean parameter, to determine if the loss function should be plotted. Defaults to FALSE.
#' @return A list containing the normalized InteractionSet, i.e. assay(normIset) will contain the normalized signal. In the
#' metadata, the number of interchromosomal contacts and intrachromosomal contacts is added, as well as the contact distances,
#' i.e. the distances between the genomic positions where two interacting bins or restriction fragments start. In addition,
#' the loss function is given for every iteration of the normalization. If plotLoss = T, the metadata will also contain a
#' plotted version of the loss function, which can be displayed with metadata(normIset)$lossPlot.
#' @examples
#' norm_Iset <- normalize_chrom(raw_Iset, numberOfIterations = 40L, plotLoss = FALSE)
#'
#' @export normalize_chrom

normalize_chrom <- function(Iset, numberOfIterations = 40L, plotLoss = FALSE){

  # input checking
  if (class(Iset)!= "InteractionSet") stop("Input parameter 'Iset' needs to be an InteractionSet object.")
  if (class(numberOfIterations)!= "integer") stop("Input parameter 'numberOfIterations' needs to be a positive integer number.")
  binSize <- S4Vectors::metadata(Iset)$binSize

  # convert Iset to LFM
  list <- Iset_to_LFM(Iset)
  LFM <- list[[1]]
  RFanno <- list[[2]]
  RFpairs <- list[[3]]

  # normalize
  normList <- normalize_LFM_iteratively(LFM, numberOfIterations, plotLoss)

  # convert normalized LFM to Iset
  if(is.null(binSize)) binned=F else binned=T
  normIset <- LFM_to_Iset(normList[[1]], binned=binned)

  # store every normalization output that does not fit in the standard InteractionSet into the metadata.
  S4Vectors::metadata(normIset) <- quality_assessment(RFpairs, RFanno)
  S4Vectors::metadata(normIset)$binSize <- binSize
  S4Vectors::metadata(normIset)$lossFunction <- normList$lossFunction
  if(plotLoss) S4Vectors::metadata(normIset)$lossPlot <- normList$lossPlot

  return(normIset)
}


#' Normalize a ligation frequency matrix with iterative proportional fitting
#'
#' A binned or unbinned ligation frequency matrix (LFM) is normalized using iterative proportional fitting (IPF).
#'
#' @param LFM A sparse ligation frequency matrix (LFM) with binned or unbinned data.
#' @param numberOfIterations The number of iterations for normalization. Usually, 20 should be more than enough.
#' @param plotLoss Optional boolean parameter, to determine if the loss function should be plotted. Defaults to FALSE.
#' @return A list containing the normalized, triangular matrix \code{finalLFM}, the dataframes with the row sums
#' \code{IPFrowsums} and the column sums \code{IPFcolsums} of every iteration, and a vector with the loss function
#' for every iteration \code{lossFunction}.
#' @examples
#' M <- matrix(c(3,0,0,0,7,0,1,9,7),3,3)
#' M <- Matrix(M, sparse = TRUE)
#' normalizedLFM <- normalize_LFM_iteratively(M, 10)
#'
#' @export normalize_LFM_iteratively
#' @import ggplot2


normalize_LFM_iteratively <- function(LFM, numberOfIterations, plotLoss = FALSE){
  if(class(LFM) != "dgCMatrix") stop("Please provide a dgcMatrix as input!")
  # we need a symmetric matrix. forceSymmetric doesn't work for large matrices.
  symLFM <- triangle_to_sym(LFM)
  ipf <- IPF_alg(symLFM, numberOfIterations)

  # rowBiases and colBiases are the matrices in which the row/col sums of each iteration are stored
  # add up the biases from all iterations to form the bias vector
  rowBias <- exp(Matrix::rowSums(log(ipf$rowBiases)))
  colBias <- exp(Matrix::rowSums(log(ipf$colBiases)))

  # smoothed huber estimator (robust estimator of the mean). $s = standard deviation, $mu = location estimator
  estimator <- smoothmest::smhuber(rowBias / colBias)$mu

  # entries from zero row sums are replaced with the estimator
  rowBias[Matrix::rowSums(symLFM) == 0] <- estimator

  # scale the x so that row and colsums are 1
  rowBias <- rowBias / sqrt(estimator)

  finalLFM <- symLFM / rowBias

  # same normalization for rows and columns since the matrix is symmetric
  finalLFM <- Matrix::t(Matrix::t(finalLFM) / rowBias)

  # Convert the symmetric matrix back to a triangle matrix
  finalLFM <- sym_to_triangle(finalLFM)

  lossFunction <- ipf$lf
  lossPlot <- NULL
  if (plotLoss){
    iterations <- c(1L:length(lossFunction))
    df <- data.frame(iterations, lossFunction)
    lossPlot <- ggplot(df, aes(iterations, lossFunction)) +
               geom_bar(stat = "identity") + xlab("Iterations") + ylab("Loss function")
  }

  return( list( finalLFM = finalLFM, lossFunction = lossFunction, lossPlot = lossPlot, rowBias = rowBias, colBias = colBias))
}



#' Iterative proportional fitting algorithm
#'
#' Iterative proportional fitting (IPF) is used to normalize a sparse matrix.
#'
#' @param symLFM A symmetric sparse ligation frequency matrix (LFM) with binned or unbinned data.
#' @param numberOfIterations The number of iterations for normalization. Usually, 10 should be more than enough.
#' @return A list containing the normalized, symmetric matrix \code{normLFM}, the dataframes with the row sums
#' \code{rowBiases} and the column sums \code{colBiases} of every iteration, and a vector with the loss function \code{lf}
#' for every iteration \code{lossFunction}.
#' @examples
#' M <- matrix(c(3,0,1,0,7,9,1,9,7),3,3)
#' M <- Matrix(M, sparse = TRUE)
#' ipf <- IPF_alg(M,10)
#'

IPF_alg <- function(symLFM, numberOfIterations){

  # initialize the matrices in which the row/col sums of each iteration are stored (rowBiases and colBiases)
  rowBiases <- matrix(1, nrow = nrow(symLFM), ncol = numberOfIterations)
  colBiases <- rowBiases

  # initialize the loss function vector
  lf <- vector("numeric", length = numberOfIterations)
  normLFM <- symLFM

  for(i in seq(numberOfIterations)){

    rowsums <- Matrix::rowSums(normLFM)
    # rowSums + 1:current colsums to make it triangular
    rowBiases[ ,i] <- ifelse(rowsums > 0, rowsums, 1)
    normLFM <- normLFM / rowBiases[ ,i]   # divide the normLFM by the rowsums (or by 1)

    colsums <- Matrix::colSums(normLFM)
    colBiases[ ,i] <- ifelse(colsums > 0, colsums, 1)
    normLFM <- Matrix::t(Matrix::t(normLFM) / colBiases[ ,i])

    # append new loss function to old loss function
    lf[i] <- loss_function(rowsums, colsums)
    message("loss function: ", lf[i], "   iteration: ", i, '/', numberOfIterations)
  }

  ipf <- list(normLFM = normLFM, rowBiases = rowBiases, colBiases = colBiases, lf = lf)  # produce named list elements
  return(ipf)
}


#' Loss function
#'
#' Loss function of the iterative proportional fitting algorithm
#'
#' @param rowsums The current row sums of the matrix
#' @param colsums The current column sums of the matrix
#' @return A decimal number giving the loss function of the current row and column sums.
#' @examples
#' M <- matrix(c(3,0,1,0,7,9,1,9,7),3,3)
#' M <- Matrix::Matrix(M, sparse = TRUE)
#' lf <- loss_function(Matrix::rowSums(M), Matrix::colSums(M))

loss_function <- function(rowsums, colsums){
  rowsums <- rowsums[rowsums > 0]
  rs <- abs(rowsums - 1)              #abs: absolute value (only positive)

  colsums <- colsums[colsums > 0]
  cs <- abs(colsums - 1)

  return(0.5 * sum(c(rs,cs)))
}



