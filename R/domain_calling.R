#' Domain calling using a variant of the TopDom algorithm.
#'
#' TODO: add description of the algorithm.
#'
#' @param Iset An interactionSet object, with the bin size stored as metadata
#' under the name "binSize".
#' @param window_size Size of the sliding window around each bin (given in number
#' of bins). The bins  max(bin - window_size  -1, 1) are considered downstream and
#' the bins max(bin + window_size, total_no_of_bins) are considered upstream of the
#' current bin
#' @param pv_thresh The p-value threshold that is used for the post-hoc filtering
#' of domain boundaries (step 3 of the algorithm.).
#' @param parallel if FALSE, no parallelization. if TRUE, parallel
#' execution using \code{BiocParallel}, see next argument \code{BPPARAM}.
#' A note on running in parallel using \code{BiocParallel}: it may be
#' advantageous to remove large, unneeded objects from your current
#' R environment before calling the function,
#' as it is possible that R's internal garbage collection
#' will copy these files while running on worker nodes.
#' @param BPPARAM an optional parameter object passed internally
#' to \code{\link{bplapply}} when \code{parallel=TRUE}.
#' If not specified, the parameters last registered with
#' \code{\link{register}} will be used.
#'
#' @examples
#' 1+1
#'
#' @referencs
#' Shin et. al. - TopDom : An efficient and Deterministic Method for identifying Topological Domains in
#' Genomes, 2016, http://dx.doi.org/10.1093/nar/gkv1505
#' @export domain_calling_TopDom
#' @import BiocParallel
#' @importFrom psych logit
#' @importFrom purrr map_dbl


domain_calling_TopDom <- function(Iset, window_size = 5, pv_thresh = 0.05,
                                 parallel=FALSE, BPPARAM=bpparam()){

  if( !("binSize" %in% names(metadata(Iset))) ){
    stop("The Iset you are providing is not binned. Please bin it first.")
  }

  lfm <- Iset_to_LFM(Iset)$LFM

  # turn window size into an index set around a bin
  idx <- c(-seq_len(window_size  + 1), 0, seq_len(window_size))

  total_no_bins <- nrow(lfm)

  bins_idx <- map(seq_len(total_no_bins), .get_bin_idx,
                  idx = idx, total_no_bins = total_no_bins)

  signal <- bplapply(bins_idx, .get_signal)

  logit_signal <- map(signal, map, ~.my_logit(.x))

  mean_logit <- map(logit_signal, map_dbl, ~mean(.x))

  mean_diamond_logit <- map_dbl(mean_logit, "diamond")

  # use the TopDom built in function to identify changepoints
  cps <- .Change.Point(seq_along(mean_diamond_logit),
                       mean_diamond_logit)$cp

  # identify local minima from the fitted curve
  pits <- .find_pits(cps, total_no_bins)

  # posthoc filter step of the minima based on testing
  pv_bins <- map_dbl(logit_signal, .get_posthoc_pvalue)
  pits <- subset(pits, pv_bins[pits] < pv_thresh)


  # assemble the return values

  # add the TAD IDs to the metadata columns of the regions
  mcols(regions(Iset))$TD_IDs <- rep(seq_along(c(pits, 1)), c(pits[1], diff(pits),
                                    total_no_bins - tail(pits, 1)))

  # add the mean signal of the bin diamond and upstream and downstream triangles
  mcols(regions(Iset))$mean_diamond_logit <- mean_diamond_logit
  mcols(regions(Iset))$mean_upstream_logit <- map_dbl(mean_logit, "upstream_triangle")
  mcols(regions(Iset))$mean_upstream_logit <- map_dbl(mean_logit, "downstream_triangle")

  # add the identified pits (local minima) to the metadata
  metadata(Iset)$TD_boundaries <- pits


  return(Iset)
}

# The index returned extracts the upper right corner from the bin onwards,
# i.e. the bin position is at the lower left end of the rectangle
# so you go up and then to the right from the diagonal!

.get_bin_idx <- function(bin, idx, total_no_bins){
  tmp <- sort(unique(pmin(pmax(1, bin + idx), total_no_bins)))
  list(row_idx = tmp[tmp <= bin], col_idx = tmp[tmp >= bin])
}



# get the signal in the diamond, corresponding to the signal between the bins in the
# upstream and downstream windows, and in the triangles corresponding
# to the interactions of the bins in the upstream and downstream windows with each other
.get_signal <- function(bin_idx, dat_mat = chr_12_ic_pooled){

  diamond <- as.vector(dat_mat[bin_idx$row_idx, bin_idx$col_idx])
  upstream_triangle <- dat_mat[bin_idx$row_idx, bin_idx$row_idx]
  downstream_triangle <- dat_mat[bin_idx$col_idx, bin_idx$col_idx]

  res <- list("diamond" = diamond,
              "upstream_triangle" = upstream_triangle,
              "downstream_triangle" = downstream_triangle)

  res[c("upstream_triangle",
        "downstream_triangle")] <- map(res[c("upstream_triangle",
                                          "downstream_triangle")],
                                                            ~.x[upper.tri(.x)])

  return(res)

}


# linear fitting function copied from TopDom 0.0.2

.Change.Point <- function( x, y )
{
  if( length(x) != length(y))
  {
    print("ERROR : The length of x and y should be the same")
    return(0)
  }

  n_bins <- length(x)
  Fv <- rep(NA, n_bins)
  Ev <- rep(NA, n_bins)
  cp <- 1

  i=1
  Fv[1]=0
  while( i < n_bins )
  {
    j=i+1
    Fv[j] = sqrt( (x[j]-x[i])^2 + (y[j] - y[i] )^2 )

    while(j<n_bins)
    {
      j=j+1
      k=(i+1):(j-1)
      Ev[j] = ( sum( abs( (y[j]-y[i] )*x[k] - (x[j] -x[i])*y[k] - (x[i]*y[j]) + (x[j]*y[i]) ) ) / sqrt( (x[j]-x[i])^2 + (y[j] - y[i] )^2 ) )
      Fv[j] = sqrt( (x[j]-x[i])^2 + (y[j] - y[i])^2 ) - ( sum( abs( (y[j]-y[i] )*x[k] - (x[j] -x[i])*y[k] - (x[i]*y[j]) + (x[j]*y[i]) ) ) / sqrt( (x[j]-x[i])^2 + (y[j] - y[i] )^2 ) )

      #################################################
      #Not Original Code
      if( is.na(Fv[j]) || is.na(Fv[j-1]) ) {
        j = j-1
        cp <- c(cp, j)
        break
      }
      ####################################################3
      if(Fv[j] < Fv[j-1] ) {
        j = j - 1
        cp <- c(cp, j )
        break
      }
    }
    i=j
  }

  cp <- c(cp, n_bins)

  return(list(cp=cp, objF=Fv, errF=Ev))
}


# check slopes at neighbouring turning points to find local minima

.find_pits <- function(cps, total_no_bins, mean_signal = mean_diamond_logit){

  x <- seq_len(total_no_bins)[cps]
  y <- mean_signal[cps]

  slopes <- diff(y) / diff(x)
  slopes_binary <- slopes
  slopes_binary[slopes_binary != 0] <- ifelse(slopes[slopes != 0] >= 0, +1L, -1L)

  cps[which(diff(slopes_binary) >= 2) + 1]

}

# logit transformation with thresholding of large and small value


.my_logit <- function(x){

  map_dbl(pmin(pmax(logistic(-10), x), logistic(10)),  logit)

}


# posthoc filter step: get p-values bases on wilcoxon that test
# whether the signal in the diamond is lower than in the up-
# and downstream triangles. A low p-value indicates a domain boundary.

.get_posthoc_pvalue <- function(bin_signal){

  pv <-  wilcox.test(x = bin_signal$diamond,
                     y = c(bin_signal$upstream_triangle,
                           bin_signal$downstream_triangle),
                     alternative = "less",
                     exact = FALSE)$p.value
}
