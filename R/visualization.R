#' Triangle visualization of Hi-C data with tracks for chipseq, genome axis, gene models, and custom data.
#'
#' An interaction matrix is plotted, rotated, and if required, smoothed. Using Gviz, tracks can be added as required.
#'
#' @param Iset An interactionSet object, with the bin size stored as metadata under the name "binSize"
#' @param chr The chromosome name to be plotted with or without "chr", depending on your interactionSet object.
#' @param from Start of the region you want to plot.
#' @param to End of the region you want to plot.
#' @param gen The genome identifier, used primarily to determine which gene model will be used.
#' Give the organism name (e.g. "hsapiens", "mmusculus") if you want to use the most recent assembly from biomart (internet
#' connection required). Give one of "mm9", "mm10", "hg19", or "hg38" to use a gene model bundled with the package.
#' Leave out if you want to use your own gene model.
#' @param geneModels Only input a geneModel GRanges if you want to use your own custom geneModel.
#' @param printTranscripts Boolean. Produce gene model track? Defaults to TRUE.
#' @param outDir Give the directory where the PDF plot will be saved. Defaults to the working directory.
#' @param name A string. Filename of the output PDF.
#' @param scaleCol The color scale threshold. All values above the threshold will ne shown in the highest intensity color.
#' @param chipseqs A named list of GRanges objects with ChIP-seq scores in the metadata column "score". Optional.
#' @param chipScale A named list corresponding to chipseqs with a vector containing the highest and lowest score for each list
#' element. Defaults to taking the non-zero 99\% quantile (rounded to the nearest ten) as the maximum and -0.5 as the minimum value.
#' @param colmapChipseqs A named character vector corresponding to chipseqs with a color name for each chipseq. Defaults to a
#' vector of twelve colors. For a larger number of chipseq tracks, please define colmapChipseqs.
#' @param customAnno A custom annotation, given as a list of BED dataframes. Optional.
#' @param smooth Which smoothing type should be applied? Options are: "none", "box", and "gaussian".
#' @param smoothing Smoothing parameter. The percentage of signal to be taken from the surrounding pixels. Smoothing is optional.
#' @param filterSize Smoothing parameter. The dimension of filter mask, e.g. filterSize = 3 means 3x3 filter. Smoothing is optional.
#' @param sigma In the case of Gaussian smoothing, how large should the sigma of the distribution be?
#' @return Output will be in PDF format.
#' @examples
#' visualize(Iset, chr = "11", from = 30000000, to = 30100000,
#'           gen = "mm9", customAnno = cpgIslands, smoothing = 80, filterSize = 3)
#'
#' @export visualize
#'

visualize <- function(Iset, chr, from = NULL, to = NULL, gen = NULL, geneModels = NULL, printTranscripts = TRUE, outDir = getwd(),
                                  name = "triangle_visualization.pdf", scaleCol = 0.02, chipseqs = NULL, chipScale = NULL,
                                  colmapChipseqs = NULL, customAnno = NULL, smooth = "none", smoothing = NULL, filterSize = NULL,
                                  sigma = NULL){

  stopifnot(xor(is.null(gen), is.null(geneModels)))
  binSize <- S4Vectors::metadata(Iset)$binSize

  max <- max(end(InteractionSet::regions(Iset)[GenomicRanges::seqnames(InteractionSet::regions(Iset)) == chr]))

  if (is.null(from))
    from <- 1
  if (is.null(to))
    to <- max

  mat <- Iset_region_to_LFM(Iset, chr, from, to)
  mat <- triangle_to_sym(mat)
  img <- as(mat, "matrix")


  ### Build the tracks ###

  ## Rotated matrix track #vary size =37
  interaction.matrix <- Gviz::CustomTrack(plottingFunction = function(GdObject, prepare){
      if(!prepare) rotated_image(img, scaleCol, smooth, smoothing, filterSize, sigma, binSize)
      return(GdObject)
    }, size = 37, name = "rotated_image", col.title = "white", background.title = NA, showTitle = F)


  ## GeneModel track

  if(printTranscripts){

    if(S4Vectors::grepl("chr", chr) == 0) chr.reg = paste0("chr", chr) else chr.reg = chr
    id <- c("mm9", "mm10", "hg19", "hg38")

    if(!is.null(gen) && gen %in% id){

      # offline: load pregenerated generegiontracks
      genomes <- data.frame(id = id, file = c("gm_mm9.RData", "gm_mm10.RData", "gm_hg19.RData", "gm_hg38.RData"))
      rdata <- as.character(genomes[genomes$id == gen,]$file)
      load(file.path("data", rdata))
    }
    if(!is.null(gen) && gen %in% id && !is.null(geneModels)){

      # offline: use custom geneModels provided by user or already in the environment
      gm <- Gviz::GeneRegionTrack(geneModels, chromosome = chr, genome = 'gen', #start = from, end = to,
                                  #transcriptAnnotation = "symbol", stacking = "squish",
                                  size = 1, name = "RefSeq", col.line = NULL, cex.group = 0.7,
                                  fill = "black", utr5 = "black", utr3 = "black", protein_coding = "black", #make all annotations black
                                  collapseTranscripts = "longest", background.title = NA#, showTitle = F
      )
    }else if(is.null(geneModels)){

      # online: get most recent gene model for any organism from biomart
      dataset <- paste0(gen, "_gene_ensembl")
      mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = dataset)
      if(!dataset %in% biomaRt::listDatasets(mart)$dataset) stop(gen, " is not a valid organism identifier.
                                                      Correct examples: mmusculus, hsapiens")
      gm <- Gviz::BiomartGeneRegionTrack(chromosome = chr.reg, genome = gen , biomart = mart, start = from, end = to,
                                         size = 1, name = "RefSeq", col.line = NULL, cex.group = 0.7,
                                         fill = "black", utr5 = "black", utr3 = "black", protein_coding = "black", #make all annotations black
                                         collapseTranscripts = "longest",  background.title = NA#, showTitle = FALSE
      )
    }
  }else{
    gm <- NULL
  }

  ## Chipseq track

  if(length(chipseqs) > 0) {
    if(is.null(colmapChipseqs) && length(chipseqs) <= 12){
      colmapChipseqs <- c("steelblue", "gray", "tomato3", "green4", "black", "peachpuff3", "steelblue", "gray", "tomato3",
                        "green4", "black", "peachpuff3")
      names(colmapChipseqs) <- names(chipseqs)
    }

    ct <- lapply(names(chipseqs), function(csl) {
      #chipseqs[[csl]] <- sort(chipseqs[[csl]])
      theRegion <- chipseqs[[csl]][chipseqs[[csl]]@seqnames == chr.reg]

      # set the minimum and maximum values of the chipseq scale
      if(is.null(chipScale)){
        # default scale is from -0.5 to the 99% quantile value of the non-zero scores
        ylim <- c(-0.5, round( stats::quantile( score(chipseqs[[csl]][score(chipseqs[[csl]]) != 0,]),
                                         probs = seq(0,1,0.01))[100], digits = -1))
      }else{
        # user-specified scale
        ylim <- c(as.numeric(chipScale[[csl]][[1]]), as.numeric(chipScale[[csl]][[2]]))
      }

      Gviz::DataTrack(theRegion, strand = "*", start = from, end = to,
                genome = gen, col.histogram = colmapChipseqs[[csl]],
                fill.histogram = colmapChipseqs[[csl]], name = csl,
                col.axis = colmapChipseqs[[csl]],  size = 4,# cex.legend = 2, #cex.axis = 0.75,
                ylim = ylim, col.title = colmapChipseqs[[csl]], separator = 2,
                background.title = "white") } )
  }else{
    ct <- NULL
  }


  ## custom annotation track

  if(length(customAnno) > 0) {
    ca <- lapply(names(customAnno), function(csl) {
      theRegion <- customAnno[[csl]]
      theRegion <- theRegion[ theRegion[,1] == chr.reg & theRegion[,3] >= from & theRegion[,2] <= to, ]
      if(nrow(theRegion) > 0){
        if(ncol(theRegion) >= 5) score <- theRegion[,5] else score <- NULL
        if(ncol(theRegion) >= 6 && theRegion[,6] %in% c('+', '-', 'NA')){
          strand <- theRegion[,6]
        } else {
            strand <- NULL
            }
        theRegion <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(chr.reg),
                             ranges = IRanges::IRanges(start = as.numeric(theRegion[,2]), end = as.numeric(theRegion[,3])),
                             score = score,
                             strand = strand
        )

      }else{
        theRegion <- NULL
      }
      Gviz::AnnotationTrack(theRegion, genome = gen, name = "", col.axis = "black", cex.axis = 1,
                      background.title = "NA", chromosome = chr.reg, shape = 'arrow', fill = 'gray50', size = 2)
    } )

  }else{
    ca <- NULL
  }

  ax <- Gviz::GenomeAxisTrack(littleTicks = FALSE, add35 = TRUE, add53 = TRUE, distFromAxis = 5, cex = 1.2)


  #plot pdf
  grDevices::pdf(file.path(outDir, name), width = 30, height = 35)
  Gviz::plotTracks(c(interaction.matrix, gm, ct, ca, ax),
             from = from, to = to, margin = 100, col = NULL, col.title="grey", showTitle = F, cex.title = 0.5,
             transcriptAnnotation = "symbol", window = "auto", type = "histogram", fontsize = 35)
  grDevices::dev.off()
}


#' Plot a rotated and, if required, smoothed interaction matrix.
#'
#' A matrix with the interactions is rotated and cut along the diagonal to form a triangle plot of the interactome. Smoothing
#' is performed utilizing the filter2 convolution filter from the EBImage package.
#'
#' @param img The symmetric, non-sparse interaction matrix (for feasibility, choose only a limited chromosomal region).
#' @param scaleCol The color scale threshold. All values above the threshold will ne shown in the highest intensity color.
#' @param smoothing Smoothing parameter. The percentage of signal to be taken from the surrounding pixels. Smoothing is optional.
#' @param smooth Which smoothing type should be applied? Options are: "none", "box", and "gaussian".
#' @param filterSize Smoothing parameter. The number of pixels surrounding the central pixel. (Size of filter mask.) Smoothing is optional.
#' @param sigma In the case of Gaussian smoothing, how large should the sigma of the distribution be?
#' @param binSize Size of bins into which the interaction matrix is binned.
#' @examples
#' myTriangle <- rotatedImageIC(img, scaleCol = 0.02, smoothing = 80, filterSize = 3, binSize = 10000)
#'
# ' @import Gviz
# ' @import EBImage
# ' @import grid
#'

rotated_image <- function(img, scaleCol, smooth, smoothing, filterSize, sigma, binSize){
  if(!smooth %in% c("none", "box", "gaussian")) stop("Invalid smooth setting. Please set to 'none', 'box', or 'gaussian'.")
  start <- as.integer(sapply(strsplit(colnames(img), "_"), function(x) x[2]))
  end <- as.integer(sapply(strsplit(colnames(img), "_"), function(x) x[3]))
  x1 <- rep(start, each = length(start))
  y1 <- rep(start, length(start))

  ncols <- 256
  lims <- c(base::min(x1) - 1, base::max(x1) + binSize - 1)
  gradient <- c(grDevices::colorRampPalette(c("white", "gray90", "blue2", "red4", "orange", "yellow"))(ncols - 1), "#A8A8A8")
  colorScale <- seq(0, scaleCol, length.out = ncols)

  ## smoothing
  if(smooth == "box"){
    if(is.null(smoothing)) smoothing <- 80
    if(is.null(filterSize)) filterSize <- 3
    stopifnot(smoothing <= 100 && smoothing > 0 && filterSize >= 3 && filterSize %% 2 == 1)
    surround <- filterSize^2 - 1

    filter <- matrix(rep(smoothing / (surround * 100), filterSize^2), nrow = filterSize)
    middle <- stats::median(filterSize)
    filter[middle, middle] <- (100 - smoothing) / 100
    img <- EBImage::filter2(img, filter)

  }else if(smooth == "gaussian"){
    if(is.null(sigma)) sigma <- 0.3
    img <- EBImage::gblur(img, sigma)
  }

  ## color.scale is set here
  img[ img > scaleCol ] <- scaleCol
  val <- matrix(cut(as.numeric(img), colorScale, labels = FALSE), nrow = nrow(img), ncol = ncol(img))
  val[img == 0] <- 1

  col <- gradient[val]

  grid::pushViewport(grid::viewport(x = lims[1], width = diff(lims), default.units = "native", clip = TRUE, just = c(0, 0.5)))
  wpwdth <- Gviz:::vpLocation()$isize["width"]/2
  cl <- sqrt(sum(rep(wpwdth^2, 2)))
  grid::pushViewport(grid::viewport(width = grid::unit(cl, "inches"), height = grid::unit(cl, "inches"), y = 0, angle = -45,
                                    xscale = lims, yscale = lims))
  grid::grid.rect(x = x1, y = y1, width = binSize - 1, height = binSize, hjust = 0, vjust = 0,
                  default.units = "native", gp = grid::gpar(col = col, fill = col))
  grid::popViewport(2)
}

