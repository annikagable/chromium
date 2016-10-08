#' Quality assessment of whole genome Hi-C data.
#'
#' Calculate number of ligation products (treated as interactions), number of inter- and intrachromosomal interactions,
#' and distances between interactions.
#'
#' here go the details
#'
#' @param RFpairs A dataframe with two columns containing restriction fragment pairs. The restriction fragment with the lower ID
#' is always in the left column, and all restriction fragments are sorted in ascending order first by the left column and then by
#' the right colunn.
#' @param RFanno A dataframe giving the genomic annotation of the restriction fragments. Columns are: chromosome, start of
#' restriction fragment, end of restriction fragment, restriction fragment ID. The table is sorted by chromosomal position
#' and restriction fragments are numbered along their genomic position.
#' @param Iset An InteractionSet object with lossFunction from normalization in its metadata.
#' @return A list with the quality assessment values interactionCount, interchromosomal, intrachromosomal, contactDistances,
#' and if it exist, the lossFunction values of each iteration of the normalization.
#' @examples example

quality_assessment = function(RFpairs, RFanno, Iset=NULL){

  # get number of ligation products
  interactionCount = nrow(RFpairs)

  # Get the number of interchromosomal and the number of intrachromosomal contacts to know if the experiment is ok
  interIntra = table(RFanno$chr[RFpairs[,1]]==RFanno$chr[RFpairs[,2]])
  inter = as.data.frame(interIntra)[1,2]
  intra = as.data.frame(interIntra)[2,2]

  # Get the distances between contacts to compare how many short or long range contacts there are.
  contactDistances = RFanno$start[RFpairs[,2]] - RFanno$start[RFpairs[,1]]

  qa = list(interactionCount=interactionCount,
            interchromosomal=inter,
            intrachromosomal=intra,
            contactDistances=contactDistances)

  if (!is.null(Iset)){
    qa$lossFunction=metadata(Iset)$lossFunction
  }
  return(qa)
}




#' Output a quality assessment report.
#'
#' Outputs a PDF with number of ligation products (treated as interactions), number of inter- and intrachromosomal interactions,
#' and plots of distances between interactions and the loss function of normalization (if data are normalized).
#'
#' @param qa The list from the quality_assessment() function.
#' @param outDir The directory where the quality report should be created.
#' @param reportName File name of the output PDF file
#'
#' @examples
#' plot_qa(qa, "./Outputs", "My_qality_report.pdf")


plot_qa = function(qa, outDir=getwd(), reportName="Quality_report.pdf"){

  # create an output directory if it doesn't exist
  if(!dir.exists(file.path(outDir))){dir.create(file.path(outDir))}

  # "plot" the quality assessment file
  pdf(file.path(outDir, reportName))
  plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")

  text(5, 9, "Quality assessment report")
  text(4, 7, "Interactioncount = ")
  text(10, 7, qa$interactionCount)
  text(4, 6, "Number of interchromosomal interactions = ")
  text(10, 6, qa$interchromosomal)
  text(4, 5, "Number of intrachromosomal interactions = ")
  text(10, 5, qa$intrachromosomal)

  distancePlot = qplot(qa$contactDistances, geom="histogram", xlab="Distance [bp]", main="Plot of contact distances")
  print(distancePlot)

  if("lossFunction" %in% names(qa)){

    iterations = c(1L:length(qa$lossFunction))
    df=data.frame(iterations, qa$lossFunction)
    lossPlot = ggplot(df, aes(iterations, qa$lossFunction)) + geom_bar(stat = "identity") +
                xlab("Iterations") + ylab("Loss function") + ggtitle("Loss function of the normalization")
    print(lossPlot)
  }
  dev.off()
}
