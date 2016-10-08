#http://www.ensembl.org/info/data/ftp/index.html

dataDirectory="./Data"
gtfFile=
#  'Mus_musculus.NCBIM37.67.gtf' #mm9
# 'Mus_musculus.GRCm38.83.gtf' #mm10
# 'Homo_sapiens.GRCh37.75.gtf' #hg19
 'Homo_sapiens.GRCh38.83.gtf' #hg38
id=
#   "mm9"
#   "mm10"
#   "hg19"
   "hg38"

gtf = import(file.path(dataDirectory, gtfFile))
seqlevelsStyle(gtf)='ucsc'
geneModels = data.frame(chromosome=as.character(chrom(gtf)),
                        start=as.numeric(start(gtf)),
                        end=as.numeric(end(gtf)),
                        width=width(gtf),
                        strand=strand(gtf),
                        feature=as.character(elementMetadata(gtf[,'source'])@listData$source),
                        gene=as.character(elementMetadata(gtf[,'gene_id'])@listData$gene_id),
                        exon=as.character(elementMetadata(gtf[,'transcript_id'])@listData$transcript_id),
                        transcript=as.character(elementMetadata(gtf[,'transcript_id'])@listData$transcript_id),
                        symbol=as.character(elementMetadata(gtf[,'gene_name'])@listData$gene_name))
geneModels = geneModels[geneModels$chromosome %in% paste0('chr',c(1:99,'X','Y')),]


save(geneModels, file=paste0("data/gm_", id, ".RData"))
