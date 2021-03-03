
## ----options, echo=FALSE-------------------------------------------------
opts_chunk$set(fig.path=paste0("figure/", sub("(.*).Rmd","\\1",basename(knitr:::knit_concord$get('infile'))), "-"))


## ------------------------------------------------------------------------
source("http://bioconductor.org/biocLite.R")
installifnot <- function (packageName){
  if (!(require(packageName, character.only=TRUE))) biocLite(packageName)
}
installifnot("pasillaBamSubset")
installifnot("TxDb.Dmelanogaster.UCSC.dm3.ensGene")
installifnot("Rsamtools")
installifnot("GenomicRanges")
installifnot("GenomicAlignments")
installifnot("biomaRt")
installifnot("GenomicFeatures")
installifnot("Gviz")
installifnot("ggbio")


## ------------------------------------------------------------------------
require (pasillaBamSubset)
fl1 <- untreated1_chr4()
fl2 <- untreated3_chr4()


## ------------------------------------------------------------------------
file.copy(from=fl1,to=basename(fl1))
file.copy(from=fl2,to=basename(fl2))
require(Rsamtools)
indexBam(basename(fl1))
indexBam(basename(fl2))


## ------------------------------------------------------------------------
require(GenomicRanges)


## ------------------------------------------------------------------------
require(GenomicAlignments)


## ------------------------------------------------------------------------
x <- readGAlignments(fl1)
xcov <- coverage(x)
z <- GRanges("chr4",IRanges(456500,466000))
# Bioconductor 2.14
xcov[z]
# Bioconductor 2.13
xcov$chr4[ranges(z)]
xnum <- as.numeric(xcov$chr4[ranges(z)])
plot(xnum)


## ------------------------------------------------------------------------
y <- readGAlignmentPairs(fl2)
ycov <- coverage(y)
ynum <- as.numeric(ycov$chr4[ranges(z)])
plot(xnum, type="l", col="blue", lwd=2)
lines(ynum, col="red", lwd=2)


## ------------------------------------------------------------------------
plot(xnum, type="l", col="blue", lwd=2, xlim=c(6200,6600))
lines(ynum, col="red", lwd=2)


## ------------------------------------------------------------------------
require(biomaRt)
m <- useMart("ensembl", dataset = "dmelanogaster_gene_ensembl")
lf <- listFilters(m)
lf[grep("name", lf$description, ignore.case=TRUE),]
map <- getBM(mart = m,
  attributes = c("ensembl_gene_id", "flybasename_gene"),
  filters = "flybasename_gene", 
  values = "lgs")
map


## ------------------------------------------------------------------------
require(GenomicFeatures)
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
grl <- exonsBy(TxDb.Dmelanogaster.UCSC.dm3.ensGene, by="gene")
gene <- grl[[map$ensembl_gene_id[1]]]
gene


## ------------------------------------------------------------------------
rg <- range(gene)
plot(c(start(rg), end(rg)), c(0,0), type="n", xlab=seqnames(gene)[1], ylab="")
arrows(start(gene),rep(0,length(gene)),
       end(gene),rep(0,length(gene)),
       lwd=3, length=.1)


## ------------------------------------------------------------------------
plot(c(start(rg), end(rg)), c(0,0), type="n", xlab=seqnames(gene)[1], ylab="")
arrows(start(gene),rep(0,length(gene)),
       end(gene),rep(0,length(gene)),
       lwd=3, length=.1, 
       code=ifelse(as.character(strand(gene)[1]) == "+", 2, 1))


## ------------------------------------------------------------------------
require(Gviz)
gtrack <- GenomeAxisTrack()
atrack <- AnnotationTrack(gene, name = "Gene Model")
plotTracks(list(gtrack, atrack))


## ------------------------------------------------------------------------
xgr <- as(xcov, "GRanges")
ygr <- as(ycov, "GRanges")
dtrack1 <- DataTrack(xgr[xgr %over% z], name = "sample 1")
dtrack2 <- DataTrack(ygr[ygr %over% z], name = "sample 2")
plotTracks(list(gtrack, atrack, dtrack1, dtrack2))
plotTracks(list(gtrack, atrack, dtrack1, dtrack2), type="polygon")


## ------------------------------------------------------------------------
require (ggbio)
autoplot(gene)
autoplot(fl1, which=z)
autoplot(fl2, which=z)


