
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


## ------------------------------------------------------------------------
require(pasillaBamSubset)
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)


## ------------------------------------------------------------------------
txdb <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
grl <- exonsBy(txdb, by="gene")
grl[100]
grl[[100]]
grl[[100]][1]


## ------------------------------------------------------------------------
fl1 <- untreated1_chr4()
fl2 <- untreated3_chr4()
fl1


## ------------------------------------------------------------------------
require (Rsamtools)
require (GenomicRanges)


## ------------------------------------------------------------------------
require(GenomicAlignments)


## ------------------------------------------------------------------------
fls <- BamFileList(c(fl1, fl2), yieldSize=5e4)
names(fls) <- c("first","second")


## ------------------------------------------------------------------------
so1 <- summarizeOverlaps(features=grl,
                         reads=fls,
                         ignore.strand=TRUE)
so1


## ------------------------------------------------------------------------
head(assay(so1))
colSums(assay(so1))


## ------------------------------------------------------------------------
rowData(so1)
colData(so1)
colData(so1)$sample <- c("one","two")
colData(so1)
metadata(rowData(so1))


## ------------------------------------------------------------------------
x <- assay(so1)[,1]
hist(x[x > 0], col="grey")
hist(x[x > 0 & x < 10000], col="grey")
plot(assay(so1) + 1, log="xy")


## ------------------------------------------------------------------------
fls <- BamFileList(fl2, yieldSize=5e4)
so2 <- summarizeOverlaps(features=grl,
                         reads=fls,
                         ignore.strand=TRUE,
                         singleEnd=FALSE, 
                         fragments=TRUE)
colSums(assay(so2))
colSums(assay(so1))
plot(assay(so1)[,2], assay(so2)[,1], xlim=c(0,5000), ylim=c(0,5000),
     xlab="single end counting", ylab="paired end counting")
abline(0,1)
abline(0,.5)


