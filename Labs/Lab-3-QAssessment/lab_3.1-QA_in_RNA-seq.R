## ----include=FALSE-------------------------------------------------------
opts_chunk$set(fig.path='images/grafic', tidy=FALSE)


## ----InstallIfNot--------------------------------------------------------
installifnot <- function (pkg){
  if (!require(pkg, character.only=T)){
  source("http://bioconductor.org/biocLite.R")
  biocLite(pkg)
}else{
  require(pkg, character.only=T)
  }
}
installifnot("ShortRead")
installifnot("Rsubread")
installifnot("Rsamtools")
installifnot("GenomicFeatures")



## ----getData-------------------------------------------------------------
require(ShortRead)
fq<-readFastq('Pool1small.fastq')


## ----checkReads----------------------------------------------------------
fq
sread(fq)
sread(fq)[1:10]
quality(fq)
# width(fq)
detail(fq)


## ----checkQuality--------------------------------------------------------
qaReads<-qa(fq,lane='Pool1small')
#qaReads<-qa('.',pattern='Pool1small.fastq',type='fastq')
show(qaReads)
qaReads[['readCounts']] # Number of reads
qaReads[['baseCalls']]  # Frequencies of each base
head(qaReads[["frequentSequences"]], n=5)  # Frequent sequences
qaReads[['perCycle']]$baseCall # Base Call per cycle
qaReads[['perCycle']]$quality[1:50,]  # Quality per cycle

# nucleotide distribution

axc<-alphabetByCycle(sread(fq))
axc[,1:10]
axc[DNA_BASES, 1:10]

axcp <- axc/ colSums(axc)
axcp[DNA_BASES, 1:10]

axc2<-prop.table(axc, margin=2)*100

axc2<-axc2[DNA_BASES,]

tables(fq)

# nucleotide with maximum frequency per cycle and overrepresented sequence

paste(DNA_BASES[apply(axc, 2, which.max)],
      collapse="")

cnt <- vcountPattern("ATTTAAA", sread(fq))
sum(cnt > 0)


## ----plotQuality---------------------------------------------------------
matplot(t(axc2), type="l", xlab="Cycle",
        ylab="Nt freq", lwd=2, lty=seq(along=DNA_BASES),
        col=seq(along=DNA_BASES))
legend("topright", DNA_BASES, lty=seq(along=DNA_BASES),
       col=seq(along=DNA_BASES), lwd=3, inset=0.01)


## ----preProc-------------------------------------------------------------
# Sum of qualities per read
qsr<-alphabetScore(quality(fq))
qsr[1:10]
# Mean quality per read
qar<-qsr/width(fq)
qar[1:10]
# select reads with mean quality greater or equal to 20
fq20q<-fq[qar>=20]
length(fq20q)

axq<-alphabetByCycle(quality(fq))
axq[,1:10]


# select reads of 72 bp
fq72b<-fq[width(fq)==72]
length(fq72b)

# mean quality per cycle
qxc <- as(quality(fq72b), "matrix")
dim(qxc)
qxc[1:5, 1:20]


## ----plotQual------------------------------------------------------------
plot(colMeans(qxc), type="b", lwd=3, xlab="Cycle",
     ylab="Quality score")


## ----filterReads---------------------------------------------------------
filt<-nFilter(threshold=0L)
fqn<-fq[filt(fq)]
length(fqn)
axc<-alphabetByCycle(sread(fqn))
axc[,1:10]


## ----removeNonNuc--------------------------------------------------------
fqnet<-clean(fq)
length(fqnet)
axc<-alphabetByCycle(sread(fqnet))
axc[,1:10]

# keep only first 30 bases
fq30p<-narrow(fq72b,start=1,end=30)
sread(fq30p)
qxc <- as(quality(fq30p), "matrix")
plot(colMeans(qxc), type="b", lwd=3, xlab="Cycle",ylab="Quality score")


# remove sequences resembling a polyA
head(qaReads[["frequentSequences"]], n=5)  # Frequent sequences
distance<-srdistance(fq,polyn('A',width(fq)[1]))[[1]]
# histogram of distances
hist(distance,xlab=paste('Distance to',polyn('A',width(fq)[1])),ylab='Number of reads')

# create a mask to select reads that resemble the sequence
polyAs<-distance<30
head(polyAs)
sum(polyAs)
fqnoPA<-fq[!polyAs]
length(fqnoPA)


## ----writeFastQ, eval=FALSE----------------------------------------------
## # save in fastq format
## writeFastq(fqnoPA,file='fqnoPA.fastq')


## ----indexAndAlign, eval=FALSE-------------------------------------------
## require(Rsubread)
## buildindex(basename='subreadIndex',reference='chr6.fa')
## align(index='subreadIndex',readfile1='Pool1small.fastq',+output_file='Pool1small.sam')
## propmapped('Pool1small.sam')


## ----SAM2BAM-------------------------------------------------------------
require(Rsamtools)
asBam(file='Pool1small.sam',destination='Pool1small')


## ----readBAM-------------------------------------------------------------
require(ShortRead)
aln<-readAligned('.','Pool1small.bam$',type='BAM')
aln


## ----getAnnotsAnnots, cache=TRUE-----------------------------------------
require(GenomicFeatures) 
supportedUCSCtables()
sgdTx<-makeTranscriptDbFromUCSC(genome='hg19',tablename='knownGene')
# saveFeatures(sgdTx,file='sgdTx.db.sqlite')
saveDb(sgdTx,file='sgdTx.db.sqlite')


## ----getAnnots-----------------------------------------------------------
require(GenomicFeatures) 
sgdTx.db<-loadDb('sgdTx.db.sqlite')
## 1st set all the sequences to be inactive:
isActiveSeq(sgdTx.db)[seqlevels(sgdTx.db)] <- FALSE ## Then set only "chr6" to be active: 
isActiveSeq(sgdTx.db)["chr6"] <- TRUE ## You can see which are active like this: isActiveSeq(sgdTx.db)
isActiveSeq(sgdTx.db)

# Extract exons grouped by genes
sgdGenes<-exonsBy(sgdTx.db,by='gene')
length(sgdGenes)
sgdGenes[1:3]

# Read Gapped alignments
require(GenomicRanges)
aln<-readGAlignments('Pool1small.bam') # GappedAlignments objects
# verify equality of names
seqlevels(aln)
seqlevels(sgdGenes)

ov<-summarizeOverlaps(sgdGenes,aln,ignore.strand=T)
head(assays(ov)$counts)

Pool1ov<-ov
TaulaPool1Chr6<-assays(Pool1ov)$counts[,1]
write.table(TaulaPool1Chr6,file='TaulaPool1.txt',quote=F,)


