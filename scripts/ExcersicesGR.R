#tutorial de Rtraclayer

dir(system.file("extdata", package = "compGenomRData"))


Rango <- GRanges(seqnames = c("chr1", "chr1", "chr2"), 
                 IRanges(start = c(1000,11100,20000), end = c(10300,11500,20030)
                        ),  strand=c("+","-","+"), score=c(10,20,15))
Rango

(Rango2 <- Rango[strand(Rango)=="+",])
Rango3 <- Rango[seqnames(Rango)=="chr1",]
Rango3
require(rtracklayer)
session <- browserSession("UCSC",url = 'http://genome-euro.ucsc.edu/cgi-bin/')
genome(session) <- "mm9"

Cpgislandquery <- ucscTableQuery(session, track="CpG Islands",table="cpgIslandExt",
                                 range=GRangesForUCSCGenome("mm9", "chr12"))
cpgi <- track(Cpgislandquery)
RefseqGenes <- ucscTableQuery(session, track="RefSeq Genes",
                                       table="refGene",
                                       range=GRangesForUCSCGenome("mm9", "chr12"))
RefseqGenes
Cpgislandquery

#track(RefseqGenes) #don't works
RefseqGenes


tab1<- getTable(RefseqGenes)
head(tab1)
ref.gr1=GRanges(seqnames=tab1[,3], ranges=IRanges(start=tab1[,5],
                                                   end=tab1[,6]),
               strand=tab1[,4], name=tab1[,2])

ref.gr1
promotores <- promoters(ref.gr1)
promotores
?promoters

subsetByOverlaps(promotores,cpgi)
counts=countOverlaps(promotores,cpgi)
counts
sum(counts)
length(promotores)
porcentaje <- sum(counts)/length(promotores)
hist(counts)

####6
filePeaksp1=system.file("extdata",
                           "wgEncodeHaibTfbsGm12878Sp1Pcr1xPkRep1.broadPeak.gz",
                           package="compGenomRData")


filePeaksp2=system.file("extdata",
                        "wgEncodeHaibTfbsGm12878Sp1Pcr1xPkRep2.broadPeak.gz",
                        package="compGenomRData")



ref.pks1=read.table(filePeaksp1, header = FALSE, stringsAsFactors = FALSE)
ref.pks1


GR.pks1=GRanges(seqnames=ref.pks1[,1], ranges=IRanges(start=ref.pks1[,2],
                                                   end=ref.pks1[,3]),
               strand="*", name=ref.pks1[,4])
GR.pks1


sp1<- genomation::readBroadPeak(filePeaksp1)
sp2 <- genomation::readBroadPeak(filePeaksp2)
sp1 <- sp1[seqnames(sp1)=="chr21",]
sp2 <- sp2[seqnames(sp2)=="chr21",]
sp2
spcanonical <- subsetByOverlaps(sp1,sp2)
spcanonical

require(rtracklayer)
session <- browserSession("UCSC",url = 'http://genome-euro.ucsc.edu/cgi-bin/')
genome(session) <- "anoGam1"
ucscGenomes()
BiocManager::install("Gviz")

library(Gviz)

