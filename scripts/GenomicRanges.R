install.packages("installr")
library("installr")
updateR()
BiocManager::install("GenomicRanges")
library(GenomicRanges)
library(ggbio)


devtools::install_github("compgenomr/compGenomRData")
options(timeout=800)

gr=GRanges(seqnames = c("chr1", "chr2", "chr2"),
           ranges=IRanges(start=c(50,150,200),
                          end=c(100,200,300)),
           strand=c("+","-","-"), 
           names=c("id1", "id3", "id2"),
           scores=c("100","90","50"))
gr

mcols(gr)=DataFrame(name2=c("pax6", "meis1", "zic4"), score2=c("1","2","3"))
gr
values(gr)
gr$name3=c("A","B", "C")
gr

elementMetadata(gr)
#Read CpGi data set
filePath=system.file("extdata", "cpgi.hg19.chr21.bed",
                     package = "compGenomRData")



cpgi.df=read.table(filePath, header=FALSE,
                   stringsAsFactors=FALSE)

cpgi.df

cpgi.df=cpgi.df[grep("-", cpgi.df[,1], invert=TRUE),]


cpgi.df





cpgi.gr=GRanges(seqnames=cpgi.df[,1], 
                ranges=IRanges(start=cpgi.df[,2],
                               end=cpgi.df[,3]))
cpgi.gr
#Read refseq file
filePathRefseq=system.file("extdata",
                           "refseq.hg19.chr21.bed",
                           package="compGenomRData")
ref.df=read.table(filePathRefseq, header = FALSE, stringsAsFactors = FALSE)

ref.df
ref.gr=GRanges(seqnames=ref.df[,1], ranges=IRanges(start=ref.df[,2],
                                                   end=ref.df[,3]),
               strand=ref.df[,6], name=ref.df[,4])
ref.gr
ref.df

tss.gr=ref.gr

end(tss.gr[strand(tss.gr)=="+",])= start(tss.gr[strand(tss.gr)=="+",])
start(tss.gr[strand(tss.gr)=="-",])= end(tss.gr[strand(tss.gr)=="-",])

tss.gr=tss.gr[!duplicated(tss.gr),]
tss.gr
######
BiocManager::install("genomation")
require(rtracklayer)
import.bed(filePathRefseq)

##UCSC
session <- browserSession("UCSC",  url="http://genome-euro.ucsc.edu/cgi-bin")
session <- browserSession("UCSC",url = 'http://genome-euro.ucsc.edu/cgi-bin/')
genome(session) <- "mm9"
query <- ucscTableQuery(session, track="CpG Islands", table="cpgIslandExt",
                       range=GRangesForUCSCGenome("mm9", "chr12"))
session
query
track(query)
genome(session)
#####
library(genomation)
filePathPeaks=system.file("extdata",
                          "wgEncodeHaibTfbsGm12878Sp1Pcr1xPkRep1.broadPeak.gz",
                          package = "compGenomRData")
pk1.gr=readBroadPeak(filePathPeaks)
pk1.gr


elementMetadata(pk1.gr)
subsetByOverlaps(pk1.gr, cpgi.gr)

counts=countOverlaps(pk1.gr, cpgi.gr)
counts
findOverlaps(pk1.gr, cpgi.gr)
n.ind=nearest(pk1.gr, cpgi.gr)
n.ind
dists=distanceToNearest(pk1.gr, cpgi.gr, select="arbitrary")
dists
dist2plot=mcols(dists)[,1]
hist(log10(dist2plot), xlab="log10(distance to nearest TSS)",
     main="Distances")





##Dealing with mapped high-troughput sequencing read

promoter.gr=tss.gr
start(promoter.gr)=start(promoter.gr)-1000
end(promoter.gr)=end(promoter.gr)+1000
promoter.gr=promoter.gr[seqnames(promoter.gr)=="chr21"]
promoter.gr


library(Rsamtools)
bamfilePath=system.file("extdata",
                        "wgEncodeHaibTfbsGm12878Sp1Pcr1xAlnRep1.chr21.bam",
                        package="compGenomRData")
param <- ScanBamParam(which = promoter.gr)
param
counts <- countBam(bamfilePath, param = param)
counts

library(GenomicAlignments)
alns <- readGAlignments(bamfilePath, param=param)
alns

covs=coverage(alns)


covs












bwFile=system.file("extdata", "wgEncodeHaibTfbsA549.chr21.bw", 
                   package = "compGenomRData")
bw.gr=import(bwFile, which=promoter.gr)
cov.bw=coverage(bw.gr, weight = "score")
cov.bw
cov.bw2=import(bwFile, which=promoter.gr, as="RleList")
cov.bw2






myViews=Views(cov.bw, as(promoter.gr, "IRangesList"))
myViews
myViews[[1]]

plot(myViews[[1]][[5]])

head(viewMeans(myViews[[1]]))
head(viewMaxs(myViews[[1]]))


#SummarizedExperiment class----
#simulate an RNA-seq read counts table
nrows <- 200
ncols <- 6
counts <- matrix(runif(nrows*ncols, 1, 1e4), nrows)
counts
1e4
#create gene location
rowRanges <- GRanges(rep(c("chr1","chr2"), c(50,150)),
                     IRanges(floor(runif(200,1e5,1e6)), width=100),
                     strand=sample(c("+","-"), 200, TRUE),
                     feature_id=paste0("gene",1:200 ))
rowRanges

colData <- DataFrame(timepoint=1:6,
                     row.names = LETTERS[1:6])
colData
se=SummarizedExperiment(assays = list(counts=counts),
                        rowRanges = rowRanges, colData = colData)
se
colData(se)
rowData(se)
assays(se)
assays(se)[[1]]
se[1:5, 1:3]
colData(se[1:5, 1:3])
rowData(se[1:5, 1:3])
assays(se[1:5, 1:3])
se[,se$timepoint==1]
roi <- GRanges(seqnames = "chr1", ranges=100000:1100000)
roi
subsetByOverlaps(se,roi)

cpgi.track=AnnotationTrack(cpgi.gr, name = "CpG")
cpgi.track

gene.track <- BiomartGeneRegionTrack(genome="hg19",
                                     chromosome="chr21",
                                     start=27698681, end=28083310, 
                                     name="ENSEMBL")
gene.track

chipseqFile=system.file("extdata", "wgEncodeHaibTfbsA549.chr21.bw",
                        package="compGenomRData")
cov.track=DataTrack(chipseqFile, type="1",
                    name="coverage")
cov.track
track.list=list(cpgi.track, gene.track, cov.track)
plotTracks(list(cpgi.track, gene.track, cov.track), from = 27698681, to=28083310)








#Summaries of genomic intervals on multiple loci
library(genomation)
transcripFile=system.file("extdata", "refseq.hg19.chr20.bed", 
                          package="compGenomRData")



feat=readTranscriptFeatures(transcripFile, remove.unusual = TRUE,
                            up.flank = 500, down.flank = 500)
feat
prom=feat$promoters
prom
H3K4me3File=system.file("extdata",
                        "H1.ESC.H3K4me3.chr20.bw",
                        package = "compGenomRData")

sm=ScoreMatrix(H3K4me3File, prom,
               type="bigwig", strand.aware=TRUE)
plotMeta(sm, profile.names = "H3K4me3", xcoords = c(-500,500),
         ylab="H3k4me enrichment", dispersion="se",
         xlab="bases around TSS")

heatMatrix(sm, order = TRUE, xcoords = c(-500,500),
           xlab="bases around TSS")


DNAseFile=system.file("extdata", "H1.ESC.dnase.chr20.bw",
                      package="compGenomRData")
sml=ScoreMatrixList(c(H3K4me3=H3K4me3File,
                      DNase=DNAseFile), prom,
                    type="bigwig", strand.aware = TRUE)
plotMeta(sml)

set.seed(1029)
multiHeatMatrix(sml, order=TRUE, xcoords=c(-500, 500),
                xlab="bases around TSS", winsorize = c(0,95),
                matrix.main=c("H3K4me3", "DNAse"),
                column.scale = TRUE,
                clustfun=function(x) kmeans(x, centers = 3)$cluster)

#ggbio
library(ggbio)
BiocManager::install("genomation")
ideoCyto
p <- autoplot(seqinfo(ideoCyto$hg19), layout="karyogram")
CpGiFile=system.file("extdata", "CpGi.hg19.table.txt", package = "compGenomRData")

cpgi.gr=genomation::readGeneric(CpGiFile, chr=1, start=2, end=3,
                               header=T, keep.all.metadata=T,
                               remove.unusual=T)
p+layout_karyogram(cpgi.gr)


p+layout_karyogram(cpgi.gr, aes(x=start, y=obsExp), geom="point",
                   ylim=c(2,50), color="red",
                   size=0.1, rect.height=1)
p <- ggplot()+layout_circle(ideoCyto$hg19, geom="ideo", fill="white", cytoband=T,
                            radius=39, trackWidth=2)
p <- p+ layout_circle(cpgi.gr, geom="point", grid=T,
                      size=0.01, aes(y=obsExp), color="red",
                      radius=42, trackWidth=10)
p <- p+layout_circle(as(seqinfo(ideoCyto$hg19), "GRanges"),
                   geom="text", aes(label=seqnames),
                   vjust=0, radius=55, trackWidth=7, size=3)

p <- p + layout_circle(as(seqinfo(ideoCyto$hg19),"GRanges"), 
                       geom = "text", aes(label = seqnames), 
                       vjust = 0, radius = 55, trackWidth = 7,
                       size=3) 

p
