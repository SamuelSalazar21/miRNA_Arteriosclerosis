if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c('qvalue','plot3D','ggplot2','pheatmap','cowplot', 
                       'cluster', 'NbClust', 'fastICA', 'NMF','matrixStats',
                       'Rtsne', 'mosaic', 'knitr', 'genomation',
                       'ggbio', 'Gviz', 'DESeq2', 'RUVSeq',
                       'gProfileR', 'ggfortify', 'corrplot',
                       'gage', 'EDASeq', 'citr', 'formatR',
                       'svglite', 'Rqc', 'ShortRead', 'QuasR',
                       'methylKit','FactoMineR', 'iClusterPlus',
                       'enrichR','caret','xgboost','glmnet',
                       'DALEX','kernlab','pROC','nnet','RANN',
                       'ranger','GenomeInfoDb', 'GenomicRanges',
                       'GenomicAlignments', 'ComplexHeatmap', 'circlize', 
                       'rtracklayer', 'BSgenome.Hsapiens.UCSC.hg38',
                       'BSgenome.Hsapiens.UCSC.hg19','tidyr',
                       'AnnotationHub', 'GenomicFeatures', 'normr',
                       'MotifDb', 'TFBSTools', 'rGADEM', 'JASPAR2018'
)) #instalación de paquetes de bioinformatica
install.packages("devtools") #instalación de paquetes devtools
devtools::install_github("compgenomr/compGenomRData") #intalaciión de paquete desde un repositorio de git hub
options(timeout=800)

enhancerFilePath=system.file("extdata",
                             "subset.enhancers.hg18.bed",
                             package="compGenomRData")
cpgiFilePath=system.file("extdata",
                         "subset.cpgi.hg18.bed",
                         package="compGenomRData")
# read enhancer marker BED file
enh.df <- read.table(enhancerFilePath, header = FALSE) 

# read CpG island BED file
cpgi.df <- read.table(cpgiFilePath, header = FALSE) 

# check first lines to see how the data looks like
head(enh.df)
write.table(cpgi.df,file="cpgiexample.txt",quote=FALSE,
            row.names=FALSE,col.names=FALSE,sep="\t")
#################Manipulate GRanges
library(GenomicRanges)
gr=GRanges(seqnames =c("chr1","chr2","chr3"),
           ranges=IRanges(start = c(50,150,200),
                          end=c(100,200,300)),
           strand = c("+","-","-"))
gr=GRanges(seqnames =c("chr1","chr2","chr2"),
           ranges=IRanges(start = c(50,150,200),
                          end=c(100,200,300)),
           names=c("id1","id3","id2"),
           svores=c(100,90,50))
gr
elementMetadata(gr)
values(gr)
gr$name3=c("A","B","C")
gr

hola=2
# read CpGi data set____________
  
  
filePath=system.file("extdata","cpgi.hg19.chr21.bed",
                      package="compGenomRData")


cpgi.df = read.table(filePath, header = FALSE,
                     stringsAsFactors=FALSE) 
# remove chr names with "_"
cpgi.df =cpgi.df [grep("_",cpgi.df[,1],invert=TRUE),]

cpgi.gr=GRanges(seqnames=cpgi.df[,1],
                ranges=IRanges(start=cpgi.df[,2],
                               end=cpgi.df[,3]))
cpgi.gr


filePathRefseq=system.file("extdata",
                           "refseq.hg19.chr21.bed",
                           package="compGenomRData")

ref.df=read.table(filePathRefseq, header = F,stringsAsFactors = F)
ref.gr=GRanges(seqnames = ref.df[,1],
               ranges=IRanges(start = ref.df[,2],end = ref.df[,3]),
               strand=ref.df[,6], name=ref.df[,4])
ref.gr
#Get TSS
tss.gr=ref.gr
tss.gr
end(tss.gr[strand(tss.gr)=="+",]) = start(tss.gr[strand(tss.gr)=="+",])
start(tss.gr[strand(tss.gr)=="-",]) = end(tss.gr[strand(tss.gr)=="-",])
tss.gr=tss.gr[!duplicated(tss.gr)]
tss.gr
######--------------------------
require(rtracklayer)
import.bed(filePathRefseq)






####UCSC------
session <- browserSession("UCSC", url="http://genome-euro.ucsc.edu/cgi-bin/")
genome(session) <- "mm9"


query <- ucscTableQuery(session, track="CpG Islands",table="cpgIslandExt", 
                        range=GRangesForUCSCGenome("mm9","chr12"))
track(query)

library(genomation)
filePathPeaks=system.file("extdata",   
                          "wgEncodeHaibTfbsGm12878Sp1Pcr1xPkRep1.broadPeak.gz",
                          package="compGenomRData")
pk1.gr=readBroadPeak(filePathPeaks)
pk1.gr
subsetByOverlaps(pk1.gr,cpgi.gr)
counts=countOverlaps(pk1.gr, cpgi.gr)
head(counts)
findOverlaps(pk1.gr,cpgi.gr)
n.ind=nearest(pk1.gr, cpgi.gr)
n.ind
dist=distanceToNearest(pk1.gr, cpgi.gr, select="arbitrary")
dist
dist2plot=mcols(dist)[,1]
hist(log10(dist2plot), xlab="Log 10 nearest", main="Distances")
###Counting mapped reads for a set of regions------
promoter.gr=tss.gr
promoter.gr

start(promoter.gr)=start(promoter.gr)-1000
end(promoter.gr)=end(promoter.gr)-1000
promoter.gr=promoter.gr[seqnames(promoter.gr)=="chr21"]
promoter.gr
###
library(Rsamtools)
bamfilePath=system.file("extdata",
                        "wgEncodeHaibTfbsGm12878Sp1Pcr1xAlnRep1.chr21.bam",
                        package="compGenomRData")
param <- ScanBamParam(which = promoter.gr)
counts=countBam(bamfilePath, param=param)
counts
library(GenomicAlignments)
alns <- readGAlignments(bamfilePath, param=param)
alns

#Rle vector------
covs=coverage(alns)
covs
library(rtracklayer)
exdata
list.dir("extdata")
bwFiles=system.file("extdata",
                    "wgEncodeHaibTfbsA549.chr21.bw", package="compGenomRData" )
bw.gr=import(bwFiles, which=promoter.gr)
bw.gr
cov.bw=coverage(bw.gr, weight = "score")
cov.bw
#Extracting subsections of Rle and Rle list objects----
myViews=Views(cov.bw, as(promoter.gr, "IRangesList"))
myViews
plot(myViews[[1]][[5]])

head(viewMeans(myViews[[1]]))
head(viewMaxs(myViews[[1]]))
#Sumarzed Experiment class----
library(SummarizedExperiment)
nrows <- 200
ncols <- 6

counts <- matrix(runif(nrows*ncols,1,1e4), nrows)
counts
length(counts)
colData <- data.frame(timepoint=1:6, row.names=LETTERS[1:6])
colData
rowRanges <- GRanges(rep(c("chr1","chr2"), c(50,150)),
                     IRanges(floor(runif(200,1e5,1e6)), width=100),
                     strand = sample(c("+","-"), 200, TRUE),
                     feature_id=paste0("gene",1:200))
se=SummarizedExperiment(assays = list(counts=counts), 
                        rowRanges = rowRanges, colData = colData)
se

colData(se)
rowData(se)
assays(se)
assays(se)[[1]]
assays(se)$counts


###Subset


assays(se[1:5, 1:3])$counts
se[,se$timepoint==1]

roi <- GRanges(seqnames = "chr1", ranges=100000:1100000)
roi
subsetByOverlaps(se, roi)
###Gviz
library(Gviz)
cpgi.track=AnnotationTrack(cpgi.gr, name="CpG")
gene.track <- BiomartGeneRegionTrack(genome="hg19",
                                     chromosome = "chr21",
                                     start=27698681, end=28083310,
                                     name="ENSEMBL")
gene.track

chipseqFile <- system.file("extdata", "wgEncodeHaibTfbsA549.chr21.bw",
                           package = "compGenomRData")

dir("extdata")

cov.track=DataTrack(chipseqFile, 
                    name="coverage")

cov.track

plotTracks(list( cpgi.track,gene.track, cov.track), from=27698681, to=28083310,
           chromosome = "chr21")



