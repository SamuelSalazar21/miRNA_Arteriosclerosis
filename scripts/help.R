queryi<- IRanges(c(1,65,66,300),c(5,80,80,500)) 
queryi
subjecti<- IRanges(c(5,30,100,105), c(10,60,104,200)) 
subjecti
subjecti; queryi
nearest(queryi,subjecti)
findOverlaps(queryi,subjecti)
countOverlaps(queryi,subjecti)
subsetByOverlaps(queryi,subjecti)
distanceToNearest(queryi,subjecti, select="arbitrary")


##########
colData2 <- data.frame(Muestra=factor(1:6),
                       condición=factor(c("A","A","B","B","C","C")),
                       tratamiento=factor(rep(0:1,3)))
colData2
exprs <- matrix(rnorm(6*10), ncol=6, nrow = 10)
exprs

BiocManager::install("EnsDb.Hsapiens.v86")
library("EnsDb.Hsapiens.v86")
txbd <- EnsDb.Hsapiens.v86

g <- genes(txbd)
g <- keepStandardChromosomes(g,pruning.mode = "coarse")
g
rowRanges2 <- g[1:10]
rowRanges2
se2 <- SummarizedExperiment(assay=list(exprs=exprs),
                            colData = colData2,
                            rowRanges = rowRanges2)
se2
###SummarizedExperiment
BiocManager::install("airway")
data(airway, package="airway")
sea <- airway
sea
assays(sea)$counts
rowRanges(sea)
colData(sea)
sea[, sea$dex=="trt"]

colData(sea[, sea$dex=="trt"])
colData(sea[, sea$cell=="N61311"])
metadata(sea)

###Create a summarized experiment
nrows <- 200
ncols <- 6
RowRanges <- GRanges(rep(c("chr1","chr2"), c(50,150)),
                  IRanges(floor(runif(200,1e5,1e6)), width = 100),
                  strand = sample(c("+","-"), 200, T),
                  feature_id=sprintf("ID%03d", 1:200))
RowRanges
counts <- matrix(runif(nrows*ncols,1,1e4), nrows)
counts
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names = LETTERS[1:6])
colData

SummarizedExperiment(assays = list(counts=counts),colData = colData, 
                     rowRanges = RowRanges)

###subsetting
#subset the first five transcripts and first three samples


sea[1:5,1:3]
sea[ ,sea$cell=="N61311"]
###Getters and setters

counts <- matrix(1:15, 5,3, dimnames = list(LETTERS[1:5], LETTERS[1:3]))
dates <- SummarizedExperiment(assays = list(counts=counts),
                              rowData = DataFrame(month=month.name[1:5], day=1:5))
dates[rowData(dates)$month=="January",]
assays(sea)
assays(sea)[[1]][1:5,1:5]
assay(sea)[1:5, 1:5]
assay(sea)

roi <- GRanges(seqnames = "1", ranges=100000:1100000)
subsetByOverlaps(sea,roi)
##GVIZ----
library(Gviz)
data("cpgIslands")
class(cpgIslands)
chr <- as.character(unique(seqnames(cpgIslands)))
chr
gen <- genome(cpgIslands)
gen
atrack <- AnnotationTrack(cpgIslands, name="CpG")
atrack
plotTracks(atrack)

gtrack <- GenomeAxisTrack()
plotTracks(list(gtrack, atrack))
itrack <- IdeogramTrack(genome=gen, chromosome = chr)
plotTracks(list(itrack,gtrack, atrack))
data("geneModels")
grtrack <- GeneRegionTrack(geneModels, genome=gen, chromosome=chr,
                           name="Gene Model")
plotTracks(list(itrack, gtrack, atrack, grtrack))

plotTracks(list(itrack,gtrack, atrack, grtrack),
           from=26700000, to= 26750000)

plotTracks(list(itrack,gtrack, atrack, grtrack),
           extend.left = 0.5, extend.right= 100000, col=NULL)
library(BSgenome.Hsapiens.UCSC.hg19)

strack <- SequenceTrack(Hsapiens, chromosome = chr)
strack


plotTracks(list(itrack,gtrack, atrack, grtrack, strack),
           from = 26591822, to= 26591852, cex=0.8)


set.seed(255)

lim <- c(26700000, 26750000)

coords <- sort(c(lim[1],
                 sample(seq(from=lim[1], to=lim[2]),99),
                 lim[2]))
coords

dat <- runif(100, min = -10, max = 10)
dtrack <- DataTrack(data=dat, start=coords[-length(coords)],
                    end=coords[-1], chromosome=chr, genome=gen,
                    name="uniform")
dtrack




plotTracks(list(itrack,gtrack, atrack, grtrack, dtrack),
           from = lim[1], to= lim[2])


plotTracks(list(itrack, atrack, grtrack, dtrack),
           from = lim[1], to= lim[2], type="heatmap")

grtrack <- GeneRegionTrack(geneModels, genome=gen, chromosome=chr,
                          name="Modelo genetico",
                          transcriptAnnotation="symbol",
                          background.title="black")
grtrack

displayPars(grtrack) <- list(background.panel="blue", col=NULL)
plotTracks(list(itrack,gtrack, atrack , grtrack))

plotTracks(list(itrack,gtrack, atrack , grtrack),
           background.pane="yellow", background.title="green")

dp <- available.genomes()
dp <- availableDisplayPars(itrack)
dp


###Gprofiler2----
gostres <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = TRUE)

names(gostres)
View(head(gostres$result, 3))

get_version_info(organism = "hsapiens")
