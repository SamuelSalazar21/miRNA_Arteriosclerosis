###https://ivanek.github.io/analysisOfGenomicsDataWithR/02_IntroToBioc_html.html
ls()
ir <- IRanges(start = c(3,10,50), end=c(17,28,104))
ir
gr <- GRanges(seqnames = c("chr1","chr1","chr2"),
              ranges=ir,
              strand = c("+","+","-"))
seqlevelsStyle(gr) <- "Ensembl"
seqlevelsStyle(gr)

gr[1:2]

gr+5
gr

seqlengths(gr)
seqlengths(gr) <- c("1"=100,"2"=200)
seqlengths(gr)
gr+5
trim(gr+5)
genome(gr) <- "conejito"
gr

shift(gr,5)

flank(gr, 2, start = T, both = F)
gr

mcols(gr)$score <- 1:3
gr
####EXcercises



url <- "C:/Users/PC/Downloads/gencode.v28.annotation.1.1.10M.gtf"
read.table(url, header=T)


gtf <- rtracklayer::import(url)
gtf
class(gtf)
str(gtf)
length(gtf)
gtf
genome(gtf)
seqnames(gtf)
gtsub <- subset(gtf, gene_name=="CHD5")
length(subset(gtsub, type=="transcript"))
###List of genomic ranges
library(GenomicRanges)

grl <- GRangesList(
  tx1.1=GRanges(seqnames="chr1",
  ranges=IRanges(start=c(50,125,188),
                 end=c(75,130,257)),
  strand="+",
  symbol="tx1.1", transcript="tx1.1",
  gene="gene1", exon=paste0("ex1.",1:3)),
  
  tx1.2=GRanges(seqnames = "chr1", 
                ranges=IRanges(start = c(50,125,200),
                               end=c(75,130,227)),
                strand = "+",
                symbol="tx1.2", transcript="tx1.2",
                gene="gene1", exon=paste0("ex1.", c(1,2,4))),
  tx2.1=GRanges(seqnames = "chr2", 
                ranges=IRanges(start = c(289,318),
                               end=c(300,350)),
                strand = "-",
                symbol="tx2.1", transcript="tx2.1",
                gene="gene2", exon=paste0("ex2.", c(1,2)))
  
  
  
)

seqlengths(grl) <- c(chr1=1000, chr2=500)
grl
library(Gviz)

g1 <- GeneRegionTrack(grl$tx1.1, showId=T, col=NULL,
                      fill="#377EB8", name="",
                      )
g2 <- GeneRegionTrack(grl$tx1.2, showId=T, col=NULL,
                      fill="#E41A1C", name="",
)



gtr <- GenomeAxisTrack()
tracks <- c(gtr,g1,g2)
plotTracks(tracks)

library(SummarizedExperiment)
se <- SummarizedExperiment(
  assays = list(exprs=matrix(10*runif(15),5,3)),
  colData = DataFrame(sample=paste0("Sample", 1:3),
                      condition=c("A","A","B")),
  rowData = DataFrame(gene=paste0("G", 1:5))
)
rownames(se) <- rowData(se)$gene
colnames(se) <- colData(se)$sample
se


assayNames(se)
assay(se, "exprs")
assay(se)
colData(se)
rowData(se)

sesub <- se[1:2,]
sesub
assay(sesub)
View(assay(se))

colData(sesub)
rowData(sesub)



se <- SummarizedExperiment(
  assays = list(exprs=matrix(10*runif(15), 5, 3)),
  colData = DataFrame(sample=paste0("Sample", 1:3),
                      condition=c("A","A","B")),
  rowRanges = GRanges(seqnames = "chr1",
                      ranges=IRanges(start=100*runif(5),
                                     width = 50),
                      strand = "+",
                      gene=paste0("G", 1:5))
)
class(se)
rowData(se)
colData(se)
rowRanges(se)
#####
library(airway)
data("airway")
airway
rowRanges(airway)
length(rowData(airway))
unique(colData(airway))

View(assay(airway))
class(airway)

assayNames(airway)
colSums(assay(airway, "counts"))

###Class inheritance
BiocManager::install("sloop")
library(SingleCellExperiment)
getClass("SingleCellExperiment")
getClass("RangedSummarizedExperiment")
showMethods(class="SingleCellExperiment")
library(sloop)
s4_methods_class("SingleCellExperiment")


library(Biostrings)
dna <- DNAString(x="AGGCATAGA")
dna
alphabet(DNAString())
IUPAC_CODE_MAP
alphabet(AAString())
AAString(x="EQI")

reverseComplement(dna)
translate(dna)
(dna_multiple <- DNAStringSet(c(x="AGGCATAGA", y="GGTATGAG")))
dna_multiple$y
dna_multiple[[2]]
dna_multiple[2]
width(dna_multiple)
Biostrings::subseq(dna_multiple, start = 1, end=3)
Biostrings::subseq(dna_multiple, start = 1, end=c(4,2))
Biostrings::subseq(dna_multiple, start=1, end=width(dna_multiple)-4)

DANA <- DNAString("AATTGGCCRGGCCAATT")
DANA


fasta <- Biostrings::readDNAStringSet("C:/Users/PC/Downloads/gencode.v28.transcripts.1.1.10M.fa")

fasta

min(width(fasta))
length(fasta)
range(width(fasta))
Biostrings::subseq(fasta, start=1, end=10)

freq<- alphabetFrequency(fasta)
fracs <- sweep(freq, MARGIN = 1, STATS = rowSums(freq), FUN = "/")
hist(fracs[,"G"])

LML <- DNAString(x="AGTGTTCGCGGGAGCGCCGCACCTACACCAGCCAACCCAGATCCCGAGGTCCGACAGCGCCCGGCCCAGATCCCCACGCCTGCCAGGAGCAAGCCGAGAGCCAGCCGGCCGGCGCACTCCGACTCCGAGCAGTCTCTGTCCTTCGACCCGAGCCCCGCGCCCTTTCCGGGACCCCTGCCC")
LML
translate(LML)
reverseComplement(LML)




###Annotation resources----
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db")
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
?TxDb.Hsapiens.UCSC.hg38.knownGene

gr <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
gr

grl <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by="gene")
grl



?org.Hs.eg.db
columns(org.Hs.eg.db)[12]
select(org.Hs.eg.db,keys = c("RAF1", "FLT3", "MAPK1"), 
       columns = c("ENTREZID", "ENSEMBL", "GENENAME"))

BiocManager::install("macrophage")       
library(tximeta)
library(macrophage)
library(SummarizedExperiment)

dir <- system.file("extdata", package="macrophage")
dir
list.files(dir)


coldata <- read.csv(file.path(dir, "coldata.csv"))[,c(1,2,3,5)]
View(coldata)
dim(coldata)

coldata$files <- file.path(dir, "quants", coldata$names, "quant.sf.gz")
head(coldata)       
all(file.exists(coldata$files))

se <- tximeta(coldata=coldata, type="salmon", dropInfReps=T)

list.files(file.path(dir, "quants", coldata$names[1]), recursive=T)

View(assay(se))


qf <- read.delim(file.path(dir, "quants", "SAMEA103884898", "quant.sf.gz"),
                 header=T, as.is=T)


all(rownames(se)==qf$Name)



colSums(assay(se, "abundance"))


gr <- GRanges(seqnames = "chr3",
              ranges=IRanges(start = 1e6,
                             end = 2e6),
                             strand="*")
gr
sesub <- IRanges::subsetByOverlaps(se, gr)  
sesub

rowData(se)

colData(se)

View(assay(se))



se2<- se[, se$condition_name %in% c("naive", "IFNg") ]

View(assay(se2))





all(colData(se)[, 4]==colData(se2)[,4])



rowRanges(se)



se2<- se[, se$condition_name %in% "naive"]
colData(se2)
se
se2 <- se[seqnames(rowRanges(se))=="chrX",
          se$condition_name %in% "naive"]
se2
dim(se2)
colData(se2)
rowRanges(se2)
##Summarizing on the gene level----
seg <- summarizeToGene(se)
seg
seg
rowData(seg)


length(unique(rowData(seg)))

rowData(seg)
rowRanges(seg)


rowRanges(IRanges::subsetByOverlaps(seg, gr))
library(org.Hs.eg.db)
seg <- addIds(seg, "SYMBOL")

head(rowData(seg))

AnnotationDbi::columns(org.Hs.eg.db)
seg2 <- addIds(se, "GO")
head(rowData(seg2))

rowData(seg2)


seg <- addIds(seg, GO, multiVals="list")
head <- (rowData(seg))
seg <- addIds(seg, "GO" )
rowData(seg)
seg <- addIds(seg, "ENTREZID")

sort(colSums(assay(seg, "counts")))
assay(seg)
sessionInfo()





library(DESeq2)
BiocManager::install("vsn2")
library(vsn2)
seg <- DESeqDataSet(seg, design = ~condition_name)
seg

seg <- seg[rowSums(assay(seg, "counts"))>0,]
seg


vst <- varianceStabilizingTransformation(seg, blind=T)
seg <- DESeq2::estimateSizeFactors(seg)
meanSdPlot(counts(seg, normalized=T))

meanSdPlot(assay(normTransform(seg)))
meanSdPlot(assay(vst))

rv <- rowVars(assay(vst))
sel <- order(rv, decreasing = T)[1:500]
pca <- prcomp(t(assay(vst)[sel,]))
pca

df <- data.frame(pca$x[, 1:2], colData(vst))
df
library(ggplot2)
ggplot(df, aes(x=PC1, y=PC2, color= condition_name))+
  geom_point(size=5) + theme_bw()

##mRNA seq 3: Differential Gene Expression with DESeq----


packages <- c("tidyverse", "DT", "tximeta", "macrophage", "SummarizedExperiment", "DESeq2")
str(packages)
class(packages)
library( "tximeta")

dir <- system.file("extdata", package = "macrophage")
coldata <- read.csv(file.path(dir, "coldata.csv")) %>% 
  dplyr::select(c(names, sample_id, line_id, condition_name))
coldata

coldata %>% group_by(condition_name) %>% 
  dplyr::summarise(n_replicates=n())



coldata <- coldata %>% mutate(files=file.path(dir, "quants", names, "quant.sf.gz"))
coldata



se <- tximeta::tximeta(coldata = coldata, type = "salmon", dropInfReps=T)

seg <- summarizeToGene(se)
rownames(seg)
rownames(seg) <- str_replace(rownames(seg), "\\.[:digit:]+$", "")
library(org.Hs.eg.db)
seg <- addIds(seg, "SYMBOL")
seg <- addIds(seg, "GO", multiVals="list")
head(rowData(seg))

seg$condition_name <- factor(seg$condition_name,
                             levels=c("naive", "IFNg", "SL1344", "IFNg_SL1344"))

seg$condition_name

dds <- DESeqDataSet(seg, design = ~condition_name)
dds
dds
row.names(dds)
nrow(dds)


keep <- rowSums(counts(dds))>=10
sum(keep)
dds <- dds[keep,]
nrow(dds)

dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)
dds <- DESeq(dds)
dds
plotDispEsts(dds)


res <- results(dds)
head(res)

mcols(res, use.names = T)[6,2]

summary(res)



table(res$padj < 0.1)

res.05 <- results(dds, alpha=0.05)

table(res.05$padj<5)


sum(res$padj<0.05, na.rm = T)
sum(res$pvalue<0.05, na.rm = T)

sum(!is.na(res$pvalue))


resSig <- subset(res, padj < 0.1)

head(resSig[order(resSig$log2FoldChange),], n=10)


head(resSig[order(resSig$log2FoldChange, decreasing = T),], n=10)


topGene <- rownames(res)[which.min(res$padj)]
topGene
plotCounts(dds, gene=topGene, intgroup = c("condition_name"))

geneCounts <- plotCounts(dds, gene=topGene, intgroup = c("condition_name"), returnData = T)
geneCounts

ggplot(geneCounts, aes(x=condition_name, y= count, color=condition_name))+
  geom_point(position=position_jitter(width = 0.1, height = 0, seed=123), size=4)+
             scale_y_log10()+
               labs(x="group", y="normalized count") +
  theme_bw(base_size = 12)

res.noshr <- results(dds, contrast = c("condition_name", "IFNg_SL1344", "naive"))
plotMA(res.noshr)

BiocManager::install("apeglm")
library("apeglm")
resultsNames(dds)

res <- lfcShrink(dds,
                 coef = "condition_name_IFNg_SL1344_vs_naive",
                 type="apeglm")

plotMA(res)
selGene <- "ENSG00000000938"
with(res[selGene,],{
  points(baseMean, log2FoldChange, col="black", cex=2, lwd=2)
  text(baseMean, log2FoldChange, selGene, pos=2, col="black")
}
     )

hist(res$pvalue[res$baseMean>1], breaks = 0:20/20,
     col="grey50", border = "white")


###Gene Clustering(heatmap)----

vsd <- varianceStabilizingTransformation(dds, blind=T)
BiocManager::install("genefilter")
library(genefilter)
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing=T), 20)
topVarGenes

mat <- assay(vsd)[topVarGenes,]
mat <- mat-rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("line_id", "condition_name")])
library(pheatmap)
library("AnnotationDbi")
library(org.Hs.eg.db)
pheatmap(mat, annotation_col = anno)


columns(org.Hs.eg.db)

res$symbol <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")

res$entrez <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = "ENTREZID",
                     keytype = "ENSEMBL",
                     multiVals = "first")

resOrdered <- res[order(res$pvalue),]
head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)
write.table(resOrderedDF, file="C:/Users/PC/Documents/results2.tsv", sep="\t")
getwd()


BiocManager::install("XML")

library(ReportingTools)
version


htmlRep <- HTMLReport(shortName = "reporte",
                      title = "Mi Reporte",
                      reportDirectory ="Documents/reporte" )


publish(resOrderedDF, htmlRep)

## Gene set enrichment analysis----

library(org.Hs.eg.db)

library(fgsea)

BiocManager::install("fgsea")
  library(reactome.db)  
library(stringr)
ens2eg <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = rownames(res),
                                columns="ENTREZID",
                                keytype = "ENSEMBL")
eg2reactome <- AnnotationDbi::select(reactome.db, keys=unique(ens2eg$ENTREZID),
                                     columns = c("PATHID", "PATHNAME"))

eg2reactome <- eg2reactome[!is.na(eg2reactome$PATHID),]

eg2reactome$PATHNAME <- str_replace(eg2reactome$PATHNAME,"^Homo sapiens: ", "" )

eg2reactome$REACTOME <- str_c(eg2reactome$PATHID, eg2reactome$PATHNAME, sep=": ")

eg2reactome <- eg2reactome[ , c("ENTREZID", "REACTOME")]


head(eg2reactome)



View(eg2reactome)

ens2reactome <- merge(ens2eg, eg2reactome)



View(ens2reactome)


pathwayList <- split(ens2reactome$ENSEMBL, ens2reactome$REACTOME)

str(pathwayList)


head(res$log2FoldChange)


set.seed(29)

fgseaRes <- fgsea(pathways=pathwayList,
                  stats=res$log2FoldChange,
                  eps= 0.0,
                  minSize=10,
                  maxSize=500)
fgseaRes







head(fgseaRes[order(pval),])



library(ggplot2)

plotEnrichment(pathwayList[["R-HSA-449147: Signaling by Interleukins"]], res$log2FoldChange)+
  labs(title = "R-HSA-449147: Signaling by Interleukins")


topPathwaysUp <- fgseaRes[ES>0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES<0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
topPathways
plotGseaTable(pathwayList[topPathways], res$log2FoldChange, fgseaRes, gseaParam = 0.5)


