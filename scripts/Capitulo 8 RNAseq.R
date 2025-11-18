#CPM

counts_file <-system.file("extdata/rna-seq/SRP029880.raw_counts.tsv",
                          package = "compGenomRData")

coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv",
                            package="compGenomRData")
counts <- as.matrix(read.table(counts_file, header=T, sep="\t"))
counts

head(counts)
length(unique(counts[,1]))

summary(counts[,1:11])
cpm <- apply(subset(counts, select=c(-width)),2,
             function(x)x/sum(as.numeric(x))*10^6)
head(cpm)
colSums(cpm)

geneLengths <- as.vector(subset(counts, select=c(width)))
geneLengths

rpkm <- apply(X=subset(counts, select=c(-width)),
              MARGIN = 2,
              FUN = function(x){
                10^9*x/geneLengths/sum(as.numeric(x))
              } )

colSums(rpkm)


##computing TPM

rpk <- apply(subset(counts, select=c(-width)), 2, 
             function(x)x/(geneLengths/1000))

tpm <- apply(rpk, 2, function(x)x/sum(as.numeric(x))*10^6)


colSums(tpm)

#8.36 Exploratory analysis of the read count table

V <- apply(tpm,1,var)

selectedGenes <- names(V[order(V, decreasing=T)][1:100])
selectedGenes
library(pheatmap)
pheatmap(tpm[selectedGenes, ], scale="row", show_rownames = F)


colData <- read.table(coldata_file, header=T, sep = "\t",
                      stringsAsFactors = T)
head(colData)
View(colData)
pheatmap(tpm[selectedGenes, ], scale="row", show_rownames = F,
         annotation_col = colData) #tutoraril heapmap
##PCA----
library(stats)
library(ggplot2)
library(ggbio)
M <- t(tpm[selectedGenes,])
?t()


M <- log2(M+1)


pcaResults <- prcomp(M)
pcaResults
autoplot( pcaResults, aes(data=colData, colour="group"))#don't works
summary(pcaResults)

#correlation plots


correlationMatrix <- cor(tpm)
View(correlationMatrix)
library(corrplot)
corrplot(correlationMatrix, order="hclust",
         addrect=2, addCoef.col = "pink",
         number.cex = 0.7)

pheatmap(correlationMatrix,
         annotation_col=colData,
         cutree_cols = 2)
###8.37 Differential Expression analysis

countData <- as.matrix(subset(counts, select=c(-width)))
colData <- read.table(coldata_file, header = T, sep = "\t",
                      stringsAsFactors = T)
colData
View(countData)

designFormula <- " ~group"
designFormula

library(DESeq2)
library(stats)

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = as.formula(designFormula))
print(dds)
                    

rownames(dds)
colnames(dds)
counts(dds)
colData(dds)

dds <- dds[rowSums(DESeq2::counts(dds))>1,]          
dds <- DESeq(dds)
print(dds)
DEresults=results(dds, contrast = c("group", "CASE", "CTRL"))
DEresults <- DEresults[order(DEresults$pvalue),]
DEresults
#Diagnostic plots
library(DESeq2)
DESeq2::plotMA(object=dds, ylim=c(-5,5))


library(ggplot2)
  ggplot(data=as.data.frame(DEresults), aes(x=pvalue))+
    geom_histogram(bins = 100)

countsNormalized <- DESeq2::counts(dds, normalized=T)
selectedGenes <- names(sort(apply(countsNormalized, 1, var),
                            decreasing = T)[1:500])

selectedGenes

###plotPCA(countsNormalized[selectedGenes,],
   #     col=as.numeric(colData$group), adj=0.5,
    #    xlim=c(-0.5,0.5), ylim=C(-0.5, 0.6)) # Don´t works


rld <- rlog(dds)

DESeq2::plotPCA(rld, ntop=500, intgroup="group")+
  ylim(-50,50)+ theme_bw()

library(EDASeq)
par(mfrow=c(1,2))
plotRLE(countData, outline=F, ylim=c(-4,4),
        col=as.numeric(colData$group), main="Raw Counts")

plotRLE(DESeq2::counts(dds, normalized=T), 
        outline=F, ylim=c(-4,4),
        col=as.numeric(colData$group), 
        main="Normalized Counts")

##Functional enrichment analysis----
install.packages("gprofiler2")
library(gprofiler2)
library(knitr)
library(DESeq2)
DEresults <- results(dds, contrast=c("group", "CASE", "CTRL"))
DE <- DEresults[!is.na(DEresults$padj),]
DE
DE <- DE[DE$padj<0.1,]
DE <- DE[abs(DE$log2FoldChange)>1,]
DE

genesofinterest <- rownames(DE)
genesofinterest

goResults <- gost(query = genesofinterest,
                       organism="hsapiens",
                  sources = "GO",
                  
                       )
names(goResults)
head(goResults$result)


View(goResults$meta)
goResults <- goResults$result
goResults
#Gene set enrichment analysis GSEA

goResults <- goResults[order(goResults$p.value),]

go <- goResults[goResults$intersection_size<100,]
go
colnames(go)

geneSet1 <- unlist(strsplit(go[1,]$intersection_size, ","))
geneSet1 <- unlist(strsplit(go[1,]$intersection, ','))

####Accounting for additional sources of variation

counts_file <- system.file("extdata/rna-seq/SRP021193.raw_counts.tsv",
                           package = "compGenomRData")



colData_file <- system.file("extdata/rna-seq/SRP021193.colData.tsv",
                            package = "compGenomRData")





counts <- read.table(counts_file)
colData <- read.table(colData_file, header=T, sep="\t",
                      stringsAsFactors = T)
View(colData)

geneLengths <- counts$width
View(geneLengths)

rpk <- apply(subset(counts, select=c(-width)),2,
             function(x)x/geneLengths/1000)
rpk

tmp <- apply(rpk, 2, function(x)x/sum(as.numeric(x))*10^6)

slectedGenes <- names(sort(apply(tpm, 1, var),
                           decreasing = T) [1:100])
tmp

pheatmap::pheatmap(tpm[slectedGenes,],
         scale="row",
         annotation_col = colData,
         show_rownames = F)

countData <- as.matrix(subset(counts, select=c(-width)))

countData

dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData = colData,
                              design= ~ LibrarySelection + group)
dds


dds <- DESeq(dds)
DEresults <- results(dds, contrast = c("group", "CASE", "CTRL"))
DEresults
###Accounting for estimated covariates using RUVSeq-----




counts_file <- system.file("extdata/rna-seq/SRP049988.raw_counts.tsv",
                           package = "compGenomRData")

colData_file <- system.file("extdata/rna-seq/SRP049988.colData.tsv",
                           package = "compGenomRData")

counts <- read.table(counts_file)
colData <- read.table(colData_file, header = T,
                      sep="\t", stringsAsFactors = T)
View(colData)

colData$source_name <- ifelse(colData$group=="CASE",
                              "EHF_OVEREXPRESSION", "EMPTY VECTOR")


rpk <- apply(subset(counts, select=c(-width)),2,
             function(x)x/geneLengths/1000)
rpk

tmp <- apply(rpk, 2, function(x)x/sum(as.numeric(x))*10^6)

slectedGenes <- names(sort(apply(tpm, 1, var),
                           decreasing = T) [1:100])
tmp

pheatmap::pheatmap(tpm[slectedGenes,],
                   scale="row",
                   annotation_col = colData,
                   show_rownames = F)

library(EDASeq)

countData <- as.matrix(subset(counts, select=c(-width)))
set <- newSeqExpressionSet(counts=countData,
                           phenoData = colData)



par(mfrow=c(1,2))


plotRLE(set, outline=F, ylim=c(-4,4), col=as.numeric(colData$group))
plotPCA(set, col=as.numeric(colData$group), adj=0.5,
        ylim=c(-0.7,0.5), xlim=c(-0.5,0.5))



par(mfrow=c(1,2))


plotRLE(tpm, outline=F, ylim=c(-4,4), col=as.numeric(colData$group))
plotPCA(tpm, col=as.numeric(colData$group), adj=0.5,
        ylim=c(-0.7,0.5), xlim=c(-0.5,0.5))





