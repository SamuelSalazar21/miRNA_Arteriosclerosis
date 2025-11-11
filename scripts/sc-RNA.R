BiocManager::install("TENxPBMCData")
BiocManager::install("SingleCellExperiment")
library(TENxPBMCData)
library(SingleCellExperiment)

sce=TENxPBMCData("pbmc5k-CITEseq")
sce

colData(sce)
rowData(sce)
assays(sce)


rownames(sce)=paste0(rowData(sce)$ENSEMBL_ID, "_",
                     rowData(sce)$Symbol_TENx)
colnames(sce)=sce$Sequence
counts_method1=assays(sce)[["counts"]]
counts_method1[1:5,1:3]


counts_method2=assay(sce, "counts")
counts_method2[1:5, 1:3]



counts_method3=counts(sce)
counts_method3[1:5, 1:3]



frac_0_per_cell=colMeans(counts(sce)==0)

frac_0_per_cell
mean(counts(sce)==0)
head(colSums(counts(sce)))

summary(colSums(counts(sce)))

class(counts(sce))
counts(sce)=as(counts(sce), "dgCMatrix")
class(counts(sce))
counts(sce)[1:10, 1:10]

plot(rowSums(counts(sce)), 1-rowMeans(counts(sce) == 0), 
     log = "x",xlab='total count per gene',
     ylab='fraction of cell in which the gene is detected')
