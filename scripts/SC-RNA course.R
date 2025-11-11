##singlecellcourse.org
BiocManager::install("DropletUtils")
tung_counts <- read.table("tung/molecules.txt", sep="\t")

tung_annotation <- read.table("tung/annotation.txt", sep="\t", header=T)

View(tung_annotation)

tung <- SingleCellExperiment(
  assays=list(counts=as.matrix(tung_counts)),
  colData=tung_annotation
)

rm(tung_annotation, tung_counts)

tung

rownames(tung)

View(assay(tung, "counts"))


library(DropletUtils)
sce <- read10xCounts("tung/pbmc_1k_raw")
sce <- read10xCounts("tung/pbmc_1k_filtered")


dir("tung")

tung
class(colData(tung))
class(assay(tung))

colData(tung)
table(colData(tung)$batch)
assay(tung, "logcounts") <- log2(counts(tung)+1)
View(assay(tung))


logcounts(tung)[1:10, 1:4]

colMeans(counts(tung))

mean(counts(tung)[,1])

colData(tung)$mean_counts <- colMeans(counts(tung))

colData(tung)$total_counts <- colSums(counts(tung))

colData(tung)
colData(tung)

assay(tung, "cpm") <- sweep(counts(tung), 2, tung$total_counts/1e6, "/")
colSums(cpm(tung))[1:10]
cpm(tung)









