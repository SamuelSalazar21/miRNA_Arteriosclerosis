if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Gviz")
install.packages("rlang")
BiocManager::install("ggbio")
packageVersion("rlang")
remove.packages("rlang")

library(Gviz)
afrom <- 2960000
ato <- 3160000
alTrack <- AlignmentsTrack(system.file(package="Gviz","extdata", "gapped.bam"),
                          isPaired=TRUE)
alTrack
plotTracks(alTrack, from=afrom, to=ato, chromosome = "chr12")
plotTracks(alTrack, from=afrom, to=ato, chromosome = "chr12", type="coverage")
plotTracks(alTrack, from=afrom, to=ato, chromosome = "chr12", type="pileup")

library(ggbio)
p.ideo <- Ideogram(genome="hg19")
p.ideo + xlim(GRanges("chr2", IRanges(1e8, 1e8+1000000)))
BiocManager::install("Homo.sapiens")
library(Homo.sapiens)
data(genesymbol, package = "biovizBase")
wh <- genesymbol[c("BRCA1", "NBR1")]
wh <- range(wh, ignore.strand=TRUE)
p.txdb <- autoplot(Homo.sapiens, which= wh)
p.txdb

autoplot(Homo.sapiens, which=wh, label.color="black", color="brown", 
         fill="white")

autoplot(Homo.sapiens, which=wh, gap.geom="segment")
autoplot(Homo.sapiens, which=wh, stat="reduce")




###TxDb object #graficar genes
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
autoplot(txdb, which=wh)
###EnsDb object
BiocManager::install("EnsDb.Hsapiens.v75")
library("EnsDb.Hsapiens.v75")
ensdb <- "EnsDb.Hsapiens.v75"
ensdb
autoplot(ensdb, GeneNameFilter("PHKG2"))
autoplot(ensdb, ~symbol=="PHKG2", names.eprx="gene_name")
gr <- GRanges(seqnames = 16, IRanges(30768000,30770000), strand="+")
gr
autoplot(ensdb, GRangesFilter(gr), names.expr="gene_name")





##Make gene model from GRangesList Object
BiocManager::install("biovizBase")
library(biovizBase)


gr.txdb <- crunch(txdb, which=wh)
gr.txdb
values(gr.txdb)
colnames(values(gr.txdb))[4] <- "model"
grl <- split(gr.txdb, gr.txdb$tx_id)
grl
names(grl) <- sample(LETTERS, size=length(grl),replace = T)
grl

autoplot(grl, aes(type=model))
ggplot()+geom_alignment(grl, type="model")
##Add a reference track

BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
library(BSgenome.Hsapiens.UCSC.hg19)
bg <- BSgenome.Hsapiens.UCSC.hg19
p.bg <- autoplot(bg, which=wh)
p.bg
p.bg+zoom(1/100)
p.bg+zoom(1/1000)
p.bg+zoom(1/2500)
p.bg <- autoplot(bg, which=resize(wh, width=width(wh)/2000), geom="segment")

p.bg


fl.bam <- system.file("extdata", "wg-brca1.sorted.bam",package = "biovizBase")

