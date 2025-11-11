b=1e3  #bin
G= 3e9 #genome size
N=G/b

f=10000/N



data_path=system.file("extdata/chip-seq", package = "compGenomRData")
chip_files=list.files(data_path, full.names=T)
chip_files
library(GenomeInfoDb)
hg_chrs=getChromInfoFromUCSC("hg38")
hg_chrs
hg_chrs=subset(hg_chrs, grepl("chr21$", chrom))
hg_chrs

seqlengths=with(hg_chrs, setNames(size, chrom))
seqlengths
help(with)


library(GenomicRanges)


tilling_window=tileGenome(seqlengths, tilewidth = 1000)


tilling_window=unlist(tilling_window)

tilling_window

library(GenomicAlignments)
bam_files=list.files(
  path=data_path,
  full.names = T,
  pattern = "bam$"
)
bam_files

so=summarizeOverlaps(tilling_window, bam_files)
so

counts=assays(so)[[1]]
counts

cpm=t(t(counts)*(10e6/colSums(counts)))
cpm=cpm[rowSums(cpm)>0,]
cpm  




colnames(cpm)=sub(".chr21.bam", "", colnames(cpm))
colnames(cpm)
colnames(cpm)=sub("GM12878_hg38_", "", colnames(cpm))
colnames(cpm)

correlation_matrix=cor(cpm, method = "pearson")
correlation_matrix

library(ComplexHeatmap)
library(circlize)

heatmap_col=circlize::colorRamp2(
  breaks = c(-1,0,1),
  colors=c("blue", "white","red")
)


Heatmap(
  matrix = correlation_matrix,
  col=heatmap_col
)



###9.3 Visualization in the genome browser


bam_files=list.files(
  path=data_path,
  full.names=T,
  pattern = "bam$"
  
)

chip_file=bam_files[1]
chip_file
reads=readGAlignments(chip_file)
reads
reads=granges(reads)
reads



reads=resize(reads, width = 200, fix="start")


reads=keepSeqlevels(reads, "chr21", pruning.mode = "coarse")
reads

cov=coverage(reads, width = seqlengths)
cov


output_file=sub(".bam", ".bigWig", chip_file)

rtracklayer::export.bw(cov, "output_file")

library(Gviz)
axis=GenomeAxisTrack(
  range=GRanges("chr21", IRanges(1, width = seqlengths))
)

gcov=as(cov, "GRanges")
dtrack=DataTrack(gcov, name="CTCF", type="1")

track_list=list(axis, dtrack)


plotTracks(
  trackList = track_list,
  sizes=c(.1,1),
  background.title="black", type=c("a","p")
)
ws2

###Jacard similarity

chip_file

reads=readGAlignments(chip_file)
reads=granges(reads)
reads=resize(reads, width=1, fix="start")
reads=keepSeqlevels(reads, "chr21", pruning.mode = "coarse")
reads



reads=split(reads, strand(reads))
reads


cov=lapply(reads, function(x){
  coverage(x, width=seqlengths)[[1]]>0
})
cov


cov=lapply(cov, as.vector)
cov


wsize=1:400
jaccard=function(x,y)sum((x&y))/sum((x|y))

cc=shiftApply(
  SHIFT = wsize,
  X= cov[["+"]],
  Y=cov[["-"]],
  FUN=jaccard)
cc

cc=data.frame(fragment_size=wsize, cross_correlations=cc)

cc

library(ggplot2)
ggplot(data = cc, aes(fragment_size, cross_correlations)) +
  geom_point() +
  geom_vline(xintercept = which.max(cc$cross_correlation), 
             size=2, color='red', linetype=2) +
  theme_bw() +
  theme(
    axis.text = element_text(size=10, face='bold'),
    axis.title = element_text(size=14,face="bold"),
    plot.title = element_text(hjust = 0.5)) +
  xlab('Shift in base pairs') +
  ylab('Jaccard similarity') 

###9.5.5 GC quantification----




library(GenomeInfoDb)
library(GenomicRanges)


hg_chrs=getChromInfoFromUCSC("hg38")
hg_chrs=subset(hg_chrs, grepl("chr21$", chrom))
hg_chrs

seqlengths=with(hg_chrs, setNames(size, chrom))
seqlengths

tilling_window=unlist(tileGenome(
  seqlengths = seqlengths,
  tilewidth = 1000
))
tilling_window


library(BSgenome.Hsapiens.UCSC.hg38)
set=getSeq(BSgenome.Hsapiens.UCSC.hg38, tilling_window)

set

nuc=oligonucleotideFrequency(set, width = 2)
nuc= as.data.frame(nuc)
nuc

nuc=round(nuc/1000,3)
View(nuc)


so=summarizeOverlaps(tilling_window, bam_files)
so

counts= assays(so)[[1]]
counts

cpm=t(t(counts))*(1000000/colSums(counts))
cpm


cpm_log=log10(cpm+1)
cpm_log

gc=cbind(data.frame(cpm_log), GC=nuc["GC"])
gc

ggplot(
  data = gc, 
  aes(
    x = GC, 
    y = GM12878_hg38_CTCF_r1.chr21.bam
  )) +
  geom_point(size=2, alpha=.3) +
  theme_bw() +
  theme(
    axis.text  = element_text(size=10, face='bold'),
    axis.title = element_text(size=14,face="bold"),
    plot.title = element_text(hjust = 0.5)) +
  xlab('GC content in one kilobase windows') +
  ylab('log10( cpm + 1 )') +
  ggtitle('CTCF Replicate 1')












library(tidyr)
gcd=pivot_longer(
  data=gc,
  cols=-GC,
  names_to="experiment", 
  values_to="cpm"
)

gcd=subset(gcd, grepl("CTCF", experiment))
gcd

gcd$experiment=sub("chr21.","", gcd$experiment)
gcd

ggplot(data = gcd, aes(GC, log10(cpm+1))) +
  geom_point(size=2, alpha=.05) +
  theme_bw() +
  facet_wrap(~experiment, nrow=1)+
  theme(
    axis.text = element_text(size=10, face='bold'),
    axis.title = element_text(size=14,face="bold"),
    plot.title = element_text(hjust = 0.5)) +
  xlab('GC content in one kilobase windows') +
  ylab('log10( cpm + 1 )') +
  ggtitle('CTCF Replicates 1 and 2')

library(AnnotationHub)
hub=AnnotationHub()

AnnotationHub::query(
  x=hub,
  pattern= c("ENSEMBL", "Homo", "GRCh38", "chr", "gtf")
)

gtf <- hub[["AH61126"]]

ensembl_seqlevels=seqlevels(gtf)

ucsc_seqlevels=paste0("chr", ensembl_seqlevels)

seqlevels(gtf, pruning.mode="coarse")=ucsc_seqlevels

gtf=gtf[seqnames(gtf)=="chr21"]

seqlevels(gtf)


annotation_list=GRangesList(
  tss=promoters(
    x=subset(gtf, type="gene"),
    upstream=1000,
    downstream=1000),
  exon=subset(gtf, type="exon"),
  intron=subset(gtf, type=="gene")
  
)
annotation_list


annotatReads=function(bam_file, annotation_list){
  library(dplyr)
  message(basename(bam_file))
  
  bam=readGAlignments(bam_file)
  result=as.data.frame(
    findOverlaps(bam, annotation_list)
    
    
  )
  annotation_name=names(annotation_list)[result$subjectHits]
  result$annotation=annotation_name
  result= result[order(result$subjectHits),]
  result=subset(result, !duplicated(queryHits))
  
  result=group_by(.data=result, annotation)
  
  result=summarise(.data = result, counts=length(annotation))
  
  result=rbind(
    result,
    
    data.frame(
      annotation="intergenic",
      counts= length(bam)-sum(result$counts)
    )
  )
  
  result$frequency=with(result, round(counts/sum(counts), 2))
  result$experiment=basename(bam_file)
  return(result)
  
  
  
}
bam_files=list.files(data_path, full.names = T, pattern = "bam$")
bam_files


annot_reads_list <- lapply(bam_files, function(x){
  annotatReads(
    bam_file = x,
    annotation_list = annotation_list
  )
})


annotation_list

annot_reads_df=dplyr::bind_rows(annot_reads_list)
annot_reads_df$experiment

View(annot_reads_df)

experiment_name=annot_reads_df$experiment

experiment_name=sub(".chr21.bam", "", experiment_name)
experiment_name=sub("GM12878_hg38_", "", experiment_name)

annot_reads_df$experiment=experiment_name
library(ggplot2)

ggplot(data = annot_reads_df, 
       aes(
         x    = experiment, 
         y    = frequency, 
         fill = annotation
       )) +
  geom_bar(stat='identity') +
  theme_bw() +
  scale_fill_brewer(palette='Set2') +
  theme(
    axis.text = element_text(size=10, face='bold'),
    axis.title = element_text(size=14,face="bold"),
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab('Sample') +
  ylab('Percentage of reads') +
  ggtitle('Percentage of reads in annotation')



#Peak Calling----
chip_files = list(
  H3K4me3  = 'GM12878_hg38_H3K4me3.chr21.bw',
  
  H3K36me3 = 'GM12878_hg38_H3K36me3.chr21.bw',
  
  POL2     = 'GM12878_hg38_POLR2A.chr21.bw'
)


chip_files=lapply(chip_files, function(x){
  file.path(data_path, x)
}
  )
chip_files

library(rtracklayer)
chip_profiles=lapply(chip_files, rtracklayer::import.bw)
chip_profiles


library(AnnotationHub)
hub=AnnotationHub()
gtf=hub[["AH61126"]]

seqlevels(gtf, pruning.mode="coarse")="21"
ensembl_seqlevels=seqlevels(gtf)

ensembl_seqlevels
ucsc_seqlevels= paste0("chr", ensembl_seqlevels)
ucsc_seqlevels

seqlevels(gtf, pruning.mode="coarse")=ucsc_seqlevels
seqlevels(gtf)


library(GenomicFeatures)
library(Gviz)




txdb=makeTxDbFromGRanges(gtf)



gene_track=Gviz::GeneRegionTrack(txdb, chr="chr21", genome="hg38")
gene_track


hg_chrs=getChromInfoFromUCSC("hg38")
hg_chrs=subset(hg_chrs, (grepl("chr21$", chrom)))

seqlengths=with(hg_chrs, setNames(size, chrom))
seqlengths

chr_track=IdeogramTrack(
  chromosome="chr21",
  genome="hg38"
)


axis=GenomeAxisTrack(
  range=GRanges("chr21", IRanges(1, width = seqlengths))
)


data_tack=lapply(names(chip_profiles), function(experiment_name){
  DataTrack(
    range=chip_profiles[[experiment_name]],
    name=experiment_name,
    type="h",
    lwd=5
  )
  
})



data_tack




start <- min(start(subset(gtf, gene_name="URB1")))
end <- max(end(subset(gtf, gene_name="URB1")))


plotTracks(
  trackList = c(chr_track, axis, gene_track, data_tack),
  sizes = c(2,2,2,2,2,2),
  background.title="black",
  collapseTranscripts="longest",
  transcriptAnnotation="symbol",
  
  from=start,
  to=end 
  
)


