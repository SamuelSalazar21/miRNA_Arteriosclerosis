#fasqc
library(Rqc)
folder=system.file(package = "ShortRead", "extdata/E-MTAB-1147")

qcRes=rqc(path=folder, pattern = ".fastq.gz", openBrowser = FALSE)
qcRes=rqc(path=folder, pattern = ".fastq.gz", openBrowser = T)

rqcCycleQualityBoxPlot(qcRes)
rqcCycleBaseCallsLinePlot(qcRes)
rqcReadFrequencyPlot(qcRes)
rqcReadFrequencyCalc(qcRes)


###fasqcr----
library(fastqcr)
BiocManager::install("fastqcr")
fastqc_install()

qc.dir <- system.file("fastqc_results", package = "fastqcr")
qc <- qc_aggregate(qc.dir)
qc
qc_fails(qc, "module")
qc_fails(qc,"sample")
qc_report(qc.dir, result.file = "multi-qc-report")


###Building one-sample QC Reports (+ Interpretation)

qc.file <- system.file("fastqc_results", "S1_fastqc.zip", package = "fastqcr")
qc_report(qc.file, result.file = "One-Sample-report", interpret = T)

##filtering and trimming reads
library(QuasR)

Fastqfiles <-system.file(package = "ShortRead",
                         "extdata/E-MTAB-1147",
                         c("ERR127302_1_subset.fastq.gz",
                           "ERR127302_2_subset.fastq.gz")) 
outfiles <- paste(tempfile(pattern = c("procesado_1_", "procesado_2_")),
                  ".fastq", sep = " ")
outfiles

preprocessReads(Fastqfiles, outfiles,
                nBases = 1,
                truncateEndBases = 3,
                Lpattern = "ACCCGGGA",
                minLength = 40)
library(ShortRead)
fastqFile <- system.file(package="ShortRead",
                         "extdata/E-MTAB-1147",
                         "ERR127302_1_subset.fastq.gz")
fq=readFastq(fastqFile)
fq
qPerBase=as(quality(fq), "matrix")
qPerBase
qcount <- rowSums(qPerBase <=20)
fq[qcount==0]

writeFastq(fq[qcount==0],
           paste(fastqFile,"Qfiltered", sep="_"))
f <- FastqStreamer(fastqFile, readerBlockSize = 1000)
while(length(fq <- yield(f))){
  qPerBase=as(quality(fq), "matrix")
  qcount=rowSums(qPerBase <=20)
  writeFastq(fq[qcount==0],
             paste(fastqFile, "Qfiltered", sep="_"), mode="a")
}


#Mapping/aligning reads to genome
library(QuasR)
file.copy(system.file(package = "QuasR","extdata"), ".", recursive = T)
genomeFile <- "extdata/hg19sub.fa"
sampleFile <- "extdata/samples_chip_single.txt"
proj <- qAlign(sampleFile,genomeFile)
proj

