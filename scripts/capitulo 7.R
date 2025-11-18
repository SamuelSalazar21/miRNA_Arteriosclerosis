####Script Short Reads



BiocManager::install("Rqc")
library(Rqc)

dir(system.file("extdata"
                package="ShortRead"))

folder=system.file(package="ShortRead", "extdata/E-MTAB-1147")
folder
qcRes=rqc(path=folder, pattern = ".fastq.gz", openBrowser = T)

rqcCycleQualityBoxPlot(qcRes)
rqcCycleBaseCallsLinePlot(qcRes)
rqcReadFrequencyPlot(qcRes)
 BiocManager::install("fastqcr")
library("fastqcr") 
fastqc_install() 
fastqc(fq.dir = folder, qc.dir = "fastqc_results")

#Filtering and trimming reads----
BiocManager::install("QuasR")
library("QuasR")
fasqfile <- system.file(package="ShortRead",
                           "extdata/E-MTAB-1147",
                           c("ERR127302_1_subset.fastq.gz",
                             "ERR127302_2_subset.fastq.gz")
)
                        
outfiles <- paste(tempfile(pattern=c("processed_1_", "processed_2_")), ".fastq",
                  sep="") 
outfiles                        


preprocessReads(fasqfile, outfiles,
                nBases=1,
                truncateEndBases = 3,
                Lpattern = "ACCCGGGA",
                minLength = 40)


fasqfile <- system.file(package="ShortRead",
                        "extdata/E-MTAB-1147",
                        "ERR127302_1_subset.fastq.gz")
fasqfile

fq=readFastq(fasqfile)
fq
qPerBase=as(quality(fq), "matrix")
qPerBase

quality(fq)

qcount <- rowSums(qPerBase<=20)
qcount
fq[qcount==0]


writeFastq(fq[qcount==0],
           paste(fasqfile, "Qfiltered", sep="-"))
f <- FastqStreamer(fasqfile, readerBlockSize = 1000)
while(length(fq <- yield(f))){
  qPerBase=as(quality(fq), "matrix")
  qcount=rowSums(qPerBase<=20)
  writeFastq(fq[qcount==0],
             paste(fasqfile,"Qfiltered",sep="-"),
             mode="a")
}


#mapping aligment
library(QuasR)
file.copy(system.file(package ="QuasR", "extdata"),".", recursive=T)
genomeFile <- "extdata/hg19sub.fa"
sampleFile <- "extdata/samples_chip_single.txt"
proj <- qAlign(sampleFile, genomeFile)
proj



###Exercises




folder=(system.file(package = "QuasR", "extdata"))
folder
dir(folder)
qcCes=rqc(path=folder, pattern = "^chip", openBrowser = T)
qcCes


rqcCycleQualityBoxPlot(qcCes)#calidad de cada base según su posición
rqcCycleBaseCallsLinePlot(qcCes) #porcentaje CG
rqcReadFrequencyPlot(qcCes)# proporcion de read repetidos





fastqFiles <- system.file(package="QuasR",
                          "extdata",
                          c("chip_1_1.fq.bz2","chip_2_1.fq.bz2"))
fastqFiles


outfiles <- paste(tempfile(pattern=c("processed_1_",
                                     "processed_2_")),".fastq",sep="")
outfiles




#
proced<- preprocessReads(fastqFiles, outfiles,
                         truncateEndBases = 4)
proced

fqp2=readFastq("C:/Users/PC/Documents/extdata/processed_2_45b82def26c8.fastq")
fqp2
qperB2=as(qualsity(fqp2), "matrix")
qcount2=rowSums(qperB2<20)
writeFastq(fqp2[qcount2==0],paste("Qfiltered", "processed_2_45b82def26c8.fastq", sep="-"))

fqp1 <- readFastq("C:/Users/PC/Documents/extdata/processed_1_45b84c213522.fastq")

qperB1=as(quality(fqp1), "matrix")
qperB1
qcount1=rowSums(qperB1<=20)
fqp1[qcount1==0]
writeFastq(fqp1[qcount1==0], paste("Qfiltered","processed_1_45b84c213522.fastq" ,sep="-"))


##3punto7
file.copy(system.file(package ="QuasR", "extdata"),".", recursive=T)


sampleFile2 <- "extdata/samples_ejercio_7.txt"
sampleFile2
fin <- qAlign(sampleFile2, genomeFile)

fin2 <- qAlign(sampleFile2, genomeFile)
fin
fin2






