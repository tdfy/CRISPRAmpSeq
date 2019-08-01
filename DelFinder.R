setwd("c:/Export/Fang/PIGA/map/")
library(ggplot2)
library(GenomicFeatures)
library(ggbio)
library(rtracklayer)
library(GenomicRanges)
library(magick)
library(dplyr)
library(Rsamtools)
library(GenomicAlignments)
library(VariantAnnotation)
library(reshape)


filenames <- list.files("c:/Export/20190705_Amp_SeqFiltered-138247127/", pattern=glob2rx("190705_B2M_*_CIGAR_summary.csv"), full.names=TRUE) # 190705_B2M_2_CIGAR_summary.csv

BAM <- 'ccs.unknown.5--5.sorted.25.bam.csv'


CIGAR_output <- read.csv(paste(BAM), header=TRUE, sep=",")

names(CIGAR_output)[4]<-"Hits"

indel<- ggplot(CIGAR_output, aes(x=Var_len, y= Hits, color=Variant,fill=Variant)) + geom_bar(stat='identity')+theme_bw()+ scale_y_continuous(label=function(Hit){abs(Hit)})


penta_hope = paste(BAM, ".png", sep="")
ggsave(penta_hope)

