setwd("c:/Export/20190508_Amplicon_Seq-126488167/")
library(ggplot2)
library(GenomicFeatures)
library(ggbio)
library(rtracklayer)
library(GenomicRanges)
library(magick)


samp_name <- c("R1C1 S1D5","R1C1 S1D7","R1C1 S1D9","R1C1 Post-H","D(-4)_1","D(-4)_2")


filenames <- list.files("c:/Export/20190508_Amplicon_Seq-126488167/", pattern=glob2rx("*pivot_summary.csv"), full.names=TRUE) # A_B2M_2_CIGAR_summary.csv

for(target in filenames) #<----- 1st & 5th samples, retains B2Mk labels and Y axis
{

  # link <- substr(target,43,45)
  link <-lapply(strsplit(target, "_"), `[`, 3)
  link2 <-lapply(strsplit(toString(link), "/"), `[`, 2)
  
  
  target <- read.csv(target, header=TRUE, sep=",")
  # colnames(target) <- c("Variant_Type","Sample_Name","Percentage")
  target$Percentage=round(target$Percentage, 1)
  target$Sample_Name <- sapply(strsplit(as.character(target$Sample_Name), ""), tail, 1)
  target$Sample_Name <-as.numeric(target$Sample_Name)
  target$Sample_Name <-samp_name[target$Sample_Name]
  target$Sample_Name <- vapply(target$Sample_Name, paste, collapse = ", ", character(1L))

  target_comp <-ggplot(target, aes(x = Sample_Name, y = Percentage, fill=Variant, label=Percentage)) +
    geom_bar(colour='black',stat = "identity")+theme_bw()+geom_text(aes(label=ifelse(abs(Percentage)>1,as.character(Percentage),'')),position = position_stack(vjust = 0.5))+
    scale_fill_manual(values=c(D="red",DX="red3", I="green2", IX="green4", NoEdit = "cornflowerblue",X="blue",DIX="purple"),breaks=c("D", "DX", "I","IX","NoEdit","X","DIX"),labels=c("Deletion", "Deletion/\nSubstitution", "Insertion","Insertion/\nSubstitution","No Edit","Substitution","Insertion/\nDeleltion/Sub."))+
    theme(axis.text.x  = element_text(size=7)) + xlab(NULL) + guides(fill = guide_legend(title = "Edit Type")) + theme(legend.position="right")+scale_y_continuous(expand = c(0, 0))+ scale_x_discrete(limits = samp_name)
  

  hope = paste(link2, ".png", sep="")
  ggsave(hope)
  
}
