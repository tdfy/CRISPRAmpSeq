setwd("c:/Export/20190508_Amplicon_Seq-126488167/")
library(ggplot2)
library(GenomicFeatures)
library(ggbio)
library(rtracklayer)
library(GenomicRanges)
library(magick)

samp_name <- list("VP-19-002.01.G1/TDN,R1C1 S1D5","VP-19-002.01.G1/TDN,R1C1 S1D7","VP-19-002.01.G1/TDN,R1C1 S1D9","VP-19-002.01.G1/TDN,R1C1 Post-H","VP-19-002.01.G1 D(-4)_1","VP-19-002.01.G1 D(-4)_2")

B2M <- image_read("B2M.png")
PIK <- image_read("PIK.png")
CIITA <- image_read("CIITA.png")
TRAC <- image_read("TRAC3.png")

#_______________________________________________________________________________________________________________#

avsINT1.granges <- GRanges(seqnames="15",
                           ranges=IRanges(start=44711473,
                                          end=44711546),
                           strand="*")

avsEXON.granges <- GRanges(seqnames="15",
                           ranges=IRanges(start=44711547,
                                          end=44711613),
                           strand="*")


avsINT2.granges <- GRanges(seqnames="15",
                           ranges=IRanges(start=44711614,
                                          end=44711726),
                           strand="*")

GUIDEstart <- 44711473 + 94
GUIDEend <- 44711473 + 115


avsB2Mpam.granges <- GRanges(seqnames="15",
                             ranges=IRanges(start=GUIDEstart,
                                            end=GUIDEend),
                             strand="*")

locus.granges <- GRanges(seqnames="15",
                         ranges=IRanges(start=44711473,
                                        end=44711726),
                         strand="*")

CIGAR = c(0) 
Hit = c(0) 
Position = start(avsINT1.granges):end(avsINT2.granges)
Variant = c("X")
Sample_Name = c("blank")
Dummy_df = data.frame(CIGAR, Hit, Position, Variant, Sample_Name)



filenames <- list.files("c:/Export/20190508_Amplicon_Seq-126488167/", pattern=glob2rx("190508_B2M_*_CIGAR_summary.csv"), full.names=TRUE) # 190508_B2M_2_CIGAR_summary.csv

for(sample in filenames[c(1,5)]) #<----- 1st & 5th samples, retains B2Mk labels and Y axis
{
  mydataB2M <- read.csv(sample, header=TRUE, sep=",")
  mydataB2M$Position= as.numeric(mydataB2M$Position)+44711473
  
  
  if (nrow(mydataB2M)==0)
  {
    mydataB2M <- Dummy_df 
  }
  
  
  link <- substr(sample,50,54)
  link_no <-substr(sample,54,54)
  har <- as.integer((link_no))
  new_title <- toString(samp_name[har])
  
  
  B2M_cigar <- ggplot(mydataB2M,aes(y = Hit, x = Position, fill = Variant)) + geom_bar(width=1.5,alpha=.8, stat="identity")+theme_bw()+scale_y_continuous(limit=c(0,50000),expand = c(.01, .01))+scale_x_continuous(expand = c(0, 0))+ylab(NULL)+
    scale_fill_manual(values=c(D="red", I="green4", X="blue"),breaks=c("D", "I", "X"),labels=c("Deletion", "Insertion", "Substitution"))+ 
    theme(legend.position= c(0.8, 0.85),axis.text.x  = element_text(size=12))+guides(fill = guide_legend(title = "Edit Type"))
  
  B2M_loc <- autoplot(avsEXON.granges, geom_rect(),aes(fill="black",color="black"))+scale_fill_identity() +scale_color_identity()+theme_bw()+
    geom_chevron(avsINT1.granges,color="black", fill=NA)+geom_chevron(avsINT2.granges,color="black", fill=NA)+geom_rect(avsB2Mpam.granges,color="black", fill="white")+
    scale_x_continuous(expand = c(0, 0)) + theme(axis.text.x  = element_text(size=8))
  
  sample <- tracks("Reads"=B2M_cigar,"B2M"=B2M_loc,heights=c(0.1,0.01),label.text.cex=c(1,1),title = new_title,label.bg.fill=('white'))
  
  hope = paste(link, ".png", sep="")
  ggsave(hope)
}

for(sample in filenames[c(4,6)]) #<----- 4th & 7th samples,  Y axis (right)
{
  mydataB2M <- read.csv(sample, header=TRUE, sep=",")
  mydataB2M$Position= as.numeric(mydataB2M$Position)+44711473
  
  if (nrow(mydataB2M)==0)
  {
    mydataB2M <- Dummy_df 
  }
  
  
  link <- substr(sample,50,54)
  link_no <-substr(sample,54,54)
  har <- as.integer((link_no))
  new_title <- toString(samp_name[har])
  
  
  B2M_cigar <- ggplot(mydataB2M,aes(y = Hit, x = Position, fill = Variant)) + geom_bar(width=1.5,alpha=.8, stat="identity")+theme_bw()+scale_y_continuous(limit=c(0,50000),position="right",expand = c(.01, .01))+
    scale_x_continuous(expand = c(0, 0))+ylab(NULL)+scale_fill_manual(values=c(D="red", I="green4", X="blue"),breaks=c("D", "I", "X"),labels=c("Deletion", "Insertion", "Substitution"))+
    theme(legend.position= "none",axis.text.x  = element_text(size=12))
  
  
  sample <- tracks("Reads"=B2M_cigar,"B2M"=B2M_loc,heights=c(0.1,0.01),label.text.cex=c(1,1),title = new_title,label.bg.fill=('white'))
  
  hope = paste(link, ".png", sep="")
  ggsave(hope)
}


for(sample in filenames[c(2,3)]) #<-- All other samples do not retain legend, y-axis, B2Mk labels (cropped later)
{
  
  mydataB2M <- read.csv(sample, header=TRUE, sep=",")
  mydataB2M$Position= as.numeric(mydataB2M$Position)+44711473
  
  if (nrow(mydataB2M)==0)
  {
    mydataB2M <- Dummy_df 
  }
  
  
  link <- substr(sample,50,54)
  link_no <-substr(sample,54,54)
  har <- as.integer((link_no))
  new_title <- toString(samp_name[har])
  
  B2M_cigar <- ggplot(mydataB2M,aes(y = Hit, x = Position, fill = Variant)) + geom_bar(width=1.5,alpha=.8, stat="identity")+theme_bw()+scale_y_continuous(limit=c(0,50000),expand = c(.01, .01))+scale_x_continuous(expand = c(0, 0))+ylab(NULL)+
    scale_fill_manual(values=c(D="red", I="green4", X="blue"),breaks=c("D", "I", "X"),labels=c("Deletion", "Insertion", "Substitution"))+ 
    theme(legend.position="none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.text.x  = element_text(size=12))
  
  
  sample <- tracks("Reads"=B2M_cigar,"B2M"=B2M_loc,heights=c(0.1,0.01),label.text.cex=c(1,1),title = new_title,label.bg.fill=('white'))
  
  
  hope = paste(link, ".png", sep="")
  ggsave(hope)
}

B2M_1 <- image_read("B2M_1.png")
B2M_2 <- image_read("B2M_2.png")
B2M_3 <- image_read("B2M_3.png")
B2M_4 <- image_read("B2M_4.png")
B2M_5 <- image_read("B2M_5.png")
B2M_6 <- image_read("B2M_6.png")
# B2M_7 <- image_read("B2M_7.png")

B2M_2 <- image_crop(B2M_2,"-150x1703+0+0")
B2M_3 <- image_crop(B2M_3,"-150x1703+0+0")
B2M_4 <- image_crop(B2M_4,"-150x1703+0+0")
# B2M_5 <- image_crop(B2M_5,"-150x1703+0+0")
B2M_6 <- image_crop(B2M_6,"-150x1703+0+0")
# B2M_7 <- image_crop(B2M_7,"-150x1703+0+0")

tiger_style <- image_append(c(B2M_1, B2M_2,B2M_3,B2M_4))
woo <- image_append(c(B2M_5, B2M_6,image_scale(B2M, "1700x1900")))
ghost <- image_append(c(tiger_style,woo),stack=TRUE)
indeed11<-image_annotate(ghost, "A", font = 'serif', size = 70,location = "+175+1560")
indeed12<-image_annotate(indeed11, "B", font = 'serif', size = 70,location = "+175+3260")
indeed13<-image_annotate(indeed12, "C", font = 'serif', size = 70,location = "+3100+3290")

image_write(indeed13,"B2M_VP-19-002.01.G1.png")



###____________________CIITA_____________________________________________________________###

avsINT1.granges <- GRanges(seqnames="16",
                           ranges=IRanges(start=10916241,
                                          end=10916366),
                           strand="*")

avsEXON.granges <- GRanges(seqnames="16",
                           ranges=IRanges(start=10916367,
                                          end=	10916459),
                           strand="*")


avsINT2.granges <- GRanges(seqnames="16",
                           ranges=IRanges(start=10916460,
                                          end=10916471),
                           strand="*")

GUIDEstart <- 10916241 + 145
GUIDEend <- 10916241 + 168 # + 3


avsCIITApam.granges <- GRanges(seqnames="16",
                               ranges=IRanges(start=GUIDEstart,
                                              end=GUIDEend),
                               strand="*")

filenames <- list.files("c:/Export/20190508_Amplicon_Seq-126488167/", pattern=glob2rx("190508_CIITA_*_CIGAR_summary.csv"), full.names=TRUE) # 190508_CIITA_2_CIGAR_summary.csv

for(sample in filenames[c(1,5)]) #<----- 1st & 5th samples, retains CIITAk labels and Y axis
{
  mydataCIITA <- read.csv(sample, header=TRUE, sep=",")
  
  mydataCIITA$Position= as.numeric(mydataCIITA$Position)+10916241
  
  if (nrow(mydataCIITA)==0)
  {
    mydataCIITA <- Dummy_df 
  }
  
  CIITA_subset <-subset(mydataCIITA, Position!=10916241+53)
  
  
  link <- substr(sample,50,56)
  link_no <-substr(sample,56,56)
  har <- as.integer((link_no))
  new_title <- toString(samp_name[har])
  
  
  CIITA_cigar <- ggplot(CIITA_subset,aes(y = Hit, x = Position, fill = Variant)) + geom_bar(width=1.5,alpha=.8, stat="identity")+theme_bw()+scale_y_continuous(limit=c(0,50000),expand = c(.01, .01))+scale_x_continuous(expand = c(0, 0))+ylab(NULL)+
    scale_fill_manual(values=c(D="red", I="green4", X="blue"),breaks=c("D", "I", "X"),labels=c("Deletion", "Insertion", "Substitution"))+ 
    theme(legend.position= c(0.2, 0.85),text = element_text(size=18))+guides(fill = guide_legend(title = "Edit Type"))
  
  
  CIITA_loc <- autoplot(avsEXON.granges, geom_rect(),aes(fill="black",color="black"))+scale_fill_identity() +scale_color_identity()+theme_bw()+
    geom_chevron(avsINT1.granges,color="black", fill=NA)+geom_chevron(avsINT2.granges,color="black", fill=NA)+geom_rect(avsCIITApam.granges,color="black", fill="white")+
    scale_x_continuous(expand = c(0, 0)) + theme(text = element_text(size=18))
  
  sample <- tracks("Reads"=CIITA_cigar,"CIITA"=CIITA_loc,heights=c(0.1,0.01),label.text.cex=c(1,1),title = new_title,label.bg.fill=('white'))
  
  hope = paste(link, ".png", sep="")
  ggsave(hope)
}

for(sample in filenames[c(4,7)]) #<----- 4th & 7th samples,  Y axis (right)
{
  mydataCIITA <- read.csv(sample, header=TRUE, sep=",")
  
  mydataCIITA$Position= as.numeric(mydataCIITA$Position)+10916241
  
  if (nrow(mydataCIITA)==0)
  {
    mydataCIITA <- Dummy_df 
  }
  
  
  CIITA_subset <-subset(mydataCIITA, Position!=10916241+53)
  
  link <- substr(sample,50,56)
  link_no <-substr(sample,56,56)
  har <- as.integer((link_no))
  new_title <- toString(samp_name[har])
  
  
  CIITA_cigar <- ggplot(CIITA_subset,aes(y = Hit, x = Position, fill = Variant)) + geom_bar(width=1.5,alpha=.8, stat="identity")+theme_bw()+scale_y_continuous(limit=c(0,50000),position="right",expand = c(.01, .01))+
    scale_x_continuous(expand = c(0, 0))+ylab(NULL)+scale_fill_manual(values=c(D="red", I="green4", X="blue"),breaks=c("D", "I", "X"),labels=c("Deletion", "Insertion", "Substitution"))+
    theme(legend.position= "none",text = element_text(size=18))
  
  
  sample <- tracks("Reads"=CIITA_cigar,"CIITA"=CIITA_loc,heights=c(0.1,0.01),label.text.cex=c(1,1),title = new_title,label.bg.fill=('white'))
  
  hope = paste(link, ".png", sep="")
  ggsave(hope)
}


for(sample in filenames[c(2,3,6)]) #<-- All other samples do not retain legend, y-axis, CIITAk labels (cropped later)
{
  
  mydataCIITA <- read.csv(sample, header=TRUE, sep=",")
  
  mydataCIITA$Position= as.numeric(mydataCIITA$Position)+10916241
  
  if (nrow(mydataCIITA)==0)
  {
    mydataCIITA <- Dummy_df 
  }
  
  # CIITA_subset <-subset(mydataCIITA, Position!=10916241+53)
  CIITA_subset <-mydataCIITA
  
  link <- substr(sample,50,56)
  link_no <-substr(sample,56,56)
  har <- as.integer((link_no))
  new_title <- toString(samp_name[har])
  
  CIITA_cigar <- ggplot(CIITA_subset,aes(y = Hit, x = Position, fill = Variant)) + geom_bar(width=1.5,alpha=.8, stat="identity")+theme_bw()+scale_y_continuous(limit=c(0,50000),expand = c(.01, .01))+scale_x_continuous(expand = c(0, 0))+ylab(NULL)+
    scale_fill_manual(values=c(D="red", I="green4", X="blue"),breaks=c("D", "I", "X"),labels=c("Deletion", "Insertion", "Substitution"))+ 
    theme(legend.position="none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),text = element_text(size=18))
  
  
  sample <- tracks("Reads"=CIITA_cigar,"CIITA"=CIITA_loc,heights=c(0.1,0.01),label.text.cex=c(1,1),title = new_title,label.bg.fill=('white'))
  
  
  hope = paste(link, ".png", sep="")
  ggsave(hope)
}

CIITA_1 <- image_read("CIITA_1.png")
CIITA_2 <- image_read("CIITA_2.png")
CIITA_3 <- image_read("CIITA_3.png")
CIITA_4 <- image_read("CIITA_4.png")
CIITA_5 <- image_read("CIITA_5.png")
CIITA_6 <- image_read("CIITA_6.png")
# CIITA_7 <- image_read("CIITA_7.png")

CIITA_2 <- image_crop(CIITA_2,"-150x1703+0+0")
CIITA_3 <- image_crop(CIITA_3,"-150x1703+0+0")
CIITA_4 <- image_crop(CIITA_4,"-150x1703+0+0")
# CIITA_5 <- image_crop(CIITA_5,"-150x1703+0+0")
CIITA_6 <- image_crop(CIITA_6,"-150x1703+0+0")
# CIITA_7 <- image_crop(CIITA_7,"-150x1703+0+0")

tiger_style <- image_append(c(CIITA_1, CIITA_2,CIITA_3,CIITA_4))
woo <- image_append(c(CIITA_5, CIITA_6,image_scale(CIITA, "1700x1900")))
ghost <- image_append(c(tiger_style,woo),stack=TRUE)
indeed11<-image_annotate(ghost, "A", font = 'serif', size = 70,location = "+175+1560")
indeed12<-image_annotate(indeed11, "B", font = 'serif', size = 70,location = "+175+3260")
indeed13<-image_annotate(indeed12, "C", font = 'serif', size = 70,location = "+3100+3290")

image_write(indeed13,"CIITA_VP-19-002.01.G1.png")



###____________________TRAC_____________________________________________________________###
avsINT1.granges <- GRanges(seqnames="14",
                           ranges=IRanges(start=22550515,
                                          end=22550556),
                           strand="*")

avsEXON.granges <- GRanges(seqnames="14",
                           ranges=IRanges(start=22550557,
                                          end=	22550664),
                           strand="*")


avsINT2.granges <- GRanges(seqnames="14",
                           ranges=IRanges(start=22550664,
                                          end=22550749),
                           strand="*")

GUIDEstart <- 22550515 + 108
GUIDEend <- 22550515 + 130



avsTRACpam.granges <- GRanges(seqnames="14",
                              ranges=IRanges(start=GUIDEstart,
                                             end=GUIDEend),
                              strand="*")


filenames <- list.files("c:/Export/20190508_Amplicon_Seq-126488167/", pattern=glob2rx("190508_TRAC3_*_CIGAR_summary.csv"), full.names=TRUE) # 190508_TRAC_2_CIGAR_summary.csv

for(sample in filenames[c(1,5)]) #<----- 1st & 5th samples, retains track labels and Y axis
{
  mydataTRAC <- read.csv(sample, header=TRUE, sep=",")
  
  mydataTRAC$Position= as.numeric(mydataTRAC$Position)+22550515
  
  if (nrow(mydataTRAC)==0)
  {
    mydataTRAC <- Dummy_df 
  }
  
  TRAC_subset <-subset(mydataTRAC, Position!=22550515+189)
  
  link <- substr(sample,50,56)
  link_no <-substr(sample,56,56)
  har <- as.integer((link_no))
  new_title <- toString(samp_name[har])
  
  
  TRAC_cigar <- ggplot(TRAC_subset,aes(y = Hit, x = Position, fill = Variant)) + geom_bar(width=1.5,alpha=.8, stat="identity")+theme_bw()+scale_y_continuous(limit=c(0,50000),expand = c(.01, .01))+scale_x_continuous(expand = c(0, 0))+ylab(NULL)+
    scale_fill_manual(values=c(D="red", I="green4", X="blue"),breaks=c("D", "I", "X"),labels=c("Deletion", "Insertion", "Substitution"))+ 
    theme(legend.position= c(0.8, 0.85),text = element_text(size=18))+guides(fill = guide_legend(title = "Edit Type"))
  
  
  TRAC_loc <- autoplot(avsEXON.granges, geom_rect(),aes(fill="black",color="black"))+scale_fill_identity() +scale_color_identity()+theme_bw()+
    geom_chevron(avsINT1.granges,color="black", fill=NA)+geom_chevron(avsINT2.granges,color="black", fill=NA)+geom_rect(avsTRACpam.granges,color="black",fill="white")+
    scale_x_continuous(expand = c(0, 0)) + theme(text = element_text(size=18))
  
  sample <- tracks("Reads"=TRAC_cigar,"TRAC"=TRAC_loc,heights=c(0.1,0.01),label.text.cex=c(1,1),title = new_title,label.bg.fill=('white'))
  
  
  print(link)
  hope = paste(link, ".png", sep="")
  ggsave(hope)
}

for(sample in filenames[c(4,6)]) #<----- 4th & 7th samples,  Y axis (right)
{
  mydataTRAC <- read.csv(sample, header=TRUE, sep=",")
  
  mydataTRAC$Position= as.numeric(mydataTRAC$Position)+22550515
  
  if (nrow(mydataTRAC)==0)
  {
    mydataTRAC <- Dummy_df 
  }
  
  TRAC_subset <-subset(mydataTRAC, Position!=22550515+189)
  
  link <- substr(sample,50,56)
  link_no <-substr(sample,56,56)
  har <- as.integer((link_no))
  new_title <- toString(samp_name[har])
  
  
  TRAC_cigar <- ggplot(TRAC_subset,aes(y = Hit, x = Position, fill = Variant)) + geom_bar(width=1.5,alpha=.8, stat="identity")+theme_bw()+scale_y_continuous(limit=c(0,50000),position="right",expand = c(.01, .01))+
    scale_x_continuous(expand = c(0, 0))+ylab(NULL)+scale_fill_manual(values=c(D="red", I="green4", X="blue"),breaks=c("D", "I", "X"),labels=c("Deletion", "Insertion", "Substitution"))+
    theme(legend.position= "none",text = element_text(size=18))
  
  sample <- tracks("Reads"=TRAC_cigar,"TRAC"=TRAC_loc,heights=c(0.1,0.01),label.text.cex=c(1,1),title = new_title,label.bg.fill=('white'))
  
  hope = paste(link, ".png", sep="")
  ggsave(hope)
}


for(sample in filenames[c(2,3)]) #<-- All other samples do not retain legend, y-axis, track labels (cropped later)
{
  
  mydataTRAC <- read.csv(sample, header=TRUE, sep=",")
  
  mydataTRAC$Position= as.numeric(mydataTRAC$Position)+22550515
  
  if (nrow(mydataTRAC)==0)
  {
    mydataTRAC <- Dummy_df 
  }
  
  
  TRAC_subset <-subset(mydataTRAC, Position!=22550515+189)
  
  link <- substr(sample,50,56)
  link_no <-substr(sample,56,56)
  har <- as.integer((link_no))
  new_title <- toString(samp_name[har])
  
  TRAC_cigar <- ggplot(TRAC_subset,aes(y = Hit, x = Position, fill = Variant)) + geom_bar(width=1.5,alpha=.8, stat="identity")+theme_bw()+scale_y_continuous(limit=c(0,50000),expand = c(.01, .01))+scale_x_continuous(expand = c(0, 0))+ylab(NULL)+
    scale_fill_manual(values=c(D="red", I="green4", X="blue"),breaks=c("D", "I", "X"),labels=c("Deletion", "Insertion", "Substitution"))+ 
    theme(legend.position="none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),text = element_text(size=18))
  
  sample <- tracks("Reads"=TRAC_cigar,"TRAC"=TRAC_loc,heights=c(0.1,0.01),label.text.cex=c(1,1),title = new_title,label.bg.fill=('white'))
  
  print(link)
  hope = paste(link, ".png", sep="")
  ggsave(hope)
}

TRAC3_1 <- image_read("TRAC3_1.png")
TRAC3_2 <- image_read("TRAC3_2.png")
TRAC3_3 <- image_read("TRAC3_3.png")
TRAC3_4 <- image_read("TRAC3_4.png")
TRAC3_5 <- image_read("TRAC3_5.png")
TRAC3_6 <- image_read("TRAC3_6.png")
# TRAC_7 <- image_read("TRAC3_7.png")

TRAC3_2 <- image_crop(TRAC3_2,"-150x1703+0+0")
TRAC3_3 <- image_crop(TRAC3_3,"-150x1703+0+0")
TRAC3_4 <- image_crop(TRAC3_4,"-150x1703+0+0")
# TRAC_5 <- image_crop(TRAC_5,"-150x1703+0+0")
TRAC3_6 <- image_crop(TRAC3_6,"-150x1703+0+0")
# TRAC_7 <- image_crop(TRAC3_7,"-150x1703+0+0")

tiger_style <- image_append(c(TRAC3_1, TRAC3_2,TRAC3_3,TRAC3_4))
woo <- image_append(c(TRAC3_5, TRAC3_6,image_scale(TRAC, "1700x1900")))
ghost <- image_append(c(tiger_style,woo),stack=TRUE)
indeed11<-image_annotate(ghost, "A", font = 'serif', size = 70,location = "+175+1560")
indeed12<-image_annotate(indeed11, "B", font = 'serif', size = 70,location = "+175+3260")
indeed13<-image_annotate(indeed12, "C", font = 'serif', size = 70,location = "+3100+3290")


image_write(indeed13,"TRAC3_VP-19-002.01.G1.png")

#__________PIK__________________________________________________________________________#

PIK.granges <- GRanges(seqnames="17",
                       ranges=IRanges(start=8901263,
                                      end=8901484),
                       strand="*")

PIK.sgRNA <- GRanges(seqnames="17",
                       ranges=IRanges(start=8901382,
                                      end=8901402),
                       strand="*")

CIGAR = c(0) 
Hit = c(0) 
Position = start(PIK.granges):end(PIK.granges)
Variant = c("X")
Sample_Name = c("blank")
PIK_Dummy_df = data.frame(CIGAR, Hit, Position, Variant, Sample_Name)


filenames <- list.files("c:/Export/20190508_Amplicon_Seq-126488167/", pattern=glob2rx("190508_PIK_*_CIGAR_summary.csv"), full.names=TRUE) # 190508_PIK_2_CIGAR_summary.csv

for(sample in filenames[c(1,5)]) #<----- 1st & 5th samples, retains PIKk labels and Y axis
{
  mydataPIK <- read.csv(sample, header=TRUE, sep=",")
  
  mydataPIK$Position= as.numeric(mydataPIK$Position)+8901263
  
  if (nrow(mydataPIK)==0)
  {
    mydataPIK <- PIK_Dummy_df 
  }
  
  link <- substr(sample,50,54)
  link_no <-substr(sample,54,54)
  har <- as.integer((link_no))
  new_title <- toString(samp_name[har])
  
  
  PIK_cigar <- ggplot(mydataPIK,aes(y = Hit, x = Position, fill = Variant)) + geom_bar(width=1.5,alpha=.8, stat="identity")+theme_bw()+scale_y_continuous(limit=c(0,50000),expand = c(.01, .01))+scale_x_continuous(expand = c(0, 0))+ylab(NULL)+
    scale_fill_manual(values=c(D="red", I="green4", X="blue"),breaks=c("D", "I", "X"),labels=c("Deletion", "Insertion", "Substitution"))+ 
    theme(legend.position= c(0.8, 0.85),text = element_text(size=12))+guides(fill = guide_legend(title = "Edit Type"))
  
  
  PIK_loc <- autoplot(PIK.granges, geom = "chevron",offset=0.3,aes(fill="black",color="black"))+scale_fill_identity() +scale_color_identity()+theme_bw()+
    geom_rect(PIK.sgRNA,color="black",fill=NA)+
    scale_x_continuous(expand = c(0, 0))+theme(text = element_text(size=12))
  
  sample <- tracks("Reads"=PIK_cigar,"PIK3"=PIK_loc,heights=c(0.1,0.01),label.text.cex=c(1,1),title = new_title,label.bg.fill=('white'))
  
  hope = paste(link, ".png", sep="")
  ggsave(hope)
}


for(sample in filenames[c(4,7)]) #<----- 4th & 7th samples,  Y axis (right)
{
  mydataPIK <- read.csv(sample, header=TRUE, sep=",")
  
  mydataPIK$Position= as.numeric(mydataPIK$Position)+8901263
  
  if (nrow(mydataPIK)==0)
  {
    mydataPIK <- PIK_Dummy_df 
  }
  
  link <- substr(sample,50,54)
  link_no <-substr(sample,54,54)
  har <- as.integer((link_no))
  new_title <- toString(samp_name[har])
  
  
  PIK_cigar <- ggplot(mydataPIK,aes(y = Hit, x = Position, fill = Variant)) + geom_bar(width=1.5,alpha=.8, stat="identity")+theme_bw()+scale_y_continuous(limit=c(0,50000),position="right",expand = c(.01, .01))+
    scale_x_continuous(expand = c(0, 0))+ylab(NULL)+scale_fill_manual(values=c(D="red", I="green4", X="blue"),breaks=c("D", "I", "X"),labels=c("Deletion", "Insertion", "Substitution"))+
    theme(legend.position= "none",text = element_text(size=12))
  
  sample <- tracks("Reads"=PIK_cigar,"PIK3"=PIK_loc,heights=c(0.1,0.01),label.text.cex=c(1,1),title = new_title,label.bg.fill=('white'))
  
  hope = paste(link, ".png", sep="")
  ggsave(hope)
}


for(sample in filenames[c(2,3,6)]) #<-- All other samples do not retain legend, y-axis, PIKk labels (cropped later)
{
  
  mydataPIK <- read.csv(sample, header=TRUE, sep=",")
  
  mydataPIK$Position= as.numeric(mydataPIK$Position)+8901263
  
  if (nrow(mydataPIK)==0)
  {
    mydataPIK <- PIK_Dummy_df 
  }
  
  link <- substr(sample,50,54)
  link_no <-substr(sample,54,54)
  har <- as.integer((link_no))
  new_title <- toString(samp_name[har])
  
  PIK_cigar <- ggplot(mydataPIK,aes(y = Hit, x = Position, fill = Variant)) + geom_bar(width=1.5,alpha=.8, stat="identity")+theme_bw()+scale_y_continuous(limit=c(0,50000),expand = c(.01, .01))+scale_x_continuous(expand = c(0, 0))+ylab(NULL)+
    scale_fill_manual(values=c(D="red", I="green4", X="blue"),breaks=c("D", "I", "X"),labels=c("Deletion", "Insertion", "Substitution"))+ 
    theme(legend.position="none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),text = element_text(size=12))
  
  
  sample <- tracks("Reads"=PIK_cigar,"PIK3"=PIK_loc,heights=c(0.1,0.01),label.text.cex=c(1,1),title = new_title,label.bg.fill=('white'))
  
  
  hope = paste(link, ".png", sep="")
  ggsave(hope)
}

PIK_1 <- image_read("PIK_1.png")
PIK_2 <- image_read("PIK_2.png")
PIK_3 <- image_read("PIK_3.png")
PIK_4 <- image_read("PIK_4.png")
PIK_5 <- image_read("PIK_5.png")
PIK_6 <- image_read("PIK_6.png")
# PIK_7 <- image_read("PIK_7.png")

PIK_2 <- image_crop(PIK_2,"-150x1703+0+0")
PIK_3 <- image_crop(PIK_3,"-150x1703+0+0")
PIK_4 <- image_crop(PIK_4,"-150x1703+0+0")
# PIK_5 <- image_crop(PIK_5,"-150x1703+0+0")
PIK_6 <- image_crop(PIK_6,"-150x1703+0+0")
# PIK_7 <- image_crop(PIK_7,"-150x1703+0+0")



tiger_style <- image_append(c(PIK_1, PIK_2,PIK_3,PIK_4))
woo <- image_append(c(PIK_5, PIK_6,image_scale(PIK, "1700x1900")))
ghost <- image_append(c(tiger_style,woo),stack=TRUE)

indeed11<-image_annotate(ghost, "A", font = 'serif', size = 70,location = "+175+1560")
indeed12<-image_annotate(indeed11, "B", font = 'serif', size = 70,location = "+175+3260")
indeed13<-image_annotate(indeed12, "C", font = 'serif', size = 70,location = "+3100+3290")

image_write(indeed13,"PIK_VP-19-002.01.G1.png")
