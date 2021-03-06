setwd("c:/Export/20191202_307_360/Seq_files/")

library(ggplot2)
library(GenomicFeatures)
library(ggbio)
library(rtracklayer)
library(GenomicRanges)
library(magick)
library(dplyr)

samp_name <- list()

B2M <- image_read("B2M_307.png") %>% image_border("grey", "2x1") 
# B2M <-image_annotate(B2M, "D", font = 'serif', size = 70,location = "+20+20")
PIK <- image_read("PIK_307.png") %>% image_border("grey", "2x1")
# PIK <-image_annotate(PIK, "D", font = 'serif', size = 70,location = "+20+20")
CIITA <- image_read("CIITA_307.png") %>% image_border("grey", "2x1")
# CIITA <-image_annotate(CIITA, "D", font = 'serif', size = 70,location = "+20+20")
TRAC <- image_read("TRAC_307.png") %>% image_border("grey", "2x1")
# TRAC <-image_annotate(TRAC, "D", font = 'serif', size = 70,location = "+20+20")

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
Variant = c("D")
Sample_Name = c("blank")
Var_len = 1
plot_pos = start(avsINT1.granges):end(avsINT2.granges)
Percentage = c(0)
Dummy_1 = data.frame(CIGAR, Hit, Position, Variant, Sample_Name, Var_len,plot_pos,Percentage)

Variant = c("I")
Dummy_2 = data.frame(CIGAR, Hit, Position, Variant, Sample_Name,Var_len,plot_pos,Percentage)
Dummy_df = rbind(Dummy_2,Dummy_1)


Var_len = c(0)
Hit = c(0)
Variant = c('I')
Dummy_ins_agg = data.frame(Var_len,Hit,Variant)

Var_len = c(0)
Hit = c(0)
Variant = c('D')
Dummy_del_agg = data.frame(Var_len,Hit,Variant)

filenames <- list.files("c:/Export/20191202_307_360/Seq_files/", pattern=glob2rx("20191209_B2M_*_CIGAR_summary.csv"), full.names=TRUE) # 190705_B2M_2_CIGAR_summary.csv

max_list <- list()


for (sample in filenames[c(1:4)])
{
  df<-read.csv(sample, header=TRUE, sep=",")
  melt <- aggregate(df$Percentage, by=list(Position=df$Position), FUN=sum)
  upper <- max(melt$x)
  max_list <- append(upper,max_list)
}

max_max <- max(unlist(max_list))

for (sample in filenames[c(1:4)])
{
  mydataB2M <- read.csv(sample, header=TRUE, sep=",")
  mydataB2M$Position= as.numeric(mydataB2M$Position)+44711473
  mydataB2M$plot_pos= mydataB2M$Position+(mydataB2M$Var_len/2)
  
  if (nrow(mydataB2M)==0)
  {
    mydataB2M <- Dummy_df
    del <- Dummy_df
    ins <- Dummy_df
  }
  
  
  del <- subset(mydataB2M, Variant=='D')
  ins <- subset(mydataB2M, Variant=='I')

  ### Penta DF Block #####
  
  alt1 <-mydataB2M

  alt <- alt1[,c("Variant","Hit","Var_len")]
  
  deletion <- subset(alt, Variant=='D')
  insertion <- subset(alt, Variant=='I')
  
  if (nrow(deletion)!=0){
    agg_del <- aggregate(Hit ~ Var_len, data=deletion, FUN=sum)
    agg_del$Hit = 0 - as.numeric(agg_del$Hit)
    agg_del$Variant <- "D"
  }else {agg_del <-Dummy_del_agg}
  
  if (nrow(insertion)!=0){
    agg_ins <- aggregate(Hit ~ Var_len, data=insertion, FUN=sum)
    agg_ins$Variant <- "I"
  }else {agg_ins <-Dummy_ins_agg }
  
  
  new_data = rbind(agg_del,agg_ins)
  

  ###_______________________________________________________#########
  
  mydataB2M <-mydataB2M[
    with(mydataB2M, order(-mydataB2M$plot_pos,-mydataB2M$Var_len)),
    ]
  
  
  del <-del[
    with(del, order(-del$plot_pos,-del$Var_len)),
    ]
  
  ins <-ins[
    with(ins, order(-ins$plot_pos,-ins$Var_len)),
    ]
  

  link <- substr(sample,47,49)
  link_no <-substr(sample,51,51)
  har <- as.integer((link_no))
  new_title <- toString(samp_name[har])
  
  B2M_cigar <- ggplot(mydataB2M,aes(y = Percentage, x = plot_pos, fill = Variant, group = Variant)) + geom_bar(data=del,width=del$Var_len,alpha=.5, stat="identity")+
    geom_bar(data=ins,width=ins$Var_len,alpha=.5, stat="identity")+
    theme_bw() + scale_x_continuous(expand = c(0, 0))+ylab(NULL) +scale_y_continuous(limits=c(0,max_max),expand = c(.01, .01))+
    scale_fill_manual(values=c(D="red", I="green4", X="blue"),breaks=c("D", "I", "X"),labels=c("Deletion", "Insertion", "Substitution"))+
    theme(legend.position= c(0.8, 0.85),axis.text.x  = element_text(size=12))+guides(fill = guide_legend(title = "Edit Type"))
  
  

    B2M_loc <- autoplot(avsEXON.granges, geom_rect(),aes(fill="black",color="black"))+scale_fill_identity() +scale_color_identity()+theme_bw()+
    geom_chevron(avsINT1.granges,color="black", fill=NA)+geom_chevron(avsINT2.granges,color="black", fill=NA)+geom_rect(avsB2Mpam.granges,color="black", fill="white")+
    scale_x_continuous(expand = c(0, 0)) + theme(axis.text.x  = element_text(size=6))
  
  sample <- tracks("Percentage"=B2M_cigar,"B2M"=B2M_loc,heights=c(0.1,0.01),label.text.cex=c(1,1),title = new_title,label.bg.fill=('white'))
  
  hope = paste(link,link_no, ".png", sep="")
  ggsave(hope)
  
  if (max(abs(new_data$Hit)) > 1000){
    
    penta <- ggplot(new_data, aes(x = Var_len, y = Hit, fill=Variant)) + geom_bar(stat="identity", position="identity")+
      theme_bw() + scale_fill_manual(values=c(D="red", I="green4"),labels=c("Deletion", "Insertion"))+theme(legend.position="none") + 
      ylab('Reads')+ xlab('Indel Length (bp)') + scale_y_continuous(label=function(Hit){abs(Hit)})+ 
      ggtitle("decoy") + theme(plot.title = element_text(color="white", size=15))
    
  }else {
    penta <- ggplot(new_data, aes(x = Var_len, y = Hit, fill=Variant)) + geom_bar(stat="identity", position="identity")+
      theme_bw() + scale_fill_manual(values=c(D="red", I="green4"),labels=c("Deletion", "Insertion"))+theme(legend.position="none") + 
      ylab('Reads')+ xlab('Indel Length (bp)') + scale_y_continuous(limits = c(-100000,100000))+ 
      ggtitle("decoy") + theme(plot.title = element_text(color="white", size=15))}
  
  # print(new_data)
  
  penta_hope = paste(link,link_no, "_penta.png", sep="")
  ggsave(penta_hope, width = 2, height = 5)
  samp_plot <- image_read(paste(link,link_no,".png",sep=""))
  pent_plot <- image_read(paste(link,link_no,"_penta.png",sep=""))

  pen_cat <- image_append(c(samp_plot,pent_plot))
  pen_cat <- image_border(pen_cat, "grey", "2x1")
  pen_cat <-image_annotate(pen_cat, "A", font = 'serif', size = 70,location = "+120+100")
  pen_cat <-image_annotate(pen_cat, "B", font = 'serif', size = 70,location = "+120+1550")
  pen_cat <-image_annotate(pen_cat, "C", font = 'serif', size = 70,location = "+1180+100")
  

  image_write(pen_cat,(paste(link,link_no,"_Cat.png",sep="")))
  
}

B2M_1 <- image_read("B2M1_Cat.png")
B2M_2 <- image_read("B2M2_Cat.png")
B2M_3 <- image_read("B2M3_Cat.png")
B2M_4 <- image_read("B2M4_Cat.png")
# B2M_5 <- image_read("B2M_5_Cat.png")


tiger_style <- image_append(c(B2M_1, B2M_2))
method_style <- image_append(c(B2M_3, B2M_4))
# woo <- image_append(c(B2M_4,B2M))
woo <- B2M
ghost <- image_append(c(tiger_style,method_style,woo),stack=TRUE)


image_write(ghost,(paste(new_title,"_B2M.png",sep="")))

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

CIGAR = c(0)
Hit = c(0)
Position = start(avsINT1.granges):end(avsINT2.granges)
Variant = c("D")
Sample_Name = c("blank")
Var_len = 1
plot_pos = start(avsINT1.granges):end(avsINT2.granges)
Percentage = c(0)
Dummy_1 = data.frame(CIGAR, Hit, Position, Variant, Sample_Name, Var_len,plot_pos,Percentage)

Variant = c("I")
Dummy_2 = data.frame(CIGAR, Hit, Position, Variant, Sample_Name,Var_len,plot_pos,Percentage)
Dummy_df = rbind(Dummy_2,Dummy_1)


filenames <- list.files("c:/Export/20191202_307_360/Seq_files/", pattern=glob2rx("20191209_CIITA_*_CIGAR_summary.csv"), full.names=TRUE) # 190705_B2M_2_CIGAR_summary.csv


max_list <- list()


for (sample in filenames[c(1:4)])
{
  df<-read.csv(sample, header=TRUE, sep=",")
  melt <- aggregate(df$Percentage, by=list(Position=df$Position), FUN=sum)
  upper <- max(melt$x)
  max_list <- append(upper,max_list)
}

max_max <- max(unlist(max_list))

for (sample in filenames[c(1:4)])
{
  mydataCIITA <- read.csv(sample, header=TRUE, sep=",")
  mydataCIITA$Position= as.numeric(mydataCIITA$Position)+10916241 ###<------------------ CHECK
  mydataCIITA$plot_pos= mydataCIITA$Position+(mydataCIITA$Var_len/2)
  
  if (nrow(mydataCIITA)==0)
  {
    mydataCIITA <- Dummy_df
    del <- Dummy_df
    ins <- Dummy_df
  }
  
  del <- subset(mydataCIITA, Variant=='D')
  ins <- subset(mydataCIITA, Variant=='I')
  
  
  ### Penta DF Block #####
  
  alt1 <-mydataCIITA
  
  alt <- alt1[,c("Variant","Hit","Var_len")]
  
  deletion <- subset(alt, Variant=='D')
  insertion <- subset(alt, Variant=='I')
  
  if (nrow(deletion)!=0){
    agg_del <- aggregate(Hit ~ Var_len, data=deletion, FUN=sum)
    agg_del$Hit = 0 - as.numeric(agg_del$Hit)
    agg_del$Variant <- "D"
  }else {agg_del <-Dummy_del_agg}
  
  if (nrow(insertion)!=0){
    agg_ins <- aggregate(Hit ~ Var_len, data=insertion, FUN=sum)
    agg_ins$Variant <- "I"
  }else {agg_ins <-Dummy_ins_agg }
  
  
  new_data = rbind(agg_del,agg_ins)
  ##_________________________________________________________##
  
  mydataCIITA <-mydataCIITA[
    with(mydataCIITA, order(-mydataCIITA$plot_pos,-mydataCIITA$Var_len)),
    ]
  
  
  del <-del[
    with(del, order(-del$plot_pos,-del$Var_len)),
    ]
  
  ins <-ins[
    with(ins, order(-ins$plot_pos,-ins$Var_len)),
    ]
  
  
    link <- substr(sample,47,51)
    link_no <-substr(sample,53,53)
    har <- as.integer((link_no))
    new_title <- toString(samp_name[har])
  
  CIITA_cigar <- ggplot(mydataCIITA,aes(y = Percentage, x = plot_pos, fill = Variant, group = Variant)) + geom_bar(data=del,width=del$Var_len,alpha=.5, stat="identity")+
    geom_bar(data=ins,width=ins$Var_len,alpha=.5, stat="identity")+
    theme_bw() + scale_x_continuous(expand = c(0, 0))+ylab(NULL) +scale_y_continuous(limits=c(0,max_max),expand = c(.01, .01))+
    scale_fill_manual(values=c(D="red", I="green4", X="blue"),breaks=c("D", "I", "X"),labels=c("Deletion", "Insertion", "Substitution"))+
    theme(legend.position= c(0.8, 0.85),axis.text.x  = element_text(size=12))+guides(fill = guide_legend(title = "Edit Type"))
  
  
  
  CIITA_loc <- autoplot(avsEXON.granges, geom_rect(),aes(fill="black",color="black"))+scale_fill_identity() +scale_color_identity()+theme_bw()+
    geom_chevron(avsINT1.granges,color="black", fill=NA)+geom_chevron(avsINT2.granges,color="black", fill=NA)+geom_rect(avsCIITApam.granges,color="black", fill="white")+
    scale_x_continuous(expand = c(0, 0)) + theme(axis.text.x  = element_text(size=6))
  
  sample <- tracks("Percentage"=CIITA_cigar,"CIITA"=CIITA_loc,heights=c(0.1,0.01),label.text.cex=c(1,1),title = new_title,label.bg.fill=('white'))
  
  hope = paste(link,link_no, ".png", sep="")
  ggsave(hope)
  
  if (max(abs(new_data$Hit)) > 1000){
    
    penta <- ggplot(new_data, aes(x = Var_len, y = Hit, fill=Variant)) + geom_bar(stat="identity", position="identity")+
      theme_bw() + scale_fill_manual(values=c(D="red", I="green4"),labels=c("Deletion", "Insertion"))+theme(legend.position="none") + 
      ylab('Reads')+ xlab('Indel Length (bp)') + scale_y_continuous(label=function(Hit){abs(Hit)})+ 
      ggtitle("decoy") + theme(plot.title = element_text(color="white", size=15))
    
  }else {
    penta <- ggplot(new_data, aes(x = Var_len, y = Hit, fill=Variant)) + geom_bar(stat="identity", position="identity")+
      theme_bw() + scale_fill_manual(values=c(D="red", I="green4"),labels=c("Deletion", "Insertion"))+theme(legend.position="none") + 
      ylab('Reads')+ xlab('Indel Length (bp)') + scale_y_continuous(limits = c(-100000,100000))+ 
      ggtitle("decoy") + theme(plot.title = element_text(color="white", size=15))}
  
  
  # quad <- ggplot(new_data, aes(x = Var_len, y = Hit, fill=Variant)) + geom_bar(stat="identity", position="identity")+
  #   theme_bw() + scale_fill_manual(values=c(D="red", I="green4"),labels=c("Deletion", "Insertion"))+theme(legend.position="none") + 
  #   ylab('Reads')+ xlab('Indel Length (bp)')
  # 
  # penta <- quad + scale_y_continuous(limit = c(min(ggplot_build(quad)$data[[1]]$y),max(abs(ggplot_build(quad)$data[[1]]$y))),label=function(Hit){abs(Hit)})+ 
    # ggtitle("decoy") + theme(plot.title = element_text(color="white", size=15)) ##<-------------- alternate penta histogram to account for low edit frequencies, breaks as fxn of 
  
  
  penta_hope = paste(link,link_no, "_penta.png", sep="")
  ggsave(penta_hope, width = 2, height = 5)
  samp_plot <- image_read(paste(link,link_no,".png",sep=""))
  pent_plot <- image_read(paste(link,link_no,"_penta.png",sep=""))
  
  pen_cat <- image_append(c(samp_plot,pent_plot))
  pen_cat <- image_border(pen_cat, "grey", "2x1")
  pen_cat <-image_annotate(pen_cat, "A", font = 'serif', size = 70,location = "+120+100")
  pen_cat <-image_annotate(pen_cat, "B", font = 'serif', size = 70,location = "+120+1550")
  pen_cat <-image_annotate(pen_cat, "C", font = 'serif', size = 70,location = "+1180+100")
  
  
  
  image_write(pen_cat,(paste(link,link_no,"_Cat.png",sep="")))
  
}

CIITA_1 <- image_read("CIITA1_Cat.png")
CIITA_2 <- image_read("CIITA2_Cat.png")
CIITA_3 <- image_read("CIITA3_Cat.png")
CIITA_4 <- image_read("CIITA4_Cat.png")
# CIITA_5 <- image_read("CIITA_5_Cat.png")


tiger_style <- image_append(c(CIITA_1, CIITA_2))
method_style <- image_append(c(CIITA_3, CIITA_4))
woo <- CIITA
# woo <- image_append(c(CIITA_4,CIITA))
ghost <- image_append(c(tiger_style,method_style,woo),stack=TRUE)


image_write(ghost,(paste(new_title,"_CIITA.png",sep="")))


# ###____________________TRAC_____________________________________________________________###
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

CIGAR = c(0)
Hit = c(0)
Position = start(avsINT1.granges):end(avsINT2.granges)
Variant = c("D")
Sample_Name = c("blank")
Var_len = 1
plot_pos = start(avsINT1.granges):end(avsINT2.granges)
Percentage = c(0)
Dummy_1 = data.frame(CIGAR, Hit, Position, Variant, Sample_Name, Var_len,plot_pos,Percentage)

Variant = c("I")
Dummy_2 = data.frame(CIGAR, Hit, Position, Variant, Sample_Name,Var_len,plot_pos,Percentage)
Dummy_df = rbind(Dummy_2,Dummy_1)



filenames <- list.files("c:/Export/20191202_307_360/Seq_files/", pattern=glob2rx("20191209_TRAC_*_CIGAR_summary.csv"), full.names=TRUE) # 190705_B2M_2_CIGAR_summary.csv

max_list <- list()


for (sample in filenames[c(1:4)])
{
  df<-read.csv(sample, header=TRUE, sep=",")
  melt <- aggregate(df$Percentage, by=list(Position=df$Position), FUN=sum)
  upper <- max(melt$x)
  max_list <- append(upper,max_list)
}

max_max <- max(unlist(max_list))

for (sample in filenames[c(1:4)])
{
  mydataTRAC <- read.csv(sample, header=TRUE, sep=",")
  mydataTRAC$Position= as.numeric(mydataTRAC$Position)+22550515
  mydataTRAC$plot_pos= mydataTRAC$Position+(mydataTRAC$Var_len/2)
  
  if (nrow(mydataTRAC)==0)
  {
    mydataTRAC <- Dummy_df
    del <- Dummy_df
    ins <- Dummy_df
  }
  
  del <- subset(mydataTRAC, Variant=='D')
  ins <- subset(mydataTRAC, Variant=='I')
  
  ### Penta DF Block #####
  
  alt1 <-mydataTRAC
  
  alt <- alt1[,c("Variant","Hit","Var_len")]
  
  deletion <- subset(alt, Variant=='D')
  insertion <- subset(alt, Variant=='I')
  
  if (nrow(deletion)!=0){
    agg_del <- aggregate(Hit ~ Var_len, data=deletion, FUN=sum)
    agg_del$Hit = 0 - as.numeric(agg_del$Hit)
    agg_del$Variant <- "D"
  }else {agg_del <-Dummy_del_agg}
  
  if (nrow(insertion)!=0){
    agg_ins <- aggregate(Hit ~ Var_len, data=insertion, FUN=sum)
    agg_ins$Variant <- "I"
  }else {agg_ins <-Dummy_ins_agg }
  
  
  new_data = rbind(agg_del,agg_ins)
  
  
  ##_________________________________________________________##
  
  mydataTRAC <-mydataTRAC[
    with(mydataTRAC, order(-mydataTRAC$plot_pos,-mydataTRAC$Var_len)),
    ]
  
  
  del <-del[
    with(del, order(-del$plot_pos,-del$Var_len)),
    ]
  
  ins <-ins[
    with(ins, order(-ins$plot_pos,-ins$Var_len)),
    ]
  
  
  link <- substr(sample,47,50)
  link_no <-substr(sample,52,52)
  har <- as.integer((link_no))
  new_title <- toString(samp_name[har])
  
  TRAC_cigar <- ggplot(mydataTRAC,aes(y = Percentage, x = plot_pos, fill = Variant, group = Variant)) + geom_bar(data=del,width=del$Var_len,alpha=.5, stat="identity")+
    geom_bar(data=ins,width=ins$Var_len,alpha=.5, stat="identity")+
    theme_bw() + scale_x_continuous(expand = c(0, 0))+ylab(NULL) +scale_y_continuous(limits=c(0,max_max),expand = c(.01, .01))+
    scale_fill_manual(values=c(D="red", I="green4", X="blue"),breaks=c("D", "I", "X"),labels=c("Deletion", "Insertion", "Substitution"))+
    theme(legend.position= c(0.8, 0.85),axis.text.x  = element_text(size=12))+guides(fill = guide_legend(title = "Edit Type"))
  
  
  
  TRAC_loc <- autoplot(avsEXON.granges, geom_rect(),aes(fill="black",color="black"))+scale_fill_identity() +scale_color_identity()+theme_bw()+
    geom_chevron(avsINT1.granges,color="black", fill=NA)+geom_chevron(avsINT2.granges,color="black", fill=NA)+geom_rect(avsTRACpam.granges,color="black", fill="white")+
    scale_x_continuous(expand = c(0, 0)) + theme(axis.text.x  = element_text(size=6))
  
  sample <- tracks("Percentage"=TRAC_cigar,"TRAC"=TRAC_loc,heights=c(0.1,0.01),label.text.cex=c(1,1),title = new_title,label.bg.fill=('white'))
  
  hope = paste(link,link_no, ".png", sep="")
  ggsave(hope)
  
  if (max(abs(new_data$Hit)) > 1000){
    
    penta <- ggplot(new_data, aes(x = Var_len, y = Hit, fill=Variant)) + geom_bar(stat="identity", position="identity")+
      theme_bw() + scale_fill_manual(values=c(D="red", I="green4"),labels=c("Deletion", "Insertion"))+theme(legend.position="none") + 
      ylab('Reads')+ xlab('Indel Length (bp)') + scale_y_continuous(label=function(Hit){abs(Hit)})+ 
      ggtitle("decoy") + theme(plot.title = element_text(color="white", size=15))
    
  }else {
    penta <- ggplot(new_data, aes(x = Var_len, y = Hit, fill=Variant)) + geom_bar(stat="identity", position="identity")+
      theme_bw() + scale_fill_manual(values=c(D="red", I="green4"),labels=c("Deletion", "Insertion"))+theme(legend.position="none") + 
      ylab('Reads')+ xlab('Indel Length (bp)') + scale_y_continuous(limits = c(-100000,100000))+ 
      ggtitle("decoy") + theme(plot.title = element_text(color="white", size=15))}
  
  penta_hope = paste(link,link_no, "_penta.png", sep="")
  ggsave(penta_hope, width = 2, height = 5)
  samp_plot <- image_read(paste(link,link_no,".png",sep=""))
  pent_plot <- image_read(paste(link,link_no,"_penta.png",sep=""))
  
  pen_cat <- image_append(c(samp_plot,pent_plot))
  pen_cat <- image_border(pen_cat, "grey", "2x1")
  pen_cat <-image_annotate(pen_cat, "A", font = 'serif', size = 70,location = "+120+100")
  pen_cat <-image_annotate(pen_cat, "B", font = 'serif', size = 70,location = "+120+1550")
  pen_cat <-image_annotate(pen_cat, "C", font = 'serif', size = 70,location = "+1180+100")
  
  
  
  image_write(pen_cat,(paste(link,link_no,"_Cat.png",sep="")))
  
}

TRAC_1 <- image_read("TRAC1_Cat.png")
TRAC_2 <- image_read("TRAC2_Cat.png")
TRAC_3 <- image_read("TRAC3_Cat.png")
TRAC_4 <- image_read("TRAC4_Cat.png")
# TRAC_5 <- image_read("TRAC_5_Cat.png")


tiger_style <- image_append(c(TRAC_1, TRAC_2))
method_style <- image_append(c(TRAC_3, TRAC_4))
# woo <- image_append(c(TRAC_4,TRAC))
woo <- TRAC
ghost <- image_append(c(tiger_style,method_style,woo),stack=TRUE)


image_write(ghost,(paste(new_title,"_TRAC.png",sep="")))




# #__________PIK__________________________________________________________________________#

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
Position = start(avsINT1.granges):end(avsINT2.granges)
Variant = c("D")
Sample_Name = c("blank")
Var_len = 1
plot_pos = start(avsINT1.granges):end(avsINT2.granges)
Percentage = c(0)
Dummy_1 = data.frame(CIGAR, Hit, Position, Variant, Sample_Name, Var_len,plot_pos,Percentage)

Variant = c("I")
Dummy_2 = data.frame(CIGAR, Hit, Position, Variant, Sample_Name,Var_len,plot_pos,Percentage)
Dummy_df = rbind(Dummy_2,Dummy_1)

filenames <- list.files("c:/Export/20191202_307_360/Seq_files/", pattern=glob2rx("20191209_PIK_*_CIGAR_summary.csv"), full.names=TRUE) # 190705_B2M_2_CIGAR_summary.csv


max_list <- list()


for (sample in filenames[c(1:4)])
{
  df<-read.csv(sample, header=TRUE, sep=",")
  melt <- aggregate(df$Percentage, by=list(Position=df$Position), FUN=sum)
  upper <- max(melt$x)
  max_list <- append(upper,max_list)
}

max_max <- max(unlist(max_list))

for (sample in filenames[c(1:4)])
{
  mydataPIK <- read.csv(sample, header=TRUE, sep=",")
  mydataPIK$Position= as.numeric(mydataPIK$Position)+8901263
  mydataPIK$plot_pos= mydataPIK$Position+(mydataPIK$Var_len/2)
  
  if (nrow(mydataPIK)==0)
  {
    mydataPIK <- Dummy_df
    del <- Dummy_df
    ins <- Dummy_df
  }
  
  del <- subset(mydataPIK, Variant=='D')
  ins <- subset(mydataPIK, Variant=='I')
  

  
  ### Penta DF Block #####
  
  
  alt1 <-mydataPIK
  
  alt <- alt1[,c("Variant","Hit","Var_len")]
  
  deletion <- subset(alt, Variant=='D')
  insertion <- subset(alt, Variant=='I')
  
  if (nrow(deletion)!=0){
    agg_del <- aggregate(Hit ~ Var_len, data=deletion, FUN=sum)
    agg_del$Hit = 0 - as.numeric(agg_del$Hit)
    agg_del$Variant <- "D"
  }else {agg_del <-Dummy_del_agg}
  
  if (nrow(insertion)!=0){
    agg_ins <- aggregate(Hit ~ Var_len, data=insertion, FUN=sum)
    agg_ins$Variant <- "I"
  }else {agg_ins <-Dummy_ins_agg }
  
  
  new_data = rbind(agg_del,agg_ins)
  
  ##_________________________________________________________##
  
  mydataPIK <-mydataPIK[
    with(mydataPIK, order(-mydataPIK$plot_pos,-mydataPIK$Var_len)),
    ]
  
  
  del <-del[
    with(del, order(-del$plot_pos,-del$Var_len)),
    ]
  
  ins <-ins[
    with(ins, order(-ins$plot_pos,-ins$Var_len)),
    ]
  
   
  link <- substr(sample,47,49)
  link_no <-substr(sample,51,51)
  har <- as.integer((link_no))
  new_title <- toString(samp_name[har])
  
  PIK_cigar <- ggplot(mydataPIK,aes(y = Percentage, x = plot_pos, fill = Variant, group = Variant)) + geom_bar(data=del,width=del$Var_len,alpha=.5, stat="identity")+
    geom_bar(data=ins,width=ins$Var_len,alpha=.5, stat="identity")+
    theme_bw() + scale_x_continuous(expand = c(0, 0))+ylab(NULL) +scale_y_continuous(limits=c(0,100),expand = c(.01, .01))+
    scale_fill_manual(values=c(D="red", I="green4", X="blue"),breaks=c("D", "I", "X"),labels=c("Deletion", "Insertion", "Substitution"))+
    theme(legend.position= c(0.8, 0.85),axis.text.x  = element_text(size=12))+guides(fill = guide_legend(title = "Edit Type"))
  
  
  PIK_loc <- autoplot(PIK.granges, geom = "chevron",offset=0.3,aes(fill="black",color="black"))+scale_fill_identity() +scale_color_identity()+theme_bw()+
    geom_rect(PIK.sgRNA,color="black",fill=NA)+
    scale_x_continuous(expand = c(0, 0))+theme(text = element_text(size=8))
  
  
  sample <- tracks("Percentage"=PIK_cigar,"PIK"=PIK_loc,heights=c(0.1,0.01),label.text.cex=c(1,1),title = new_title,label.bg.fill=('white'))
  
  hope = paste(link,link_no, ".png", sep="")
  ggsave(hope)
  
  if (max(abs(new_data$Hit)) > 1000){
    
    penta <- ggplot(new_data, aes(x = Var_len, y = Hit, fill=Variant)) + geom_bar(stat="identity", position="identity")+
      theme_bw() + scale_fill_manual(values=c(D="red", I="green4"),labels=c("Deletion", "Insertion"))+theme(legend.position="none") + 
      ylab('Reads')+ xlab('Indel Length (bp)') + scale_y_continuous(label=function(Hit){abs(Hit)})+ 
      ggtitle("decoy") + theme(plot.title = element_text(color="white", size=15))
    
  }else {
    penta <- ggplot(new_data, aes(x = Var_len, y = Hit, fill=Variant)) + geom_bar(stat="identity", position="identity")+
      theme_bw() + scale_fill_manual(values=c(D="red", I="green4"),labels=c("Deletion", "Insertion"))+theme(legend.position="none") + 
      ylab('Reads')+ xlab('Indel Length (bp)') + scale_y_continuous(limits = c(-100000,100000))+ 
      ggtitle("decoy") + theme(plot.title = element_text(color="white", size=15))}
  
  
  
  penta_hope = paste(link,link_no, "_penta.png", sep="")
  ggsave(penta_hope, width = 2, height = 5)
  samp_plot <- image_read(paste(link,link_no,".png",sep=""))
  pent_plot <- image_read(paste(link,link_no,"_penta.png",sep=""))
  
  pen_cat <- image_append(c(samp_plot,pent_plot))
  pen_cat <- image_border(pen_cat, "grey", "2x1")
  pen_cat <-image_annotate(pen_cat, "A", font = 'serif', size = 70,location = "+120+100")
  pen_cat <-image_annotate(pen_cat, "B", font = 'serif', size = 70,location = "+120+1550")
  pen_cat <-image_annotate(pen_cat, "C", font = 'serif', size = 70,location = "+1180+100")
  
  
  
  image_write(pen_cat,(paste(link,link_no,"_Cat.png",sep="")))
  
}

PIK_1 <- image_read("PIK1_Cat.png")
PIK_2 <- image_read("PIK2_Cat.png")
PIK_3 <- image_read("PIK3_Cat.png")
PIK_4 <- image_read("PIK4_Cat.png")
# PIK_5 <- image_read("PIK5_Cat.png")


tiger_style <- image_append(c(PIK_1, PIK_2))
method_style <- image_append(c(PIK_3, PIK_4))
# woo <- image_append(c(PIK_4,PIK))
woo <- PIK
ghost <- image_append(c(tiger_style,method_style,woo),stack=TRUE)


image_write(ghost,(paste(new_title,"_PIK.png",sep="")))

