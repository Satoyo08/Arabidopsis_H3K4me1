# what this script do: make plots in Figure.1
source('custom_functions.r',echo=FALSE)
library(tidyverse)


# ------Fig1.a structure------------------------
# read data : domain annotations were downloaded from UniProt as jsonfiles, and converted into csv. ATX1_len; length of ATX1 protein (AA)
ATX1_st<-read.csv("data/Figure1/ATX1.csv",header=TRUE);ATX1_len<-1062
ATX2_st<-read.csv("data/Figure1/ATX2.csv",header=TRUE);ATX2_len<-1083
ATX3_st<-read.csv("data/Figure1/ATX3.csv",header=TRUE);ATX3_len<-1018
ATX4_st<-read.csv("data/Figure1/ATX4.csv",header=TRUE);ATX4_len<-1027
ATX5_st<-read.csv("data/Figure1/ATX5.csv",header=TRUE);ATX5_len<-1043
ATXR7_st<-read.csv("data/Figure1/ATXR7.csv",header=TRUE);ATXR7_len<-1423
ATXR3_st<-read.csv("data/Figure1/ATXR3.csv",header=TRUE);ATXR3_len<-2335

STRUCTURE<-list(ATX1_st,ATX2_st,ATX3_st,ATX4_st,ATX5_st,ATXR7_st,ATXR3_st);length<-c(ATX1_len,ATX2_len,ATX3_len,ATX4_len,ATX5_len,ATXR7_len,ATXR3_len)
All_st<-rbind(ATX1_st,ATX2_st,ATX3_st,ATX4_st,ATX5_st,ATXR7_st,ATXR3_st)
# set colors
pallets<-as.data.frame(c(brewer.pal(8, "Set1"),brewer.pal(11, "Paired")[c(9,11,1)]))
rownames(pallets)<-unique(All_st$text)[c(7,1,2,3,4,9,5,6,8,10,11)];colnames(pallets)<-c("colors")
#plot the structure
par(mar=c(0,0,0,0))
plot(1,1,xlim=c(0,max(length)),ylim=c(0,length(length)*10),cex=0,type='n',axes=F)
for(i in 1:length(length)){
  rect(0,4+(i-1)*10,length[i],6+(i-1)*10,col="gray",border="gray")
  rect(STRUCTURE[[i]]$start,3+(i-1)*10,STRUCTURE[[i]]$end,7+(i-1)*10,col=as.character(pallets[as.character(STRUCTURE[[i]]$text),]),border=as.character(pallets[as.character(STRUCTURE[[i]]$text),]))
  text(STRUCTURE[[i]]$start,3+(i-1)*10,STRUCTURE[[i]]$text)
}
par(mar=c(1,1,1,1));legend("top",legend=rownames(pallets),fill=as.character(pallets$colors))
# ------------ Fig.1b --------------
# the metaplots were generated with ngs.plot.r with the following command
# bams/lists are not available in the Github (if you need them contact o.satoyo@bs.s.u-tokyo.ac.jp)
# ngs.plot.r -G Tair10 -R genebody -C revgen7mutans_config.txt -O revgen7_top3000 -GO none -L 500
# revgen7mutans_config.txt ---
#atx1_K4me1_1_sorted.bam:Col_K4me1_1_sorted.bam ~/refs/atx1_top3000_gene_list.txt "atx1"
#atx2_K4me1_1_sorted.bam:Col_K4me1_1_sorted.bam ~/refs/atx2_top3000_gene_list.txt "atx2"
#atx3_K4me1_1_sorted.bam:Col_K4me1_1_sorted.bam ~/refs/atx3_top3000_gene_list.txt "atx3"
#atx4_K4me1_1_sorted.bam:Col_K4me1_1_sorted.bam ~/refs/atx4_top3000_gene_list.txt "atx4"
#atx5_K4me1_1_sorted.bam:Col_K4me1_1_sorted.bam ~/refs/atx5_top3000_gene_list.txt "atx5"
#atxr7_K4me1_1_sorted.bam:Col_K4me1_1_sorted.bam ~/refs/atxr7_top3000_gene_list.txt "atxr7"

# -----------Fig.1 c heatmap
# the heatmap was generated with ngs.plot.r,
# ngs.plot.r -G Tair10 -R genebody -C atx1_2_12_r7_12r7_all_config.txt -O atx1_2_12_r7_12r7_all -GO none -L 500 -SC global -RR 250 -CO gold:springgreen3:midnightblue
# atx1_2_12_r7_12r7_all_config.txt ---
#atx1_K4me1_1_sorted.bam:Col_K4me1_1_sorted.bam ~/refs/gene_list/sortedIDlist_K4_12r7-WT.txt "atx1"
#atx2_K4me1_1_sorted.bam:Col_K4me1_1_sorted.bam ~/refs/gene_list/sortedIDlist_K4_12r7-WT.txt "atx2"
#/homeB/satoyo08/ChIP2/x12_me1.bam:/homeB/satoyo08/ChIP2/WT_me1.bam ~/refs/gene_list/sortedIDlist_K4_12r7-WT.txt "atx12"
#atxr7_K4me1_1_sorted.bam:Col_K4me1_1_sorted.bam ~/refs/gene_list/sortedIDlist_K4_12r7-WT.txt "atxr7"
#/homeB/satoyo08/ChIP1/x127_me1.bam:/homeB/satoyo08/ChIP1/WT_me1_2.bam  ~/refs/gene_list/sortedIDlist_K4_12r7-WT.txt "atx12r7"

# -----------Fig1.d scatter plots -------------
# read data : use RPKM normalized count for H3K4me1 ChIP-seq, RPM normalized count for H3K4me2,3.
RPKM<-read.table("data/Figure1/ChIP1_RPKM.txt",header=TRUE)
RPM<-read.table("data/Figure1/ChIP1_RPM.txt",header=TRUE)

#define **-marked genes
DF<-gonly(RPKM);S<-DF[order(DF$x127_me1-DF$WT_me1),];new_red<-S[1:3000,1]
DF<-gonly(RPM);S<-DF[order(DF$x345_me2-DF$WT_me2),];new_green<-S[1:3000,1]
DF<-gonly(RPM);S<-DF[order(DF$r3_me3-DF$WT_me3),];new_blue<-S[1:3000,1]

#--3x3 plot
AL<-merge(RPKM,RPM,by='ID');AL[,ncol(AL)+1]<-0
AL[AL$ID %in% new_red,]$V34<-AL[AL$ID %in% new_red,]$V34+1
AL[AL$ID %in% new_green,]$V34<-AL[AL$ID %in% new_green,]$V34+2
AL[AL$ID %in% new_blue,]$V34<-AL[AL$ID %in% new_blue,]$V34+4


wt_col<-c(7,7,7,24,24,24,25,25,25)
mutant_col<-c(11,15,3,28,32,20,29,33,21)
par(mfrow=c(3,3),mar=c(2,2,1,1),cex.lab=1.2)
for(i in 1:length(wt_col)){
  plot_shu(AL[,wt_col[i]],AL[,mutant_col[i]],col=rgb_pal[AL[,ncol(AL)] + 1],pch=20,lwd=0.01) 
}


# -----------Fig1.f venn-diagram -------------
## venn diagram was drawn with illustrator
# hyper geometric test
q<-length(intersect(new_red,new_blue)) # intersect of the **-marked genes defined above
m<-3000
n<-nrow(gonly(RPKM))-m
k<-3000
qh<-phyper(q, m, n, k,lower.tail = TRUE,log.p = TRUE) 
nrow(gonly(RPKM))*(3000/nrow(gonly(RPKM)))^2 #Expectation


