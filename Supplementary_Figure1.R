# what this script do: make plots in Supplementary Fig.1
source('custom_functions.r',echo=FALSE)
library(tidyverse)
library(multcomp)
# ------------S Fig.1 a --------------
Ch02<-read.table("data/Figure1/RPKM.txt",header=TRUE)
cairo_ps(file = "singe_scatter.eps", onefile = FALSE, fallback_resolution = 300,width=6,height=1)
par(mfrow=c(1,6),mar=c(2,2,2,2),cex=0.2,cex.axis=1.5)
for(i in 5:10){plot_shu(Ch02[,11],Ch02[,i]) };dev.off()

# ------S Fig1.b clustering of mutants------------------------
# read data sets; RPM-normalized H3K4me1 ChIP-seq counts in genebodies. Genebodies were split into six bin(V2~V7). 
atx1<-read.table("data/Figure1/atx1_K4me1_1hexa_RPM.txt",header=F)
atx2<-read.table("data/Figure1/atx2_K4me1_1hexa_RPM.txt",header=F)
atx3<-read.table("data/Figure1/atx3_K4me1_1hexa_RPM.txt",header=F)
atx4<-read.table("data/Figure1/atx4_K4me1_1hexa_RPM.txt",header=F)
atx5<-read.table("data/Figure1/atx5_K4me1_1hexa_RPM.txt",header=F)
atxr7<-read.table("data/Figure1/atxr7_K4me1_1hexa_RPM.txt",header=F)
Col<-read.table("data/Figure1/Col_K4me1_1hexa_RPM.txt",header=F)
RAW<-list(Col,atx1,atx2,atx3,atx4,atx5,atxr7)
q<-0.07865964
# for each mutant, the top 3'000 affected genes are extracted and normalized.
m<-c();sd<-c();SA<-list();TOPSA<-list();for(i in 1:7){ 
  RAW[[i]][,8]<-apply(RAW[[i]][,2:7],1,sum) 
  SA[[i]]<-data.frame(RAW[[i]][,1],RAW[[i]][,2:8]-RAW[[1]][,2:8]) 
  TOPSA[[i]]<-subset(SA[[i]],SA[[i]][,8]< quantile(SA[[i]][,8],q))
  m<-cbind(m,apply(TOPSA[[i]][,2:7],2,mean)) 
  sd<-cbind(sd,apply(TOPSA[[i]][,2:7],2,sd))
}

d<-c();for(i in 2:7){d<-cbind(d,as.vector(as.matrix(TOPSA[[i]][,2:7])))} 

# scale the 6 bins within each gene
sc<-list()
msc=c()
for(i in 2:7){
  sc[[i]]<-t(scale(t(TOPSA[[i]][,2:7])))
  msc<-cbind(msc,apply(sc[[i]],2,mean))
}
d<-c();for(i in 2:7){d<-cbind(d,as.vector(as.matrix(sc[[i]])))} 
Mr<-d;names<-c("atx1", "atx2"  ,"atx3" ,"atx4", "atx5","atxr7");

# plot the dendrogram
Mr<-msc[,c(6,1,2,3,4,5)];colnames(Mr)<-names[c(6,1,2,3,4,5)]
Mr<-scale(Mr)
h<-dist(t(Mr),method="euclidean");hc<-hclust(h,method="centroid");hcd <- as.dendrogram(hc)
par(mar=c(0,0,0,3),cex=2)
plot(hcd,horiz=TRUE)


# -----------S Fig1.c ---------
ChIP2_RPKM<-read.table("data/Figure1/ChIP2_RPKM.txt",header=TRUE)
ChIP2_RPM<-read.table("data/Figure1/ChIP2_RPM.txt",header=TRUE)
wt<-c(3,3,3,4,4,4,5,5,5);muts<-c(c(7,11,15),c(7,11,15)+1,c(7,11,15)+2)
jpeg('doubles_scatter.jpeg',height=600*a,width=600*a,res=300)
par(mfrow=c(3,3),mar=c(3,3,0.1,0.1),cex=0.4,cex.axis=2)
for(i in 1:3){plot_shu(ChIP2_RPKM[,wt[i]],ChIP2_RPKM[,muts[i]])}
for(i in 4:length(wt)){
  print(colnames(ChIP2_RPM)[c(wt[i],muts[i])])
  plot_shu(ChIP2_RPM[,wt[i]],ChIP2_RPM[,muts[i]])
}
dev.off()

# -----------S Fig.d ---------
# for the source data please see Supplementalry Data 2 : Summary of Spike-in analysis

# -------- S Fig.1 e -------------
arapo<-read.table("data/refs/araport11_all_sorted.bed",header=F)[,2:4];colnames(arapo)<-c('start','end','ID')
RNA8<-read.table("data/Figure1/PerGene_RNA8_araport.txt",header=TRUE)
cDNA=read.table("data/refs/cDNA_length_araport.txt",header=F);colnames(cDNA)<-c('ID','cDNAlength') # protein coding,with primary transcript only
# normalize the mRNA-seq count
sum<-apply(RNA8[,2:ncol(RNA8)],2,sum); RPM<-data.frame(RNA8$ID,t(t(RNA8[,2:ncol(RNA8)])*10^6/sum));colnames(RPM)[1]<-"ID"
AL<-merge(RPM,cDNA,by="ID");FPKM<-data.frame(AL[,1],AL[,2:(ncol(AL)-1)]*1000/AL$cDNAlength);colnames(FPKM)[1]<-"ID"
AL<-merge(FPKM,arapo,by='ID')
# prep the dataframe, group by **-marked genes. 
PC<-gonly(AL)[,c(1:ncol(AL))];R<-AL[AL$ID %in% new_red,];G<-AL[AL$ID %in% new_green,];B<-AL[AL$ID %in% new_blue,]
k<-ncol(AL);PC[,k+1]<-'0';R[,k+1]<-'1';G[,k+1]<-'2';B[,k+1]<-'3'
ALL<-rbind(PC,R,G,B);ALL[,ncol(PC)+1]<-log2(apply(ALL[,8:10],1,mean));ALL[,ncol(PC)+2]<-log2(ALL$end-ALL$start);

# violinplot of log2 length,
ggplot(data=ALL)+
  aes(x=V13,y=V15,fill=V13)+
  geom_violin(bw="nrd0",color='gray70')+
  theme_classic()+
  scale_fill_manual(values=c('gray70',rgb_pal[c(2,3,5)]))+
  theme(panel.background = element_rect(fill = "white",color='black'),axis.text.y=element_text(angle=90,hjust=0.5),
        axis.title=element_blank(),legend.position = 'none',axis.line = element_blank()
  )+
  geom_boxplot(width=0.1,color='gray50',fill='black',outlier.size=0.2,outlier.color = 'gray50')

# Dunnett test (multiple comparison of means)
fx=factor(ALL$V13)
vx=ALL$V15
summary(glht(aov(vx~fx),linfct=mcp(fx="Dunnett")))

# log2 FPKM
ggplot(data=ALL)+
  aes(x=V13,y=V14,fill=V13)+
  geom_violin(bw="nrd0",color='gray70')+
  theme_classic()+
  scale_fill_manual(values=c('gray70',rgb_pal[c(2,3,5)]))+
  theme(panel.background = element_rect(fill = "white",color='black'),axis.text.y=element_text(angle=90,hjust=0.5),
        axis.title=element_blank(),legend.position = 'none',axis.line = element_blank()
  )+
  geom_boxplot(width=0.1,color='gray50',fill='black',outlier.size=0.2,outlier.color = 'gray50')

# Dunnett test, removing log2(0)=-Inf 
fx=factor(ALL[!ALL$V14== -Inf,]$V13)
vx=ALL[!ALL$V14== -Inf,]$V14
summary(glht(aov(vx~fx),linfct=mcp(fx="Dunnett")))

# ------------S Fig. f,g,h ------------------
# the metaplots and heatmaps were generated with ngs.plot.r with the following commands
# ngs.plot.r -G Tair10 -R genebody -C the_blue.txt -O the_blue -GO none -L 500 &
# /homeB/satoyo08/ChIP1/WT_me1_2.bam ~/refs/gene_list/the_new_blue.txt "WT_H3K4me1"
# /homeB/satoyo08/ChIP1/WT_me2.bam ~/refs/gene_list/the_new_blue.txt "WT_H3K4me2"
# /homeB/satoyo08/ChIP1/WT_me3.bam ~/refs/gene_list/the_new_blue.txt "WT_H3K4me3"
# /homeB/satoyo08/ChIP1/r3_me1.bam ~/refs/gene_list/the_new_blue.txt "atxr3_H3K4me1"
# /homeB/satoyo08/ChIP1/r3_me2.bam ~/refs/gene_list/the_new_blue.txt "atxr3_H3K4me2"
# /homeB/satoyo08/ChIP1/r3_me3.bam ~/refs/gene_list/the_new_blue.txt "atxr3_H3K4me3"
# 
# ngs.plot.r -G Tair10 -R genebody -C the_red.txt -O the_red -GO none -L 500 &
# /homeB/satoyo08/ChIP1/WT_me1_2.bam ~/refs/gene_list/the_new_red.txt "WT_H3K4me1"
# /homeB/satoyo08/ChIP1/WT_me2.bam ~/refs/gene_list/the_new_red.txt "WT_H3K4me2"
# /homeB/satoyo08/ChIP1/WT_me3.bam ~/refs/gene_list/the_new_red.txt "WT_H3K4me3"
# /homeB/satoyo08/ChIP1/x127_me1.bam ~/refs/gene_list/the_new_red.txt "atx12r7_H3K4me1"
# /homeB/satoyo08/ChIP1/x127_me2.bam ~/refs/gene_list/the_new_red.txt "atx12r7_H3K4me2"
# /homeB/satoyo08/ChIP1/x127_me3.bam ~/refs/gene_list/the_new_red.txt "atx12r7_H3K4me3"
# 
# ngs.plot.r -G Tair10 -R genebody -C the_green.txt -O the_green -GO none -L 500 &
# /homeB/satoyo08/ChIP1/WT_me1_2.bam ~/refs/gene_list/the_new_green.txt "WT_H3K4me1"
# /homeB/satoyo08/ChIP1/WT_me2.bam ~/refs/gene_list/the_new_green.txt "WT_H3K4me2"
# /homeB/satoyo08/ChIP1/WT_me3.bam ~/refs/gene_list/the_new_green.txt "WT_H3K4me3"
# /homeB/satoyo08/ChIP1/x345_me1.bam ~/refs/gene_list/the_new_green.txt "atx345_H3K4me1"
# /homeB/satoyo08/ChIP1/x345_me2.bam ~/refs/gene_list/the_new_green.txt "atx345_H3K4me2"
# /homeB/satoyo08/ChIP1/x345_me3.bam ~/refs/gene_list/the_new_green.txt "atx345_H3K4me3"

