source('custom_functions.r');library(randomForest)

# SFig.16ab----
# scatter plots
RPKM<-read.table('data/Figure1/ChIP1_RPKM.txt',header=T)
ChIP3_RPKM<-read.table("data/refs/ChIP3_RPKM.txt",header=T)
ChIP11_K36_RPKM<-read.table('data/refs/ChIP11_ashhetc_RPKM.txt',header=T)
ChIP11_RPKM2<-read.table('data/refs/ChIP11_triple_RPKM.txt',header=T) 

D13<-gonly(merge(RPKM,ChIP3_RPKM,by="ID"))[,1:36]
DF<-gonly(merge(merge(ChIP11_K36_RPKM,D13,by='ID'),ChIP11_RPKM2,by='ID'))
plot_shu(D13$x127_me1-D13$WT_me1,D13$x127_K36me3-D13$Col_K36me3,limmaxx=0.9995,limmaxy=0.9995,limminx=0.0005,limminy=0.0005)
cor.test(D13$x127_me1-D13$WT_me1,D13$x127_K36me3-D13$Col_K36me3,method='spearman')

plot_shu(DF$x12r7_me1-DF$WT_K4me1,DF$x12r7_K36me3-DF$WT_H3K36me3,limmaxx=0.9995,limmaxy=0.9995,limminx=0.0005,limminy=0.0005)
cor.test(DF$x12r7_me1-DF$WT_K4me1,DF$x12r7_K36me3-DF$WT_H3K36me3,method='spearman')

plot_shu(DF$sdg8_H3K4me1-DF$WT_SDG8_H3K4me1,DF$ashh2_H3K36me3-DF$WT_H3K36me3,limmaxx=0.9995,limmaxy=0.9995,limminx=0.0005,limminy=0.0005)
cor.test(DF$sdg8_H3K4me1-DF$WT_SDG8_H3K4me1,DF$ashh2_H3K36me3-DF$WT_H3K36me3,method='spearman')

# SFig.16ab heatmaps
# calculated with ngs.plot.r
