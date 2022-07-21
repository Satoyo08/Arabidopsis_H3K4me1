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

# SFig.16c
#overlap between ATX1/2/R7-marked genes and hypo-H3K36me3 genes in ashh mutatns
intersects<-list.files('/Users/Satoyo/data/ChIP11/ashh_etc/intersect',pattern='ash');setwd('/Users/Satoyo/data/ChIP11/ashh_etc/intersect')
intersect_files<-intersects[grep('reduction',intersects)];intersect_files
target<-new_red
HYPG<-as.data.frame(matrix(nrow=length(intersect_files),ncol=7));INT<-HYPG;enr<-HYPG;K<-HYPG
for(i in 1:length(intersect_files)){
  ids<-unique(read.table(as.character(intersect_files[i]))$V4)
  q<-length(intersect(ids,target)) #number of white balls drawn without replacement from an urn which contains both black and white balls
  m<-length(target)# the number of white balls in the urn.
  n<-27408-m# the number of black balls in the urn.
  k<-length(ids) # Yの数　k	the number of balls drawn from the urn.
  e<-k*m/(n+m)#Expected
  q/e #enrichment
  pv<-phyper(q, m, n, k,lower.tail = FALSE,log.p = FALSE) # probability lower.tail=TRUE だと、今の引きよりも少なくひく確率
  HYPG[i,]<-c(as.character(intersect_files[i]),k,m,q,e,q/e,pv)
}
HYPG

# plot-bubble
ds<-HYPG
xax<-rep(c(1,2),5)
yax<-c(1,1,2,2,3,3,4,4,5,5)

par(mar=c(0.5,0.5,0.5,0.5))
colf<-colorRampPalette(rev(c(pallets[1],pallets[5],'gray80')))
pval<- -log10(as.numeric(ds$V7))
pval[pval>5]<-5 # cutoff min p-val
plot(1,1,type='n',ylim=c(0,6),xlim=c(0,3),xaxt='n',yaxt='n')
points(xax,rev(yax),cex=as.numeric(ds$V6)*1.3,col=colf(1002)[floor(pval*1000/5)+2],pch=16)
text(xax,rev(yax),ds$V2,cex=0.6)
text(xax,rev(yax)-0.15,paste('(',ds$V4,')'),cex=0.6)

#plot legends
plot(1,1,type='n',ylim=c(0,6),xlim=c(0,3),xaxt='n',yaxt='n')
points(rep(1,6),6:1,cex=c(0.25,0.5,1,2,4,8)*1.3,col='gray',pch=16)
n<-100;bin<-5/n
rect(rep(2,n),(1:n)*bin-bin,rep(2,n)+0.5,(1:n)*bin,
     col=colf(1002)[(n:0)*bin*1000/5+2],# range -log(p) 5~0
     border=colf(1002)[(n:0)*bin*1000/5+2])
rect(2,1*bin-bin,2.5,n*bin,border=1)
text(rep(1,6),6:1,c('0.25','0.5','1','2','4','8'))
text(rep(2,6)+0.7,0:5,c('1e-05','1e-04','1e-03','0.01','0.1','1'))
segments(rep(2,6)+0.5,0:5,rep(2,6)+0.55,0:5)

