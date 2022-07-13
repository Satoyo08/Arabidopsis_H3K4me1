
source('custom_functions.r');source('RF_functions.R')
library(ROCR);library(tidyverse);library(eulerr)

# ------------ SFig.8a ------------
# Random Forest models were trained and saved by RandomForest_arabidopsis.Rmd
load("data/Figure3/ATXR7_with_K36me3_ChIP8_rep_FLD") #  load randomforest models. 'ChIP8'=Fig.3, 'ChIP7'=biological replicates shown in Supplementry Fig.5
model<-ATXR7_with_K36me3_ChIP8_rep_FLD

# sort features 
sort<-c(1,2,24,3,66,45,25,4,67,46,26,5,68,47,36,15,78,57,34,13,76,55,33,12,75,54,35,14,77,56,29,8,71,50,31,10,73,52,30,9,72,51,27,6,69,48,42,21,84,63,43,22,85,64,44,23,86,65,37,16,79,58,28,7,70,49,40,19,82,61,41,20,83,62,32,11,74,53,87,93,90,88,94,91,89,95,92,38,17,80,59,39,18,81,60,96)
col_key<-c(10,9,rep(c(3:6),19),rep(c(3,4,6),3),rep(c(3:6),2),4)
space<-c(0,2,c(2,0,0,0),rep(c(0.5,0,0,0),18),rep(c(0.5,0,0),3),rep(c(0.5,0,0,0),2),2)

# plot bars
par(mar=c(0.5,2,0.5,0))
IMPs<-apply(model$Importance,1,mean)[sort];sds<-apply(model$Importance,1,sd)[sort]
bp<-barplot(IMPs,col=pallets[col_key],space=space,names.arg=F,ylim=c(0,max(IMPs+max(sds))));lims<-par('usr')
plot(0,type='n',axes=FALSE,ann=FALSE,xlim=lims[1:2]+c(4,0),ylim=lims[3:4])
bp<-barplot(IMPs,col=pallets[col_key],space=space,names.arg=F,ylim=c(0,max(IMPs+max(sds))),add=T)
segments(bp,IMPs-sds,bp,IMPs+sds) # plot error bars
points(rep(bp,5),as.vector(model$Importance[sort,]),pch=16,col='gray40',cex=0.5) # plot each data points


# SFig.8b ---------------
FLD_sa<-read.table('RF_data_fld.txt',header=T) # https://github.com/soinagak/FLD2021/blob/main/RF_data_fld.txt
S<-TAGg[order(-TAGg$ATXR7_ChIP8TES),];posID<-S[1:3000,1]
col=pallets[6]
AL<-FLD_sa
AL[,ncol(AL)+1]<-'A'
AL[AL$ID %in% posID,ncol(AL)]<-'B'
k<-which(colnames(AL)=='TTS_fld')
DF<-data.frame(AL[,k],AL[,ncol(AL)])
colnames(DF)<-c('data','negpos')
vio<-ggplot(data=DF)+ 
  aes(x=negpos,y= -data,fill=negpos)+
  geom_violin(bw="nrd0",color='gray70',scale='width')+
  theme_classic()+
  scale_fill_manual(values=c('gray70',col))+
  theme(panel.background = element_rect(fill = "white",color='black'),axis.text.y=element_text(angle=90,hjust=0.5, family ="Helvetica"),
        axis.text.x=element_blank(),axis.title=element_blank(),legend.position = 'none',axis.line = element_blank())+
  geom_boxplot(width=0.1,color='gray50',fill='black',outlier.size=0.2,outlier.color = 'gray50')
ggsave(file="R7_bound_FLD_vio.pdf",vio,width=1.2,height=2.5)


# SFig.8c-----------
Ch02<-read.table("data/Figure1/RPKM.txt",header=TRUE)
r7<-gonly(Ch02[Ch02$Col-Ch02$atxr7>6.5,])$ID #638
fld<-read.table('data/refs/fld_K4me1_up_diff_sqrt1_281.txt',header=F)$V1
gl<-list(r7=r7,fld=fld)
plot(euler(gl),shape = "euler",quantities = TRUE) 

q<-length(intersect(r7,fld));q
m<-length(r7)
n<-27409-m
k<-length(fld)
phyper(q, m, n, k,lower.tail = FALSE,log.p = FALSE)
