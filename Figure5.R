# what this script do: plot 2D density heatmaps
# Corresponding Figures: Fig. 5, Fig. 6e-h and such

# functions --------------------------------------------------------------------------------------------------------------
source('custom_functions.r')
library(MASS);library(rgl);library(viridis)
return_rhos<-function(rnas,K4s,AL){
  rhos<-c()
  for(i in 1:length(K4s)){
    rna<-rnas[i];k4<-K4s[i]
    rho<-cor.test(AL[,rna],AL[,k4],method="spearman")
    rhos<-rbind(rhos,rho$estimate)
  }
  return(rhos)
}
plot_2d<-function(rnas,K4s,AL){
  avails<-c();for(i in rnas){print(colnames(AL)[i]);D<-AL[AL[,i]>0,];print(nrow(D));avails=c(avails,nrow(D))};avails
  min_sample<-min(avails)
  
  # Get the maximum and minimum rd densities of the data sets to be compared and store them in zmax and zmin (for ranges in heatmap)
  zmins<-c();zmaxs<-c()
  for(i in 1:length(rnas)){
    rna<-rnas[i];k4<-K4s[i]
    Dp<-AL[AL[,rna]>0,];D<-Dp[sample(1:nrow(Dp),min_sample),];rd=kde2d(x=rank(D[,k4]),y=rank(D[,rna]),h=c(3300,3300),n=70)# making the density data rd, using MASS library. n and h tunes the smoothness of the contour.
    zmins<-c(zmins,min(rd$z));zmaxs<-c(zmaxs,max(rd$z))
  }
  zmin<-min(zmins);zmax<-max(zmaxs)
  bins=20;hmp_pallets=viridis(bins)
  for(i in 1:length(rnas)){
    rna<-rnas[i];k4<-K4s[i]
    #title<-paste('AVG_',colnames(AL)[rna],colnames(AL)[k4],sep='_') ;
    #jpeg(paste(title,'_color_bar.jpeg',sep=''),height=1000,width=1180,res=300);par(mar=c(0.1,0.1,0.1,0.1),cex=1) 
    Dp<-AL[AL[,rna]>0,];D<-Dp[sample(1:nrow(Dp),min_sample),];rd=kde2d(x=rank(D[,k4]),y=rank(D[,rna]),h=c(3300,3300),n=70)
    filled.contour(rd$x,rd$y,rd$z,col=hmp_pallets,evels=seq(from=zmin,to=zmax,by=(zmax-zmin)/(bins-1)),cex=0.2);
    #dev.off()
  }
  
}
#--------------------------------------------------------------------------------------------------------------


## for atx1/2,atxr7 (Fig.5) -------------------------------------------------------------------------------------------------------
RNA8<-read.table("data/Figure1/PerGene_RNA8_araport.txt",header=TRUE)
atx_RPM<-read.table('data/refs/atx12r7_me3H3_RPM.txt',header=T)
cDNA=read.table("refs/cDNA_length_araport.txt",header=F);colnames(cDNA)<-c('ID','cDNAlength') # protein coding,with primary transcript only, same 行数as arapo$V6==1

## prep the dataframe: merge and scale
sum<-apply(RNA8[,2:ncol(RNA8)],2,sum)
RPM<-data.frame(RNA8$ID,t(t(RNA8[,2:ncol(RNA8)])*10^6/sum));colnames(RPM)[1]<-"ID"
AL<-merge(RPM,cDNA,by="ID")
FPKM<-data.frame(AL[,1],AL[,2:(ncol(AL)-1)]*1000/AL$cDNAlength);colnames(FPKM)[1]<-"ID"
AL<-gonly(merge(FPKM,atx_RPM,by="ID"))#;AL_2<-merge(AL_1,Ch02,by="ID");AL_3<-merge(AL_2,RPKM,by="ID");AL<-gonly(AL_3)
## choose columns of interest
rnas<-c(2,3,4,5,6,7,8,9,10);colnames(AL)[rnas]
K4s<-c(rep(which(colnames(AL)=="atx12_H3K4me1"),3),rep(which(colnames(AL)=="atxr7_H3K4me1"),3),rep(which(colnames(AL)=="WT_H3K4me1"),3));colnames(AL)[K4s]
## average three rna-seq replicates
AVG<-data.frame(AL$ID,apply(AL[,c(2,3,4)],1,mean),apply(AL[,c(5,6,7)],1,mean),apply(AL[,c(8,9,10)],1,mean),AL[,c(15,14,11)]);colnames(AVG)[1:4]<-c('ID','atx12','atxr7','WT')
## calculate rho
rhos<-return_rhos(rnas,K4s,AL) # rho of mRNA-H3K4me1 in {atx12,atxr7,wt} x 3 (three replicates of mRNA-seq).
rhos<-return_rhos(rep(c(8,9,10),2),K4s[1:6],AL) 
## plot 2d density
plot_2d(c(2,3,4,4,4),c(5,6,7,5,6),AVG)
t.test(rhos[7:9],rhos[1:3]);t.test(rhos[7:9],rhos[4:6]);


#for atx1/2,atxr7, replicates (Sup.Fig 12) -------------------------------------------------------------------------------------------------------
RNA8<-read.table("data/Figure1/PerGene_RNA8_araport.txt",header=TRUE)
C2<-read.table("data/Figure1/ChIP2_RPKM.txt",header=TRUE)
Ch02<-read.table("data/Figure1/RPKM.txt",header=TRUE)
RPKM<-read.table("data/ChIP1_RPKM.txt",header=TRUE)

## prep the dataframe: merge and scale
sum<-apply(RNA8[,2:ncol(RNA8)],2,sum)
RPM<-data.frame(RNA8$ID,t(t(RNA8[,2:ncol(RNA8)])*10^6/sum));colnames(RPM)[1]<-"ID"
AL<-merge(RPM,cDNA,by="ID")
FPKM<-data.frame(AL[,1],AL[,2:(ncol(AL)-1)]*1000/AL$cDNAlength);colnames(FPKM)[1]<-"ID"
AL_1<-merge(FPKM,C2,by="ID");AL_2<-merge(AL_1,Ch02,by="ID");AL_3<-merge(AL_2,RPKM,by="ID");AL<-gonly(AL_3)
## choose columns of interest
rnas<-c(2,3,4,5,6,7,8,9,10,8,9,10);colnames(AL)[rnas]
K4s<-c(rep(which(colnames(AL)=="x12_me1"),3),rep(which(colnames(AL)=="atxr7"),3),rep(which(colnames(AL)=="Col"),3),rep(which(colnames(AL)=="WT_me1.x"),3));colnames(AL)[K4s]
## average three rna-seq replicates
AVG<-data.frame(AL$ID,apply(AL[,c(2,3,4)],1,mean),apply(AL[,c(5,6,7)],1,mean),apply(AL[,c(8,9,10)],1,mean),AL[,c(16,35,36,12)]);colnames(AVG)[1:4]<-c('ID','atx12','atxr7','WT')
## calculate rho
rhos<-return_rhos(rnas,K4s,AL) # rho of mRNA-H3K4me1 in {atx12,atxr7,wt} x 3 (three replicates of mRNA-seq).
rhos<-return_rhos(rep(c(8,9,10),2),K4s_12r7[1:6],AL) 
## plot 2d density
plot_2d(c(2,3,4,4,4,4),c(5,6,7,8,5,6),AVG)


# for Chen 2017's atx3-1/4/5 -------------------------------------------------------------------------------------------------------
Chen_RNA<-read.table("data/Figure5/PerGene_Chen2017_arapo.txt",header=T)
Chen_IP<-read.table("data/Figure5/RPM_Chen2017.txt",header=T) 
sum<-apply(Chen_RNA[,2:ncol(Chen_RNA)],2,sum)
Chen_RNA<-data.frame(Chen_RNA$ID,t(t(Chen_RNA[,2:ncol(Chen_RNA)])*10^6/sum));colnames(Chen_RNA)[1]<-"ID"
AL<-merge(Chen_RNA,cDNA,by="ID");FPKM<-data.frame(AL[,1],AL[,2:(ncol(AL)-1)]*1000/AL$cDNAlength);colnames(FPKM)[1]<-"ID"
AL_1<-merge(FPKM,Chen_IP,by="ID");AL<-gonly(AL_1)
AVG<-data.frame(AL$ID,apply(AL[,c(2,3)],1,mean),apply(AL[,c(4,5)],1,mean),AL[6:13]);colnames(AVG)[1:3]<-c('ID','WT','atx345')
## calculate rho
rnas<-c(2,3,4,5,2,3);K4s<-c(11,11,13,13,13,13);colnames(AL)[rnas];colnames(AL)[K4s]
rho_chen3<-return_rhos(rnas,K4s,AL);rho_chen3
rnas<-c(2,3,4,5,2,3);K4s<-c(10,10,12,12,12,12);colnames(AL)[rnas];colnames(AL)[K4s]
rho_chen2<-return_rhos(rnas,K4s,AL);rho_chen2
## plot 2d density
rnas<-c(2,3,2);K4s<-c(8,10,10);plot_2d(rnas,K4s,AVG) #wt,atx345,Chen,replace wt K4me2
rnas<-c(2,3,2);K4s<-c(9,11,11);plot_2d(rnas,K4s,AVG) #wt,atx345,Chen,replace wt K4me3



# for atxr3 and atx3-2/4/5 my data -------------------------------------------------------------------------------------------------------
RPM<-read.table("data/Figure1/ChIP1_RPM.txt",header=TRUE)
RNA_R3<-read.table("data/Figure5/PerGene_ChIP1RNA_araport.txt",header=TRUE)
sum<-apply(RNA_R3[,2:ncol(RNA_R3)],2,sum)
RNA_R3<-data.frame(RNA_R3$ID,t(t(RNA_R3[,2:ncol(RNA_R3)])*10^6/sum));colnames(RNA_R3)[1]<-"ID"
AL<-merge(RNA_R3,cDNA,by="ID");FPKM<-data.frame(AL[,1],AL[,2:(ncol(AL)-1)]*1000/AL$cDNAlength);colnames(FPKM)[1]<-"ID"
AL_1<-merge(AL,RPM,by="ID");AL<-gonly(AL_1)
AVG<-data.frame(AL$ID,apply(AL[,c(2:4)],1,mean),apply(AL[,c(5:7)],1,mean),apply(AL[,c(11:13)],1,mean),AL[,15:30]);colnames(AVG)[1:4]<-c('ID','WT','atxr3','atx345')
## calculate rho and plot 2d density H3K4me2
rnas<-c(2,3,4,11,12,13,2,3,4);colnames(AL)[rnas];K4s<-c(rep(which(colnames(AL)=="WT_me2"),3),rep(which(colnames(AL)=="x345_me2"),6)) ;colnames(AL)[K4s];
rhos<-return_rhos(rnas,K4s,AL);plot_2d(c(2,4,2),c(11,19,19),AVG)
t.test(rhos[1:3],rhos[4:6]);#WT-atx345
t.test(rhos[1:3],rhos[7:9]);#WT-atx345wt

## calculate rho and plot 2d density H3K4me3
rnas<-c(2,3,4,11,12,13,2,3,4,5,6,7,2,3,4);colnames(AL)[rnas];K4s<-c(rep(which(colnames(AL)=="WT_me3"),3),rep(which(colnames(AL)=="x345_me3"),6),rep(which(colnames(AL)=="r3_me3"),6));colnames(AL)[K4s];
rhos<-return_rhos(rnas,K4s,AL);plot_2d(c(2,3,4,2,2),c(12,8,20,20,8),AVG)
t.test(rhos[1:3],rhos[4:6]);#WT-atx345, me3
t.test(rhos[10:12],rhos[1:3]);#WT-atxr3_me3
t.test(rhos[13:15],rhos[1:3]);#WT-atxr3wt_me3


