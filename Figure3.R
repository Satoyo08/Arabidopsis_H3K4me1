# what this script do: visualize random forest models for arabidopsis; make plots in Figure.3, Supplementary Fig.5.
# For the training of the random forest, see RandomForest_arabidopsis.Rmd

source('custom_functions.r');source('RF_functions.R')
library(ROCR);library(tidyverse)

# ------------ Fig.3a-c ------------
# Random Forest models were trained and saved by RandomForest_arabidopsis.Rmd
load("data/Figure3/ATXR7_with_K36me3_ChIP8_rep_FLD") #  load randomforest models. 'ChIP8'=Fig.3, 'ChIP7'=biological replicates shown in Supplementry Fig.5
model<-ATXR7_with_K36me3_ChIP8_rep_FLD

# sort features  plot pane= 916x285
sort<-c(1,2,24,3,66,45,25,4,67,46,26,5,68,47,36,15,78,57,34,13,76,55,33,12,75,54,35,14,77,56,29,8,71,50,31,10,73,52,30,9,72,51,27,6,69,48,42,21,84,63,43,22,85,64,44,23,86,65,37,16,79,58,28,7,70,49,40,19,82,61,41,20,83,62,32,11,74,53,87,93,90,88,94,91,89,95,92,38,17,80,59,39,18,81,60)
col_key<-c(10,9,rep(c(3:6),19),rep(c(3,4,6),3),rep(c(3:6),2))
space<-c(0,2,c(2,0,0,0),rep(c(0.5,0,0,0),18),rep(c(0.5,0,0),3),rep(c(0.5,0,0,0),2))

# plot bars
par(mar=c(0.5,2,0.5,0))
IMPs<-apply(model$Importance,1,mean)[sort];sds<-apply(model$Importance,1,sd)[sort]
bp<-barplot(IMPs,col=pallets[col_key],space=space,names.arg=F,ylim=c(0,max(IMPs+max(sds))));lims<-par('usr')
plot(0,type='n',axes=FALSE,ann=FALSE,xlim=lims[1:2]+c(4,0),ylim=lims[3:4])
bp<-barplot(IMPs,col=pallets[col_key],space=space,names.arg=F,ylim=c(0,max(IMPs+max(sds))),add=T)
segments(bp,IMPs-sds,bp,IMPs+sds) # plot error bars
points(rep(bp,5),as.vector(model$Importance[sort,]),pch=16,col='gray40',cex=0.5) # plot each data points


# ---------- Fig. 3d-f ---------------
# plot ROC for the model loaded above

par(mar=c(2,2,0.5,0.5)) 
error_FPR_position=c(1:9)/10
ref<-data.frame(c((1:1000)/1000));head(ref);colnames(ref)<-'FPRref';nrow(ref);auc<-rep(0,5);TPRs<-c() #initiation

for(i in c(1:5)*2){
  pred<-prediction(model$probabilities[,i],model$groundtruth)
  perf <- performance(pred, "tpr", "fpr")
  auct<- performance(pred, "auc");auc[i/2]<-auct@y.values[[1]]
  FTEN<-data.frame(round(perf@x.values[[1]],3),round(perf@y.values[[1]],3))
  temp<-c()  # error bar 
  for(erF in error_FPR_position){
    temp<-c(temp,FTEN[which.min(abs(FTEN[,1]-erF)),2])
  }
  TPRs<-rbind(TPRs,temp)
  colnames(FTEN)<-c('FPRref',paste('tpr',i,sep='_'));ref<-merge(ref,FTEN,by='FPRref',all.x=TRUE)
  ref<-ref[!duplicated(ref$FPRref),]
}

head(ref);nrow(ref);ref<-na.omit(ref)
plot(c(0,ref$FPRref),c(0,apply(ref[,2:6],1,mean)),type='l',xlim=c(0,1),ylim=c(0,1),xlab='',ylab='')
text(0.8,0.1,paste("AUC = ",round(mean(auc),digits=3)))

sds<-apply(TPRs,2,sd);means<-apply(TPRs,2,mean)
segments(error_FPR_position,means-sds,error_FPR_position,means+sds) #add error bars
segments(error_FPR_position-0.01,means-sds,error_FPR_position+0.01,means-sds)
segments(error_FPR_position-0.01,means+sds,error_FPR_position+0.01,means+sds) 


# ----------- Fig. 3 ghi --------------
# make violin plot for selected features
# for Fig.3g
S<-TAGg[order(-TAGg$ATXR7_ChIP8TES),]
posID<-S[1:3000,1]
feature<-c("TES_CMA602","TES_CMA603");col=pallets[6]
AL[,ncol(AL)+1]<-'A'
AL[AL$ID %in% posID,ncol(AL)]<-'B'
dataset_name<-'ATXR7'
i<-1;vio<-vio_pairs_log();ggsave(paste('ggplot_violin',dataset_name,feature[i],'.pdf',sep=''),vio,width=1.2,height=2.5)
i<-2;vio<-vio_pairs_log();ggsave(paste('ggplot_violin_log',dataset_name,feature[i],'.pdf',sep=''),vio,width=1.2,height=2.5)

# for Fig.3 h
S<-TAGg[order(-TAGg$ATX1_ChIP8TSS),] # or -TAGg$ATX2_ChIP8TSS for Fig.3 i
posID<-S[1:3000,1]
feature<-c("TSS_H4K16ac","TSS_H2Bub");col=pallets[3]
AL[,ncol(AL)+1]<-'A'
AL[AL$ID %in% posID,ncol(AL)]<-'B'
dataset_name<-'ATX1'
i<-1;vio<-vio_pairs();ggsave(paste('ggplot_violin',dataset_name,feature[i],'.pdf',sep=''),vio,width=1.2,height=2.5)
i<-2;vio<-vio_pairs_log();ggsave(paste('ggplot_violin_log',dataset_name,feature[i],'.pdf',sep=''),vio,width=1.2,height=2.5)





