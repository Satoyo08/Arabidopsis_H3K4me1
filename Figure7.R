source('custom_functions.r')
library(randomForest) ;library(ROCR);library(wesanderson)

# Figure 7 b-j,SFig16 c-h Random Forest training and output ----
# these are the data sets ---
TSS_feature<-read.table("data/Figure7/TSS_feature_DF.bed",header=T)
enh_feature<-read.table("data/Figure7/enh_V65_feature_DF.bed",header=T)
set1A<-read.table('data/Figure7/Set1A_RPMs.txt',header=T) 
Mll2<-read.table('data/Figure7/Mll2_RPMs.txt',header=T) 
Mll34<-read.table('data/Figure7/Mll34_RPMs.txt',header=T) 

system('cp /Users/Satoyo/Desktop/mouse_Mll/coverage_beds/Mll34_RPMs.txt /Users/Satoyo/Desktop/論文作業/論文figs/remale_later/data/Figure7/.')
# load functions ----
return_train_CV_mm<-function(DF,posID){ 
  cv_num=750
  c<-c((ncol(DF)+1),2:ncol(DF))
  negDF<-DF[!(DF$ID %in% posID),];negDF[ncol(negDF)+1]<-as.factor(0)
  posDF<-DF[DF$ID %in% posID,];posDF[ncol(posDF)+1]<-as.factor(1)
  
  n_idx_a<-sample(1:nrow(negDF),length(posID))
  n_cv_idx<-sample(n_idx_a,cv_num)
  p_cv_idx<-sample(1:length(posID),cv_num)
  n_train_idx<-n_idx_a[!n_idx_a %in% p_cv_idx]
  p_train_idx<-c(1:length(posID))[-p_cv_idx]
  Train<-rbind(posDF[p_train_idx,c],negDF[n_train_idx,c])  
  CrVal<-rbind(posDF[p_cv_idx,c],negDF[p_cv_idx,c])  
  return(list(Train=Train,CrVal=CrVal))
}

repeat_train_mm<-function(DF,posID){
  for(i in 1:5){
    D<-return_train_CV_mm(DF,posID);print(i);print(Sys.time())
    model<-randomForest(D$Train[,2:ncol(D$Train)],y=D$Train[,1],ntree=1000) 
    if(i==1){
      imps<-model$importance
      pred<-as.numeric(predict(model,D$CrVal[,2:ncol(D$Train)]))-1
      prob<-predict(model,D$CrVal[,2:ncol(D$Train)],type='prob')
    }else{
      imps<-cbind(imps,model$importance)
      pred<-cbind(pred,as.numeric(predict(model,D$CrVal[,2:ncol(D$Train)]))-1)
      prob<-cbind(prob,predict(model,D$CrVal[,2:ncol(D$Train)],type='prob'))
    }
  }
  return(list(Importance=imps,prediction=pred,probabilities=prob,groundtruth=D$CrVal[,1]))
}

bar_colors<-wes_palette('Darjeeling2')[2:3]
plotRF<-function(model,AL,posID,outfilename,bar_colors,sort){
  P<-AL[AL$ID %in% posID,];N<-AL[-(AL$ID %in% posID),];signs<-(sign(apply(P[,2:ncol(P)],2,mean)-apply(N[,2:ncol(N)],2,mean))+1)/2+1 #if the feature is pos>neg 2, or else 1, to set the bar color
  pdf(file = outfilename,width=5,height=1.91)
  par(mar=c(0.5,2,0.5,0));IMPs<-apply(model$Importance,1,mean)[sort];sds<-apply(model$Importance,1,sd)[sort];
  bp<-barplot(IMPs,ylim=c(0,max(model$Importance)+max(sds)),col=bar_colors[signs[sort]],border='gray20');lims<-par('usr');segments(bp,IMPs-sds,bp,IMPs+sds,col='gray20')
  points(rep(bp,5),as.vector(model$Importance[sort,]),pch=16,col='gray40',cex=0.5) # plot each data points
  dev.off()
}


plot_ROC<-function(predicted,groundtruth){
  pred<-prediction(as.numeric(predicted),as.numeric(groundtruth))
  perf <- performance(pred, "tpr", "fpr")
  auct<- performance(pred, "auc");auc<-auct@y.values[[1]]
  plot(perf);text(0.95,0.1,paste("AUC = ",round(auc,digits=3)),pos=2)
  return(invisible(list(perf=perf,auc=auc)))
}

F_score<-function(predicted,groundtruth){ #Fscore, Presision,Recall. predicted calssification と ground truth classificationのベクタを受け取る。
  TBL<-table(predicted,groundtruth) #true_false cross table
  Precision<-TBL[1,1]/sum(TBL[1,])
  Recall<-TBL[1,1]/sum(TBL[,1])
  F.score<-2*Precision*Recall/(Precision+Recall)
  FPC<-data.frame(F.score,Precision,Recall)
  return(FPC)
}


# Excute training  ----
# Set1A TSS----
set1A<-read.table('data/Figure7/Set1A_RPMs.txt',header=T) 
AL<-TSS_feature[,c(4,7:ncol(TSS_feature))]
DF<-set1A;S<-DF[order(-DF$WT_ChIP_Set1A_1-DF$WT_ChIP_Set1A_2),]
posID<-S[1:3000,4];
#Set1A_RF1<-repeat_train_mm(AL,posID) 
#savefilename="Set1A_RF1"
#save(Set1A_RF1,file=savefilename)
load('data/Figure7/Set1A_RF1')
model<-Set1A_RF1
sort<-c(12,13,17,8,11,18,7,9,10,5,14,3,6,19,1,15,16,4,2)
plotRF(model,AL,posID,'Set1A_RF1_TSS.pdf',bar_colors,sort)

#plot ROC with POl2
pdf(file = 'ROC_Pol2_Set1A_TSS.pdf',width=2,height=2)
groundtruth<-rep(0,nrow(AL));groundtruth[AL$ID %in% posID]<-1
k<-grep('Pol',colnames(AL))
if(length(k)>1){print('there are multiple Pol2 columns');predicted<-'error';break}else{predicted<-AL[,k]}
par(mar=c(1,1,0.5,0.5),cex=1)
RC<-plot_ROC(predicted,groundtruth)
dev.off()



# Set1A enh ----
set1A<-read.table('data/Figure7/Set1A_enhV65_RPMs.txt',header=T) 
AL<-data.frame(paste(enh_feature$Chr,enh_feature$St,sep='_'),enh_feature[,7:ncol(enh_feature)]);colnames(AL)[1]<-"ID"
DF<-set1A;S<-DF[order(-DF$WT_ChIP_Set1A_1-DF$WT_ChIP_Set1A_2),]
posID<-paste(S[1:3000,1],S[1:3000,2],sep='_')
#Set1A_RF_enh<-repeat_train_mm(AL,posID)
#savefilename="Set1A_RF1"
#save(Set1A_RF1,file=savefilename)
load('data/Figure7/Set1A_RF_enh')
model<-Set1A_RF_enh
sort<-c(12,13,17,8,11,18,7,9,10,5,14,3,6,19,1,15,16,4,2)
plotRF(model,AL,posID,'Set1A_RF_enh.pdf',bar_colors,sort)




# Mll2 enh ----
Mll2<-read.table('data/Figure7/Mll2_RPMs.txt',header=T) 
AL<-TSS_feature[,c(4,7:ncol(TSS_feature))]
DF<-Mll2;S<-DF[order(DF$Mll2_Mll2KO-DF$Mll2_ab_CT2),]
posID<-S[1:3000,4];
#Mll2_RF_rev<-repeat_train_mm(AL,posID)
#savefilename="Mll2_RF1"
#save(Mll2_RF1,file=savefilename)
load('data/Figure7/Mll2_RF1')
model<-Mll2_RF1
sort<-c(12,13,17,8,11,18,7,9,10,5,14,3,6,19,1,15,16,4,2)
plotRF(model,AL,posID,'Mll2_RF1.pdf',bar_colors,sort)

pdf(file = 'ROC_POl2_Mll2_TSS.pdf',width=2,height=2)
groundtruth<-rep(0,nrow(AL));groundtruth[AL$ID %in% posID]<-1
k<-grep('Pol',colnames(AL))
if(length(k)>1){print('there are multiple Pol2 columns');predicted<-'error';break}else{predicted<-AL[,k]}
par(mar=c(1,1,0.5,0.5),cex=1)
RC<-plot_ROC(predicted,groundtruth)
dev.off()

# Mll2 enh----
Mll2<-read.table('data/Figure7/Mll2_enhV65_RPMs.txt',header=T) 
AL<-data.frame(paste(enh_feature$Chr,enh_feature$St,sep='_'),enh_feature[,7:ncol(enh_feature)]);colnames(AL)[1]<-"ID"
DF<-Mll2;S<-DF[order(DF$Mll2_Mll2KO-DF$Mll2_ab_CT2),]
posID<-paste(S[1:3000,1],S[1:3000,2],sep='_')
#Mll2_RF_enh<-repeat_train_mm(AL,posID)
#savefilename="Mll2_enh_RF1_2000"
#save(Mll2_RFMll2_RF_enh_2000_enh,file=savefilename)
load('data/Figure7/Mll2_enh_Rf1')
model<-Mll2_RF_enh
sort<-c(12,13,17,8,11,18,7,9,10,5,14,3,6,19,1,15,16,4,2)
plotRF(model,AL,posID,'Mll2_RF_enh.pdf',bar_colors,sort)

# ** 3.3.3 Mll34 ----
Mll34<-read.table('data/Figure7/Mll34_RPMs.txt',header=T) 
enh_feature<-read.table("data/Figure7/enh_feature_DF.bed",header=T)
AL<-data.frame(paste(enh_feature$Chr,enh_feature$St,sep='_'),enh_feature[,7:ncol(enh_feature)]);colnames(AL)[1]<-"ID"
DF<-Mll34;S<-DF[order(c(DF$Mll34_ChIP.seq_dKO_rep1+DF$Mll334_ChIP.seq_dKO_rep2)-c(DF$Mll34_ChIP.seq_WT_rep1+DF$Mll34_ChIP.seq_WT_rep1)),]
posID<-paste(S[1:3000,1],S[1:3000,2],sep='_')
#Mlll34_RF1<-repeat_train_mm(AL,posID)
#savefilename="Mlll34_RF1"
#save(Mlll34_RF1,file=savefilename)
load('data/Figure7/Mlll34_RF1')
model<-Mlll34_RF1
sort<-c(12,13,17,8,11,18,7,9,10,5,14,3,6,19,1,15,16,4,2)
plotRF(model,AL,posID,'Mlll34_RF1.pdf',bar_colors,sort)

pdf(file = 'ROC_POl2_Mll34_enh.pdf',width=2,height=2)
groundtruth<-rep(0,nrow(AL));groundtruth[AL$ID %in% posID]<-1
k<-grep('Pol',colnames(AL))
if(length(k)>1){print('there are multiple Pol2 columns');predicted<-'error';break}else{predicted<-AL[,k]}
par(mar=c(1,1,0.5,0.5),cex=1)
RC<-plot_ROC(predicted,groundtruth)
dev.off()

# Mll34_TSS ----
Mll34<-read.table('data/Figure7/Mll34_TSS_RPMs.txt',header=T) 
AL<-TSS_feature[,c(4,7:ncol(TSS_feature))]
DF<-Mll34;S<-DF[order(c(DF$Mll34_ChIP.seq_dKO_rep1+DF$Mll334_ChIP.seq_dKO_rep2)-c(DF$Mll34_ChIP.seq_WT_rep1+DF$Mll34_ChIP.seq_WT_rep1)),]
posID<-S[1:3000,4];
#Mll34_RF_TSS<-repeat_train_mm(AL,posID)
#savefilename="Mll34_RF_TSS"
#save(Mll34_RF_TSS,file=savefilename)
load('data/Figure7/Mll34_RF_TSS')
model<-Mll34_RF_TSS
sort<-c(12,13,17,8,11,18,7,9,10,5,14,3,6,19,1,15,16,4,2)
plotRF(model,AL,posID,'Mlll34_RF_TSS.pdf',bar_colors,sort)




# ROC with error bars ---------
model<-Mlll34_RF1 # input models
error_FPR_position=c(1:9)/10
par(mar=c(1,1,0.5,0.5))
ref<-data.frame(c((1:1000)/1000));head(ref);colnames(ref)<-'FPRref';nrow(ref);auc<-rep(0,5);TPRs<-c() #initiation
for(i in c(1:5)*2-1){
  pred<-prediction(model$probabilities[,i],model$groundtruth)
  perf <- performance(pred, "tpr", "fpr")
  auct<- performance(pred, "auc");auc[(i/2)+1]<-auct@y.values[[1]]
  print(auc)
  FTEN<-data.frame(round(perf@x.values[[1]],3),round(perf@y.values[[1]],3))
  temp<-c()  # error bar 
  for(erF in error_FPR_position){
    temp<-c(temp,FTEN[which.min(abs(FTEN[,1]-erF)),2])
  }
  TPRs<-rbind(TPRs,temp)
  colnames(FTEN)<-c('FPRref',paste('tpr',i,sep='_'));ref<-merge(ref,FTEN,by='FPRref',all.x=TRUE);print(nrow(ref))
  ref<-ref[!duplicated(ref$FPRref),];print(nrow(ref))
}

head(ref);nrow(ref);ref<-na.omit(ref);nrow(ref)
plot(c(0,ref$FPRref),c(0,apply(ref[,2:6],1,mean)),type='l',xlim=c(0,1),ylim=c(0,1),xlab='',ylab='')
text(0.95,0.1,paste("AUC = ",round(mean(auc),digits=3)),pos=2)
sds<-apply(TPRs,2,sd);means<-apply(TPRs,2,mean)
segments(error_FPR_position,means-sds,error_FPR_position,means+sds) #add error bars
segments(error_FPR_position-0.01,means-sds,error_FPR_position+0.01,means-sds)
segments(error_FPR_position-0.01,means+sds,error_FPR_position+0.01,means+sds) 
dev.off()
