
TAGg<-read.table('data/RF_for_arabisopsis/TAGg.txt')　#TAG-Input. the higher the 
TAG_R3<-read.table('data/RF_for_arabisopsis/TAGg_R3.txt')
AL<-read.table('data/RF_for_arabisopsis/AL.txt')
AL_R<-read.table('data/RF_for_arabisopsis/AL_R.txt')

getChr5<-function(AtID){
  C1to4<-AtID[c(grep("AT5G",AtID,invert=T))]
  C5<-AtID[c(grep("AT5G",AtID))]
  C<-list(C1to4,C5)
  names(C)<-c("C1to4","C5")
  return(C)
}


# 説明変数のDF,見分けたいpositive遺伝子のAtIDリストposID,DFから除きたい列名に部分マッチする文字列omitsを与えると、train dataとCV data を返す関数 除きたい列がなければ列名に存在しないを入れとく
return_train_CV<-function(DF,posID,omits){ 
  DATA<-DF[,c(grep(omits,colnames(DF),invert=T))]
  negID<-DATA[!(DATA$ID %in% posID),]$ID
  
  posIDlist=getChr5(posID)
  negIDlist=getChr5(negID)
  
  Tr_neg<-sample(negIDlist[[1]],length(posIDlist[[1]]))
  CV_neg<-sample(negIDlist[[2]],length(posIDlist[[2]]))
  
  Train<-rbind(DATA[DATA$ID %in% Tr_neg,],DATA[DATA$ID %in% posIDlist[[1]],])
  CrVal<-rbind(DATA[DATA$ID %in% CV_neg,],DATA[DATA$ID %in% posIDlist[[2]],])
  Train[,1]<-as.factor(c(rep(0,nrow(DATA[DATA$ID %in% Tr_neg,])),rep(1,nrow(DATA[DATA$ID %in% posIDlist[[1]],]))))
  CrVal[,1]<-as.factor(c(rep(0,nrow(DATA[DATA$ID %in% CV_neg,])),rep(1,nrow(DATA[DATA$ID %in% posIDlist[[2]],]))))
  return(list(Train=Train,CrVal=CrVal))
}

F_score<-function(predicted,groundtruth){ #Fscore, Presision,Recall. predicted calssification と ground truth classificationのベクタを受け取る。
  TBL<-table(predicted,groundtruth) #true_false cross table
  Precision<-TBL[1,1]/sum(TBL[1,])
  Recall<-TBL[1,1]/sum(TBL[,1])
  F.score<-2*Precision*Recall/(Precision+Recall)
  FPC<-data.frame(F.score,Precision,Recall)
  return(FPC)
}

repeat_train<-function(DF,posID,omits){
  for(i in 1:5){
    D<-return_train_CV(DF,posID,omits);print(i);print(Sys.time())
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


repeat_train_2<-function(DF,posID,omits){
  mods<-list()
  for(i in 1:5){
    D<-return_train_CV(DF,posID,omits);print(i);print(Sys.time())
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
    mods<-c(mods,list(model))
  }
  return(list(Importance=imps,prediction=pred,probabilities=prob,groundtruth=D$CrVal[,1],models=mods))
}
