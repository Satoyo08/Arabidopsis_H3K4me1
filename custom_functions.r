setwd('/Users/Satoyo/Desktop/論文作業/論文figs/remale_later')
# colors 
library(RColorBrewer)
pallets<-brewer.pal(11, "Spectral")
makeTransparent<-function(someColor, alpha=100){
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
rgb_pal<-c(makeTransparent('gray80',150),'#FF5E5B','#FFED66','#FFA661','#00CECB','#809693','#80DE99','#C09E7A') 


# dataframe preps
add_mean<-function(W){W[,7]<-apply(W[,2:6],1,mean);return(W)}

gonly<-function(X){ # for a given dataframe, match AtID with araport 11 annotation and return a subdataframe containing "protein coding genes"
  arapo<-read.table("data/refs/araport11_all_sorted.bed",header=F) #V6; 1=pc,2=pseudo,3=TE,4=noncoding,5=noveltranscrived
  colnames(arapo)<-c("Chrr","st","en","ID","d","type","V7")
  X1<-merge(X,arapo,by="ID")
  XX<-(X1[!duplicated(X1$ID),])
  X1<-subset(XX,XX$Chrr!="C"&XX$Chrr!="M"&XX$type==1) 
  return(X1)
}


labelDF<-function(DATA,cutoff){
  DATA[,1][DATA[,1]>=cutoff]<-0 ;DATA[,1][DATA[,1]< cutoff]<-1
  print(table(DATA[,1]))
  DATA[,1]<-as.factor(DATA[,1])
  return(DATA)
}

# makind graphs
plot_shu<-function(x,y,col=rgb(0,0,0,0.3),pch=16,lwd=1,type='p',limmaxx=0.98,limmaxy=0.98,limminx=0,limminy=0){#quantile 0.98で軸を絞ってplot
  plot(x,y,xlim=c(quantile(x,limminx),quantile(x,limmaxx)),ylim=c(quantile(y,limminy),quantile(y,limmaxy)),pch=pch,col=col,xlab="",ylab="",lwd=lwd,type=type)
}



#---------for SVM qstring----------------------
# list of DNAstrings DNAS, RC_DNAS required as global variables
eq_revcon<-function(i,j,cost_rc,cost){
  if(DNAS[[i]]==RC_DNAS[[j]]){val<-cost_rc}else{val<-cost}
  return(val)
}

offset1_hits_dir<-function(i,j,cost_off_same,cost_off_rc,cost_nomatch){
  ofmn1<-function(DNAStr){subseq(DNAStr,1,length(DNAStr)-1)}
  ofmp1<-function(DNAStr){subseq(DNAStr,2,length(DNAStr))}
  if(ofmn1(DNAS[[i]])==ofmp1(DNAS[[j]])){
    val<-cost_off_same
  }else if(ofmp1(DNAS[[i]])==ofmn1(DNAS[[j]])){
    val<-cost_off_same
  }else if(ofmn1(DNAS[[i]])==ofmp1(RC_DNAS[[j]])){
    val<-cost_off_rc
  }else if(ofmp1(DNAS[[i]])==ofmn1(RC_DNAS[[j]])){
    val<-cost_off_rc
  }else{
    val<-cost_nomatch
  }
  return(val)
}

my_distance_dir<-function(i,j,cost_rc,cost_off_same,cost_off_rc,cost_nomatch){
  if(i==j){val=0}else{
    val=eq_revcon(i,j,cost_rc,offset1_hits_dir(i,j,cost_off_same,cost_off_rc,cost_nomatch))
  }
  return(val)
}

kmer_list_to_matrix_dir<-function(feature_list){
  DM<-data.frame();k=0
  for(i in 1:length(feature_list)){
    for(j in i:length(feature_list)){
      dist<-my_distance_dir(i,j,1,2,9,10) #距離の定義部。i,j,cost_off_same,cost_off_rc,cost_nomatch
      DM[i,j]<-dist
      if(i<j){DM[j,i]<-dist}
    }
  }
  
  X<-DM
  dist_mi <- 1/as.dist(X) # one over, as qgraph takes similarity matrices as input
  dist_mi[dist_mi==0.1]<-0
  return(list(dist_mi,X))  
}

replace_string<-function(x,dict){
  for(i in 1:nrow(dict)){
    x<-replace(x,which(x==dict[i,1]),as.character(dict[i,2]))
  }
  return(x)
}

# violin plot
vio_pairs_log<-function(){
  k<-which(colnames(AL)==feature[i])
  DF<-data.frame(log(AL[,k]+0.5),AL[,ncol(AL)])
  colnames(DF)<-c('data','negpos')
  ggplot(data=DF)+
    aes(x=negpos,y=data,fill=negpos)+
    geom_violin(bw="nrd0",color='gray70')+
    theme_classic()+
    scale_fill_manual(values=c('gray70',col))+
    theme(panel.background = element_rect(fill = "white",color='black'),axis.text.y=element_text(angle=90,hjust=0.5),
          axis.text.x=element_blank(),axis.title=element_blank(),legend.position = 'none',axis.line = element_blank()
    )+
    geom_boxplot(width=0.1,color='gray50',fill='black',outlier.size=0.2,outlier.color = 'gray50')
}
vio_pairs<-function(){
  k<-which(colnames(AL)==feature[i])
  DF<-data.frame(AL[,k],AL[,ncol(AL)])
  colnames(DF)<-c('data','negpos')
  ggplot(data=DF)+
    aes(x=negpos,y=data,fill=negpos)+
    geom_violin(bw="nrd0",color='gray70')+
    theme_classic()+
    scale_fill_manual(values=c('gray70',col))+
    theme(panel.background = element_rect(fill = "white",color='black'),axis.text.y=element_text(angle=90,hjust=0.5),
          axis.text.x=element_blank(),axis.title=element_blank(),legend.position = 'none',axis.line = element_blank()
    )+
    geom_boxplot(width=0.1,color='gray50',fill='black',outlier.size=0.2,outlier.color = 'gray50')
}

plot_ROC_fromTPRFPR<-function(filename){
  options(scipen=2)
  aucscore=read.table(filename,nrow=1);FPRTPR=read.table(filename,skip=1,header=T)
  error_FPR_position=c(1:9)/10
  TPRs<-c();for(erF in error_FPR_position){
    temp<-erF
    for(r in c(1,3,5,7,9)){temp<-c(temp,FPRTPR[which.min(abs(FPRTPR[,r]-erF)),r+1])}
    TPRs<-rbind(TPRs,temp)
  }
  sds<-apply(TPRs[,2:6],1,sd);means<-apply(TPRs[,2:6],1,mean)
  ref<-data.frame(c((1:100000)/100000));head(ref);colnames(ref)<-'FPRref';nrow(ref)
  for(i in 1:5){
    FTEN<-FPRTPR[,c(2*i-1,2*i)];colnames(FTEN)<-c('FPRref',paste('tpr',i));ref<-merge(ref,FTEN,by='FPRref',all.x=TRUE)
  }
  ref<-na.omit(ref);ref<-ref[!duplicated(ref$FPRref),] # cleanup. this step lose data, but its effect on apperance for ROC plot seems to be ignorable.
  ref[ncol(ref)+1]<-apply(ref[,2:ncol(ref)],1,mean)
  par(mar=c(2,2,0.1,0.1))
  plot(ref[,1],ref[,7],type='l') # draw ROC
  segments(error_FPR_position,means-sds,error_FPR_position,means+sds) #add error bars
  segments(error_FPR_position-0.01,means-sds,error_FPR_position+0.01,means-sds)
  segments(error_FPR_position-0.01,means+sds,error_FPR_position+0.01,means+sds)
  text(0.8,0.1,aucscore$V1)
  return(invisible(list(DF=ref,auc_score=aucscore,errorbars=TPRs)))
}

# common variables and data
## definitions of ATX12R7-marked, ATX345-marked, ATXR3-marked genes
RPKM<-read.table("data/Figure1/ChIP1_RPKM.txt",header=TRUE)
RPM<-read.table("data/Figure1/ChIP1_RPM.txt",header=TRUE)
DF<-gonly(RPKM);S<-DF[order(DF$x127_me1-DF$WT_me1),];new_red<-S[1:3000,1]
DF<-gonly(RPM);S<-DF[order(DF$x345_me2-DF$WT_me2),];new_green<-S[1:3000,1]
DF<-gonly(RPM);S<-DF[order(DF$r3_me3-DF$WT_me3),];new_blue<-S[1:3000,1]

# araport11 annotation table
arapo<-read.table("data/refs/araport11_all_sorted.bed",header=F);colnames(arapo)<-c("Chr","st","en","ID","d","X6","V7");len<-data.frame(arapo$ID,arapo$en-arapo$st);colnames(len)<-c("ID","len")
# rep2 of H3K4me-atx12r7, atx12, atxr7, atx1, atx2 series
atx12r7_rpm=read.table('data/refs/atx12r7_me3H3_RPM.txt',header=T)

