
source('custom_functions.r');source("RF_functions.R");library(randomForest)
# ------- in silico hub simulation ------------
# --- training of WT ATX2 model
S<-TAGg[order(-TAGg$ATX2_ChIP8TSS),];posID<-S[1:3000,1];ATX2_bound<-posID
ATX2_ChIP8TSS_fm<-repeat_train_2(AL,posID,'remove nothing') # This model was trained in the same way as Figure.3c('ATX2_with_K36me3_ChIP8_rep), and behave very similary to Figure.3c. The reason we use this model here instead of ATX2_with_K36me3_ChIP8_rep model is because we missed saving 'model' object for ATX2_with_K36me3_ChIP8_rep.  
#save(ATX2_ChIP8TSS_fm2,file='ATX2_ChIP8TSS_fm2')
load('data/supplemental/ATX2_ChIP8TSS_fm2')  
model<-ATX2_ChIP8TSS_fm2
AL_R<-AL;AL0<-AL_R
AL0$TES_H2Bub<-0;AL0$TSS_H2Bub<-0;AL0$mid_gb_k_H2Bub<-0;AL0$middle_of_genebody_H2Bub<-0 #replacing the values of H2Bub features to 0, thereby simulating hub

# get list of the genethat lost ATX2 in in silico hub
cands<-list()
for(i in 1:5){
  pAL<-predict(model$models[[i]],AL_R[,2:ncol(AL_R)]);table(pAL)
  p0AL<-predict(model$models[[i]],AL0[,2:ncol(AL_R)]);table(p0AL)
  cands<-c(cands,list(intersect(posID[pAL==1],posID[p0AL==0]))) # the genes, according to the model's prediction, to which ATX will not bind in the absence of H2Bub
}
cand2<-as.data.frame(table(unlist(cands)))[(as.data.frame(table(unlist(cands)))$Freq)>0,1];length(cand2) # this list is not always identical and a little vary time to time, depending on the prediction above (thus may deviate from n=265 in the figure)

# --- training of WT ATX1 model
S<-TAGg[order(-TAGg$ATX1_ChIP8TSS),];posID<-S[1:3000,1];ATX1_bound<-posID
ATX1_ChIP8TSS_fm<-repeat_train_2(AL,posID,'remove nothing') # This model was trained in the same way as Figure.3b('ATX1_with_K36me3_ChIP8_rep), and behave very similary to Figure.3c. The reason we use this model here instead of ATX1_with_K36me3_ChIP8_rep model is because we missed saving 'model' object for ATX2_with_K36me3_ChIP8_rep.  
#save(ATX1_ChIP8TSS_fm,file='ATX1_ChIP8TSS_fm2')
load('data/supplemental/ATX1_ChIP8TSS_fm2')  
model<-ATX1_ChIP8TSS_fm2
AL_R<-AL;AL0<-AL_R
AL0$TES_H2Bub<-0;AL0$TSS_H2Bub<-0;AL0$mid_gb_k_H2Bub<-0;AL0$middle_of_genebody_H2Bub<-0 

# get list of the genethat lost ATX1 in in silico hub
cands<-list()
for(i in 1:5){
  pAL<-predict(model$models[[i]],AL_R[,2:ncol(AL_R)]);table(pAL)
  p0AL<-predict(model$models[[i]],AL0[,2:ncol(AL_R)]);table(p0AL)
  cands<-c(cands,list(intersect(posID[pAL==1],posID[p0AL==0]))) 
}
cand1<-as.data.frame(table(unlist(cands)))[(as.data.frame(table(unlist(cands)))$Freq)>0,1] # this list is not always identical and a little vary time to time, dep ending on the prediction above
cand12<-union(cand1,cand2) # genes that lost either ATX1 or ATX2 in in silico hub

#---- SFig 7b -----
ChIP10_FLG<-gonly(read.table('data/supplemental/ChIP10_RPM_TSS.txt',header=T))
ChIP11_FLG<-gonly(read.table('data/supplemental/ChIP11_RPM_TSS.txt',header=T))


pChIP10<-ChIP10_FLG[ChIP10_FLG$ID %in% posID,]
pChIP11<-ChIP11_FLG[ChIP11_FLG$ID %in% posID,]
c2ChIP10<-ChIP10_FLG[ChIP10_FLG$ID %in% cand2,]
c2ChIP11<-ChIP11_FLG[ChIP11_FLG$ID %in% cand2,]

boxplot(ChIP10_FLG$ATX2_tag,ChIP10_FLG$ATX2_tag_hub,
        ChIP11_FLG$ATX2,ChIP11_FLG$ATX2_hub,
        c2ChIP10$ATX2_tag,c2ChIP10$ATX2_tag_hub,
        c2ChIP11$ATX2,c2ChIP11$ATX2_hub,
        outline = F,notch=T,col=c(rep(c('gray70','darkolivegreen3'),2),rep(c('gray90','khaki2'),2)))


# test for the difference 
wilcox.test(ChIP10_FLG$ATX2_tag/ChIP10_FLG$ATX2_tag_hub,c2ChIP10$ATX2_tag/c2ChIP10$ATX2_tag_hub)
wilcox.test(ChIP11_FLG$ATX2/ChIP11_FLG$ATX2_hub,c2ChIP11$ATX2/c2ChIP11$ATX2_hub)

# calculate the overlaps
K4_files<-list.files('/Users/Satoyo/data/ChIP11/hubs/',pattern='intersect');K4_files_r<-K4_files[grep('reduction',K4_files)]
tag_files<-list.files('/Users/Satoyo/data/ChIP11/FLAG/',pattern='intersect')

setwd('/Users/Satoyo/data/ChIP11/hubs/')
intersect_files<-K4_files_r
HYPG<-as.data.frame(matrix(nrow=length(intersect_files),ncol=7))
for(i in 1:length(intersect_files)){
  ids<-unique(read.table(as.character(intersect_files[i]))$V4)
  q<-length(intersect(ids,cand)) #number of white balls drawn without replacement from an urn which contains both black and white balls
  m<-length(cand)# the number of white balls in the urn.
  n<-27408-m# the number of black balls in the urn.
  k<-length(ids) # Yの数　k	the number of balls drawn from the urn.
  e<-k*m/(n+m)#Expected
  q/e #enrichment
  pv<-phyper(q, m, n, k,lower.tail = FALSE,log.p = FALSE) # probability lower.tail=TRUE だと、今の引きよりも少なくひく確率
  HYPG[i,]<-c(as.character(intersect_files[i]),k,m,q,e,q/e,pv)
}
HYPG

# overlaps between Hypo-H3K4me1 regions in hub and ATX2-lost regions in hub; SFig7d, two left column 
HYPG<-as.data.frame(matrix(nrow=length(K4_files),ncol=length(tag_files)+1));INT<-HYPG;enr<-HYPG;K<-HYPG;KK<-HYPG
setwd('/Users/Satoyo/data/ChIP11/hubs/')
for(i in 1:length(K4_files)){
  K4_ids<-unique(read.table(as.character(K4_files[i]))$V4)
  print(length(K4_ids));print(K4_files[i])
  HYPG[i,1]<-as.character(K4_files[i])
  for(j in 1:length(tag_files)){
    tag_ids<-unique(read.table(paste('/Users/Satoyo/data/ChIP11/FLAG/',as.character(tag_files[j]),sep=''))$V4)
    q<-length(intersect(tag_ids,K4_ids))
    m<-length(K4_ids)
    n<-27408-m
    k<-length(tag_ids) 
    e<-k*m/(n+m)#Expected
    q/e #enrichment
    pv<-phyper(q, m, n, k,lower.tail = FALSE,log.p = FALSE) 
    HYPG[i,j+1]<-pv;INT[i,j+1]<-q;enr[i,j+1]<-q/e;K[i,j+1]<-m;KK[i,j+1]<-k
  }
}
# p-val, enrichment, intersect
HYPG[c(4,10,6,12,2,8),c(1,5,9)]
enr[c(4,10,6,12,2,8),c(1,5,9)]
INT[c(4,10,6,12,2,8),c(1,5,9)]
KK[c(4,10,6,12,2,8),c(1,5,9)]


#scores for bound,unbound

#intersect_files<-tag_files
#par(mfrow=c(1,length(intersect_files)+1),mar=c(0.5,2,0.5,0))
# 
# NRM_prob[,7]<-apply(NRM_prob[,2:6],1,mean);UB0_prob[,7]<-apply(UB0_prob[,2:6],1,mean)
# NRM_pos<-NRM_prob[AL_R$ID %in% posID,]$V7
# UB0_pos<-UB0_prob[AL_R$ID %in% posID,]$V7
# boxplot(NRM_pos,UB0_pos,fill=c('gray70','olivedrab3'),ylim=c(0,1))
# for(i in 1:length(intersect_files)){
#   ids<-unique(read.table(as.character(intersect_files[i]))$V4)
#   NRM_goi<-NRM_prob[AL_R$ID %in% ids,]$V7
#   UB0_goi<-UB0_prob[AL_R$ID %in% ids,]$V7
#   boxplot(NRM_goi,UB0_goi,fill=c('gray70','olivedrab3'),ylim=c(0,1))
#   #boxplot(NRM_pos-UB0_pos,NRM_goi-UB0_goi,fill=c('gray70','olivedrab3'))
# }
# 
# ratio<-list()
# ratio<-c(ratio,list(UB0_pos/(NRM_pos+0.00001)))
# for(i in 1:length(intersect_files)){
#   ids<-unique(read.table(as.character(intersect_files[i]))$V4)
#   NRM_goi<-NRM_prob[AL_R$ID %in% ids,]$V7
#   UB0_goi<-UB0_prob[AL_R$ID %in% ids,]$V7
#   ratio<-c(ratio,list(UB0_goi/(NRM_goi+0.00001)))
# }
# dev.off()
# boxplot(ratio)

### bubble_plot_hub
par(mar=c(0.5,0.5,0.5,0.5))
colf<-colorRampPalette(rev(c(pallets[1],pallets[5],'gray80')))
pval<- -log10(p_val)

# hypo_H3K4me1 in hubs and ATX2 rep1,rep2, either ATX1 or ATX2 in silico
# calculated by L54 
enr<-c(1.41,1.40,1.43,1.45,0.52,0.52,1.73,1.87,1.83,1.67,0.00,0.00,1.73,1.83,2.03,2.61,0.00,0.00)
p_val<-c(8.7E-11,2.9E-13,1.8E-13,1.1E-10,0.61,0.61,4.3E-04,5.1E-06,2.7E-05,2.5E-03,0.20,0.20,4.7E-03,5.5E-04,9.8E-05,2.3E-06,0.12,0.12)
intersect_num<-c(279,359,331,235,1,1,38,53,47,30,0,0,22,30,30,27,0,0)
yax<-rep(c(1:6),3)
xax<-rep(1:3,each=6)
pval<- -log10(p_val)
pval[pval>5]<-5
plot(1,1,type='n',ylim=c(0,max(yax)+1),xlim=c(0,max(xax)+1),xaxt='n',yaxt='n')
points(xax,rev(yax),cex=as.numeric(enr)*1.3,col=colf(1002)[floor(pval*1000/5)+2],pch=16)
text(xax,rev(yax),intersect_num,cex=0.6,col='black')
# genes that lost ATX2 in planta-in silico
enr<-c(1.53,2.13)
p_val<-c(8.81E-06,3.54E-03)
pval<- -log10(p_val);pval[pval>5]<-5
intersect_num<-c(84,13)
yax<-c(1,2)
xax<-c(1,1)

# hypo-H2Bub genes and ATX12 bound or ATX12r7-marked
enr<-c(1.30,2.23,3.05,6.06)
p_val<-c(0.19,4.0E-28,3.6E-03,4.2E-208);pval<- -log10(p_val);pval[pval>5]<-5
intersect_num<-c(3,181,5,348)
xax<-rep(c(1,2),each=2)

# hypo-H4K16ac genes
enr<-c(0.99,0.32,0.92,0.70)
p_val<-c(0.70,1.00,0.67,0.42);pval<- -log10(p_val);pval[pval>5]<-5
intersect_num<-c(43,2,40,1)




#legends
plot(1,1,type='n',ylim=c(0,6),xlim=c(0,3))
points(rep(1,6),6:1,cex=c(0.25,0.5,1,2,4,8)*1.3,col='gray',pch=16)
points(rep(2,6),1:6,cex=3*0.8,col=colf(1002)[-log10(c(0.00001,0.0001,0.001,0.01,0.1,1))*1000/5+2],pch=16)
text(rep(1,6),6:1,c('0.25','0.5','1','2','4','8'),cex=0.6)
text(rep(2,6)+0.5,1:6,c('1e05','1e04','1e03','0.01','0.1','1'),cex=0.6)







# make bubble plots
bubblle_plots<-function(sizeDF,colorDF,scaler=1)
  
  par(mar=c(0.5,0.5,0.5,0.5))
colf<-colorRampPalette(c(pallets[5],pallets[1]))
plot(1,1,type='n',ylim=c(0,17),xlim=c(0,17))
for(i in 1:16){
  points(1:LL,rep(17-i,16),cex=as.numeric(enr[,i])*0.8,col=colf(1002)[floor(HYPGl[,i]*1000/max(HYPGl))+2],pch=16)
}

plot(1,1,type='n',ylim=c(0,17),xlim=c(0,17))
points(1:10,rep(17-i,10),cex=c(0:10)*0.8,col='gray',pch=16)
points(1:10,rep(17-3,10),cex=3*0.8,col=colf(1002)[-log10(c(0.00001,0.0001,0.001,0.01,0.1,1))*1000/5+2],pch=16)





#--------ashh_etc----------
setwd('/Users/Satoyo/data/ChIP11/ashh_etc/intersect')

# 重なるか？ATX1,2 bound genes, atx12r7-marked genes, ATX12-marked genes とH2Bub,H4K16ac, H3K36me3 が低下したgenes
DF<-gonly(ChIP2_RPKM);S<-DF[order(DF$x12_me1-DF$WT_me1),];x12_marked<-S[1:1500,1]

sdgK4<-read.table('/Users/Satoyo/data/ChIP3/sdg8_k4me1_nomodel_peaks.narrowPeak_intersect')
intersects<-list.files('/Users/Satoyo/data/ChIP11/ashh_etc/intersect',pattern='atx12r7')
intersect_files<-intersects[grep('reduction',intersects)];intersect_files
target<-sdgK4$V4
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
cpy(HYPG)





#重なるか？ATX1,2 bound genes, atx12r7-marked genes, atx12r7でK36me3さがったgenes とashh targets
intersects<-list.files('/Users/Satoyo/data/ChIP11/ashh_etc/intersect',pattern='ash');setwd('/Users/Satoyo/data/ChIP11/ashh_etc/intersect')
intersect_files<-intersects[grep('reduction',intersects)];intersect_files

DF<-gonly(ChIP2_RPKM);S<-DF[order(DF$x12_me1-DF$WT_me1),];x12_marked<-S[1:3000,1]
target<-intersect(ATX1_bound,ATX2_bound)
target<-x12_marked
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
cpy(HYPG)

# ashhs-bubble
ds<-ctrlV()
rep1<-ds[c(1,3,5,7,9),]
rep2<-ds[c(1,3,5,7,9)+1,]
xax<-rep(c(1,2),5)
yax<-c(1,1,2,2,3,3,4,4,5,5)

par(mar=c(0.5,0.5,0.5,0.5))
colf<-colorRampPalette(rev(c(pallets[1],pallets[5],'gray80')))
pval<- -log10(ds$V7)
pval[pval>5]<-5
plot(1,1,type='n',ylim=c(0,6),xlim=c(0,3),xaxt='n',yaxt='n')
points(xax,rev(yax),cex=as.numeric(ds$V6)*1.3,col=colf(1002)[floor(pval*1000/5)+2],pch=16)
text(xax,rev(yax),ds$V2,cex=0.6)
text(xax,rev(yax),paste('(',ds$V4,')'),cex=0.6)

plot(1,1,type='n',ylim=c(0,6),xlim=c(0,3))
points(rep(1,6),6:1,cex=c(0.25,0.5,1,2,4,8)*1.3,col='gray',pch=16)
points(rep(2,6),1:6,cex=3*0.8,col=colf(1002)[-log10(c(0.00001,0.0001,0.001,0.01,0.1,1))*1000/5+2],pch=16)
text(rep(1,6),6:1,c('0.25','0.5','1','2','4','8'))
text(rep(2,6)+0.5,1:6,c('1e05','1e04','1e03','0.01','0.1','1'))


plot_shu(ALm$ashh4_H3,ALm$len)
plot_shu(ALrpkm$ashh4_H3,ALm$len)
head(ALm)

ChIP11_K36_rpm<-read.table('/Users/Satoyo/data/ChIP11/ashh_etc/intersect/ChIP11_ashhetc_RPM.txt',header=T)
ALm<-merge(ChIP11_K36_rpm,len,by='ID');ALrpkm<-data.frame(ALm$ID,ALm[,2:(ncol(ALm)-1)]*1000/ALm$len);colnames(ALrpkm)[1]<-'ID'
write.table(ALrpkm,file='/Users/Satoyo/data/ChIP11/ashh_etc/intersect/ChIP11_ashhetc_RPKM.txt',row.names = F,quote=F)
