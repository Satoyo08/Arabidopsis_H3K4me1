# ------------ Supplementary Fig3.a,b ------------
# plotted using ngs.plot.r, similary to Fig.1c

# ------------ Supplementary Fig3.c ------------
source('custom_functions.r');source('RF_functions.R')

# ATX(R)s bound genes
S<-TAGg[order(-TAGg$ATX1_ChIP7TSS),];ATX1_bound_CHIP7<-S[1:3000,1];
S<-TAGg[order(-TAGg$ATX2_ChIP7TSS),];ATX2_bound_CHIP7<-S[1:3000,1];
S<-TAGg[order(-TAGg$ATXR7_soiTES),];ATXR7_bound_CHIP7<-S[1:3000,1];
S<-TAGg[order(-TAGg$ATX1_ChIP8TSS),];ATX1_bound<-S[1:3000,1];
S<-TAGg[order(-TAGg$ATX2_ChIP8TSS),];ATX2_bound<-S[1:3000,1];
S<-TAGg[order(-TAGg$ATXR7_ChIP8TES),];ATXR7_bound<-S[1:3000,1];

write.table(ATX1_bound_CHIP7,file='data/list_of_genes/ATX1_bound_genes_rep2.txt',col.names = F,row.names=FALSE, quote=FALSE,sep = "\t")
write.table(ATX2_bound_CHIP7,file='data/list_of_genes/ATX2_bound_genes_rep2.txt',col.names = F,row.names=FALSE, quote=FALSE,sep = "\t")
write.table(ATXR7_bound_CHIP7,file='data/list_of_genes/ATXR7_bound_genes_rep2.txt',col.names = F,row.names=FALSE, quote=FALSE,sep = "\t")
write.table(ATX1_bound,file='data/list_of_genes/ATX1_bound_genes_rep1.txt',col.names = F,row.names=FALSE, quote=FALSE,sep = "\t")
write.table(ATX2_bound,file='data/list_of_genes/ATX2_bound_genes_rep1.txt',col.names = F,row.names=FALSE, quote=FALSE,sep = "\t")
write.table(ATXR7_bound,file='data/list_of_genes/ATXR7_bound_genes_rep1.txt',col.names = F,row.names=FALSE, quote=FALSE,sep = "\t")


# atx12r7-marked genes
ALm<-merge(atx12r7_rpm,len,by='ID');ALrpkm<-data.frame(ALm$ID,ALm[,2:17]*1000/ALm$len);colnames(ALrpkm)[1]<-'ID'
DF<-gonly(ALrpkm)
S<-DF[order(DF$x12r7_H3K4me1-DF$WT_H3K4me1),];atx12r7_r2<-S[1:3000,1]
the_list<-list(new_red,atx12r7_r2,ATX1_bound,ATX1_bound_CHIP7,ATX2_bound,ATX2_bound_CHIP7,ATXR7_bound,ATXR7_bound_CHIP7)

HYPG<-as.data.frame(matrix(nrow=length(the_list),ncol=length(the_list)));INT<-HYPG;enr<-HYPG;K<-HYPG
for(i in 1:length(the_list)){
  K4_ids<-the_list[[i]]
  for(j in 1:length(the_list)){
    tag_ids<-the_list[[j]]
    q<-length(intersect(tag_ids,K4_ids)) 
    m<-length(K4_ids)
    n<-27408-m
    k<-length(tag_ids)
    e<-k*m/(n+m) #Expected
    q/e #enrichment
    pv<-phyper(q, m, n, k,lower.tail = FALSE,log.p = FALSE)
    HYPG[i,j]<-pv;INT[i,j]<-q;enr[i,j]<-q/e;K[i,j]<-m
  }
}
HYPGl<--log10(HYPG+0.00001)
enr[enr>9]<-0 # remove diagonal
LL<-length(the_list)
par(mar=c(0.5,0.5,0.5,0.5))
colf<-colorRampPalette(c(pallets[5],pallets[1]))
plot(1,1,type='n',ylim=c(0,9),xlim=c(0,9))
for(i in 1:8){
  points(1:LL,rep(9-i,8),cex=as.numeric(enr[,i])*0.8,col=colf(1002)[floor(HYPGl[,i]*1000/max(HYPGl))+2],pch=16)
}

#legends
plot(1,1,type='n',ylim=c(0,17),xlim=c(0,17))
points(1:10,rep(17-i,10),cex=c(0:10)*0.8,col='gray',pch=16)
