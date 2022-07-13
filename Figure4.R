# what this script do: visualize lSVM models. make plots in Figure.4, Supplementary Fig 9,10,11,18
# make Supplementary Table.1
# For the training of the SVM models, see SVM_training.ipyenv


source('custom_functions.r',echo=FALSE)
library(Biostrings);library(qgraph);library(eulerr)

# ------------ Fig.4 a-c, SFig.9ab, SFig.18 ; ROC plots  ------------
# plot ROC for SVM models. SVM models were trained and saved by SVM.py
ROC_6mer<-plot_ROC_fromTPRFPR('data/Figure4/ROC/ChIP8_ATX1_re_6mer_FPRTPR') # loadFPR/TPR data sets of interest. 'ChIP8'=Fig.4, 'ChIP7'=biological replicates shown in Supplementry Fig.4. 
ROC_hi<-plot_ROC_fromTPRFPR('data/Figure4/ROC/ATX1_ChIP8_highmer_FPRTPR') # 'selected'. models trained only with highly weighted 120 6-mers 
error_x<-rep(c(1:9)/10,3)
par(mar=c(2,2,0.2,0.2))
RC<-ROC_6mer;sds<-apply(RC$errorbars[,2:6],1,sd);means<-apply(RC$errorbars[,2:6],1,mean)
plot(RC$DF[,1],RC$DF[,7],type='l',lty="solid",xlim=c(0,1),ylim=c(0,1),labels=F);segments(rep(RC$errorbars[,1],3)+c(rep(0,9),rep(-0.01,18)),c(means-sds,means-sds,means+sds),rep(RC$errorbars[,1],3)+c(rep(0,9),rep(0.01,18)),c(means+sds,means-sds,means+sds))
text(1,0.25,RC$auc_score,pos=2);par(new=T)
RC<-ROC_hi;sds<-apply(RC$errorbars[,2:6],1,sd);means<-apply(RC$errorbars[,2:6],1,mean)
plot(RC$DF[,1],RC$DF[,7],type='l',lty="dashed",xlim=c(0,1),ylim=c(0,1),labels=F);segments(rep(RC$errorbars[,1],3)+c(rep(0,9),rep(-0.01,18)),c(means-sds,means-sds,means+sds),rep(RC$errorbars[,1],3)+c(rep(0,9),rep(0.01,18)),c(means+sds,means-sds,means+sds))
text(1,0.3,RC$auc_score,pos=2)


# Supplementary Fig.18
ROC_5mer<-plot_ROC_fromTPRFPR('data/Figure4/ROC/ChIP8_ATX1_re_5mer_FPRTPR')
ROC_7mer<-plot_ROC_fromTPRFPR('data/Figure4/ROC/ChIP8_ATX1_7mer_FPRTPR')
ROC_kern<-plot_ROC_fromTPRFPR('data/Figure4/ROC/ChIP8_ATX1_kernel_FPRTPR')
ROC_6mer<-plot_ROC_fromTPRFPR(xc)
error_x<-rep(c(1:9)/10,3)
par(mar=c(2,2,0.2,0.2))
RC<-ROC_6mer;sds<-apply(RC$errorbars[,2:6],1,sd);means<-apply(RC$errorbars[,2:6],1,mean)
plot(RC$DF[,1],RC$DF[,7],type='l',lty="solid",xlim=c(0,1),ylim=c(0,1),labels=F);segments(rep(RC$errorbars[,1],3)+c(rep(0,9),rep(-0.01,18)),c(means-sds,means-sds,means+sds),rep(RC$errorbars[,1],3)+c(rep(0,9),rep(0.01,18)),c(means+sds,means-sds,means+sds))
text(1,0.25,paste('linear 6 mer (',RC$auc_score,')'),pos=2);par(new=T)
RC<-ROC_5mer;sds<-apply(RC$errorbars[,2:6],1,sd);means<-apply(RC$errorbars[,2:6],1,mean)
plot(RC$DF[,1],RC$DF[,7],type='l',lty="dashed",col=pallets[2],xlim=c(0,1),ylim=c(0,1),labels=F);segments(rep(RC$errorbars[,1],3)+c(rep(0,9),rep(-0.01,18)),c(means-sds,means-sds,means+sds),rep(RC$errorbars[,1],3)+c(rep(0,9),rep(0.01,18)),c(means+sds,means-sds,means+sds),col=pallets[2])
text(1,0.3,paste('linear 5 mer (',RC$auc_score,')'),pos=2);par(new=T)
RC<-ROC_7mer;sds<-apply(RC$errorbars[,2:6],1,sd);means<-apply(RC$errorbars[,2:6],1,mean)
plot(RC$DF[,1],RC$DF[,7],type='l',lty="dotted",col=pallets[4],xlim=c(0,1),ylim=c(0,1));segments(rep(RC$errorbars[,1],3)+c(rep(0,9),rep(-0.01,18)),c(means-sds,means-sds,means+sds),rep(RC$errorbars[,1],3)+c(rep(0,9),rep(0.01,18)),c(means+sds,means-sds,means+sds),col=pallets[4])
text(1,0.2,paste('linear 7 mer (',RC$auc_score,')'),pos=2);par(new=T)
RC<-ROC_kern;sds<-apply(RC$errorbars[,2:6],1,sd);means<-apply(RC$errorbars[,2:6],1,mean)
plot(RC$DF[,1],RC$DF[,7],type='l',lty="twodash",col=pallets[10],xlim=c(0,1),ylim=c(0,1));segments(rep(RC$errorbars[,1],3)+c(rep(0,9),rep(-0.01,18)),c(means-sds,means-sds,means+sds),rep(RC$errorbars[,1],3)+c(rep(0,9),rep(0.01,18)),c(means+sds,means-sds,means+sds),col=pallets[10])
text(1,0.15,paste('kernel 6 mer (',RC$auc_score,')'),pos=2)



# ------------ Fig.4 d, Supplementary Fig.9 cd  ------------
weights1<-read.table("data/Figure4/weights/ChIP8_ATX1_weights",header=F);weights1<-add_mean(weights1)
weights2<-read.table("data/Figure4/weights/ChIP8_ATX2_weights",header=F);weights2<-add_mean(weights2)
weights2r<-read.table("data/Figure4/weights/ChIP7_ATX2_weights",header=F);weights2r<-add_mean(weights2r)
weights1r<-read.table("data/Figure4/weights/ChIP7_ATX1_weights",header=F);weights1r<-add_mean(weights1r)
par(mar=c(2.5,2.5,0.2,0.2))
jpeg('SFig4d.jpeg',height=2160/3,width=2160/3,res=300);par(mar=c(2.5,2.5,0.2,0.2),cex=0.4,cex.axis=2)
plot(weights1$V7,weights2$V7,pch=16);cor.test(weights1$V7,weights2$V7) # Fig 4d
plot(weights1$V7,weights1r$V7,pch=16,cex.axis=2,xlab='',ylab='',lwd=2);cor.test(weights1$V7,weights1r$V7) # Supplementary Fig 9c
plot(weights2$V7,weights2r$V7,pch=16,cex.axis=2);cor.test(weights2$V7,weights2r$V7) # Supplementary Fig 9d
dev.off()
# make Supplementary Table.1 
ST1<-data.frame(weights1[,c(1,7)],weights2$V7,weights1r$V7,weights2r$V7);colnames(ST1)<-c('6mer','ATX1_model','ATX2_model','ATX1_model_replicate','ATX2_model_replicate')
#write.csv(ST1,file='Supplementary_Table1.csv')

# -------------- Fig.4 e-f, SFig10a-c, SFig11 string/cluster visualization of weighted 6-mers  -------------------------
# 6-mer list was organized like below; 
W<-weights1
Sorted_W<-W[order(-W[,7]),];ZSorted_W<-W[order(abs(W[,7])),]
pos<-data.frame(Sorted_W[1:60,c(1,7)],rownames(Sorted_W[1:60,c(1,7)]));colnames(pos)<-c('kmer','weights','index')
neg<-data.frame(Sorted_W[nrow(W):(nrow(W)-59),c(1,7)],rownames(Sorted_W[nrow(W):(nrow(W)-59),c(1,7)]));colnames(neg)<-c('kmer','weights','index')
rid<-sample(c(1:nrow(W)),60);rand<-data.frame(Sorted_W[rid,c(1,7)],rownames(Sorted_W[rid,c(1,7)]));colnames(rand)<-c('kmer','weights','index')
zero<-data.frame(ZSorted_W[1:60,c(1,7)],rownames(ZSorted_W[1:60,c(1,7)]));colnames(zero)<-c('kmer','weights','index')
WRITE<-data.frame(pos$kmer,neg$kmer,rand$kmer,zero$kmer)
# write.table(file="data/Figure4/SVM_string/ChIP8_ATX1_KMOIS.txt",t(WRITE),quote=F,row.names = F) 
# 6-mers in this K-mer Of Interest table (KMOI) were then passed to SVM_motif.ipyenb, where nucleotide frequency of 3 base flanking each 6-mers were calcutlated and converted to meme format. 
# the resulting meme format file representing each 6mers were searched for matching motifs with TOMTOM. see motif_search.md
# tomtom $query all_a_thal_uniq.meme -o "all_thale_withflank/"$query  # all_a_thal_uniq.meme were made by concatenating all the Arabidopsis databases included in the MEME Suite's motif_databases.12.19 and removind duplicates.

#-----INPUT FILES: change readfile names for ATX1 or ATX2. set k.
k<-1 #1,2,3,4=positive 6-mers,negative 6-mers,random 6-mers, nearly zero-weighted 6-mers
#read TOMTOM annotation hits. the first rows of these files corresponds to the Supplementary Table.2 
fATX_n<-read.csv("data/Figure4/SVM_string/ATX2_neg0.1_.csv") 
fATX_p<-read.csv("data/Figure4/SVM_string/ATX2_pos0.1_.csv")
fATX_r<-read.csv("data/Figure4/SVM_string/ATX2_rand0.1_.csv")
fATX_z<-read.csv("data/Figure4/SVM_string/ATX2_zeros0.1_.csv")
DB<-list(fATX_p,fATX_n,fATX_r,fATX_z)
#read weighted 6-mers 
kmsl<-read.table("data/Figure4/SVM_string/ATX2_KMOIS.txt",header=TRUE)


#The symbolic notation of CISBP is replaced by the following string
CISBPs<-c('M0256_1.02','M0581_1.02','M0588_1.02','M0591_1.02','M0768_1.02','M1276_1.02','M1694_1.02','M2371_1.02','M4266_1.02','M0844_1.02','M0155_1.02','M1345_1.02','M1337_1.02','M1562_1.02','M1657_1.02','M1303_1.02','M1309_1.02','M2370_1.02','M1684_1.02','M1329_1.02','M1333_1.02')
CISBPn<-c('ATAREB1_CISBP','CAMTA3_CISBP','CxC_AT2G20110_CISBP','CxC_AT4G29000_CISBP','GATA9_CISBP','GT-1_CISBP','WRKY18_CISBP','AtSPL8_CISBP','scTBP_CISBP','ATHB22_CISBP','bHLH_HBI1_CISBP','ARR18_CISBP','homeodomain_like_AT5G05090_CISBP','SPL11_CISBP','TCP21_CISBP','RVE6_CISBP','MYBD_CISBP','AtSPL3_CISBP','WRKY60_CISBP','MYB94_CISBP','MYB98MYB98')
patterns=c('WRKY','MYB','TCP','bHLH','C2H2','SPL','C3H','NAC','C2C2','Trihelix','ABI3VP','TBP','miR')

#fill colors
colorkeys=c(pallets[1:11],'steelblue3','gray80')
naColor='gray90';others='gray50' 
dict=data.frame(patterns,colorkeys);sorter=data.frame(c(patterns,'NA'),1:(length(patterns)+1));colnames(sorter)<-c("Var1","sortbythis")
all_colors_list=c(colorkeys,naColor)

#qgraphs-------------------
feature_list<-as.character(unlist(kmsl[k,]));feature_list_s<-sort(feature_list)
DNAS<-list();for(i in 1:length(feature_list)){DNAS[i]<-DNAString(feature_list_s[i])}
RC_DNAS<-list();for(i in 1:length(feature_list)){RC_DNAS[i]<-reverseComplement(DNAS[[i]])}
M<-kmer_list_to_matrix_dir(feature_list);dist_mi<-M[[1]];X<-M[[2]] 

#make filling color string by TOMTOM hits
CSV<-t(DB[[k]]);a<-CSV[,1] 
for(i in 1:length(CISBPn)){a[grep(CISBPs[i],a)]<-CISBPn[i]}
for(p in 1:length(patterns)){a[grep(patterns[p],a)]<-patterns[p]}
a<-replace_string(a,dict);a[is.na(a)]<-naColor;
a=replace(a,!(a %in% all_colors_list),others) 

# make border color string by grep 
border_colors=rep(NA,60) 
pal<-c('gray30',pallets[7],pallets[4],'purple')
border_colors<-replace(border_colors,grep('AAGAGA|GAGAGA|AAGAGG|GAGAGG|AGAGAA|GGAGAA|AGAGAG|GGAGAG',feature_list),pal[1]);border_colors  # RAGAGR or RGAGAR
border_colors<-replace(border_colors,grep('AAACCC|AACCCT|ACCCAA|ACCCTA',feature_list),pal[2]) ;border_colors  #telobox
border_colors<-replace(border_colors,grep('AAGCCC|AGGCCC|GGCCCA|AGCCCA|GCCCAA|GCCCAT|GGGCTT|GGGCCT|TGGGCC|TGGGCT|TTGGGC|ATGGGC|AATGGG|CCCATT',feature_list),pal[3]);border_colors #my own consensus for ARGCCCAWT 
border_colors<-replace(border_colors,grep('TATAAA|TATATA|ATATAA|ATAAAT|TAAATA|ATATAT|TTATAA|TTATAT',feature_list),pal[4]);border_colors  # TATA box from yamamoto 

shapes<-rep('circle',60)

#draw qgraph!
dev.off()
par(mar=c(0,0,0,0))
Q<-qgraph(dist_mi, layout='spring',shape=shapes[order(feature_list)],vsize=3,color=a,border.color=border_colors[order(feature_list)],border.width=2,labels=feature_list_s,node.label.position=4,label.cex=2.5,edge.color='forestgreen')
if(k<3){text(Q$layout[,1],Q$layout[,2],order(feature_list),cex=0.8)}




# ---------- Fig. 4h -----------------------------------
# overlap with sppRNAs

TAGg<-read.table('data/RF_for_arabisopsis/TAGg.txt')
load(file='data/Figure4/sppRNA_uniq_gene_list_for_R')
S<-TAGg[order(-TAGg$ATX1_ChIP8TSS),];ATX1_bound<-S[1:3000,1]
S<-TAGg[order(-TAGg$ATX2_ChIP8TSS),];ATX2_bound<-S[1:3000,1]

q<-length(intersect(ATX1_bound,sppRNA_LIST$hen2))
m<-3000
n<-27409-m
k<-length(sppRNA_LIST$hen2)
phyper(q, m, n, k,lower.tail = FALSE,log.p = FALSE) # probability
gLIST<-list(bound=ATX2_bound,sppRNA=sppRNA_LIST$hen2) 
plot(euler(gLIST),shape = "euler",quantities = TRUE,fill=pallets[c(4,8)],col=pallets[c(4,8)]) 
