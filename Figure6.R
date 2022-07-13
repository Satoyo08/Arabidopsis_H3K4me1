source('custom_functions.r');source('RF_functions.R')
library(ROCR);library(tidyverse)


# Figure 6a
ip<-read.table('data/Figure6/ATXR3_nomodel_summits.bed_closest_TSS_hist',head=F)
plot(ip[,2],ip[,1]-1,col='gray40',xlim=c(-600,1000),type='l',ylim=c(0,max(ip[,1])))

# Figure 6b-c ; see Figure3.r
# Figure 6d
S<-TAG_R3[order(-TAG_R3$SDG2_TSS2),];posID<-S[1:3000,1];
posID<-S[1:3000,1]
feature<-c("mid_gb_k_CMA603","TSS2_CMA601")
AL<-AL_R
AL[,ncol(AL)+1]<-'A'
AL[AL$ID %in% posID,ncol(AL)]<-'B'
col=pallets[5];i<-1;vio<-vio_pairs_log();vio;ggsave(paste('ggplot_violin_log','ATXR3',feature[i],'.pdf',sep=''),vio,width=1.2,height=2.5)
col=pallets[3];i<-2;vio<-vio_pairs_log();vio;ggsave(paste('ggplot_violin_log','ATXR3',feature[i],'.pdf',sep=''),vio,width=1.2,height=2.5)
# Figure 6 e-h ; see Figure5.r
