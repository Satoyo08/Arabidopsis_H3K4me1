
WB<-read.table('data/refs/WB_signal_counts.csv',sep=',',skip=1)
means<-rev(apply(WB[2:13,],1,mean,na.rm=T))#;means[1]<-NAs
sds<-rev(apply(WB[,2:13],1,sd,na.rm=T))
plot(1,1,ylim=c(0,1.2),xlim=c(0.5,3.5)
     ,ylab='mutant/WT(H3K4me1/H4,H3)',xaxt='n',xlab='',type='n')
rect(c(1,2,3)-0.1,means,c(1,2,3)+0.3,0,col='gray80')
points(rep(nrow(WB):1,ncol(WB))+runif(nrow(WB)*ncol(WB))/5,unlist(WB[1:nrow(WB),1:ncol(WB)]),pch=16,cex=0.7)
segments(c(1,2,3)+0.1,means-sds,c(1,2,3)+0.1,means+sds)
l<-0.03
segments(c(1,2,3)+0.1-l,means-sds,c(1,2,3)+0.1+l,means-sds)
segments(c(1,2,3)+0.1-l,means+sds,c(1,2,3)+0.1+l,means+sds)
