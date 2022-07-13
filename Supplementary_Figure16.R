# Figure.16ab -----------------------------------------
# relative distance of mouse H3K4methyltransferase peaks - TSS, enhancer
hist_list<-list.files('data/Figure7') 
ROI<-c('TSS','enh_R1','enh_V'); k=1 # 1 for plotting TSS, 2 for R1 enhancer, 3 for V6.5 enhancer
keywords<-c('Set1A_rep1','Set1A_rep2','Mll2','Mll34_rep1','Mll34_rep2')
ylimm<-c(40,55,15,15,15)
in_list<-hist_list[grep(ROI[k],hist_list)]
pdf(file = 'peak_summits_R1.pdf',width=1.5*2,height=2.7*2)
par(mar=c(2,2,0.2,0.5),mfrow=c(length(keywords),1),cex=0.7);i<-1 
for(key in keywords){
  openthis<-in_list[grep(key,in_list)];print(openthis)
  ip<-read.table(as.character(paste('data/Figure7/',openthis,sep='')),head=F)
  if(k==1){ #TSS
    plot(ip[,2],ip[,1]-1,col=rgb(0,0,0,0),xlim=c(-600,600),type='l',ylim=c(0,ylimm[i]),axes=F)
    pos<-par('usr')
    rect(-150,pos[3],300,pos[4],col=pallets[6],border=F)
    par(new=T);plot(ip[,2],ip[,1]-1,col='gray40',xlim=c(-600,600),type='l',ylim=c(0,ylimm[i]),xaxt='n')
    axis(side=1,at=c(-600,-150,0,300,600),line=0)
    
  }else{ # enhencer
    plot(ip[,2],ip[,1]-1,col=rgb(0,0,0,0),xlim=c(-3000,3000),type='l',ylim=c(0,ylimm[i]),axes=F)
    pos<-par('usr')
    rect(-900,pos[3],900,pos[4],col=pallets[6],border=F)
    par(new=T);plot(ip[,2],ip[,1]-1,col='gray40',xlim=c(-3000,3000),type='l',ylim=c(0,ylimm[i]),xaxt='n')
    axis(side=1,at=c(-3000,-2000,-900,0,900,2000,3000),line=0)
  }
  abline(v=0,lty='dashed')
  i<-i+1
}
dev.off()

# for SFig16 c-h, see Figure7.r