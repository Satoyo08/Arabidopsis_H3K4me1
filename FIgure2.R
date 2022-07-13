# What this script do: make plots in Figure2, as well as Supplementary Fig.5c, Fig.6a etc

source('custom_functions.r',echo=FALSE)
setwd('data/Figure2')

# --------- Fig.2b,c -------
# The heatmap and metaplot were generated using ngs.plot.r with default option

# ---------- Fig.2d, Supplementary Fig.5c, Fig.6a --------
# first, peaks were called using MACS2;
# macs2 callpeak -t $sample_bam -c $control_bam -n $output_file -f BAM -g 1.3e8 -q 0.3 --nomodel
# then the distance between peak summits and nearest TSS/TES was calculated using bedtools closest;
# bedtools closest -a $output_file'_summits.bed' -b sorted_araport11_all_exact_TSS.bed -D b -t all | cut -f 1,2,3,9,11,12 |grep -v -e 'Pt' -e 'Mt' -e 'M' -e 'C'> $output_file'_summits.bed_closest_TSS_prehist'
# cut -f 6 $output_file'_summits.bed_closest_TSS_prehist' |sort -n |uniq -c > $output_file'_summits.bed_closest_TSS_hist'

# load peak position data generated as described above
in_list<-read.table('input_list',head=F)
keywords<-c('ChIP7_ATX1_nomodel_summits.bed_closest_TSS_hist','ChIP7_ATX2_nomodel_summits.bed_closest_TSS_hist','soinagak_ATXR7_nomodel_summits.bed_closest_TES_hist') # For Supplementary Fig.5c
keywords<-c('covaris_re_ATX1_nomodel_summits.bed_closest_TSS_hist','covaris_re_ATX2_nomodel_summits.bed_closest_TSS_hist','covaris_re_ATXR7_nomodel_summits.bed_closest_TES_hist') # For Fig.2d
keywords<-c('ATXR3_nomodel_summits.bed_closest_TSS_hist') # For Fig.6

par(mar=c(3,3,2,0.2),mfrow=c(length(keywords),1));i<-1
for(key in keywords){
  openthis<-in_list[grep(key,in_list$V1),];print(openthis)
  ip<-read.table(as.character(openthis),head=F)
  plot(ip[,2],ip[,1]-1,col=rgb(0,0,0,0),xlim=c(-600,600),type='l',ylim=c(0,max(ip[,1])),axes=F)
  pos<-par('usr')
  par(new=T);plot(ip[,2],ip[,1]-1,col='gray40',xlim=c(-600,600),type='l',ylim=c(0,max(ip[,1])))
  if(i<=2){
    axis(side=1,at=c(-150,200,600))  
  }
  abline(v=0,lty='dashed')
  i<-i+1
}
