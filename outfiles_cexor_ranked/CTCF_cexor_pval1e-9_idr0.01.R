
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
list.of.packages <- c("CexoR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite("CexoR")

library("CexoR")
library('GenomicRanges')


for (i in c('Y',"X",1:22)  ){    
  
  rep1 <- paste("rep1_chr",i,".bam", sep='')
  rep2 <- paste("rep2_chr",i,".bam", sep='')
  rep3 <- paste("rep3_chr",i,".bam", sep='')

  chrom <- read.table('hg19.chrom.sizes.txt', head=FALSE)
  chrSel  <- which( as.character(chrom$V1)  == paste("chr",i,sep='') )

  chipexo <- cexor(bam=c(rep1,rep2,rep3), chrN=as.character(chrom$V1[chrSel]), chrL=chrom$V2[chrSel] , idr=0.01, N=1e6, 
                    p=1e-9, dpeaks=c(0,150), dpairs=200, bedfile=FALSE)  


  #binding events
  be <-  as.data.frame(chipexo$bindingEvents)[,c(1:4,6:9)]
  be$Mean.neg.log10pvalue <- rowMeans( be[,6:8] )
  #head(be)
  
  write.table(x=be, file = paste(paste("chr",i,sep=''),"_cexor_binding_events.bed",sep="" ) , append = FALSE, quote = FALSE, sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  rm(be)
  
  #peak centres
  pc <-  as.data.frame(chipexo$bindingCentres )  [,c(1:3,6:9)]
  pc$Mean.neg.log10pvalue <- rowMeans( pc[,5:7] )
  #head(pc)
  
  write.table(x=pc, file = paste(paste("chr",i,sep=''),"_cexor_peak_centres.bed",sep="" ) , append = FALSE, quote = FALSE, sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  rm(pc)
  
  rm(rep1,rep2,rep3,chrom,chrSel,chipexo)
}


system('head -n 1 chr1_cexor_binding_events.bed > head_cexor_binding_events.bed' )
system('tail -q -n +2 *_cexor_binding_events.bed >  temp.bed')
system('cat head_cexor_binding_events.bed temp.bed >  CTCF_cexor_binding_events_pval1e-9_idr0.01.bed')
system('rm temp.bed')

system('head -n 1 chr1_cexor_peak_centres.bed > head_cexor_peak_centres.bed' )
system('tail -q -n +2 *_cexor_peak_centres.bed >  temp.bed')
system('cat head_cexor_peak_centres.bed temp.bed >  CTCF_cexor_peak_centres_pval1e-9_idr0.01.bed')
system('rm temp.bed')


#sort
be <- read.table('CTCF_cexor_binding_events_pval1e-9_idr0.01.bed', head=TRUE)
write.table(x=be[order(-be$Mean.neg.log10pvalue),], file = 'CTCF_cexor_binding_events_pval1e-9_idr0.01_ranked.bed' , 
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

pc <- read.table('CTCF_cexor_peak_centres_pval1e-9_idr0.01.bed', head=TRUE)
write.table(x=pc[order(-pc$Mean.neg.log10pvalue),], file = 'CTCF_cexor_peak_centres_pval1e-9_idr0.01_ranked.bed' , 
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = TRUE)
