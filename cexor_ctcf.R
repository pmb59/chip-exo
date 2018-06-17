## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
list.of.packages <- c("CexoR")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) biocLite("CexoR")

library("CexoR")
library('GenomicRanges')
W = 50 
counter = 0
CHR = c()
CEXOR = c()
RHEE_PUGH = c()
PROP = c()
PROP2 = c()
OVL= c()
OVL2= c()

for (i in c('Y',"X",1:22)  ){
  
  counter <- counter + 1

  rep1 <- paste("rep1_chr",i,".bam", sep='')
  rep2 <- paste("rep2_chr",i,".bam", sep='')
  rep3 <- paste("rep3_chr",i,".bam", sep='')

  chrom <- read.table('hg19.chrom.sizes.txt', head=FALSE)
  chrSel  <- which( as.character(chrom$V1)  == paste("chr",i,sep='') )

  chipexo <- cexor(bam=c(rep1,rep2,rep3), chrN=as.character(chrom$V1[chrSel]), chrL=chrom$V2[chrSel] , idr=0.01, N=1e6, 
                    p=1e-9, dpeaks=c(0,150), dpairs=200)  

  system(paste("mv cexor_peak_centres_IDR_0.01.bed ", paste("chr",i,sep=''),"_cexor_peak_centres.bed",sep="" )  )
  system(paste("mv cexor_peak_events_IDR_0.01.bed ",  paste("chr",i,sep=''),"_cexor_peak_events.bed",sep="" )  )
  
  pdf(file=paste("chr",i,'.pdf',sep='') , height=7.5, width=4)
  plotcexor(bam=c(rep1,rep2,rep3), peaks=chipexo, EXT=500)
  dev.off()
  
  #Get overlaps
  reference <- read.table('CTCF_hg19.bed', head=F)
  reference$V3 <- reference$V2 + W
  reference$V2 <- reference$V2 - W 

  gr <- GRanges(
    seqnames=Rle(as.character(reference$V1)) ,
    ranges=IRanges(reference$V2, reference$V3 )
  )

  ovl <- countOverlaps(query=chipexo$bindingEvents, subject=gr,
                maxgap=-1L, minoverlap=0L,
                type=c("any"),  #"any", "start", "end", "within", "equal"
                ignore.strand=FALSE)

  OVL[counter] <- length(which(ovl!=0))
  
  CHR[counter] <- paste("chr",i,sep='')
  CEXOR[counter] <- length( chipexo$bindingEvents )
  PROP[counter] <- 100* (length(which(ovl!=0))  / CEXOR[counter] )
  RHEE_PUGH[counter] <- length( which(reference$V1 == paste("chr",i,sep='')) )
  rm(ovl)
  
  ovl <- countOverlaps(query=gr, subject=chipexo$bindingEvents,
                       maxgap=-1L, minoverlap=0L,
                       type=c("any"),  #"any", "start", "end", "within", "equal"
                       ignore.strand=FALSE)
  PROP2[counter] <- 100* (length(which(ovl!=0))  / RHEE_PUGH[counter] )
  OVL2[counter] <- length(which(ovl!=0))
  #write.table(x=data.frame(chr=CHR, rhee_pugh=RHEE_PUGH,CexoR=CEXOR, ovl_rhee_pugh=OVL2, prop_rhee_pugh_cov=PROP2, ovl_CexoR=OVL, prop_CexoR_cov=PROP), file = "proportions.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE,col.names = TRUE)

}
#save R object
x=data.frame(chr=CHR, rhee_pugh=RHEE_PUGH,CexoR=CEXOR, ovl_rhee_pugh=OVL2, prop_rhee_pugh_cov=PROP2, ovl_CexoR=OVL, prop_CexoR_cov=PROP)
save(x,file='cexor_ctcf.Rdata')



