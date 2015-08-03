##################################################
## RoadMap + ENOCDE Enrichment
## Alexander W Drong 15/07/15 
##################################################

## SETUP:
## 1) Download all files from http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/ into ../IN
## 2) Run http://www.broadinstitute.org/mpg/snpsnap/match_snps.html on test SNPs and unzip contents to ../SNPsnap_test_run
## 3) Specify Parameters or leave as default
## 4) Run this Script from scripts directory 

## PARAMETERS

## histone marks to test
histone.marks <- c("H3K4me1.narrowPeak","H3K4me3.narrowPeak","H3K27ac.narrowPeak", "DNase.macs2.narrowPeak")

## cis - interval around SNPs
cis.interval <- 500000
project.name <- "snps_23andMe_clumped_1e5"

##################################################

require(RColorBrewer)
require(foreach)
require(GenomicRanges)
require(doMC)

registerDoMC(cores=8)

matched <- read.table(paste("../SNPsnap_",project.name,"/matched_snps.txt",sep=""), header=T, as.is=T)

anno.granges <- foreach(i=1:10001) %do% {
  print(i)
  chr <- as.numeric(unlist(lapply(strsplit(matched[,i],":"),'[',1)))
  pos <- as.numeric(unlist(lapply(strsplit(matched[,i],":"),'[',2)))
  GRanges(Rle(chr), IRanges(pos - cis.interval, pos + cis.interval), strand="*")
}


anno.granges.list <- GRangesList(anno.granges)


for(hmark in histone.marks){

  files <- dir("../IN")
  files <- files[grep(hmark, files)]

  peak.granges <- foreach(file=files, .combine=rbind) %dopar% {

    print(file)

    dataf <- read.table(paste("../IN/",file,sep=""), as.is=T)
    dataf[,1] <- gsub("chr","",dataf[,1])
    dataf[dataf[,1]==23, 1] <- "X"
    dataf[dataf[,1]==24, 1] <- "Y"

    data.granges <- GRanges(Rle(dataf[, 1]), IRanges(dataf[, 2], dataf[, 3]), strand="*")

    ovs <- countOverlaps(anno.granges.list,data.granges)


    obs <- ovs[1]
    null <- ovs[2:10001]
    p.vals <- sum(null>=obs)/10000
        
    c(p.vals,obs,mean(null))
  }

  peak.granges[peak.granges[,1]==0,1] <- 1e-4 

  rownames(peak.granges) <- substr(files,0,4)

  save(hmark,files,peak.granges, file=paste("../results/",project.name,"results-",hmark,"-",cis.interval,".RData",sep=""))

  samples <- read.table("../roadmap_samples.txt",as.is=T, sep="\t",comment.char="")
  rownames(samples) <- samples[,1]

  samples.id <- intersect(samples[,1], rownames(peak.granges))

  peak.granges <- peak.granges[samples.id,]
  samples <- samples[samples.id,]

  legend <- samples[!duplicated(samples[,2:3]),2:3]

  bonf <- 0.05 / nrow(peak.granges)

  pdf(paste("../plots/",project.name,"_",hmark,"-fold-", cis.interval, ".pdf",sep=""), width=20, height=10)
  layout(rbind(1,2), heights=c(7,1))
  par(mar=c(20, 4, 4, 2))
  y <- (peak.granges[,2]/peak.granges[,3])
  plot(1:nrow(peak.granges),y, col=samples[,3],xaxt="n", xlab="", ylab="Observed/Expected", main=hmark, cex.axis=1.2, cex.lab=1.2,ylim=c(min(y),max(y)+1), pch=c(1,19)[as.numeric(peak.granges[,1]<bonf)+1])
  grid(nx=NA,ny=NULL)
  axis(side=1,at=1:nrow(peak.granges),labels=F, las=2)
  mtext(text=samples[,4], side=1,at=1:nrow(samples),col=samples[,3],las=2,line=1)
  abline(h=bonf)
  vlines <- cumsum(table(factor(samples[,2],levels=unique(samples[,2])))) + .5
  abline(v=vlines,lty=2, col="grey")
  text(x=vlines+0.0,y=max(y)+1,labels=legend[,1], srt=90, col=legend[,2], pos=2)
  par(mar=c(0,0,0,0))
  dev.off()

  pdf(paste("../plots/",project_name,"_",hmark,"-p-", cis.interval, ".pdf",sep=""), width=20, height=10)
  layout(rbind(1,2), heights=c(7,1))
  par(mar=c(20, 4, 4, 2))
  y <- -log10(peak.granges[,1])
  plot(1:nrow(peak.granges),y, col=samples[,3],xaxt="n", xlab="", ylab="-log10(p)", main=hmark, cex.axis=1.2, cex.lab=1.2,ylim=c(min(y),max(y)+2), pch=c(1,19)[as.numeric(peak.granges[,1]<bonf)+1])
  grid(nx=NA,ny=NULL)
  axis(side=1,at=1:nrow(peak.granges),labels=F, las=2)
  mtext(text=samples[,4], side=1,at=1:nrow(samples),col=samples[,3],las=2,line=1)
  abline(h=1)
  vlines <- cumsum(table(factor(samples[,2],levels=unique(samples[,2])))) + .5
  abline(v=vlines,lty=2,col="grey")
  text(x=vlines+0.0,y=max(y)+2,labels=legend[,1], srt=90, col=legend[,2], pos=2)
  par(mar=c(0,0,0,0))
  dev.off()
}
