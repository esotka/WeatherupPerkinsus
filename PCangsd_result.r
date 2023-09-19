
library(colorRamps)
library(scales)
rm(list=ls())
C <- as.matrix(read.table("pmarinus_36ind_final_393loci.cov"))
e <- eigen(C)
meta = read.csv("../meta_36ind_final.csv")
inds = readLines("../bamlist36_final")
inds = gsub(".bwa.bam","",inds)
inds = gsub("_downsampled","",inds)
meta = meta[match(inds,meta$file),]
st_va = factor(paste(meta$state,meta$VA_reg))
pdf("PCangsd_result.pdf")
plot(e$vectors[,1],e$vectors[,2],xlab="PC1",ylab="PC2",
   col=alpha(blue2red(4)[st_va],.5),pch=20,cex=3)
legend("topleft",fill=alpha(blue2red(4),.5),legend=levels(st_va),cex=.7)
for (i in 1:4)
{
  tmp <- e$vectors[st_va==levels(st_va)[i],1:2]
  xbar.x <- mean(tmp[,1])
  xbar.y <- mean(tmp[,2])
  segments(x0 = xbar.x,y0 = xbar.y,x1 = tmp[,1],y1 = tmp[,2],col=alpha(blue2red(4)[i]))
}


dev.off()

