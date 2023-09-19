
rm(list=ls())
library(colorRamps)
library(scales)
library(dplyr)

#### PCA plots (perkinsus) ####

meta = read.csv("data/meta_36ind_final.csv")
C <- as.matrix(read.table("data/pmarinus_36ind_final_393loci.cov"))
e <- eigen(C)
inds = meta$file
st_va = factor(paste(meta$state,meta$VA_reg))
reg = factor(meta$region)
pdf("output/Pmar_virg_PC.pdf",width=4,height=8)
par(mfrow=c(2,1),mar=c(2,2,2,2))
plot(e$vectors[,1],e$vectors[,2],xlab="PC1",ylab="PC2",
   col=alpha(blue2red(4)[st_va],.5),pch=c(18,20)[reg],cex=2)
legend("topleft",col = alpha(blue2red(4),.5),legend=levels(st_va),cex=.5,pt.cex=1.5,pch=c(20,20,18,18))

mtext(substitute(italic("A. Perkinsus")),3,line=-1)
for (i in 1:4)
{
  tmp <- e$vectors[st_va==levels(st_va)[i],1:2]
  xbar.x <- mean(tmp[,1])
  xbar.y <- mean(tmp[,2])
  segments(x0 = xbar.x,y0 = xbar.y,x1 = tmp[,1],y1 = tmp[,2],col=alpha(blue2red(4)[i],0.5))
  }

### pc virginica
out = read.table("data/cvir_36ind.52100snp.gt.noNAs.txt",skip=1)
rownames(out) = out[,1]; out = out[,-1]
colnames(out) = meta$file
#### remove outlier 7388-6
out_sub = out[!colnames(out)=="7388-6"]
meta_sub = meta[!meta$file=="7388-6",]
st_va_sub = factor(paste(meta_sub$state,meta_sub$VA_reg))
reg_sub = factor(meta_sub$region)
pc <- prcomp(t(out_sub))
summary(pc) #to get my % explained
plot(pc$x[,1],pc$x[,2],xlab="PC1",ylab="PC2",
   col=alpha(blue2red(4)[st_va_sub],.5),pch=c(18,20)[reg_sub],cex=2)
mtext(substitute(italic("B. Crassostrea")),3,line=-1)
for (i in 1:4)
{
  tmp <- pc$x[st_va_sub==levels(st_va_sub)[i],1:2]
  xbar.x <- mean(tmp[,1])
  xbar.y <- mean(tmp[,2])
  segments(x0 = xbar.x,y0 = xbar.y,x1 = tmp[,1],y1 = tmp[,2],col=alpha(blue2red(4)[i],0.5))
  }

dev.off()
