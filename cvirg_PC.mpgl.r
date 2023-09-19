
## MPGL taken from samtools genotype likelihoods (probabilites)
### pulled by "ProcessGL_fromVCF.R"


rm(list=ls())
library(vcfR)
library(readxl)
library(ggplot2)
library(ade4)
library(colorRamps)

geno = read.table("./cvirg_36ind_52458snp.mpgl.txt") # col = loci; row = ind
dim(geno)
meta = read.csv("../meta_36ind_final.csv")
spp <- meta$state
va = meta$VA_reg
spp_va = factor(paste(spp,va,sep=" "))
pc = prcomp(geno)
plot(pc$x[,1],pc$x[,2],xlab="PC1",ylab="PC2",
   col=alpha(blue2red(4)[spp_va],.5),pch=20,cex=3)
legend("bottomleft",fill=alpha(blue2red(4),.5),legend=levels(spp_va),cex=.7)
for (i in 1:4)
{
  tmp <- pc$x[spp_va==levels(spp_va)[i],1:2]
  xbar.x <- mean(tmp[,1])
  xbar.y <- mean(tmp[,2])
  segments(x0 = xbar.x,y0 = xbar.y,x1 = tmp[,1],y1 = tmp[,2],col=alpha(blue2red(4)[i]))
}


