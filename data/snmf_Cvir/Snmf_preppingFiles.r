### make snmf
library(vcfR)
library(poppr)
library(hierfstat)
library(LEA)
rm(list=ls())
dat = read.table("./cvir_36ind.52100snp.gt.noNAs.txt",skip=1)
rownames(dat) = dat[,1]; dat = dat[,-1]
meta = read.csv("../meta_36ind_final.csv")
colnames(dat) = meta$file
rownames(dat) = paste("loc",1:nrow(dat),sep="")
#### remove outlier 7388-6
dat_sub = dat[,!colnames(dat)=="7388-6"]
meta_sub = meta[!meta$file=="7388-6",]
st_va_sub = factor(paste(meta_sub$state,meta_sub$VA_reg))
reg_sub = factor(meta_sub$region)
snp.hf <- df2genind(t(dat_sub),ncode=1)
snp.hf2 <- genind2hierfstat(dat=snp.hf,pop = st_va_sub)
write.struct(snp.hf2,fname="virg.str",MARKERNAMES = F,MISSING=-9)
### remove NAs in first column

#struct2geno("virg.str",ploidy=2,FORMAT = 2) #run this once. it takes awhile.

project = NULL
project = snmf("virg.str.geno", 
               K=2:10,
               entropy = TRUE,
               repetitions = 10, # need to include more
               project = "new") #19 is the max populations sampled
project = load.snmfProject("virg.str.snmfProject")
pdf("snmf_cross-val.pdf")
plot(project, col = "blue", pch = 19, cex = 1.2)
dev.off()

pdf("snmf.pdf",width=10,height=4)
for (k in 2:10)
{
best = which.min(cross.entropy(project, K = k))
qmatrix <- Q(project,K=k,run=best)
pop = st_va_sub
pop_order = factor(pop) # by pop

qmatrix_sort = qmatrix[order(pop_order),]
pop_order2 = pop_order[order(pop_order)]
barplot(t(qmatrix_sort),border=NA,space=0,xlab="Individuals",ylab="Admixture",main=k)
endLine <- as.vector(tapply((1:length(pop_order2)),pop_order2,max))
segments(x0=endLine,y0=-1,x1=endLine,y1=1,col="red",lwd=2)
meanPop <- as.vector(tapply((1:length(pop_order2)),pop_order2,mean))
mtext(levels(pop_order2),at=meanPop,cex=0.5)
}
dev.off()
