### process GL ### 
rm(list=ls())
library(vcfR)
library(parallel)
dat = read.vcfR("./cvirg_36ind_52100snp.vcf.gz")
gt <- extract.gt(dat,element="PL",as.numeric = F)
colnames(gt) = gsub(".bwa.bam","",colnames(gt))

meta = read.csv("../meta_36ind_final.csv")
#colnames(gt)==meta$file ALL ARE TRUE
ninds = ncol(gt)
nloci = nrow(gt)

gp.out <- mclapply(1:nloci, function(i)# by row
{
  tmppl <- gt[i,]
  tmppl.list <- strsplit(tmppl,",")
  tmppl.unlist <- unlist(tmppl.list)
  tmp <- data.frame(aa = as.numeric(tmppl.unlist[seq(1,ninds*3,3)]),#homozygote A
                    ab = as.numeric(tmppl.unlist[seq(2,ninds*3,3)]),#het
                    bb = as.numeric(tmppl.unlist[seq(3,ninds*3,3)]))#homozygote B
  tmp2 <- 10^(tmp/-10) # convert likelihoods
  tsum <- rowSums(tmp2) 
  tstd <- tmp2/tsum #normalize
  tgenest <- c()
  for(j in 0:2){tgenest <- cbind(tgenest,j*tstd[,j+1])} # convert to 0:2 score.
  tgenestSum <- round(rowSums(tgenest),5)
},mc.cores = 4)

gp.out2 <- as.data.frame(matrix(as.numeric(unlist(gp.out)),nrow=ninds)) # rows = ind cols = loci

rownames(gp.out2) = meta$file
colnames(gp.out2) = rownames(gt)
### write mpgl file
options(scipen=999) ### removes scientific notation
write.table(gp.out2,"./cvirg_36ind_52458snp.mpgl.txt",
            sep = " ",quote=F)


