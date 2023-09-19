### convert vcf of genotype calls from bcftools view (CC, CT, etc..) into single letter IUPAC codes and then make a fasta file
rm(list=ls())
library(vcfR) # read.vcfR()
library(seqinr) # bma() Computing an IUPAC nucleotide symbol
library(ape)
dat <- read.vcfR("../virginica_36ind/cvirg_36ind_52100snp.vcf.gz")
gt2 <- extract.gt(dat,return.alleles=T)
loci <- rownames(gt2)
#loci.tmp <- paste(lapply(strsplit(loci,"_"),"[[",1),lapply(strsplit(loci,"_"),"[[",2),sep="_")
#rownames(gt2) <- loci.tmp
colnames(gt2) = gsub(".bwa.bam","",colnames(gt2))


### convert to IAPAC code

library(parallel)

## by locus.

out <- mclapply(1:nrow(gt2), function(i) #
  {
tmp1 <- paste(substr(gt2[i,],1,1),substr(gt2[i,],3,3),sep="")
tmp1 <- ifelse(nchar(tmp1)==1,"N",tmp1)
tmp2 <- mclapply(as.list(tmp1),FUN=function(j) bma(s2c(j)))
})

out2 <- as.data.frame(matrix(unlist(out),nrow=ncol(gt2))) # individuals
out2.noNAs <- out2[,colSums(is.na(out2))==0] ### get rid of all NAs

write.table(x = out2.noNAs,"Cvir_36ind_52100snp.txt",quote=F,row.names = F,col.names = F,sep = "")

dna <- read.table("./Cvir_36ind_52100snp.txt")

# now make a fasta file
prop <- c()
for (i in 1:dim(dna)[1])
{prop <- c(prop,length(gregexpr("n",dna[i,])[[1]])/nchar(dna[i,]))} # prop of missing data
#> table(prop<=0.2)
#TRUE 
# 36 

#ind25 <- colnames(gt2)[prop<=0.25]

#dna25 <- dna[prop<=0.25,]
#write.table(x = dna25,"Perk_allLoci_ind25prop.txt",quote=F,row.names = F,col.names = F,sep = "")
ind = colnames(gt2)
meta = read.csv("../meta_36ind_final.csv")
meta$state_vaReg = paste(meta$state,meta$VA_reg)
ind2 = paste(meta$state_vaReg[match(ind,meta$file)],ind)
ind2 = gsub(" ","",ind2)

for (i in 1:length(ind)) 
{
  tmp <- read.table("Cvir_36ind_52100snp.txt",skip=i-1,nrows=1)
  if (i == 1)
  {
    write.table(x = paste(">",ind2[i],"\n",tmp,"\n",sep=""),"Cvir_36ind_52100snp.fas",quote=F,row.names = F,col.names = F)
  }
  else
  {
    write.table(x = paste(">",ind2[i],"\n",tmp,"\n",sep=""),"Cvir_36ind_52100snp.fas",quote=F,row.names = F,col.names = F,append = T)
  }
}

