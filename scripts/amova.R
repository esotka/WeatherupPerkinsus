### AMOVA on mpgl

library(pegas)
rm(list=ls())
meta <- read.csv('data/meta_36ind_final.csv')
# check meta$file==rownames(gprob) ALL TRUE
gprob<- t(read.table('data/cvir_36ind.52100snp.gt.noNAs.txt'))
rownames(gprob) = meta$file
ids <- rownames(gprob) 
pop <- factor(paste(meta$state,meta$VA_reg))
reg <- factor(meta$reg)
st.d <- dist(gprob)## distances between rows
print(m1 <- amova(st.d ~ pop,nperm=1000))
sig2 <- setNames(m1$varcomp$sigma2,rownames(m1$varcomp))
print(getPhi(sig2))

### pairwise PhiST (AL, LA, VA_easternShore, VA_Bay)
pairs_to_use = combn(levels(pop),2)
outStats = c()
for (i in 1:ncol(pairs_to_use))
{
gprob2 <- gprob[pop%in%c(pairs_to_use[1,i],pairs_to_use[2,i]),]
ids2 <- rownames(gprob2)
pop2 <- factor(paste(meta$state,meta$VA_reg)[match(ids2,meta$file)])
st.d2 <- dist(gprob2)## distances between rows
m1 <- amova(st.d2 ~ pop2,nperm=1000)
p = round(m1$varcomp[1,2],3)
sig2 <- setNames(m1$varcomp$sigma2,rownames(m1$varcomp))
phist = round(getPhi(sig2)[1],3)
outStats = rbind(outStats,data.frame(pop1=pairs_to_use[1,i],pop2=pairs_to_use[2,i],phist,p))
}

### pairwise PhiST (Gulf, VA_easternShore, VA_Bay)
pop_v2 = as.character(pop)
pop_v2[pop_v2%in%c("LA ","AL ")] = "Gulf"; pop_v2 = factor(pop_v2)
print(m1 <- amova(st.d ~ pop_v2,nperm=1000))
sig2 <- setNames(m1$varcomp$sigma2,rownames(m1$varcomp))
print(getPhi(sig2))

pairs_to_use = combn(levels(pop_v2),2)
for (i in 1:ncol(pairs_to_use))
{
gprob2 <- gprob[pop_v2%in%c(pairs_to_use[1,i],pairs_to_use[2,i]),]
ids2 <- rownames(gprob2)
pop2 <- factor(paste(meta$state,meta$VA_reg)[match(ids2,meta$file)])
st.d2 <- dist(gprob2)## distances between rows
m1 <- amova(st.d2 ~ pop2,nperm=1000)
p = round(m1$varcomp[1,2],3)
sig2 <- setNames(m1$varcomp$sigma2,rownames(m1$varcomp))
phist = round(getPhi(sig2)[1],3)
outStats = rbind(outStats,data.frame(pop1=pairs_to_use[1,i],pop2=pairs_to_use[2,i],phist,p))
}

print(outStats)
