#### genetic difference within pop'ns (use all popns)
library(phytools)
library(ape)
rm(list=ls())

pdf("output/CoPhylogeny_vir_vs_perk.pdf",height=7,width=10)
obj<-cophylo(read.tree("data/phylogenies/Perk_36ind-PhyML_tree.tr"),read.tree("data/phylogenies/Cvir_36ind_52100snp.fas.contree"),rotate=T)
state = factor(substr(obj$assoc,1,2))
plot(obj,link.col=c("black","black","red")[state])
text(x=c(.25,-.25),y=c(1.03,1.03),c("Oyster","Perkinsus"))
dev.off()
