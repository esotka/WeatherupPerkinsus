### admix
rm(list=ls())
ks <- c("K2","K3","K4","K5","K6","K7","K8","K9","K10")
meta = read.csv("../meta_36ind_final.csv")
state_reg = factor(paste(meta$state,meta$VA_reg))

state_reg_ordered = state_reg[order(state_reg)]

pdf("admix.pdf",width=5,height=3)
par(mar=c(4,4,1,1))
for (k in 1:length(ks))
{
  dat <- read.delim(paste(ks[k],"run1.qopt",sep=""),sep=" ",header = F)
  dat <- dat[order(state_reg),]
  dat <- dat[,-(dim(dat)[2])]
  barplot(t(dat),col=1:nrow(dat),space=0,border=NA,xlab="",ylab="admixture",
          names.arg = rep("",nrow(dat)),main=ks[k])
  # labels
    x <- 1:dim(dat)[1]
  for(i in 1:length(levels(state_reg)))
    {
  xtmp <- x[state_reg_ordered==levels(state_reg)[i]]
  segments(max(xtmp),1,max(xtmp),-.05,col="white",lwd=1)
  mtext(side=1,at=mean(xtmp),levels(state_reg)[i],cex=.5,line=.5)
  }
}
dev.off()



