rm(list=ls())
ks <- c("K2","K3","K4","K5","K6")#,"K7","K8","K9","K10")

############# ENV DATA ##############
meta = read.csv("data/meta_36ind_final.csv")
state_reg = factor(paste(meta$state,meta$VA_reg))
state_reg_ordered = state_reg[order(state_reg)]

kcol <- c("black", #1
          "red",#2
          "yellow",#3
          "gainsboro",#4
          "brown",#5
          "deepskyblue",#6
          "darkgreen",#7
          "dodgerblue4",#8
          "darkorchid",#9
          "burlywood")##10
colorder <- list(
  c(1,2), #ks=2
  c(1,2,3), #ks=3
  c(2,1,3,4), #ks=4
  c(4,2,1,3,5), #ks=5
  c(3,2,1,4,5,6), #ks=6
  c(2,5,3,7,1,6,4), #ks=7
  c(4,7,6,1,2,3,5,8), #ks=8
  c(4,8,3,9,6,5,2,1,7), #ks=9
  c(9,1,2,8,4,6,5,7,10,3)) #ks=10

pdf('output/ngsAdmix.pretty.pdf',width=5,height=4)
par(mfrow=c(11,1),mar=c(.2,0,0,0),xpd = TRUE)

for (k in 1:length(ks))
{
  dat <- read.delim(paste("data/ngsadmix_perk/",ks[k],"run1.qopt",sep=""),sep=" ",header = F)
  dat <- dat[,-(dim(dat)[2])] # remove extra column
  dat <- dat[,order(colorder[[k]])]
  dat <- dat[order(state_reg),]
  fig <- barplot(t(dat),col=kcol[1:length(colorder[[k]])],space=0,border=NA,xlab="",ylab="",
          names.arg = rep("",nrow(dat)),horiz=F)#,main=substr(ks[k],1,3))
  mtext(substr(ks[k],1,3),side=2,line=-1.5,at=.5,cex=0.7)
  x <- 1:dim(dat)[1]
  for(i in 1:length(levels(state_reg_ordered)))
    {
  xtmp <- x[state_reg_ordered==levels(state_reg_ordered)[i]]
  segments(max(xtmp),1,max(xtmp),-.05,col="white",lwd=2)
  mtext(side=1,at=mean(xtmp)-0.5,levels(state_reg)[i],cex=.5,line=.5)
  }
}
xtick<-c(0,3,12,30,36)
axis(side=1, at=xtick, labels = FALSE)
 
    
dev.off()

