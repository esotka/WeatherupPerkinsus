rm(list=ls())
library(LEA)
############# ENV DATA ##############
meta = read.csv("data/meta_36ind_final.csv")
meta_sub = meta[!meta$file=="7388-6",]
st_va_sub = factor(paste(meta_sub$state,meta_sub$VA_reg))
st_va_sub_ordered = st_va_sub[order(st_va_sub)]
reg_sub = factor(meta_sub$region)

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
  c(3,1,2), #ks=3
  c(2,1,4,3), #ks=4
  c(5,4,2,3,1), #ks=5
  c(3,1,6,2,5,4), #ks=6
  c(2,5,3,7,1,6,4), #ks=7
  c(4,7,6,1,2,3,5,8), #ks=8
  c(4,8,3,9,6,5,2,1,7), #ks=9
  c(9,1,2,8,4,6,5,7,10,3)) #ks=10

pdf('output/SNMF.pretty_Cvir.pdf',width=5,height=4)
par(mfrow=c(11,1),mar=c(0.2,0,0,0),xpd = TRUE)

project = load.snmfProject("data/snmf_Cvir/virg.str.snmfProject")
for (k in 1:5)
{
best = which.min(cross.entropy(project, K = k+1))
qmatrix <- Q(project,K=k+1,run=best)
pop = st_va_sub
pop_order = factor(pop) # by pop

qmatrix_sort = qmatrix[order(pop_order),]
pop_order2 = pop_order[order(pop_order)]
fig <- barplot(t(qmatrix_sort),col=kcol[colorder[[k]]],space=0,border=NA,xlab="",ylab="",
          names.arg = rep("",nrow(qmatrix_sort)),horiz=F)#,main=substr(ks[k],1,3))
 mtext(paste("K",k+1,sep=""),side=2,line=-1.5,at=.5,cex=0.7)
  x <- 1:dim(qmatrix)[1]
  for(i in 1:length(levels(st_va_sub_ordered)))
    {
  xtmp <- x[st_va_sub_ordered==levels(st_va_sub_ordered)[i]]
  segments(max(xtmp),1,max(xtmp),-.05,col="white",lwd=2)
  mtext(side=1,at=mean(xtmp)-0.5,levels(st_va_sub)[i],cex=.5,line=.5)
  }
}
xtick<-c(0,3,12,29,35)
axis(side=1, at=xtick, labels = FALSE)
  
dev.off()


