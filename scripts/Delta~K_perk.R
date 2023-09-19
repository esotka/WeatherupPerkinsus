
### make k plot
rm(list=ls())
ks <- paste(c("K1","K2","K3","K4","K5","K6","K7","K8","K9","K10"),"run",sep="")
filename <- list.files(path="data/ngsadmix_perk/",pattern=".log")
lnl <- c()
for (k in 1:length(ks))
{
  filename.k <- filename[grep(pattern=ks[k],filename)]
  for (i in 1:length(filename.k))
  {
    liketmp <- readLines(paste("data/ngsadmix_perk/",filename.k[i],sep=""))
    liketmp <- liketmp[grep(pattern = "best like=",liketmp)]
    lnl <- rbind(lnl,
                 data.frame(k=ks[k],run=i,
                 lnl=substr(liketmp,regexpr(pattern="=",liketmp)[1]+1,regexpr(pattern="=",liketmp)[1]+10)))
    
  }
}
lnl$lnl <- as.numeric(as.character(lnl$lnl))
lnl$k = factor(lnl$k); lnl$k = factor(lnl$k,levels(lnl$k)[c(2:10,1)])
plot(lnl~k,data=lnl)
options(scipen=999)
xbar <- tapply(lnl$lnl,lnl$k,mean,na.rm=T)
std <- tapply(lnl$lnl,lnl$k,sd,na.rm=T)+0.000001
out <- data.frame(xbar,std)
# Evanno et al. 2005 ∆K = m(|L(K + 1) − 2 L(K ) + L(K − 1)|)/s[L(K )]
out$L.prime.k <- c(NA,xbar[-1]-xbar[-length(xbar)])
out$L.dblprime.k <- c(NA,out$L.prime.k[-c(1,length(xbar))]-(out$L.prime.k)[-(1:2)],NA)
out$delta <- out$L.dblprime.k/out$std
print(out)
pdf("output/EvannoOutput~K_perk.pdf",width=5,height=3)
par(mar=c(4,4,1,1))
plot(out$L.dblprime.k,xaxt="n",xlab="k",type="b",ylab="L.dblprime.k")
mtext(at=1:length(xbar),rownames(out),side = 1,line=1,cex=.7)
dev.off()
