#!/usr/bin/env Rscript
#Name: unpaired_cov_hist.R
#Author: Gang Wu <Gang.Wu@stjude.org>
#Requirement: there must be disease [wg|exon].coverage.counts.txt in the current directory
#             for unpaired, use unpaired_cov_hist.R instead
#Note: the ceiling cov (covs) is subject to change depending on sequencing depths
#
# $1 = case sample, e.g., SJLGG042_D
# $2 = coverage ceiling (for histogram plot), e.g. 101

#parse the arguments
args=(commandArgs(TRUE))

id=args[1]
covs=as.numeric(args[2])

##++++++++++ define functions ++++++++++++++++++++++
myprod<-function(x){as.numeric(x[1])*as.numeric(x[2])}
mysq<-function(x){((as.numeric(x[1])-as.numeric(x[3]))^2)*as.numeric(x[2])}
mymsn <- function(x,p) {  
  cof<-qnorm(1-p/2)
  subtotal <- apply(x,1,myprod)
  total <- sum(subtotal) #sum(k(i))
  n <- sum(as.numeric(x[,2]))
  x$m<-total/n
  m<-x$m[1]
  this.sq<-sqrt(sum(apply(x,1,mysq))/n)
  out<-data.frame(n=round(n,digits=0),
		  mean=round(m,digits=1),                  
		  sd=round(this.sq,digits=1),
		  lower=round(m-cof*this.sq,digits=1),
		  upper=round(m+cof*this.sq,digits=1),
		  pct25=round(qnorm(0.25,mean=m,sd=this.sq),digits=1),
		  median=round(qnorm(0.5,mean=m,sd=this.sq),digits=1),
		  pct75=round(qnorm(0.75,mean=m,sd=this.sq),digits=1))
  out
}
my.rcumsum <- function(x){
 y=100-cumsum(x)
 y1=c(100,y)
 y1[1:length(y1)-1]
}
chrs<-c(1:22,"X","Y")

##++++++++ process the sample ++++++++++++++++++++++
wg<-read.table(paste(id,".wg.coverage.counts.txt",sep=""),header=T,sep="\t")

png(paste(id,".coverage.distribution.chr-by-chr.png",sep=""),height=2000,width=3000,type="cairo")
par(mfrow=c(6,4))
for(j in 1:length(chrs)){
	chr=paste("Chr",chrs[j],sep="")
	wg.chr<-wg[,j+1]
	wg.chr.prop<-100*(wg.chr/sum(wg.chr))
	barplot(wg.chr.prop[1:covs],names.arg=wg$Coverage[1:covs],ylim=c(0,5),
		main=chr,cex.main=5,xlab="Coverage",ylab="Percent of Total (%)",
		  axis.lty=1,cex.names=0.7,legend.text=id,beside=T,col="grey")
}
dev.off()

exon<-read.table(paste(id,".exon.coverage.counts.txt",sep=""),header=T,sep="\t")

we.total<-cbind(wg$Total[1:covs],exon$Total[1:covs])
we.prop<-as.data.frame(100*prop.table(we.total/100,2))

colnames(we.prop)<-c("wg","ex")
we.prop$wgg.cum<-my.rcumsum(we.prop$wg)
we.prop$exg.cum<-my.rcumsum(we.prop$ex)
we.prop$cov<- seq(0,(covs-1),1)
png(paste(id,"cumulative_plot","png",sep="."),width=400,height=400,type="cairo")
plot(wgg.cum ~ cov,data=we.prop, xlab="Coverage",ylab="Percent of Bases (%)",
     main=paste(id," cumulative coverage",sep=""),col="grey",pch=22)
points(exg.cum ~ cov,data=we.prop,col="black",pch=19)
legend("topright",c("Genomic","Exonic"),pch=c(22,19),col=c("grey","black"),bty="n",cex=0.9)
dev.off()

png(paste(id,"coverage.distribution.all.chrs","png",sep="."),height=600,width=1600,type="cairo")
barplot(as.matrix(t(we.prop[,1:2])),names.arg=wg$Coverage[1:covs],ylim=c(0,5),
	main=id,xlab="Coverage",ylab="Percent of Total (%)",
	  axis.lty=1,cex.names=0.7,legend.text=c("Genomic","Exonic"),beside=T,col=c("grey","black"))
dev.off()


msns<-rbind(
	mymsn(exon[2:covs,c("Coverage","Total")],0.05),
	mymsn(wg[2:covs,c("Coverage","Total")],0.05)
)
dimnames(msns)[[1]]<-c("exon","wg")
write.table(msns,paste(id,"coverage.summary.table","txt",sep="."),quote=F,row.names=T,col.names=T,sep="\t")

