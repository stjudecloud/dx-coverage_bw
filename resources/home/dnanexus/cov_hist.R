#!/usr/bin/Rscript
#Name: cov_hist.R
#Author: Gang Wu <Gang.Wu@stjude.org>
#Requirement: there must be paired 'control' and 'case' [wg|exon].coverage.counts.txt in the current directory
#             for unpaired, use unpaired_cov_hist.R instead

# $1 = control sample, e.g., SJLGG042_G
# $2 = case sample, e.g., SJLGG042_D
# $3 = coverage ceiling (for histogram plot), e.g. 101

#parse the arguments
args=(commandArgs(TRUE))
control=args[1]
case=args[2]
covs=as.numeric(args[3])

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
## genomic coverage####
wg.g<-read.table(paste(control,".wg.coverage.counts.txt",sep=""),header=T,sep="\t")
wg.d<-read.table(paste(case,".wg.coverage.counts.txt",sep=""),header=T,sep="\t")
wg.gd.total <- cbind(wg.g$Total[1:covs],wg.d$Total[1:covs])
wg.gd.prop<-100*prop.table(wg.gd.total/100,2)
imin<-min(dim(wg.g)[[1]],dim(wg.d)[[1]])

all.chrs=grep("chr",colnames(wg.g))
png(paste(case,".coverage.distribution.chr-by-chr.png",sep=""),height=2000,width=3000,type="cairo")
par(mfrow=c(6,4))
for(j in 1:length(all.chrs)){
	chr=paste(case," ",colnames(wg.g)[all.chrs[j]],sep="")
	wg.chr<-cbind(wg.g[1:imin,all.chrs[j]],wg.d[1:imin,all.chrs[j]])
	wg.chr.prop<-100*prop.table(wg.chr/100,2)
	barplot(as.matrix(t(wg.chr.prop[1:covs,])),names.arg=wg.g$Coverage[1:covs],ylim=c(0,5),
		main=chr,cex.main=3,xlab="Coverage",ylab="Percent of Total (%)",args.legend=list(cex=2),
		  axis.lty=1,cex.names=0.8,legend.text=c("Normal","Disease"),beside=T,col=c("green","red"))
}
dev.off()

### exonic coverage ####
exon.g<-read.table(paste(control,".exon.coverage.counts.txt",sep=""),header=T,sep="\t")
exon.d<-read.table(paste(case,".exon.coverage.counts.txt",sep=""),header=T,sep="\t")
xmin<-min(dim(exon.g)[[1]],dim(exon.d)[[1]])
exon.gd.total <- cbind(exon.g$Total[1:xmin],exon.d$Total[1:xmin])
exon.gd.prop <- as.data.frame(prop.table(exon.gd.total/100,2))
colnames(exon.gd.prop)<-c("ex.g","ex.d")
rownames(exon.gd.prop)<-exon.g$Coverage[1:xmin]
exon.gd.prop$germline<-cumsum(exon.gd.prop[,1])
exon.gd.prop$diagnosis<-cumsum(exon.gd.prop[,2])

### combine genomic and exonic coverage ####
we.gd.total<-cbind(wg.g$Total[1:covs],wg.d$Total[1:covs],exon.g$Total[1:covs],exon.d$Total[1:covs])
we.gd.prop<-as.data.frame(100*prop.table(we.gd.total/100,2))
colnames(we.gd.prop)<-c("wg.g","wg.d","ex.g","ex.d")
we.gd.prop$wgg.cum<-my.rcumsum(we.gd.prop$wg.g)
we.gd.prop$wgd.cum<-my.rcumsum(we.gd.prop$wg.d)
we.gd.prop$exg.cum<-my.rcumsum(we.gd.prop$ex.g)
we.gd.prop$exd.cum<-my.rcumsum(we.gd.prop$ex.d)
we.gd.prop$cov<- seq(0,(covs-1),1)
png(paste(case,"cumulative_plot","png",sep="."),width=400,height=400,type="cairo")
plot(wgg.cum ~ cov,data=we.gd.prop, xlab="Coverage",ylab="Percent of Bases (%)",
     main=paste(case," cumulative coverage",sep=""),col="green",pch=22)
points(wgd.cum ~ cov,data=we.gd.prop,col="red",pch=22)
points(exg.cum ~ cov,data=we.gd.prop,col="green",pch=19)
points(exd.cum ~ cov,data=we.gd.prop,col="red",pch=19)
legend("topright",c("Genomic Normal","Genomic Disease","Exonic Normal","Exonic Disease"),
      pch=c(22,22,19,19),col=c("green","red"),bty="n",cex=0.9)
dev.off()

png(paste(case,"coverage.distribution.all.chrs","png",sep="."),height=1200,width=1600,type="cairo")
par(mfrow=c(2,1))
barplot(as.matrix(t(we.gd.prop[,c(1,3)])),names.arg=wg.g$Coverage[1:covs],ylim=c(0,5),
	main="Normal Coverage",xlab="Coverage",ylab="Percent of Total (%)",
	  axis.lty=1,cex.names=0.7,legend.text=c("Genomic","Exonic"),beside=T,col=c("green",rgb(0,30/255,0)))

barplot(as.matrix(t(we.gd.prop[,c(2,4)])),names.arg=exon.g$Coverage[1:covs],ylim=c(0,5),
	main="Disease Coverage",xlab="Coverage",ylab="Percent of Total (%)",
	  axis.lty=1,cex.names=0.7,legend.text=c("Genomic","Exonic"),beside=T,col=c("red",rgb(30/255,0,0)))
dev.off()


msns<-rbind(
	mymsn(exon.g[2:covs,c("Coverage","Total")],0.05),
	mymsn(exon.d[2:covs,c("Coverage","Total")],0.05),
	mymsn(wg.g[2:covs,c("Coverage","Total")],0.05),
	mymsn(wg.d[2:covs,c("Coverage","Total")],0.05)
)
dimnames(msns)[[1]]<-c("exon.germline","exon.diagnosis","wg.germline","wg.diagnosis")
write.table(msns,paste(case,"coverage.summary.table","txt",sep="."),quote=F,row.names=T,col.names=T,sep="\t")

