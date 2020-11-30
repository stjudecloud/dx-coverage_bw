#!/usr/bin/Rscript
#Name:     cov_boxplots.R
#Author: Gang Wu <Gang.Wu@stjude.org>
# Make coverage boxplots, based on GC content of the exons.
# Requirements: ${id}.{$seg}.coverage.means.with.GC.txt files under the current directory

# $1 = case sample, e.g., SJLGG042_D

#parse the arguments
args=(commandArgs(TRUE))
id=args[1]
seg=args[2]


dat<-read.table(paste(id,seg,"coverage.means.with.GC.txt",sep="."),header=F,sep="\t")
dat$V5=dat$V4/median(dat$V4)
gc25<-quantile(dat[,3],0.25)
gc75<-quantile(dat[,3],0.75)

y.cov=4*median(dat$V4,na.rm=T)

m25<-dat[dat$V3<=gc25 & dat$V2<=30,]
m75<-dat[dat$V3>=gc75 & dat$V2<=30,]
m50<-dat[dat$V3>gc25 & dat$V3<gc75  & dat$V2<=30,]
png(paste(id,seg,"coverage.boxplots.png",sep="."),height=250,width=1000,type="cairo")
par(mfrow=c(1,3))
boxplot(V4~V2,data=m25,ylim=c(0,y.cov),col="green",xaxt="n",
        main=id,xlab=seg,ylab="Average Coverage")
axis(1,1:30,1:30,las=2)
abline(h=mean(dat$V4),lty=2)
legend("topright","Low GC",col="green")
boxplot(V4~V2,data=m50,ylim=c(0,y.cov),col="yellow",xaxt="n",
       main=id,xlab=seg,ylab="Average Coverage")
axis(1,1:30,1:30,las=2)
abline(h=mean(dat$V4),lty=2)
legend("topright","Moderate GC",col="green")
boxplot(V4~V2,data=m75,ylim=c(0,y.cov),col="red",xaxt="n",
        main=id,xlab=seg,ylab="Average Coverage")
axis(1,1:30,1:30,las=2)
abline(h=mean(dat$V4),lty=2)
legend("topright","High GC",col="green")
dev.off()

png(paste(id,seg,"relative.coverage.boxplots.png",sep="."),height=250,width=1000,type="cairo")
par(mfrow=c(1,3))
boxplot(V5~V2,data=m25,ylim=c(0,4),col="green",xaxt="n",
        main=id,xlab=seg,ylab="Relative Average Coverage")
axis(1,1:30,1:30,las=2)
abline(h=mean(dat$V5),lty=2)
legend("topright","Low GC",col="green")
boxplot(V5~V2,data=m50,ylim=c(0,4),col="yellow",xaxt="n",
       main=id,xlab=seg,ylab="Relative Average Coverage")
axis(1,1:30,1:30,las=2)
abline(h=mean(dat$V5),lty=2)
legend("topright","Moderate GC",col="green")
boxplot(V5~V2,data=m75,ylim=c(0,4),col="red",xaxt="n",
        main=id,xlab=seg,ylab="Relative Average Coverage")
axis(1,1:30,1:30,las=2)
abline(h=mean(dat$V5),lty=2)
legend("topright","High GC",col="green")
dev.off()

