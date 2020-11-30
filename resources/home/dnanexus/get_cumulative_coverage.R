#!/usr/bin/env Rscript
# Get the cumulative coverage stats
#
# input: *.coverage.counts.txt from coverage summary pipeline, for wg, exon, codingexon, etc

# $1 = prefix, e.g., "SJHYPOALL012_G.codingexon", or SJTALL004_D.wg

#parse the arguments
args=(commandArgs(TRUE))
prefix=args[1]

exon.d<-read.table(paste(prefix,"coverage.counts.txt",sep="."),header=T,sep="\t")
rownames(exon.d)<-exon.d$Coverage
exon.d$pct<-100*exon.d$Total/sum(as.numeric(exon.d$Total))
exon.d$cumsum<-cumsum(exon.d$pct)
exon.d$rcumsum=100-exon.d$cumsum
write.table(round(t(exon.d[c(1,2,5,10,20,30,40,45),c("pct","cumsum","rcumsum")]),digits=3),
            paste(prefix,"cumulative.txt",sep="."),quote=F,row.names=T,col.names=T,sep="\t")


