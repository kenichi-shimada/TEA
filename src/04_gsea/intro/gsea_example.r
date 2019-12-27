library(RColorBrewer)
library(gplots)

setwd(src_dir);setwd("04_gsea")
source("GSEA.EnrichmentScore.r")

setwd(rda_dir)
eff <- readRDS(file="MS-eff.rds")
cmpds <- sub("X","",rownames(eff))
x <- load(file="slps_tgts.rda") 
x <- load(file="sig_tgts.rda") # dmso.pos,dmso.neg,atoc.pos,atoc.neg
x <- load(file="gsea_data.rda") # cmpds,eff,dicts,list.cmpds
dicts <- dicts[-1]

slp.thres <- 3

## DMSO.pos
sig.1 <- apply(slp.dmso[dmso.pos,],1,function(x){
	idx <- which(!is.na(x))
	i.sig <- which.max(x[idx])
	return(idx[i.sig])
})
names(sig.1) <- dmso.pos

d.eff <- eff$dmso
names(d.eff) <- cmpds

#this.dmso.pos <- dmso.pos
this.dmso.pos <- c("DHB1")#,"ADA2B","ADA2C")
this.i <- 1
d.all <- unique(unlist(dicts[[this.i]]))
sorted.d.all <- names(sort(d.eff[d.all],decreasing=T))
this.d <- unique(unlist(dicts[[this.i]][this.dmso.pos]))
sorted.this.d <- sorted.d.all[sorted.d.all %in% this.d]
gs <- GSEA.EnrichmentScore(sorted.d.all,this.d,weighted.score.type=0)

##
sd <- sort(d.eff[d.all],decreasing=T)
sd <- (abs(sd))^(1/2)*sign(sd)

max.sd <- max(abs(sd))
th <- 100
nd <- as.character(round((sd/max.sd)*th))
uniq.cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(th*2+1)
names(uniq.cols) <- -th:th

##
idx <- which(gs$indicator==1)
cols <- uniq.cols[nd]

setwd(plot_dir)
pdf("res.pdf",width=5,height=4)
plot(gs$RES,main="DHB1",lwd=2,col="green2",
	ylim=range(gs$RES,-0.04),axes=F,xlab="# compound",
	ylab="Running ES",type="n")
abline(h=0)
maxi <- gs$arg.ES
lines(rep(maxi,2),c(0,gs$RES[maxi]),lty=2)
lines(gs$RES,lwd=2,col="green2")
axis(1)
axis(2)

for(i in seq(gs$RES)){
	lines(rep(i,2),c(par()$usr[3],-0.025),col=cols[i])
}

for(i in idx){
	lines(rep(i,2),c(par()$usr[3],-0.025),col=1)
}
abline(h=-0.025)
box()
dev.off()




