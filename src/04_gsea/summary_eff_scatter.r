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
	i.sig <- which(x[idx] > slp.thres)
	return(max(i.sig))
})
names(sig.1) <- dmso.pos

d.eff <- eff$dmso
names(d.eff) <- cmpds

dicts.dmso.pos <- lapply(dmso.pos,function(n){
	this.i <- sig.1[n]
	d.all <- unique(unlist(dicts[[this.i]]))
	sorted.d.all <- names(sort(d.eff[d.all],decreasing=T))
	this.d <- dicts[[this.i]][[n]]
	sorted.this.d <- sorted.d.all[sorted.d.all %in% this.d]
	gs <- GSEA.EnrichmentScore(sorted.d.all,this.d,weighted.score.type=0)

	mem <- which(gs$indicator==1)
	le <- sorted.this.d[seq(mem[mem <= gs$arg.ES])]
	return(le)
})
names(dicts.dmso.pos) <- dmso.pos

## selective & essnetial
ois1 <- sapply(dicts.dmso.pos,function(x){
	sapply(dicts.dmso.pos,function(y){
		n1 <- length(x)
		n2 <- length(y)
		n3 <- length(intersect(x,y))
		oi <- n3/min(n1,n2)
		return(oi)
	})
})
rownames(ois1) <- colnames(ois1) <- dmso.pos

## ATOC.pos
sig2 <- apply(slp.atoc[atoc.pos,],1,function(x){
	idx <- which(!is.na(x))
	i.sig <- which(x[idx] > slp.thres)
	return(max(i.sig))
})
names(sig2) <- atoc.pos

a.eff <- eff$atoc
names(a.eff) <- cmpds

dicts.atoc.pos <- lapply(atoc.pos,function(n){
	this.i <- sig2[n]
	a.all <- unique(unlist(dicts[[this.i]]))
	sorted.a.all <- names(sort(a.eff[a.all],decreasing=T))
	this.d <- dicts[[this.i]][[n]]
	sorted.this.d <- sorted.a.all[sorted.a.all %in% this.d]
	gs <- GSEA.EnrichmentScore(sorted.a.all,this.d,weighted.score.type=0)

	mem <- which(gs$indicator==1)
	le <- sorted.this.d[seq(mem[mem <= gs$arg.ES])]
	return(le)
})
names(dicts.atoc.pos) <- atoc.pos

## selective & essnetial
ois2 <- sapply(dicts.atoc.pos,function(x){
	sapply(dicts.atoc.pos,function(y){
		n1 <- length(x)
		n2 <- length(y)
		n3 <- length(intersect(x,y))
		oi <- n3/min(n1,n2)
		return(oi)
	})
})
rownames(ois2) <- colnames(ois2) <- atoc.pos

n.sc1 <- sapply(dicts.dmso.pos,length)
n.sc2 <- sapply(dicts.atoc.pos,length)

sc1 <- colorRampPalette(brewer.pal(9,"YlOrBr"))((28/2))[round((n.sc1)/2)])
sc2 <- colorRampPalette(brewer.pal(9,"YlOrBr"))((28/2))[round((n.sc2)/2)])
cols <- colorRampPalette(brewer.pal(9,"Blues"))(11)

##
setwd(plot_dir);setwd("tea")
pdf("dmso_pos.pdf",width=5,height=5) 
heatmap.2(ois1,trace="none",margins=c(12,12),cexRow=.4,cexCol=.4,
	col=cols,RowSideColors=sc1,ColSideColors=sc1,main="essential & selective")
dev.off()

pdf("atoc_pos.pdf",width=5,height=5) 
heatmap.2(ois2,trace="none",margins=c(12,12),cexRow=.4,cexCol=.4,
	col=cols,RowSideColors=sc2,ColSideColors=sc2,main="essential & selective")
dev.off()

pdf("index_sc.pdf",width=5,height=5) 
par(mar=c(5,5,2,2))
par(tcl=0)
image(1:14,1,matrix(1:14,ncol=1),
	col=colorRampPalette(brewer.pal(9,"YlOrBr"))(14),
	main="",axes=F,xlab="",ylab="")
box()
axis(1,at=c(1,7.5,14),labels=c(0,14,28))
dev.off()

##
setwd(rda_dir)
save(dicts.dmso.pos,dicts.atoc.pos,file="dicts_sig-tgts_le.rda")

dmso.pos.pro <- list(d1=c("ADA2A","ADA2B","ADA2C"),
	d2=c("AK1C3","STS","DHB1", "ESR1"),
	d3="NR1H3",
	d4="CP51A")
atoc.pos.pro <- list(a1="ALDH2",a2="C11B1",a3="CP51A",a4="PGTB1")

dmso.pos.grp <- lapply(dmso.pos.pro,function(x)unique(unlist(dicts.dmso.pos[x])))
atoc.pos.grp <- lapply(atoc.pos.pro,function(x)unique(unlist(dicts.atoc.pos[x])))


setwd(rda_dir)
load(file="dicts_sig_tgts.rda") # dicts.dmso.pos,dicts.atoc.pos,
dicts <- readRDS("dicts_MS-ecfp.rds")[-1]

setwd(rda_dir)
eff <- readRDS(file="MS-eff.rds")
cmpds <- sub("X","",rownames(eff))

##
dmso.pos.pro <- list(d1=c("ADA2A","ADA2B","ADA2C"),
	d2=c("AK1C3","STS","DHB1", "ESR1"),
	d3="NR1H3",
	d4="CP51A")
atoc.pos.pro <- list(a1="ALDH2",a2="C11B1",a3="CP51A",a4="PGTB1")

dmso.pos.grp <- lapply(dmso.pos.pro,function(x)unique(unlist(dicts.dmso.pos[x])))
atoc.pos.grp <- lapply(atoc.pos.pro,function(x)unique(unlist(dicts.atoc.pos[x])))

##
dc <- densCols(eff)

dicts.cmpds.all <- unique(unlist(dicts))
dc[!cmpds %in% dicts.cmpds.all] <- "grey80"

## DMSO.pos
setwd(plot_dir);setwd("eff_sig_cmpds")
for(i in 1:4){
	idx <- cmpds %in% dmso.pos.grp[[i]]
	pdf(paste0("CIL56_screening_dmso_",i,".pdf"),width=5,height=5)
	par(mar=c(4,4,2,2))
	plot(eff,pch=20,col=dc,asp=1,
		xlab="dAUC (in DMSO)",ylab="dAUC (in ATOC)",
		main="Summary of CIL56-enhancer/suppressor screening")
	abline(h=0,v=0)
	abline(a=0,b=1,lty=2)
	points(eff$d[dc=="grey80"],eff$a[dc=="grey80"],pch=20,col="grey80")
	points(eff$d[dc!="grey80" & !idx],eff$a[dc!="grey80" & !idx],pch=20,
		col=dc[dc!="grey80" & !idx])
	points(eff$d[idx],eff$a[idx],pch=20,col=2)
	dev.off()
}

## ATOC.pos
setwd(plot_dir);setwd("eff_sig_cmpds")
for(i in 1:4){
	idx <- cmpds %in% atoc.pos.grp[[i]]
	pdf(paste0("CIL56_screening_atoc_",i,".pdf"),width=5,height=5)
	par(mar=c(4,4,2,2))
	plot(eff,pch=20,col=dc,asp=1,
		xlab="dAUC (in DMSO)",ylab="dAUC (in ATOC)",
		main="Summary of CIL56-enhancer/suppressor screening")
	abline(h=0,v=0)
	abline(a=0,b=1,lty=2)
	points(eff$d[dc=="grey80"],eff$a[dc=="grey80"],pch=20,col="grey80")
	points(eff$d[dc!="grey80" & !idx],eff$a[dc!="grey80" & !idx],pch=20,
		col=dc[dc!="grey80" & !idx])
	points(eff$d[idx],eff$a[idx],pch=20,col=2)
	dev.off()
}
