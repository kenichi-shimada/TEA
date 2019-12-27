library(RColorBrewer)
library(gplots)

setwd(rda_dir)
x <- load(file="gsea_data.rda") # cmpds,eff,dicts,list.cmpds

## in O2
setwd(rda_dir);setwd("fgsea")
fns <- paste0("fgsea_1e7_",1:24,".rds")[c(2:12,14:24)]

setting <- expand.grid(thres=2:12,mod=c("dmso","atoc"),stringsAsFactors=F)

##
slps.all  <- list()
for(i in 1:22){
	fn <- fns[[i]]
	res <- readRDS(fn)
	slps.all[[i]] <- sign(res$ES) * -log10(res$pval)
}

dicts <- dicts[-1]

tgts <- unique(unlist(sapply(dicts,names))) # 502
slp.thres <- 3
slp.dmso <- sapply(1:11,function(i){
	x <- slps.all[[i]]
	this.tgts <- names(dicts[[i]])
	this.slp <- rep(NA,length(tgts))
	names(this.slp) <- tgts
	this.slp[this.tgts] <- x
	return(this.slp)
}) 
dmso.pos <- names(which(apply(slp.dmso[,],1,max,na.rm=T) > slp.thres))
dmso.neg <- names(which(apply(slp.dmso[,],1,min,na.rm=T) < -slp.thres))
# dmso.pos <- "CALM"

##
slp.atoc <- sapply(1:11,function(i){
	x <- slps.all[[i+11]]
	this.tgts <- names(dicts[[i]])
	this.slp <- rep(NA,length(tgts))
	names(this.slp) <- tgts
	this.slp[this.tgts] <- x
	return(this.slp)
}) 
atoc.pos <- names(which(apply(slp.atoc[,],1,max,na.rm=T) > slp.thres))
atoc.neg <- names(which(apply(slp.atoc[,],1,min,na.rm=T) < -slp.thres))


##
setwd(rda_dir)
save(slp.dmso,slp.atoc,file="slps_tgts.rda")
save(dmso.pos,dmso.neg,atoc.pos,atoc.neg,file="sig_tgts.rda")

##
labs <- c("1",paste0("1e-",seq(10,100,10)))

##
setwd(plot_dir);setwd("tea")
# pdf("slp_dmso_scan.pdf",width=5,height=4)
pdf("slp_dmso_scan_1.pdf",width=5,height=4)
par(mar=c(4,4,2,2))
plot(c(1,11),range(slp.dmso,na.rm=T),type="n",main="Ferroptosis",
	xlab="threshold (E < )",ylab="Signed log P-value",axes=F)
box()
axis(1,at=1:11,labels=labs,las=2)
axis(2)
uniq.cols <- c("grey70",2,4)
cols <- rep(uniq.cols[1],nrow(slp.dmso))
cols[tgts %in% dmso.pos] <- uniq.cols[2]
cols[tgts %in% dmso.neg] <- uniq.cols[3]

for(j in uniq.cols){
	idx <- which(cols==j)
	for(i in idx){
		non.na <- which(!is.na(slp.dmso[i,]))
		lines(non.na,slp.dmso[i,non.na],col=j)
		if(length(non.na)==1){
			points(non.na,slp.dmso[i,non.na],col=j,pch=20,cex=.5)
		}
		if(j != uniq.cols[1]){
			text(max(non.na),slp.dmso[i,max(non.na)],tgts[i],pch=20,cex=.5)
		}
	}
}
abline(h=c(3,-3),lty=2)
abline(h=0)
dev.off()


# pdf("slp_atoc_scan.pdf",width=5,height=4)
pdf("slp_atoc_scan_1.pdf",width=5,height=4)
par(mar=c(4,4,2,2))
plot(c(1,11),range(slp.atoc,na.rm=T),type="n",main="Necrosis",
	xlab="threshold (E < )",ylab="Signed log P-value",axes=F)
box()
axis(1,at=1:11,labels=labs,las=2)
axis(2)
uniq.cols <- c("grey70",2,4)
cols <- rep(uniq.cols[1],nrow(slp.atoc))
cols[tgts %in% atoc.pos] <- uniq.cols[2]
cols[tgts %in% atoc.neg] <- uniq.cols[3]

for(j in uniq.cols){
	idx <- which(cols==j)
	for(i in idx){
		non.na <- which(!is.na(slp.atoc[i,]))
		lines(non.na,slp.atoc[i,non.na],col=j)
		if(length(non.na)==1){
			points(non.na,slp.atoc[i,non.na],col=j,pch=20,cex=.5)
		}
		if(j != uniq.cols[1]){
			text(max(non.na),slp.atoc[i,max(non.na)],tgts[i],pch=20,cex=.5)
		}
	}
}
abline(h=c(3,-3),lty=2)
abline(h=0)
dev.off()

#
xd <- range(slp.dmso,na.rm=T)
xa <- range(slp.atoc,na.rm=T)

titles <- c("E < 1",paste0("E < 1e-",seq(10,100,10)))
for(i in c(1,4,7,10)){
	setwd(plot_dir);setwd("tea")
	pdf(paste0("scatter_",i,".pdf"),width=5,height=5)
	par(mar=c(4,4,2,2))
	slps <- data.frame(do.call(cbind,slps.all[c(0,11)+i]))
	names(slps) <- c("DMSO","ATOC")
	rownames(slps) <- tgts <- names(dicts[[i]])
	# slps <- slps[tgts,]

	par(mar=c(4,4,1,1))
	dc <- densCols(slps)
	plot(slps$DMSO,slps$ATOC,pch=20,
		xlab="Signed log-P (ferroptosis)",
		ylab="Signed log-P (necrosis) ",asp=1,
		xlim=xd,ylim=xa,
		main=titles[i],
		col=dc)
	# points(slps$DMSO,slps$ATOC,pch=20,cex=.3)
	# idx <- grep("^LOX",tgts,ignore.case=T,value=F)
	idx <- slps$DMSO >= 3 | slps$ATOC >= 3
	points(slps$DMSO[idx],slps$ATOC[idx],pch=20,cex=1,col=2)
	text(slps$DMSO[idx],slps$ATOC[idx],tgts[idx],pos=3)
	abline(h=0,v=0)
	abline(h=3,v=3,lty=2)
	dev.off()
}




