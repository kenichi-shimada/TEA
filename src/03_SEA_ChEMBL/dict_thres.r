library(dplyr)
setwd("~/Dropbox (HMS)/projects/02 TEA/204_SEA_ChEMBL/dicts/")
load(file="MS-ecfp-all.rda") ##93535 relationships are e.value <= 1

dict <- dict %>% 
	filter(grepl("^cb_.+_HUMAN",code)) %>%
	filter(grepl("^cb_(.+)_HUMAN",code)) %>%
	mutate(target=sub("^cb_(.+)_HUMAN","\\1",code))

evalue <- dict$E.value
maxtc <- dict$MaxTC

dict$cmpd <- sub("ms_","",dict$Compound.ID) 
evalue[evalue==0] <- min(evalue[evalue!=0])
range(evalue)

x <- -log10(evalue)
y <- maxtc

idx <- -log10(evalue) >= -Inf

x <- x[idx]
y <- y[idx]

xrange <- range(x)
yrange <- range(y)

xhist <- hist(x, breaks = 100,plot=F)
yhist <- hist(y, breaks = 100,plot=F)
top <- max(c(xhist$counts, yhist$counts))

##
setwd(plot_dir)
png("evalue_maxtc.png",width=700,height=700)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(4,1), c(1,4), TRUE)
# layout.show(nf)
cols <- densCols(x,y)

par(mar = c(4,4,1,1))
plot(x, y, xlim = xrange, ylim = yrange, col=cols,pch=20,cex=.2,
	xlab = "-log10(E-value)", ylab = "MaxTC")

par(mar = c(0,4,1,1))
barplot(xhist$counts, axes = FALSE, ylim = c(0, top), space = 0)

par(mar = c(4,0,1,1))
barplot(yhist$counts, axes = FALSE, xlim = c(0, top), space = 0, horiz = TRUE)

dev.off()

##
setwd(plot_dir)
pdf("ecdf_evalue.pdf",width=5,height=5)
par(mar=c(4,4,2,2))
plot.ecdf(-log10(evalue),xlab="-log10(E.value)",
	ylab="empirical CDF",main="")
abline(v=seq(0,100,10),col="grey70")

nl <- length(evalue)
xs <- seq(0,100,10)
ys <- 1-sapply(xs,function(x)sum(-log10(evalue)>x))/nl
points(xs,ys,pch=20)

dev.off()

ns <- t(sapply(xs,function(x){
	idx <- -log10(evalue)>x
	n.all <- sum(idx)
	n.cmpd <- length(unique(dict$Compound.ID[idx]))
	n.tgt <- length(unique(dict$code[idx]))
	return(c(n.all=n.all,n.cmpd=n.cmpd,n.tgt=n.tgt))
}))
write.csv(ns,file="stats_dict_thres.csv")


##
setwd(rda_dir <- "~/Dropbox (HMS)/projects/02 TEA/TEA/rdas")
eff <- readRDS("MS-eff.rds")
cmpds <- sub("X","",rownames(eff))

dict <- dict[dict$cmpd %in% cmpds,] # 284680

library(dplyr)
dicts <- lapply(c(-Inf,seq(0,100,10)),function(th){
	dict.1 <- dict %>% filter(-log10(E.value) > th)
	grps <- tapply(dict.1$Compound.ID,dict.1$target,function(x)sub("ms_","",x))
	is.good.size <- sapply(grps,length) >= 10 & sapply(grps,length) <= 200
	return(grps[is.good.size])
})
this.n <- sapply(dicts[[1]],length)
sapply(dicts,function(x)grep("LOX",names(x),value=T))

setwd(plot_dir)
pdf("n.targets.pdf",width=4,height=5)
par(mar=c(7,5,4,1))
ns <- sapply(dicts,length)
names(ns) <- c("all","e < 1",paste0("e < 1e-",seq(10,100,10)))
barplot(ns,las=2,ylab="# targets",main="")
dev.off()
##
setwd(rda_dir)
saveRDS(dicts,file="dicts_MS-ecfp.rds")

