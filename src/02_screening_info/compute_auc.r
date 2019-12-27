library(dplyr)

setwd(rda_dir)
x <- load("MS-experimental-data.RData")
norms <- readRDS("normalized_adjusted.plates.rds")

## finally, compute dAUC
AUC <- function(x,na.rm=T,normalize=T){
  n <- length(x)
  stopifnot(n>1)
  if(normalize){
    return(sum(x,x[-c(1,n)],na.rm=T)/(2*(n-1),
  }else{
    return(sum(x,x[-c(1,n)],na.rm=T))
  }
}

##
# compound -> 
#	mother, well (layout) -> 
#	barcodes for each cil, ms, modulator ->
#	get the normalized viability

# cond 1: CIL56 == 0 (MS only) -> filter out mods that are lethal (0.5,1,2,4uM)
# cond 2: CIL56 + MS divided by modulators, same conc (cell)
# cond 3: with or without MS - is it positive or not

anno <- anno %>% filter(anno$barcode %in% names(norms))

MS.tox <- function(cmpd,mod=c("DMSO","ATOC")){
	mo <- layout$plate.384[layout$cmpd==cmpd]
	row <- match(layout$row.384[layout$cmpd==cmpd],LETTERS[1:16])
	col <- layout$column.384[layout$cmpd==cmpd]
	mod <- mod[1]

	tanno1 <- anno %>% 
		filter(mother==mo & modulator==mod & cell=="HT1080" & CIL56 ==0)
	cmpd.conc <- unique(tanno1$modulator)
	bcs.list <- tapply(tanno1$barcode,tanno1$MS,function(x){
		vis <- sapply(norms[x],function(n)n[row,col])
		med.vi <- median(vis)
		return(med.vi)
	})
	return(bcs.list)
}

cmpds <- (layout %>% 
	filter(!cmpd %in% "DMSO" & !grepl("unknown",cmpd) & !grepl("new",cmpd)))$cmpd # 1941
mos <- as.numeric(factor((layout %>% filter(cmpd %in% cmpds))$plate.384))

MS.dmso <- lapply(cmpds,MS.tox,mod="DMSO")
MS.atoc <- lapply(cmpds,MS.tox,mod="ATOC")

names(MS.dmso) <- names(MS.atoc) <- cmpds
table(sapply(MS.dmso,length)) # 160 have 1, 1814 have 4
table(sapply(MS.atoc,length)) # 1974 have 4

correct.format <- sapply(MS.dmso,length)==4

AUC.dmso <- sapply(MS.dmso[correct.format],AUC)
AUC.atoc <- sapply(MS.atoc[correct.format],AUC)

##
dc <- densCols(AUC.dmso,AUC.atoc)
idx <- AUC.dmso > .8 & AUC.atoc > .8 & AUC.dmso + AUC.atoc > 1.8
dc[!idx] <- "grey70"

pdf("AUC_MS-only.pdf",width=5,height=5)

par(mar=c(4,4,2,2))
plot(AUC.dmso,AUC.atoc,col=dc,pch=20,asp=1,
	xlab="AUC (with DMSO)",ylab="AUC (with ATOC)")
abline(h=0,v=0,lty=1)
abline(h=c(1),v=c(1),lty=2)
abline(a=1.8,b=-1,col=4)
points(AUC.dmso[fer],AUC.atoc[fer],pch=20,col=2)
points(AUC.dmso["1504233"],AUC.atoc["1504233"],pch=20,col=2)
text(AUC.dmso[fer],AUC.atoc[fer],fer,pos=3)

dev.off()

##
xy <- locator(10,type="l")

library(sp)
fer <- names(AUC.dmso)[which(point.in.polygon(AUC.dmso,AUC.atoc,xy$x,xy$y)==1)]

fers <- lapply(fer,function(cmpd){
	lst <- list(d=MS.dmso[[cmpd]],a=MS.atoc[[cmpd]])
})
names(fers) <- c("MS00300007", "MS00201522", "MS001503913",)

x.mains <- c("euparin","ganbogic acid amide","nigericin sodium")
names(x.mains) <- fer
xs <- c(0.5,1,2,4)

##
setwd(plot_dir)
for(x in fer){
	pdf(paste0(x,".pdf"),width=4,height=3)
	par(mar=c(4,4,2,2))
	tmp <- fers[[x]]
	plot(c(0.5,4),c(0,1.2),log="x",type="n",
		main=x.mains[x],xlab=paste0("[",x,"] (uM)"),
		ylab="Normalized AlamarBlue",axes=F)
	box()
	axis(1,at=xs)
	axis(2)
	abline(h=0:1,lty=1:2)
	lines(xs,tmp$d,col=1)
	points(xs,tmp$d,pch=20,col=1)
	lines(xs,tmp$a,col=2)
	points(xs,tmp$a,pch=20,col=2)
	dev.off()
}

## 
mos <- unique(layout$plate.384)
vi.vehs <- t(sapply(mos,function(mo){
	bcs <- anno$barcode[anno$mother==mo]
	wells <- layout[layout$plate.384==mo & layout$cmpd=="DMSO",3:4]
	rows <- match(wells$r,LETTERS[1:16])
	cols <- wells$c
	vi.veh <- sapply(c("DMSO","ATOC"),function(mod){
		tanno1 <- anno %>% filter(mother==mo & modulator==mod)
		tmp <- tapply(tanno1$barcode,tanno1$CIL56,function(bcs){
			median(sapply(norms[bcs],function(n)median(n[rows,cols])))
		})
		tmp <- tmp[-1]/tmp[1]
	})
	names(vi.veh) <- c("DMSO","ATOC")
	aucs <- sapply(vi.veh,AUC)
	return(aucs)
}))

mod.eff <- function(cmpd,mod=c("DMSO","ATOC"),normalize=T,std=F){
	mo <- layout$plate.384[layout$cmpd==cmpd]
	row <- match(layout$row.384[layout$cmpd==cmpd],LETTERS[1:16])
	col <- layout$column.384[layout$cmpd==cmpd]
	mod <- mod[1]

	##
	tanno1 <- anno[anno$mother==mo & anno$mod==mod & anno$cell=="HT1080",]

	vis <- tapply(tanno1$barcode,list(tanno1$CIL56,tanno1$MS),function(bcs){
		median(sapply(norms[bcs],function(x)x[row,col]))
	})

	names(dimnames(vis)) <- c("CIL56","MS")
	n.vis <- (vis/rep(vis[1,],each=nrow(vis)))[-1,]
	aucs <- apply(n.vis,2,AUC)
	d.aucs <- aucs-vi.vehs[mo,mod]
	return(d.aucs)
}

used.cmpds <- cmpds[correct.format][idx]

eff.dmso <- t(sapply(used.cmpds,mod.eff,mod="DMSO"))
eff.atoc <- t(sapply(used.cmpds,mod.eff,mod="ATOC"))

eff <- cbind(dmso=eff.dmso[,4],atoc=eff.atoc[,4])
boxplot(lapply(1:4,function(i)eff.dmso[,i]))

dc <- densCols(eff)

##
eff <- data.frame(eff)
dap <- which(eff$d > .6 & eff$a > .6)
dp <- which(eff$d > .6 & eff$a < .1)
ap <- which(eff$d < 0 & eff$a > .35)
dan <- which(eff$d < -.3 & eff$a < -.3)
dn <- which(eff$d < -.4 & eff$a > -.1)
an <- which(eff$d > -.1 & eff$d < 0 & eff$a < -.33)

keys.i <- c(dap,dp,ap,dan,dn,an)
keys <- lapply(list(dap,dp,ap,dan,dn,an),function(x)used.cmpds[x])

##
setwd(plot_dir)
pdf("CIL56_screening.pdf",width=5,height=5)
par(mar=c(4,4,2,2))
plot(eff,pch=20,col=dc,asp=1,
	xlab="dAUC (in DMSO)",ylab="dAUC (in ATOC)",
	main="Summary of CIL56-enhancer/suppressor screening")
abline(h=0,v=0)
abline(a=0,b=1,lty=2)
points(eff$d[keys.i],eff$a[keys.i],pch=20,col=2)
text(eff$d[keys.i],eff$a[keys.i],sub("X","",rownames(eff))[keys.i],pos=3)
dev.off()

setwd(rda_dir)
saveRDS(eff,file="MS-eff.rds")

##
setwd(plot_dir)
dc <- densCols(eff)
  
x <- eff$d
y <- eff$a

rng <- range(x,y)
xhist <- hist(x, breaks = seq(rng[1],rng[2],length=100), plot = FALSE)
yhist <- hist(y, breaks = seq(rng[1],rng[2],length=100), plot = FALSE)
top <- max(c(xhist$counts, yhist$counts))
xrange <- yrange <- rng

pdf("CIL56_screening_histo.pdf",width=5,height=5)
nf <- layout(matrix(c(2,0,1,3),2,2,byrow = TRUE), c(5,1), c(1,5), TRUE)
# layout.show(nf)

par(mar = c(4,4,1,1))
plot(eff,pch=20,col=dc,asp=1,
	xlab="dAUC (in DMSO)",ylab="dAUC (in ATOC)",
	xlim=xrange,ylim=yrange,
	main="")
abline(h=0,v=0)
abline(a=0,b=1,lty=2)
par(mar = c(0,4,1,1))
barplot(xhist$counts, axes = FALSE, ylim = c(0, top), space = 0,)
par(mar = c(4,0,1,1))
barplot(yhist$counts, axes = FALSE, xlim = c(0, top), space = 0, horiz = TRUE)
     
dev.off()

setwd(rda_dir)
saveRDS(eff,file="MS-eff.rds")