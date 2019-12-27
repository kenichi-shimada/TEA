setwd(rda_dir)
load(file="dicts_sig-tgts_le.rda") # dicts.dmso.pos,dicts.atoc.pos,setwd(rda_dir)

##
dmso.pos.pro <- list(d1="CP51A",
	d2=c("AK1C3","STS","DHB1", "ESR1"),
	d3=c("ADA2A","ADA2B","ADA2C"),
	d4="NR1H3")
atoc.pos.pro <- list(a1="ALDH2",a2=c("C11B1","CP51A","PGTB1"))

dmso.pos.grp <- lapply(dmso.pos.pro,function(x){
	tab <- table(unlist(dicts.dmso.pos[x]))
	n <- names(tab[tab==max(tab)])
	return(n)
})
atoc.pos.grp <- lapply(atoc.pos.pro,function(x){
	tab <- table(unlist(dicts.atoc.pos[x]))
	n <- names(tab[tab==max(tab)])
	return(n)
})

dmso.titles <- c("CP51A","AK1C3/STS/DHB1/ESR1","ADA2A/B/C","NR1H3")
atoc.titles <- c("ALDH2","C11B1/CP51A/PGTB1")
all.titles <- c("CP51A/C11B1/PGTB1","AK1C3/STS/DHB1/ESR1","ADA2A/B/C","NR1H3","ALDH2")

## 
eff <- readRDS(file="MS-eff.rds")
cmpds <- sub("X","",rownames(eff))

eff.d <- eff$dmso
names(eff.d) <- cmpds

lsts <- lapply(dmso.pos.grp,function(xs){
	names(head(sort(eff.d[xs],decreasing=T),2))
})

dmso.sig <- unlist(lsts)

##
eff.a <- eff$atoc
names(eff.a) <- cmpds

lsts <- lapply(atoc.pos.grp,function(xs){
	names(head(sort(eff.a[xs],decreasing=T),2))
})

atoc.sig <- unlist(lsts)

sig.cmpds <- c(dmso.sig[1:2],atoc.sig[4],
	dmso.sig[3:8],atoc.sig[1:2])

sig.pchs <- rep(21:25,c(3,2,2,2,2))


## DMSO.pos
setwd(plot_dir);setwd("eff_sig_cmpds")
pdf("CIL56_screening_targets.pdf",width=5,height=5)
par(mar=c(4,4,2,2))

dc <- densCols(eff)
idx <- cmpds %in% sig.cmpds

plot(eff,pch=20,col=dc,asp=1,
	xlab="dAUC (in DMSO)",ylab="dAUC (in ATOC)",
	main="Modulators - predicted targets")
abline(h=0,v=0)
abline(a=0,b=1,lty=2)
points(eff$d[!idx],eff$a[!idx],pch=20,col=dc[!idx])
points(eff$d[match(sig.cmpds,cmpds)],eff$a[match(sig.cmpds,cmpds)],
	col=NA,pch=sig.pchs,bg=2)
legend("topleft",all.titles,pch=21:25,col=NA,pt.bg=2,bg="white",cex=.8)
dev.off()

##
tapply(sig.cmpds,rep(1:5,c(3,2,2,2,2)),identity)


xy1 <- locator(10,type="l")
xy2 <- locator(15,type="l")
xy3 <- locator(10,type="l")


eff <- readRDS(file="MS-eff.rds")
x <- load("MS-experimental-data.RData")

setwd(plot_dir)
pdf("sum_hist_eff.pdf",width=5,height=5)
hist(eff$d,breaks=100,border=NA,col="grey70",
	main="Effect of modulator compounds",
	xlab="dAUC",ylab="# compounds")
abline(h=0,v=0)
dev.off()


mats <- data.frame(array(NA,c(0,4))) ##not normalized by 0%
for(mo in unique(layout$plate.384)){
  ind <- match(mo,unique(layout$plate.384))
  thisa <- anno[anno$mother==mo & anno$cell=="HT1080",]
  dmso.c <- unique(thisa$CIL[thisa$mod=="DMSO"])
  atoc.c <- unique(thisa$CIL[thisa$mod=="ATOC"])
  tl <- layout[layout$plate.384==mo & layout$cmpd=="DMSO",]
  rows <- match(tl$row.384,LETTERS[1:16])
  cols <- tl$column.384
  vi <- sapply(thisa$b,function(pl){
    sapply(seq(nrow(tl)),function(i)norms[[pl]][rows[i],cols[i]])
  })
  factor.k1 <- factor(rep(
    apply(thisa[match(colnames(vi),thisa$b),c(2,4)],1,paste,collapse="."),each=nrow(vi)))
  levels(factor.k1) <- paste0("k",c(5,1,2,3,6,4,7))
  # levels(factor.k1) <- paste0("k",1:7)

  factor.k2 <- factor(rep(thisa$modulator[match(colnames(vi),thisa$b)],each=nrow(vi)))
  mod <- sub("(.).+","\\1",factor.k2)
  factor.t <- factor(paste0(paste("t",1:13,sep="")[as.numeric(sub(".+(..)$","\\1",mo))],mod))
  df <- data.frame(vi=as.vector(vi),k1=factor.k1,k2=factor.k2,t=factor.t)
  mats <-rbind(mats,df)
}

##
if(0){
  vi1s <- sapply(levels(mats$t),function(fa)quantile(mats$vi[mats$t==fa&mats$k1=="k1"],1))## max is t12
  ##c0.max <- quantile(mats$vi[mats$t=="t12"&mats$k1=="k1"],seq(0,1,.1))
  c0.max <- 60000 ## from here, set as 60000
}

##
mats.dmso <- mats %>% filter(k2=="DMSO") %>% 
  mutate(t=factor(t)) %>%
  group_by(k1,t) %>% summarize(vi=median(vi)) %>%
  arrange(t,k1)
mats.atoc <- mats %>% filter(k2=="ATOC") %>% 
  mutate(t=factor(t)) %>%
  group_by(k1,t) %>% summarize(vi=median(vi)) %>%
  arrange(t,k1)

dmso.vi <- tapply(mats.dmso$vi,mats.dmso$t,identity)
atoc.vi <- tapply(mats.atoc$vi,mats.atoc$t,identity)

##
setwd(plot_dir)
pdf("CIL56_batch-eff.pdf",width=7,height=5)
par(mar=c(4,4,2,7))
plot(c(.31,2.5),0:1,type="n",
  xlab="[CIL56] (ug/mL)",ylab= "Normalized viability",log="x",
  axes=F,main="CIL56 efficacy on different days")
box()
axis(1,at=2.5*(1/2)^(3:0))
axis(2)
abline(v=1.25,lty=2,col="grey70")
cols <- colorRampPalette(brewer.pal(11,"Spectral"))(13)
for(i in 1:13){
  lines(1.25*(1/2)^(2:0),dmso.vi[[i]][-1],
    col=cols[i],
    lwd=1)
  points(1.25*(1/2)^(2:0),dmso.vi[[i]][-1],
    col=cols[i],
    pch=20,cex=1)
}
for(i in 1:13){
  lines(2.5*(1/2)^(1:0),atoc.vi[[i]][-1],
    col=cols[i],
    lwd=1)
  points(2.5*(1/2)^(1:0),atoc.vi[[i]][-1],
    col=cols[i],
    pch=20,cex=1)
}
par(xpd=T)
legend(10^par()$usr[2],par()$usr[4],paste("t",1:13),
  lty=1,pch=20,cex=1,col=cols)
par(xpd=F)
dev.off()
