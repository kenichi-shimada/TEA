library(RColorBrewer)
library(VICTOR)

setwd(rda_dir);
x <- load("MS-experimental-data.RData")

## visualize plates
setwd(plot_dir)
pdf("738-plates.pdf",width=11,height=8.5)
a <- sapply(names(plates$r),function(x)plate.image(plates$r[[x]],main=x))
dev.off()

bcs <- tapply(anno$barcode,anno$mother,identity)
summed.pl <- lapply(bcs,function(bc){
	Reduce('+', plates$r[bc])
})

setwd(plot_dir)
pdf("summary_13-mothers.pdf",width=11,height=8.5)
x <- sapply(names(summed.pl),function(x)plate.image(summed.pl[[x]],main=x))
dev.off()

png("summary_13mothers.png",width=800,height=600)
par(mfrow=c(4,4))
par(mar=c(3,3,2,1))
x <- sapply(names(summed.pl),function(x)plate.image(summed.pl[[x]],main=x))
dev.off()

## assay plate layouts
plate <- array(NA,c(16,24)) # yellow
plate[,c(2,23)] <- 4 # media
plate[3:14,3:22] <- 3 # cells
plate[c(4,13),c(3,8,13,18,22)] <- 2 # lethal
plate[5:12,3:22] <- 1 # lethal+mod

##
setwd(plot_dir)
cols <- brewer.pal(4,"Accent")
png("assay_plate_layout.png",width=600,height=300)
par(mar = c(3, 3, 2, 12))
pl <- plate.image(plate,main="Assay plate layout",col=cols)
par(xpd=T)
legend(par()$usr[2]+.1,16,
	rev(c("media","cells (no lethal)","lethal only","lethal + modulator")),
	fill=cols)
par(xpd=F)
dev.off()

### checking raw values
coords <- list(
	media=which(plate==4),
	cells=which(plate==3),
	lethal=which(plate==2),
	lm=which(plate==1))

range.raws <-  t(sapply(raws,function(r){
	sapply(coords,function(i)quantile(r[i],.95))
})
plot(range(seq(nrow(range.raws))),range(range.raws),type="n",
	xlab="range of raw values",ylab=)

order.m <- unique(anno$mother[match(rownames(range.raws),anno$barcode)])
n.m <- table(anno$mother[match(rownames(range.raws),anno$barcode)])
positions <- c(0.5,cumsum(n.m[match(order.m,names(n.m))])+.5)
	abline(h=0,v=positions,col="grey50")
	mid.pos <- (positions[-1]+positions[-14])/2
text(mid.pos,10000,order.m,srt=90)
points(range.raws[,1],pch=".",col=4)
points(range.raws[,2],pch=".",col=2)


