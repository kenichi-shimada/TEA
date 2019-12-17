plate.image <- function(plate,main="Plate Heatmap",col=heat.colors(256),...){
  par(tcl=0)
  image(1:24,1:16,t(plate[16:1,]),col=col,axes=F,xlab="",ylab="",
        main=main,asp=1)
  sapply(2:16,function(x)lines(c(0.5,24.5),rep(x-.5,2),lty=1,col="grey50"))
  sapply(2:24,function(x)lines(rep(x-.5,2),c(0.5,16.5),lty=1,col="grey50"))
  sapply(c(0.5,8.5,16.5),function(x)lines(c(0.5,24.5),rep(x,2),lty=1,col="grey20"))
  sapply(c(0.5,12.5,24.5),function(x)lines(rep(x,2),c(0.5,16.5),lty=1,col="grey20"))
  par(xpd=T)
  text(-0.5,16:1,LETTERS[seq(16)],cex.axis=.6,las=1)
  text(1:24,-0.5,seq(24),cex.axis=.5,las=1)
  par(xpd=F)
  invisible()
}

