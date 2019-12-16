plot.dc <-
function(x,cmpd,conc=20*(1/2)^(seq(nrow(x))-1),
           labels=unique(colnames(x)),unit="ug/ml",cond=NA,main="",
           xlab=paste("[",cmpd,"](",unit,")",sep=""),ylab="Normalized AlamarBlue",
           leg.x=NA,leg.y=NA,pcol=seq(length(labels)),lcol=seq(length(labels)),...){
  old.par <- par(no.readonly = TRUE); on.exit(old.par)
  par(mar=c(10,6,1,1)+.1)
  plot(range(conc),c(0,1.1),type="n",axes=F,log="x",xlab="",ylab="")
  box()
  axis(1,at=round(conc,3),las=2,cex.axis=2.5)
  axis(2,at=seq(0,1,.2),cex.axis=2.5)
  mtext(xlab,1,line=8,cex=3)
  mtext(ylab,2,line=4,cex=3)
  abline(h=0,col="grey50")
  abline(h=1,col="grey50")
  ##points
  for(j in seq(ncol(x))){
    this <- x[,j]
    color <- pcol[match(colnames(x)[j],labels)]
    points(conc,this,col=color,pch=20)
  }
  ##lines
  for(k in labels){
    this <- apply(x[,which(!is.na(match(colnames(x),k)))],1,median,na.rm=T)
    color <- lcol[match(k,labels)]
    lines(conc,this,col=color,lwd=3)
  }
  if(is.na(leg.x))leg.x <- 10^par()$usr[1]
  if(is.na(leg.y))leg.y <- 0.4
  legend(leg.x,leg.y,labels,lty=1,lwd=3,pch=20,col=lcol,cex=2,...)
}

