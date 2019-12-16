interp.mat <-
function(data,rn,cn,col=rev(brewer.pal(11,"Spectral")),row.log=TRUE,cmpd){
  require(akima)
  require(RColorBrewer)
  stopifnot(is.matrix(data)|is.data.frame(data))
  if(is.data.frame(data)) data <- as.matrix(data)
  conv <- t(data[rev(seq(nrow(data))),])
  if(exists("cn")){ ##col -> x-axis
    x <- rep(cn,times=ncol(conv))
  }else{
    x <- rep(seq(nrow(conv)),times=ncol(conv))
  }
  if(exists("rn")){
    lrn <- if(row.log==TRUE){log10(rn)}else{rn}
    y <- rep(rev(lrn),each=nrow(conv))
  }else{
    lrn <- rn <- rev(seq(ncol(conv)))
    y <- rep(rev(lrn),each=nrow(conv))
  }
  ip <- interp(x,y,as.vector(conv))
  if(exists("cmpd")){
    ylab <- paste("[",cmpd,"](ug/ml)",sep="")
  }else{
    ylab <- "[cmpd](ug/ml)"
  }
  par(mar=c(5,8,3,1)+.1)
  ###levels
  std <- data[nrow(data),1]
  levels <- c(seq(0,std,length.out=6),seq(std,max(data),length.out=7)[-1])
  norm.levels <- round(levels/std *100)
  filled.contour(ip$x,ip$y,ip$z,col=col,levels=levels,
                 plot.axes=c(axis(1,at=time,labels=round(time),cex.axis=1),
                   axis(2,at=lrn,labels=round(rn,4),las=1,cex.axis=1)),
                 key.title="viability\n(%)",
                 key.axes=axis(4,at=levels,labels=norm.levels),
                 xlab=mtext("lapse(hour)",1,line=3.5))
  mtext(ylab,2,line=5,cex=1)
  return()
}

