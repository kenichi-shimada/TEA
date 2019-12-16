simple.heatmap <-
function(x,col.na="grey",row.names=T,col.names=T,
         row.clust=T,col.clust=T,clust.method="average",
         dist.method=c("pearson","spearman","euclidean","manhattan"),
         cols=c("RdGn","YlBl"),cex.x=.5,cex.y=.5,asp=1,
         x.adj=0,y.adj=0,index=F){
  ##
  if(any(row.clust,col.clust)){
    dist.method <- dist.method[1]
    if(dist.method %in% c("spearman","pearson")){
      if(row.clust){
        hc.row <- hclust(as.dist(1-cor(t(x),method=dist.method,use="na.or.complete")),"average")
        or <- hc.row$labels[hc.row$order]
      }else{
        or <- rownames(x)
      }
      if(col.clust){
        hc.col <- hclust(as.dist(1-cor(x,method=dist.method,use="na.or.complete")),"average")
        oc <- hc.col$labels[hc.col$order]
      }else{
        oc <- colnames(x)
      }
    }else if(dist.method %in% c("euclidean","manhattan")){
      if(row.clust){
        hc.row <- hclust(dist(x,method=dist.method),"average")
        or <- hc.row$labels[hc.row$order]
      }else{
        or <- rownames(x)
      }
      if(col.clust){
        hc.col <- hclust(dist(t(x),method=dist.method),"average")
        oc <- hc.col$labels[hc.col$order]
      }else{
        oc <- colnames(x)
      }
    }
  }else{
    or <- rownames(x)
    oc <- colnames(x)
  }
  ##
  cols <-  cols[1]
  if(cols=="YlBl"){
    colors <- c(rgb(0,0,256:1,maxColorValue=256),rgb(0:256,0:256,0,maxColorValue=256))
  }else if(cols=="RdGn"){
    colors <- c(rgb(0,256:1,0,maxColorValue=256),rgb(0:256,0,0,maxColorValue=256))
  }
  tmp <- x
  nr <- nrow(tmp); nc <- ncol(tmp)
  tmp[is.na(x)] <- min(x,na.rm=T)-(max(x,na.rm=T)-min(x,na.rm=T))/(length(colors)-1)
  if(any(is.na(x)))colors <- c(col.na,colors)
  if(!index){
    image(seq(nc),seq(nr),t(tmp[or,oc][rev(seq(nr)),]), axes=F,xlab="",ylab="",col=colors,add=F,asp=asp)
    ##
    if(any(row.names,col.names)){
      par(tcl=0)
      par(xpd=T)
      if(row.names){
        text(0,seq(nr)+x.adj,rev(or),cex=cex.x,pos=2)
      }
      if(col.names){
        text(seq(nc)+y.adj,0,oc,cex=cex.y,pos=2,srt=90)
      }
      par(xpd=F)
    }
  }else{
    image(seq(colors),1,array(seq(colors),c(length(colors),1)),col=colors,axes=F,xlab="",ylab="",asp=10)
    par(tcl=-0.5)
    box()
    axis(1,at=c(min(seq(colors)),median(seq(colors)),max(seq(colors))),
         labels=c(paste("< ",min(x,na.rm=T),sep=""),
           (min(x,na.rm=T)+max(x,na.rm=T))/2,
           c(paste("> ",max(x,na.rm=T),sep=""))))
  }
}

