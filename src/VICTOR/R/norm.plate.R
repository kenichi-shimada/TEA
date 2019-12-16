norm.plate <-
function(x,row=2:15,bg.col=c(2,23),pos.col=c(3,22),method=c("median")){
  if(is.list(x)){
    norm <- as.list(NULL)
    for(i in seq(x)){
      tmp <- x[[i]]
      bg <- median(tmp[row,bg.col])
      pos <- median(tmp[row,pos.col])
      tmp.norm <- (tmp-bg)/(pos-bg)
      norm[[names(x)[i]]] <- tmp.norm
    }
  }else if(is.matrix(x)){
    tmp <- x
    bg <- median(tmp[row,bg.col])
    pos <- median(tmp[row,pos.col])
    tmp.norm <- (tmp-bg)/(pos-bg)
    norm <- tmp.norm
  }
  return(norm)
}

