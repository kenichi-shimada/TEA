read.victor <-
function(file,path=getwd(),time=TRUE,delimiter="\t",keyword=NULL,read=TRUE,
                        na.action=c("error","na.rm","ignore")){
  extfile <- paste(path,file,sep="/")
  tmp <- scan(extfile,what=character(0),sep=delimiter)
  conds <- c("Plate","Repeat","End time","Start temp.","End temp.","BarCode")

  bcs <- grep(conds[6],tmp) ##number can be different every time dependent on the choice of keyword/na.action
  na.action <- na.action[1] ##by default
  if(exists("read") && !read)return(tmp[bcs+24])
  if(!is.null(keyword)){
    bcs <- bcs[grep(keyword,tmp[bcs+24])] ##have to be modified
    names(bcs) <- tmp[bcs+24]
  }else{
    if(na.action=="error" && length(grep("N/A",tmp[bcs+24])) > 0){
      print("barcodes can't be N/A as they are used as identifiers.\n") ##check
      return(tmp[bcs+24])
    } else if(grep("na.rm",na.action) && length(grep("N/A",tmp[bcs+24])) > 0){
      bcs <- bcs[-grep("N/A",tmp[bcs+24])]
      names(bcs) <- tmp[bcs+24]
    } else{
      names(bcs) <- tmp[bcs+24]
    }
  }
  n.pl <- length(bcs)
  params <- data.frame(t(sapply(bcs,function(x)tmp[x+(19:24)])),stringsAsFactors=FALSE)
  params[,c(1,2,4,5)] <- sapply(c(1,2,4,5),function(x){as.numeric(params[,x])})
  colnames(params) <- conds
  names <- if(length(levels(params$BarCode))==length(bcs)){
    params$BarCode
  }else{
    seq(length(bcs))
  }
  params[,c(1,2,4,5)] <- sapply(c(1,2,4,5),function(x){as.numeric(params[,x])})
  if(time){
    time <- sub("^12:","0:",params$"End time")
    time <- sub(" AM",":0",time)
    time <- sub(" PM",":12",time)
    time <- t(sapply(time,function(x)as.numeric(strsplit(x,"[ :]")[[1]])))
    mins <- time %*% c(60,1,1/60,60)
    params$minutes <- mins
    rownames(time) <- c()
  }
  mats <- lapply(bcs,function(x)
                 matrix(as.numeric(tmp[x+(139:522)]),nrow=16,ncol=24,byrow=T))
  if(any(sapply(mats,function(x)any(is.na(x))))){
    mats <- lapply(bcs,function(x)
                   matrix(as.numeric(tmp[x+(139:538)]),nrow=16,ncol=25,byrow=T)[,-25])
  }
  class(mats) <- "VictorReadout"
  rownames(params) <- NULL
  return(list(params=params,readouts=mats))
}

