read.multi.victor <-
function(pattern="\\.txt",path=getwd(),...){
  files <- dir(path,pattern=pattern)
  names(files) <- files
  all <- lapply(files,read.victor,path=path,...)

  ##parameters
  all.params <- lapply(files,function(x)all[[x]]$params)
  params <- do.call("rbind",all.params)
  params$file <- sub("\\.[0-9]+$","",rownames(params))
  rownames(params) <- NULL 

  ##readouts
  all.ro <- mapply(function(x,y)all[[x]]$readouts[[y]],params$f,params$B,SIMPLIFY
                   =FALSE)
  names(all.ro) <- paste(params$f,params$B,sep=":")
  (list(params=params,readouts=all.ro))
}

