### compiling dictionary file from spreadsheet (exported to txt) to R objects.
setwd("~/Dropbox (HMS)/projects/02 TEA/204_SEA_ChEMBL/dicts/txts_MS")
# txts <- dir(pattern="path\\.txt$")
txts <- dir(pattern="ecfp\\.txt$")
types <- sub("_ecfp.txt","",txts)
dict <- data.frame(row.names=c("Compound.ID","code","Target.name","e.value","MaxTC"))

for(ind in 1:5){
  txt <- txts[ind]
  type <- types[ind]
  db <- t(sapply(readLines(txt),function(x)strsplit(x,split="\t")[[1]]))
  rownames(db) <- c()
  colnames(db) <- db[1,]
  db <- db[-1,]
  df <- data.frame(db,stringsAsFactors=F)
  df$E.value <- as.numeric(df$E.value)
  df$MaxTC <- as.numeric(df$MaxTC)
  if(ind %in% 1:3){
    df$code <- paste(type,"_",df$Uniprot.code,sep="")
    df$Uniprot.code <- NULL
  }else if(ind %in% 4:5){
    df$code <- paste(type,"_",df$Target.code,sep="")
    df$Target.name <- paste(df$Target.name,";",df$Uniprot.code,sep="")
    df$Target.code <- NULL
    df$Uniprot.code <- NULL
  }
  # thres.eval <- 
  # tmp <- df[df$E.value <= thres.eval,]
  tmp <- df
  ##
  if(ind==2){
    concs <- log10(as.numeric(sub(".+ \\[([01\\.]+) .+","\\1",tmp$Target)))
    tmp$code <- paste(sub("E$","",tmp$code),concs,sep="")
  }
  print(paste(ind,":",txt))
  print(table(sapply(unique(tmp$code),function(u)length(unique(tmp$T[tmp$code==u])))))
  print(table(sapply(unique(tmp$T),function(u)length(unique(tmp$code[tmp$T==u])))))
  ##
  dict <- rbind(dict,tmp[,c("Compound.ID","code","Target.name","E.value","MaxTC")])
}
## all 'code's are unique.
setwd("~/Dropbox (HMS)/projects/02 TEA/204_SEA_ChEMBL/dicts/")
save(dict,file="MS-ecfp-all.rda") ##93535 relationships are e.value <= 1

