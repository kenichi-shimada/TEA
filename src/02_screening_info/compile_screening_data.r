library(open.xlsx)
library(VICTOR)

##########################################################################################
##    loading data
##########################################################################################
##loading plate readouts
setwd(data_dir);setwd("readout")
plates <- read.multi.victor()

setwd(data_dir)
## rename plates
bcs <- names(plates$r)
names(plates$r) <- sub("^.+:(.+)$","\\1",bcs) ##rename: each element can be called by its barcode
table(table(names(plates$r))) ## everything is unique, 738 plates

## layout and annotation
layout <- read.csv("layout.csv",stringsAsFactors=F);layout <- layout[layout$a!="x",]
anno <- read.csv("annotation.csv",stringsAsFactors=F);anno <- anno[anno$a!="x"&anno$a!="b",]


##########################################################################################
##    bg subtracted raws
##########################################################################################
### bg.sub.raws
raws <- list()
mothers <- unique(anno$mother)
for(m in mothers){
  thisa <- anno[anno$mother==m,]
  tbcs <- thisa$b
  tmpr <- plates$r[tbcs]
  bgbc <- thisa$b[thisa$cell=="media"]

  bg <- tmpr[[bgbc]]
  tmpr <- tmpr[-match(bgbc,names(tmpr))]
  bg <- bg - median(bg[2:15,c(2,23)])
  thisn <- lapply(tmpr,function(raw){
    bg.col <- median(raw[2:15,c(2,23)])
    k <- 1 ##to be changed
    bg.sub <- raw - bg * k - bg.col
    if(m=="MS384010"){
      pos.col <- c(4:6,19:21)
      pos.row <- c(3,14)
    }else if(m=="MS384002"){
      pos.col <- c(4:7,18:21)
      pos.row <- c(4,13)
    }else{
      pos.col <- c(4:7,14:17)
      pos.row <- c(4,13)
    }
    ##pos <- median(bg.sub[pos.row,pos.col]) ##removed for the first time
    ##norm <- bg.sub/pos ##removed for the first time
    norm <- bg.sub ##used for the first time
  })
  raws <- c(raws,thisn)
}

## ori.norms
ori.norms <- list()
for(m in mothers){
  thisa <- anno[anno$mother==m,]
  tbcs <- thisa$b
  tmpr <- plates$r[tbcs]
  bgbc <- thisa$b[thisa$cell=="media"]

  bg <- tmpr[[bgbc]]
  tmpr <- tmpr[-match(bgbc,names(tmpr))]
  bg <- bg - median(bg[2:15,c(2,23)])
  thisn <- lapply(tmpr,function(raw){
    bg.col <- median(raw[2:15,c(2,23)])
    k <- 1 ##to be changed
    bg.sub <- raw - bg * k - bg.col
    if(m=="MS384010"){
      pos.col <- c(4:6,19:21)
      pos.row <- c(3,14)
    }else if(m=="MS384002"){
      pos.col <- c(4:7,18:21)
      pos.row <- c(4,13)
    }else{
      pos.col <- c(4:7,14:17)
      pos.row <- c(4,13)
    }
    pos <- median(bg.sub[pos.row,pos.col]) ##removed for the first time
    norm <- bg.sub/pos ##removed for the first time
  })
  ori.norms <- c(ori.norms,thisn)
}

##normalization (within plate): somehow this should be treated later.
bc100 <- anno$b[anno$cell=="HT1080" & anno$CIL==0]
pos.plate <- array(1,c(16,24))
pos.plate[3:14,3:22] <-
  round(array(apply(sapply(ori.norms[bc100],function(x)as.vector(x)),1,median),c(16,24)),2)[3:14,3:22]
filter <- 1/pos.plate
norms <- lapply(ori.norms,function(x)x * filter)

####
setwd(rda_dir)
save(plates,layout,anno,raws,mothers,norms,file="MS-experimental-data.RData")
