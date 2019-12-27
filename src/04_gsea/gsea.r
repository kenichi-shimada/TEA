library(fgsea)
library(BiocParallel)
register(SerialParam())

setwd("/n/groups/mitchison/Kenichi/projects/TEA/gsea/")
load(file="gsea_data.rda") # cmpds,eff,dicts,list.cmpds

i <- as.numeric(commandArgs(TRUE))[1] ## 1-24

setting <- expand.grid(thres=1:12,mod=c("dmso","atoc"),stringsAsFactors=F)

thres <- setting$thres[i]
mod <- setting$mod[i]

dicts <- dicts[[thres]]
cmpds <- list.cmpds[[thres]]

eff.1 <- eff[paste0("X",cmpds),mod]
names(eff.1) <- cmpds
sorted.cmpds <- sort(eff.1,decreasing=T)

## no filter
set.seed(12345)
fgseaRes <- fgsea(pathways = dicts,
              stats = sorted.cmpds,
              minSize=10,
              maxSize=200,
              nperm=1e7)

setwd("/n/groups/mitchison/Kenichi/projects/TEA/gsea/fgsea")
fn <- paste0("fgsea_1e7_",i,".rds")
saveRDS(fgseaRes,file=fn)
