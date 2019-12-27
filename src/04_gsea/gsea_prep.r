setwd(rda_dir)
dicts <- readRDS("dicts_MS-ecfp.rds")
eff <- readRDS("MS-eff.rds")
cmpds <- sub("X","",rownames(eff))

##
list.cmpds <- lapply(dicts,function(x){
	nl <- sapply(x,length)
	this.cmpds <- unique(unlist(x))
	# table(this.cmpds %in% cmpds)
	return(this.cmpds)
})

##
setwd(rda_dir)
save(cmpds,eff,dicts,list.cmpds,file="gsea_data.rda")
