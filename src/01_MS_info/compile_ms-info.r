library(gdata)
library(dplyr)

setwd(data_dir);setwd("microsource/dbinfo")
db.info <- read.xls("MS library 96 plate map.xls",1)[1:4]

plate.names <- grep("MS384",dir(),value=T)

tmps <- lapply(1:13,function(i){
	tmp <- read.xls(plate.names[i],1,stringsAsFactors=F)
	tmp <- tmp[6:29,4:23]
	tmp <- data.frame(matrix(unlist(tmp),ncol=3,byrow=T)) %>%
		`colnames<-`(c("P1","add1","ID")) %>% 
		mutate(r1=sub("[0-9]{2}$","",add1)) %>%
		mutate(c1=sub("^[A-Z]","",add1)) %>%
		mutate(r2=rep(LETTERS[5:12],times=20)) %>%
		mutate(c2=rep(sprintf("%02s",3:22),each=8)) %>%
		mutate(P2=sprintf("MS384%03s",i))
	return(tmp)
})

plate.map <- do.call(rbind,tmps) %>% filter(!is.na(P1))
dim(plate.map) # 2000 x 8

setwd(rda_dir)
save(db.info,plate.map,file="MS_info.rda")

##
