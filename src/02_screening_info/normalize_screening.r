library(open.xlsx)
library(VICTOR)

##########################################################################################
##    loading data
##########################################################################################
setwd(rda_dir);
load("MS-experimental-data.RData")

##check quality of the plates, by eye. ## first experiment showed some ATOC treatment didn't work. should remove them explicitly.

##########################################################################################
##    plate normalization
##########################################################################################
## normalize between the plate
## 1. within plate normalization (as usual)
## 2. normalize for AB incubation time (challenging)

## normalize for inubation time is really necessary.
## 1. get k's (first guess of how much cells are from different concs of CIL56 treatment) from bg.sub.raws.
## 2. based on these values, get t's (incubation time on each day)
## 3. get k's for each modulator treatment from bg.sub.raws?? or use ff and get the params from norms?
## 4. normalize (if not yet)
## 5. compute differences of area under curve between modulator and DMSO treamtnets (dAUC)

f <- function(k,t=5,vi.max=1)(vi.max*(1-exp(-k*t)))
ff <- function(k.cmpd,k.dmso,t,vi.max=1)f(k.cmpd,t,vi.max)/f(k.dmso,t,vi.max)

curve(f,0,1)

##########################################################################################
##    get initial k's
##########################################################################################

##bg subtracted raws
if(0){
  raws ##defined
  norms ##defined
}

##why such systematic differences?: many outer (top/bottom) wells have higher signals.
##plot section row-wise.
##########################################################################################
##    mats: data.frame for non-linear regression
##########################################################################################
mats <- data.frame(array(NA,c(0,3))) ##not normalized by 0%
for(mo in unique(layout$plate.384)){
  ind <- match(mo,unique(layout$plate.384))
  thisa <- anno[anno$mother==mo & anno$cell=="HT1080",]
  dmso.c <- unique(thisa$CIL[thisa$mod=="DMSO"])
  atoc.c <- unique(thisa$CIL[thisa$mod=="ATOC"])
  tl <- layout[layout$plate.384==mo & layout$cmpd=="DMSO",]
  rows <- match(tl$row.384,LETTERS[1:16])
  cols <- tl$column.384
  vi <- sapply(thisa$b,function(pl){
    sapply(seq(nrow(tl)),function(i)norms[[pl]][rows[i],cols[i]])
  })
  factor.k1 <- factor(rep(
    apply(thisa[match(colnames(vi),thisa$b),c(2,4)],1,paste,collapse="."),each=nrow(vi)))
  levels(factor.k1) <- paste0("k",c(5,1,2,3,6,4,7))

  factor.k2 <- factor(rep(thisa$modulator[match(colnames(vi),thisa$b)],each=nrow(vi)))
  factor.t <- factor(rep(paste("t",1:13,sep="")[as.numeric(sub(".+(..)$","\\1",mo))],each=length(vi)))
  df <- data.frame(vi=as.vector(vi),k1=factor.k1,k2=factor.k2,t=factor.t)
  mats <-rbind(mats,df)
}

##
if(0){
  vi1s <- sapply(levels(mats$t),function(fa)quantile(mats$vi[mats$t==fa&mats$k1=="k1"],1))## max is t12
  ##vi.max <- quantile(mats$vi[mats$t=="t12"&mats$k1=="k1"],seq(0,1,.1))
  vi.max <- 60000 ## from here, set as 60000
}

##########################################################################################
##   compute initial values
##########################################################################################
inits <- array(NA,c(7,13))
rownames(inits) <- paste("k",1:7,sep="")
colnames(inits) <- paste("t",1:13,sep="")
uparams <- as.matrix(unique(mats[,c(2,4)]))
for(i in seq(nrow(uparams))){
  v0 <- (mats[mats$k1==uparams[i,1]&mats$t==uparams[i,2],])
  thisv <- v0$vi
  options(show.error.messages=FALSE)
  init <- try(nls(thisv ~ f(k,1,vi.max=1),start=list(k=1))$m$getPars())
  if(class(init)!="numeric")init <- NA
  inits[uparams[i,1],uparams[i,2]] <- init
  options(show.error.messages=TRUE)
}

##plot to check the global pattern
par(mar=c(7,7,4,1)+.1)
image(1:13,1:7,t(inits[7:1,]),axes=F,xlab="",ylab="")
box()
axis(1,at=1:13,paste("t",1:13,sep=""),las=2)
axis(3,at=1:13,paste("t",1:13,sep=""),las=2)
axis(2,at=1:7,rev(paste(rep(c("DMSO","ATOC"),c(4,3)),
  c("0.00","0.63","1.25","2.50","0.00","2.50","5.00"))),las=1)
## initial values => init for t and k
inits.sub <- inits[-c(1,5),]
init.k <- apply(inits.sub,1,function(x)prod(x)^(1/length(x)))
##pos.init <- apply(inits,1,median,na.rm=T)[c(5,1)];names(pos.init) <- c("ATOC","DMSO")
##pos.init <- init.k[c(5,1)];names(pos.init) <- c("ATOC","DMSO")
##init.k2 <- apply(inits,1,function(x)prod(x,na.rm=T)^(1/length(x)))
##init.k3 <- apply(inits,1,median,na.rm=T)
init.t <- qr.solve(init.k,inits.sub) ##matrix operation

## use these initial values and new mats, to compute new values
inds <- as.character(mats$k1) %in% c("k2","k3","k4","k6","k7")
kin1 <- init.k[match(as.character(mats$k1),names(init.k))]
times <- paste(mats$t,".",mats$k2,sep="")
mats.1 <- data.frame(vi=mats$vi,kin1=kin1,time=times)[inds,]

## using final.t, fit nonlinear regression to all of the samples.
tmp.t <- nls(vi ~ f(kin1,t[time]),data=mats.1,
               start=list(t=rep(init.t,each=2)))$m$getPars()
reform.t <- matrix(tmp.t,nrow=2)[2:1,order(sub("^t([0-9])\\.","t0\\1.",levels(mats.1$t))[seq(2,26,2)])]
rownames(reform.t) <- c("DMSO","ATOC")
colnames(reform.t) <- paste("t",1:13,sep="")
final.t <- reform.t/apply(reform.t,1,median)
barplot(final.t,beside=T,main="estimated incubation time",col=c("lightblue","lightpink"))
abline(h=0:1,col="grey20")
legend("topright",c("DMSO","ATOC"),fill=c("lightblue","lightpink"))
##
tfs <- list()
barcodes <- names(norms)
##bc <- barcodes[13]
for(bc in barcodes){
  ti <- final.t[anno$modulator[anno$barcode==bc],
                paste("t",as.numeric(sub("MS3840","",anno$mother[anno$barcode==bc])),sep="")]
  ## make data frame
  tf <- array(sapply(norms[[bc]],function(x){
    if(x>=1)return(x)
    k <- -log(1-x)/ti
    new.x <- 1-exp(-k*1)
    return(new.x)
  }),c(16,24))
  tfs[[bc]] <- tf
}
tfs.old <-tfs
tfs <- norms

##
if(0){
  t1 <- seq(0,5,length.out=501)
  y1 <- f(1,t1)
  y2 <- f(1/2,t1)
  y3 <- f(1/3,t1)
  ##
  plot(range(t1),c(0,1),main="AlamarBlue incubation",lwd=2,type="n",
       xlab="AB incubation time (a.u.)",
       ylab="Normalized AlamarBlue")
  abline(h=0,v=0,col="grey50")
  abline(h=1,lty=2,col="grey50")
  lines(t1,y1,col=1)
  lines(t1,y2,col=2)
  lines(t1,y3,col=3)
  abline(v=c(.1,1,2,4),col="grey80")
  legend(4.2,0.4,c("100%","50%","33%"),col=1:3,lty=1)
  ##
  v <- array(NA,c(3,5))
  rownames(v) <- c("100%","50%","30%")
  colnames(v) <- c("true value","t=0.1","t=1","t=2","t=4")
  v[,1] <- c(1,1/2,1/3)
  v[,2] <- sapply(list(y1,y2,y3),function(x)x[t1==.1])
  v[,3] <- sapply(list(y1,y2,y3),function(x)x[t1==1])
  v[,4] <- sapply(list(y1,y2,y3),function(x)x[t1==2])
  v[,5] <- sapply(list(y1,y2,y3),function(x)x[t1==4])
  barplot(v,beside=T,col=c("lightblue","lightpink","lightgreen"))
  barplot(v/rep(v[1,],each=3),beside=T,col=c("lightblue","lightpink","lightgreen"))
}

## finally, compute dAUC
tf.norm <- function(x,ti){
  is.lt1 <- x < 1
  tmp <- x[is.lt1]
  k <- -log(1-tmp)/ti
  new.x <- 1-exp(-k*1)
  x[is.lt1] <- new.x
  return(x)
}
auc <- function(x,na.rm=T,normalize=T){
  n <- length(x)
  stopifnot(n>1)
  if(normalize){
    return(sum(x,x[-c(1,n)],na.rm=T)/n)
  }else{
    return(sum(x,x[-c(1,n)],na.rm=T))
  }
}
mod.auc <- function(mother,mod,normalize=T){
  tanno1 <- anno[anno$mother==mother & anno$modulator==mod & anno$cell=="HT1080",]
  ##
  cils <- as.character(unique(anno$CIL56[anno$modulator==mod]))
  mss <- as.character(unique(anno$MS[anno$modulator==mod]))
  ##
  mat <- array(NA,c(length(cils),length(mss)))
  rownames(mat) <- cils
  colnames(mat) <- mss
  dmso.well <- layout[layout$plate.384==mother&layout$cmpd=="DMSO",c("row.384","column.384")]
  dmso.i <- seq(nrow(dmso.well))
  dmso.r <- match(dmso.well[[1]],LETTERS[1:16])
  dmso.c <- dmso.well[[2]]
  for(ms in mss){
    for(cil in cils){
      thisbc <- tanno1$barcode[tanno1$CIL56==cil&tanno1$MS==ms]
      if(length(thisbc)==0){mat[cil,ms]<-NA}else{
        mat[cil,ms] <-  median(sapply(tfs[thisbc],function(x)
                                      median(sapply(dmso.i,function(i)x[dmso.r[i],dmso.c[i]]))))
      }
    }
  }
  if(mother=="MS384001"&mod=="DMSO"){
    if(normalize){
      n.mat <- diag(mat)/diag(mat)[1]
    }else{
      n.mat <- diag(mat)
    }
    dmso.auc <- auc(n.mat[-1])
  }else{
    if(normalize){
      n.mat <- mat/rep(mat[1,],each=length(cils))
    }else{
      n.mat <- mat
    }
    dmso.auc <- apply(n.mat[-1,],2,auc)
  }
  return(dmso.auc)
}
##
cmpd.auc <- function(ind,mod,normalize=T,std=F){
  cmpd <- layout$cmpd[ind]
  mother <- layout$plate.384[ind]
  stopifnot(length(layout$cmpd[layout$plate.384==mother&layout$cmpd==cmpd])==1)
  trow <- match(layout$row.384[ind],LETTERS[1:16]);
  tcol <- layout$column.384[ind]
  ##
  tanno1 <- anno[anno$mother==mother & anno$modulator==mod & anno$cell=="HT1080",]
  ##
  cils <- as.character(unique(anno$CIL56[anno$modulator==mod]))
  mss <- as.character(unique(anno$MS[anno$modulator==mod]))
  mat <- array(NA,c(length(cils),length(mss)))
  rownames(mat) <- cils
  colnames(mat) <- mss
  for(ms in mss){
    for(cil in cils){
      thisbc <- tanno1$barcode[tanno1$CIL56==cil&tanno1$MS==ms]
      if(length(thisbc)==0)mat[cil,ms] <- NA
      if(length(thisbc)==1)mat[cil,ms] <- tfs[[thisbc]][trow,tcol]
      if(length(thisbc)>1)mat[cil,ms] <-  median(sapply(tfs[thisbc],function(x)x[trow,tcol]))
    }
  }
  if(mother=="MS384001"&mod=="DMSO"){
    std.val <- diag(mat)[1]
    if(normalize){
      n.mat <- diag(mat)/std.val
    }else{
      n.mat <- diag(mat)
    }
    inh.auc <- auc(n.mat[-1])
  }else{
    std.val <- mat[1,]
    if(normalize){
      n.mat <- mat/rep(std.val,each=length(cils))
    }else{
      n.mat <- mat
    }
    inh.auc <- apply(n.mat[-1,],2,auc)
  }
  if(std){
    return(std.val)
  }else{
    return(inh.auc)
  }
}

##
dAUC <- array(NA,c(nrow(layout),2))
rownames(dAUC) <- layout$cmpd
colnames(dAUC) <- c("DMSO","ATOC")
stds <- aucs <- dAUC
##auc1 <- list();auc1$DMSO <- auc2$ATOC <- list()
##auc2 <- list();auc2$DMSO <- auc2$ATOC <- list()

##
thres <- 0.2
for(nr in seq(nrow(layout))){
  cmpd <- layout$cmpd[nr];if(cmpd=="DMSO")next;
  mo <- layout$plate.384[nr]
  trow <- match(layout$row.384[nr],LETTERS[1:16])
  tcol <- layout$column.384[nr]
  ## compute dAUC
  for(mod in c("DMSO","ATOC")){
    ##auc1[[mod]][[cmpd]] <- auc.noinh <- mod.auc(mo,mod,normalize=T)
    ##auc2[[mod]][[cmpd]] <- auc.inh <- cmpd.auc(ind,mod,normalize=T)
    normalize <- TRUE
    auc.noinh <- mod.auc(mo,mod,normalize=normalize)
    auc.inh <- cmpd.auc(nr,mod,normalize=normalize)
    ##auc.raw <- cmpd.auc(nr,mod,normalize=F)
    std <- cmpd.auc(nr,mod,std=T)
    d.auc <- auc.inh - auc.noinh
    ##d.auc.raw <- auc.raw - auc.noinh
    if(length(auc.noinh)==4 & length(auc.inh)==4& length(std)==4){
      max.ind <- which(abs(d.auc)==max(abs(d.auc)) & std > thres)
      if(length(max.ind)==0){
        std <- d.auc.raw <- d.auc <- NA
      }else if(length(max.ind)>0){
        d.auc <- d.auc[max.ind]
        ##d.auc.raw <- d.auc.raw[max.ind]
        std <- std[max.ind]
      }
    }else if(length(auc.noinh)!=length(auc.inh))stop();
    dAUC[nr,mod] <- d.auc
    ##aucs[nr,mod] <- d.auc.raw
    stds[nr,mod] <- std
    ##rm(auc.noinh,auc.inh)
  }
}
if(0){
  auc.atoc <- dAUC[,2]/2
  auc.dmso <- dAUC[,1]/4
  aucs.atoc <- aucs[,2]/2
  aucs.dmso <- aucs[,1]/4
}else{
  auc.atoc <- dAUC[,2]/sd(dAUC[,2],na.rm=T)
  auc.dmso <- dAUC[,1]/sd(dAUC[,1],na.rm=T)
  aucs.atoc <- aucs[,2]
  aucs.dmso <- aucs[,1]
}
## comparison of dose curves again

##bop
par(mfrow=c(1,1))##par(mfrow=c(3,5))
inds <- layout$plate.384 %in% unique(layout$plate.384)
r.auc <- range(auc.dmso,auc.atoc,na.rm=T)
par(mar=c(4,4,4,4)+.1)
col.mo <- rainbow(15)[as.numeric(factor(layout$plate.384))+1]
sum(col.mo)
plot(auc.dmso[inds],auc.atoc[inds],pch=20,col=col.mo,
     ##col=rainbow(13)[round(stds[,1]*10,0)[inds]],
     ##col=rainbow(13)[as.numeric(sub("MS3840","",layout$plate.384))[inds]],
     xlab="dAUC(DMSO)",ylab="dAUC(ATOC)",asp=1,
     xlim=r.auc,ylim=r.auc,main="all plates",cex=.3)
abline(h=0,v=0,a=0,b=1,col="grey80")

##
mm <- array(NA,c(16,24))
colnames(mm) <- as.character(seq(24))
rownames(mm) <- LETTERS[1:16]
##pdf("positional_effect_dAUC.pdf")
par(mfrow=c(3,5))
for(j in 1:13){
  mm1 <- mm2 <- mm
  col.mo <- as.numeric(factor(layout$plate.384)) %in% j
  for(ind in which(col.mo)){
    tr <- layout$row.384[ind]
    tc <- layout$column.384[ind]
    mm1[tr,tc] <- auc.dmso[ind]
    mm2[tr,tc] <- auc.atoc[ind]
  }
  thres <- 2
  mm1[is.na(mm1)] <- mm2[is.na(mm2)] <- 0
  mm1[mm1 > thres] <- thres; mm1[mm1 < -thres] <- -thres
  mm2[mm2 > thres] <- thres; mm2[mm2 < -thres] <- -thres
  par(mar=c(0,0,2,0)+.1)
  if(1){
    image(1:24,1:16,t(mm1[16:1,]),zlim=c(-1,1) * thres,
          col=rev(c(rgb(255:0,0,0,maxColorValue=255),rgb(0,1:255,0,maxColorValue=255))),
          asp=1,main=paste("DMSO",j),axes=F,xlab="",ylab="")
  }else{
    image(1:24,1:16,t(mm2[16:1,]),zlim=c(-1,1) * thres,
          col=rev(c(rgb(255:0,0,0,maxColorValue=255),rgb(0,1:255,0,maxColorValue=255))),
          asp=1,main=paste("ATOC",j),axes=F,xlab="",ylab="")
  }
}
##
if(0){
  par(xpd=T)
  legend(par()$usr[2],par()$usr[4],as.character((2:12)/10),pch=20,col=rainbow(13)[2:12],ncol=1)
  par(xpd=F)
}

##
for(p384 in unique(layout$plate.384)){
  inds <-layout$plate.384 %in% p384
  par(mar=c(4,4,4,4)+.1)
  plot(auc.dmso[inds],auc.atoc[inds],pch=20,
       col=rainbow(13)[round(stds[,1]*10,0)[inds]],
       ##col=rainbow(13)[as.numeric(sub("MS3840","",layout$plate.384))[inds]],
       xlab="dAUC(DMSO)",ylab="dAUC(ATOC)",asp=1,
       xlim=r.auc,ylim=r.auc,main=p384)
  abline(h=0,v=0,a=0,b=1,col="grey80")
  ##legend("topleft",as.character(1:13),pch=20,col=rainbow(13),ncol=3)
  par(xpd=T)
  legend(par()$usr[2],par()$usr[4],as.character((2:12)/10),pch=20,col=rainbow(13)[2:12],ncol=1)
  par(xpd=F)
  ##
  if(0){
    inds <-layout$plate.384 %in% "MS384006"
    plot(aucs.dmso[inds],auc.dmso[inds],#ylim=c(-3,3),xlim=c(-3,3),
         pch=20,col=rainbow(13)[as.numeric(sub("MS3840","",layout$plate.384))[inds]],
         xlab="raw intensity change",ylab="normalized intensity change")
    abline(h=0,v=0,a=0,b=1,col="grey80")
    legend("bottomright",as.character(1:13),pch=20,col=rainbow(13))
  }
  ##
}
##eop

## CIL56 only - effect of normalization
cil.only <- array(NA,c(13,5))
rownames(cil.only) <- unique(layout$plate.384)
cms <- unique(cbind(anno$C,anno$mod))
mothers <- unique(anno$mother)
##plot
par(mfrow=c(1,1))
par(mar=c(8,6,4,8)+.1)
plot(c(1,7),c(0,1.2),type="n",xlab="",ylab="",axes=F)
par(xpd=T)
legend(par()$usr[2],par()$usr[4],sort(mothers),lwd=2,col=1:8,lty=rep(1:2,c(8,5)))
par(xpd=F)
box()
abline(h=seq(0,1.2,.2),col=c(1,"grey80")[c(1,2,2,2,2,1,2)])
abline(v=1:7,col="grey80")
mtext("[CIL56](ug/mL),modulator",1,line=6)
mtext("Normalized/Transformed AlamarBlue",2,line=4)
axis(2,at=seq(0,1.2,.2))
axis(1,at=1:7,labels=apply(cms[,2:1],1,paste,collapse="_"),las=2)
tmp <- rep(NA,7)

for(mo in sort(mothers)){
  tl <- layout[layout$plate.384==mo&layout$cmpd=="DMSO",3:4]
  dmso.r <- match(tl$r, LETTERS[1:16])
  dmso.c <- tl$c
  i <- match(mo,sort(mothers))
  for(nr in seq(nrow(cms))){
    bcs <- anno$barcode[anno$mother==mo & anno$CIL==cms[nr,1] & anno$mod==cms[nr,2] & anno$cell=="HT1080"]
    vals <- sapply(tfs[bcs],function(mat)median(sapply(seq(nrow(tl)),function(i){
      row <- match(tl[i,1],LETTERS[1:16])
      col <- tl[i,2]
      mat[row,col]
    })))
    points(rep(nr,length(vals)),vals,col=rep(1:8,length.out=13)[i],pch=20)
    tmp[nr] <- median(vals)
  }
  lines(1:4,tmp[1:4],lwd=2,col=rep(1:8,length.out=13)[i],lty=rep(1:2,c(8,5))[i],pch=20)
  lines(5:7,tmp[5:7],lwd=2,col=rep(1:8,length.out=13)[i],lty=rep(1:2,c(8,5))[i],pch=20)
}
abline(v=4.5)
mo <- sort(mothers)[2]
## why no 'no change' points??
## remained to be answered for now.

