# Run a test of aggregated vs individual. Aggregate the data, run the model, load the individual results and take the ratio
# Initial setup
#attach("C:/Users/simonh/Documents/I/Pop dy modeling/MFCL/R functions/R4MFCL.RData",pos=2)
#load("C:/Users/simonh/Documents/L/yft/2009/Data Preparation/CPUE/compare_methods/operat.RData")
require(R4MFCL)
#subset
a <- sample(1:nrow(llop),size=floor(0.4*nrow(llop)))
op <- llop[a,]
# rename fields
a <- match(c("bet_n","yft_n","latlong5","hook","boat_id"),names(op))
colnames(op)[a] <-  c("bet","yft","latlong","hooks","vessid")
# Test with SPC-held data
#- load SPC data
# map the data
library(maps)
library(mapproj)
library(mapdata)
plot.pacific.species(eez_file="C:/Users/simonh/Documents/I/Pop dy modeling/MFCL/R functions/EZNEW2.TXT")
points(op$lond,op$latd,cex=0.5)

#- aggregate SPC data

# Now with JPLL data
setwd("D:/SimonHOYLE2011/weight_methods/")

a <- tapply(op$bet!=0,list(op$yrqtr,op$reg,op$newfishingcat),mean)
write.csv(a,file="BET_pos_catches.csv")
a <- tapply(op$yft!=0,list(op$yrqtr,op$reg,op$newfishingcat),mean)
write.csv(a,file="YFT_pos_catches.csv")
a <- tapply(op$bet!=0,list(op$yrqtr,op$reg,op$newfishingcat),length)
write.csv(a,file="BET_ncatches.csv")
a <- tapply(op$yft!=0,list(op$yrqtr,op$reg,op$newfishingcat),length)
write.csv(a,file="YFT_ncatches.csv")

# test sampling and its effect on results
table(op$yrqtr,useNA="always")
mt <- "logn"
glmdat <- select_data(op,runreg=3,runsp="bet",mt,minqtrs,maxqtrs,addmain,addother=addother,fc="OS")
n <- make_strat(glmdat)
d10 <- samp_data(glmdat,n,10)
d20 <- samp_data(glmdat,n,20)
d30 <- samp_data(glmdat,n,30)
d50 <- samp_data(glmdat,n,50)
d10.wts <- mk_wts(d10,wttype="area")
d20.wts <- mk_wts(d20,wttype="area")
d30.wts <- mk_wts(d30,wttype="area")
d50.wts <- mk_wts(d50,wttype="area")
fmla.boat <- make_formula(runsp="bet",mt,addboat=T,splitboat=F,addmain=addmain,addbranch=addbranch)
model.opboat.d10  <- glm(fmla.boat,data=d10,weights=d10.wts);gc()
model.opboat.d20  <- glm(fmla.boat,data=d20,weights=d20.wts);gc()
model.opboat.d30  <- glm(fmla.boat,data=d30,weights=d30.wts);gc()
model.opboat.d50  <- glm(fmla.boat,data=d50,weights=d50.wts);gc()
dyr <- sort(unique(d50$yrqtr))
a <- length(unique(d10$yrqtr))
coefs.d10 <- get.coefs(model.opboat.d10,a)
a <- length(unique(d20$yrqtr))
coefs.d20 <- get.coefs(model.opboat.d20,a)
a <- length(unique(d30$yrqtr))
coefs.d30 <- get.coefs(model.opboat.d30,a)
a <- length(unique(d50$yrqtr))
coefs.d50 <- get.coefs(model.opboat.d50,a)
plot.agg.slope.ratio(coefs.d10,coefs.d50,dyr,dyr, titl=paste("Region",runreg,fishlab,methlab,"d10 vs d50"),lab1="d10",lab2="d50",fname="Compare subsampling")  # plot results
plot.agg.slope.ratio(coefs.d20,coefs.d50,dyr,dyr, titl=paste("Region",runreg,fishlab,methlab,"d20 vs d50"),lab1="d20",lab2="d50",fname="Compare subsampling")  # plot results
plot.agg.slope.ratio(coefs.d30,coefs.d50,dyr,dyr, titl=paste("Region",runreg,fishlab,methlab,"d30 vs d50"),lab1="d30",lab2="d50",fname="Compare subsampling")  # plot results
rm(model.opboat.d20,model.opboat.d30,model.opboat.d10,model.opboat.d50); gc()

mt <- "qp"
glmdat <- select_data(op,runreg=3,runsp="bet",mt,minqtrs,maxqtrs,addmain,addother=addother,fc="OS")
n <- make_strat(glmdat)
d10 <- samp_data(glmdat,n,10)
d20 <- samp_data(glmdat,n,20)
d30 <- samp_data(glmdat,n,30)
d50 <- samp_data(glmdat,n,50)
d10.wts <- mk_wts(d10,wttype="area")
d20.wts <- mk_wts(d20,wttype="area")
d30.wts <- mk_wts(d30,wttype="area")
d50.wts <- mk_wts(d50,wttype="area")
fmla.boat <- make_formula(runsp="bet",modtype="qp",addboat=T,splitboat=F,addmain=addmain,addbranch=addbranch)
model.opboat.d10.qp  <- glm(fmla.boat,data=d10,weights=d10.wts,family=quasipoisson);gc()
model.opboat.d20.qp  <- glm(fmla.boat,data=d20,weights=d20.wts,family=quasipoisson);gc()
model.opboat.d30.qp  <- glm(fmla.boat,data=d30,weights=d30.wts,family=quasipoisson);gc()
model.opboat.d50.qp  <- glm(fmla.boat,data=d50,weights=d50.wts,family=quasipoisson);gc()
a <- length(unique(d10$yrqtr))
coefs.d10.qp <- get.coefs(model.opboat.d10.qp,a)
a <- length(unique(d20$yrqtr))
coefs.d20.qp <- get.coefs(model.opboat.d20.qp,a)
a <- length(unique(d30$yrqtr))
coefs.d30.qp <- get.coefs(model.opboat.d30.qp,a)
a <- length(unique(d50$yrqtr))
coefs.d50.qp <- get.coefs(model.opboat.d50.qp,a)
plot.agg.slope.ratio(coefs.d10.qp,coefs.d50.qp,dyr,dyr, titl=paste("Region",runreg,fishlab,methlab,"d10 vs d50 qp"),lab1="d10",lab2="d50",fname="Compare subsampling qp")  # plot results
plot.agg.slope.ratio(coefs.d20.qp,coefs.d50.qp,dyr,dyr, titl=paste("Region",runreg,fishlab,methlab,"d20 vs d50 qp"),lab1="d20",lab2="d50",fname="Compare subsampling qp")  # plot results
plot.agg.slope.ratio(coefs.d30.qp,coefs.d50.qp,dyr,dyr, titl=paste("Region",runreg,fishlab,methlab,"d30 vs d50 qp"),lab1="d30",lab2="d50",fname="Compare subsampling qp")  # plot results
save(file="Compare weighting qp.RData","model.opboat.d10.qp","model.opboat.d20.qp","model.opboat.d30.qp","model.opboat.d50.qp")
rm("model.opboat.d10.qp","model.opboat.d20.qp","model.opboat.d30.qp","model.opboat.d50.qp"); gc() 

# rerun_with_delta_lognormal_no_swo
attach("D:/SimonHOYLE2010/Rfiles/R4mfcl 2010_04_15/R4MFCL.RData",pos=2)
setwd("D:/SimonHOYLE2011/")
source("./Rfiles/support_functions.r")
require(splines)
require(boot)
load("op.RData")
setwd("D:/SimonHOYLE2011/weight_methods_delta")
minqtrs=2; maxqtrs=200; runreg <- 2;addmain<-F; addbranch<-F;addother=F;addalb=F;fc="both";runsp="bet";fc="OS"
#maxqtrs=200; minqtrs_byreg = c(2,2,8,8,2,2,2,8,2)
maxqtrs=200; minqtrs_byreg = c(2,2,2,2,2,2,2,2,2)
for(runreg in c(2:5,9,1,6)) {
#for(runreg in c(9,1,6)) {
#for(runreg in c(6)) {
  for(fc in c("OS","DW"))
#  for(fc in c("DW"))
    {
    for(runsp in c("bet","yft"))
      {
      minqtrs <- minqtrs_byreg[runreg]
      glmdat.bin <- select_data(op,runreg,runsp,mt="deltabin",minqtrs,maxqtrs,fc=fc)
      if(nrow(glmdat.bin)>250000) {
        n <- make_strat(glmdat.bin)
        glmdat.bin <- samp_data(glmdat.bin,n,50)
        }
      glmdat.pos <- glmdat.bin[glmdat.bin[,3]>0,]
      aggdat <- aggregate_data(glmdat.bin,sp=runsp)
      wtt.agg.equal <- mk_wts(aggdat,wttype="equal")
      wtt.agg.area  <- mk_wts(aggdat,wttype="area")
      wtt.opbin.equal  <- mk_wts(glmdat.bin,wttype="equal")
      wtt.opbin.area   <- mk_wts(glmdat.bin,wttype="area")
      wtt.oppos.equal  <- mk_wts(glmdat.pos,wttype="equal")
      wtt.oppos.area   <- mk_wts(glmdat.pos,wttype="area")
      fname <- paste("Compare weightings ",runsp,"_R",runreg,"eq ",minqtrs,"-",maxqtrs,"_qtrs ",fc,sep="")
      fmla.agg <- make_formula(runsp,modtype="logn",addboat=F)
      fmla.opbin <- make_formula(runsp,modtype="deltabin",addboat=F)
      fmla.oppos <- make_formula(runsp,modtype="deltapos",addboat=F)
      fmla.oplogn <- make_formula(runsp,modtype="logn",addboat=F)
      fmla.boatbin <- make_formula(runsp,modtype="deltabin",addboat=T)
      fmla.boatpos <- make_formula(runsp,modtype="deltapos",addboat=T)
      fmla.boatlogn <- make_formula(runsp,modtype="logn",addboat=T)
      model.agg.equal <- glm(fmla.agg,data=aggdat,weights=wtt.agg.equal);gc()
      a <- length(unique(aggdat$yrqtr))
      coefs.agg.equal <- get.coefs(model.agg.equal,a)
      save(file=paste(fname," model.agg.equal.RData",sep=""),model.agg.equal);rm(model.agg.equal);gc()
      model.agg.area  <- glm(fmla.agg,data=aggdat,weights=wtt.agg.area);gc()
      coefs.agg.area <- get.coefs(model.agg.area,a)
      save(file=paste(fname," model.agg.area.RData",sep=""),model.agg.area);rm(model.agg.area);gc()
      model.opbin.equal  <- glm(fmla.opbin,data=glmdat.bin,family="binomial");gc()
      abin <- length(unique(glmdat.bin$yrqtr))
      coefs.opbin.equal <- get.bin.coefs(model.opbin.equal,abin,glmdat.bin)
      save(file=paste(fname," model.opbin.equal.RData",sep=""),model.opbin.equal);rm(model.opbin.equal);gc()
      model.oppos.equal  <- glm(fmla.oppos,data=glmdat.pos,weights=wtt.oppos.equal,family="gaussian");gc()
      apos <- length(unique(glmdat.pos$yrqtr))
      coefs.oppos.equal <- get.coefs(model.oppos.equal,apos)
      save(file=paste(fname," model.oppos.equal.RData",sep=""),model.oppos.equal);rm(model.oppos.equal);gc()
#      model.opbin.area   <- glm(fmla.opbin,data=glmdat.bin,weights=wtt.opbin.area,family="binomial");gc()
#      coefs.opbin.area <- get.coefs(model.opbin.area,abin)
#      save(file=paste(fname," model.opbin.area.RData",sep=""),model.opbin.area);rm(model.opbin.area);gc()
      model.oppos.area   <- glm(fmla.oppos,data=glmdat.pos,weights=wtt.oppos.area,family="gaussian");gc()
      coefs.oppos.area <- get.coefs(model.oppos.area,apos)
      save(file=paste(fname," model.oppos.area.RData",sep=""),model.oppos.area);rm(model.oppos.area);gc()
      model.oplogn.area   <- glm(fmla.oplogn,data=glmdat.bin,weights=wtt.opbin.area,family="gaussian");gc()
      coefs.oplogn.area <- get.coefs(model.oplogn.area,abin)
      save(file=paste(fname," model.oplogn.area.RData",sep=""),model.oplogn.area);rm(model.oplogn.area);gc()
      model.opboatbin.equal  <- glm(fmla.boatbin,data=glmdat.bin,weights=wtt.opbin.equal,family="binomial");gc()
      coefs.opboatbin.equal <- get.bin.coefs(model.opboatbin.equal,abin,glmdat.bin)
      plot_effects(model.opboatbin.equal,indat=glmdat.bin); savePlot(paste(fname,"bin_effects.png",sep=""),type="png")
      save(file=paste(fname," model.opboatbin.equal.RData",sep=""),model.opboatbin.equal);rm(model.opboatbin.equal);gc()
      model.opboatpos.equal  <- glm(fmla.boatpos,data=glmdat.pos,weights=wtt.oppos.equal,family="gaussian");gc()
      coefs.opboatpos.equal <- get.coefs(model.opboatpos.equal,apos)
      save(file=paste(fname," model.opboatpos.equal.RData",sep=""),model.opboatpos.equal);rm(model.opboatpos.equal);gc()
#      model.opboatbin.area   <- glm(fmla.boatbin,data=glmdat.bin,weights=wtt.opbin.area,family="binomial");gc()
#      coefs.opboatbin.area <- get.coefs(model.opboatbin.area,abin)
#      plot_effects(model.opboatbin.area,indat=glmdat.bin); savePlot(paste(fname,"bin_effects.png",sep=""),type="png")
#      save(file=paste(fname," model.opboatbin.area.RData",sep=""),model.opboatbin.area);rm(model.opboatbin.area);gc()
      model.opboatpos.area   <- glm(fmla.boatpos,data=glmdat.pos,weights=wtt.oppos.area,family="gaussian");gc()
      coefs.opboatpos.area <- get.coefs(model.opboatpos.area,apos)  
      windows(); par(mfrow=c(2,1),mar=c(4,4,3,1))                           # plot diagnostics
      plotdiags(model.opboatpos.area$residuals,ti="Model including vessel")
      savePlot(paste(fname,"_diags.png",sep=""),type="png")
      plot_effects(model.opboatpos.area,indat=glmdat.pos); savePlot(paste(fname,"pos_effects.png",sep=""),type="png")
      save(file=paste(fname," model.opboatpos.area.RData",sep=""),model.opboatpos.area);rm(model.opboatpos.area);gc()
      model.opboatlogn.area   <- glm(fmla.boatlogn,data=glmdat.bin,weights=wtt.opbin.area,family="gaussian");gc()
      coefs.opboatlogn.area <- get.coefs(model.opboatlogn.area,abin)  
      windows(); par(mfrow=c(2,1),mar=c(4,4,3,1))                           # plot diagnostics
      plotdiags(model.opboatlogn.area$residuals,ti="Model including vessel")
      savePlot(paste(fname,"_diags_logn.png",sep=""),type="png")
      plot_effects(model.opboatlogn.area,indat=glmdat.bin); savePlot(paste(fname,"logn_effects_.png",sep=""),type="png")
      save(file=paste(fname," model.opboatlogn.area.RData",sep=""),model.opboatlogn.area);rm(model.opboatlogn.area);gc()
      
      aggyr <- sort(unique(aggdat$yrqtr))
      opyr <- sort(unique(glmdat.pos$yrqtr))
      myr <- aggyr[match(opyr,aggyr)]
      coefs2.agg.equal <- coefs.agg.equal[match(opyr,aggyr)]
      coefs2.agg.area <- coefs.agg.area[match(opyr,aggyr)]
      fishlab <- switch(runsp,yft="Yellowfin",bet="Bigeye"); methlab <- switch(mt,deltabin="Delta-binomial",deltapos="Delta-positive",logl="Lognormal(+0.5)")
      a <- match(names(coefs.oppos.equal),names(coefs.opbin.equal))             # make op.equal indices
      coefs.op.equal <- coefs.opbin.equal[a] * coefs.oppos.equal
#      a <- match(names(coefs.oppos.area),names(coefs.opbin.area))             # make op.area indices
#      coefs.op.area <- coefs.opbin.area[a] * coefs.oppos.area
      a <- match(names(coefs.oppos.area),names(coefs.opbin.equal))             # make op.area indices
      coefs.op.area2 <- coefs.opbin.equal[a] * coefs.oppos.area
      a <- match(names(coefs.opboatpos.equal),names(coefs.opboatbin.equal))             # make opboat.equal indices
      coefs.opboat.equal <- coefs.opboatbin.equal[a] * coefs.opboatpos.equal
#      a <- match(names(coefs.opboatpos.area),names(coefs.opboatbin.area))             # make opboat.area indices
#      coefs.opboat.area <- coefs.opboatbin.area[a] * coefs.opboatpos.area
      a <- match(names(coefs.opboatpos.area),names(coefs.opboatbin.equal))             # make opboat.area indices
      coefs.opboat.area2 <- coefs.opboatbin.equal[a] * coefs.opboatpos.area
      save(file=paste(fname,"indices.RData",sep=""),list=ls(pattern="coefs."))
      
      plot.agg.slope.ratio(coefs.agg.equal,coefs.agg.area,aggyr,aggyr, titl=paste("Region",runreg,fishlab,fc,"agg equal vs agg area"),lab1="agg equal",lab2="agg area",fname)  # plot results
      plot.agg.slope.ratio(coefs.op.equal,  coefs.op.area2,  opyr, opyr,  titl=paste("Region",runreg,fishlab,fc,"op equal vs op area2 d_logn"),lab1="op equal",lab2="op area2",fname)  # plot results
      plot.agg.slope.ratio(coefs.agg.equal,coefs.op.equal,opyr,opyr,titl=paste("Region",runreg,fishlab,fc,"agg equal vs op equal d_logn"),lab1="agg equal",lab2="op equal",fname)  # plot results
      plot.agg.slope.ratio(coefs.agg.area, coefs.op.area2, aggyr,opyr,titl=paste("Region",runreg,fishlab,fc,"agg area vs op area2 d_logn"),  lab1="agg area", lab2="op area2",fname)  # plot results
      plot.agg.slope.ratio(coefs.agg.equal,coefs.opboat.equal,aggyr,opyr,titl=paste("Region",runreg,fishlab,fc,"agg equal vs opboat equal d_logn"),lab1="agg equal",lab2="opboat equal",fname)  # plot results
      plot.agg.slope.ratio(coefs.agg.area, coefs.opboat.area2, aggyr,opyr,titl=paste("Region",runreg,fishlab,fc,"agg area vs opboat area2 d_logn"),  lab1="agg area", lab2="opboat area2",fname)  # plot results
      plot.agg.slope.ratio(coefs.op.equal,coefs.opboat.equal,opyr,opyr,titl=paste("Region",runreg,fishlab,fc,"op equal vs opboat equal d_logn"),lab1="op equal",lab2="opboat equal",fname)  # plot results
      plot.agg.slope.ratio(coefs.op.area2, coefs.opboat.area2, opyr,opyr,titl=paste("Region",runreg,fishlab,fc,"op area vs opboat area2 d_logn"),  lab1="op area", lab2="opboat area2",fname)  # plot results
      plot.agg.slope.ratio(coefs.opboat.equal, coefs.opboat.area2, opyr,opyr,titl=paste("Region",runreg,fishlab,fc,"opboat equal vs opboat area2 d_logn"),  lab1="opboat equal", lab2="opboat area2",fname)  # plot results
      plot.agg.slope.ratio(coefs.agg.equal, coefs.opboat.area2, aggyr,opyr,titl=paste("Region",runreg,fishlab,fc,"agg equal vs opboat area2 d_logn"),  lab1="agg equal", lab2="opboat area2",fname)  # plot results
#      plot.agg.slope.ratio(coefs.op.area2, coefs.opboatlogn.area, opyr,opyr,titl=paste("Region",runreg,fishlab,fc,"op area2 vs opboatlogn area d_logn"),  lab1="op area2", lab2="opboatlogn area",fname)  # plot results
#      plot.agg.slope.ratio(coefs.opboat.area2, coefs.opboatlogn.area, opyr,opyr,titl=paste("Region",runreg,fishlab,fc,"opboat area2 vs opboatlogn area d_logn"),  lab1="opboat area2", lab2="opboatlogn area",fname)  # plot results
      rm(list=ls(pattern="coefs."))
      graphics.off()
     }
    }
  }
  
attach("D:/SimonHOYLE2010/Rfiles/R4mfcl 2010_04_15/R4MFCL.RData",pos=2)
setwd("D:/SimonHOYLE2011/")
source("./Rfiles/support_functions.r")
require(splines)
require(boot)
load("op.RData")
setwd("D:/SimonHOYLE2011/weight_methods_delta_allbait")
minqtrs=2; maxqtrs=200; runreg <- 2;addmain<-F; addbranch<-F;addother=F;addalb=F;fc="both";runsp="bet";fc="OS"
maxqtrs=200; minqtrs_byreg = c(2,2,2,2,2,2,2,2,2)
for(runreg in c(2:5,9,1,6)) {
  for(fc in c("OS","DW"))
    {
    if(runreg!=6 | fc == "DW") {
      for(runsp in c("bet","yft"))
        {
        minqtrs <- minqtrs_byreg[runreg]
        glmdat.bin <- select_data(op,runreg,runsp,mt="deltabin",minqtrs,maxqtrs,fc=fc,bait="allbait")
        if(nrow(glmdat.bin)>2000) {
          n <- make_strat(glmdat.bin)
          glmdat.bin <- samp_data(glmdat.bin,n,2)
          }
        glmdat.pos <- glmdat.bin[glmdat.bin[,3]>0,]
        aggdat <- aggregate_data(glmdat.bin,sp=runsp)
        wtt.agg.equal <- mk_wts(aggdat,wttype="equal")
        wtt.agg.area  <- mk_wts(aggdat,wttype="area")
        wtt.opbin.equal  <- mk_wts(glmdat.bin,wttype="equal")
        wtt.opbin.area   <- mk_wts(glmdat.bin,wttype="area")
        wtt.oppos.equal  <- mk_wts(glmdat.pos,wttype="equal")
        wtt.oppos.area   <- mk_wts(glmdat.pos,wttype="area")
        fname <- paste("Compare weightings ",runsp,"_R",runreg,"eq ",minqtrs,"-",maxqtrs,"_qtrs ",fc,"_allbait",sep="")
        fmla.agg <- make_formula(runsp,modtype="logn",addboat=F)
        fmla.opbin <- make_formula(runsp,modtype="deltabin",addboat=F)
        fmla.oppos <- make_formula(runsp,modtype="deltapos",addboat=F)
        fmla.oplogn <- make_formula(runsp,modtype="logn",addboat=F)
        fmla.boatbin <- make_formula(runsp,modtype="deltabin",addboat=T)
        fmla.boatpos <- make_formula(runsp,modtype="deltapos",addboat=T)
        fmla.boatlogn <- make_formula(runsp,modtype="logn",addboat=T)
        model.agg.equal <- glm(fmla.agg,data=aggdat,weights=wtt.agg.equal);gc()
        a <- length(unique(aggdat$yrqtr))
        coefs.agg.equal <- get.coefs(model.agg.equal,a)
        save(file=paste(fname," model.agg.equal.RData",sep=""),model.agg.equal);rm(model.agg.equal);gc()
        model.agg.area  <- glm(fmla.agg,data=aggdat,weights=wtt.agg.area);gc()
        coefs.agg.area <- get.coefs(model.agg.area,a)
        save(file=paste(fname," model.agg.area.RData",sep=""),model.agg.area);rm(model.agg.area);gc()
        model.opbin.equal  <- glm(fmla.opbin,data=glmdat.bin,family="binomial");gc()
        abin <- length(unique(glmdat.bin$yrqtr))
        coefs.opbin.equal <- get.bin.coefs(model.opbin.equal,abin,glmdat.bin)
        save(file=paste(fname," model.opbin.equal.RData",sep=""),model.opbin.equal);rm(model.opbin.equal);gc()
        model.oppos.equal  <- glm(fmla.oppos,data=glmdat.pos,weights=wtt.oppos.equal,family="gaussian");gc()
        apos <- length(unique(glmdat.pos$yrqtr))
        coefs.oppos.equal <- get.coefs(model.oppos.equal,apos)
        save(file=paste(fname," model.oppos.equal.RData",sep=""),model.oppos.equal);rm(model.oppos.equal);gc()
        model.oppos.area   <- glm(fmla.oppos,data=glmdat.pos,weights=wtt.oppos.area,family="gaussian");gc()
        coefs.oppos.area <- get.coefs(model.oppos.area,apos)
        save(file=paste(fname," model.oppos.area.RData",sep=""),model.oppos.area);rm(model.oppos.area);gc()
        model.oplogn.area   <- glm(fmla.oplogn,data=glmdat.bin,weights=wtt.opbin.area,family="gaussian");gc()
        coefs.oplogn.area <- get.coefs(model.oplogn.area,abin)
        save(file=paste(fname," model.oplogn.area.RData",sep=""),model.oplogn.area);rm(model.oplogn.area);gc()
        model.opboatbin.equal  <- glm(fmla.boatbin,data=glmdat.bin,weights=wtt.opbin.equal,family="binomial");gc()
        coefs.opboatbin.equal <- get.bin.coefs(model.opboatbin.equal,abin,glmdat.bin)
        plot_effects(model.opboatbin.equal,indat=glmdat.bin); savePlot(paste(fname,"bin_effects.png",sep=""),type="png")
        save(file=paste(fname," model.opboatbin.equal.RData",sep=""),model.opboatbin.equal);rm(model.opboatbin.equal);gc()
        model.opboatpos.equal  <- glm(fmla.boatpos,data=glmdat.pos,weights=wtt.oppos.equal,family="gaussian");gc()
        coefs.opboatpos.equal <- get.coefs(model.opboatpos.equal,apos)
        save(file=paste(fname," model.opboatpos.equal.RData",sep=""),model.opboatpos.equal);rm(model.opboatpos.equal);gc()
        model.opboatpos.area   <- glm(fmla.boatpos,data=glmdat.pos,weights=wtt.oppos.area,family="gaussian");gc()
        coefs.opboatpos.area <- get.coefs(model.opboatpos.area,apos)  
        windows(); par(mfrow=c(2,1),mar=c(4,4,3,1))                           # plot diagnostics
        plotdiags(model.opboatpos.area$residuals,ti="Model including vessel")
        savePlot(paste(fname,"_diags.png",sep=""),type="png")
        plot_effects(model.opboatpos.area,indat=glmdat.pos); savePlot(paste(fname,"pos_effects.png",sep=""),type="png")
        save(file=paste(fname," model.opboatpos.area.RData",sep=""),model.opboatpos.area);rm(model.opboatpos.area);gc()
        model.opboatlogn.area   <- glm(fmla.boatlogn,data=glmdat.bin,weights=wtt.opbin.area,family="gaussian");gc()
        coefs.opboatlogn.area <- get.coefs(model.opboatlogn.area,abin)  
        windows(); par(mfrow=c(2,1),mar=c(4,4,3,1))                           # plot diagnostics
        plotdiags(model.opboatlogn.area$residuals,ti="Model including vessel")
        savePlot(paste(fname,"_diags_logn.png",sep=""),type="png")
        plot_effects(model.opboatlogn.area,indat=glmdat.bin); savePlot(paste(fname,"logn_effects_.png",sep=""),type="png")
        save(file=paste(fname," model.opboatlogn.area.RData",sep=""),model.opboatlogn.area);rm(model.opboatlogn.area);gc()
        
        aggyr <- sort(unique(aggdat$yrqtr))
        opyr <- sort(unique(glmdat.pos$yrqtr))
        myr <- aggyr[match(opyr,aggyr)]
        coefs2.agg.equal <- coefs.agg.equal[match(opyr,aggyr)]
        coefs2.agg.area <- coefs.agg.area[match(opyr,aggyr)]
        fishlab <- switch(runsp,yft="Yellowfin",bet="Bigeye"); methlab <- switch(mt,deltabin="Delta-binomial",deltapos="Delta-positive",logl="Lognormal(+0.5)")
        a <- match(names(coefs.oppos.equal),names(coefs.opbin.equal))             # make op.equal indices
        coefs.op.equal <- coefs.opbin.equal[a] * coefs.oppos.equal
        a <- match(names(coefs.oppos.area),names(coefs.opbin.equal))             # make op.area indices
        coefs.op.area2 <- coefs.opbin.equal[a] * coefs.oppos.area
        a <- match(names(coefs.opboatpos.equal),names(coefs.opboatbin.equal))             # make opboat.equal indices
        coefs.opboat.equal <- coefs.opboatbin.equal[a] * coefs.opboatpos.equal
        a <- match(names(coefs.opboatpos.area),names(coefs.opboatbin.equal))             # make opboat.area indices
        coefs.opboat.area2 <- coefs.opboatbin.equal[a] * coefs.opboatpos.area
        save(file=paste(fname,"indices.RData",sep=""),list=ls(pattern="coefs."))
        
        plot.agg.slope.ratio(coefs.agg.equal,coefs.agg.area,aggyr,aggyr, titl=paste("Region",runreg,fishlab,fc,"agg equal vs agg area"),lab1="agg equal",lab2="agg area",fname)  # plot results
        plot.agg.slope.ratio(coefs.op.equal,  coefs.op.area2,  opyr, opyr,  titl=paste("Region",runreg,fishlab,fc,"op equal vs op area2 d_logn"),lab1="op equal",lab2="op area2",fname)  # plot results
        plot.agg.slope.ratio(coefs.agg.equal,coefs.op.equal,opyr,opyr,titl=paste("Region",runreg,fishlab,fc,"agg equal vs op equal d_logn"),lab1="agg equal",lab2="op equal",fname)  # plot results
        plot.agg.slope.ratio(coefs.agg.area, coefs.op.area2, aggyr,opyr,titl=paste("Region",runreg,fishlab,fc,"agg area vs op area2 d_logn"),  lab1="agg area", lab2="op area2",fname)  # plot results
        plot.agg.slope.ratio(coefs.agg.equal,coefs.opboat.equal,aggyr,opyr,titl=paste("Region",runreg,fishlab,fc,"agg equal vs opboat equal d_logn"),lab1="agg equal",lab2="opboat equal",fname)  # plot results
        plot.agg.slope.ratio(coefs.agg.area, coefs.opboat.area2, aggyr,opyr,titl=paste("Region",runreg,fishlab,fc,"agg area vs opboat area2 d_logn"),  lab1="agg area", lab2="opboat area2",fname)  # plot results
        plot.agg.slope.ratio(coefs.op.equal,coefs.opboat.equal,opyr,opyr,titl=paste("Region",runreg,fishlab,fc,"op equal vs opboat equal d_logn"),lab1="op equal",lab2="opboat equal",fname)  # plot results
        plot.agg.slope.ratio(coefs.op.area2, coefs.opboat.area2, opyr,opyr,titl=paste("Region",runreg,fishlab,fc,"op area vs opboat area2 d_logn"),  lab1="op area", lab2="opboat area2",fname)  # plot results
        plot.agg.slope.ratio(coefs.opboat.equal, coefs.opboat.area2, opyr,opyr,titl=paste("Region",runreg,fishlab,fc,"opboat equal vs opboat area2 d_logn"),  lab1="opboat equal", lab2="opboat area2",fname)  # plot results
        plot.agg.slope.ratio(coefs.agg.equal, coefs.opboat.area2, aggyr,opyr,titl=paste("Region",runreg,fishlab,fc,"agg equal vs opboat area2 d_logn"),  lab1="agg equal", lab2="opboat area2",fname)  # plot results
        rm(list=ls(pattern="coefs."))
        graphics.off()
        }
      }
    }
  }
  

# Rerun quasipoisson with SWO targeting removed (two methods)
attach("D:/SimonHOYLE2010/Rfiles/R4mfcl 2010_04_15/R4MFCL.RData",pos=2)
setwd("D:/SimonHOYLE2011/")
source("./Rfiles/support_functions.r")
require(splines)
require(boot)
load("op.RData")
setwd("D:/SimonHOYLE2011/weight_methods_qp_noswo/")
minqtrs=2; maxqtrs=200; runreg <- 2;addmain<-F; addbranch<-F;addother=F;addalb=F;fc="both"
maxqtrs=200; minqtrs_byreg = c(2,2,2,2,2,2,2,2,2)
runreg=1;runsp="bet";fc="OS";mt="logn"
for(runreg in c(1:5,9,6)) {
  for(fc in c("OS","DW"))
    {
    if(runreg!=6 | fc == "DW") {
      for(bait in c("no23","allbait")) {
        for(runsp in c("bet","yft"))
          {
          minqtrs <- minqtrs_byreg[runreg]
          glmdat <- select_data(op,runreg,runsp,mt,minqtrs,maxqtrs,fc=fc)
          if(nrow(glmdat)>250000) {
            n <- make_strat(glmdat)
            glmdat <- samp_data(glmdat,n,50)
            }
          aggdat <- aggregate_data(glmdat,sp=runsp)
          wtt.agg.equal <- mk_wts(aggdat,wttype="equal")
          wtt.agg.area  <- mk_wts(aggdat,wttype="area")
          wtt.op.equal  <- mk_wts(glmdat,wttype="equal")
          wtt.op.area   <- mk_wts(glmdat,wttype="area")
          fname <- paste("Compare weightings ",runsp,"_R",runreg,"eq ",minqtrs,"-",maxqtrs,"_qtrs ",fc,"_",bait,sep="")
          fmla.agg <- make_formula(runsp,modtype="logn",addboat=F)
          fmla.op <- make_formula(runsp,modtype="qp",addboat=F)
          fmla.boat <- make_formula(runsp,modtype="qp",addboat=T)
          model.agg.equal <- glm(fmla.agg,data=aggdat,weights=wtt.agg.equal);gc()
          a <- length(unique(aggdat$yrqtr))
          coefs.agg.equal <- get.coefs(model.agg.equal,a)
          save(file=paste(fname," model.agg.equal.RData",sep=""),model.agg.equal);rm(model.agg.equal);gc()
          model.agg.area  <- glm(fmla.agg,data=aggdat,weights=wtt.agg.area);gc()
          coefs.agg.area <- get.coefs(model.agg.area,a)
          save(file=paste(fname," model.agg.area.RData",sep=""),model.agg.area);rm(model.agg.area);gc()
          model.op.equal  <- glm(fmla.op,data=glmdat,weights=wtt.op.equal,family="quasipoisson");gc()
          a <- length(unique(glmdat$yrqtr))
          coefs.op.equal <- get.coefs(model.op.equal,a)
          save(file=paste(fname," model.op.equal.RData",sep=""),model.op.equal);rm(model.op.equal);gc()
          model.op.area   <- glm(fmla.op,data=glmdat,weights=wtt.op.area,family="quasipoisson");gc()
          coefs.op.area <- get.coefs(model.op.area,a)
          save(file=paste(fname," model.op.area.RData",sep=""),model.op.area);rm(model.op.area);gc()
          model.opboat.equal  <- glm(fmla.boat,data=glmdat,weights=wtt.op.equal,family="quasipoisson");gc()
          coefs.opboat.equal <- get.coefs(model.opboat.equal,a)
          save(file=paste(fname," model.opboat.equal.RData",sep=""),model.opboat.equal);rm(model.opboat.equal);gc()
          model.opboat.area   <- glm(fmla.boat,data=glmdat,weights=wtt.op.area,family="quasipoisson");gc()
          coefs.opboat.area <- get.coefs(model.opboat.area,a)
          save(file=paste(fname," model.opboat.area.RData",sep=""),model.opboat.area);rm(model.opboat.area);gc()
          
          aggyr <- sort(unique(aggdat$yrqtr))
          opyr <- sort(unique(glmdat$yrqtr))
          myr <- aggyr[match(opyr,aggyr)]
          coefs2.agg.equal <- coefs.agg.equal[match(opyr,aggyr)]
          coefs2.agg.area <- coefs.agg.area[match(opyr,aggyr)]
          fishlab <- switch(runsp,yft="Yellowfin",bet="Bigeye"); methlab <- switch(mt,deltabin="Delta-binomial",deltapos="Delta-positive",logl="Lognormal(+0.5)")
          save(file=paste(fname,"indices.RData",sep=""),list=ls(pattern="coefs."))
          
          plot.agg.slope.ratio(coefs2.agg.equal,coefs2.agg.area,aggyr,aggyr, titl=paste("Region",runreg,fishlab,fc,"agg equal vs agg area"),lab1="agg equal",lab2="agg area",fname)  # plot results
          plot.agg.slope.ratio(coefs.op.equal,  coefs.op.area,  opyr, opyr,  titl=paste("Region",runreg,fishlab,fc,"op equal vs op area qp"),lab1="op equal",lab2="op area",fname)  # plot results
          plot.agg.slope.ratio(coefs2.agg.equal,coefs.op.equal,aggyr,opyr,titl=paste("Region",runreg,fishlab,fc,"agg equal vs op equal qp"),lab1="agg equal",lab2="op equal",fname)  # plot results
          plot.agg.slope.ratio(coefs2.agg.area, coefs.op.area, aggyr,opyr,titl=paste("Region",runreg,fishlab,fc,"agg area vs op area qp"),  lab1="agg area", lab2="op area",fname)  # plot results
          plot.agg.slope.ratio(coefs2.agg.equal,coefs.opboat.equal,aggyr,opyr,titl=paste("Region",runreg,fishlab,fc,"agg equal vs opboat equal qp"),lab1="agg equal",lab2="opboat equal",fname)  # plot results
          plot.agg.slope.ratio(coefs2.agg.area, coefs.opboat.area, aggyr,opyr,titl=paste("Region",runreg,fishlab,fc,"agg area vs opboat area qp"),  lab1="agg area", lab2="opboat area",fname)  # plot results
          plot.agg.slope.ratio(coefs.op.equal,coefs.opboat.equal,opyr,opyr,titl=paste("Region",runreg,fishlab,fc,"op equal qp vs opboat equal qp"),lab1="op equal",lab2="opboat equal",fname)  # plot results
          plot.agg.slope.ratio(coefs.op.area, coefs.opboat.area, opyr,opyr,titl=paste("Region",runreg,fishlab,fc,"op area qp vs opboat area qp"),  lab1="op area", lab2="opboat area",fname)  # plot results
          plot.agg.slope.ratio(coefs.agg.equal, coefs.opboat.area, aggyr,opyr,titl=paste("Region",runreg,fishlab,fc,"agg equal vs opboat area qp"),  lab1="agg equal", lab2="opboat area",fname)  # plot results
          save(file=paste(fname,"indices.RData",sep=""),list=ls(pattern="coefs."))
          rm(list=ls(pattern="coefs."))
          graphics.off()
          }
        }
      }
    }
  }


mt <- "qp" 
fmla <- make_formula(runsp,mt,addboat=F,splitboat=F,addmain=addmain,addbranch=addbranch)
fmla.boat <- make_formula(runsp,mt,addboat=T,splitboat=F,addmain=addmain,addbranch=addbranch)
model.op.area.qp <- glm(fmla,data=glmdat,weights=wtt.op.area,family=quasipoisson);gc()
model.opboat.area.qp <- glm(fmla.boat,data=glmdat,weights=wtt.op.area,family=quasipoisson);gc()
a <- length(unique(glmdat$yrqtr))
coefs.op.area.qp <- get.coefs(model.op.area.qp,a)
coefs.opboat.area.qp <- get.coefs(model.opboat.area.qp,a)

plot.agg.slope.ratio(coefs.op.area,  coefs.op.area.qp,  opyr, opyr,  titl=paste("Region",runreg,fishlab,methlab,"op area vs op area qp"),lab1="op area",lab2="op area qp")  # plot results
plot.agg.slope.ratio_2000(coefs.op.area.qp,  coefs.opboat.area.qp,  opyr, opyr,  titl=paste("Region",runreg,fishlab,methlab,"op area qp vs opboat area qp"),lab1="op area qp",lab2="opboat area qp")  # plot results


savePlot(paste(opname,"_agg results.png",sep=""),type="png")
write.csv(cbind(myr,coefs2.agg,cv2.agg,coefs.base,cv.base,coefs.boat,cv.boat),file=paste(fname,"agg.csv",sep=""))  # save summaries
rm(model.agg);gc()



  
#- test analysis with equal weighting per data pt
#- tst with equal weighting per area
#- test with weigting proportional to total catch in last 20 years

savewd <- getwd()
setwd("D:/SimonHOYLE2011/testwt")
load("D:/SimonHOYLE/op.RData")

newpd <- op[op$vessid!=2,]
compdir <- "D:/SimonHOYLE2011/testwt/"

aggregate_data <- function(dat) {
  tab1 <- aggregate(dat$bet, list(dat$yrqtr, dat$latlong, dat$hbf), sum)
  tab2 <- aggregate(dat$hooks, list(dat$yrqtr, dat$latlong, dat$hbf), sum)  
  taball <- cbind(tab1[,1:4],tab2[,4])
  names(taball) <-  c("yrqtr", "latlong","hbf","bet","hooks")
  taball$yrqtr <- as.numeric(as.character(taball$yrqtr))
  taball$latlong <- as.character(taball$latlong)
  taball$hbf <- as.numeric(as.character(taball$hbf))
  taball$bet <- as.numeric(as.character(taball$bet))
  taball$hooks <- as.numeric(as.character(taball$hooks))
  return(taball)
  }
  
mt="deltapos";runsp="bet"; minqtrs=2; maxqtrs=200; runreg <- 3;addmain<-F; addbranch<-F;addother=F;addalb=F
maxqtrs=200; minqtrs_byreg = c(2,2,8,8,2,2,2,8)
for(runreg in c(2)) {
  mt <- "logn"
  minqtrs <- minqtrs_byreg[runreg]
  for(runsp in c("bet")) {
    fname <- paste(runsp,"_R",runreg,"eq ",minqtrs,"-",maxqtrs,"_qtrs ",mt,ifelse(addother,"_othersp",""),ifelse(addalb,"_albcat",""),sep="")
    glmdat <- select_data(newpd,runreg,runsp,mt,minqtrs,maxqtrs,addmain,addother=addother)
    aggdat <- aggregate_data(glmdat)
    fmla.agg <- make_formula(runsp,mt,addboat=F,splitboat=F,addmain=addmain,addbranch=addbranch,addother=addother,addalb=addalb)
    model.agg <- glm(fmla.agg,data=aggdat);gc()
    a <- length(unique(glmdat$yrqtr))
    coefs.agg <- get.coefs(model.agg,a)
    cv.agg <- summary(model.agg)$coefficients[1:a,2]
    aggyr <- sort(unique(aggdat$yrqtr))
    fmla <- make_formula(runsp,mt,addboat=F,splitboat=F,addmain=addmain,addbranch=addbranch,addother=addother,addalb=addalb)
    model <- glm(fmla.agg,data=glmdat);gc()
    a <- length(unique(glmdat$yrqtr))
    coefs.op <- get.coefs(model,a)
    cv.op <- summary(model)$coefficients[1:a,2]
    opyr <- sort(unique(glmdat$yrqtr))
    myr <- aggyr[match(opyr,aggyr)]
    coefs2.agg <- coefs.agg[match(opyr,aggyr)]
    cv2.agg <- cv.agg[match(opyr,aggyr)]
    fishlab <- switch(runsp,yft="Yellowfin",bet="Bigeye"); methlab <- switch(mt,deltabin="Delta-binomial",deltapos="Delta-positive",logl="Lognormal(+0.5)")
    plot.agg.slope.ratio(coefs2.agg,coefs.op,aggyr,opyr,titl=paste("Region",runreg,fishlab,methlab,"op vs agg"))  # plot results
    savePlot(paste(opname,"_agg results.png",sep=""),type="png")
    write.csv(cbind(myr,coefs2.agg,cv2.agg,coefs.base,cv.base,coefs.boat,cv.boat),file=paste(fname,"agg.csv",sep=""))  # save summaries
    storesumm(model.agg,paste(fname,"_summ_base",sep=""))
    rm(model.agg);gc()

    fmla <- make_formula(runsp,mt,addboat=F,splitboat=F,addmain=addmain,addbranch=addbranch,addother=addother,addalb=addalb)
    model <- glm(fmla.agg,data=glmdat);gc()
    a <- length(unique(glmdat$yrqtr))
    coefs.op <- get.coefs(model,a)
    cv.op <- summary(model)$coefficients[1:a,2]
    opyr <- sort(unique(glmdat$yrqtr))
    myr <- aggyr[match(opyr,aggyr)]
    coefs2.agg <- coefs.agg[match(opyr,aggyr)]
    cv2.agg <- cv.agg[match(opyr,aggyr)]
    fishlab <- switch(runsp,yft="Yellowfin",bet="Bigeye"); methlab <- switch(mt,deltabin="Delta-binomial",deltapos="Delta-positive",logl="Lognormal(+0.5)")
    plot.agg.slope.ratio(coefs2.agg,coefs.op,aggyr,opyr,titl=paste("Region",runreg,fishlab,methlab,"op vs agg"))  # plot results


    graphics.off()
    }
  }

setwd(savewd)


