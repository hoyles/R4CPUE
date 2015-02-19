# Data cleaning
dataclean <- function(dat,checktg=F,allHBF=F) {
  dat$op_yr <- as.numeric(dat$op_yr)
  dat$op_mon <- as.numeric(dat$op_mon)
  dat$op_day <- as.numeric(dat$op_day)
  dat <- dat[dat$op_day < 32,]
  dat <- dat[!is.na(dat$op_day),]
  dat$lat <- as.numeric(dat$lat)
  dat$latcode <- as.numeric(dat$latcode)
  dat$lon <- as.numeric(dat$lon)
  dat$loncode <- as.numeric(dat$loncode)
  dat <- dat[dat$loncode %in% c(1,2),]
  dat$hbf <- as.numeric(dat$hbf)
  dat$tonnage <- as.numeric(dat$tonnage)
  dat$hooks <- as.numeric(dat$hooks)
  dat$alb <- as.numeric(dat$alb)
  dat$bet <- as.numeric(dat$bet)
  dat$yft <- as.numeric(dat$yft)
  dat$swo <- as.numeric(dat$swo)
  if(sum(is.na(dat$lat))>0) dat[is.na(dat$lat),]$lat <- 0
  if(sum(is.na(dat$alb))>0) dat[is.na(dat$alb),]$alb <- 0
  if(sum(is.na(dat$bet))>0) dat[is.na(dat$bet),]$bet <- 0
  if(sum(is.na(dat$yft))>0) dat[is.na(dat$yft),]$yft <- 0
  if(sum(is.na(dat$swo))>0) dat[is.na(dat$swo),]$swo <- 0
  dat <- dat[!is.na(dat$hooks),]
  dat <- dat[dat$hooks<10000,] # clean up outliers
  dat <- dat[dat$hooks>200,] 
  dat <- dat[dat$yft<250,]  
  dat <- dat[dat$bet<250,]
  dat <- dat[dat$alb<250,]
  dat <- dat[dat$tonnage<30000 | is.na(dat$tonnage),]
#  dat[dat$fishingcat =="0",]
#  dat <- dat[dat$fishingcat !=".",]
  dat <- dat[dat$fishingcat !="0",]
#  dat <- dat[dat$hbf != "  .",]
  dat <- dat[is.na(dat$hbf)==F | dat$op_yr < 1976,]
  dat <- dat[dat$hbf < 26 | is.na(dat$hbf)==T,]
#  if(allHBF==F) {
#    dat[dat$hbf>22,]$hbf <- 22     # pool hbf > 22 into 22
#    dat <- dat[dat$hbf > 4,]   # remove swordfish targeting in R1 and R2
#    }
#  dat$ncrew <- as.numeric(dat$ncrew)
  if(checktg) dat <- dat[dat$target == 3 | is.na(dat$target),] # tuna target  (remove to avoid a change in 1994 - but recent trend is more important)
  return(dat)
  }

# Data preparation

dataprep <- function(dat,alldat=F) {
  dat$lat_raw <- dat$lat
  dat$lon_raw <- dat$lon
  dat$lat[dat$latcode==2] <- (dat$lat_raw[dat$latcode==2]+1) * -1
  dat$lon[dat$loncode==2] <- 360 - (dat$lon_raw[dat$loncode==2] + 1)
  dat <- dat[dat$lon >= 0,]
  
  dat$lat5 <- 5 * floor(dat$lat/5)
  dat$lon5 <- 5 * floor(dat$lon/5)
  
  dat$reg <- 0
  dat[dat$lat <  40 & dat$lat >=  20 & dat$lon >= 110 & dat$lon < 170,]$reg <- 1
  dat[dat$lat <  40 & dat$lat >=  20 & dat$lon >= 170 & dat$lon < 210,]$reg <- 2
  dat[dat$lat <  20 & dat$lat >= -10 & dat$lon >= 110 & dat$lon < 170,]$reg <- 3
  dat[dat$lat <  20 & dat$lat >= -10 & dat$lon >= 170 & dat$lon < 210,]$reg <- 4
  dat[dat$lat < -10 & dat$lat >= -40 & dat$lon >= 140 & dat$lon < 170,]$reg <- 5
  dat[dat$lat < -10 & dat$lat >= -40 & dat$lon >= 170 & dat$lon < 210,]$reg <- 6
  dat[dat$lat <  40 & dat$lat >=  20 & dat$lon >= 210,]$reg <- 7
  dat[dat$lat <  20 & dat$lat >= -40 & dat$lon >= 210,]$reg <- 8

  dat$subreg <- 0
  dat[dat$lat <  40 & dat$lat >=  20 & dat$lon >= 110 & dat$lon < 170,]$subreg <- 1
  dat[dat$lat <  40 & dat$lat >=  20 & dat$lon >= 170 & dat$lon < 210,]$subreg <- 2
  dat[dat$lat <  20 & dat$lat >=   0 & dat$lon >= 110 & dat$lon < 150,]$subreg <- 3.1
  dat[dat$lat <  20 & dat$lat >=   0 & dat$lon >= 150 & dat$lon < 170,]$subreg <- 3.2
  dat[dat$lat <   0 & dat$lat >= -10 & dat$lon >= 110 & dat$lon < 150,]$subreg <- 3.3
  dat[dat$lat <   0 & dat$lat >= -10 & dat$lon >= 150 & dat$lon < 170,]$subreg <- 3.4
  dat[dat$lat <  20 & dat$lat >= -10 & dat$lon >= 170 & dat$lon < 180,]$subreg <- 4.1
  dat[dat$lat <  20 & dat$lat >= -10 & dat$lon >= 180 & dat$lon < 210,]$subreg <- 4.2
  dat[dat$lat < -10 & dat$lat >= -40 & dat$lon >= 140 & dat$lon < 170,]$subreg <- 5
  dat[dat$lat < -10 & dat$lat >= -40 & dat$lon >= 170 & dat$lon < 210,]$subreg <- 6
  dat[dat$lat <  40 & dat$lat >=  20 & dat$lon >= 210,]$subreg <- 7
  dat[dat$lat <  20 & dat$lat >= -40 & dat$lon >= 210,]$subreg <- 8

  
  dat$callsign[dat$callsign == "      "] <- ".     "
  dat$vessid <- as.numeric(as.factor(paste(dat$callsign)))
  if (alldat==F) dat <- dat[dat$vessid != 1,]
  dat$vessid <- as.numeric(as.factor(dat$vessid))
  
  dat$yrqtr <- dat$op_yr + floor((dat$op_mon)/3)/4 + 0.125
  dat$latlong <- paste(dat$lat5,dat$lon5,sep=".")
  dat <- dat[dat$yrqtr < 2012,]
  #dat <- dat[dat$reg %in% 1:6,]

  dat$newfishingcat <- NA
  dat <- dat[dat$fishingcat<=3,]
  dat <- dat[dat$op_yr < 1967 | dat$op_yr > 1970 | dat$fishingcat < 3,]
  dat <- dat[dat$op_yr <= 1957 | dat$fishingcat != ".",]
  dat[dat$op_yr <=1957 & dat$reg %in% c(1),]$newfishingcat <- 1
  dat[dat$op_yr <=1957 & dat$reg %in% c(0,2:8),]$newfishingcat <- 2
  dat[dat$op_yr >1957 & dat$op_yr<=1993 & dat$fishingcat==1,]$newfishingcat <- 1
  dat[dat$op_yr >1993 & dat$fishingcat==3,]$newfishingcat <- 1
  dat[dat$op_yr >1957 & dat$op_yr<=1966 & dat$fishingcat %in% c(2,3),]$newfishingcat <- 2
  dat[dat$op_yr >1966 & dat$op_yr<=1970 & dat$fishingcat %in% c(2),]$newfishingcat <- 2
  dat[dat$op_yr >=1971 & dat$op_yr<=1993 & dat$fishingcat %in% c(2,3),]$newfishingcat <- 2
  dat[dat$op_yr >1993 & dat$fishingcat %in% c(1,2),]$newfishingcat <- 2
  
  dat <- dat[dat$yrqtr > 1945,]
  
  dat$trip_yr <- as.numeric(substr(as.character(dat$trip_st),1,4))
  dat <- dat[dat$trip_yr > 1945 | is.na(dat$trip_yr),]

  dat$tripid <- paste(dat$vessid,dat$trip_st,sep="_")
  dat$tripid[dat$vessid == 1] <- NA
  dat$tripid[dat$trip_st == "       0"] <- NA

  a <- table(dat$tripid)
  dat$sets_per_trip <- NA
  dat$sets_per_trip <- a[match(dat$tripid,names(a))]
  
  return(dat)
  }


get.coefs <- function(model,nyrs) {
  coefs <- exp(c(0,summary(model)$coefficients[2:nyrs,1]))
  coefs <- coefs/mean(coefs)
  return(coefs)
  }

get.bin.coefs <- function(model,nyrs,dat) {
  mn <- logit(mean(dat[,3]!=0))
  a <- c(0,summary(model)$coefficients[2:nyrs,1])
  a <- a - mean(a) + mn 
  coefs <- inv.logit(a)
#  coefs <- coefs/mean(coefs)
  return(coefs)
  }

#get.coefs <- function(model,yr) {
#  a <- grep("yrqtr",names(model$coefficients))
#  estyrs <- as.numeric(gsub("as.factor(yrqtr)","",names(model$coefficients)[a],fixed=T))
#  coefs <- exp(c(0,summary(model)$coefficients[2:nyrs,1]))
#  coefs <- coefs/mean(coefs)
#  }

plot.slope.ratio <- function(coefs1,coefs2,yr,titl) {
  # base goes first, boat goes second
  windows(height=14,width=12)
  par(mfrow=c(2,1),mar=c(4,4,3,1))
  plot(yr,coefs1/coefs2,xlim=c(1976,2010),ylim=c(0,2),xlab="Year",ylab="Ratio of coefficients")
  title(main=titl,cex.main=1.2)
  lin <-lm(coefs1/coefs2 ~ yr)
  logl <-lm(log(coefs1/coefs2) ~ yr)
  lines(yr,exp(logl$coefficients[1] + logl$coefficients[2]*yr),lty=3)
  tt <- paste(prettyNum(100*logl$coefficients[2],digits=2,format="f"),
  "% \261 ",prettyNum(100*summary(logl)$coefficients[2,2],digits=2,format="f"),", p = ",
  prettyNum(summary(logl)$coefficients[2,4],
  digits=2,format="f"),sep="")
  text(min(yr)+5, 1.9, tt,font=2,col="red",cex=1.1)
  par(mar=c(5,4,1,1))
  plot(yr,coefs1,type="l",ylab="Relative abundance estimate",xlab="Year",ylim=c(0,2.5))
  lines(yr,coefs2,col="red")
  }

plot.agg.slope.ratio <- function(coefs1,coefs2,aggyr,opyr,titl,lab1="coefs1",lab2="coefs2",fname=NULL) {
  # agg goes first, op goes second
  windows(height=14,width=12)
  par(mfrow=c(2,1),mar=c(4,4,3,1))
  myr <- aggyr[match(opyr,aggyr)]
  coefs1 <- coefs1[match(opyr,aggyr)]
  plot(myr,coefs1/coefs2,ylim=c(0,2),xlab="Year",ylab="Ratio of coefficients")
  title(main=titl,cex.main=1.2)
  legend("bottomright",legend=paste(lab1,"/",lab2),col=1,lty=3)
  lin <-lm(coefs1/coefs2 ~ myr)
  logl <-lm(log(coefs1/coefs2) ~ myr)
  lines(myr,exp(logl$coefficients[1] + logl$coefficients[2]*myr),lty=3)
  tt <- paste(prettyNum(100*logl$coefficients[2],digits=2,format="f"),
  "% \261 ",prettyNum(100*summary(logl)$coefficients[2,2],digits=2,format="f"),", p = ",
  prettyNum(summary(logl)$coefficients[2,4],
  digits=2,format="f"),sep="")
  text(min(myr)+10, 1.9, tt,font=2,col="red",cex=1.1)
  par(mar=c(5,4,1,1))
  plot(myr,coefs1,type="l",ylab="Relative abundance estimate",xlab="Year",ylim=c(0,3.5))
  lines(myr,coefs2,col="red")
  legend("topright",legend=c(lab1,lab2),col=1:2,lty=c(1,1))
  if(is.null(fname)==F) savePlot(paste(fname,lab1,lab2),type="png")
  }

plot.agg.slope.ratio_2000 <- function(coefs1,coefs2,aggyr,opyr,titl,lab1="coefs1",lab2="coefs2") {
  # agg goes first, op goes second
  windows(height=14,width=12)
  par(mfrow=c(2,1),mar=c(4,4,3,1))
  myr <- aggyr[match(opyr,aggyr)]
  coefs1 <- coefs1[match(opyr,aggyr)]
  plot(myr,coefs1/coefs2,xlim=c(1976,2010),ylim=c(0,2),xlab="Year",ylab="Ratio of coefficients")
  title(main=titl,cex.main=1.2)
  legend("bottomright",legend=paste(lab1,"/",lab2),col=1,lty=3)
  
  lin <-lm(coefs1/coefs2 ~ myr)
  logl <-lm(log(coefs1/coefs2) ~ myr)
  logl2000 <-lm(log(coefs1[myr<=2000]/coefs2[myr<=2000]) ~ myr[myr<=2000])
  lines(myr[myr<=2000],exp(logl2000$coefficients[1] + logl2000$coefficients[2]*myr[myr<=2000]),lty=3)
  tt <- paste(prettyNum(100*logl2000$coefficients[2],digits=2,format="f"),
  "% \261 ",prettyNum(100*summary(logl2000)$coefficients[2,2],digits=2,format="f"),", p = ",
  prettyNum(summary(logl2000)$coefficients[2,4],
  digits=2,format="f"),sep="")
  text(min(myr)+5, 1.9, tt,font=2,col="red",cex=1.1)
  par(mar=c(5,4,1,1))
  plot(myr,coefs1,type="l",ylab="Relative abundance estimate",xlab="Year",ylim=c(0,2.5))
  lines(myr,coefs2,col="red")
  legend("topleft",legend=c(lab1,lab2),col=1:2,lty=c(1,1))
  }

plot.res <- function(summ1,summ2,name1,name2) {
  nyrs <- length(grep("yrqtr",rownames(summ1$coefficients)))+1
  coefs1 <- get.summ.coefs(summ1,nyrs)
  coefs2 <- get.summ.coefs(summ2,nyrs)
  yrs <- rownames(summ1$coefficients)[grep("yrqtr",rownames(summ1$coefficients))]
  yrs <- c(1980.25,as.numeric(substring(yrs,17)))
  windows(height=13,width=11)
  par(mfrow=c(2,1),mar=c(5,4,1,2))
  plot.slope.ratio(coefs1,coefs2,yrs,"")
  plot(yrs,coefs1,type="l",ylab="Relative abundance estimate",xlab="Year",ylim=c(0,2.5))
  lines(yrs,coefs2,col="red")
  legend("topleft",legend=c(name1, name2),lty=c(1,1),col=c("black","red"))
}

get.summ.coefs <- function(model,nyrs) {
  coefs <- exp(c(0,model$coefficients[2:nyrs,1]))
  coefs <- coefs/mean(coefs)
  }

make_formula <- function(runsp,modtype,addboat,splitboat=F,addmain=F,addbranch=F,addother=F,addalb=F) {
  fmla <- "~ as.factor(yrqtr) + as.factor(latlong) + poly(hbf,6)"
  modhead <- switch(modtype,logn=paste("log(",runsp,"/hooks+0.01)"),deltabin=paste(runsp,"!= 0"),deltapos=paste("log((",runsp,")/hooks)"),qp=runsp,propn=runsp)
  fmla <- paste(modhead,fmla)
  if(modtype %in% c("deltabin","qp")) fmla <- paste(fmla,"+ ns(hooks,6)")
  if(addboat & !splitboat) fmla <- paste(fmla,"+ as.factor(vessid)")
  if(addboat & splitboat) fmla <- paste(fmla,"+ as.factor(splitvess)")
  if(addmain) fmla <- paste(fmla,"+ as.factor(mainline) + as.factor(mainline):ns(hbf,6)")
  if(addbranch) fmla <- paste(fmla,"+ as.factor(branchline)")
  if(addother) fmla <- paste(fmla,"+ as.factor(other)")
  if(addalb) fmla <- paste(fmla,"+ as.factor(alb_cat)")
  return(fmla)
  }
  
aggregate_data <- function(dat,sp) {
  tab1 <- aggregate(dat[,3], list(dat$yrqtr, dat$latlong, dat$hbf), sum)
  tab2 <- aggregate(dat$hooks, list(dat$yrqtr, dat$latlong, dat$hbf), sum)  
  taball <- cbind(tab1[,1:4],tab2[,4])
  names(taball) <-  c("yrqtr", "latlong","hbf",sp,"hooks")
  taball$yrqtr <- as.numeric(as.character(taball$yrqtr))
  taball$latlong <- as.character(taball$latlong)
  taball$hbf <- as.numeric(as.character(taball$hbf))
  taball[,4] <- as.numeric(as.character(taball[,4]))
  taball$hooks <- as.numeric(as.character(taball$hooks))
  return(taball)
  }
  
mk_wts <- function(dat,wttype,catch=NULL) {
  if(wttype=="equal") wts <- NULL
  if(wttype=="propn") wts <- catch
  if(wttype=="area") {
    a <- tapply(dat$latlong,list(dat$latlong,dat$yrqtr),length)
    i <- match(dat$latlong,rownames(a))
    j <- match(dat$yrqtr,colnames(a))
    n <- mapply("[", list(a), i, j)
    wts <- 1/n
    }
  if(wttype=="catch") {
    if(is.null(catch)) catch <- tapply(dat$bet,list(dat$latlong),sum)
    a <- tapply(dat$latlong,list(dat$latlong,dat$yrqtr),length)
    i <- match(dat$latlong,rownames(a))
    j <- match(dat$yrqtr,colnames(a))
    n <- mapply("[", list(a), i, j)
    cwts <- mapply("[", list(catch), i)/sum(catch)
    wts <- cwts/n
    }
  return(wts)  
  }
  
mk_wts_integer <- function(dat,wttype,catch=NULL) {
  if(wttype=="equal") wts <- NULL
  if(wttype=="area") {
    a <- tapply(dat$latlong,list(dat$latlong,dat$yrqtr),length)
    i <- match(dat$latlong,rownames(a))
    j <- match(dat$yrqtr,colnames(a))
    n <- mapply("[", list(a), i, j)
    wts <- 1/n
    wts <- floor(10 * max(n) * wts)
    }
  if(wttype=="catch") {
    if(is.null(catch)) catch <- tapply(dat$bet,list(dat$latlong),sum)
    a <- tapply(dat$latlong,list(dat$latlong,dat$yrqtr),length)
    i <- match(dat$latlong,rownames(a))
    j <- match(dat$yrqtr,colnames(a))
    n <- mapply("[", list(a), i, j)
    cwts <- mapply("[", list(catch), i)/sum(catch)
    wts <- cwts/n
    }
  return(wts)  
  }
  
make_strat <- function(dat) {
  a <- tapply(dat$latlong,list(dat$latlong,dat$yrqtr),length)
  i <- match(dat$latlong,rownames(a))
  j <- match(dat$yrqtr,colnames(a))
  n <- mapply("[", list(a), i, j)
  return(n)
  }

samp_data <- function(dat,n,nsamp) {
  p <- nsamp / n
  r <- runif(length(n))
  d2 <- dat[r<p,]
  return(d2)
  }
  
# p(bet !=0) ~ yrqtr + latlong + poly(hbf,6) + ns(hooks,6) + vessel id
# log(bet) ~ yrqtr + latlong + poly(hbf,6) + vessel id
  
# select the dataset
select_data <- function(indat,runreg,runsp,mt,minqtrs=2,maxqtrs=500,addmain=F,addother=F,addalb=F,fc="both",bait="no23",llstrat=5,doboth=F) {
  gdat <- indat[indat$reg==runreg,]
  if(runreg==9) gdat <- indat[indat$reg==3 & indat$lat >= -5 & indat$lat <10,]
  if(llstrat!=5) gdat$latlong <- paste(llstrat*floor(gdat$lat/llstrat),llstrat*floor(gdat$lon/llstrat),sep=".")
  if(mt=="deltapos") gdat <- gdat[gdat[,runsp] > 0,]
  if(bait=="no23") gdat <- gdat[(gdat$bait!=2 & gdat$bait !=3) | is.na(gdat$bait),]
  a <- table(gdat$vessid,gdat$yrqtr)
  a <- apply(a>0,1,sum)
  table(a)
  a <- a[a >= minqtrs & a <= maxqtrs]
  gdat <- gdat[gdat$vessid %in% names(a),]
  a <- table(gdat$yrqtr);a
  a <- a[a>=100]
  gdat <- gdat[gdat$yrqtr %in% names(a),]
  a <- table(gdat$vessid);a
  a <- a[a>=100]
  gdat <- gdat[gdat$vessid %in% names(a),]
  if(sum(is.na(gdat$hbf)) > 0) gdat[is.na(gdat$hbf),]$hbf <- 5
  gdat <- gdat[gdat$hbf >= 5,]
  fcnum <- switch(fc,"both"=0,"OS"=1,"DW"=2)
  if (fc!="both") gdat <- gdat[gdat$newfishingcat==fcnum,]
  if(addmain) {
    gdat <- gdat[gdat$target==3,c("vessid","hooks","yft", "bet", "alb", "hbf", "yrqtr","latlong","mainline","branchline")]
    gdat <- gdat[gdat$mainline %in% c(1,2),]
    gdat <- gdat[gdat$branchline %in% c(1,2),]
    } else {
    gdat <- gdat[,c("vessid","hooks","yft", "bet", "alb", "hbf", "yrqtr","latlong")]
    } 
  if(addother) {
      a <- (gdat[,switch(runsp,"bet"="yft","yft"="bet")]+1)/gdat$hooks
      divs <- c(0,0.1, 0.5, 0.9,1) 
      b <- quantile(a,divs)
      gdat$other <- findInterval(a,b)
      }
  if(addalb) {
      a <- (gdat[,"alb"]+1)/gdat$hooks
      divs <- c(0,0.1, 0.5, 0.9,1) 
      b <- quantile(a,divs)
      gdat$alb_cat <- findInterval(a,b)
      }
  a <- grep(switch(runsp,"bet"="yft","yft"="bet"),names(gdat))
  if(!doboth) gdat <- gdat[,-a]
  a <- grep("alb",names(gdat),fixed=T)[1]
  gdat <- gdat[,-a]
  return(gdat)
  }

qqDist <- function (x, standardise = F, add.median = F, ...) 
{
    n <- length(x)
    seq.length <- min(1000, n)
    if (standardise) {
        SEQ <- seq(1, 2 * n + 1, length = seq.length)/2
        U <- qnorm(qbeta(0.975, SEQ, rev(SEQ)))
        L <- qnorm(qbeta(0.025, SEQ, rev(SEQ)))
        if (add.median) 
            M <- qnorm(qbeta(0.5, SEQ, rev(SEQ)))
    }
    else {
        SD <- sqrt(var(x) * (n + 1)/n)
        SEQ <- seq(1, 2 * n + 1, length = seq.length)/2
        U <- mean(x) + SD * qt(qbeta(0.975, SEQ, rev(SEQ)), n - 
            1)
        L <- mean(x) + SD * qt(qbeta(0.025, SEQ, rev(SEQ)), n - 
            1)
        if (add.median) 
            M <- mean(x) + SD * qt(qbeta(0.5, SEQ, rev(SEQ)), 
                n - 1)
    }
    X <- qnorm((SEQ - 0.25)/(n + 0.5))
    qqnorm(x, main = "", ...)
    lines(X, U, type = "l",col=2)
    lines(X, L, type = "l",col=2)
    if (add.median) 
        lines(X, M, type = "l",col=2)
    invisible()
}

plotdiags <- function(res,ti="") {
  hist(res,nclass=200,freq=F,xlab="Residuals",main=ti)
  lines((-30:30)/10,dnorm((-30:30)/10,sd=sd(res)),col=2)
  sdres <- res/sd(res)
  qqDist(sdres,add.median=T)
  }

splitvessels <- function(indat,period) {
  vess <- unique(indat$vessid)
  indat$oldvess <- indat$vessid
  indat$vessid <- ""
  minyr <- maxyr <- vess
  for (i in 1:length(vess)) {
    a <- grep(vess[i],indat$oldvess)
    minyr <- min(indat[a,]$yrqtr)
    indat[a,]$vessid <- paste(indat[a,]$oldvess,floor((indat[a,]$yrqtr - minyr)/period))
    }
  return(indat)
  }
  
mt <- "deltabin" 
mt <- "deltapos"  

plot_effects <- function(model,indat,addmain=F,addbranch=F,addalb=F,addother=F,ti="") {
  cf <- model$coefficients
  pred <- predict(model,data=indat,type="terms",se.fit=T)
  fishlab <- switch(runsp,yft="Yellowfin",bet="Bigeye"); methlab <- switch(mt,deltabin="Delta-binomial",deltapos="Delta-positive",logl="Lognormal(+0.01)",propn="Proportion Bigeye")
  nfigs <- 6 + addmain + addbranch + addalb + addother
  mf <- switch(nfigs-5,c(2,3),c(3,3),c(3,3),c(3,3),c(3,4))
  hw <- c(14,19)
  windows(height=hw[1],width=hw[2])
  par(mfrow=mf,mar=c(5,4,2,1),oma=c(0,0,3,0)) 
  pr <- pred$fit ; prse <- pred$se.fit
  termlist <- dimnames(pred$fit)[[2]]
  llpos <- grep("latlong",termlist) 
  hbfpos  <- grep("hbf",termlist)
  mainpos <- grep("mainline",termlist)
  vesspos <- grep("vessid",termlist)
  branchpos <- grep("branchline",termlist)
  if(length(grep("hooks",termlist))>0) db <- T else db <- F
  
  index <- sort(unique(indat$yrqtr))
  b <- match(index,indat$yrqtr)
  out <- exp(pr[b,1])
  se1 <- exp(pr[b,1] - 1.96 * prse[b,1])
  se2 <- exp(pr[b,1] + 1.96 * prse[b,1])
  plot(as.numeric(as.character(index)),out,ylim=c(0,3),xlab="Year",ylab="Effect size",main="Year effects")
  segments(as.numeric(as.character(index)), se1,  as.numeric(as.character(index)), se2, lty=1, col="slate grey")
  points(as.numeric(as.character(index)), out, pch=16)

  index <- sort(unique(indat$latlong))
  b <- match(index,indat$latlong)
  out <- exp(pr[b,llpos])
  se1 <- exp(pr[b,llpos] - 1.96 * prse[b,llpos])
  se2 <- exp(pr[b,llpos] + 1.96 * prse[b,llpos])
  plot(as.numeric(as.character(index)),out,ylim=c(0,3),xlab="Latlong",ylab="Effect size",main="Spatial effects")
  segments(as.numeric(as.character(index)), se1,  as.numeric(as.character(index)), se2, lty=1, col="slate grey")
  points(as.numeric(as.character(index)), out, pch=16)
#  plot(indat$latlong,pr[,2],ylim=c(-.75,.75),xlab="Latitude",ylab="Effect size")

  ##plot coefficients
  library(maps)
  ll <- as.numeric(indat$latlong)
  index <- sort(unique(ll))
  lats <- trunc(ll)
  lons <- 1000*abs((ll - trunc(ll)))
  alat <- round(lats[match(index, ll)]*2)/2 + 2.5  
  alon <- round(lons[match(index, ll)]*2)/2 + 2.5
# lat <- trunc((dat$newlat[match(index, dat$latlong)] + 50)/5) * 5 + 2.5 - 50  
# long <- trunc((dat$newlong[match(index, dat$latlong)])/5) * 5 + 2.5
  coefs <- exp(pr[match(index, ll),llpos])
  coefs2 <- tapply(coefs, list(alon, alat), mean)
  image(sort(as.numeric(unique(alon))), sort(unique(alat)), coefs2, zlim=c(0.5, 2.5), ylab="Lat", xlab="Long")
  contour(sort(unique(alon)), sort(unique(alat)), coefs2, levels= c(0,1,2,2.5), add=TRUE,col=4)
  map(add=TRUE)

  index <- sort(unique(indat$vessid))
  b <- match(index,indat$vessid)
  out <- exp(pr[b,vesspos])
  se1 <- exp(pr[b,vesspos] - 1.96 * prse[b,vesspos])
  se2 <- exp(pr[b,vesspos] + 1.96 * prse[b,vesspos])
  plot(as.factor(index),out,type="n",ylim=c(0,3),xlab="Vessel ID",ylab="Effect size",main="")
  segments(match(index,index), se1,  match(index,index), se2, lty=1, col="slate grey")
  points(as.factor(index), out, pch=16)

  vts <- tapply(pr[,vesspos],list(indat$yrqtr,indat$vessid),mean)
  y<-exp(as.vector(vts))
  x<-rep(as.numeric(dimnames(vts)[[1]]),dim(vts)[[2]])
  plot(x,y,type="p",ylim=c(0,3),xlab="Vessel ID",ylab="Effect size",main="Vessel effects",pch=16,cex=0.9)
  mn <- tapply(exp(pr[,vesspos]),indat$yrqtr,mean)
  lines(as.numeric(dimnames(vts)[[1]]), mn, col=2,lwd=2)
#  plot(as.factor(indat$vessid),pr[,vesspos],ylim=c(-.75,.75),xlab="Vessel",ylab="Effect size")

  rsp <- pr[,hbfpos[1]]
  rsp.se <- prse[,hbfpos[1]]
  if(addmain) {
    a <- indat$mainline==1
    rsp2 <- pr[,hbfpos[1]] + pr[,hbfpos[2]] + pr[,mainpos[1]]
    rsp2.se <- exp(prse[,hbfpos[1]] + prse[,hbfpos[2]] + prse[,mainpos[1]])
    b <- match(sort(unique(indat[a,]$hbf)),indat[a,]$hbf)
    plot(indat[a,]$hbf[b],exp(rsp[a][b]),ylim=c(0,10),xlab="Nylon mainline HBF",ylab="Effect size",main=paste("Nylon HBF"))
    lines(indat[a,]$hbf[b],exp(rsp[a][b]+1.96*rsp.se[a][b]),lty=2)
    lines(indat[a,]$hbf[b],exp(rsp[a][b]-1.96*rsp.se[a][b]),lty=2)
    b <- match(sort(unique(indat[a==F,]$hbf)),indat[a==F,]$hbf)
    plot(indat[a==F,]$hbf[b],exp(rsp2[a==F][b]),ylim=c(0,3),xlab="Other mainline HBF",ylab="Effect size",main = "Other HBF")
    lines(indat[a==F,]$hbf[b],exp(rsp2[a==F][b] + 1.96*rsp2.se[a==F][b]),lty=2)
    lines(indat[a==F,]$hbf[b],exp(rsp2[a==F][b] - 1.96*rsp2.se[a==F][b]),lty=2)
    } else {
    b <- match(sort(unique(indat$hbf)),indat$hbf)
    plot(indat$hbf[b],exp(pr[,hbfpos][b]),ylim=c(0,3),xlab="HBF",ylab="Effect size",main="HBF effects")
    lines(indat$hbf[b],exp(rsp[b]+2*rsp.se[b]),lty=2)
    lines(indat$hbf[b],exp(rsp[b]-2*rsp.se[b]),lty=2)
    }

  if(addbranch) {
    index <- sort(unique(indat$branchline))
    b <- match(index,indat$branchline)
    out <- exp(pr[b,1])
    se1 <- exp(pr[b,1] - 1.96 * prse[b,1])
    se2 <- exp(pr[b,1] + 1.96 * prse[b,1])
    plot(as.numeric(as.character(index)),out,ylim=c(0,3),xlab="Branchline",ylab="Effect size",main="Branchline material effects")
    segments(as.numeric(as.character(index)), se1,  as.numeric(as.character(index)), se2, lty=1, col="slate grey")
    points(as.numeric(as.character(index)), out, pch=16)
    }

  if(addother) {
    index <- sort(unique(indat$other))
    b <- match(index,indat$other)
    out <- exp(pr[b,1])
    se1 <- exp(pr[b,1] - 1.96 * prse[b,1])
    se2 <- exp(pr[b,1] + 1.96 * prse[b,1])
    plot(as.numeric(as.character(index)),out,ylim=c(0,3),xlab="Other species category",ylab="Effect size",main="Other species effects")
    segments(as.numeric(as.character(index)), se1,  as.numeric(as.character(index)), se2, lty=1, col="slate grey")
    points(as.numeric(as.character(index)), out, pch=16)
    }

  if(addalb) {
    index <- sort(unique(indat$alb_cat))
    b <- match(index,indat$alb_cat)
    out <- exp(pr[b,1])
    se1 <- exp(pr[b,1] - 1.96 * prse[b,1])
    se2 <- exp(pr[b,1] + 1.96 * prse[b,1])
    plot(as.numeric(as.character(index)),out,ylim=c(0,3),xlab="Albacore catch category",ylab="Effect size",main="Albacore effects")
    segments(as.numeric(as.character(index)), se1,  as.numeric(as.character(index)), se2, lty=1, col="slate grey")
    points(as.numeric(as.character(index)), out, pch=16)
    }
  title(main=ti,cex.main=1.5,outer=T)
  }

plot.pacific <- function(plot_title="",lims=c(100,260,-45,45)) {
  plot(1,1, yaxt="n", xaxt="n", type="n", xlim=c(lims[1]+10,lims[2]-10), ylim=c(lims[3]+5,lims[4]-5), ylab="", xlab="", bg="lightblue")
  polygon(c(lims[1],lims[2],lims[2],lims[1]), c(lims[3],lims[3],lims[4],lims[4]), col="lightblue")
#  polygon(eez[,1], eez[,2], lwd=1, col="white")
#  lines(eez[,1], eez[,2], lwd=1, col="slate grey")
  map('world2Hires',  yaxt="n", xaxt="n", add=T, resolution=1)
  map('world2Hires',  region = c("USA","Hawaii","Mexico","Japan","China","South Korea","North Korea","Philippines","Vietnam","Laos","Taiwan","Fiji", "Vanuatu", "Malaysia", "Australia", "New Zealand", "Indonesia", "New Caledonia", "Papua New Guinea", "Solomon Islands"), fill=T, add=T, yaxt="n", xaxt="n", col="black", density=50)
  box(lwd=3)
  lines(c(210, 210, 230, 230), c(45, -2.5, -2.5, -45), lwd=2, lty=2)
  lines(c(170, 170), c(-35, 40), lwd=2, lty=1)
  lines(c(120, 210), c(20, 20), lwd=2, lty=1)
  lines(c(120, 230), c(-10, -10), lwd=2, lty=1)
  lines(c(120, 230), c(-35, -35), lwd=2, lty=1)
  axis(1, at=seq(lims[1],lims[2],by=10), labels=F)
  axis(2, at=seq(lims[3],lims[4],by=5), labels=F)
  latseq <- seq(lims[3]+10,lims[4]-10,by=10) ;latseq2 <- as.character(latseq) 
  lonseq <- seq(lims[1]+20,lims[2]-20,by=20) ;lonseq2 <- as.character(lonseq) 
  latseq2[latseq < 0] <- paste(abs(latseq[latseq < 0]),"S",sep="")
  latseq2[latseq > 0] <- paste(latseq[latseq > 0],"N",sep="")
  lonseq2[lonseq < 180] <- paste(lonseq2[lonseq < 180],"E",sep="")
  lonseq2[lonseq > 180] <- paste(360-lonseq[lonseq > 180],"W",sep="")
  axis(2, at=latseq, labels=latseq2, cex.axis=0.75)
  axis(1, at=lonseq, labels=lonseq2, cex.axis=0.75)
  mtext(side=3, line=0.5, plot_title)
}

storesumm <- function(mod,fname) {
  summry <- summary(mod)
  save(summry,file=paste(fname,".RData"))
  a <- capture.output(summry); 
  cat(a,file=paste(fname,".txt"),sep="\n",append=F)
  }
  
combine_delta_xl <- function(fnamedelta,fnamepos) {
    xl_delta <- read.csv(paste(fnamedelta,".csv",sep=""))
    xl_pos <- read.csv(paste(fnamepos,".csv",sep=""))
    pos <- match(xl_pos[,2],xl_delta[,2])
    coefs.boat <- xl_delta[pos,5] * xl_pos[,5]
    coefs.base <- xl_delta[pos,3] * xl_pos[,3]
    yrpos <- xl_delta[pos,2]
    fishlab <- switch(runsp,yft="Yellowfin",bet="Bigeye"); 
    plot.slope.ratio(coefs.base,coefs.boat,yrpos,titl=paste("Region",runreg,fishlab,"Delta lognormal combined"))
#    par(mar=c(5,4,1,1))
#    plot(yrpos,coefs.base,type="l",ylab="Relative abundance estimate",xlab="Year",ylim=c(0,2.5))
#    lines(yrpos,coefs.boat,col="red")
    fname2 <- gsub("deltabin","deltacomb",fnamedelta)
    savePlot(paste(fname2,".png",sep=""),type="png")
    write.csv(cbind(yrpos,coefs.base,coefs.boat),file=paste(fname2,".csv",sep=""))
    graphics.off()
    }
    
combine_delta_xl_indices <- function(fnamedelta,fnamepos) {
    xl_delta <- read.csv(paste(fnamedelta,".csv",sep=""))
    xl_pos <- read.csv(paste(fnamedelta,".csv",sep=""))
    coefs.boat <- xl_delta[,3] * xl_pos[,3]
    fishlab <- switch(runsp,yft="Yellowfin",bet="Bigeye"); 
    yrpos <- xl_delta[,2]
    windows();par(mar=c(5,4,1,1))
    plot(yrpos,coefs.boat,type="l",ylab="Relative abundance estimate",xlab="Year",ylim=c(0,2.5))
    fname2 <- gsub("deltabin","deltacomb",fnamedelta)
    savePlot(paste(fname2,".png",sep=""),type="png")
    write.csv(cbind(yrpos,coefs.boat),file=paste(fname2,".csv",sep=""))
    graphics.off()
    }
    
combine_delta <- function(fnamedelta,fnamepos) {
    load(paste("model.",fnamedelta,".base.RData",sep=""))
    load(paste("model.",fnamedelta,".boat.RData",sep=""))
    yrbin <- as.numeric(model.base$xlevels[[1]])
    a <- length(yrbin)
    coefsdeltabin.base <- get.coefs(model.base,a)
    coefsdeltabin.boat <- get.coefs(model.boat,a)
    rm(model.base,model.boat); gc()
    load(paste("model.",fnamepos,".base.RData",sep=""))
    load(paste("model.",fnamepos,".boat.RData",sep=""))
    yrpos <- as.numeric(model.base$xlevels[[1]])
    a <- length(yrpos)
    coefsdeltapos.base <- get.coefs(model.base,a)
    coefsdeltapos.boat <- get.coefs(model.boat,a)
    rm(model.base,model.boat); gc()
    a <- match(names(coefsdeltapos.base),names(coefsdeltabin.base))
    coefs.base <- coefsdeltabin.base[a] * coefsdeltapos.base
    coefs.boat <- coefsdeltabin.boat[a] * coefsdeltapos.boat
    fishlab <- switch(runsp,yft="Yellowfin",bet="Bigeye"); 
    plot.slope.ratio(coefs.base,coefs.boat,yrpos,titl=paste("Region",runreg,fishlab,"Delta lognormal combined"))
    par(mar=c(5,4,1,1))
    plot(yrpos,coefs.base,type="l",ylab="Relative abundance estimate",xlab="Year",ylim=c(0,2.5))
    lines(yrpos,coefs.boat,col="red")
    fname2 <- paste(fname," deltacomb",sep="")
    savePlot(paste(fname2,".png",sep=""),type="png")
    write.csv(cbind(yr,coefs.base,coefs.boat),file=paste(fname,".csv",sep=""))
    graphics.off()
    }
    


