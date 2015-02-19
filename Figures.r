#####################################
#Prepare figures
# Logsheets per year by region
windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,2,2,1))
for (r in 1:6) {
  llv <- op[op$reg==r,]
  a <- table(llv$vessid,llv$op_yr)
  a <- apply(a>0,1,sum)
  print(table(a))

  a <- tapply(llv$vessid,llv$op_yr,length)
  plot(names(a),a,xlab="Year",ylab="Logsheet records",main=paste("Region",r))
  }
savePlot(filename=paste("Number of records",sep=""),type="png")

lu <- function(x) {
  length(unique(x))
  }

# Vessels per year by region
windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,2,2,1))
for (r in 1:6) {
  llv <- op[op$reg==r,]
  a <- tapply(llv$vessid,llv$op_yr,lu)
  plot(names(a),a,xlab="Year",ylab="",main=paste("Region",r))
  }
savePlot(filename=paste("Unique vessels by year",sep=""),type="png")

# Time distribution of vessels
vy <- list()
windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,2,2,1))
for (r in 1:6) {
  llv <- op[op$reg==r,]
  vess <- unique(llv$vessid)
  minyr <- maxyr <- vess
  for (i in 1:length(vess)) {
    minyr[i] <- min(llv[vess[i]==llv$vessid,]$op_yr)
    maxyr[i] <- max(llv[vess[i]==llv$vessid,]$op_yr)
    }
  vessyrs <- data.frame(vess=vess,minyr=as.numeric(minyr),maxyr=as.numeric(maxyr),stringsAsFactors=F)
  vessyrs <- vessyrs[order(-floor(vessyrs$minyr),-floor(vessyrs$maxyr)),]
  vy[[r]] <- vessyrs
  plot(1:length(vess),1:length(vess),xlim=c(1978,2010),type="n",xlab="Years",ylab="Vessel",main=paste("Region",r))
  for (i in 1:length(vess)) {
    lines(c(floor(vessyrs[i,]$minyr),floor(vessyrs[i,]$maxyr)),c(i,i))
    }
  }
savePlot(filename=paste("Time distribution of vessels 1",sep=""),type="png")

windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,2,2,1))
for (r in 1:6) {
  vessyrs <- vy[[r]]
  plot(1:400,1:400,xlim=c(1952,2010),type="n",xlab="Years",ylab="Vessel",main=paste("Region",r))
  for (i in 1:length(vess)) {
    lines(c(floor(vessyrs[i,]$minyr),floor(vessyrs[i,]$maxyr)),c(i,i))
    }
  }
savePlot(filename=paste("Time distribution of vessels 1b",sep=""),type="png")

windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,2,2,1))
for (r in 1:6) {
  vessyrs <- vy[[r]]
  vessyrs <- vessyrs[order(-floor(vessyrs$maxyr)),]
  llv <- op[op$reg==r,]
  vess <- unique(llv$vessid)
  plot(1:length(vess),1:length(vess),xlim=c(1978,2010),type="n",xlab="Years",ylab="Vessel",main=paste("Region",r))
  for (i in 1:length(vess)) {
    lines(c(floor(vessyrs[i,]$minyr),floor(vessyrs[i,]$maxyr)),c(i,i))
    }
  }
savePlot(filename=paste("Time distribution of vessels 2",sep=""),type="png")


windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,2,2,1))
for (r in 1:6) {
  vessyrs <- vy[[r]]
  llv <- op[op$reg==r,]
  vess <- unique(llv$vessid)
  minyr <- vess
  plot(1:length(vess),1:length(vess),xlim=c(1978,2010),type="n",xlab="Years",ylab="Vessel",main=paste("Region",r))
  for (i in 1:length(vess)) {
    a <- floor(llv[vessyrs[i,1]==llv$vessid,]$yr)
    pp <- unique(a)
    pp2 <- tapply(a,a,length)
#    points(pp,rep(i,length(pp)),cex=0.6,pch=3)
    symbols(pp,rep(i,length(pp)),sqrt(pp2)/40, add = T, inches =FALSE)
    }
  }
savePlot(filename=paste("Time distribution of vessels 4",sep=""),type="png")

# Effort by region
windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1))
for (r in 1:6) {
  llv <- op[op$reg==r,]
  a <- tapply(llv$hooks,llv$yrqtr,sum)
  plot(names(a),a,xlab="Year",ylab="Hooks",main=paste("Region",r))
  }
savePlot(filename="Effort by region by yrqtr",type="png")

# Effort by region by yrqtr by fishingcat
windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1))
for (r in 1:6) {
  llv <- op[op$reg==r,]
  a <- tapply(llv$hooks,llv$yrqtr,sum)
  plot(names(a),a,xlab="Year",ylab="Hooks",main=paste("Region",r),col=1,pch=16)
  llv <- op[op$reg==r & op$newfishingcat==1,]  # offshore
  a <- tapply(llv$hooks,llv$yrqtr,sum)
  points(names(a),a,col=2,pch=2)
  llv <- op[op$reg==r & op$newfishingcat==2,]  # DW
  a <- tapply(llv$hooks,llv$yrqtr,sum)
  points(names(a),a,col=4,pch=4)
  }
legend("topright",legend=c("Total","Offshore","Distant water"),col=c(1,2,4),pch=c(16,2,4))
savePlot(filename="Effort by region by yrqtr by fishingcat",type="png")

# Sets by region by yrqtr by fishingcat
windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1))
for (r in 1:6) {
  llv <- op[op$reg==r,]
  a <- tapply(llv$hooks,llv$yrqtr,length)
  plot(names(a),a,xlab="Year",ylab="Sets",main=paste("Region",r),col=1,pch=16)
  llv <- op[op$reg==r & op$newfishingcat==1,]
  a <- tapply(llv$hooks,llv$yrqtr,length)
  points(names(a),a,col=2,pch=2)
  llv <- op[op$reg==r & op$newfishingcat==2,]
  a <- tapply(llv$hooks,llv$yrqtr,length)
  points(names(a),a,col=4,pch=4)
  }
legend("topright",legend=c("Total","Offshore","Distant water"),col=c(1,2,4),pch=c(16,2,4))
savePlot(filename="Sets by region by yrqtr by fishingcat",type="png")

# total effort by region and yearqtr
e <- tapply(op$hooks,list(op$yrqtr,op$reg),sum)
write.table(file="effort by region and yrqtr.csv",e,sep=",")

# Tonnage - should be a stacked histogram, but no time
a <- unique(op$vessid)
a <- op[match(a,op$vessid),]
windows(height=14,width=12); par(mfrow=c(2,1),mar=c(3,4,2,1))
hist(a[a$newfishingcat==1,]$tonnage,breaks=seq(0,500,20),ylim=c(0,300),xlab="GRT (metric tonnes)",ylab="Number of vessels",main="Vessels",col=4,density=15,angle=135)
hist(a[a$newfishingcat==2,]$tonnage,breaks=seq(0,500,20),col=1,density=25,angle=45,add=T)
legend("topright",legend=c("Distant water","Offshore"),col=c(4,1),density=c(25,15),angle=c(45,135))

hist(op[op$newfishingcat==1,]$tonnage,breaks=seq(0,500,20),xlab="GRT (metric tonnes)",,ylab="Number of sets",main="Sets",col=4,density=15,angle=135)
hist(op[op$newfishingcat==2,]$tonnage,breaks=seq(0,500,20),xlab="GRT (metric tonnes)",main="",col=1,density=25,angle=45,add=T)
legend("topright",legend=c("Distant water","Offshore"),col=c(4,1),density=c(25,15),angle=c(45,135))
savePlot(file="plot tonnage by fishingcat",type="png")

# Effort coverage by region
# Load JPLLagg from SPC database
ll <- read.table("D:/SPC_NRIFSF_data2011/QUERY091.csv",sep=",",header=T)
ll$region <- rep(0, length(ll$hhooks))
ll$region <- ifelse(ll$latd > 20 & ll$latd < 40 & ll$lond > 110 & ll$lond < 170, 1, ll$region)
ll$region <- ifelse(ll$latd > 20 & ll$latd < 40 & ll$lond > 170 & ll$lond < 210, 2, ll$region)
#ll$region <- ifelse(ll$latd > -10 & ll$latd < 20 & ll$lond > 135 & ll$lond < 170, 3, ll$region)
ll$region <- ifelse(ll$latd > -10 & ll$latd < 20 & ll$lond > 110 & ll$lond < 170, 3, ll$region)
ll$region <- ifelse(ll$latd > -10 & ll$latd < 20 & ll$lond > 170 & ll$lond < 210, 4, ll$region)
ll$region <- ifelse(ll$latd > -35 & ll$latd < -10 & ll$lond > 140 & ll$lond < 170, 5, ll$region)
ll$region <- ifelse(ll$latd > -35 & ll$latd < -10 & ll$lond > 170 & ll$lond < 210, 6, ll$region)
##INDO/PH region
##ll$region <- ifelse(ll$LATD > -10 & ll$LATD < 20 & ll$LOND > 110 & ll$LOND < 135, 7, ll$region)
ll <- ll[ll$region > 0,]
ll$yrqtr <- ll$yy + ll$qtr/4 - 0.125
effort <- tapply(ll$hhooks, list(ll$yrqtr, ll$region), sum)*100
albc <- tapply(ll$alb_no, list(ll$yrqtr, ll$region), sum)
betc <- tapply(ll$bet_no, list(ll$yrqtr, ll$region), sum)
yftc <- tapply(ll$yft_no, list(ll$yrqtr, ll$region), sum)
swoc <- tapply(ll$swo_no, list(ll$yrqtr, ll$region), sum)

windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1))
for (r in 1:6) {
  llv <- op[op$reg==r,]
  e <- tapply(llv$hooks,llv$yrqtr,sum)
  agge <- effort[,r] 
  plot(names(e),e/agge[match(names(e),names(agge))],xlab="Year",ylab="Hooks",main=paste("Region",r),ylim=c(0,3))
  abline(h=1,lty=2)
  }
savePlot(filename="Effort coverage by region by yrqtr",type="png")
windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1))
for (r in 1:6) {
  llv <- op[op$reg==r,]
  alb<- tapply(llv$alb,llv$yrqtr,sum)
  agge <- albc[,r] 
  plot(names(alb),alb/agge[match(names(alb),names(agge))],xlab="Year",ylab="ALB catch",main=paste("Region",r),ylim=c(0,3))
  abline(h=1,lty=2)
  }
savePlot(filename="ALB coverage by region by yrqtr",type="png")
windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1))
for (r in 1:6) {
  llv <- op[op$reg==r,]
  bet <- tapply(llv$bet,llv$yrqtr,sum)
  agge <- betc[,r] 
  plot(names(bet),bet/agge[match(names(bet),names(agge))],xlab="Year",ylab="BET catch",main=paste("Region",r),ylim=c(0,3))
  abline(h=1,lty=2)
  }
savePlot(filename="BET coverage by region by yrqtr",type="png")
windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1))
for (r in 1:6) {
  llv <- op[op$reg==r,]
  yft <- tapply(llv$yft,llv$yrqtr,sum)
  agge <- yftc[,r] 
  plot(names(yft),yft/agge[match(names(yft),names(agge))],xlab="Year",ylab="YFT catch",main=paste("Region",r),ylim=c(0,3))
  abline(h=1,lty=2)
  }
savePlot(filename="YFT coverage by region by yrqtr",type="png")
windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1))
for (r in 1:6) {
  llv <- op[op$reg==r,]
  swo <- tapply(llv$swo,llv$yrqtr,sum)
  agge <- swoc[,r] 
  plot(names(swo),swo/agge[match(names(swo),names(agge))],xlab="Year",ylab="SWO catch",main=paste("Region",r),ylim=c(0,3))
  abline(h=1,lty=2)
  }
savePlot(filename="SWO coverage by region by yrqtr",type="png")

load("alldatraw.RData")
dat$lat_raw <- as.numeric(dat$lat)
dat$lon_raw <- dat$lon
dat$lat[dat$latcode==2] <- (dat$lat_raw[dat$latcode==2]+1) * -1
dat$lon[dat$loncode==2] <- 360 - (dat$lon_raw[dat$loncode==2] + 1)
dat$reg <- 0
dat[dat$lat <  40 & dat$lat >=  20 & dat$lon >= 110 & dat$lon < 170,]$reg <- 1
dat[dat$lat <  40 & dat$lat >=  20 & dat$lon >= 170 & dat$lon < 210,]$reg <- 2
dat[dat$lat <  20 & dat$lat >= -10 & dat$lon >= 110 & dat$lon < 170,]$reg <- 3
dat[dat$lat <  20 & dat$lat >= -10 & dat$lon >= 170 & dat$lon < 210,]$reg <- 4
dat[dat$lat < -10 & dat$lat >= -35 & dat$lon >= 140 & dat$lon < 170,]$reg <- 5
dat[dat$lat < -10 & dat$lat >= -35 & dat$lon >= 170 & dat$lon < 210,]$reg <- 6
dat$yrqtr <- dat$op_yr + floor((dat$op_mon)/3)/4 + 0.125
dat <- dat[dat$reg >0 & dat$reg <7,]
windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1))
for (r in 1:6) {
  cld <- op[op$reg==r,]
  rwd <- dat[dat$reg==r,]
  cle <- tapply(cld$hooks,cld$yrqtr,sum)
  rwe <- tapply(rwd$hooks,rwd$yrqtr,sum)
  plot(names(cle),cle/rwe[match(names(cle),names(rwe))],xlab="Year",ylab="Proportion of hooks",main=paste("Region",r),ylim=c(0,1.1))
  }
savePlot(filename="Cleaned effort proportion by region by yrqtr",type="png")
rm(dat)

# Target by region through time
for (fc in c("OS","DW")) {
  windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1),oma=c(0,0,3,0))
  b <- op[op$yrqtr >1994 & op$newfishingcat==switch(fc,"OS"=1,"DW"=2),]
  for (r in 1:6) {
    llv <- b[b$reg==r,]
    a <- tapply(llv$hooks,list(llv$yrqtr,factor(llv$target,levels=1:3)),sum)
    tota <- apply(a,1,sum,na.rm=T)
    tota[tota==0]<-1 
    a[is.na(a)==T]<-0 
    plot(names(tota),a[,1]/tota,type="b",col=1,xlim=c(1994,2010),xlab="Year",ylab="Proportion of effort",main=paste("Region",r),ylim=c(0,1))
    lines(names(tota),a[,2]/tota,type="b",col=2,pch=2,lty=1)
    lines(names(tota),a[,3]/tota,type="b",col=3,pch=3,lty=1)
    }
  legend("bottomright",legend=c("Swordfish","Sharks","Other (Tuna)"),col=c(1,2,3),pch=c(1,2,3),lty=c(1,1,1))
  title(switch(fc,"DW"="Distant water","OS"="Offshore"),outer=T)
  savePlot(filename=paste("Target by region by yrqtr",fc),type="png")
}

# Bait type by region through time
for (fc in c("OS","DW")) {
  windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1),oma=c(0,0,3,0))
  b <- op[op$newfishingcat==switch(fc,"OS"=1,"DW"=2),]
  for (r in 1:6) {
    llv <- b[b$reg==r,]
    a <- tapply(llv$hooks,list(llv$yrqtr,factor(llv$bait,levels=1:4)),sum)
    tota <- apply(a,1,sum,na.rm=T)
    tota[tota==0]<-1 
    a[is.na(a)==T]<-0 
    plot(names(tota),a[,1]/tota,type="p",col=1,xlim=c(1952,2010),xlab="Year",ylab="Proportion of effort",main=paste("Region",r),ylim=c(0,1))
    lines(names(tota),a[,2]/tota,type="p",col=2,pch=2,lty=1)
    lines(names(tota),a[,3]/tota,type="p",col=3,pch=3,lty=1)
    lines(names(tota),a[,4]/tota,type="p",col=4,pch=4,lty=1)
    }
  legend("bottomright",legend=c("Pacific saury", "Squid", "Live bait", "Other"),col=c(1,2,3,4),pch=c(1,2,3,4),lty=c(1,1,1,1))
  title(switch(fc,"DW"="Distant water","OS"="Offshore"),outer=T)
  savePlot(filename=paste("Bait type by region by yrqtr",fc),type="png")
}



# Catch by region
windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1))
for (r in 1:6) {
  llv <- op[op$reg==r,]
  yft <- tapply(llv$yft,llv$yrqtr,sum)
  bet <- tapply(llv$bet,llv$yrqtr,sum)
  alb <- tapply(llv$alb,llv$yrqtr,sum)
  swo <- tapply(llv$swo,llv$yrqtr,sum)
  maxy <- max(c(yft,bet,alb,swo))
  plot(names(yft),yft,ylim=c(0,maxy),xlab="Year",ylab="Catch",main=paste("Region",r))
  points(names(bet),bet,col=2,pch=2)
  points(names(alb),alb,col=3,pch=3)
  points(names(swo),swo,col=4,pch=4)
  }
legend("topright",legend=c("Yellowfin","Bigeye","Albacore","Swordfish"),col=c(1,2,3,4),pch=c(1,2,3,4))
savePlot(filename="Catch by region allsp by yrqtr",type="png")

# CPUE by region
windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1))
for (r in 1:6) {
  llv <- op[op$reg==r,]
  yft <- tapply(llv$yft,llv$yrqtr,sum)
  eff <- tapply(llv$hooks,llv$yrqtr,sum)
  bet <- tapply(llv$bet,llv$yrqtr,sum)
  alb <- tapply(llv$alb,llv$yrqtr,sum)
  swo <- tapply(llv$swo,llv$yrqtr,sum)
  yft <- 100*yft/eff
  bet <- 100*bet/eff
  alb <- 100*alb/eff
  swo <- 100*swo/eff
  maxy <- max(c(yft,bet,alb,swo))
  plot(names(yft),yft,ylim=c(0,maxy),xlab="Year",ylab="Catch per hundred hooks",main=paste("Region",r))
  points(names(bet),bet,col=2,pch=2)
  points(names(alb),alb,col=3,pch=3)
  points(names(swo),swo,col=4,pch=4)
  }
legend("topright",legend=c("Yellowfin","Bigeye","Albacore","Swordfish"),col=c(1,2,3,4),pch=c(1,2,3,4))
savePlot(filename="CPUE nominal allsp by region by yrqtr",type="png")

# CPUE by region by fishingcat
windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1),oma=c(0,0,3,0))
for (fc in 1:2) {
  for (r in 1:6) {
    llv <- op[op$reg==r & op$newfishingcat==fc,]
    yft <- tapply(llv$yft,llv$yrqtr,sum)
    eff <- tapply(llv$hooks,llv$yrqtr,sum)
    bet <- tapply(llv$bet,llv$yrqtr,sum)
    alb <- tapply(llv$alb,llv$yrqtr,sum)
    swo <- tapply(llv$swo,llv$yrqtr,sum)
    yft <- 100*yft/eff
    bet <- 100*bet/eff
    alb <- 100*alb/eff
    swo <- 100*swo/eff
    maxy <- max(c(yft,bet,alb,swo))
    plot(names(yft),yft,ylim=c(0,maxy),xlab="Year",ylab="Catch per hundred hooks",main=paste("Region",r))
    points(names(bet),bet,col=2,pch=2)
    points(names(alb),alb,col=3,pch=3)
    points(names(swo),swo,col=4,pch=4)
    }
  legend("topright",legend=c("Yellowfin","Bigeye","Albacore","Swordfish"),col=c(1,2,3,4),pch=c(1,2,3,4))
  title(switch(fc,"Offshore","Distant water"),outer=T)
  savePlot(filename=paste("CPUE nominal allsp by region by yrqtr",fc),type="png")
  }

# 5 degree squares fished
dimu <- function(x) { length(unique(x)) }
windows(height=14,width=12)
par (mfrow=c(3,2),mar=c(3,4,2,1))
for (r in 1:6) {
  llv <- op[op$reg == r,]
  yq <- sort(unique(llv$yrqtr))
  strats <- tapply(paste(llv$lat5,llv$lon5),llv$yrqtr,dimu)
  plot(yq, strats, type="p", xlim=c(1952, 2010),pch=1,col=1,ylim=c(0,max(strats)),cex=1,ylab="5 x 5 spatial strata with reported effort",main=paste("Region",r))
#  mtext(side=3, paste("Region", r),line=0.5)
  }
savePlot("Number of spatial strata",type="png")

# Proportion sets with zero catches
windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1))
for (r in 1:6) {
  llv <- op[op$reg==r,]
  ay <- tapply(llv$yft==0,llv$yrqtr,sum)
  a <- tapply(llv$hooks,llv$yrqtr,length)
  ab <- tapply(llv$bet==0,llv$yrqtr,sum)
  ay <- 100*ay/a
  ab <- 100*ab/a
  maxy <- max(c(ay,ab))
  plot(names(ay),ay,ylim=c(0,maxy),xlab="Year",ylab="Proportion of zero catches",main=paste("Region",r))
  points(names(ab),ab,col=2,pch=2)
  }
legend("topright",legend=c("Yellowfin","Bigeye"),col=c(1,2),pch=c(1,2))
savePlot(filename="Proportion zeroes by region by yrqtr",type="png")

# Proportion sets with zero catches allspp
windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1))
for (r in 1:6) {
  llv <- op[op$reg==r,]
  ay <- tapply(llv$yft==0,llv$yrqtr,sum)
  a <- tapply(llv$hooks,llv$yrqtr,length)
  ab <- tapply(llv$bet==0,llv$yrqtr,sum)
  aalb <- tapply(llv$alb==0,llv$yrqtr,sum)
  aswo <- tapply(llv$swo==0,llv$yrqtr,sum)
  ay <- 100*ay/a
  ab <- 100*ab/a
  aalb <- 100*aalb/a
  aswo <- 100*aswo/a
  maxy <- max(c(ay,ab,aswo,aalb))
  plot(names(ay),ay,ylim=c(0,maxy),xlab="Year",ylab="Proportion of zero catches",main=paste("Region",r))
  points(names(ab),ab,col=2,pch=2)
  points(names(aalb),aalb,col=3,pch=3)
  points(names(aswo),aswo,col=4,pch=4)
  }
legend("topright",legend=c("Yellowfin","Bigeye","Albacore","Swordfish"),col=c(1,2,3,4),pch=c(1,2,3,4))
savePlot(filename="Proportion zeroes by region by yrqtr allspp",type="png")

# Proportion sets with zero catches allspp by fishingcat
windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1),oma=c(0,0,3,0))
for (fc in 1:2) {
  for (r in 1:6) {
    llv <- op[op$reg==r & op$newfishingcat==fc,]
    ay <- tapply(llv$yft==0,llv$yrqtr,sum)
    a <- tapply(llv$hooks,llv$yrqtr,length)
    ab <- tapply(llv$bet==0,llv$yrqtr,sum)
    aalb <- tapply(llv$alb==0,llv$yrqtr,sum)
    aswo <- tapply(llv$swo==0,llv$yrqtr,sum)
    ay <- 100*ay/a
    ab <- 100*ab/a
    aalb <- 100*aalb/a
    aswo <- 100*aswo/a
    maxy <- max(c(ay,ab,aswo,aalb))
    plot(names(ay),ay,ylim=c(0,maxy),xlab="Year",ylab="Proportion of zero catches",main=paste("Region",r))
    points(names(ab),ab,col=2,pch=2)
    points(names(aalb),aalb,col=3,pch=3)
    points(names(aswo),aswo,col=4,pch=4)
    }
  legend("topright",legend=c("Yellowfin","Bigeye","Albacore","Swordfish"),col=c(1,2,3,4),pch=c(1,2,3,4))
  title(switch(fc,"Offshore","Distant water"),outer=T)
  savePlot(filename=paste("Proportion zeroes by region by yrqtr allspp",fc),type="png")
}

# Maps of effort through time
library(maps)
library(mapproj)
library(mapdata)
windows(width=14,height=12); par(mfrow=c(3,2),mar=c(3,4,2,1))
for(yr in seq(1950,2009,by=10)) {
  plot.pacific(plot_title=yr)
  a <- op[op$yrqtr > yr & op$yrqtr < yr+10,]
  lats <- sort(unique(a$lat5))
  lons <- sort(unique(a$lon5))
  a <- table(a$lat5,a$lon5)
  for (i in 1:15) {
    symbols(lons+2.5,rep(lats[i],length(lons))+2.5,circles=sqrt(a[i,])/40,col=1,add = T, inches =F)
    }
  }
savePlot("Plot map sets", type="png")

for (fc in 1:2) {
  windows(width=14,height=12); par(mfrow=c(3,2),mar=c(3,4,2,1),oma=c(0,0,3,0))
  a2 <- op[op$newfishingcat==fc,]
  for(yr in seq(1950,2009,by=10)) {
    plot.pacific(plot_title=yr)
    a <- a2[a2$yrqtr > yr & a2$yrqtr < yr+10,]
    lats <- sort(unique(a$lat5))
    lons <- sort(unique(a$lon5))
    a <- table(a$lat5,a$lon5)
    for (i in 1:length(lats)) {
      symbols(lons+2.5,rep(lats[i],length(lons))+2.5,circles=sqrt(a[i,])/30,col=1,add = T, inches =F)
      }
    }
  title(switch(fc,"Offshore","Distant water"),outer=T)
  savePlot(paste("Plot map sets",fc), type="png")
  }
 #   symbols(pp,rep(i,length(pp)),sqrt(pp2)/40, add = T, inches =FALSE)

# Maps of mean HPB through time
windows(width=14,height=12); par(mfrow=c(3,2),mar=c(3,4,2,1))
for(yr in seq(1950,2002,by=10)) {
  plot.pacific(plot_title=yr)
  a <- op[op$yrqtr > yr & op$yrqtr < yr+10,]
  lats <- sort(unique(a$lat5))
  lons <- sort(unique(a$lon5))
  a <- tapply(a$hbf,list(a$lat5,a$lon5),mean,na.rm=T)
  for (i in 1:15) {
    text(lons+2.5,rep(lats[i],length(lons))+2.5,floor(a[i,]),col=1,add = T, inches =F)
    }
  }
savePlot("Plot map mean HBF", type="png")

# Maps of mean HPB through time by FC
for (fc in 1:2) {
  windows(width=14,height=12); par(mfrow=c(3,2),mar=c(3,4,2,1),oma=c(0,0,3,0))
  a2 <- op[op$newfishingcat==fc,]
  for(yr in seq(1950,2002,by=10)) {
    plot.pacific(plot_title=yr)
    a <- a2[a2$yrqtr > yr & a2$yrqtr < yr+10,]
    lats <- sort(unique(a$lat5))
    lons <- sort(unique(a$lon5))
    a <- tapply(a$hbf,list(a$lat5,a$lon5),mean,na.rm=T)
    for (i in 1:length(lats)) {
      text(lons+2.5,rep(lats[i],length(lons))+2.5,floor(a[i,]),col=1,add = T, inches =F)
      }
    }
  title(switch(fc,"Offshore","Distant water"),outer=T)
  savePlot(paste("Plot map mean HBF",fc), type="png")
  }

# Maps of median HPB through time 
windows(width=14,height=12); par(mfrow=c(3,2),mar=c(3,4,2,1))
for(yr in seq(1950,2002,by=10)) {
  plot.pacific(plot_title=yr)
  a <- op[op$yrqtr > yr & op$yrqtr < yr+10,]
  lats <- sort(unique(a$lat5))
  lons <- sort(unique(a$lon5))
  a <- tapply(a$hbf,list(a$lat5,a$lon5),median,na.rm=T)
  for (i in 1:length(lats)) {
    text(lons+2.5,rep(lats[i],length(lons))+2.5,(a[i,]),col=1)
    }
  }
savePlot("Plot map median HBF", type="png")

for (fc in 1:2) {
  windows(width=14,height=12); par(mfrow=c(3,2),mar=c(3,4,2,1),oma=c(0,0,3,0))
  a2 <- op[op$newfishingcat==fc,]
  for(yr in seq(1950,2002,by=10)) {
    plot.pacific(plot_title=yr)
    a <- a2[a2$yrqtr > yr & a2$yrqtr < yr+10,]
    lats <- sort(unique(a$lat5))
    lons <- sort(unique(a$lon5))
    a <- tapply(a$hbf,list(a$lat5,a$lon5),median,na.rm=T)
    for (i in 1:length(lats)) {
      text(lons+2.5,rep(lats[i],length(lons))+2.5,(a[i,]),col=1)
      }
    }
  title(switch(fc,"Offshore","Distant water"),outer=T)
  savePlot(paste("Plot map median HBF",fc), type="png")
  }

# Plot hbf by region by year
for (fc in c(1,2)){
  windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1),oma=c(1,1,3,1)); r<-3; hh <- 18
  b<- op[op$newfishingcat==fc,]
  for (r in 1:6) {
    a <- b[b$reg==r,]
    yrs <- sort(unique(a$op_yr))
    a <- table(a$hbf,a$op_yr)
    plot(1,1, type="n", xlab="Year", ylab="HBF", ylim = c(1,25), xlim=c(1950, 2010),main=paste("Region",r))
    ilist <- as.numeric(row.names(a)) 
    for(i in 1:length(ilist)){
      symbols(yrs , rep(ilist[i], length(yrs)), sqrt(a[i,]) / 60, add = T, inches =FALSE)
      if(trunc(i/5) == i/5){
      lines(c(1950,2009), rep(i,2), lty=3)}
      }
    }
  title(switch(fc,"Offshore","Distant water"),outer=T,cex=1)
  savePlot(paste("plot hbf by region by yr",switch(fc,"OS","DW")),type="png")
  }

# Plot sets proportion by line type by hbf by yr by region by fishing category
for (fc in c(1,2)){
  windows(height=14,width=12); par(mfrow=c(3,2),mar=c(3,4,2,1),oma=c(1,1,3,1)); r<-3; hh <- 18
  hbf_grps <- list(c(5:9),10:13, 14:17,18:19,20:25)
  a<- op[op$newfishingcat==fc,]
  for (r in 1:6) {
    plot(1,1,xlim=c(1994,2010),ylim=c(0,1),xlab="Year",ylab="Proportion Nylon",type="n",main=paste("Region",r))
    a1 <- a[a$reg==r,]
    for (hh in 1:5) {
      a2 <- a1[a1$hbf %in% hbf_grps[[hh]],]
      ayr <- floor(a2$op_yr/2)*2
      a2 <- table(a2$mainline,ayr)
      if(sum(a2)>1000) {
        prop <- a2[1,]/(a2[1,] + a2[2,])
        lines(as.numeric(names(prop)),prop,col=hh)
        }
      }
    }
  legend("bottomright",legend=hbf_grps,col=1:5,lty=1)
  title(switch(fc,"Offshore","Distant water"),outer=T,cex=1)
  savePlot(paste("plot prop sets by hbf region and yr",switch(fc,"OS","DW")),type="png")
  }


# Proportion zero by HBF by region by year
countzero <- function(a) {
  z <- sum(a==0)
  tot <- length(a)
  pz <- z/tot ; return(pz)
  }

windows(height=14,width=12); par(mfrow=c(3,2),mar=c(4,4,2,2))
for (r in 1:6) {
  a <- op[op$reg==r,]
  yrs <- sort(unique(a$op_yr))
  a <- tapply(a$bet,list(a$hbf,a$op_yr),countzero)
  plot(1,1, type="n", xlab="Year", ylab="HBF", ylim = c(1,25), xlim=c(1952, 2010),main=paste("Region",r))
  ilist <- as.numeric(row.names(a)) 
  for(i in 1:length(ilist)){
    symbols(yrs , rep(ilist[i], length(yrs)), sqrt(a[i,]), add = T, inches =FALSE)
    }
  }
savePlot("plot pzero bet by hbf by region by yr",type="png")

windows(height=14,width=12); par(mfrow=c(3,2),mar=c(4,4,2,2),oma=c(0,0,3,0))
for(fc in 1:2) {
  for (r in 1:6) {
    a <- op[op$reg==r & op$newfishingcat == fc,]
    yrs <- sort(unique(a$op_yr))
    a <- tapply(a$bet,list(a$hbf,a$op_yr),countzero)
    plot(1,1, type="n", xlab="Year", ylab="HBF", ylim = c(1,25), xlim=c(1952, 2010),main=paste("Region",r))
    ilist <- as.numeric(row.names(a)) 
    for(i in 1:length(ilist)){
      symbols(yrs , rep(ilist[i], length(yrs)), sqrt(a[i,]), add = T, inches =FALSE)
      }
    }
  title(switch(fc,"Offshore","Distant water"),outer=T,cex=1)
  savePlot(paste("plot pzero bet by hbf by region by yr",fc),type="png")
}

windows(height=14,width=12); par(mfrow=c(3,2),mar=c(4,4,2,2))
for (r in 1:6) {
  a <- op[op$reg==r,]
  yrs <- sort(unique(a$op_yr))
  a <- tapply(a$yft,list(a$hbf,a$op_yr),countzero)
  plot(1,1, type="n", xlab="Year", ylab="HBF", ylim = c(1,25), xlim=c(1952, 2010),main=paste("Region",r))
  ilist <- as.numeric(row.names(a)) 
  for(i in 1:length(ilist)){
    symbols(yrs , rep(ilist[i], length(yrs)), sqrt(a[i,]), add = T, inches =FALSE)
    }
  }
savePlot("plot pzero yft by hbf by region by yr",type="png")

windows(height=14,width=12); par(mfrow=c(3,2),mar=c(4,4,2,2),oma=c(0,0,3,0))
for(fc in 1:2) {
  for (r in 1:6) {
    a <- op[op$reg==r & op$newfishingcat == fc,]
    yrs <- sort(unique(a$op_yr))
    a <- tapply(a$yft,list(a$hbf,a$op_yr),countzero)
    plot(1,1, type="n", xlab="Year", ylab="HBF", ylim = c(1,25), xlim=c(1952, 2010),main=paste("Region",r))
    ilist <- as.numeric(row.names(a)) 
    for(i in 1:length(ilist)){
      symbols(yrs , rep(ilist[i], length(yrs)), sqrt(a[i,]), add = T, inches =FALSE)
      }
    }
  title(switch(fc,"Offshore","Distant water"),outer=T,cex=1)
  savePlot(paste("plot pzero yft by hbf by region by yr",fc),type="png")
}

windows(height=14,width=12); par(mfrow=c(3,2),mar=c(4,4,2,2))
for (r in 1:6) {
  a <- op[op$reg==r,]
  yrs <- sort(unique(a$op_yr))
  a <- tapply(a$alb,list(a$hbf,a$op_yr),countzero)
  plot(1,1, type="n", xlab="Year", ylab="HBF", ylim = c(1,25), xlim=c(1952, 2010),main=paste("Region",r))
  ilist <- as.numeric(row.names(a)) 
  for(i in 1:length(ilist)){
    symbols(yrs , rep(ilist[i], length(yrs)), sqrt(a[i,]), add = T, inches =FALSE)
    }
  }
savePlot("plot pzero alb by hbf by region by yr",type="png")

windows(height=14,width=12); par(mfrow=c(3,2),mar=c(4,4,2,2),oma=c(0,0,3,0))
for(fc in 1:2) {
  for (r in 1:6) {
    a <- op[op$reg==r & op$newfishingcat == fc,]
    yrs <- sort(unique(a$op_yr))
    a <- tapply(a$alb,list(a$hbf,a$op_yr),countzero)
    plot(1,1, type="n", xlab="Year", ylab="HBF", ylim = c(1,25), xlim=c(1952, 2010),main=paste("Region",r))
    ilist <- as.numeric(row.names(a)) 
    for(i in 1:length(ilist)){
      symbols(yrs , rep(ilist[i], length(yrs)), sqrt(a[i,]), add = T, inches =FALSE)
      }
    }
  title(switch(fc,"Offshore","Distant water"),outer=T,cex=1)
  savePlot(paste("plot pzero alb by hbf by region by yr",fc),type="png")
}

windows(height=14,width=12); par(mfrow=c(3,2),mar=c(4,4,2,2))
for (r in 1:6) {
  a <- op[op$reg==r,]
  yrs <- sort(unique(a$op_yr))
  a <- tapply(a$swo,list(a$hbf,a$op_yr),countzero)
  plot(1,1, type="n", xlab="Year", ylab="HBF", ylim = c(1,25), xlim=c(1952, 2010),main=paste("Region",r))
  ilist <- as.numeric(row.names(a)) 
  for(i in 1:length(ilist)){
    symbols(yrs , rep(ilist[i], length(yrs)), sqrt(a[i,]), add = T, inches =FALSE)
    }
  }
savePlot("plot pzero swo by hbf by region by yr",type="png")

windows(height=14,width=12); par(mfrow=c(3,2),mar=c(4,4,2,2),oma=c(0,0,3,0))
for(fc in 1:2) {
  for (r in 1:6) {
    a <- op[op$reg==r & op$newfishingcat == fc,]
    yrs <- sort(unique(a$op_yr))
    a <- tapply(a$swo,list(a$hbf,a$op_yr),countzero)
    plot(1,1, type="n", xlab="Year", ylab="HBF", ylim = c(1,25), xlim=c(1952, 2010),main=paste("Region",r))
    ilist <- as.numeric(row.names(a)) 
    for(i in 1:length(ilist)){
      symbols(yrs , rep(ilist[i], length(yrs)), sqrt(a[i,]), add = T, inches =FALSE)
      }
    }
  title(switch(fc,"Offshore","Distant water"),outer=T,cex=1)
  savePlot(paste("plot pzero swo by hbf by region by yr",fc),type="png")
}

# Proportion zero by HBF by subregion R3 by year
plot_pzero <- function(a,sp) {
  yrs <- sort(unique(a$op_yr))
  a <- tapply(a[,sp],list(a$hbf,a$op_yr),countzero)
  plot(1,1, type="n", xlab="Year", ylab="HBF", ylim = c(1,25), xlim=c(1950, 2010),main=paste("Subregion",r))
  ilist <- as.numeric(row.names(a)) 
  for(i in 1:length(ilist)){
    symbols(yrs , rep(ilist[i], length(yrs)), sqrt(a[i,]), add = T, inches =FALSE)
    }
  }

windows(height=14,width=12); par(mfrow=c(2,2),mar=c(4,4,2,2))
for (r in c(3.1, 3.2, 3.3, 3.4)) {
  a <- op[op$subreg==r,]
  plot_pzero(a,"bet")
  }
savePlot("plot pzero bet by hbf by subreg R3 by yr",type="png")

windows(height=14,width=12); par(mfrow=c(2,2),mar=c(4,4,2,2))
for (r in c(3.1, 3.2, 3.3, 3.4)) {
  a <- op[op$subreg==r,]
  plot_pzero(a,"yft")
  }
savePlot("plot pzero yft by hbf by subreg R3 by yr",type="png")

# Proportion of zero catches by HBF and year
plot_pzero_ll <- function(a,sp,la,lo) {
 # if(la == -30) browser()
  yrs <- sort(unique(a$op_yr))
  a <- tapply(a[,sp],list(a$hbf,a$op_yr),countzero)
  plot(1,1, type="n", xlab="Year", ylab="HBF", ylim = c(1,23), xlim=c(1955, 2010),main=paste(la+2.5,lo+2.5,sep=", "))
  ilist <- as.numeric(row.names(a)) 
  if(length(ilist) > 0) {
    for(i in 1:length(ilist)){
      symbols(yrs , rep(ilist[i], length(yrs)), sqrt(a[i,]), add = T, inches =FALSE)
      }
    }
  }
windows(height=14,width=24); par(mfrow=c(6,8),mar=c(2,2,1,0))
for (la in seq(15,-10,by=-5)) {
  for (lo in seq(130,165,by=5)) {
    a <- op[op$lat5==la & op$lon5 %in% c(lo),]
    if(dim(a)[[1]] >0)  plot_pzero_ll(a,"bet",la,lo) else plot(1,1,type="n")
    }
  }
savePlot("plot pzero bet by hbf by R3 latlong2",type="png")

windows(height=14,width=24); par(mfrow=c(6,8),mar=c(2,2,1,0))
for (la in seq(15,-10,by=-5)) {
  for (lo in seq(130,165,by=5)) {
    a <- op[op$lat5==la & op$lon5 %in% c(lo),]
    if(dim(a)[[1]] >0)  plot_pzero_ll(a,"yft",la,lo) else plot(1,1,type="n")
    }
  }
savePlot("plot pzero yft by hbf by R3 latlong2",type="png")

windows(height=14,width=24); par(mfrow=c(6,8),mar=c(2,2,1,0))
for (la in seq(15,-10,by=-5)) {
  for (lo in seq(130,165,by=5)) {
    a <- op[op$lat5==la & op$lon5 %in% c(lo),]
    if(dim(a)[[1]] >0)  plot_pzero_ll(a,"alb",la,lo) else plot(1,1,type="n")
    }
  }
savePlot("plot pzero alb by hbf by R3 latlong2",type="png")
#R4
windows(height=14,width=24); par(mfrow=c(6,8),mar=c(2,2,1,0))
for (la in seq(15,-10,by=-5)) {
  for (lo in seq(170,205,by=5)) {
    a <- op[op$lat5==la & op$lon5 %in% c(lo),]
    if(dim(a)[[1]] >0)  plot_pzero_ll(a,"bet",la,lo) else plot(1,1,type="n")
    }
  }
savePlot("plot pzero bet by hbf by R4 latlong2",type="png")

windows(height=14,width=24); par(mfrow=c(6,8),mar=c(2,2,1,0))
for (la in seq(15,-10,by=-5)) {
  for (lo in seq(170,205,by=5)) {
    a <- op[op$lat5==la & op$lon5 %in% c(lo),]
    if(dim(a)[[1]] >0)  plot_pzero_ll(a,"yft",la,lo) else plot(1,1,type="n")
    }
  }
savePlot("plot pzero yft by hbf by R4 latlong2",type="png")

windows(height=14,width=24); par(mfrow=c(6,8),mar=c(2,2,1,0))
for (la in seq(15,-10,by=-5)) {
  for (lo in seq(170,205,by=5)) {
    a <- op[op$lat5==la & op$lon5 %in% c(lo),]
    if(dim(a)[[1]] >0)  plot_pzero_ll(a,"alb",la,lo) else plot(1,1,type="n")
    }
  }
savePlot("plot pzero alb by hbf by R4 latlong2",type="png")

#R 1,2,5,6
pzero_hbf_sp_yq <- function(sp,laseq,loseq,fname,ti="") {
  windows(height=14,width=24); par(mfrow=c(length(laseq),length(loseq)),mar=c(2,2,1,0),oma=c(1,1,4,1))
  b <- op[op$lat5 %in% laseq & op$lon5 %in% loseq,]
  for (la in laseq) {
    for (lo in loseq) {
      a <- b[b$lat5==la & b$lon5 %in% c(lo),]
      if(dim(a)[[1]] >0)  plot_pzero_ll(a,sp,la,lo) else plot(1,1,type="n")
      }
    }
  title(ti,outer=T)
  savePlot(fname,type="png")
  }
pzero_hbf_sp_yq("bet",laseq=seq(35,20,-5),loseq=seq(130,165,5),fname="plot pzero bet by hbf by R1 latlong",ti="Bigeye Region 1")
pzero_hbf_sp_yq("yft",laseq=seq(35,20,-5),loseq=seq(130,165,5),fname="plot pzero yft by hbf by R1 latlong",ti="Yellowfin Region 1")
pzero_hbf_sp_yq("alb",laseq=seq(35,20,-5),loseq=seq(130,165,5),fname="plot pzero alb by hbf by R1 latlong",ti="Albacore Region 1")
  
pzero_hbf_sp_yq("bet",laseq=seq(35,20,-5),loseq=seq(170,205,5),fname="plot pzero bet by hbf by R2 latlong",ti="Bigeye Region 2")
pzero_hbf_sp_yq("yft",laseq=seq(35,20,-5),loseq=seq(170,205,5),fname="plot pzero yft by hbf by R2 latlong",ti="Yellowfin Region 2")
pzero_hbf_sp_yq("alb",laseq=seq(35,20,-5),loseq=seq(170,205,5),fname="plot pzero alb by hbf by R2 latlong",ti="Albacore Region 2")
  
pzero_hbf_sp_yq("bet",laseq=seq(-15,-35,-5),loseq=seq(140,165,5),fname="plot pzero bet by hbf by R5 latlong",ti="Bigeye Region 5")
pzero_hbf_sp_yq("yft",laseq=seq(-15,-35,-5),loseq=seq(140,165,5),fname="plot pzero yft by hbf by R5 latlong",ti="Yellowfin Region 5")
pzero_hbf_sp_yq("alb",laseq=seq(-15,-35,-5),loseq=seq(140,165,5),fname="plot pzero alb by hbf by R5 latlong",ti="Albacore Region 5")
  
pzero_hbf_sp_yq("bet",laseq=seq(-15,-35,-5),loseq=seq(170,205,5),fname="plot pzero bet by hbf by R6 latlong",ti="Bigeye Region 6")
pzero_hbf_sp_yq("yft",laseq=seq(-15,-35,-5),loseq=seq(170,205,5),fname="plot pzero yft by hbf by R6 latlong",ti="Yellowfin Region 6")
pzero_hbf_sp_yq("alb",laseq=seq(-15,-35,-5),loseq=seq(170,205,5),fname="plot pzero alb by hbf by R6 latlong",ti="Albacore Region 6")
  
#R 1,2,3,4,5,6
load("alldatraw.RData")
dat <- dataclean(dat,allHBF=T)
dat <- dataprep(dat)
op <- dat[dat$reg > 0 & dat$reg <=6,c("tonnage","fishingcat","ncrew","target","mainline","branchline","op_yr","op_mon","op_day","lat","lon","hbf",
          "hooks","bet","yft","alb","lat5","lon5","reg","subreg","vessid","yrqtr","latlong","cstart_yr","cstart_mon","cstart_day")]

plot_hbf_5x5 <- function(a,la,lo) {
  yrs <- sort(unique(a$op_yr))
  a <- tapply(a$hooks,list(a$hbf,a$op_yr),sum)
  plot(1,1, type="n", xlab="Year", ylab="HBF", ylim = c(1,25), xlim=c(1955, 2010),main=paste(la+2.5,lo+2.5,sep=", "))
  ilist <- as.numeric(row.names(a)) 
  if(length(ilist) > 0) {
    for(i in 1:length(ilist)){
      symbols(yrs , rep(ilist[i], length(yrs)), sqrt(a[i,])/400, add = T, inches =F)
    }
  }
}
hbf_sp_yq <- function(sp,laseq,loseq,fname,ti="") {
  windows(height=14,width=24); par(mfrow=c(length(laseq),length(loseq)),mar=c(2,2,1,0),oma=c(1,1,4,1))
  b <- op[op$lat5 %in% laseq & op$lon5 %in% loseq,]
  for (la in laseq) {
    for (lo in loseq) {
      a <- b[b$lat5==la & b$lon5 %in% c(lo),]
      if(dim(a)[[1]] >0)  plot_hbf_5x5(a,la,lo) else plot(1,1,type="n")
      }
    }
  title(ti,outer=T)
  savePlot(fname,type="png")
  }
hbf_sp_yq(laseq=seq(35,20,-5),loseq=seq(130,165,5),fname="plot hbf by R1 latlong",ti="HBF Region 1")
hbf_sp_yq(laseq=seq(35,20,-5),loseq=seq(170,205,5),fname="plot hbf by R2 latlong",ti="HBF Region 2")
hbf_sp_yq(laseq=seq(15,-10,-5),loseq=seq(130,165,5),fname="plot hbf by R3 latlong",ti="HBF Region 3")
hbf_sp_yq(laseq=seq(15,-10,-5),loseq=seq(170,205,5),fname="plot hbf by R4 latlong",ti="HBF Region 4")
hbf_sp_yq(laseq=seq(-15,-35,-5),loseq=seq(140,165,5),fname="plot hbf by R5 latlong",ti="HBF Region 5")
hbf_sp_yq(laseq=seq(-15,-35,-5),loseq=seq(170,205,5),fname="plot hbf by R6 latlong",ti="HBF Region 6")

op <- op[op$hbf>4,]
hbf_sp_yq(laseq=seq(35,20,-5),loseq=seq(130,165,5),fname="plot hbf noHBF4 R1 latlong noHBF4",ti="HBF Region 1")
hbf_sp_yq(laseq=seq(35,20,-5),loseq=seq(170,205,5),fname="plot hbf noHBF4 R2 latlong",ti="HBF Region 2")
hbf_sp_yq(laseq=seq(15,-10,-5),loseq=seq(130,165,5),fname="plot hbf noHBF4 R3 latlong",ti="HBF Region 3")
hbf_sp_yq(laseq=seq(15,-10,-5),loseq=seq(170,205,5),fname="plot hbf noHBF4 R4 latlong",ti="HBF Region 4")
hbf_sp_yq(laseq=seq(-15,-35,-5),loseq=seq(140,165,5),fname="plot hbf noHBF4 R5 latlong",ti="HBF Region 5")
hbf_sp_yq(laseq=seq(-15,-35,-5),loseq=seq(170,205,5),fname="plot hbf noHBF4 R6 latlong",ti="HBF Region 6")
  
load("op.RData")

# Proportion of zero catches
plot_pzero_both <- function(a,la,lo) {
  cx <- 0.8
  yrs <- sort(unique(a$yrqtr))
  alb <- tapply(a$alb,list(a$yrqtr),countzero)
  yft <- tapply(a$yft,list(a$yrqtr),countzero)
  bet <- tapply(a$bet,list(a$yrqtr),countzero)
#  swo <- tapply(a$swo,list(a$yrqtr),countzero)
  plot(1,1, type="n", xlab="Year", ylab="p(zero catch)", ylim = c(0,1), xlim=c(1952, 2010),main=paste(la+2.5,", ",lo+2.5,sep=""))
  points(yrs, yft,cex=cx)
  points(yrs,bet,col=2,pch=2,cex=cx)
  points(yrs,alb,col=3,pch=3,cex=cx)
#  points(yrs,swo,col=4,pch=4,cex=cx)
  }
pzero_all_sp_yq <- function(laseq,loseq,fname,ti="") {
  windows(height=14,width=24); par(mfrow=c(length(laseq),length(loseq)),mar=c(2,2,1,0),oma=c(1,1,4,1))
  b <- op[op$lat5 %in% laseq & op$lon5 %in% loseq,]
  for (la in laseq) {
    for (lo in loseq) {
      a <- b[b$lat5==la & b$lon5 %in% c(lo),]
      if(dim(a)[[1]] >0)  plot_pzero_both(a,la,lo) else plot(1,1,type="n")
      }
    }
  title(ti,outer=T)
  savePlot(fname,type="png")
  }
  
pzero_all_sp_yq(laseq=seq(35,20,-5),loseq=seq(130,165,5),fname="plot pzero allsp by R1 latlong",ti="Probability of zero catch Region 1")
pzero_all_sp_yq(laseq=seq(35,20,-5),loseq=seq(170,205,5),fname="plot pzero allsp by R2 latlong",ti="Probability of zero catch Region 2")
pzero_all_sp_yq(laseq=seq(15,-10,-5),loseq=seq(130,165,5),fname="plot pzero allsp by R3 latlong",ti="Probability of zero catch Region 3")
pzero_all_sp_yq(laseq=seq(15,-10,-5),loseq=seq(170,205,5),fname="plot pzero allsp by R4 latlong",ti="Probability of zero catch Region 4")
pzero_all_sp_yq(laseq=seq(-15,-35,-5),loseq=seq(140,165,5),fname="plot pzero allsp by R5 latlong",ti="Probability of zero catch Region 5")
pzero_all_sp_yq(laseq=seq(-15,-35,-5),loseq=seq(170,205,5),fname="plot pzero allsp by R6 latlong",ti="Probability of zero catch Region 6")


# Albacore catch distribution histograms
write.csv(table(op$hbf,floor(op$yrqtr/5)*5,op$reg),file="hbf by region by 5 years.csv")
windows();par(mfrow=c(3,2))
for(i in 1:6) {
  a <- op[op$reg==i,]
  hist(a$alb,xlab="Albacore catch",xlim=c(0,250),main=paste("Region",i))
}
savePlot("hist alb by reg",type="png")

windows(height=14,width=24); par(mfrow=c(6,8),mar=c(2,2,1,0))
for (la in seq(15,-10,by=-5)) {
  for (lo in seq(130,165,by=5)) {
    a <- op[op$lat5==la & op$lon5 %in% c(lo,lo+5),]
    hist(a$alb,xlab="Albacore catch",xlim=c(0,250),main=paste(la,lo,sep="."))
    }
  }
savePlot("hist alb R3 latlong",type="png")
windows(height=14,width=24); par(mfrow=c(6,8),mar=c(2,2,1,0))
for (la in seq(15,-10,by=-5)) {
  for (lo in seq(170,205,by=5)) {
    a <- op[op$lat5==la & op$lon5 %in% c(lo,lo+5),]
    hist(a$alb,xlab="Albacore catch",xlim=c(0,250),main=paste(la,lo,sep="."))
    }
  }
savePlot("hist alb R4 latlong",type="png")

# Proportions of CPUE over 20 per 1000 hooks
count_over_n <- function(a,n) {
  z <- sum(a>n)
  tot <- length(a)
  pz <- z/tot ; return(pz)
  }
plot_pcpue <- function(a,la,lo,cpue) {
  cx<-0.9
  yrs <- sort(unique(a$yrqtr))
  alb <- tapply(a$alb/a$hooks,list(a$yrqtr),count_over_n,cpue/1000)
  yft <- tapply(a$yft/a$hooks,list(a$yrqtr),count_over_n,cpue/1000)
  bet <- tapply(a$bet/a$hooks,list(a$yrqtr),count_over_n,cpue/1000)
  plot(1,1, type="n", xlab="Year", ylab="p(zero catch)", ylim = c(0,1), xlim=c(1952, 2010),main=paste(la+2.5,lo+2.5,sep=", "))
  points(yrs,yft,col=1,pch=1,cex=cx)
  points(yrs,bet,col=2,pch=2,cex=cx)
  points(yrs, alb,col=3,pch=3,cex=cx)
  }
pcpue_all_sp_yq <- function(cpue,laseq,loseq,fname,ti="") {
  windows(height=14,width=24); par(mfrow=c(length(laseq),length(loseq)),mar=c(2,2,1,0),oma=c(1,1,4,1))
  b <- op[op$lat5 %in% laseq & op$lon5 %in% loseq,]
  for (la in laseq) {
    for (lo in loseq) {
      a <- b[b$lat5==la & b$lon5 %in% c(lo),]
      if(dim(a)[[1]] >0)  plot_pcpue(a,la,lo,cpue) else plot(1,1,type="n")
      }
    }
  title(ti,outer=T)
  savePlot(fname,type="png")
  }
  
pcpue_all_sp_yq(cpue=10,laseq=seq(35,20,-5),loseq=seq(130,165,5),fname="plot pcpue allsp by R1 latlong",ti="Probability of cpue > 1/100 hooks Region 1")
pcpue_all_sp_yq(cpue=10,laseq=seq(35,20,-5),loseq=seq(170,205,5),fname="plot pcpue allsp by R2 latlong",ti="Probability of cpue > 1/100 hooks Region 2")
pcpue_all_sp_yq(cpue=10,laseq=seq(15,-10,-5),loseq=seq(130,165,5),fname="plot pcpue allsp by R3 latlong",ti="Probability of cpue > 1/100 hooks Region 3")
pcpue_all_sp_yq(cpue=10,laseq=seq(15,-10,-5),loseq=seq(170,205,5),fname="plot pcpue allsp by R4 latlong",ti="Probability of cpue > 1/100 hooks Region 4")
pcpue_all_sp_yq(cpue=10,laseq=seq(-15,-35,-5),loseq=seq(140,165,5),fname="plot pcpue allsp by R5 latlong",ti="Probability of cpue > 1/100 hooks Region 5")
pcpue_all_sp_yq(cpue=10,laseq=seq(-15,-35,-5),loseq=seq(170,205,5),fname="plot pcpue allsp by R6 latlong",ti="Probability of cpue > 1/100 hooks Region 6")


# Catch 
plot_catch <- function(a,la,lo) {
  cx <- 0.9
  yrs <- sort(unique(a$op_yr))
  yft <- tapply(a$yft,list(a$op_yr),sum)
  bet <- tapply(a$bet,list(a$op_yr),sum)
  alb <- tapply(a$alb,list(a$op_yr),sum)
  plot(1,1, type="n", xlab="Year", ylab="Catch (numbers)", ylim = c(0,100000), xlim=c(1952, 2010),main=paste(la+2.5,lo+2.5,sep=", "))
  lines(yrs,yft,col=1,cex=cx)
  lines(yrs,bet,col=2,lty=2,cex=cx)
  lines(yrs,alb,col=3,lty=3,cex=cx)
  }
catch_all_sp_yq <- function(laseq,loseq,fname,ti="") {
  windows(height=14,width=24); par(mfrow=c(length(laseq),length(loseq)),mar=c(2,2,1,0),oma=c(1,1,4,1))
  b <- op[op$lat5 %in% laseq & op$lon5 %in% loseq,]
  for (la in laseq) {
    for (lo in loseq) {
      a <- b[b$lat5==la & b$lon5 %in% c(lo),]
      if(dim(a)[[1]] >0)  plot_catch(a,la,lo) else plot(1,1,type="n")
      }
    }
  legend("topright",legend=c("Yellowfin","Bigeye","Albacore"),lty=c(1:3),col=c(1:3))
  title(ti,outer=T)
  savePlot(fname,type="png")
  }
  
catch_all_sp_yq(laseq=seq(35,20,-5),loseq=seq(130,165,5),fname="plot catch allsp by R1 latlong",ti="Catch Region 1")
catch_all_sp_yq(laseq=seq(35,20,-5),loseq=seq(170,205,5),fname="plot catch allsp by R2 latlong",ti="Catch Region 2")
catch_all_sp_yq(laseq=seq(15,-10,-5),loseq=seq(130,165,5),fname="plot catch allsp by R3 latlong",ti="Catch Region 3")
catch_all_sp_yq(laseq=seq(15,-10,-5),loseq=seq(170,205,5),fname="plot catch allsp by R4 latlong",ti="Catch Region 4")
catch_all_sp_yq(laseq=seq(-15,-35,-5),loseq=seq(140,165,5),fname="plot catch allsp by R5 latlong",ti="Catch Region 5")
catch_all_sp_yq(laseq=seq(-15,-35,-5),loseq=seq(170,205,5),fname="plot catch allsp by R6 latlong",ti="Catch Region 6")

# Effort by fishing category  
plot_effort_fc <- function(a,la,lo) {
  yrs <- sort(unique(a$op_yr))
  eff <- tapply(a$hooks,list(factor(a$newfishingcat,levels=c(1:2)),a$op_yr),sum,na.rm=T)
  os_eff <- eff[1,]
  dw_eff <- eff[2,]
  plot(1,1, type="n", xlab="Year", ylab="Effort (hooks)", ylim = c(0,5e6), xlim=c(1952, 2010),main=paste(la+2.5,lo+2.5,sep=", "))
  lines(yrs,os_eff,col=1,lty=1)
  lines(yrs,dw_eff,col=2,lty=2)
  }
effort_all_fc_yq <- function(laseq,loseq,fname,ti="") {
  windows(height=14,width=24); par(mfrow=c(length(laseq),length(loseq)),mar=c(2,2,1,0),oma=c(1,1,4,1))
  b <- op[op$lat5 %in% laseq & op$lon5 %in% loseq,]
  for (la in laseq) {
    for (lo in loseq) {
      a <- b[b$lat5==la & b$lon5 %in% c(lo),]
      if(dim(a)[[1]] >0)  plot_effort_fc(a,la,lo) else plot(1,1,type="n")
      }
    }
  legend("topright",legend=c("Offshore","Distant water"),lty=c(1,2),col=c(1,2))
  title(ti,outer=T)
  savePlot(fname,type="png")
  }
  
effort_all_fc_yq(laseq=seq(35,20,-5),loseq=seq(130,165,5),fname="plot effort by fc R1 latlong",ti="Effort Region 1")
effort_all_fc_yq(laseq=seq(35,20,-5),loseq=seq(170,205,5),fname="plot effort by fc R2 latlong",ti="Effort Region 2")
effort_all_fc_yq(laseq=seq(15,-10,-5),loseq=seq(130,165,5),fname="plot effort by fc R3 latlong",ti="Effort Region 3")
effort_all_fc_yq(laseq=seq(15,-10,-5),loseq=seq(170,205,5),fname="plot effort by fc R4 latlong",ti="Effort Region 4")
effort_all_fc_yq(laseq=seq(-15,-35,-5),loseq=seq(140,165,5),fname="plot effort by fc R5 latlong",ti="Effort Region 5")
effort_all_fc_yq(laseq=seq(-15,-35,-5),loseq=seq(170,205,5),fname="plot effort by fc R6 latlong",ti="Effort Region 6")


# Median CPUE
plot_median_cpue <- function(a,la,lo) {
  cx=0.9
  yrs <- sort(unique(a$yrqtr))
  alb <- tapply(a$alb/a$hooks,list(a$yrqtr),median)
  yft <- tapply(a$yft/a$hooks,list(a$yrqtr),median)
  bet <- tapply(a$bet/a$hooks,list(a$yrqtr),median)
  plot(1,1, type="n", xlab="Year", ylab="p(zero catch)", ylim = c(0,0.04), xlim=c(1952, 2010),main=paste(la+2.5,lo+2.5,sep=", "))
  points(yrs,yft,col=1,pch=1,cex=cx)
  points(yrs,bet,col=2,pch=2,cex=cx)
  points(yrs, alb,col=3,pch=3,cex=cx)
  }
median_all_sp_yq_fc <- function(cpue,laseq,loseq,fname,ti="",fc="both") {
  windows(height=14,width=24); par(mfrow=c(length(laseq),length(loseq)),mar=c(2,2,1,0),oma=c(1,1,4,1))
  b <- op[op$lat5 %in% laseq & op$lon5 %in% loseq,]
  b <- b[b$newfishingcat %in% switch(fc,"both"=c(1,2),"OS"=1,"DW"=2),]
  for (la in laseq) {
    for (lo in loseq) {
      a <- b[b$lat5==la & b$lon5 %in% c(lo),]
#  browser()
      if(dim(a)[[1]] >0)  plot_median_cpue(a,la,lo) else plot(1,1,type="n")
      }
    }
  title(ti,outer=T)
  savePlot(fname,type="png")
  }
  
for(fcc in c("OS","DW","both")) {
  median_all_sp_yq_fc(laseq=seq(35,20,-5),loseq=seq(130,165,5),fname=paste("plot median allsp by R1 latlong",fcc),ti=paste("Median cpue Region 1",fcc),fc=fcc)
  median_all_sp_yq_fc(laseq=seq(35,20,-5),loseq=seq(170,205,5),fname=paste("plot median allsp by R2 latlong",fcc),ti=paste("Median cpue Region 2",fcc),fc=fcc)
  median_all_sp_yq_fc(laseq=seq(15,-10,-5),loseq=seq(130,165,5),fname=paste("plot median allsp by R3 latlong",fcc),ti=paste("Median cpue Region 3",fcc),fc=fcc)
  median_all_sp_yq_fc(laseq=seq(15,-10,-5),loseq=seq(170,205,5),fname=paste("plot median allsp by R4 latlong",fcc),ti=paste("Median cpue Region 4",fcc),fc=fcc)
  median_all_sp_yq_fc(laseq=seq(-15,-35,-5),loseq=seq(140,165,5),fname=paste("plot median allsp by R5 latlong",fcc),ti=paste("Median cpue Region 5",fcc),fc=fcc)
  median_all_sp_yq_fc(laseq=seq(-15,-35,-5),loseq=seq(170,205,5),fname=paste("plot median allsp by R6 latlong",fcc),ti=paste("Median cpue Region 6",fcc),fc=fcc)
}

# Mean CPUE
plot_mean_cpue <- function(a,la,lo) {
  cx=0.9
  yrs <- sort(unique(a$yrqtr))
  alb <- tapply(a$alb/a$hooks,list(a$yrqtr),mean)
  yft <- tapply(a$yft/a$hooks,list(a$yrqtr),mean)
  bet <- tapply(a$bet/a$hooks,list(a$yrqtr),mean)
  plot(1,1, type="n", xlab="Year", ylab="p(zero catch)", ylim = c(0,0.04), xlim=c(1952, 2010),main=paste(la+2.5,lo+2.5,sep=", "))
  points(yrs,yft,col=1,pch=1,cex=cx)
  points(yrs,bet,col=2,pch=2,cex=cx)
  points(yrs, alb,col=3,pch=3,cex=cx)
  }
mean_all_sp_yq_fc <- function(cpue,laseq,loseq,fname,ti="",fc="both") {
  windows(height=14,width=24); par(mfrow=c(length(laseq),length(loseq)),mar=c(2,2,1,0),oma=c(1,1,4,1))
  b <- op[op$lat5 %in% laseq & op$lon5 %in% loseq,]
  b <- b[b$newfishingcat %in% switch(fc,"both"=c(1,2),"OS"=1,"DW"=2),]
  for (la in laseq) {
    for (lo in loseq) {
      a <- b[b$lat5==la & b$lon5 %in% c(lo),]
      if(dim(a)[[1]] >0)  plot_mean_cpue(a,la,lo) else plot(1,1,type="n")
      }
    }
  title(ti,outer=T)
  savePlot(fname,type="png")
  }
  
for(fcc in c("OS","DW","both")) {
  mean_all_sp_yq_fc(laseq=seq(35,20,-5),loseq=seq(130,165,5),fname=paste("plot mean allsp by R1 latlong",fcc),ti=paste("Mean cpue Region 1",fcc),fc=fcc)
  mean_all_sp_yq_fc(laseq=seq(35,20,-5),loseq=seq(170,205,5),fname=paste("plot mean allsp by R2 latlong",fcc),ti=paste("Mean cpue Region 2",fcc),fc=fcc)
  mean_all_sp_yq_fc(laseq=seq(15,-10,-5),loseq=seq(130,165,5),fname=paste("plot mean allsp by R3 latlong",fcc),ti=paste("Mean cpue Region 3",fcc),fc=fcc)
  mean_all_sp_yq_fc(laseq=seq(15,-10,-5),loseq=seq(170,205,5),fname=paste("plot mean allsp by R4 latlong",fcc),ti=paste("Mean cpue Region 4",fcc),fc=fcc)
  mean_all_sp_yq_fc(laseq=seq(-15,-35,-5),loseq=seq(140,165,5),fname=paste("plot mean allsp by R5 latlong",fcc),ti=paste("Mean cpue Region 5",fcc),fc=fcc)
  mean_all_sp_yq_fc(laseq=seq(-15,-35,-5),loseq=seq(170,205,5),fname=paste("plot mean allsp by R6 latlong",fcc),ti=paste("Mean cpue Region 6",fcc),fc=fcc)
  graphics.off()
}

# SWO and SWO targeting
plot_median_cpue_swo <- function(a,la,lo) {
  cx=0.9
  yrs <- sort(unique(a$yrqtr))
  alb <- tapply(a$alb/a$hooks,list(a$yrqtr),median)
  yft <- tapply(a$yft/a$hooks,list(a$yrqtr),median)
  bet <- tapply(a$bet/a$hooks,list(a$yrqtr),median)
  swo <- tapply(a$swo/a$hooks,list(a$yrqtr),median)
  plot(1,1, type="n", xlab="Year", ylab="p(zero catch)",xlim=c(1952, 2010),ylim=c(0,.02),main=paste(la+2.5,lo+2.5,sep=", "))
  points(yrs,yft,col=1,pch=1,cex=cx)
  points(yrs,bet,col=2,pch=2,cex=cx)
#  points(yrs,alb,col=3,pch=3,cex=cx)
  points(yrs,swo,col=4,pch=4,cex=cx)
  }
median_all_spswo_yq_fc <- function(cpue,laseq,loseq,fname,ti="",fc="both") {
  windows(height=14,width=24); par(mfrow=c(length(laseq),length(loseq)),mar=c(2,2,1,0),oma=c(1,1,4,1))
  b <- op[op$lat5 %in% laseq & op$lon5 %in% loseq,]
  b <- b[b$newfishingcat %in% switch(fc,"both"=c(1,2),"OS"=1,"DW"=2),]
  for (la in laseq) {
    for (lo in loseq) {
      a <- b[b$lat5==la & b$lon5 %in% c(lo),]
      if(dim(a)[[1]] >0)  plot_median_cpue_swo(a,la,lo) else plot(1,1,type="n")
      }
    }
  title(ti,outer=T)
  savePlot(fname,type="png")
  }
  
for(fcc in c("OS","DW","both")) {
  median_all_spswo_yq_fc(laseq=seq(35,20,-5),loseq=seq(130,165,5),fname=paste("plot median allspswo by R1 latlong",fcc),ti=paste("Median cpue Region 1",fcc),fc=fcc)
  median_all_spswo_yq_fc(laseq=seq(35,20,-5),loseq=seq(170,205,5),fname=paste("plot median allspswo by R2 latlong",fcc),ti=paste("Median cpue Region 2",fcc),fc=fcc)
  median_all_spswo_yq_fc(laseq=seq(-15,-35,-5),loseq=seq(140,165,5),fname=paste("plot median allspswo by R5 latlong",fcc),ti=paste("Median cpue Region 5",fcc),fc=fcc)
  median_all_spswo_yq_fc(laseq=seq(-15,-35,-5),loseq=seq(170,205,5),fname=paste("plot median allspswo by R6 latlong",fcc),ti=paste("Median cpue Region 6",fcc),fc=fcc)
  graphics.off()
}

a <- op[op$lat5 %in% c(25,30) & op$lon5 %in% seq(130,165,5),]
windows(height=14,width=24); par(mfrow=c(3,4))
for (y in seq(1950,2009,5)) {
  b <- a[a$op_yr %in% seq(y,y+4,1),]
  plot(b$bet,b$swo,xlim=c(0,150),ylim=c(0,60),cex=0.6,xlab="Bigeye",ylab="Swordfish",main=paste(y,"-",y+4))
  text(140,50,cor(b$bet,b$swo))
  }
savePlot("plot R1 bet v swo sets",type="png")
a <- op[op$lat5 %in% c(25,30) & op$lon5 %in% seq(130,165,5),]
windows(height=14,width=24); par(mfrow=c(3,4))
for (y in seq(1950,2009,5)) {
  b <- a[a$op_yr %in% seq(y,y+4,1) & a$hbf %in% c(3,4),]
  plot(b$bet,b$swo,xlim=c(0,150),ylim=c(0,60),cex=0.6,xlab="Bigeye",ylab="Swordfish",main=paste(y,"-",y+4))
  text(140,50,cor(b$bet,b$swo))
  }
savePlot("plot R1 bet v swo sets HBF3_4",type="png")
a <- op[op$lat5 %in% c(25,30) & op$lon5 %in% seq(130,165,5),]
windows(height=14,width=24); par(mfrow=c(3,4))
for (y in seq(1950,2009,5)) {
  b <- a[a$op_yr %in% seq(y,y+4,1) & a$hbf >=5,]
  plot(b$bet,b$swo,xlim=c(0,150),ylim=c(0,60),cex=0.6,xlab="Bigeye",ylab="Swordfish",main=paste(y,"-",y+4))
  text(140,50,cor(b$bet,b$swo))
  }
savePlot("plot R1 bet v swo sets HBF5_up",type="png")
a <- op[op$lat5 %in% c(25,30) & op$lon5 %in% seq(180,205,5),]
windows(height=14,width=24); par(mfrow=c(3,4))
for (y in seq(1950,2009,5)) {
  b <- a[a$op_yr %in% seq(y,y+4,1),]
  plot(b$bet,b$swo,xlim=c(0,150),ylim=c(0,60),cex=0.6,xlab="Bigeye",ylab="Swordfish",main=paste(y,"-",y+4))
  text(140,50,cor(b$bet,b$swo))
  }
savePlot("plot R2east bet v swo sets",type="png")
windows(height=14,width=24); par(mfrow=c(3,4))
for (y in seq(1950,2009,5)) {
  b <- a[a$op_yr %in% seq(y,y+4,1),]
  plot(b$bet,b$alb,xlim=c(0,150),ylim=c(0,250),cex=0.6,xlab="Bigeye",ylab="Albacore",main=paste(y,"-",y+4))
  text(140,240,cor(b$bet,b$alb))
  }
savePlot("plot R2east bet v alb sets",type="png")
windows(height=14,width=24); par(mfrow=c(3,4))
for (y in seq(1950,2009,5)) {
  b <- a[a$op_yr %in% seq(y,y+4,1),]
  plot(b$alb,b$swo,xlim=c(0,250),ylim=c(0,60),cex=0.6,xlab="Albacore",ylab="Swordfish",main=paste(y,"-",y+4))
  text(200,50,cor(b$alb,b$swo))
  }
savePlot("plot R2east alb v swo sets",type="png")
windows(height=14,width=24); par(mfrow=c(3,4))
for (y in seq(1950,2009,5)) {
  b <- a[a$op_yr %in% seq(y,y+4,1),]
  hist(b$swo,breaks=c(0:9,seq(10,1000,5)),ylim=c(0,1),xlim=c(0,50),main=paste(y,"-",y+4))
  }
savePlot("plot R2east swo hist sets",type="png")
windows(height=14,width=24); par(mfrow=c(3,4))
for (y in seq(1950,2009,5)) {
  b <- a[a$op_yr %in% seq(y,y+4,1) & a$hbf<5,]
  plot(b$bet,b$swo,xlim=c(0,150),ylim=c(0,60),cex=0.6,xlab="Bigeye",ylab="Swordfish",main=paste(y,"-",y+4))
  }
savePlot("plot R2east bet v swo sets hbf3-4",type="png")
windows(height=14,width=24); par(mfrow=c(3,4))
for (y in seq(1950,2009,5)) {
  b <- a[a$op_yr %in% seq(y,y+4,1) & a$hbf>5,]
  plot(b$bet,b$swo,xlim=c(0,150),ylim=c(0,60),cex=0.6,xlab="Bigeye",ylab="Swordfish",main=paste(y,"-",y+4))
  }
savePlot("plot R2east bet v swo sets hbf 6-up",type="png")
windows(height=14,width=24); par(mfrow=c(3,4))
for (y in seq(1950,2009,5)) {
  b <- a[a$op_yr %in% seq(y,y+4,1),]
  hist(b$hooks,xlim=c(0,6000),main=paste(y,"-",y+4))
  }
savePlot("plot R2east hooks 5yr",type="png")
windows(height=14,width=24); par(mfrow=c(3,4))
for (y in seq(1950,2009,5)) {
  b <- a[a$op_yr %in% seq(y,y+4,1),]
  plot(b$hooks,b$swo,xlim=c(0,6000),cex=0.6,main=paste(y,"-",y+4))
  b2 <- tapply(b$swo,250*(floor(b$hooks/250)),mean)
  lines(as.numeric(names(b2)),b2,lwd=2,col=2)
  }
savePlot("plot R2east hooks v swo sets",type="png")
for(r in 1:6) {
  a <- op[op$reg==r,]
  windows(height=14,width=24); par(mfrow=c(3,4),oma=c(0,0,1,0))
  for (y in seq(1950,2009,5)) {
    b <- a[a$op_yr %in% seq(y,y+4,1),]
    hist(b$hooks,xlim=c(0,6000),breaks=seq(0,8000,100),main=paste(y,"-",y+4))
    }
  title(paste("Region",r),outer=T)
  savePlot(paste("plot hooks R",r,sep=""),type="png")
  }
for(fc in 1:2) {
  for(r in 1:6) {
    a <- op[op$reg==r & op$newfishingcat==fc,]
    windows(height=14,width=24); par(mfrow=c(3,4),oma=c(0,0,1,0))
    for (y in seq(1950,2009,5)) {
      b <- a[a$op_yr %in% seq(y,y+4,1),]
      hist(b$hooks,xlim=c(0,6000),breaks=seq(0,8000,100),main=paste(y,"-",y+4))
      }
    title(paste("Region",r,switch(fc,"OS","DW")),outer=T)
    savePlot(paste("plot hooks R",r,switch(fc," OS"," DW"),sep=""),type="png")
    }
  }
a <- op[op$lat5 %in% c(25,30) & op$lon5 %in% seq(185,205,5),]
a <- a[a$op_yr < 1970,]
table(a$latlong,100*a$op_yr + a$op_mon + a$op_day/100)


# Proportion of zero catches
plot_pzero_swo <- function(a,la,lo) {
  cx <- 0.8
  yrs <- sort(unique(a$yrqtr))
#  alb <- tapply(a$alb,list(a$yrqtr),countzero)
#  yft <- tapply(a$yft,list(a$yrqtr),countzero)
  bet <- tapply(a$bet,list(a$yrqtr),countzero)
  swo <- tapply(a$swo,list(a$yrqtr),countzero)
  plot(1,1, type="n", xlab="Year", ylab="p(zero catch)", ylim = c(0,1), xlim=c(1952, 2010),main=paste(la+2.5,", ",lo+2.5,sep=""))
#  points(yrs, yft,cex=cx)
  points(yrs,bet,col=2,pch=2,cex=cx)
#  points(yrs,alb,col=3,pch=3,cex=cx)
  points(yrs,swo,col=4,pch=4,cex=cx)
  }
pzero_fq_yq_swo <- function(laseq,loseq,fname,ti="",fc="both") {
  windows(height=14,width=24); par(mfrow=c(length(laseq),length(loseq)),mar=c(2,2,1,0),oma=c(1,1,4,1))
  b <- op[op$lat5 %in% laseq & op$lon5 %in% loseq,]
  b <- b[b$newfishingcat %in% switch(fc,"both"=c(1,2),"OS"=1,"DW"=2),]
  for (la in laseq) {
    for (lo in loseq) {
      a <- b[b$lat5==la & b$lon5 %in% c(lo),]
      if(dim(a)[[1]] >0)  plot_pzero_swo(a,la,lo) else plot(1,1,type="n")
      }
    }
  title(ti,outer=T)
  savePlot(fname,type="png")
  }
  
for(fcc in c("OS","DW","both")) {
  pzero_fq_yq_swo(laseq=seq(35,20,-5),loseq=seq(130,165,5),fname=paste("plot pzero allspswo by R1 latlong",fcc),ti=paste("Probability of zero catch Region 1",fcc),fc=fcc)
  pzero_fq_yq_swo(laseq=seq(35,20,-5),loseq=seq(170,205,5),fname=paste("plot pzero allspswo by R2 latlong",fcc),ti=paste("Probability of zero catch Region 2",fcc),fc=fcc)
#  pzero_fq_yq_swo(laseq=seq(15,-10,-5),loseq=seq(130,165,5),fname=paste("plot pzero allspswo by R3 latlong",fcc),ti=paste("Probability of zero catch Region 3",fcc),fc=fcc)
#  pzero_fq_yq_swo(laseq=seq(15,-10,-5),loseq=seq(170,205,5),fname=paste("plot pzero allspswo by R4 latlong",fcc),ti=paste("Probability of zero catch Region 4",fcc),fc=fcc)
  pzero_fq_yq_swo(laseq=seq(-15,-35,-5),loseq=seq(140,165,5),fname=paste("plot pzero allspswo by R5 latlong",fcc),ti=paste("Probability of zero catch Region 5",fcc),fc=fcc)
  pzero_fq_yq_swo(laseq=seq(-15,-35,-5),loseq=seq(170,205,5),fname=paste("plot pzero allspswo by R6 latlong",fcc),ti=paste("Probability of zero catch Region 6",fcc),fc=fcc)
  graphics.off()
}



setwd("D:/SimonHOYLE/AddDataforAlbacore/indices"); runsp="bet"; maxqtrs=200; addmain<-F; addbranch<-F;addother=F;addalb=T ;minqtrs_byreg = c(4,4,10,10,4,2)
for(runreg in c(1,2,3,4,5,6)) {
  minqtrs <- minqtrs_byreg[runreg]
  for(runsp in c("bet")) {
    fnamedelta <- paste(runsp,"_R",runreg," indices HBF10 ",minqtrs,"-",maxqtrs,"_qtrs ","deltabin",ifelse(addother,"_othersp",""),ifelse(addalb,"_albcat",""),sep="")
    fnamepos <-   paste(runsp,"_R",runreg," indices HBF10 ",minqtrs,"-",maxqtrs,"_qtrs ","deltapos",ifelse(addother,"_othersp",""),ifelse(addalb,"_albcat",""),sep="")
    combine_delta_xl_indices(fnamedelta,fnamepos)
    }
  }

setwd("D:/SimonHOYLE/AddDataforAlbacore/albcat"); runsp="bet"; maxqtrs=200; addmain<-F; addbranch<-F;addother=F;addalb=T ;maxqtrs=200; minqtrs_byreg = c(2,2,8,8,2,2,2,8)
for(runreg in c(1,2,3,4,5,6)) {
  minqtrs <- minqtrs_byreg[runreg]
  for(runsp in c("bet")) {
    fnamedelta <- paste(runsp,"_R",runreg,"eq ",minqtrs,"-",maxqtrs,"_qtrs ","deltabin",ifelse(addother,"_othersp",""),ifelse(addalb,"_albcat",""),sep="")
    fnamepos <-   paste(runsp,"_R",runreg,"eq ",minqtrs,"-",maxqtrs,"_qtrs ","deltapos",ifelse(addother,"_othersp",""),ifelse(addalb,"_albcat",""),sep="")
    combine_delta_xl(fnamedelta,fnamepos)
    }
  }
for(runreg in c(3)) {
  minqtrs <- 9
  for(runsp in c("bet")) {
    fnamedelta <- paste(runsp,"_R",runreg,"all ",minqtrs,"-",maxqtrs,"_qtrs ","deltabin",ifelse(addother,"_othersp",""),ifelse(addalb,"_albcat",""),sep="")
    fnamepos <-   paste(runsp,"_R",runreg,"all ",minqtrs,"-",maxqtrs,"_qtrs ","deltapos",ifelse(addother,"_othersp",""),ifelse(addalb,"_albcat",""),sep="")
    combine_delta_xl(fnamedelta,fnamepos)
    }
  }


setwd("D:/SimonHOYLE/AddDataforAlbacore/HBF10"); runsp="bet"; maxqtrs=200; addmain<-F; addbranch<-F;addother=F;addalb=F ;maxqtrs=200; minqtrs_byreg = c(2,2,8,8,2,2,2,8)
for(runreg in c(1,2,3,4,5,6)) {
  minqtrs <- minqtrs_byreg[runreg]
  for(runsp in c("bet")) {
    fnamedelta <- paste(runsp,"_R",runreg,"eq HBF10 ",minqtrs,"-",maxqtrs,"_qtrs ","deltabin",ifelse(addother,"_othersp",""),ifelse(addalb,"_albcat","")," (2)",sep="")
    fnamepos <-   paste(runsp,"_R",runreg,"eq HBF10 ",minqtrs,"-",maxqtrs,"_qtrs ","deltapos",ifelse(addother,"_othersp",""),ifelse(addalb,"_albcat","")," (2)",sep="")
    combine_delta_xl(fnamedelta,fnamepos)
    }
  }
for(runreg in c(3)) {
  minqtrs <- 9
  for(runsp in c("bet")) {
    fnamedelta <- paste(runsp,"_R",runreg,"all HBF10 ",minqtrs,"-",maxqtrs,"_qtrs ","deltabin",ifelse(addother,"_othersp",""),ifelse(addalb,"_albcat",""),sep="")
    fnamepos <-   paste(runsp,"_R",runreg,"all HBF10 ",minqtrs,"-",maxqtrs,"_qtrs ","deltapos",ifelse(addother,"_othersp",""),ifelse(addalb,"_albcat",""),sep="")
    combine_delta_xl(fnamedelta,fnamepos)
    }
  }


setwd("D:/SimonHOYLE/AddDataforAlbacore/std_eq_and_all"); runsp="bet"; maxqtrs=200; addmain<-F; addbranch<-F;addother=F;addalb=F ;maxqtrs=200; minqtrs_byreg = c(2,2,8,8,2,2,2,8)
for(runreg in c(1,2,3,4,5,6)) {
  minqtrs <- minqtrs_byreg[runreg]
  for(runsp in c("bet")) {
    fnamedelta <- paste(runsp,"_R",runreg,"eq ",minqtrs,"-",maxqtrs,"_qtrs ","deltabin",ifelse(addother,"_othersp",""),ifelse(addalb,"_albcat",""),sep="")
    fnamepos <-   paste(runsp,"_R",runreg,"eq ",minqtrs,"-",maxqtrs,"_qtrs ","deltapos",ifelse(addother,"_othersp",""),ifelse(addalb,"_albcat",""),sep="")
    combine_delta_xl(fnamedelta,fnamepos)
    }
  }
for(runreg in c(3)) {
  minqtrs <- 9
  for(runsp in c("bet")) {
    fnamedelta <- paste(runsp,"_R",runreg,"all ",minqtrs,"-",maxqtrs,"_qtrs ","deltabin",ifelse(addother,"_othersp",""),ifelse(addalb,"_albcat",""),sep="")
    fnamepos <-   paste(runsp,"_R",runreg,"all ",minqtrs,"-",maxqtrs,"_qtrs ","deltapos",ifelse(addother,"_othersp",""),ifelse(addalb,"_albcat",""),sep="")
    combine_delta_xl(fnamedelta,fnamepos)
    }
  }

