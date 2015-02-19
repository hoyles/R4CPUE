# working_2012.r
#setwd("D:/CPUE/WCPO2012/")
#setwd("C:/Users/simonh/Dropbox/Analyses/JapanCPUE/")
require(splines)
require(R4MFCL)
source("./Rfiles/support_functions.r")

nms <- c("op_yr","op_mon","op_day","lat","latcode","lon","loncode","jvname","callsign","tonnage","fishingcat",
      "licnse","target","mainline","branchline","bait","hbf","hooks","alb","bet","yft","swo","trip_st")
wdths <- c(4,2,2,2,1,3,1,20,6,7,1,5,1,1,1,1,3,6,3,3,3,3,8)
cc <- c("integer","integer","integer","integer","integer","integer","integer","character","character","real",
        "integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","date")
cbind(nms,wdths,cc)

#(logical, integer, numeric, complex, character, raw), or "factor", "Date" or "POSIXct"
read.fwf("./data/WCPFC5293op+cruiseID.dat",widths=wdths,col.names=nms,buffersize=200000,stringsAsFactors=F,n=20);gc()
d1 <- read.fwf("./data/WCPFC5293op+cruiseID.dat",widths=wdths,col.names=nms,buffersize=200000,stringsAsFactors=F);gc()
read.fwf("./data/WCPFC9411op+cruiseID.dat",widths=wdths,col.names=nms,buffersize=200000,stringsAsFactors=F,n=20);gc()
d2 <- read.fwf("./data/WCPFC9411op+cruiseID.dat",widths=wdths,col.names=nms,buffersize=200000,stringsAsFactors=F);gc()
dat <- rbind(d1,d2)

save(d1,file="raw5293.RData")
save(d2,file="raw9411.RData")
save(dat,file="alldatraw.RData")
load("raw5293.RData")
load("raw9411.RData")
load("alldatraw.RData")
dat$lat <- as.numeric(dat$lat)
dat$hbf <- as.numeric(dat$hbf)
dat$hooks <- as.numeric(dat$hooks)
dat$alb <- as.numeric(dat$alb)
dat$bet <- as.numeric(dat$bet)
dat$yft <- as.numeric(dat$yft)
dat$swo <- as.numeric(dat$swo)
table(dat$lat,useNA="always")
table(dat$hbf,useNA="always")
table(dat$hooks,useNA="always")
table(dat$alb,useNA="always")
table(dat$bet,useNA="always")
table(dat$yft,useNA="always")
table(dat$swo,useNA="always")
table(dat[dat$alb=="  .",]$bet)
a <- dat[dat$alb=="  .",]
summary(a)
table(a$bet)
table(a$swo)
table(a$yft)
dat[is.na(dat$alb)==T,]$alb <- 0
dat[is.na(dat$bet)==T,]$bet <- 0
dat[is.na(dat$yft)==T,]$yft <- 0
dat[is.na(dat$swo)==T,]$swo <- 0
rm(d1,d2)
gc()
a <- table(dat$op_yr)
write.table(a,"table sets per year.txt")
save(dat,file="alldat2.RData")

# check tonnage - many sets have NA for tonnage
tapply(dat$tonnage,list(dat$op_yr,dat$fishingcat),mean,na.rm=T)

#####################  START HERE
# Prepare data
#memory.limit(size=4000)
require(R4MFCL)
setwd("D:/CPUE/WCPO2012/")
require(splines)
source("./Rfiles/support_functions.r")
load("alldatraw.RData")
dim(dat)
dat <- dataclean(dat)
dim(dat)
summary(dat)
dat <- dataprep(dat,alldat=T)
dim(dat)
save(dat,file="dat_prepared.RData")
#load("dat_prepared.RData")

op <- dat[,c("tonnage","newfishingcat","jvname","target","mainline","branchline","op_yr","op_mon","op_day","lat","lon","hbf","bait",
 "hooks","bet","yft","alb","swo","lat5","lon5","reg","subreg","vessid","yrqtr","latlong","tripid","trip_yr","sets_per_trip")]
save(op,file="op.RData")

#load(file="op.RData")
rm(d1,d2,dat,a,cc,wdths,nms);gc()

windows(13,13); par(mfrow=c(3,2))
windows(13,13); par(mfrow=c(3,2))
windows(13,13); par(mfrow=c(3,2))
devlist <- dev.list()
for (i in 1:6) {
  a <- op[op$reg==i,]
  bet <- tapply(a$bet,a$yrqtr,sum)
  yft <- tapply(a$yft,a$yrqtr,sum)
  alb <- tapply(a$alb,a$yrqtr,sum)
  swo <- tapply(a$swo,a$yrqtr,sum)
  yq <- as.numeric(names(bet))
  effort <- tapply(a$hooks,a$yrqtr,sum)
  dev.set(2)
  plot(yq,bet,col=2,ylim=c(0,max(c(bet,yft,alb,swo))))
  points(yq,yft,col=1)
  points(yq,alb,col=3)
  points(yq,swo,col=4)
  dev.set(3)
  plot(yq,effort)
  dev.set(4)
  plot(yq,bet/effort,col=2,ylim=c(0,max(c(bet/effort,yft/effort,alb/effort,swo/effort))))
  points(yq,yft/effort,col=1)
  points(yq,alb/effort,col=3)
  points(yq,swo/effort,col=4)
  }         

cbind(table(dat$op_yr))
cbind(table(d1$op_yr))
dat[dat$vessid==1,][1:100,]
unique(dat$callsign[dat$vessid==1])
unique(dat$callsign[dat$vessid==2])
table(dat$trip_st)[1:100]
unique(dat$trip_st)

#####################################################################
# Plot catch by species 1x1 x decade
a5 <- aggregate(cbind(bet,yft,alb,swo,hooks) ~ lon + lat + eval(10*floor((op_yr+5)/10)-5),data=op,FUN=sum)
a <- aggregate(cbind(bet,yft,alb,swo,hooks) ~ lon + lat + eval(10*floor((op_yr)/10)),data=op,FUN=sum)
names(a)[3] <- names(a5)[3] <- "decade"

plot_catchmap <- function(indat=a,vbl,dcd,latlim=c(-40,40),lonlim=c(120,210),brk=seq(0,1,.05),brk2=seq(0,1,.1),ti="") {
  plot(1:5,1:5,ylim=latlim,xlim=lonlim,type="n",xlab="Longitude",ylab="Latitude")
  indat <- cbind(indat,vbl)
  indat <- indat[indat$lon <= 210,]
  a1 <- with(indat[indat$decade==dcd,],tapply(vbl,list(lon,lat),sum))
  image(as.numeric(rownames(a1)),as.numeric(colnames(a1)),a1,add=T,col=heat.colors(length(brk)-1),breaks=brk)
  contour(as.numeric(rownames(a1)),as.numeric(colnames(a1)),a1,add=T,levels=brk2)
  map(database="world2Hires",add=T)
  title(paste(dcd,ti))
  }
plot_cpuemap <- function(indat=a,vb1,vb2,dcd,latlim=c(-40,40),lonlim=c(120,210),brk=seq(0,1,.05),brk2=seq(0,1,.1),ti="") {
  plot(1:5,1:5,ylim=latlim,xlim=lonlim,type="n",xlab="Longitude",ylab="Latitude")
  indat <- cbind(indat,vb1,vb2)
  indat <- indat[indat$lon <= 210,]
  a1 <- with(indat[indat$decade==dcd,],tapply(vb1,list(lon,lat),sum)/tapply(vb2,list(lon,lat),sum))
  image(as.numeric(rownames(a1)),as.numeric(colnames(a1)),a1,add=T,col=heat.colors(length(brk)-1),breaks=brk)
  contour(as.numeric(rownames(a1)),as.numeric(colnames(a1)),a1,add=T,levels=brk2)
  map(database="world2Hires",add=T)
  title(paste(dcd,ti))
  }

windows(width=20,height=15);par(mfrow=c(2,3))
for(d in seq(1955,2005,10)) plot_catchmap(indat=a,vbl=a$bet/(a$bet+a$yft),dcd=d)
savePlot("PropBET in YFT_BET5",type="png")
for(d in seq(1955,2005,10)) plot_catchmap(indat=a,vbl=a$bet/(a$bet+a$yft+a$alb+a$swo),dcd=d)
savePlot("PropBET in YFT_BET_ALB_SWO5",type="png")
for(d in seq(1955,2005,10)) plot_catchmap(indat=a,vbl=a$bet/(a$bet+a$alb),dcd=d)
savePlot("PropBET in ALB_BET5",type="png")

for(d in seq(1955,2005,10)) with(a5,plot_catchmap(indat=a5,vbl=alb/(bet+yft+alb+swo),dcd=d,latlim=c(-10,20)))
savePlot("PropALB in YFT_BET_ALB_SWO core5",type="png")
for(d in seq(1950,2000,10)) plot_catchmap(indat=a,vbl=a$alb/(a$bet+a$yft+a$alb+a$swo),dcd=d,latlim=c(-10,20))
savePlot("PropALB in YFT_BET_ALB_SWO core",type="png")
for(d in seq(1955,2005,10)) with(a5,plot_catchmap(indat=a5,vbl=alb/(bet+yft+alb+swo),dcd=d))
savePlot("PropALB in YFT_BET_ALB_SWO5",type="png")
for(d in seq(1950,2000,10)) plot_catchmap(indat=a,vbl=a$alb/(a$bet+a$yft+a$alb+a$swo),dcd=d)
savePlot("PropALB in YFT_BET_ALB_SWO",type="png")

for(d in seq(1955,2005,10)) with(a5,plot_catchmap(indat=a5,vbl=alb/(yft+alb),dcd=d,latlim=c(-10,20)))
savePlot("PropALB in YFT_ALB core5",type="png")
for(d in seq(1950,2000,10)) plot_catchmap(indat=a,vbl=a$alb/(a$yft+a$alb),dcd=d,latlim=c(-10,20))
savePlot("PropALB in YFT_ALB core",type="png")
for(d in seq(1955,2005,10)) with(a5,plot_catchmap(indat=a5,vbl=alb/(yft+alb),dcd=d))
savePlot("PropALB in YFT_ALB5",type="png")
for(d in seq(1950,2000,10)) plot_catchmap(indat=a,vbl=a$alb/(a$yft+a$alb),dcd=d)
savePlot("PropALB in YFT_ALB",type="png")

for(d in seq(1955,2005,10)) plot_catchmap(indat=a,vbl=a$swo/(a$bet+a$yft+a$alb+a$swo),dcd=d)
savePlot("PropSWO in YFT_BET_ALB_SWO5",type="png")

# windows(width=20,height=15);par(mfrow=c(2,3))
for(d in seq(1955,2005,10)) with(a5,plot_catchmap(indat=a5,vbl=yft/(bet+yft+alb+swo),dcd=d))
savePlot("PropYFT in YFT_BET_ALB_SWO5",type="png")
for(d in seq(1950,2000,10)) with(a,plot_catchmap(indat=a,vbl=yft/(bet+yft+alb+swo),dcd=d))
savePlot("PropYFT in YFT_BET_ALB_SWO",type="png")
for(d in seq(1955,2005,10)) with(a5,plot_catchmap(indat=a5,vbl=yft/(bet+yft+alb+swo),dcd=d,latlim=c(-10,20)))
savePlot("PropYFT in YFT_BET_ALB_SWO core5",type="png")
for(d in seq(1950,2000,10)) with(a,plot_catchmap(indat=a,vbl=yft/(bet+yft+alb+swo),dcd=d,latlim=c(-10,20)))
savePlot("PropYFT in YFT_BET_ALB_SWO core",type="png")
for(d in seq(1955,2005,10)) with(a5,plot_catchmap(indat=a5,vbl=yft/(yft+alb),dcd=d))
savePlot("PropYFT in YFT_ALB5",type="png")
for(d in seq(1950,2000,10)) with(a,plot_catchmap(indat=a,vbl=yft/(yft+alb),dcd=d))
savePlot("PropYFT in YFT_ALB",type="png")
for(d in seq(1955,2005,10)) with(a5,plot_catchmap(indat=a5,vbl=yft/(yft+alb),dcd=d,latlim=c(-10,20)))
savePlot("PropYFT in YFT_ALB core5",type="png")
for(d in seq(1950,2000,10)) with(a,plot_catchmap(indat=a,vbl=yft/(yft+alb),dcd=d,latlim=c(-10,20)))
savePlot("PropYFT in YFT_ALB core",type="png")
for(d in seq(1955,2005,10)) with(a5,plot_catchmap(indat=a5,vbl=yft/(yft+bet),dcd=d))
savePlot("PropYFT in YFT_BET5",type="png")
for(d in seq(1950,2000,10)) with(a,plot_catchmap(indat=a,vbl=yft/(yft+bet),dcd=d))
savePlot("PropYFT in YFT_BET",type="png")
for(d in seq(1955,2005,10)) with(a5,plot_catchmap(indat=a5,vbl=yft/(yft+bet),dcd=d,latlim=c(-10,20)))
savePlot("PropYFT in YFT_BET core5",type="png")
for(d in seq(1950,2000,10)) with(a,plot_catchmap(indat=a,vbl=yft/(yft+bet),dcd=d,latlim=c(-10,20)))
savePlot("PropYFT in YFT_BET core",type="png")

for(d in seq(1955,2005,10)) plot_catchmap(indat=a,vbl=a$bet/(a$bet+a$yft),dcd=d,latlim=c(-10,20))
savePlot("PropBET in YFT_BET core5",type="png")

# windows(width=20,height=15);par(mfrow=c(2,3))
for(d in seq(1955,2005,10)) plot_catchmap(indat=a,vbl=a$bet/(a$bet+a$alb),dcd=d,latlim=c(-10,20))
savePlot("PropBET in ALB_BET core5",type="png")

# windows(width=20,height=15);par(mfrow=c(2,3))
for(d in seq(1955,2005,10)) plot_catchmap(indat=a,vbl=a$bet/(a$bet+a$yft+a$alb+a$swo),dcd=d,latlim=c(-10,20))
savePlot("PropBET in YFT_BET_ALB_SWO core5",type="png")

#####
#Catch & effort
# windows(width=20,height=15);par(mfrow=c(2,3))
brk <- 10^(seq(1,6,.5)); brk2 <- floor(10^(seq(1,6,.5))) 
for(d in seq(1955,2005,10)) plot_catchmap(indat=a,vbl=a$yft,dcd=d,latlim=c(-40,40),brk=brk,brk2=brk2)
savePlot("CatchYFT5",type="png")

# windows(width=20,height=15);par(mfrow=c(2,3))
for(d in seq(1955,2005,10)) plot_catchmap(indat=a,vbl=a$bet,dcd=d,latlim=c(-40,40),brk=brk,brk2=brk2)
savePlot("CatchBET5",type="png")

# windows(width=20,height=15);par(mfrow=c(2,3))
for(d in seq(1955,2005,10)) plot_catchmap(indat=a,vbl=a$alb,dcd=d,latlim=c(-40,40),brk=brk,brk2=brk2)
savePlot("CatchALB5",type="png")

# windows(width=20,height=15);par(mfrow=c(2,3))
for(d in seq(1955,2005,10)) plot_catchmap(indat=a,vbl=a$swo,dcd=d,latlim=c(-40,40),brk=brk,brk2=brk2)
savePlot("CatchSWO5",type="png")

# Core
# windows(width=20,height=15);par(mfrow=c(2,3))
for(d in seq(1955,2005,10)) plot_catchmap(indat=a,vbl=a$yft,dcd=d,latlim=c(-10,20),brk=brk,brk2=brk2)
savePlot("CatchYFT core5",type="png")

# windows(width=20,height=15);par(mfrow=c(2,3))
for(d in seq(1955,2005,10)) plot_catchmap(indat=a,vbl=a$bet,dcd=d,latlim=c(-10,20),brk=brk,brk2=brk2)
savePlot("CatchBET core5",type="png")

# windows(width=20,height=15);par(mfrow=c(2,3))
for(d in seq(1955,2005,10)) plot_catchmap(indat=a,vbl=a$alb,dcd=d,latlim=c(-10,20),brk=brk,brk2=brk2)
savePlot("CatchALB core5",type="png")

ebrk <- 10^(seq(1,7,.5)); ebrk2 <- floor(10^(seq(1,7,.5))) 
for(d in seq(1955,2005,10)) with(a5,plot_catchmap(indat=a5,vbl=hooks,dcd=d,latlim=c(-10,20),brk=ebrk,brk2=ebrk2))
savePlot("map effort core5",type="png")
for(d in seq(1950,2000,10)) with(a,plot_catchmap(indat=a,vbl=hooks,dcd=d,latlim=c(-10,20),brk=ebrk,brk2=ebrk2))
savePlot("map effort core",type="png")
for(d in seq(1955,2005,10)) with(a5,plot_catchmap(indat=a5,vbl=hooks,dcd=d,brk=ebrk,brk2=ebrk2))
savePlot("map effort5",type="png")
for(d in seq(1950,2000,10)) with(a,plot_catchmap(indat=a,vbl=hooks,dcd=d,brk=ebrk,brk2=ebrk2))
savePlot("map effort",type="png")

# nominal_CPUE
cpbrk <- seq(0,9,.5); cpbrk2 <- seq(0,20,1) 
for(d in seq(1955,2005,10)) with(a5,plot_cpuemap(indat=a5,vb1=alb,vb2=hooks/100,dcd=d,brk=cpbrk,brk2=cpbrk2))
savePlot("map_ALB_CPUE5",type="png")
cpbrk <- seq(0,3,.1); cpbrk2 <- seq(0,4,.2) 
for(d in seq(1955,2005,10)) with(a5,plot_cpuemap(indat=a5,vb1=bet,vb2=hooks/100,dcd=d,brk=cpbrk,brk2=cpbrk2))
savePlot("map_BET_CPUE5",type="png")
cpbrk <- seq(0,6,.1); cpbrk2 <- seq(0,4,.5) 
for(d in seq(1955,2005,10)) with(a5,plot_cpuemap(indat=a5,vb1=yft,vb2=hooks/100,dcd=d,brk=cpbrk,brk2=cpbrk2))
savePlot("map_YFT_CPUE5",type="png")
cpbrk <- seq(0,9,.5); cpbrk2 <- seq(0,20,1) 
for(d in seq(1950,2000,10)) with(a,plot_cpuemap(indat=a,vb1=alb,vb2=hooks/100,dcd=d,brk=cpbrk,brk2=cpbrk2))
savePlot("map_ALB_CPUE",type="png")
cpbrk <- seq(0,3,.1); cpbrk2 <- seq(0,4,.2) 
for(d in seq(1950,2000,10)) with(a,plot_cpuemap(indat=a,vb1=bet,vb2=hooks/100,dcd=d,brk=cpbrk,brk2=cpbrk2))
savePlot("map_BET_CPUE",type="png")
cpbrk <- seq(0,6,.1); cpbrk2 <- seq(0,4,.5) 
for(d in seq(1950,2000,10)) with(a,plot_cpuemap(indat=a,vb1=yft,vb2=hooks/100,dcd=d,brk=cpbrk,brk2=cpbrk2))
savePlot("map_YFT_CPUE",type="png")
cpbrk <- seq(0,9,.5); cpbrk2 <- seq(0,20,1)
 
for(d in seq(1955,2005,10)) with(a5,plot_cpuemap(indat=a5,vb1=alb,vb2=hooks/100,dcd=d,latlim=c(-10,20),brk=cpbrk,brk2=cpbrk2))
savePlot("map_ALB_CPUE_core5",type="png")
cpbrk <- seq(0,3,.1); cpbrk2 <- seq(0,4,.2) 
for(d in seq(1955,2005,10)) with(a5,plot_cpuemap(indat=a5,vb1=bet,vb2=hooks/100,latlim=c(-10,20),dcd=d,brk=cpbrk,brk2=cpbrk2))
savePlot("map_BET_CPUE_core5",type="png")
cpbrk <- seq(0,6,.1); cpbrk2 <- seq(0,4,.5) 
for(d in seq(1955,2005,10)) with(a5,plot_cpuemap(indat=a5,vb1=yft,vb2=hooks/100,latlim=c(-10,20),dcd=d,brk=cpbrk,brk2=cpbrk2))
savePlot("map_YFT_CPUE_core5",type="png")
cpbrk <- seq(0,9,.5); cpbrk2 <- seq(0,20,1) 
for(d in seq(1950,2000,10)) with(a,plot_cpuemap(indat=a,vb1=alb,vb2=hooks/100,dcd=d,latlim=c(-10,20),brk=cpbrk,brk2=cpbrk2))
savePlot("map_ALB_CPUE_core",type="png")
cpbrk <- seq(0,3,.1); cpbrk2 <- seq(0,4,.2) 
for(d in seq(1950,2000,10)) with(a,plot_cpuemap(indat=a,vb1=bet,vb2=hooks/100,dcd=d,latlim=c(-10,20),brk=cpbrk,brk2=cpbrk2))
savePlot("map_BET_CPUE_core",type="png")
cpbrk <- seq(0,6,.1); cpbrk2 <- seq(0,4,.5) 
for(d in seq(1950,2000,10)) with(a,plot_cpuemap(indat=a,vb1=yft,vb2=hooks/100,dcd=d,latlim=c(-10,20),brk=cpbrk,brk2=cpbrk2))
savePlot("map_YFT_CPUE_core",type="png")

# by qtr
a5q <- aggregate(cbind(bet,yft,alb,swo,hooks) ~ lon + lat + eval(10*floor((op_yr+5)/10)-5)
        + eval(floor((op_mon-1)/3)+1),data=op,FUN=sum)
aq <- aggregate(cbind(bet,yft,alb,swo,hooks) ~ lon + lat + eval(10*floor((op_yr)/10)) 
        + eval(floor((op_mon-1)/3)+1),data=op,FUN=sum)
names(aq)[3:4] <- names(a5q)[3:4] <- c("decade","qtr")

windows(width=20,height=15);par(mfrow=c(2,3))
for (dat in c("aq","a5q")) {
  for(sp in c("alb","bet","yft","swo")) {
    x <- get(dat)
    for(qtr in 1:4) {
      a <- x[x$qtr==qtr,]
      p5 <- switch(dat,"aq"=0,"a5q"=5)
      for(d in seq(1950,2000,10)+p5) plot_catchmap(indat=a,vbl=a[,sp]/(a$bet+a$yft+a$alb+a$swo),dcd=d,ti=paste("q",qtr,sep=""))
      savePlot(paste("map_",sp,"prop_qtr",qtr,"d",switch(dat,aq="",a5q=5),sep=""),type="png")
      }
    }
  }
      
# by HBF
a5q <- aggregate(cbind(bet,yft,alb,swo,hooks) ~ lon + lat + eval(10*floor((op_yr+5)/10)-5)
        + cut(hbf,breaks=c(0,13,30),labels=c(1,2)),data=op,FUN=sum)
aq <- aggregate(cbind(bet,yft,alb,swo,hooks) ~ lon + lat + eval(10*floor((op_yr)/10)) 
        + cut(hbf,breaks=c(0,13,30),labels=c(1,2)),data=op,FUN=sum)
names(aq)[3:4] <- names(a5q)[3:4] <- c("decade","hbfclass")

windows(width=20,height=15);par(mfrow=c(2,3))
for (dat in c("aq")) {
  for(sp in c("alb","bet","yft","swo")) {
    x <- get(dat)
    for(hb in c("1","2")) {
      a <- x[x$hbfclass==hb,]
      p5 <- switch(dat,"aq"=0,"a5q"=5)
      for(d in seq(1950,2000,10)+p5) plot_catchmap(indat=a,vbl=a[,sp]/(a$bet+a$yft+a$alb+a$swo),dcd=d,
        ti=paste("hbf",switch(hb,"1"="1-13","2"="14+")))
      savePlot(paste("map_",sp,"prop_hbf",hb,"d",switch(dat,aq="",a5q=5),sep=""),type="png")
      }
    }
  }
      


#####################################################################
# start cluster stuff
require(R4MFCL)
setwd("D:/CPUE/WCPO2012/")
require(splines)
source("./Rfiles/support_functions.r")
load(file="op.RData")
length(unique(op$vessid))
ls()
dim(op)
library(rpart)
library(cluster)
require(pkpkg)
require(R4MFCL)
require(MASS)

windows()
hist(op$op_yr)
head(op)
table(op$vessid)
unique(paste(dat$vessid,dat$callsign))
unique(dat$callsign)
unique(dat$callsign[dat$vessid==1])

op[op$vessid==1,][1:10,]
op[op$vessid==2,][1:10,]
dat[dat$vessid==1,][1:10,]
dat[dat$vessid==1000,][,c("yrqtr","vessid","callsign","trip_st","tripid")]
windows();par(mfrow=c(6,6))
hist(op)
hist(dat)
table(dat$trip_yr)
head(op)

windows()# get the first set in every trip
boxplot(op$sets_per_trip ~ op$trip_yr,na.rm=T)
windows(width=10,height=15);par(mfrow=c(2,1))
with(op[op$newfishingcat==1,],boxplot(sets_per_trip ~ op_yr))
with(op[op$newfishingcat==2,],boxplot(sets_per_trip ~ op_yr))
windows(width=10,height=15);par(mfrow=c(2,1))
with(op[op$newfishingcat==1 & op$op_yr==2004,],hist(sets_per_trip,nclass=530))
with(op[op$newfishingcat==1 & op$op_yr==2005,],hist(sets_per_trip,nclass=530))
aggregate(sets_per_trip ~ op_yr,op[op$newfishingcat==1,],mean,na.rm=T)
aggregate(sets_per_trip ~ op_yr,op[op$newfishingcat==2,],mean,na.rm=T)


windows();par(mfrow=c(2,1))
boxplot(op$sets_per_trip~op$newfishingcat)
barplot(table(op$newfishingcat))

tp <- aggregate(cbind(trip_yr,sets_per_trip) ~ tripid + newfishingcat,data=op,max)
windows(16,16);par(mfrow=c(2,2),mar=c(3,4,2,1))
for(f in sort(unique(tp$newfishingcat))) {
  with(tp[tp$newfishingcat==f,],boxplot(sets_per_trip ~ trip_yr,main=switch(f,"Offshore","Distant water"),ylim=c(0,590),ylab = "Sets per trip"))
  }
for(f in sort(unique(tp$newfishingcat))) {
  a <- with(tp[tp$newfishingcat==f,],tapply(sets_per_trip,trip_yr,mean,na.rm=T))
  plot(as.numeric(names(a)),a,ylim=c(0,max(a)),ylab = "Mean sets per trip")
  }
savePlot("sets_per_trip_by_fishingcat.png",type="png")

# Some plotting to see what's going on. 
windows(width=20,height=18); par(mfrow=c(3,2),mar=c(2,2,1,1))
for(rg in 1:6) {
  reg <- op[op$reg==rg,]
  yq <- sort(unique(reg$op_yr))
  yq1 <- sort(unique(reg[reg$newfishingcat==1,]$op_yr)) 
  yq2 <- sort(unique(reg[reg$newfishingcat==2,]$op_yr))
  xl <- c(1,length(yq)); yl=c(0,6000) 
  boxplot(hooks ~ op_yr,data=reg[reg$newfishingcat==1,],axes=F,xlim=xl,ylim=yl,main=paste("Region",rg),at=match(yq1,yq)-0.2,boxwex=0.4)
  axis(1,at=1:length(yq),labels=yq)
  axis(2)
  boxplot(hooks ~ op_yr,data=reg[reg$newfishingcat==2,],axes=F,main=paste("Region",rg),at=match(yq2,yq)+0.2,add=T,boxwex=0.4,border=2)
  }
savePlot("hooks_per_set_by_fishingcat_reg.png",type="png")

windows(width=20,height=18); par(mfrow=c(3,2),mar=c(3,2,2,1))
a <- op[op$op_yr >= 2000,]
for(rg in 1:6) {
  reg <- a[a$reg==rg,]
  hk <- sort(unique(floor(reg$hooks/200)*200))
  xl <- c(1,length(hk)); yl=c(0,.1) 
  boxplot(eval(yft/hooks) ~ eval(floor(hooks/200)*200),data=reg,main=paste("Region",rg),ylim=yl,at=match(hk,hk)-.3,boxwex=.2)
  boxplot(eval(bet/hooks) ~ eval(floor(hooks/200)*200),data=reg,main=paste("Region",rg),at=match(hk,hk),boxwex=.2,border=2,add=T,xaxt="n")
  boxplot(eval(alb/hooks) ~ eval(floor(hooks/200)*200),data=reg,main=paste("Region",rg),at=match(hk,hk)+.3,boxwex=.2,border=3,add=T,xaxt="n")
  }
savePlot("CPUE_by_hooks boxplot",type="png")
  
windows(width=20,height=18); par(mfrow=c(3,2),mar=c(3,2,2,1))
a <- op[op$op_yr >= 1995,]
for(rg in 1:6) {
  reg <- a[a$reg==rg & a$newfishingcat==2,]
  a1 <- tapply(reg$yft/reg$hooks,eval(floor(reg$hooks/200)*200),mean,na.rm=T)/mean(reg$yft/reg$hooks)
  a2 <- tapply(reg$bet/reg$hooks,eval(floor(reg$hooks/200)*200),mean,na.rm=T)/mean(reg$bet/reg$hooks)
  a3 <- tapply(reg$alb/reg$hooks,eval(floor(reg$hooks/200)*200),mean,na.rm=T)/mean(reg$alb/reg$hooks)
  a4 <- tapply(reg$swo/reg$hooks,eval(floor(reg$hooks/200)*200),mean,na.rm=T)/mean(reg$swo/reg$hooks)
  plot(as.numeric(names(a1)),a1,ylim=c(0,4),type="l",ylab = "Mean sets per trip")
  lines(as.numeric(names(a2)),a2,col=2)
  lines(as.numeric(names(a3)),a3,col=3)
#  lines(as.numeric(names(a4)),a4,col=4)
  }
savePlot("CPUE_by_hooks_DW",type="png")

windows(width=20,height=18); par(mfrow=c(3,2),mar=c(3,2,2,1))
a <- op[op$op_yr >= 2000,]
for(rg in 1:6) {
  reg <- a[a$reg==rg & a$newfishingcat==1,]
  a1 <- tapply(reg$yft/reg$hooks,eval(floor(reg$hbf/2)*2),mean,na.rm=T)/mean(reg$yft/reg$hooks)
  a2 <- tapply(reg$bet/reg$hooks,eval(floor(reg$hbf/2)*2),mean,na.rm=T)/mean(reg$bet/reg$hooks)
  a3 <- tapply(reg$alb/reg$hooks,eval(floor(reg$hbf/2)*2),mean,na.rm=T)/mean(reg$alb/reg$hooks)
  a4 <- tapply(reg$swo/reg$hooks,eval(floor(reg$hbf/2)*2),mean,na.rm=T)/mean(reg$swo/reg$hooks)
  plot(as.numeric(names(a1)),a1,ylim=c(0,4),type="l",ylab = "Mean sets per trip")
  lines(as.numeric(names(a2)),a2,col=2)
  lines(as.numeric(names(a3)),a3,col=3)
  }
savePlot("CPUE_by_hooks_OS",type="png")
  
reg <- a[a$reg==3 & a$newfishingcat==2,]
tapply(floor(reg$hooks/200)*200,floor(reg$hooks/200)*200,length)
  
  

windows(width=20,height=18); par(mfrow=c(3,2),mar=c(2,2,1,1))
for(rg in 1:6) {
  reg <- op[op$reg==rg,]
#  hist(reg$alb,main=paste("Region",rg))
  a <- table(reg$alb)
  plot(as.numeric(names(a)),log(a))
  }
for(rg in 1:6) {
  reg <- op[op$reg==rg,]
  print(quantile(reg$alb,c(.001,.01,.05,.5,.95,.99,.999)))
  }

# remove sets without many hooks
a <- op[op$hooks > 1000,]
length(unique(op$vessid))
length(unique(a$vessid))
barplot(table(a$newfishingcat)/table(op$newfishingcat),xlab="Offshore / DWFN",ylab="Proportion of sets remaining",ylim=c(0,1))
savePlot("Effect of 1000 hooks limit on sets",type="png")
barplot(tapply(a$hooks,a$newfishingcat,sum)/tapply(op$hooks,op$newfishingcat,sum),xlab="Offshore / DWFN",ylab="Proportion of effort remaining",ylim=c(0,1))
savePlot("Effect of 1000 hook limit on hooks",type="png")

op2 <- op[op$hooks > 1000,]
save(op2,file="op2.RData")    #
#load(file="op2.RData")    #

#### Set up variables for trip-level clustering
tpall <- aggregate(cbind(hooks,hbf,alb,yft,bet,lat,lon,swo) ~ tripid + newfishingcat + vessid + trip_yr + reg,data=op2,mean)
tpall$bet_cpue <- 1000*tpall$bet/tpall$hooks
tpall$yft_cpue <- 1000*tpall$yft/tpall$hooks
tpall$alb_cpue <- 1000*tpall$alb/tpall$hooks
tpall$swo_cpue <- 1000*tpall$swo/tpall$hooks
tp <- aggregate(cbind(bet,yft,alb,swo,hbf) ~ tripid + newfishingcat + vessid + trip_yr + reg,data=op2,mean)
a <- aggregate(bet ~ tripid + newfishingcat + vessid + trip_yr + reg,data=op2,length)
dima <- dim(tp)[2]
tp <- cbind(tp,a[,6])
names(tp)[dima+1] <- "nset"
tp$tot <- tp$alb + tp$yft + tp$bet + tp$swo
tp$tuna <- tp$alb + tp$yft + tp$bet
tp$troptuna <- tp$yft + tp$bet
tp$bet_perc <- tp$bet/tp$tuna
tp$alb_perc <- tp$alb/tp$tuna
tp$yft_perc <- tp$yft/tp$tuna
tp$bet_totpc <- tp$bet/tp$tot
tp$alb_totpc <- tp$alb/tp$tot
tp$yft_totpc <- tp$yft/tp$tot
tp$swo_totpc <- tp$swo/tp$tot
tp$bet_troppc <- tp$bet/tp$troptuna
tp$yft_troppc <- tp$yft/tp$troptuna
gc()

str(tp)
for(r in 1:6) {
  windows(height=20,width=30); par(mfrow=c(2,2))
  hist(tp[tp$reg==r,15:17])
  title(paste("Region",r),outer=T,line=-1)
  savePlot(paste("hists_by_sp_reg_",rg,sep=""),type="png")
  }
  
# Regression trees to see which factors affect alb cpue
windows(width=20,height=18); par(mfrow=c(3,2),mar=c(2,2,2,2))
for(rg in 1:6) {
  a <- tpall[tpall$reg==rg,]
  tree1 <- rpart(alb_cpue ~ lat + lon + trip_yr + newfishingcat + bet_cpue + yft_cpue + swo_cpue + hooks, data=a)
  plot(tree1,main=paste("Region",rg),margin=0.1)
  text(tree1)
  }
savePlot("Rtrees_by_reg_by_trip",type="png")

windows(width=20,height=18); par(mfrow=c(3,2),mar=c(2,2,2,2))
for(rg in 1:4) {
  a <- op2[op2$reg==rg,]
  a$alb_cpue <- a$alb/a$hooks
  a$yft_cpue <- a$yft/a$hooks
  a$bet_cpue <- a$bet/a$hooks
  tree1 <- rpart(alb_cpue ~ lat + lon + trip_yr + newfishingcat + bet_cpue + yft_cpue + hooks, data=a)
  plot(tree1,main=paste("Region",rg),margin=0.1)
  text(tree1)
  }
savePlot("Rtrees_by_reg_by_set",type="png")

# Regression trees to see which factors affect alb percent
windows(width=22,height=18); par(mfrow=c(3,2),mar=c(2,2,2,2))
for(rg in 1:4) {
  a <- tpall[tpall$reg==rg,]
  a$alb_perc <- with(a,alb / (alb + yft + bet + oth))
  tree1 <- rpart(alb_perc ~ lat + lon + trip_yr + newfishingcat + hooks, data=a)
  plot(tree1,main=paste("Region",rg),margin=0.1)
  text(tree1)
  }
savePlot("Rtrees_albperc_by_reg_by_trip",type="png")

windows(width=22,height=18); par(mfrow=c(3,2),mar=c(2,2,2,2))
for(rg in 1:4) {
  a <- op2[op2$reg==rg,]
  a$alb_perc <- with(a,alb / (alb + yft + bet))
  tree1 <- rpart(alb_cpue ~ lat + lon + trip_yr + newfishingcat + hooks, data=a)
  plot(tree1,main=paste("Region",rg),margin=0.1)
  text(tree1)
  }
savePlot("Rtrees_albperc_by_reg_by_set",type="png")

# Look at numbers of clusters
tpx <- tp[,c("bet_totpc","yft_totpc","alb_totpc","swo_totpc","trip_yr","tripid","newfishingcat","hbf","reg")]
triptots <- tp[,c("bet","yft","alb","swo","trip_yr","tripid", "hbf","reg")]
tpx <- na.omit(tpx)
for(f in 1:2) {
  windows(width=20,height=18); par(mfrow=c(3,2),mar=c(2,4,4,2))
  nsp=4
  for(rg in 1:5) {
    a <- tpx[tpx$reg==rg & tpx$newfishingcat==f,]
    a[,1:nsp] <- scale(a[,1:nsp])
    wss <- (nrow(a)-1)*sum(apply(a[,1:nsp],2,var))
    for (i in 2:15) wss[i] <- sum(kmeans(a[,1:nsp],centers=i,iter.max = 40)$withinss)
    plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares",main=paste("Region",rg))
    }
  title(switch(f,"Offshore","Distant water"),outer=T,line=-1)
  savePlot(paste("cluster_num_v_SS_by_reg_",nsp,"_spp_fc",f,sep=""),type="png")
  }

# identify catch rates by trip
for(rg in 1:6) {
  a <- tp[tp$reg==rg & tp$newfishingcat==1,]
  windows(width=20,height=18); par(mfcol=c(3,2))
  boxplot(a$bet_perc~a$trip_yr, xlab="Year", ylab="bet_perc",main="Offshore")
  boxplot(a$alb_perc~a$trip_yr, xlab="Year", ylab="alb_perc")
  boxplot(a$yft_perc~a$trip_yr, xlab="Year", ylab="yft_perc")
  title(paste("Region",rg),outer=T,line=-1)
  a <- tp[tp$reg==rg & tp$newfishingcat==2,]
  boxplot(a$bet_perc~a$trip_yr, xlab="Year", ylab="bet_perc",main="Distant water")
  boxplot(a$alb_perc~a$trip_yr, xlab="Year", ylab="alb_perc")
  boxplot(a$yft_perc~a$trip_yr, xlab="Year", ylab="yft_perc")
  title(paste("Region",rg),outer=T,line=-1)
  savePlot(paste("boxplot_sp_percent_reg_",rg,sep=""),type="png")
}

## GMP: start here
table(tp$tuna==0)
tp <- tp[tp$tuna!=0, ]   # removing 82 trips


###################################
## GMP note - choose which of these periods you want to look at - incl post 1999, 1990-1999, or the whole kaboosh (not done here)
## Start with post 1999
## SDH change later to make it pre-1984, 1984-2002, and 2003-2011
## SDH change again to do all at once
###################################
#junk5 <- junk3[junk3$ret_year>=1999,]
#dim(junk5)         #[1] 469   10

for(f in 1:2) {
for(reg in 1:6) { # This takes a long time to run for R1 and R2
  gc()
  ax2 <- tp[tp$reg==reg & tp$newfishingcat==f,c("bet_perc","yft_perc","alb_perc","trip_yr","tripid", "hbf")]
  ax <- tp[tp$reg==reg & tp$newfishingcat==f & tp$nset >= 5,c("bet_perc","yft_perc","alb_perc","trip_yr","tripid", "hbf")]
  triptots2 <- tp[tp$reg==reg & tp$newfishingcat==f,c("bet","yft","alb","trip_yr","tripid", "hbf")]
  triptots <- tp[tp$reg==reg & tp$newfishingcat==f & tp$nset >= 5,c("bet","yft","alb","trip_yr","tripid", "hbf")]
  a <- na.omit(ax)
  a2 <- na.omit(ax2)
  a2[,1:3] <- scale(a2[,1:3])
  triptots[,1:3] <- scale(triptots[,1:3])
  triptots2[,1:3] <- scale(triptots2[,1:3])
  summary(triptots[,1:3])
  
  ##########################
  ## GMP: first dendrogram to id k cluster number
  ## These dendrograms take forever and use a lot of memory
  ##########################
  # other cluster method of all trips in entire series there are x total bet trips
  d <- dist(a[1:3], method = "euclidean") # distance matrix
  fit <- hclust(d, method="ward")
  windows(10,10);par(mfrow=c(2,2))
  plot(fit, labels = FALSE, hang=-1,  main = "dendrogram (whole period)") # display dendogram  #looks like 3 (or 4) main branches
  groups <- cutree(fit, k=3) # cut tree into 3 clusters
  # draw dendogram with red borders around the 3 clusters
  rect.hclust(fit, k=3, border="red")
  table(groups)
  clarax <- clara(a[1:3],3)             #clustering based upon the percent of spp in total catch of tuna     NOTE: using clusters from hclust k above(?? think it's actually independent clustering SDH)
  save(a,d,fit,clarax,file=paste("reg_",reg,"_fishcat_",f,"_fitfiles.RData",sep=""))
  plot(clarax)
  savePlot(file=paste("cluster_allyrs_dendrogram_R",reg,"_fishcat_",f,".png",sep=""),type="png")
} }

for(f in 1:2) {
for(reg in 1:6) { 
  load(file=paste("reg_",reg,"_fishcat_",f,"_fitfiles.RData",sep=""))
  print(paste(f,reg))
#  print(aggregate(a[,1:3],by=list(clarax$clustering),FUN=mean))
  print(aggregate(clarax$data,by=list(clarax$clustering),FUN=mean))
  print(table(clarax$clustering))

  junk4b <- table(clarax$clustering)
  junk4c <- cbind(a, clarax$clustering)                                                                    #appends cluster to dataframe, then names (in row below)
  names(junk4c) <- c("bet_perc","yft_perc","alb_perc","ret_year","tripid", "hbf", "cluster")
  assign(paste("clust_R",reg,"f",f,sep=""),junk4c)
  }
  }
#[1] "1 1"
#  Group.1  bet_perc   yft_perc  alb_perc
#1       1 0.1114975 0.07122494 0.8172775
#2       2 0.6928946 0.20268766 0.1044177
#3       3 0.4229064 0.08354069 0.4935530
#
#   1    2    3 
#5764 3394 4296 
#[1] "1 2"
#  Group.1  bet_perc   yft_perc  alb_perc
#1       1 0.1535551 0.03630696 0.8101379
#2       2 0.3768418 0.07870511 0.5444531
#3       3 0.7275501 0.13295710 0.1394928
#
#   1    2    3 
#1506  951  746 
#[1] "1 3"
#  Group.1   bet_perc  yft_perc   alb_perc
#1       1 0.22148606 0.7611870 0.01732689
#2       2 0.54531966 0.4328258 0.02185450
#3       3 0.07602534 0.1801956 0.74377909
#
#   1    2    3 
#7975 4652 1563 
#[1] "1 4"
#  Group.1  bet_perc  yft_perc   alb_perc
#1       1 0.2546610 0.7104589 0.03488004
#2       2 0.4506096 0.4465342 0.10285617
#3       3 0.7046596 0.2320966 0.06324385
#
#   1    2    3 
# 889 1332 1038 
#[1] "1 5"
#  Group.1   bet_perc  yft_perc   alb_perc
#1       1 0.14801464 0.8373417 0.01464371
#2       2 0.13159492 0.6098960 0.25850905
#3       3 0.08569698 0.3766817 0.53762129
#
#  1   2   3 
#295 125  73 
#[1] "1 6"
#  Group.1   bet_perc    yft_perc  alb_perc
#1       1 0.42577576 0.563617042 0.0106072
#2       2 0.05655172 0.943448276 0.0000000
#3       3 0.34889868 0.003524229 0.6475771
#
#1 2 3 
#2 1 1 
#[1] "2 1"
#  Group.1  bet_perc   yft_perc  alb_perc
#1       1 0.1373024 0.04386102 0.8188365
#2       2 0.6758974 0.15228243 0.1718202
#3       3 0.3686458 0.07383632 0.5575178
#
#  1   2   3 
#316 391 444 
#[1] "2 2"
#  Group.1  bet_perc   yft_perc  alb_perc
#1       1 0.3128962 0.07532644 0.6117774
#2       2 0.7196390 0.09911243 0.1812486
#3       3 0.1333888 0.04208825 0.8245229
#
#  1   2   3 
#875 436 941 
#[1] "2 3"
#  Group.1  bet_perc  yft_perc   alb_perc
#1       1 0.6082017 0.3506175 0.04118075
#2       2 0.2869452 0.5943994 0.11865539
#3       3 0.1297898 0.8412907 0.02891950
#
#  1   2   3 
#361 808 684 
#[1] "2 4"
#  Group.1  bet_perc  yft_perc   alb_perc
#1       1 0.3739027 0.5672119 0.05888539
#2       2 0.8370212 0.1335001 0.02947878
#3       3 0.6378820 0.3034594 0.05865863
#
#   1    2    3 
#1564 1701 1493 
#[1] "2 5"
#  Group.1   bet_perc  yft_perc  alb_perc
#1       1 0.11318090 0.6895983 0.1972208
#2       2 0.08808055 0.1277118 0.7842077
#3       3 0.09039036 0.3753404 0.5342693
#
#   1    2    3 
# 746  877 1013 
#[1] "2 6"
#  Group.1  bet_perc   yft_perc  alb_perc
#1       1 0.5086544 0.08812804 0.4032176
#2       2 0.1243172 0.41181090 0.4638719
#3       3 0.1774858 0.05099151 0.7715226
#
#  1   2   3 
#214  99 300 
clust_R1f1$reg <- 1;clust_R2f1$reg <- 2;clust_R3f1$reg <- 3;clust_R4f1$reg <- 4;clust_R5f1$reg <- 5;clust_R6f1$reg <- 6
clust_R1f2$reg <- 1;clust_R2f2$reg <- 2;clust_R3f2$reg <- 3;clust_R4f2$reg <- 4;clust_R5f2$reg <- 5;clust_R6f2$reg <- 6
clust_R1f1$fc <- 1;clust_R2f1$fc <- 1;clust_R3f1$fc <- 1;clust_R4f1$fc <- 1;clust_R5f1$fc <- 1;clust_R6f1$fc <- 1
clust_R1f2$fc <- 2;clust_R2f2$fc <- 2;clust_R3f2$fc <- 2;clust_R4f2$fc <- 2;clust_R5f2$fc <- 2;clust_R6f2$fc <- 2
hclusts <- rbind(clust_R1f1,clust_R2f1,clust_R3f1,clust_R4f1,clust_R5f1,clust_R6f1,
                 clust_R1f2,clust_R2f2,clust_R3f2,clust_R4f2,clust_R5f2,clust_R6f2)
save(hclusts,file="hclusts.RData")
load(file="hclusts.RData")

aggregate(hclusts[,1:3],by=list(hclusts$cluster,hclusts$reg,hclusts$fc),FUN=mean)
  
# Following up with a few alternative cluster methods - I think these were done just for comparison
# Model Based Clustering
library(mclust)
fit <- Mclust(triptots[,1:3])
windows(15,15);par(mfrow=c(2,2))
plot(fit, triptots[,1:3]) # plot results
print(fit) # display the best model
plot.Mclust

# Ward Hierarchical Clustering with Bootstrapped p values
library(pvclust)
iv <- as.matrix((a[,1:3]))
summary(iv)
result <- pvclust((lung), method.dist="cor", method.hclust="average", nboot=200)
plot(result) # dendogram with p values
pvrect(result, alpha=.95)

fit <- pvclust(t(iv), method.dist="cor", method.hclust="average", use.cor="all.obs",nboot=200)
fit <- pvclust(t(iv), method.hclust="ward", method.dist="euclidean",nboot=500)
plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
  pvrect(fit, alpha=.95)


##################################
tpx <- tp[,c("bet_perc","yft_perc","alb_perc","trip_yr","tripid", "hbf","reg","newfishingcat")]
tpx <- na.omit(tpx)
nsp=3
clx <- claracl <- list()
ncl <- c(3,3,3,3,3,3)
tpx$clust <- tpx$claracl <- NA
for (f in 1:2) {
  for(rg in 1:4) {
    a <- tpx[tpx$reg==rg & tpx$newfishingcat==f,]
    a[,1:nsp] <- scale(a[,1:nsp])
    clx[[rg + 6*(f-1)]] <- kmeans(a[,1:nsp],centers=ncl[rg],iter.max = 40)
    claracl[[rg + 6*(f-1)]] <- clara(a[,1:nsp],k=ncl[rg])
    tpx[tpx$reg==rg & tpx$newfishingcat==f,]$clust <- clx[[rg + 6*(f-1)]]$cluster
    tpx[tpx$reg==rg & tpx$newfishingcat==f,]$claracl <- claracl[[rg + 6*(f-1)]]$clustering
    }
  }
  
# Apply clusters to sets
op2$clust <- as.factor(tpx[match(paste(op2$tripid,op2$reg),paste(tpx$tripid,tpx$reg)),]$clust)
op2$claracl <- as.factor(tpx[match(paste(op2$tripid,op2$reg),paste(tpx$tripid,tpx$reg)),]$claracl)
op2$hclustcl <- as.factor(hclusts[match(paste(op2$tripid,op2$reg),paste(hclusts$tripid,hclusts$reg)),]$cluster)
save(op2,file="pl2_clustered.RData")
load(file="pl2_clustered.RData")

length(unique(op2$vessid[!is.na(op2$hclustcl)]))
length(unique(op2$trip_id[!is.na(op2$hclustcl)]))
length(unique(op2$set_id[!is.na(op2$hclustcl)]))

table(op2$hclustcl,useNA="always")
a <- op2[is.na(op2$hclustcl)==F,]
with(a,aggregate(cbind(alb_perc,bet_perc,yft_perc) ~ hclustcl + reg,FUN=mean))
with(a,aggregate(cbind(alb_perc,bet_perc,yft_perc) ~ clust + reg,FUN=mean))
with(a,aggregate(cbind(alb_perc,bet_perc,yft_perc) ~ claracl + reg,FUN=mean))


#windows(width=20,height=18); par(mfrow=c(3,2))
# Plot clusters on a map. Mostly SDH code from here down. 
require(maps); require(maptools); require(mapdata)
windows(width=20,height=18); par(mfrow=c(3,2))
for(rg in 1:4) {
  a <- op2[op2$reg==rg & op2$op_yr < 1980,]
  plot(1:10,1:10,type="n",xlab="Longitude",ylab="Latitude",xlim=range(a$lon),ylim=range(a$lat))
  for(cl in 1:ncl[rg]) {
    b <- a[a$clust==cl,]
    points(b$lon,b$lat,cex=0.5,col=cl)
    }
  lines(ezpac$lon, ezpac$lat)
  map(add=T)
  }
savePlot("map_clusters_kmeans_pre1980",type="png")
windows(width=20,height=18); par(mfrow=c(3,2))
for(rg in 1:4) {
  a <- op2[op2$reg==rg & op2$op_yr < 1980,]
  plot(1:10,1:10,type="n",xlab="Longitude",ylab="Latitude",xlim=range(a$lon),ylim=range(a$lat))
  for(cl in 1:ncl[rg]) {
    b <- a[a$claracl==cl,]
    points(b$lon,b$lat,cex=0.5,col=cl)
    }
  lines(ezpac$lon, ezpac$lat)
  map(add=T)
  }
savePlot("map_clusters_clara_pre1980",type="png")
windows(width=20,height=18); par(mfrow=c(3,2))
for(rg in 1:4) {
  a <- op2[op2$reg==rg & op2$op_yr < 1980,]
  plot(1:10,1:10,type="n",xlab="Longitude",ylab="Latitude",xlim=range(a$lon),ylim=range(a$lat))
  for(cl in 1:ncl[rg]) {
    b <- a[a$hclustcl==cl,]
    points(b$lon,b$lat,cex=0.5,col=cl)
    }
  lines(ezpac$lon, ezpac$lat)
  map(add=T)
  }
savePlot("map_clusters_hclust_pre1980",type="png")

# Make 2x2 maps
postyr <- 1980
preyr <- 1980
agglocs_hc <- with(op2[op2$op_yr <= preyr,],tapply(hooks,list(2*floor(lon/2),2*floor(lat/2),hclustcl,reg),sum,na.rm=T))
agg <- with(op2[op2$op_yr <= preyr,],tapply(hooks,list(2*floor(lon/2),2*floor(lat/2),reg),sum,na.rm=T))
#agglocs_hc <- tapply(op2$hooks,list(1*floor(op2$lon/1),1*floor(op2$lat/1),op2$hclustcl,op2$reg),sum,na.rm=T)
#agg <- tapply(op2$hooks,list(1*floor(op2$lon/1),2*floor(op2$lat/1),op2$reg),sum,na.rm=T)
for(cl in 1:3) {
  windows(12,12)
  plot(1:10,1:10,type="n",xlim=c(140,250),ylim=c(-45,0),xlab="Longitude",ylab="Latitude",main=paste("Cluster",cl))
  for(reg in 1:4) {
    a <- agglocs_hc[,,cl,reg]/agg[,,reg]
    image(as.numeric(dimnames(a)[[1]]),as.numeric(dimnames(a)[[2]]),a,add=T)
    contour(as.numeric(dimnames(a)[[1]]),as.numeric(dimnames(a)[[2]]),a,add=T)
   }
  map(add=T)
  lines(c(180,180),c(-65,10),lwd=2)
  lines(c(100,300),c(-25,-25),lwd=2,col=1)
  savePlot(paste("Cluster_map_hcl_2x2_pre",preyr,cl,sep="_"),type="png")
  }

agglocs_hc[1:40,16:30,1,1]
props[1:40,16:30,1,1]
props[1:40,16:30,2,1]
props[1:40,16:30,3,1]


table(tpx$reg,tpx$clust)
table(tpx$reg,tpx$claracl)
aggregate(cbind(alb_perc,bet_perc,yft_perc) ~ clust + reg,data=op2,mean)
aggregate(cbind(alb_perc,bet_perc,yft_perc) ~ claracl + reg,data=op2,mean)
aggregate(cbind(claracl==1,claracl==2,claracl==3) ~ reg,data=op2,mean)
#   claracl reg  alb_perc   bet_perc   yft_perc
#1        1   1 0.9190255 0.01937280 0.06160166   *  40%
#2        2   1 0.7499081 0.05723085 0.19286104      36%
#3        3   1 0.4138501 0.15995933 0.42619061      24%
#4        1   2 0.7125166 0.13201548 0.15546788      37%
#5        2   2 0.9056770 0.04079331 0.05352972   *  35%
#6        3   2 0.3106432 0.32911140 0.36024536      28%
#7        1   3 0.8905033 0.05205180 0.05744493      17%
#8        2   3 0.9727803 0.01270842 0.01451126   *  82%
#9        3   3 0.4460348 0.20826673 0.34569844      0.5%
#10       1   4 0.9010575 0.04170756 0.05723498      26%
#11       2   4 0.9745593 0.01273401 0.01270669   *  73%
#12       3   4 0.5805059 0.17887636 0.24061776       1%

tapply(cbind(tp$bet_perc),list(tpx$reg,tpx$clust),mean)
table(op2$reg,op2$clust)
table(op2$reg,op2$claracl)
  
windows(width=15,height=18); par(mfrow=c(3,2))
for(r in 1:4) { 
  plot(1:10,1:10,xlim=c(1960,2011),ylim=c(0,1),xlab="Year",ylab="",main="kmeans")
  a <- with(op2[op2$reg==r,],tapply(op_yr,list(op_yr,clust),length))
  a[is.na(a)] <- 0 
  for(cl in 1:3) {
    points(as.numeric(rownames(a)),a[,cl]/apply(a,1,sum),col=cl,pch=cl)
    }
  }
legend(1990,0.5,legend=paste("cluster",1:3),pch=1:3,col=1:3)
savePlot("clusters_by_yr_kmeans",type="png")

windows(width=15,height=18); par(mfrow=c(3,2))
for(r in 1:4) { 
  plot(1:10,1:10,xlim=c(1960,2011),ylim=c(0,1),xlab="Year",ylab="",main="clara")
  a <- with(op2[op2$reg==r,],tapply(op_yr,list(op_yr,claracl),length))
  a[is.na(a)] <- 0 
  for(cl in 1:3) {
    points(as.numeric(rownames(a)),a[,cl]/apply(a,1,sum),col=cl,pch=cl)
    }
  }
legend(1990,0.5,legend=paste("cluster",1:3),pch=1:3,col=1:3)
savePlot("clusters_by_yr_clara",type="png")

windows(width=15,height=18); par(mfrow=c(3,2))
for(r in 1:4) { 
  plot(1:10,1:10,xlim=c(1960,2011),ylim=c(0,1),xlab="Year",ylab="",main="hclust")
  a <- with(op2[op2$reg==r,],tapply(op_yr,list(op_yr,hclustcl),length))
  a[is.na(a)] <- 0 
  for(cl in 1:3) {
    points(as.numeric(rownames(a)),a[,cl]/apply(a,1,sum),col=cl,pch=cl)
    }
  }
legend(1990,0.5,legend=paste("cluster",1:3),pch=1:3,col=1:3)
savePlot("clusters_by_yr_hclust",type="png")

usecl_h <- c(paste(1,c(1,2)),paste(2,c(1,2)),paste(3,c(1,2)),paste(4,c(1,2)))
usecl_c <- c(paste(1,c(1,2)),paste(2,c(1,2)),paste(3,c(1,2)),paste(4,c(1,2)))
usecl_k <- c(paste(1,c(2)),paste(2,c(2)),paste(3,c(1,2,3)),paste(4,c(1,3)))
op2$vessid <- as.factor(op2$vessid)
albpl <- op2[paste(op2$reg,op2$claracl) %in% usecl_c,]
save(albpl,file="albacore_doubleclusters_clara.RData")
albpl <- op2[paste(op2$reg,op2$hclustcl) %in% usecl_h,]
save(albpl,file="albacore_doubleclusters_hclust.RData")
albpl <- op2[paste(op2$reg,op2$hclustcl) %in% usecl_k,]
save(albpl,file="albacore_doubleclusters_kmeans.RData")

ls()
rm(a,a1,a2,ax,b,b2,channel,cl,claracl,clx,d,dima,f,fit,groups,i,l_bycat,ncl,nsp,op,pl5,r,reg,reg1,rg,tp,tp1,tpall,tpx,triptots,usecl,wss)

table(albop$reg,albop$newfishingcat)

####################################################
# CPUE analysis
####################################################

load(file="albacore_doubleclusters_hclust.RData")
albpl_h <- albpl
load(file="albacore_doubleclusters_clara.RData")
albpl_c <- albpl
load(file="albacore_doubleclusters_kmeans.RData")
albpl_k <- albpl
albpl_h$hclustcl <- as.factor(albpl_h$hclustcl)
albpl_h$clust <- as.factor(albpl_h$clust)
albpl_h$claracl <- as.factor(albpl_h$claracl)
albpl_c$hclustcl <- as.factor(albpl_c$hclustcl)
albpl_c$clust <- as.factor(albpl_c$clust)
albpl_c$claracl <- as.factor(albpl_c$claracl)
albpl_k$hclustcl <- as.factor(albpl_k$hclustcl)
albpl_k$clust <- as.factor(albpl_k$clust)
albpl_k$claracl <- as.factor(albpl_k$claracl)

albpl_c$cl <- albpl_c$claracl
albpl_h$cl <- albpl_h$hclustcl
albpl_k$cl <- albpl_k$clust

windows(15,15); par(mfrow=c(3,2))
for(r in 1:4) {
  dd <- albpl_h[albpl_h$reg==r,]
  catch  <- tapply(dd$alb, dd$set_yrqtr, sum)
  effort <- tapply(dd$hooks, dd$set_yrqtr, sum)
  cpue <- 1000*catch/effort
  ytot <- sort(unique(dd$set_yrqtr))
  timeseries <- seq(1960.125, 2012, by=0.25)
  timeNA <- data.frame(ytot[match(timeseries,ytot)])
  names(timeNA) <- c("time")
  cpue <- data.frame(cpue)
  cpue$time <- row.names(cpue)
  final <- data.frame(cpue[match(timeNA[,1],cpue[,2]),1])
  names(final) <- c("cpue")
  plot(timeNA$time, final$cpue, type="n", col="red", xlim=c(1960,2012),ylim=c(0,max(cpue$cpue)))
  lines(timeNA$time,final$cpue, col="blue")
  
  dd <- albpl_c[albpl_c$reg==r,]
  catch  <- tapply(dd$alb, dd$set_yrqtr, sum)
  effort <- tapply(dd$hooks, dd$set_yrqtr, sum)
  cpue <- 1000*catch/effort
  ytot <- sort(unique(dd$set_yrqtr))
  timeseries <- seq(1960.125, 2012, by=0.25)
  timeNA <- data.frame(ytot[match(timeseries,ytot)])
  names(timeNA) <- c("time")
  cpue <- data.frame(cpue)
  cpue$time <- row.names(cpue)
  final <- data.frame(cpue[match(timeNA[,1],cpue[,2]),1])
  names(final) <- c("cpue")
  lines(timeNA$time,final$cpue, col="red")
  
  dd <- albpl_k[albpl_k$reg==r,]
  catch  <- tapply(dd$alb, dd$set_yrqtr, sum)
  effort <- tapply(dd$hooks, dd$set_yrqtr, sum)
  cpue <- 1000*catch/effort
  ytot <- sort(unique(dd$set_yrqtr))
  timeseries <- seq(1960.125, 2012, by=0.25)
  timeNA <- data.frame(ytot[match(timeseries,ytot)])
  names(timeNA) <- c("time")
  cpue <- data.frame(cpue)
  cpue$time <- row.names(cpue)
  final <- data.frame(cpue[match(timeNA[,1],cpue[,2]),1])
  names(final) <- c("cpue")
  lines(timeNA$time,final$cpue, col="green")
  }
savePlot("CPUE_raw2_by_reg",type="png")
  
    # filter for active boats and post-2000

lenun <- function(x) length(unique(x))

for(r in 2:4) {
  print(r)
  dd <- albpl_c[albpl_c$reg==r,]
#  vessyrs <- c(5,5,1,1)[r]
  vessyrs <- c(5,10,1,1)[r]
  a <- tapply(dd$vessid,list(dd$vessid,dd$set_yrqtr),length)
  a1 <- apply(is.na(a[,])=="FALSE",1,sum)
  a1 <- a1[a1>vessyrs]
  rb <- tapply(dd$set_yrqtr,dd$vessid,max)
  rb1 <- rb[!is.na(rb) & rb>2000]
#  activeboats <- c(as.numeric(names(a1)), as.numeric(names(rb1)))
  activeboats <- as.numeric(names(a1)) # plenty of recent activity so no need for recent boats
  activeboats <- unique(activeboats)
  a <- dd[dd$vessid %in% activeboats,]

  modtxt <- "alb ~ as.factor(set_yrqtr) + offset(log(hooks))"
  if(length(unique(a$cl)) > 1) modtxt <- "alb ~ as.factor(set_yrqtr) + cl + offset(log(hooks))"
  glm_nb1 <- glm.nb(modtxt, data=a)
  save(glm_nb1,file=paste("clara_glm_nb1_",r,".RData",sep="")); rm(glm_nb1); gc()
  modtxt <- "alb ~ as.factor(set_yrqtr) + latlong + as.factor(vessid)+offset(log(hooks))"
  if(length(unique(a$cl)) > 1) modtxt <- "alb ~ as.factor(set_yrqtr) + latlong + cl + as.factor(vessid)+offset(log(hooks))"
  glm_nb2 <- glm.nb(modtxt, data=a)
  save(glm_nb2,file=paste("clara_glm_nb2_",r,".RData",sep="")); 
  
  ymax <- 3
  yr <- sort(unique(a$set_yrqtr))
  y <- length(yr)
  load(file=paste("clara_glm_nb1_",r,".RData",sep=""))
  coefs1 <- exp(c(0,summary(glm_nb1)$coefficients[2:y,1]))
  coefs1 <- coefs1/mean(coefs1)
  coefs2 <- exp(c(0,summary(glm_nb2)$coefficients[2:y,1]))
  coefs2 <- coefs2/mean(coefs2)

  catch  <- tapply(a$alb, a$set_yrqtr, sum)
  effort <- tapply(a$hooks, a$set_yrqtr, sum)
  cpue <- 1000*catch/effort
  cpue <- cpue / mean(cpue)
  timeseries <- seq(1960.125, 2012, by=0.25)
  timeNA <- data.frame(yr[match(timeseries,yr)])
  names(timeNA) <- c("time")

  coefs3 <- cbind(yr,cpue, coefs1, coefs2)
  final <- data.frame(coefs3[match(timeNA[,1],coefs3[,1]),2])
  final2 <- data.frame(coefs3[match(timeNA[,1],coefs3[,1]),3])
  final3 <- data.frame(coefs3[match(timeNA[,1],coefs3[,1]),4])
  names(final) <- c("cpue")
  names(final2) <- c("coefs1")
  names(final3) <- c("coefs2")
  final4 <- cbind(final, final2, final3)

  plot(timeNA$time, final4$coefs1, type="n", xlab=c("Time"), ylab=c("Relative abundance"), xlim=c(1960,2012),ylim=c(0,ymax))
  lines(timeNA$time,final4$cpue, col="black")
  lines(timeNA$time,final4$coefs1, col="blue")
  lines(timeNA$time,final4$coefs2, col="red")
  fname <- paste("spalb_clara", r,".png",sep="")
  savePlot(fname,type="png")

  coefs3 <- summary(glm_nb2)$coefficients[1:y,2]
  junk <- data.frame(cbind(yr, coefs2, coefs3))
  names(junk) <- c("yrqtr", "coef", "std.error")
  
  newdat <- with(a,expand.grid(set_yrqtr=sort(unique(a$set_yrqtr)),latlong=latlong[1],cl=cl[1],vessid=vessid[1],hooks=median(hooks)))
  a2 <- predict(glm_nb2,newdata=newdat,type="terms",se.fit=T)
  junk$coef2 <- exp(a2$fit[,1]) / mean(exp(a2$fit[,1]))
  junk$std_error2 <- a2$se.fit[,1]
  
  fname <- paste("spalb_clara", r,".dat",sep="")
  write.table(junk, fname, sep=c(" "), row.names=F)
  }
  

# consider activity filter to remove some vessels

lenun <- function(x) length(unique(x))
with(albpl,tapply(vessid,reg,lenun))

atot <- albpl[albop$reg==2,]

#for (vessyrs in c(10)) {
for (vessyrs in c(10)) { # for area 2
    # filter for active boats and post-2000
    a <- tapply(atot$vessid,list(atot$vessid,atot$set_yrqtr),length)
    a1 <- apply(is.na(a[,])=="FALSE",1,sum)
    a1 <- a1[a1>vessyrs]
    rb <- tapply(atot$set_yrqtr,atot$vessid,max)
    rb1 <- rb[rb>2000]
    activeboats <- c(as.numeric(names(a1)), as.numeric(names(rb1)))
    activeboats <- unique(activeboats)
    a <- atot[atot$vessid %in% activeboats,]
    rb <- tapply(a$set_yrqtr,a$vessid,max)
    rb1 <- rb[rb>2000]
    recentboats <- as.numeric(names(rb1))
    a$era <- 0
    a$era[a$vessid %in% recentboats] <- 1
}
with(a,tapply(vessid,reg,lenun))

dim(a)
unique(a$vessid) # reg 1 548
unique(a$vessid) # reg 2 1006 filter of 40 
unique(a$vessid) # reg 3 277 activity filter not used 
unique(a$vessid) # reg 4 368 
table(a$newfishingcat)



setwd("D:/L/alb/2012/assessment/Data_preparation/CPUE/myStd")
for(r in 1:4) {
  a <- read.table(paste("cpue_reg",r,".dat",sep=""),header=T)
  a <- a[!is.na(a[,4]),c(1,4,6)]
  names(a) <- c("yrqtr","coef","std.error")
  write.table(a,file=paste("spalb_",r,".dat",sep=""),row.names=F)
}

######################################
# compare latlong at a finer scale
minqtrs=8; maxqtrs=200; runreg <- 9;addmain<-F; addbranch<-F;addother=F;addalb=F;fc="both";runsp="bet";fc="OS";mt="deltabin"
glmdat5 <- select_data(opx,runreg=9,minqtrs=minqtrs,runsp=runsp,mt="deltapos")
glmdat.bin5 <- select_data(opx,runreg=9,minqtrs=minqtrs,runsp=runsp,mt="deltabin")
wtt.opbin.area5   <- mk_wts(glmdat.bin5,wttype="area")
fmla.boatbin <- make_formula(runsp,modtype="deltabin",addboat=T)
model.boatbin.area5   <- glm(fmla.boatbin,data=glmdat.bin5,weights=wtt.opbin.area5,family="binomial");gc()
abin <- length(unique(glmdat.bin5$yrqtr))
coefs.boatbin.area5 <- get.coefs(model.boatbin.area5,abin)

glmdat2 <- select_data(opx,runreg=9,minqtrs=minqtrs,runsp=runsp,mt="deltapos",llstrat=2)
glmdat.bin2 <- select_data(opx,runreg=9,minqtrs=minqtrs,runsp=runsp,mt="deltabin",llstrat=2)
wtt.opbin.area2   <- mk_wts(glmdat.bin2,wttype="area")
fmla.boatbin <- make_formula(runsp,modtype="deltabin",addboat=T)
model.boatbin.area2   <- glm(fmla.boatbin,data=glmdat.bin2,weights=wtt.opbin.area2,family="binomial");gc()
abin <- length(unique(glmdat.bin2$yrqtr))
coefs.boatbin.area2 <- get.coefs(model.boatbin.area2,abin)

wtt.oppos.area5   <- mk_wts(glmdat5,wttype="area")
fmla.boatpos <- make_formula(runsp,modtype="deltapos",addboat=T)
model.boatpos.area5   <- glm(fmla.boatpos,data=glmdat5,weights=wtt.oppos.area5,family="gaussian");gc()
apos <- length(unique(glmdat5$yrqtr))
coefs.boatpos.area5 <- get.coefs(model.boatpos.area5,apos)

wtt.oppos.area2   <- mk_wts(glmdat2,wttype="area")
fmla.boatpos <- make_formula(runsp,modtype="deltapos",addboat=T)
model.boatpos.area2   <- glm(fmla.boatpos,data=glmdat2,weights=wtt.oppos.area2,family="gaussian");gc()
apos <- length(unique(glmdat2$yrqtr))
coefs.boatpos.area2 <- get.coefs(model.boatpos.area2,apos)
opyr <- unique(glmdat2$yrqtr)

glmdat1 <- select_data(opx,runreg=9,minqtrs=minqtrs,runsp=runsp,mt="deltapos",llstrat=1)
wtt.oppos.area1   <- mk_wts(glmdat1,wttype="area")
fmla.boatpos <- make_formula(runsp,modtype="deltapos",addboat=T)
model.boatpos.area1   <- glm(fmla.boatpos,data=glmdat1,weights=wtt.oppos.area1,family="gaussian");gc()
apos <- length(unique(glmdat1$yrqtr))
coefs.boatpos.area1 <- get.coefs(model.boatpos.area1,apos)

plot.agg.slope.ratio(coefs.boatpos.area5,coefs.boatpos.area2,opyr,opyr, titl=paste("Region",runreg,fc,"oppos 5deg vs 2deg"),lab1="oppos 5deg",lab2="oppos 2deg",fname="oppos_5deg_vs_2deg")  # plot results
plot.agg.slope.ratio(coefs.boatpos.area5,coefs.boatpos.area1,opyr,opyr, titl=paste("Region",runreg,fc,"oppos 5deg vs 1deg"),lab1="oppos 5deg",lab2="oppos 1deg",fname="oppos_5deg_vs_1deg")  # plot results

windows()
windows();par(mfrow=c(3,2))
plotdiags(model.boatpos.area5$res,ti="5 degree squares")
plotdiags(model.boatpos.area2$res,ti="2 degree squares")
plotdiags(model.boatpos.area1$res,ti="1 degree squares")
savePlot("compare_latlongstrat_diagnostics",type="png")

####################################################################
# Movements between sets
# Set up data
opcore <- op2[op2$lat > -5 & op2$lat < 10 & op2$reg %in% 1:6 & !is.na(op2$tripid),]
table(opcore$op_yr)
a <- with(opcore,opcore[order(paste(vessid,op_yr,formatC(op_mon,width=2,flag=0),formatC(op_day,width=2,flag=0))),])
n <- dim(a)[1]
a$latmove <- a$lonmove <- NA
a[2:n,"latmove"] <- a[2:n,"lat"] - a[1:(n-1),"lat"]
a[2:n,"lonmove"] <- a[2:n,"lon"] - a[1:(n-1),"lon"]
a$sdate <- with(a,as.Date(paste(op_yr,formatC(op_mon,width=2,flag=0),formatC(op_day,width=2,flag=0),sep="-"),format="%Y-%m-%d"))
a$days <- NA
lstv <- c(NA,a$tripid[1:(n-1)])
a[a$tripid != lstv | is.na(lstv),]$days <- NA
a$days[2:n] <- a[2:n,"sdate"] - a[1:(n-1),"sdate"]
a[a$days > 10 | a$days < 1 | is.na(a$days),]$latmove <- NA
a[a$days > 10 | a$days < 1 | is.na(a$days),]$lonmove <- NA
a$gdays <- a$days
a[a$days > 10 | a$days < 1 | is.na(a$days),]$gdays <- NA
hist(a$gdays) 
hist(a$latmove)
hist(a$lonmove)
a$mv <- sqrt(a$latmove^2 + a$lonmove^2)
hist(a$mv,nclass=50)
table(a$mv)
a[a$mv>10 & !is.na(a$mv),]
dat <- a
dat$nxtmv <- NA
n <- dim(dat)[1]
dat$nxtmv[1:(n-1)] <- dat$mv[2:n]
# clean dataset
windows()
hist(dat$sets_per_trip,nclass=50)
dat <- dat[dat$sets_per_trip <= 150,]
dat <- dat[dat$sets_per_trip >= 20,]
dat[100:400,c("op_yr","op_mon","op_day","days","lat","lon","lonmove","latmove","mv","nxtmv")]

# mean distance between sets within a trip
dat$trop <- dat$bet + dat$yft
with(dat,plot(tapply(mv,op_yr,mean,na.rm=T)))
windows()
hist(dat$mv)
table(dat$mv>100)
x <- dat[dat$mv < 15 & !is.na(dat$mv),] 
a <- tapply(x$mv,x$op_yr,mean,na.rm=T)
plot(as.numeric(names(a)),a,ylim=c(0,1),type="l",xlab="Year",ylab="Movement")
savePlot("Movement between sets",type="png")
library(lattice)
histogram(~ mv | op_yr,data=x)
savePlot("movement_distrib_by_yr",type="png")
windows()
plot(tapply(dat$mv==0,dat$op_yr,mean,na.rm=T),type="l",ylim=c(0,1))
savePlot("Movements_prop_of_zeroes",type="png")
dat$yq <- as.factor(dat$yrqtr)
dat$ll <- as.factor(dat$latlong)
model <- glm((mv > 1.5) ~ as.factor(yrqtr) + latlong + bet + yft,data=dat,family="binomial")
vif(model) # this was okay yrqtr 5.79 and latlong 5.55 but bet and yft about 1.2 each
model_1.5 <- glm((mv > 1.5) ~ yq + ll + ns(bet,df=4) + ns(yft,df=4),data=dat,family="binomial")
model_1 <- glm((mv > 1) ~ yq + ll + ns(bet,df=4) + ns(yft,df=4),data=dat,family="binomial")
model_1trop <- glm((mv > 1) ~ yq + ll + ns(trop,df=4),data=dat,family="binomial")
summary(model_1)
#Anova(model)
allyq <- sort(unique(dat$yq))
newdat <- expand.grid(mv=F,yq=dat$yq[1],ll=dat$ll[10],bet=seq(0,quantile(dat$bet,0.95),1),yft=mean(dat$yft))
#a <- predict(model,newdata=newdat,type="response",se.fit=T)
plot(1:40,1:40,xlim=c(0,40),ylim=c(0,.3),xlab="Number of species caught",ylab="Probability of moving",type="n")
a <- predict(model_1,newdata=newdat,type="terms",se.fit=T)
res <- apply(a$fit,1,sum) +  attr(a$fit,"constant")
points(newdat$bet,inv.logit(res),col=1)
lines(newdat$bet,inv.logit(res+2*a$se.fit[,3]),col=1)
lines(newdat$bet,inv.logit(res-2*a$se.fit[,3]),col=1)
newdat <- expand.grid(mv=F,yq=dat$yq[1],ll=dat$ll[10],bet=mean(dat$bet),yft=seq(0,quantile(dat$yft,0.95),1))
a <- predict(model_1,newdata=newdat,type="terms",se.fit=T)
res <- apply(a$fit,1,sum) +  attr(a$fit,"constant")
points(newdat$yft,inv.logit(res),col=2)
lines(newdat$yft,inv.logit(res+2*a$se.fit[,4]),col=2)
lines(newdat$yft,inv.logit(res-2*a$se.fit[,4]),col=2)
legend("bottomright",legend=c("Bigeye","Yellowfin","Bigeye + Yellowfin"),col=c(1,2,3),pch=1)
newdat <- expand.grid(mv=F,yq=dat$yq[1],ll=dat$ll[10],trop=seq(0,quantile(dat$trop,0.95),1))
a <- predict(model_1trop,newdata=newdat,type="terms",se.fit=T)
res <- apply(a$fit,1,sum) +  attr(a$fit,"constant")
points(newdat$trop,inv.logit(res),col=3)
lines(newdat$trop,inv.logit(res+2*a$se.fit[,3]),col=3)
lines(newdat$trop,inv.logit(res-2*a$se.fit[,3]),col=3)
legend("bottomright",legend=c("Bigeye","Yellowfin","Bigeye + Yellowfin"),col=c(1,2,3),pch=1)
savePlot("Movements_GLM_BET_YFT_trop",type="png")

dat$decade <- 10 * floor(dat$op_yr / 10)
windows(10,10);par(mfrow=c(2,2))
for(dc in seq(1980,2010,10)) {
  dat2 <- dat[dat$decade==dc,]
  model_1 <- glm((mv > 1) ~ yq + ll + ns(bet,df=4) + ns(yft,df=4),data=dat2,family="binomial")
  model_1trop <- glm((mv > 1) ~ yq + ll + ns(trop,df=4),data=dat2,family="binomial")
  summary(model_1)
  #Anova(model)
  allyq <- sort(unique(dat2$yq))
  newdat <- expand.grid(mv=F,yq=dat2$yq[1],ll=dat2$ll[1],bet=seq(0,quantile(dat2$bet,0.95),1),yft=mean(dat2$yft))
  #a <- predict(model,newdata=newdat,type="response",se.fit=T)
  plot(1:40,1:40,xlim=c(0,40),ylim=c(0,.3),xlab="Number caught per set",ylab="Probability of moving > 1 cell",type="n",main=dc)
  a <- predict(model_1,newdata=newdat,type="terms",se.fit=T)
  res <- apply(a$fit,1,sum) +  attr(a$fit,"constant")
  points(newdat$bet,inv.logit(res),col=1)
  lines(newdat$bet,inv.logit(res+2*a$se.fit[,3]),col=1)
  lines(newdat$bet,inv.logit(res-2*a$se.fit[,3]),col=1)
  newdat <- expand.grid(mv=F,yq=dat2$yq[1],ll=dat2$ll[1],bet=mean(dat2$bet),yft=seq(0,quantile(dat2$yft,0.95),1))
  a <- predict(model_1,newdata=newdat,type="terms",se.fit=T)
  res <- apply(a$fit,1,sum) +  attr(a$fit,"constant")
  points(newdat$yft,inv.logit(res),col=2)
  lines(newdat$yft,inv.logit(res+2*a$se.fit[,4]),col=2)
  lines(newdat$yft,inv.logit(res-2*a$se.fit[,4]),col=2)
  legend("bottomright",legend=c("Bigeye","Yellowfin","Bigeye + Yellowfin"),col=c(1,2,3),pch=1)
  newdat <- expand.grid(mv=F,yq=dat2$yq[1],ll=dat2$ll[1],trop=seq(0,quantile(dat2$trop,0.95),1))
  a <- predict(model_1trop,newdata=newdat,type="terms",se.fit=T)
  res <- apply(a$fit,1,sum) +  attr(a$fit,"constant")
  points(newdat$trop,inv.logit(res),col=3)
  lines(newdat$trop,inv.logit(res+2*a$se.fit[,3]),col=3)
  lines(newdat$trop,inv.logit(res-2*a$se.fit[,3]),col=3)
  legend("bottomright",legend=c("Bigeye","Yellowfin","Bigeye + Yellowfin"),col=c(1,2,3),pch=1)
  }
savePlot("Movements_GLM_BET_YFT_trop_by_decade",type="png")

str(dat)
model_x <- glm((mv > 1) ~ ns(yrqtr,df=5) + ll + ns(yrqtr,df=5)*bet + ns(yrqtr,df=5)*yft,data=dat,family="binomial")
model_xtrop <- glm((mv > 1) ~ ns(yrqtr,df=5) + ll + ns(yrqtr,df=5)*trop,data=dat,family="binomial")
summary(model_x)
allyq <- sort(unique(dat$yq))
newdat <- expand.grid(mv=F,yrqtr=seq(1980.125,2011.875,.25),ll=dat$ll[1],bet=0,yft=mean(dat$yft))
#a <- predict(model,newdata=newdat,type="response",se.fit=T)
plot(1:40,1:40,xlim=c(1980,2011),ylim=c(0,.3),xlab="Year",ylab="Probability of moving after 0 catch",type="n")
a <- predict(model_x,newdata=newdat,type="terms",se.fit=T)
res <- apply(a$fit,1,sum) +  attr(a$fit,"constant")
lines(newdat$yrqtr,inv.logit(res),col=1)
lines(newdat$yrqtr,inv.logit(res+2*a$se.fit[,3]),col=1,lty=2)
lines(newdat$yrqtr,inv.logit(res-2*a$se.fit[,3]),col=1,lty=2)
newdat <- expand.grid(mv=F,yrqtr=seq(1980.125,2011.875,.25),ll=dat$ll[1],bet=mean(dat$bet),yft=0)
a <- predict(model_x,newdata=newdat,type="terms",se.fit=T)
res <- apply(a$fit,1,sum) +  attr(a$fit,"constant")
lines(newdat$yrqtr,inv.logit(res),col=2)
lines(newdat$yrqtr,inv.logit(res+2*a$se.fit[,4]),col=2,lty=2)
lines(newdat$yrqtr,inv.logit(res-2*a$se.fit[,4]),col=2,lty=2)
newdat <- expand.grid(mv=F,yrqtr=seq(1980.125,2011.875,.25),ll=dat$ll[1],trop=0)
a <- predict(model_xtrop,newdata=newdat,type="terms",se.fit=T)
res <- apply(a$fit,1,sum) +  attr(a$fit,"constant")
lines(newdat$yrqtr,inv.logit(res),col=3)
lines(newdat$yrqtr,inv.logit(res+2*a$se.fit[,3]),col=3,lty=2)
lines(newdat$yrqtr,inv.logit(res-2*a$se.fit[,3]),col=3,lty=2)
legend("bottomright",legend=c("Bigeye","Yellowfin","Bigeye + Yellowfin"),col=c(1,2,3),lty=1)
savePlot("Probability of moving after zero catch",type="png")

plot(1:40,1:40,xlim=c(1980,2011),ylim=c(0,1.5),xlab="Year",ylab="Odds of moving after zero vs mean catch",type="n")
newdat0 <- expand.grid(mv=F,yrqtr=seq(1980.125,2011.875,.25),ll=dat$ll[1],bet=0,yft=mean(dat$yft))
newdat1 <- expand.grid(mv=F,yrqtr=seq(1980.125,2011.875,.25),ll=dat$ll[1],bet=mean(dat$bet),yft=mean(dat$yft))
a0 <- predict(model_x,newdata=newdat0,type="response",se.fit=T)
a1 <- predict(model_x,newdata=newdat1,type="response",se.fit=T)
lines(newdat$yrqtr,a0$fit/a1$fit,col=1)
newdat0 <- expand.grid(mv=F,yrqtr=seq(1980.125,2011.875,.25),ll=dat$ll[1],bet=mean(dat$yft),yft=0)
newdat1 <- expand.grid(mv=F,yrqtr=seq(1980.125,2011.875,.25),ll=dat$ll[1],bet=mean(dat$bet),yft=mean(dat$yft))
a0 <- predict(model_x,newdata=newdat0,type="response",se.fit=T)
a1 <- predict(model_x,newdata=newdat1,type="response",se.fit=T)
lines(newdat$yrqtr,a0$fit/a1$fit,col=2)
newdat0 <- expand.grid(mv=F,yrqtr=seq(1980.125,2011.875,.25),ll=dat$ll[1],trop=0)
newdat1 <- expand.grid(mv=F,yrqtr=seq(1980.125,2011.875,.25),ll=dat$ll[1],trop=mean(dat$trop))
a0 <- predict(model_xtrop,newdata=newdat0,type="response",se.fit=T)
a1 <- predict(model_xtrop,newdata=newdat1,type="response",se.fit=T)
lines(newdat$yrqtr,a0$fit/a1$fit,col=3)
legend("bottomright",legend=c("Bigeye","Yellowfin","Bigeye + Yellowfin"),col=c(1,2,3),lty=1)

savePlot("Probability of moving after zero catch relative",type="png")

model_x <- glm((mv > 1) ~ yq + ll + yq*bet + yq*yft,data=dat,family="binomial")
model_xtrop <- glm((mv > 1) ~ yq + ll + yq*trop,data=dat,family="binomial")
plot(1:40,1:40,xaxt="n",xlim=c(0,140),ylim=c(0,2),xlab="Year",ylab="Odds of moving after zero vs mean catch",type="n")
axloc <- seq(6,132,20)
axis(1,at=axloc,labels=floor(as.numeric(as.character(sort(unique(dat$yq)))))[axloc])
newdat0 <- expand.grid(mv=F,yq=levels(dat$yq),ll=dat$ll[1],bet=0,yft=mean(dat$yft))
newdat1 <- expand.grid(mv=F,yq=levels(dat$yq),ll=dat$ll[1],bet=mean(dat$bet),yft=mean(dat$yft))
a0 <- predict(model_x,newdata=newdat0,type="response",se.fit=T)
a1 <- predict(model_x,newdata=newdat1,type="response",se.fit=T)
points(newdat0$yq,a0$fit/a1$fit,col=1,ylim=c(0,2),pch=1,cex=0.7)
lines(smooth.spline(a0$fit/a1$fit ~ newdat0$yq,n=4),col=1,lty=1,lwd=2)
newdat0 <- expand.grid(mv=F,yq=levels(dat$yq),ll=dat$ll[1],bet=mean(dat$yft),yft=0)
newdat1 <- expand.grid(mv=F,yq=levels(dat$yq),ll=dat$ll[1],bet=mean(dat$bet),yft=mean(dat$yft))
a0 <- predict(model_x,newdata=newdat0,type="response",se.fit=T)
a1 <- predict(model_x,newdata=newdat1,type="response",se.fit=T)
lines(smooth.spline(a0$fit/a1$fit ~ newdat0$yq,n=4),col=2,lty=1,lwd=2)
points(newdat0$yq,a0$fit/a1$fit,pch=1,col=2,cex=0.7)
newdat0 <- expand.grid(mv=F,yq=levels(dat$yq),ll=dat$ll[1],trop=0)
newdat1 <- expand.grid(mv=F,yq=levels(dat$yq),ll=dat$ll[1],trop=mean(dat$trop))
a0 <- predict(model_xtrop,newdata=newdat0,type="response",se.fit=T)
a1 <- predict(model_xtrop,newdata=newdat1,type="response",se.fit=T)
points(newdat0$yq,a0$fit/a1$fit,col=3,cex=0.7)
lines(smooth.spline(a0$fit/a1$fit ~ newdat0$yq,n=4),col=3,lty=1,lwd=2)
legend("bottomright",legend=c("Bigeye","Yellowfin","Bigeye + Yellowfin"),col=c(1,2,3),lty=1)
savePlot("Probability of moving after zero catch relative yq",type="png")


# median CPUE just after a move
plot_cpue_for_movements <- function(a01,a2,a0,a1,vbl,ti="") {
  plot(a01$yr,a01[,c(vbl)],type="l",xlab="Year",ylab="CPUE",ylim=c(0,max(c(a0[,c(vbl)],a1[,c(vbl)],a01[,c(vbl)],a2[,c(vbl)]))),main=ti)
  lines(a2$yr,a2[,vbl],col=2)
  lines(a0$yr,a0[,vbl],col=1,lty=2)
  lines(a1$yr,a1[,vbl],col=3)
  legend("bottomleft",legend=c("Move 0 squares","Move 0 or 1 square","Move 1 square","Move > 1 square"),col=c(1,1,3,2),lty=c(2,1,1,1))
  }
plot_cpue_ratios_for_movements <- function(a01,a2,a0,a1,vbl,ti="") {  
  plot(a0$yr,a2[,vbl]/a0[,vbl],type="l",xlab="Year",ylab="CPUE",ylim=c(0,1.6),main=ti)
  lines(a01$yr,a2[,vbl]/a01[,vbl],lty=2)
  legend("topleft",legend=c("Move > 1 square / Move 0 squares","Move > 1 square / Move 0 or 1 square"),col=c(1,1),lty=c(1,2))
  }

a2 <- aggregate(cbind(bet/hooks,yft/hooks,trop/hooks) ~ op_yr,dat[dat$mv > 1,],median,na.rm=T)
a1 <- aggregate(cbind(bet/hooks,yft/hooks,trop/hooks) ~ op_yr,dat[dat$mv == 1,],median,na.rm=T)
a0 <- aggregate(cbind(bet/hooks,yft/hooks,trop/hooks) ~ op_yr,dat[dat$mv == 0,],median,na.rm=T)
a01 <- aggregate(cbind(bet/hooks,yft/hooks,trop/hooks) ~ op_yr,dat[dat$mv <= 1,],median,na.rm=T)
names(a0) <- names(a1) <- names(a01) <- names(a2) <- c("yr","bet_cpue","yft_cpue","bet_yft_cpue")
windows(10,10);par(mfrow=c(3,2))
plot_cpue_for_movements(a01,a2,a0,a1,vbl="bet_cpue",ti="Bigeye CPUE")
plot_cpue_ratios_for_movements(a01,a2,a0,a1,vbl="bet_cpue",ti="Bigeye CPUE")
plot_cpue_for_movements(a01,a2,a0,a1,vbl="yft_cpue",ti="Yellowfin CPUE")
plot_cpue_ratios_for_movements(a01,a2,a0,a1,vbl="yft_cpue",ti="Yellowfin CPUE")
plot_cpue_for_movements(a01,a2,a0,a1,vbl="bet_yft_cpue",ti="Tuna CPUE")
plot_cpue_ratios_for_movements(a01,a2,a0,a1,vbl="bet_yft_cpue",ti="Tuna CPUE")
savePlot("just_after_move",type="png")
a2 <- aggregate(cbind(bet,yft,trop) ~ op_yr,dat[dat$mv > 1,],median,na.rm=T)
a1 <- aggregate(cbind(bet,yft,trop) ~ op_yr,dat[dat$mv == 1,],median,na.rm=T)
a0 <- aggregate(cbind(bet,yft,trop) ~ op_yr,dat[dat$mv == 0,],median,na.rm=T)
a01 <- aggregate(cbind(bet,yft,trop) ~ op_yr,dat[dat$mv <= 1,],median,na.rm=T)
names(a0) <- names(a1) <- names(a01) <- names(a2) <- c("yr","bet","yft","bet_yft")
windows(10,10);par(mfrow=c(3,2))
plot_cpue_for_movements(a01,a2,a0,a1,vbl="bet",ti="Bigeye CPUE")
plot_cpue_ratios_for_movements(a01,a2,a0,a1,vbl="bet",ti="Bigeye CPUE")
plot_cpue_for_movements(a01,a2,a0,a1,vbl="yft_cpue",ti="Yellowfin CPUE")
plot_cpue_ratios_for_movements(a01,a2,a0,a1,vbl="yft",ti="Yellowfin CPUE")
plot_cpue_for_movements(a01,a2,a0,a1,vbl="bet_yft",ti="Tuna CPUE")
plot_cpue_ratios_for_movements(a01,a2,a0,a1,vbl="bet_yft",ti="Tuna CPUE")
savePlot("just_after_move_catchperset",type="png")

# median CPUE just before a move
a2 <- aggregate(cbind(bet/hooks,yft/hooks,trop/hooks) ~ op_yr,dat[dat$nxtmv > 1,],median,na.rm=T)
a1 <- aggregate(cbind(bet/hooks,yft/hooks,trop/hooks) ~ op_yr,dat[dat$nxtmv == 1,],median,na.rm=T)
a0 <- aggregate(cbind(bet/hooks,yft/hooks,trop/hooks) ~ op_yr,dat[dat$nxtmv == 0,],median,na.rm=T)
a01 <- aggregate(cbind(bet/hooks,yft/hooks,trop/hooks) ~ op_yr,dat[dat$nxtmv <= 1,],median,na.rm=T)
names(a0) <- names(a1) <- names(a01) <- names(a2) <- c("yr","bet_cpue","yft_cpue","bet_yft_cpue")
windows(10,10);par(mfrow=c(3,2))
plot_cpue_for_movements(a01,a2,a0,a1,vbl="bet_cpue",ti="Bigeye CPUE")
plot_cpue_ratios_for_movements(a01,a2,a0,a1,vbl="bet_cpue",ti="Bigeye CPUE")
plot_cpue_for_movements(a01,a2,a0,a1,vbl="yft_cpue",ti="Yellowfin CPUE")
plot_cpue_ratios_for_movements(a01,a2,a0,a1,vbl="yft_cpue",ti="Yellowfin CPUE")
plot_cpue_for_movements(a01,a2,a0,a1,vbl="bet_yft_cpue",ti="Tuna CPUE")
plot_cpue_ratios_for_movements(a01,a2,a0,a1,vbl="bet_yft_cpue",ti="Tuna CPUE")
savePlot("just_before_move",type="png")

a2 <- aggregate(cbind(bet,yft,trop) ~ op_yr,dat[dat$nxtmv > 1,],median,na.rm=T)
a1 <- aggregate(cbind(bet,yft,trop) ~ op_yr,dat[dat$nxtmv == 1,],median,na.rm=T)
a0 <- aggregate(cbind(bet,yft,trop) ~ op_yr,dat[dat$nxtmv == 0,],median,na.rm=T)
a01 <- aggregate(cbind(bet,yft,trop) ~ op_yr,dat[dat$nxtmv <= 1,],median,na.rm=T)
names(a0) <- names(a1) <- names(a01) <- names(a2) <- c("yr","bet","yft","bet_yft")
windows(10,10);par(mfrow=c(3,2))
plot_cpue_for_movements(a01,a2,a0,a1,vbl="bet",ti="Bigeye CPUE")
plot_cpue_ratios_for_movements(a01,a2,a0,a1,vbl="bet",ti="Bigeye CPUE")
plot_cpue_for_movements(a01,a2,a0,a1,vbl="yft",ti="Yellowfin CPUE")
plot_cpue_ratios_for_movements(a01,a2,a0,a1,vbl="yft",ti="Yellowfin CPUE")
plot_cpue_for_movements(a01,a2,a0,a1,vbl="bet_yft",ti="Tuna CPUE")
plot_cpue_ratios_for_movements(a01,a2,a0,a1,vbl="bet_yft",ti="Tuna CPUE")
savePlot("just_before_move_catchperset",type="png")

# patchiness in set distribution across fleet
windows(10,10);par(mfrow=c(2,2))
cu <- function(a) length(unique(a))
uv <- aggregate(vessid ~ sdate + yrqtr,data=dat,FUN=cu)
usq <- aggregate(paste(lat,lon) ~ sdate + yrqtr,data=dat,FUN=cu)
uvm <- tapply(uv$vessid,uv$yrqtr,median)
usqm <- tapply(usq[,3],uv$yrqtr,median)
tm <- as.numeric(names(usqm))
plot(tm,usqm/uvm,type="l",ylim=c(0,1.1),xlab="Month",ylab="1degsq per vessel")
lines(tm,uvm/max(uvm),col=2)
savePlot("plot_1deg_squares_per_vessel",type="png")

usq <- aggregate(latlong ~ sdate + yrqtr,data=dat,FUN=cu)
uvm <- tapply(uv$vessid,uv$yrqtr,median)
usqm <- tapply(usq[,3],uv$yrqtr,median)
tm <- as.numeric(names(usqm))
plot(tm,usqm/uvm,type="l",ylim=c(0,1.1),xlab="Month",ylab="5degsq per vessel")
lines(tm,uvm/max(uvm),col=2)
savePlot("plot_5deg_squares_per_vessel",type="png")

# look at fishing in same square in same month
uv <- aggregate(vessid ~ yrqtr,data=dat,FUN=cu)
usq <- aggregate(paste(lat,lon) ~ yrqtr,data=dat,FUN=cu)
tm <- usq$yrqtr
plot(tm,usq[,2]/uv[,2],type="l",xlab="Month",ylab="1degsq per vessel",ylim=c(0,max(usq[,2]/uv[,2])))
par(new=T)
plot(tm,uvm,axes=F,col=2,xlab="",ylab="",ylim=c(0,max(uvm)))
axis(4,col=2)
savePlot("plot_1deg_squares_per_vessel_month",type="png")

usq <- aggregate(latlong ~ yrqtr,data=dat,FUN=cu)
uv <- aggregate(vessid ~ yrqtr,data=dat,FUN=cu)
tm <- uv$yrqtr 
plot(tm,usqm/uvm,type="l",ylim=c(0,1.1),xlab="Month",ylab="5degsq per vessel")
lines(tm,uvm/max(uvm),col=2)
savePlot("plot_5deg_squares_per_vessel_month",type="png")

# look at fishing in same square in same month
uv <- aggregate(vessid ~ yrqtr,data=dat,FUN=cu); 
uvsq <- aggregate(paste(vessid,lat,lon) ~ yrqtr,data=dat,FUN=cu);names(uvsq)[2] <- "vess_sq"
usq <- aggregate(paste(lat,lon) ~ yrqtr,data=dat,FUN=cu)
vsq <- aggregate(paste(vessid,lat,lon) ~ yrqtr,data=dat,FUN=length);names(vsq)[2] <- "sets_per_qtr"
tm <- usq$yrqtr
windows(10,10);par(mfrow=c(2,2))
plot(tm,usq[,2]/uv[,2],ylim=c(0,7),type="l",xlab="Month",ylab="1degsq per vessel")
lines(tm,uvm/max(uvm),col=2)
plot(tm,uvsq$vess_sq/uv$vessid,type="l",xlab="Month",ylab="1deg sq / vessel / qtr",ylim=c(0,max(uvsq$vess_sq/uv$vessid)))
par(add=T)
plot(tm,vsq$sets_per_qtr/uvsq$vess_sq,type="l",xlab="Month",ylab="Sets per vessel per square",ylim=c(0,4))
savePlot("plot_1deg_squares_per_vessel_month",type="png")

usq <- aggregate(latlong ~ yrqtr,data=dat,FUN=cu)
uv <- aggregate(vessid ~ yrqtr,data=dat,FUN=cu)
tm <- as.Date(paste(names(usqm),"-01",sep=""),"%Y-%m-%d") 
plot(tm,usqm/uvm,type="l",ylim=c(0,1.1),xlab="Month",ylab="5degsq per vessel")
lines(tm,uvm/max(uvm),col=2)
savePlot("plot_5deg_squares_per_vessel_month",type="png")


# model GLM catch ratios
glmdat <- select_data(opx,runreg=9,minqtrs=minqtrs,runsp="bet",mt="propn",llstrat=3,addmain=F,fc="both",doboth=T)
wtt   <- mk_wts(glmdat,wttype="propn",catch=glmdat$yft + glmdat$bet)
fmla <- make_formula(runsp,modtype="propn",addboat=T,addmain=F,addbranch=F)
glmdat$bet <- glmdat$bet/wtt
model.propn   <- glm(fmla,data=glmdat,weights=wtt,family="binomial");gc()
glmdat <- glmdat[!is.na(glmdat[,4]),]
abin <- length(unique(glmdat$yrqtr))
coefs <- inv.logit(get.coefs(model.propn,abin))
plot(coefs,type="l")
plot_effects(model=model.propn,indat=glmdat,addmain=F,addbranch=F,addalb=F,addother=F,ti="")
savePlot("Proportion_bigeye_effectplots_R3",type="png")
summ_pbetmod <- summary(model.propn)
save(summ_pbetmod,file="summod_prop_bet_R3.RData")

glmdat <- select_data(opx79,runreg=9,minqtrs=minqtrs,runsp="bet",mt="propn",llstrat=3,addmain=F,fc="OS",doboth=T)
wtt   <- mk_wts(glmdat,wttype="propn",catch=glmdat$yft + glmdat$bet)
fmla <- make_formula(runsp,modtype="propn",addboat=T,addmain=F,addbranch=F)
glmdat$bet <- glmdat$bet/wtt
model.propn   <- glm(fmla,data=glmdat,weights=wtt,family="binomial");gc()
glmdat <- glmdat[!is.na(glmdat[,4]),]
abin <- length(unique(glmdat$yrqtr))
coefs <- inv.logit(get.coefs(model.propn,abin))
plot(coefs,type="l")
plot_effects(model=model.propn,indat=glmdat,addmain=F,addbranch=F,addalb=F,addother=F,ti="")
savePlot("Proportion_bigeye_effectplots_R3_79",type="png")
summ_pbetmod <- summary(model.propn)
save(summ_pbetmod,file="summod_prop_bet_R3_79.RData")

glmdat <- select_data(opx79,runreg=9,minqtrs=minqtrs,runsp="bet",mt="propn",llstrat=3,addmain=T,fc="both",doboth=T)
wtt   <- mk_wts(glmdat,wttype="propn",catch=glmdat$yft + glmdat$bet)
fmla <- make_formula(runsp,modtype="propn",addboat=T,addmain=T,addbranch=T)
glmdat$bet <- glmdat$bet/wtt
model.propn   <- glm(fmla,data=glmdat,weights=wtt,family="binomial");gc()
glmdat <- glmdat[!is.na(glmdat[,4]),]
abin <- length(unique(glmdat$yrqtr))
coefs <- inv.logit(get.coefs(model.propn,abin))
plot(coefs,type="l")
plot_effects(model=model.propn,indat=glmdat,addmain=T,addbranch=T,addalb=F,addother=F,ti="")
savePlot("Proportion_bigeye_effectplots_R3_79_mainbrnch",type="png")
summ_pbetmod <- summary(model.propn)
save(summ_pbetmod,file="summod_prop_bet_R3_79_mainbrnch.RData")

glmdat <- select_data(opx,runreg=4,minqtrs=minqtrs,runsp="bet",mt="propn",llstrat=3,addmain=F,fc="both",doboth=T)
wtt   <- mk_wts(glmdat,wttype="propn",catch=glmdat$yft + glmdat$bet)
fmla <- make_formula(runsp,modtype="propn",addboat=T,addmain=F,addbranch=F)
glmdat$bet <- glmdat$bet/wtt
model.propn   <- glm(fmla,data=glmdat,weights=wtt,family="binomial");gc()
glmdat <- glmdat[!is.na(glmdat[,4]),]
abin <- length(unique(glmdat$yrqtr))
coefs <- inv.logit(get.coefs(model.propn,abin))
plot(coefs,type="l")
plot_effects(model=model.propn,indat=glmdat,addmain=F,addbranch=F,addalb=F,addother=F,ti="")
savePlot("Proportion_bigeye_effectplots_R4",type="png")
summ_pbetmod <- summary(model.propn)
save(summ_pbetmod,file="summod_prop_bet_R4.RData")

glmdat <- select_data(opx79,runreg=4,minqtrs=minqtrs,runsp="bet",mt="propn",llstrat=3,addmain=F,fc="both",doboth=T)
wtt   <- mk_wts(glmdat,wttype="propn",catch=glmdat$yft + glmdat$bet)
fmla <- make_formula(runsp,modtype="propn",addboat=T,addmain=F,addbranch=F)
glmdat$bet <- glmdat$bet/wtt
model.propn   <- glm(fmla,data=glmdat,weights=wtt,family="binomial");gc()
glmdat <- glmdat[!is.na(glmdat[,4]),]
abin <- length(unique(glmdat$yrqtr))
coefs <- inv.logit(get.coefs(model.propn,abin))
plot(coefs,type="l")
plot_effects(model=model.propn,indat=glmdat,addmain=F,addbranch=F,addalb=F,addother=F,ti="")
savePlot("Proportion_bigeye_effectplots_R4_79",type="png")
summ_pbetmod <- summary(model.propn)
save(summ_pbetmod,file="summod_prop_bet_R4_79.RData")

glmdat <- select_data(opx79,runreg=4,minqtrs=minqtrs,runsp="bet",mt="propn",llstrat=3,addmain=T,fc="both",doboth=T)
wtt   <- mk_wts(glmdat,wttype="propn",catch=glmdat$yft + glmdat$bet)
fmla <- make_formula(runsp,modtype="propn",addboat=T,addmain=T,addbranch=T)
glmdat$bet <- glmdat$bet/wtt
model.propn   <- glm(fmla,data=glmdat,weights=wtt,family="binomial");gc()
glmdat <- glmdat[!is.na(glmdat[,4]),]
abin <- length(unique(glmdat$yrqtr))
coefs <- inv.logit(get.coefs(model.propn,abin))
plot(coefs,type="l")
plot_effects(model=model.propn,indat=glmdat,addmain=T,addbranch=T,addalb=F,addother=F,ti="")
savePlot("Proportion_bigeye_effectplots_R4_79_mainbrnch",type="png")
summ_pbetmod <- summary(model.propn)
save(summ_pbetmod,file="summod_prop_bet_R4_79_mainbrnch.RData")
