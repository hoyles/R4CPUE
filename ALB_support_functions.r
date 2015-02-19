symbol_size_data <- function(a,symsize=0.2,colour=4) {
  if(dim(a)[1]>0) {
    a <- aggregate(a$setwt,list(a$centlat,a$centlong),sum)
#    a <- a[a$x>0,]
    symbols(a[,2],a[,1],circles=sqrt(a[,3]),add=T,inches=symsize,bg=colour)
  }
}  

clean_logbook_data <- function(dat) {
  dat <- dat[,c("trip_id","set_id","flag_id","logdate","lat_long","lon_long","ez_id",
                "in_arch","spc","sch_id","days_sch","set_start","set_end","skj_c","yft_c","bet_c",
                "skj_w","yft_w","bet_w","yft_c_adj","bet_c_adj")]
  return(dat)
}
   
prepare_logbook_data <- function(dat) {
  dat$latd <- as.numeric(substring(dat$lat, 1, 4))/100
  dat$latd <- floor(dat$latd) + (dat$latd-floor(dat$latd))/0.6
  dat$latd <- ifelse(substring(dat$lat, 5, 5) == "S", -1 * dat$latd, dat$latd)
  dat <- dat[dat$latd <= 0,] 
  
  dat$lond <- as.numeric(substring(dat$lon, 1, 5))/100
  dat$lond <- floor(dat$lond) + (dat$lond-floor(dat$lond))/0.6
  dat$lond <- ifelse(substring(dat$lon, 6, 6) == "W", 360-dat$lond, dat$lond)
  dat <- dat[dat$lond>=140 & dat$lond<250,]

  dat$lat5 <- 5*floor((dat$latd+5)/5)-5
  dat$lon5 <- 5*floor(dat$lond/5)
  dat$latlong <- as.factor(paste(dat$lat5,dat$lon5))
  dat$centlat <- dat$lat5 + 2.5
  dat$centlong <- dat$lon5 + 2.5
  dat <- dat[dat$centlat > -50,]
  dat <- dat[dat$centlat < 50,]
  dat <- dat[dat$centlong > 100,]
  dat <- dat[dat$centlong < 310,]
  
  dat$reg <- rep(0, length(dat$centlat))
  dat$reg <- ifelse(dat$centlat >= -25 & dat$centlat <   0 & dat$centlong >= 140 & dat$centlong < 180, 1, dat$reg)
  dat$reg <- ifelse(dat$centlat >= -25 & dat$centlat <   0 & dat$centlong >= 180 & dat$centlong < 250, 2, dat$reg)
  dat$reg <- ifelse(dat$centlat >= -40 & dat$centlat < -25 & dat$centlong >= 140 & dat$centlong < 180, 3, dat$reg)
  dat$reg <- ifelse(dat$centlat >= -40 & dat$centlat < -25 & dat$centlong >= 180 & dat$centlong < 250, 4, dat$reg)
  dat$reg <- ifelse(dat$centlat >= -25 & dat$centlat <   0 & dat$centlong >= 250 & dat$centlong < 330, 5, dat$reg)
  dat$reg <- ifelse(dat$centlat >= -40 & dat$centlat < -25 & dat$centlong >= 250 & dat$centlong < 330, 6, dat$reg)

  dat <- dat[dat$reg %in% 1:4,]
  
  return(dat)
}

prepare_logbook_data_2 <- function(a) {
  a$trip_date <- as.Date(a$trip_date,format="%d/%m/%Y")
  a$trip_yr <- as.numeric(format(a$trip_date,"%Y"))
  a$trip_mon <- as.numeric(format(a$trip_date,"%m"))
  a$set_date <- as.Date(a$logdate,format="%d/%m/%Y")
  a$set_yr <- as.numeric(format(a$set_date,"%Y"))
  a$set_mon <- as.numeric(format(a$set_date,"%m"))
  a$set_qtr <- floor((a$set_mon+2)/3)
  a$set_yrqtr <- a$set_yr + a$set_qtr/4 - 0.125
  return(a)
  }


spatial_overlap_index <- function(samp,catch){
  latlevs <- seq(-12.5,12.5,5);lonlevs <- seq(127.5,207.5,5);
  t1 <- with(samp,tapply(setwt,list(factor(centlat,levels=latlevs),factor(centlong,levels=lonlevs)),sum))
  t2 <- with(catch,tapply(setwt,list(factor(centlat,levels=latlevs),factor(centlong,levels=lonlevs)),sum))
  t1[is.na(t1)] <- 0
  t2[is.na(t2)] <- 0
  r1 <- t1 / sum(t1)
  r2 <- t2 / sum(t2)
  d1 <- sum(abs(r1-r2))
  d2 <- sum(log(((r1[r2>0]-r2[r2>0])^2)/r2[r2>0]))
  d3 <- (sum(t2)/sum(t1))/(sum(t2[t1>0]/(t1[t1>0]))/sum(t1>0))
  return(data.frame(absdiff=d1,propLik=d2,Gulland=d3))
}

symbol_catch_data <- function(a,symsize=0.2,colour=4) {
  if(dim(a)[1]>0) {
    a <- aggregate(a$setwt,list(a$centlat,a$centlong),sum)
    symbols(a[,2],a[,1],circles=sqrt(a[,3]),add=T,inches=symsize,bg=colour)
  }
}  
