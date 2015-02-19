# CPUE_analysis_2012
# SDH starting 2012-06-10
# 1. Extract data using 2012_data_prep.r
# 2. A bit of pre-processing using intro_processing_pl.r
# 2b. Bit more cleaning
# 3. Split up by region
# 4. Cluster a region by trip using kmeans
# 5. Link each set with a cluster
# 6. Plot the sets on a map, coloured by cluster
# 7. Standardize all the clusters for each region


setwd("D:/L/alb/2012/assessment/Data_preparation/CPUE/")
# setwd("X:/L/alb/2012/assessment/Data_preparation/CPUE")
load("2012_allfleets2.RData")
length(unique(pl$boat_id))
ls()
dim(pl)
library(rpart)
library(cluster)
require(pkpkg)
require(R4MFCL)
require(MASS)

ezpac <- read.table("I:/assessments/Pop dy modeling/MFCL/R functions/EZNEW2.TXT")
names(ezpac) <- c("lon","lat")

windows()
hist(pl$set_yr)
boxplot(pl$sets_per_trip~pl$set_yr) # high in 1960-61 because of two vessels with sets in both yrs, need to assign to start(trip_date) and graph    #GMP: also increased across time to 2005...
windows();par(mfrow=c(2,1))
boxplot(pl$sets_per_trip~pl$flag_id)
barplot(table(pl$flag_id))

tp <- aggregate(cbind(trip_yr,sets_per_trip) ~ tripest_id + flag_id,data=pl,max)
windows(16,16);par(mfrow=c(4,4),mar=c(2,2,2,0))
boxplot(tp$sets_per_trip~tp$trip_yr)
for(f in sort(unique(tp$flag_id))) {
  with(tp[tp$flag_id ==f,],boxplot(sets_per_trip ~ trip_yr,main=f,ylim=c(0,250)))
  }
savePlot("sets_per_trip_by_flag.png",type="png")

# Remove redundant fields to clean up the data frame
pl <- pl[,c("trip_id","tripest_id","set_id","flag_id","ez_id","hook","hk_bt_flt","alb_n","yft_n","bet_n","oth_n","latd","lond","lat5","lon5","latlong","centlat","centlong",
          "reg","boat_id","src_id","fleet_id","swo","bum","mls","blm","wah","dol","skj","bsh","sbf","shk","uns",
          "trip_yr","trip_mon","set_yr","set_mon","set_qtr","set_yrqtr","tuna_n","bet_cpue","alb_cpue","yft_cpue","bet_perc","alb_perc","yft_perc","wt")]
table(pl$reg)

# Some plotting to see what's going on. 
windows(width=20,height=18); par(mfrow=c(2,2),mar=c(2,2,1,1))
for(rg in 1:4) {
  reg <- pl[pl$reg==rg,]
  boxplot(alb_cpue ~ eval(floor(hook/200)*200),data=reg,main=paste("Region",rg),ylim=c(0,400))
  }
savePlot("alb_cpue_by_hooks",type="png")

windows(width=20,height=18); par(mfrow=c(2,2),mar=c(2,2,1,1))
for(rg in 1:4) {
  reg <- pl[pl$reg==rg,]
  boxplot(alb_n ~ eval(floor(hook/200)*200),data=reg,main=paste("Region",rg))
  }
savePlot("alb_catch_by_hooks",type="png")

windows(width=20,height=18); par(mfrow=c(2,2),mar=c(2,2,1,1))
for(rg in 1:4) {
  reg <- pl[pl$reg==rg,]
  boxplot(hook ~ set_yr,data=reg,main=paste("Region",rg))
  }

a <- pl[pl$set_yr >= 2000,]
for(rg in 1:4) {
  reg <- a[a$reg==rg,]
  boxplot(alb_cpue ~ eval(floor(hook/200)*200),data=reg,main=paste("Region",rg))
  }

windows(width=20,height=18); par(mfrow=c(2,2),mar=c(2,2,1,1))
a <- pl[pl$set_yr > 1980,]
for(rg in 1:4) {
  reg <- a[a$reg==rg,]
  boxplot(alb_cpue ~ eval(floor(hook/200)*200),data=reg,main=paste("Region",rg),ylim=c(0,400))
  }

windows(width=20,height=18); par(mfrow=c(2,2),mar=c(2,2,1,1))
for(rg in 1:4) {
  reg <- pl[pl$reg==rg,]
  boxplot(alb_cpue ~ set_yr,data=reg,main=paste("Region",rg),ylim=c(0,400))
  }
windows(width=20,height=18); par(mfrow=c(2,2),mar=c(2,2,1,1))
for(rg in 1:4) {
  reg <- pl[pl$reg==rg,]
  hist(reg$hook,main=paste("Region",rg))
  }

windows(width=20,height=18); par(mfrow=c(2,2),mar=c(2,2,1,1))
for(rg in 1:4) {
  reg <- pl[pl$reg==rg,]
#  hist(reg$alb_n,main=paste("Region",rg))
  a <- table(reg$alb_n)
  plot(as.numeric(names(a)),log(a))
  }
for(rg in 1:4) {
  reg <- pl[pl$reg==rg,]
  print(quantile(reg$alb_n,c(.001,.01,.05,.5,.95,.99,.999)))
  }
a <- pl[pl$alb_n > 700,]
windows()
a <- pl[pl$set_yr >2000,]
a <- pl[pl$alb_n >400,]
plot(a$alb_n,a$alb_c)   # alb_n seems to match alb_c. Not sure about alb_w
table(ceiling(a$alb_w/10),ceiling(a$alb_n/50))

# remove sets without many hooks
a <- pl[pl$hook > 1000,]
length(unique(a$boat_id))
table(pl$flag_id)
table(a$flag_id)
barplot(table(a$flag_id)/table(pl$flag_id),xlab="Flag",ylab="Proportion of sets remaining",ylim=c(0,1))
savePlot("Effect of 1000 hook limit on sets",type="png")
barplot(tapply(a$hook,a$flag_id,sum)/tapply(pl$hook,pl$flag_id,sum),xlab="Flag",ylab="Proportion of effort remaining",ylim=c(0,1))
savePlot("Effect of 1000 hook limit on hooks",type="png")

pl2 <- pl[pl$hook > 1000,]    #

#### Set up variables for trip-level clustering
tpall <- aggregate(cbind(hook,hk_bt_flt,alb_n,yft_n,bet_n,oth_n,latd,lond,swo,bum,mls,blm,wah,dol,skj,bsh,sbf,shk,uns) ~ tripest_id + flag_id + boat_id + trip_yr + reg,data=pl2,mean)
tpall$bet_cpue <- 1000*tpall$bet_n/tpall$hook
tpall$yft_cpue <- 1000*tpall$yft_n/tpall$hook
tpall$alb_cpue <- 1000*tpall$alb_n/tpall$hook
tpall$oth_cpue <- 1000*tpall$oth_n/tpall$hook
tpall$swo_cpue <- 1000*tpall$swo/tpall$hook
tpall$bum_cpue <- 1000*tpall$bum/tpall$hook
tpall$mls_cpue <- 1000*tpall$mls/tpall$hook
tpall$blm_cpue <- 1000*tpall$blm/tpall$hook
tpall$wah_cpue <- 1000*tpall$wah/tpall$hook
tpall$dol_cpue <- 1000*tpall$dol/tpall$hook
tpall$skj_cpue <- 1000*tpall$skj/tpall$hook
tpall$bsh_cpue <- 1000*tpall$bsh/tpall$hook
tpall$sbf_cpue <- 1000*tpall$sbf/tpall$hook
tpall$shk_cpue <- 1000*tpall$shk/tpall$hook
tpall$uns_cpue <- 1000*tpall$uns/tpall$hook

tp <- aggregate(cbind(bet_n,yft_n,alb_n,oth_n,swo,bum,mls,blm,wah,dol,skj,bsh,sbf,shk,uns,hk_bt_flt) ~ tripest_id + flag_id + boat_id + trip_yr + reg,data=pl2,mean)
a <- aggregate(bet_n ~ tripest_id + flag_id + boat_id + trip_yr + reg,data=pl2,length)
dima <- dim(tp)[2]
tp <- cbind(tp,a[,6])
names(tp)[dima+1] <- "nset"
tp$tot_n <- tp$alb_n + tp$yft_n + tp$bet_n + tp$oth_n
tp$tuna_n <- tp$alb_n + tp$yft_n + tp$bet_n
tp$bet_perc <- tp$bet_n/tp$tuna_n
tp$alb_perc <- tp$alb_n/tp$tuna_n
tp$yft_perc <- tp$yft_n/tp$tuna_n
tp$bet_totpc <- tp$bet_n/tp$tot_n
tp$alb_totpc <- tp$alb_n/tp$tot_n
tp$yft_totpc <- tp$yft_n/tp$tot_n
tp$oth_totpc <- tp$oth_n/tp$tot_n
tp$swo_totpc <- tp$swo/tp$tot_n
tp$bum_totpc <- tp$bum/tp$tot_n
tp$mls_totpc <- tp$mls/tp$tot_n
tp$blm_totpc <- tp$blm/tp$tot_n
tp$wah_totpc <- tp$wah/tp$tot_n
tp$dol_totpc <- tp$dol/tp$tot_n
tp$skj_totpc <- tp$skj/tp$tot_n
tp$bsh_totpc <- tp$bsh/tp$tot_n
tp$sbf_totpc <- tp$sbf/tp$tot_n
tp$shk_totpc <- tp$shk/tp$tot_n
tp$uns_totpc <- tp$uns/tp$tot_n
gc()

str(tp)
pairs(tp[,21:32])
for(r in 1:4) {
  windows(height=20,width=30); par(mfrow=c(4,4))
  hist(tp[tp$reg==r,25:39])
  title(paste("Region",r),outer=T,line=-1)
  savePlot(paste("hists_by_sp_reg_",rg,sep=""),type="png")
  }
  
# Regression trees to see which factors affect alb cpue
windows(width=20,height=18); par(mfrow=c(2,2),mar=c(2,2,2,2))
for(rg in 1:4) {
  a <- tpall[tpall$reg==rg,]
  tree1 <- rpart(alb_cpue ~ latd + lond + trip_yr + flag_id + bet_cpue + yft_cpue + swo_cpue + bsh_cpue + hook, data=a)
  plot(tree1,main=paste("Region",rg))
  text(tree1)
  }
savePlot("Rtrees_by_reg_by_trip",type="png")

windows(width=20,height=18); par(mfrow=c(2,2),mar=c(2,2,2,2))
for(rg in 1:4) {
  a <- pl2[pl2$reg==rg,]
  tree1 <- rpart(alb_cpue ~ latd + lond + trip_yr + flag_id + bet_cpue + yft_cpue + hook, data=a)
  plot(tree1,main=paste("Region",rg))
  text(tree1)
  }
savePlot("Rtrees_by_reg_by_set",type="png")

# Regression trees to see which factors affect alb percent
windows(width=22,height=18); par(mfrow=c(2,2),mar=c(2,2,2,2))
for(rg in 1:4) {
  a <- tpall[tpall$reg==rg,]
  a$alb_perc <- with(a,alb_n / (alb_n + yft_n + bet_n + oth_n))
  tree1 <- rpart(alb_perc ~ latd + lond + trip_yr + flag_id + hook, data=a)
  plot(tree1,main=paste("Region",rg))
  text(tree1)
  }
savePlot("Rtrees_albperc_by_reg_by_trip",type="png")

windows(width=22,height=18); par(mfrow=c(2,2),mar=c(2,2,2,2))
for(rg in 1:4) {
  a <- pl2[pl2$reg==rg,]
  a$alb_perc <- with(a,alb_n / (alb_n + yft_n + bet_n + oth_n))
  tree1 <- rpart(alb_cpue ~ latd + lond + trip_yr + flag_id + hook, data=a)
  plot(tree1,main=paste("Region",rg))
  text(tree1)
  }
savePlot("Rtrees_albperc_by_reg_by_set",type="png")

# Look at numbers of clusters
tpx <- tp[,c("bet_perc","yft_perc","alb_perc","oth_totpc","trip_yr","tripest_id", "hk_bt_flt","reg")]
triptots <- tp[,c("bet_n","yft_n","alb_n","oth_n","trip_yr","tripest_id", "hk_bt_flt","reg")]
tpx <- na.omit(tpx)
windows(width=20,height=18); par(mfrow=c(2,2),mar=c(2,2,2,2))
nsp=3
for(rg in 1:4) {
  a <- tpx[tpx$reg==rg,]
  a[,1:nsp] <- scale(a[,1:nsp])
  wss <- (nrow(a)-1)*sum(apply(a[,1:nsp],2,var))
  for (i in 2:15) wss[i] <- sum(kmeans(a[,1:nsp],centers=i,iter.max = 40)$withinss)
  plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares",main=paste("Region",rg))
  }
savePlot(paste("cluster_num_v_SS_by_reg_",nsp,"_spp",sep=""),type="png")

# identify catch rates by trip
for(rg in 1:4) {
  a <- tp[tp$reg==rg,]
  windows(width=20,height=18); par(mfrow=c(2,2))
  boxplot(a$bet_perc~a$trip_yr, xlab="Year", ylab="bet_perc")   #looks like high BET C prop prior to ~1984, and post 2002
  boxplot(a$alb_perc~a$trip_yr, xlab="Year", ylab="bet_perc")
  boxplot(a$yft_perc~a$trip_yr, xlab="Year", ylab="bet_perc")
  boxplot(a$oth_totpc~a$trip_yr, xlab="Year", ylab="bet_perc")
  title(paste("Region",rg),outer=T,line=-1)
  savePlot(paste("boxplot_sp_percent_reg_",rg,sep=""),type="png")
}

## GMP: start here
table(tp$tuna_n==0)
tp <- tp[tp$tuna_n!=0, ]   # removing 82 trips


###################################
## GMP note - choose which of these periods you want to look at - incl post 1999, 1990-1999, or the whole kaboosh (not done here)
## Start with post 1999
## SDH change later to make it pre-1984, 1984-2002, and 2003-2011
## SDH change again to do all at once
###################################
#junk5 <- junk3[junk3$ret_year>=1999,]
#dim(junk5)         #[1] 469   10

for(reg in 1:4) { # This takes a long time to run for R1 and R2
  gc()
  ax2 <- tp[tp$reg==reg,c("bet_perc","yft_perc","alb_perc","trip_yr","tripest_id", "hk_bt_flt")]
  ax <- tp[tp$reg==reg & tp$nset >= 5,c("bet_perc","yft_perc","alb_perc","trip_yr","tripest_id", "hk_bt_flt")]
  triptots2 <- tp[tp$reg==reg,c("bet_n","yft_n","alb_n","trip_yr","tripest_id", "hk_bt_flt")]
  triptots <- tp[tp$reg==reg & tp$nset >= 5,c("bet_n","yft_n","alb_n","trip_yr","tripest_id", "hk_bt_flt")]
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
  save(a,d,fit,clarax,file=paste("reg_",reg,"_fitfiles.RData",sep=""))
  plot(clarax)
  savePlot(file=paste("cluster_allyrs_dendrogram_R",reg,".png",sep=""),type="png")
}

for(reg in 1:4) { 
  load(file=paste("reg_",reg,"_fitfiles.RData",sep=""))
  print(reg)
  print(aggregate(a[,1:3],by=list(clarax$clustering),FUN=mean))
  print(aggregate(clarax$data,by=list(clarax$clustering),FUN=mean))
    
#[1] 1
#  Group.1   bet_perc   yft_perc  alb_perc
#1       1 0.02377515 0.05786726 0.9183576
#2       2 0.07744380 0.19064415 0.7319121
#3       3 0.12282213 0.49362036 0.3835575
#  Group.1   bet_perc   yft_perc  alb_perc
#1       1 0.02377515 0.05786726 0.9183576
#2       2 0.07744380 0.19064415 0.7319121
#3       3 0.12282213 0.49362036 0.3835575
#
#   1    2    3 
#7807 6329 4215 
#[1] 2
#  Group.1   bet_perc  yft_perc  alb_perc
#1       1 0.11961166 0.1497490 0.7306394
#2       2 0.03894412 0.0479472 0.9131087
#3       3 0.23057253 0.4102599 0.3591676
#  Group.1   bet_perc  yft_perc  alb_perc
#1       1 0.11961166 0.1497490 0.7306394
#2       2 0.03894412 0.0479472 0.9131087
#3       3 0.23057253 0.4102599 0.3591676
#
#   1    2    3 
#7617 7698 5839 
#[1] 3
#  Group.1    bet_perc    yft_perc  alb_perc
#1       1 0.026824129 0.026789908 0.9463860
#2       2 0.008043066 0.005833925 0.9861230
#3       3 0.051543031 0.101308324 0.8471486
#  Group.1    bet_perc    yft_perc  alb_perc
#1       1 0.026824129 0.026789908 0.9463860
#2       2 0.008043066 0.005833925 0.9861230
#3       3 0.051543031 0.101308324 0.8471486
#
#  1   2   3 
#396 680 163 
#[1] 4
#  Group.1    bet_perc   yft_perc  alb_perc
#1       1 0.032734083 0.04254883 0.9247171
#2       2 0.008987234 0.00744052 0.9835722
#3       3 0.115038781 0.18533184 0.6996294
#  Group.1    bet_perc   yft_perc  alb_perc
#1       1 0.032734083 0.04254883 0.9247171
#2       2 0.008987234 0.00744052 0.9835722
#3       3 0.115038781 0.18533184 0.6996294
#
#   1    2    3 
#1107 1487  165 
#

  print(table(clarax$clustering))

  junk4b <- table(clarax$clustering)
  junk4c <- cbind(a, clarax$clustering)                                                                    #appends cluster to dataframe, then names (in row below)
  names(junk4c) <- c("bet_perc","yft_perc","alb_perc","ret_year","tripest_id", "hk_bt_flt", "cluster")
  assign(paste("clust_R",reg,sep=""),junk4c)
  }

clust_R1$reg <- 1
clust_R2$reg <- 2
clust_R3$reg <- 3
clust_R4$reg <- 4
hclusts <- rbind(clust_R1,clust_R2,clust_R3,clust_R4)
save(hclusts,file="hclusts.RData")
load(file="hclusts.RData")

aggregate(hclusts[,1:3],by=list(hclusts$cluster,hclusts$reg),FUN=mean)
  
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
# Apparently the fit clustering takes too long. Might redo later since it's feasible on supercalc. 
# Better apply other species in the clustering? SWO? No maybe not. Just use tunas across the whole period. 
# This is where the trip clusters are identified, and later applied to the set-by-set data
tpx <- tp[,c("bet_perc","yft_perc","alb_perc","oth_totpc","trip_yr","tripest_id", "hk_bt_flt","reg")]
tpx <- na.omit(tpx)
nsp=3
clx <- claracl <- list()
ncl <- c(3,3,3,3)
tpx$clust <- tpx$claracl <- NA
for(rg in 1:4) {
  a <- tpx[tpx$reg==rg,]
  a[,1:nsp] <- scale(a[,1:nsp])
  clx[[rg]] <- kmeans(a[,1:nsp],centers=ncl[rg],iter.max = 40)
  claracl[[rg]] <- clara(a[,1:nsp],k=ncl[rg])
  tpx[tpx$reg==rg,]$clust <- clx[[rg]]$cluster
  tpx[tpx$reg==rg,]$claracl <- claracl[[rg]]$clustering
  }
  
# Apply clusters to sets
pl2$clust <- as.factor(tpx[match(paste(pl2$tripest_id,pl2$reg),paste(tpx$tripest_id,tpx$reg)),]$clust)
pl2$claracl <- as.factor(tpx[match(paste(pl2$tripest_id,pl2$reg),paste(tpx$tripest_id,tpx$reg)),]$claracl)
pl2$hclustcl <- as.factor(hclusts[match(paste(pl2$tripest_id,pl2$reg),paste(hclusts$tripest_id,hclusts$reg)),]$cluster)
save(pl2,file="pl2_clustered.RData")
load(file="pl2_clustered.RData")

length(unique(pl2$boat_id[!is.na(pl2$hclustcl)]))
length(unique(pl2$trip_id[!is.na(pl2$hclustcl)]))
length(unique(pl2$set_id[!is.na(pl2$hclustcl)]))

table(pl2$hclustcl,useNA="always")
a <- pl2[is.na(pl2$hclustcl)==F,]
with(a,aggregate(cbind(alb_perc,bet_perc,yft_perc) ~ hclustcl + reg,FUN=mean))
with(a,aggregate(cbind(alb_perc,bet_perc,yft_perc) ~ clust + reg,FUN=mean))
with(a,aggregate(cbind(alb_perc,bet_perc,yft_perc) ~ claracl + reg,FUN=mean))


#windows(width=20,height=18); par(mfrow=c(2,2))
# Plot clusters on a map. Mostly SDH code from here down. 
require(maps); require(maptools); require(mapdata)
windows(width=20,height=18); par(mfrow=c(2,2))
for(rg in 1:4) {
  a <- pl2[pl2$reg==rg & pl2$set_yr < 1980,]
  plot(1:10,1:10,type="n",xlab="Longitude",ylab="Latitude",xlim=range(a$lond),ylim=range(a$latd))
  for(cl in 1:ncl[rg]) {
    b <- a[a$clust==cl,]
    points(b$lond,b$latd,cex=0.5,col=cl)
    }
  lines(ezpac$lon, ezpac$lat)
  map(add=T)
  }
savePlot("map_clusters_kmeans_pre1980",type="png")
windows(width=20,height=18); par(mfrow=c(2,2))
for(rg in 1:4) {
  a <- pl2[pl2$reg==rg & pl2$set_yr < 1980,]
  plot(1:10,1:10,type="n",xlab="Longitude",ylab="Latitude",xlim=range(a$lond),ylim=range(a$latd))
  for(cl in 1:ncl[rg]) {
    b <- a[a$claracl==cl,]
    points(b$lond,b$latd,cex=0.5,col=cl)
    }
  lines(ezpac$lon, ezpac$lat)
  map(add=T)
  }
savePlot("map_clusters_clara_pre1980",type="png")
windows(width=20,height=18); par(mfrow=c(2,2))
for(rg in 1:4) {
  a <- pl2[pl2$reg==rg & pl2$set_yr < 1980,]
  plot(1:10,1:10,type="n",xlab="Longitude",ylab="Latitude",xlim=range(a$lond),ylim=range(a$latd))
  for(cl in 1:ncl[rg]) {
    b <- a[a$hclustcl==cl,]
    points(b$lond,b$latd,cex=0.5,col=cl)
    }
  lines(ezpac$lon, ezpac$lat)
  map(add=T)
  }
savePlot("map_clusters_hclust_pre1980",type="png")

# Make 2x2 maps
postyr <- 1980
preyr <- 1980
agglocs_hc <- with(pl2[pl2$set_yr <= preyr,],tapply(hook,list(2*floor(lond/2),2*floor(latd/2),hclustcl,reg),sum,na.rm=T))
agg <- with(pl2[pl2$set_yr <= preyr,],tapply(hook,list(2*floor(lond/2),2*floor(latd/2),reg),sum,na.rm=T))
#agglocs_hc <- tapply(pl2$hook,list(1*floor(pl2$lond/1),1*floor(pl2$latd/1),pl2$hclustcl,pl2$reg),sum,na.rm=T)
#agg <- tapply(pl2$hook,list(1*floor(pl2$lond/1),2*floor(pl2$latd/1),pl2$reg),sum,na.rm=T)
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
aggregate(cbind(alb_perc,bet_perc,yft_perc) ~ clust + reg,data=pl2,mean)
aggregate(cbind(alb_perc,bet_perc,yft_perc) ~ claracl + reg,data=pl2,mean)
aggregate(cbind(claracl==1,claracl==2,claracl==3) ~ reg,data=pl2,mean)
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
table(pl2$reg,pl2$clust)
table(pl2$reg,pl2$claracl)
  
windows(width=15,height=18); par(mfrow=c(2,2))
for(r in 1:4) { 
  plot(1:10,1:10,xlim=c(1960,2011),ylim=c(0,1),xlab="Year",ylab="",main="kmeans")
  a <- with(pl2[pl2$reg==r,],tapply(set_yr,list(set_yr,clust),length))
  a[is.na(a)] <- 0 
  for(cl in 1:3) {
    points(as.numeric(rownames(a)),a[,cl]/apply(a,1,sum),col=cl,pch=cl)
    }
  }
legend(1990,0.5,legend=paste("cluster",1:3),pch=1:3,col=1:3)
savePlot("clusters_by_yr_kmeans",type="png")

windows(width=15,height=18); par(mfrow=c(2,2))
for(r in 1:4) { 
  plot(1:10,1:10,xlim=c(1960,2011),ylim=c(0,1),xlab="Year",ylab="",main="clara")
  a <- with(pl2[pl2$reg==r,],tapply(set_yr,list(set_yr,claracl),length))
  a[is.na(a)] <- 0 
  for(cl in 1:3) {
    points(as.numeric(rownames(a)),a[,cl]/apply(a,1,sum),col=cl,pch=cl)
    }
  }
legend(1990,0.5,legend=paste("cluster",1:3),pch=1:3,col=1:3)
savePlot("clusters_by_yr_clara",type="png")

windows(width=15,height=18); par(mfrow=c(2,2))
for(r in 1:4) { 
  plot(1:10,1:10,xlim=c(1960,2011),ylim=c(0,1),xlab="Year",ylab="",main="hclust")
  a <- with(pl2[pl2$reg==r,],tapply(set_yr,list(set_yr,hclustcl),length))
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
pl2$boat_id <- as.factor(pl2$boat_id)
albpl <- pl2[paste(pl2$reg,pl2$claracl) %in% usecl_c,]
save(albpl,file="albacore_doubleclusters_clara.RData")
albpl <- pl2[paste(pl2$reg,pl2$hclustcl) %in% usecl_h,]
save(albpl,file="albacore_doubleclusters_hclust.RData")
albpl <- pl2[paste(pl2$reg,pl2$hclustcl) %in% usecl_k,]
save(albpl,file="albacore_doubleclusters_kmeans.RData")

ls()
rm(a,a1,a2,ax,b,b2,channel,cl,claracl,clx,d,dima,f,fit,groups,i,l_bycat,ncl,nsp,pl,pl5,r,reg,reg1,rg,tp,tp1,tpall,tpx,triptots,usecl,wss)

table(albpl$reg,albpl$flag_id)

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

windows(15,15); par(mfrow=c(2,2))
for(r in 1:4) {
  dd <- albpl_h[albpl_h$reg==r,]
  catch  <- tapply(dd$alb_n, dd$set_yrqtr, sum)
  effort <- tapply(dd$hook, dd$set_yrqtr, sum)
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
  catch  <- tapply(dd$alb_n, dd$set_yrqtr, sum)
  effort <- tapply(dd$hook, dd$set_yrqtr, sum)
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
  catch  <- tapply(dd$alb_n, dd$set_yrqtr, sum)
  effort <- tapply(dd$hook, dd$set_yrqtr, sum)
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
  a <- tapply(dd$boat_id,list(dd$boat_id,dd$set_yrqtr),length)
  a1 <- apply(is.na(a[,])=="FALSE",1,sum)
  a1 <- a1[a1>vessyrs]
  rb <- tapply(dd$set_yrqtr,dd$boat_id,max)
  rb1 <- rb[!is.na(rb) & rb>2000]
#  activeboats <- c(as.numeric(names(a1)), as.numeric(names(rb1)))
  activeboats <- as.numeric(names(a1)) # plenty of recent activity so no need for recent boats
  activeboats <- unique(activeboats)
  a <- dd[dd$boat_id %in% activeboats,]

  modtxt <- "alb_n ~ as.factor(set_yrqtr) + offset(log(hook))"
  if(length(unique(a$cl)) > 1) modtxt <- "alb_n ~ as.factor(set_yrqtr) + cl + offset(log(hook))"
  glm_nb1 <- glm.nb(modtxt, data=a)
  save(glm_nb1,file=paste("clara_glm_nb1_",r,".RData",sep="")); rm(glm_nb1); gc()
  modtxt <- "alb_n ~ as.factor(set_yrqtr) + latlong + as.factor(boat_id)+offset(log(hook))"
  if(length(unique(a$cl)) > 1) modtxt <- "alb_n ~ as.factor(set_yrqtr) + latlong + cl + as.factor(boat_id)+offset(log(hook))"
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

  catch  <- tapply(a$alb_n, a$set_yrqtr, sum)
  effort <- tapply(a$hook, a$set_yrqtr, sum)
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
  
  newdat <- with(a,expand.grid(set_yrqtr=sort(unique(a$set_yrqtr)),latlong=latlong[1],cl=cl[1],boat_id=boat_id[1],hook=median(hook)))
  a2 <- predict(glm_nb2,newdata=newdat,type="terms",se.fit=T)
  junk$coef2 <- exp(a2$fit[,1]) / mean(exp(a2$fit[,1]))
  junk$std_error2 <- a2$se.fit[,1]
  
  fname <- paste("spalb_clara", r,".dat",sep="")
  write.table(junk, fname, sep=c(" "), row.names=F)
  }
  

# consider activity filter to remove some vessels

lenun <- function(x) length(unique(x))
with(albpl,tapply(boat_id,reg,lenun))

atot <- albpl[albpl$reg==2,]

#for (vessyrs in c(10)) {
for (vessyrs in c(10)) { # for area 2
    # filter for active boats and post-2000
    a <- tapply(atot$boat_id,list(atot$boat_id,atot$set_yrqtr),length)
    a1 <- apply(is.na(a[,])=="FALSE",1,sum)
    a1 <- a1[a1>vessyrs]
    rb <- tapply(atot$set_yrqtr,atot$boat_id,max)
    rb1 <- rb[rb>2000]
    activeboats <- c(as.numeric(names(a1)), as.numeric(names(rb1)))
    activeboats <- unique(activeboats)
    a <- atot[atot$boat_id %in% activeboats,]
    rb <- tapply(a$set_yrqtr,a$boat_id,max)
    rb1 <- rb[rb>2000]
    recentboats <- as.numeric(names(rb1))
    a$era <- 0
    a$era[a$boat_id %in% recentboats] <- 1
}
with(a,tapply(boat_id,reg,lenun))

dim(a)
unique(a$boat_id) # reg 1 548
unique(a$boat_id) # reg 2 1006 filter of 40 
unique(a$boat_id) # reg 3 277 activity filter not used 
unique(a$boat_id) # reg 4 368 
table(a$flag_id)



setwd("D:/L/alb/2012/assessment/Data_preparation/CPUE/myStd")
for(r in 1:4) {
  a <- read.table(paste("cpue_reg",r,".dat",sep=""),header=T)
  a <- a[!is.na(a[,4]),c(1,4,6)]
  names(a) <- c("yrqtr","coef","std.error")
  write.table(a,file=paste("spalb_",r,".dat",sep=""),row.names=F)
}



# older stuff --------------------
pl$date <- as.Date(pl$end_date, format="%Y-%m-%d")
pl$ret_year <- as.numeric(format(pl$date, "%Y"))
pl$ret_month <- as.numeric(format(pl$date, "%m"))
pl$ret_day <- as.numeric(format(pl$date, "%d"))
table(pl$ret_year)
# 1961  1962  1963  1964  1965  1966  1967  1968  1969  1970  1971  1972  1973  1974  1975  1976  1977  1978  1979  1980  1981  1982  1983  1984  1985  1986  1987  1988  1989  1990  1991  1992  1993  1994  1995  1996  1997  1998  1999 
#  555    61  5919  5912  9556 19336 28271 22897 21132 24206 22816 21949 27777 26629 19146 14871 18804 15194 12874  8395 14634 16571  7622 11814 13697 14264 15519 13060  8502  5754  5831  5572  7273  6144  4713  5763  5632  4074  4912 
# 2000  2001  2002  2003  2004  2005  2006  2007  2008  2009  2010 
# 7827  2459  6220  5672  4571  3039  1873  2876   209   978    53 

# estimate no_sets per trip for GLM weighting
pl$tmp <- (rep(1, as.numeric(dim(pl)[1])))
junk2 <- aggregate(pl$tmp,list(pl$tripest_id),sum)       #i.e. take the trip id list, lots of 1's from the tmp column, and sum 'em to sum the number of rows (sets) in each trip 
pl$sets_per_trip <- junk2[match(pl[,2],junk2[,1]),2]     #then match columns tripest_id and 'Group.1' (tripest_id) from junk and return the sum by trip
pl$wt <- 1/sqrt(pl$sets_per_trip)
boxplot(pl$sets_per_trip~pl$set_yr) # high in 1960-61 because of two vessels with sets in both yrs, need to assign to start(trip_date) and graph    #GMP: also increased across time to 2005...
pl3 <- pl


