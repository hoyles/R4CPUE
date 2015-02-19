# Cluster function
require(rpart)
require(cluster)

#### Set up variables for trip-level clustering
opx <- op2[op2$lat > -5 & op2$lat < 10 & op2$reg %in% 1:6 & op2$trip_yr >= 1990 & !is.na(op2$tripid),]
opx79 <- op2[op2$lat > -5 & op2$lat < 10 & op2$reg %in% 1:6 & !is.na(op2$tripid),]
op3 <- op2[!is.na(op2$tripid),]

table(is.na(opx$reg))
table(is.na(op2$reg))


aggregate_for_cluster <- function(dat=opx) {
  tp <- aggregate(cbind(bet,yft,alb,swo,hbf,hooks,lat,lon) ~ tripid + newfishingcat + vessid + trip_yr + reg,data=dat,mean)
  a <- aggregate(bet ~ tripid + newfishingcat + vessid + trip_yr + reg,data=dat,length)
  tp <- cbind(tp,a[,6])
  dima <- dim(tp)[2]
  names(tp)[dima] <- "nset"
  tp$bet_cpue <- 1000*tp$bet/tp$hooks
  tp$yft_cpue <- 1000*tp$yft/tp$hooks
  tp$alb_cpue <- 1000*tp$alb/tp$hooks
  tp$swo_cpue <- 1000*tp$swo/tp$hooks
  tp$tot <- tp$alb + tp$yft + tp$bet + tp$swo
  tp <- tp[tp$tot > 0,]
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
  return(tp)
  }
  
prune_tree_and_plot <- function(tree) {
  minxerror <- tree$cptable[,4][which.min(tree$cptable[,4])] + tree$cptable[which.min(tree$cptable[,4]),5]
  if(sum(tree$cptable[,4] > minxerror)!=0) { 
    findxerror <- as.numeric(names(tree$cptable[,4][tree$cptable[,4] < minxerror][1]))
    treepr <- prune(tree,cp=tree$cptable[findxerror,1]+0.00001)
    printcp(treepr)
    if(dim(tree1$cptable)[2] > 1) {
      plot(treepr,main=paste("Region",rg),margin=0.1)
      text(treepr)
      } else plot(1:5,1:5,type="n",axes=F)
    }  else plot(1:5,1:5,type="n",axes=F)
  }

tp3 <- aggregate_for_cluster(dat=op3) 
tpx <- aggregate_for_cluster(dat=opx)

########################################################
# Regression trees on cpue by species & region
windows(width=18,height=18); par(mfrow=c(3,2),mar=c(2,2,2,2))
for(rg in 1:6) {
  a <- tp3[tp3$reg==rg,]
  tree1 <- rpart(bet_cpue ~ lat+lon+trip_yr+newfishingcat + alb_cpue + yft_cpue + swo_cpue + hooks + hbf, data=a)
  prune_tree_and_plot(tree1)
  }
savePlot("Rtrees_betcpue_by_reg_by_trip",type="png")
for(rg in 1:6) {
  a <- tp3[tp3$reg==rg,]
  tree1 <- rpart(alb_cpue ~ lat+lon+trip_yr+newfishingcat + bet_cpue + yft_cpue + swo_cpue + hooks + hbf, data=a)
  prune_tree_and_plot(tree1)
  }
savePlot("Rtrees_albcpue_by_reg_by_trip",type="png")
for(rg in 1:6) {
  a <- tp3[tp3$reg==rg,]
  tree1 <- rpart(yft_cpue ~ lat+lon+trip_yr+newfishingcat + bet_cpue + alb_cpue + swo_cpue + hooks + hbf, data=a)
  prune_tree_and_plot(tree1)
  }
savePlot("Rtrees_yftcpue_by_reg_by_trip",type="png")

# Regression trees by trip in core area on cpue by species & region
windows(width=10,height=10); par(mfrow=c(2,1),mar=c(2,2,2,2))
for(rg in 3:4) {
  a <- tpx[tpx$reg==rg,]
  tree1 <- rpart(bet_cpue ~ lat+lon+trip_yr+newfishingcat + alb_cpue + yft_cpue + swo_cpue + hooks + hbf, data=a)
  prune_tree_and_plot(tree1)
  }
savePlot("Rtrees_betcpue_by_reg_by_trip_corepost90",type="png")
for(rg in 3:4) {
  a <- tpx[tpx$reg==rg,]
  tree1 <- rpart(alb_cpue ~ lat+lon+trip_yr+newfishingcat + bet_cpue + yft_cpue + swo_cpue + hooks + hbf, data=a)
  prune_tree_and_plot(tree1)
  }
savePlot("Rtrees_albcpue_by_reg_by_trip_corepost90",type="png")
for(rg in 3:4) {
  a <- tpx[tpx$reg==rg,]
  tree1 <- rpart(yft_cpue ~ lat+lon+trip_yr+newfishingcat + bet_cpue + alb_cpue + swo_cpue + hooks + hbf, data=a)
  prune_tree_and_plot(tree1)
  }
savePlot("Rtrees_yftcpue_by_reg_by_trip_corepost90",type="png")

# Regression trees by trip in core area on proportion by species & region
windows(width=10,height=10); par(mfrow=c(2,1),mar=c(2,2,2,2))
for(rg in 3:4) {
  a <- tpx[tpx$reg==rg,]
  tree1 <- rpart(bet_troppc ~ lat+lon+trip_yr+newfishingcat + alb_cpue + swo_cpue + hooks + hbf, data=a)
  prune_tree_and_plot(tree1)
  }
savePlot("Rtrees_betprop_by_reg_by_trip_corepost90",type="png")
for(rg in 3:4) {
  a <- tpx[tpx$reg==rg,]
  tree1 <- rpart(yft_troppc ~ lat+lon+trip_yr+newfishingcat + alb_cpue + swo_cpue + hooks + hbf, data=a)
  prune_tree_and_plot(tree1)
  }
savePlot("Rtrees_yftprop_by_reg_by_trip_corepost90",type="png")

# Regression trees by set in core area on cpue
windows(width=10,height=10); par(mfrow=c(2,1),mar=c(2,2,2,2))
for(rg in 3:4) {
  a <- opx[opx$reg==rg,]
  a$alb_cpue <- a$alb/a$hooks; a$yft_cpue <- a$yft/a$hooks; a$bet_cpue <- a$bet/a$hooks; a$swo_cpue <- a$swo/a$hooks
  tree1 <- rpart(yft_cpue ~ lat+lon + trip_yr + newfishingcat + bet_cpue + alb_cpue + swo_cpue + hooks + hbf, data=a)
  prune_tree_and_plot(tree1)  
  }
savePlot("Rtrees_yftcpue_by_reg_by_set",type="png")
for(rg in 3:4) {
  a <- opx[opx$reg==rg,]
  a$alb_cpue <- a$alb/a$hooks; a$yft_cpue <- a$yft/a$hooks; a$bet_cpue <- a$bet/a$hooks; a$swo_cpue <- a$swo/a$hooks
  tree1 <- rpart(bet_cpue ~ lat+lon + trip_yr + newfishingcat + yft_cpue + alb_cpue + swo_cpue + hooks + hbf, data=a)
  prune_tree_and_plot(tree1)  
  }
savePlot("Rtrees_betcpue_by_reg_by_set",type="png")

# Regression trees by set in core area on species prop
windows(width=10,height=10); par(mfrow=c(2,1),mar=c(2,2,2,2))
for(rg in 3:4) {
  a <- opx[opx$reg==rg,]
  a$alb_cpue <- a$alb/a$hooks; a$swo_cpue <- a$swo/a$hooks
  a$yft_troppc <- a$yft/(a$yft + a$bet); a$bet_troppc <- a$bet/(a$yft + a$bet)
  tree1 <- rpart(yft_troppc ~ lat+lon + trip_yr + newfishingcat + alb_cpue + swo_cpue + hooks + hbf + mainline + branchline, data=a)
  prune_tree_and_plot(tree1)  
  }
savePlot("Rtrees_yftprop_by_reg_by_set_core",type="png")
for(rg in 3:4) {
  a <- opx[opx$reg==rg,]
  a$alb_cpue <- a$alb/a$hooks; a$swo_cpue <- a$swo/a$hooks
  a$yft_troppc <- a$yft/(a$yft + a$bet); a$bet_troppc <- a$bet/(a$yft + a$bet)
  tree1 <- rpart(bet_troppc ~ lat+lon + trip_yr + newfishingcat + alb_cpue + swo_cpue + hooks + hbf + mainline + branchline, data=a)
  prune_tree_and_plot(tree1)  
  }
savePlot("Rtrees_betprop_by_reg_by_set_core",type="png")

# Regression trees by set in core area on species prop 1976-1993
windows(width=10,height=10); par(mfrow=c(2,1),mar=c(2,2,2,2))
for(rg in 3:4) {
  a <- op2[op2$reg==rg & op2$op_yr >=1976 & op2$op_yr <= 1993 & op2$lat > -5 & op2$lat < 10,]
  a$alb_cpue <- a$alb/a$hooks; a$swo_cpue <- a$swo/a$hooks
  a$yft_troppc <- a$yft/(a$yft + a$bet); a$bet_troppc <- a$bet/(a$yft + a$bet)
  tree1 <- rpart(yft_troppc ~ lat+lon + trip_yr + newfishingcat + alb_cpue + swo_cpue + hooks + hbf + bait, data=a)
  prune_tree_and_plot(tree1)  
  }
savePlot("Rtrees_yftprop_by_reg_by_set_core_76_93bait",type="png")
for(rg in 3:4) {
  a <- op2[op2$reg==rg & op2$op_yr >=1976 & op2$op_yr <= 1993 & op2$lat > -5 & op2$lat < 10,]
  a$alb_cpue <- a$alb/a$hooks; a$swo_cpue <- a$swo/a$hooks
  a$yft_troppc <- a$yft/(a$yft + a$bet); a$bet_troppc <- a$bet/(a$yft + a$bet)
  tree1 <- rpart(bet_troppc ~ lat+lon + trip_yr + newfishingcat + alb_cpue + swo_cpue + hooks + hbf + bait, data=a)
  prune_tree_and_plot(tree1)  
  }
savePlot("Rtrees_betprop_by_reg_by_set_core_76_93bait",type="png")



########################################################
# Cluster analysis

# Look at numbers of clusters
# by region 
tpx <- tp3[,c("bet_totpc","yft_totpc","alb_totpc","swo_totpc","bet_troppc","yft_troppc","trip_yr","tripid","newfishingcat","hbf","reg")]
#triptots <- tp3[,c("bet","yft","alb","swo","trip_yr","tripid", "hbf","reg")]

tpx <- na.omit(tpx)
windows(width=20,height=18); par(mfrow=c(3,2),mar=c(2,4,4,2))
nsp=4
for(latlim in seq(10,0,-5)) {
  for(lonlim in seq(120,200,20)) {
    a <- tpx[tpx$lon >= lonlim & tpx$lon <(lonlim+20) & tpx$lat <= latlim & tpx$lat >(latlim-5),]
    a[,1:nsp] <- scale(a[,1:nsp])
    wss <- (nrow(a)-1)*sum(apply(a[,1:nsp],2,var))
    for (i in 2:15) wss[i] <- sum(kmeans(a[,1:nsp],centers=i,iter.max = 40)$withinss)
    plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares",main=paste("Region",rg))
    } }
savePlot(paste("cluster_num_v_SS_by_reg_",nsp,"_spp_fc",f,sep=""),type="png")

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

# for a single sample - run clustering
# 1. define the sample
# 2. set up the input format
# 3. run the cluster
# 4. output, check and save the results
# 5. apply the cluster back to the data


make_clusters <- function(setdat=opx,cldat=tpx,nsp=3,ncl=3,titx="") { 
  a <- na.omit(cldat)
  a <- scale(a[,1:nsp])
  d <- dist(a[,1:nsp], method = "euclidean") # distance matrix
  fit <- hclust(d, method="ward")
  plot(fit, labels = FALSE, hang=-1,  main = titx) # display dendogram  #looks like 3 (or 4)  
  groups <- cutree(fit, k=ncl) # cut tree into ncl clusters
  print(table(groups))
  rect.hclust(fit, k=ncl, border="red")
  clarax <- clara(a[,1:nsp],ncl)             #clustering based upon the percent of spp in total catch of tuna
  kmclust <- kmeans(a[,1:nsp],centers=ncl,iter.max = 40)
  setdat$kmclust <- kmclust$cluster[match(setdat$tripid,cldat$tripid)]
  setdat$claracl <- clarax$clustering[match(setdat$tripid,cldat$tripid)]
  setdat$hclustcl <- groups[match(setdat$tripid,cldat$tripid)]
  return(list(d=d,fit=fit,clarax=clarax,setdat=setdat))
  }

windows(18,18);par(mfrow=c(3,2))
for (r in 1:6) {
  a <- make_clusters(setdat=op3[op3$reg==r & !is.na(op3$tripid),],
      cldat=tp3[tp3$reg==r,c("bet_cpue","alb_cpue","yft_cpue","swo_cpue","tripid","lat","lon")],
      nsp=4,ncl=4,titx=paste("Region",r))
  assign(paste("clusters_4spp_cpue_R",r,sep=""),a)
  }
savePlot("cluster_4spp_cpue_hclust",type="png")

windows(18,18);par(mfrow=c(3,2))
for (r in 1:6) {
  a <- make_clusters(setdat=op3[op3$reg==r & !is.na(op3$tripid),],
      cldat=tp3[tp3$reg==r,c("bet_totpc","alb_totpc","yft_totpc","swo_totpc","tripid","lat","lon")],
      nsp=4,ncl=4,titx=paste("Region",r))
  assign(paste("clusters_4spp_spcomp_R",r,sep=""),a)
  }
savePlot("cluster_4spp_spcomp_hclust",type="png")


for(r in 1:6) { 
  a <- get(paste("clusters_4spp_cpue_R",r,sep=""))
  print(paste("Region",r))
  print(aggregate(cbind(bet/hooks,alb/hooks,yft/hooks,swo/hooks) ~ kmclust,data=a$setdat,FUN=mean))
  print(aggregate(cbind(bet/hooks,alb/hooks,yft/hooks,swo/hooks) ~ claracl,data=a$setdat,FUN=mean))
  print(aggregate(cbind(bet/hooks,alb/hooks,yft/hooks,swo/hooks) ~ hclustcl,data=a$setdat,FUN=mean))
  }

for(r in 1:6) { 
  a <- get(paste("clusters_4spp_spcomp_R",r,sep=""))
  a$setdat$bet_totpc <- with(a$setdat,bet/(bet+alb+yft+swo))
  a$setdat$yft_totpc <- with(a$setdat,yft/(bet+alb+yft+swo))
  a$setdat$alb_totpc <- with(a$setdat,alb/(bet+alb+yft+swo))
  a$setdat$swo_totpc <- with(a$setdat,swo/(bet+alb+yft+swo))
  print(paste("Region",r))
  print(aggregate(cbind(bet_totpc,alb_totpc,yft_totpc,swo_totpc) ~ kmclust,data=a$setdat,FUN=mean))
  print(aggregate(cbind(bet_totpc,alb_totpc,yft_totpc,swo_totpc) ~ claracl,data=a$setdat,FUN=mean))
  print(aggregate(cbind(bet_totpc,alb_totpc,yft_totpc,swo_totpc) ~ hclustcl,data=a$setdat,FUN=mean))
  }

#[1] "Region 1"
#  kmclust bet_totpc  alb_totpc  yft_totpc  swo_totpc
#1       1 0.1536289 0.06833749 0.04782119 0.73021246
#2       2 0.1851300 0.71447261 0.06030514 0.04009229
#3       3 0.5571895 0.29233113 0.05059215 0.09988726
#4       4 0.1445576 0.30283987 0.48052112 0.07208143
#  claracl bet_totpc  alb_totpc  yft_totpc  swo_totpc
#1       1 0.3385728 0.54153313 0.04572553 0.07416853
#2       2 0.1620779 0.07030975 0.05066386 0.71694852
#3       3 0.1104852 0.75441851 0.10929538 0.02580094
#4       4 0.6352546 0.22549886 0.05996981 0.07927675
#  hclustcl bet_totpc  alb_totpc  yft_totpc  swo_totpc
#1        1 0.4719572 0.37834385 0.05344760 0.09625136 BET/ALB
#2        2 0.1354442 0.77878282 0.06164946 0.02412351 ALB
#3        3 0.1629319 0.06293028 0.05101672 0.72312109 SWO
#4        4 0.1190933 0.34907554 0.50237503 0.02945614 YFT/ALB
#[1] "Region 2"
#  kmclust bet_totpc alb_totpc  yft_totpc  swo_totpc
#1       1 0.6567195 0.2040884 0.06594782 0.07324435
#2       2 0.1382498 0.1189409 0.04543195 0.69737733
#3       3 0.4570409 0.1882735 0.30256792 0.05211775
#4       4 0.2433549 0.6406329 0.04981352 0.06619861
#  claracl bet_totpc alb_totpc  yft_totpc  swo_totpc
#1       1 0.2051379 0.6914872 0.04287006 0.06050492
#2       2 0.4014663 0.4209617 0.09263449 0.08493748
#3       3 0.7450416 0.1216713 0.06996902 0.06331810
#4       4 0.1355003 0.1149961 0.04534030 0.70416330
#  hclustcl bet_totpc alb_totpc  yft_totpc  swo_totpc
#1        1 0.1954361 0.7135398 0.04431391 0.04671016  ALB
#2        2 0.6607419 0.1862159 0.10315734 0.04988481  BET
#3        3 0.3466165 0.5054282 0.06266327 0.08529206  BET/ALB
#4        4 0.1762068 0.1486632 0.04768199 0.62744807  SWO
#[1] "Region 3"
#  kmclust  bet_totpc  alb_totpc yft_totpc   swo_totpc
#1       1 0.39613578 0.02949170 0.5644550 0.009917549
#2       2 0.09277553 0.70798482 0.1943410 0.004898642
#3       3 0.64352576 0.02113801 0.3211552 0.014181035
#4       4 0.19271117 0.02376420 0.7766829 0.006841681
#  claracl  bet_totpc  alb_totpc yft_totpc   swo_totpc
#1       1 0.55207030 0.02419010 0.4125719 0.011167734
#2       2 0.31186934 0.02928189 0.6488266 0.010022153
#3       3 0.15805610 0.02222317 0.8135628 0.006157948
#4       4 0.09419902 0.70721399 0.1936928 0.004894148
#  hclustcl  bet_totpc  alb_totpc yft_totpc   swo_totpc
#1        1 0.42406791 0.01444774 0.5506179 0.010866409 BET/YFT
#2        2 0.20452259 0.03643860 0.7519490 0.007089786 YFT
#3        3 0.67787842 0.02317629 0.2857903 0.013155032 BET
#4        4 0.09200346 0.71869447 0.1847973 0.004504779 ALB
#[1] "Region 4"
#  kmclust bet_totpc  alb_totpc yft_totpc  swo_totpc
#1       1 0.3122835 0.47550627 0.2016801 0.01053017
#2       2 0.7189410 0.04047501 0.2233555 0.01722847
#3       3 0.2331442 0.01792848 0.2974835 0.45144385
#4       4 0.4169692 0.03201743 0.5395719 0.01144145
#  claracl bet_totpc  alb_totpc yft_totpc   swo_totpc
#1       1 0.3654746 0.02608587 0.5981864 0.010253178
#2       2 0.7832261 0.03055475 0.1678441 0.018375086
#3       3 0.5853077 0.05695255 0.3424368 0.015302973
#4       4 0.2343827 0.59826744 0.1591314 0.008218422
#  hclustcl bet_totpc  alb_totpc yft_totpc   swo_totpc
#1        1 0.2988523 0.01665656 0.6759633 0.008527862 YFT/BET
#2        2 0.7685379 0.04353982 0.1705537 0.017368554 BET
#3        3 0.5457279 0.04331624 0.3959485 0.015007316 BET/YFT
#4        4 0.2090365 0.61655172 0.1685203 0.005891574 ALB
#[1] "Region 5"
#  kmclust  bet_totpc alb_totpc yft_totpc  swo_totpc
#1       1 0.09245959 0.2561926 0.6221855 0.02916232
#2       2 0.07369510 0.6937159 0.1742617 0.05832723
#3       3 0.34216264 0.2129954 0.4097289 0.03511307
#4       4 0.10165820 0.4291136 0.2410797 0.22814856
#  claracl  bet_totpc alb_totpc  yft_totpc  swo_totpc
#1       1 0.09155720 0.2005735 0.68355239 0.02431686
#2       2 0.09208670 0.5329998 0.27939635 0.09551711
#3       3 0.33799688 0.1200607 0.50423843 0.03770403
#4       4 0.05955261 0.8052625 0.08723216 0.04795277
#  hclustcl  bet_totpc alb_totpc yft_totpc   swo_totpc
#1        1 0.10575738 0.5498555 0.2221687 0.122218434 ALB/YFT
#2        2 0.19107696 0.0265153 0.7756121 0.006795665 YFT
#3        3 0.04980333 0.7825375 0.1339286 0.033730618 ALB
#4        4 0.06815181 0.3621662 0.5527010 0.016980980 YFT/ALB
#[1] "Region 6"
#  kmclust  bet_totpc alb_totpc  yft_totpc swo_totpc
#1       1 0.16421188 0.3701292 0.38976075 0.0758982
#2       2 0.10164377 0.2760276 0.03049396 0.5918346
#3       3 0.30073900 0.5296506 0.04058380 0.1290266
#4       4 0.08215414 0.7438033 0.02400502 0.1500376
#  claracl  bet_totpc alb_totpc  yft_totpc  swo_totpc
#1       1 0.08388614 0.7366375 0.01794020 0.16153615
#2       2 0.04930062 0.1202030 0.02075805 0.80973833
#3       3 0.30765827 0.5175160 0.03877165 0.13605407
#4       4 0.15029794 0.4487137 0.31346368 0.08752467
#  hclustcl  bet_totpc alb_totpc  yft_totpc  swo_totpc
#1        1 0.06284170 0.7912712 0.01401011 0.13187704 ALB
#2        2 0.25396833 0.5829080 0.03181261 0.13131106 ALB/BET
#3        3 0.07476451 0.4660218 0.02225141 0.43696228 SWO/ALB
#4        4 0.16024543 0.4540296 0.30626019 0.07946474 ALB/YFT
#
for(r in 1:6) { 
  a <- get(paste("clusters_4spp_spcomp_R",r,sep=""))
  if(r==1) op_clust <- a$setdat else op_clust <- rbind(op_clust,a$setdat)
  }
save(op_clust,file="op_clust.RData")

op_clust <- op_clust[op_clust$bet+op_clust$yft+op_clust$alb+op_clust$swo > 0,]
write.csv(table(op_clust$hclustcl,op_clust$reg),file="hclust_n_reg.csv")
write.csv(with(op_clust,aggregate(cbind(bet,alb,swo,yft)/(bet+alb+swo+yft),by=list(hclustcl,reg),mean)),file="hclust_pc_reg.csv")
write.csv(table(op_clust$kmclust,op_clust$reg),file="kmeans_n_reg.csv")
write.csv(with(op_clust,aggregate(cbind(bet,alb,swo,yft)/(bet+alb+swo+yft),by=list(kmclust,reg),mean)),file="kmeans_pc_reg.csv")
write.csv(table(op_clust$claracl,op_clust$reg),file="clara_n_reg.csv")
write.csv(with(op_clust,aggregate(cbind(bet,alb,swo,yft)/(bet+alb+swo+yft),by=list(claracl,reg),mean)),file="clara_pc_reg.csv")

# histograms
for(r in 1:6) {
  windows(12,12);par(mfrow=c(6,6),mar=c(2,2,2,2),oma=c(0,0,2,0))
  for(y in 1978:2010) {
    hist(tp3[tp3$reg==r & tp3$trip_yr==y,]$bet_totpc,main=y,nclass=20,xlim=c(0,1))
    }
  title(paste("Region",r),outer=T,line=1)
  savePlot(paste("PropBET_per_trip_allR",r),type="png")
  }
for(r in 1:6) {
  windows(12,12);par(mfrow=c(6,6),mar=c(2,2,2,2),oma=c(0,0,2,0))
  for(y in 1978:2010) {
    hist(tp3[tp3$reg==r & tp3$trip_yr==y,]$alb_totpc,main=y,nclass=20,xlim=c(0,1))
    }
  title(paste("Region",r),outer=T,line=1)
  savePlot(paste("PropALB_per_trip_allR",r),type="png")
  }
for(r in 1:6) {
  windows(12,12);par(mfrow=c(6,6),mar=c(2,2,2,2),oma=c(0,0,2,0))
  for(y in 1978:2010) {
    hist(tp3[tp3$reg==r & tp3$trip_yr==y,]$yft_totpc,main=y,nclass=20,xlim=c(0,1))
    }
  title(paste("Region",r),outer=T,line=1)
  savePlot(paste("PropYFT_per_trip_allR",r),type="png")
  }
for(r in 1:6) {
  windows(12,12);par(mfrow=c(6,6),mar=c(2,2,2,2),oma=c(0,0,2,0))
  for(y in 1978:2010) {
    hist(tp3[tp3$reg==r & tp3$trip_yr==y,]$swo_totpc,main=y,nclass=20,xlim=c(0,1))
    }
  title(paste("Region",r),outer=T,line=1)
  savePlot(paste("PropSWO_per_trip_allR",r),type="png")
  }


  
#######################################################################################################
# cluster the core area
#######################################################################################################
# set up data
tpxa <- aggregate_for_cluster(dat=opx)
tpx <- tpxa[,c("bet_totpc","yft_totpc","alb_totpc","swo_totpc","bet_troppc","yft_troppc","bet_cpue","yft_cpue","alb_cpue","swo_cpue","trip_yr","tripid","newfishingcat","hbf","reg","lat","lon")]

tpx <- na.omit(tpx)
windows(width=20,height=15); par(mfrow=c(3,5),mar=c(2,4,4,2))
nsp=4
for(latlim in seq(10,0,-5)) {
  for(lonlim in seq(120,200,20)) {
    a <- tpx[tpx$lon >= lonlim & tpx$lon <(lonlim+20) & tpx$lat <= latlim & tpx$lat >(latlim-5),]
    a <- na.omit(a)
    a <- scale(a[,1:nsp])
    a[,1:nsp] <- scale(a[,1:nsp])
    if(nrow(a) > 10) {
      wss <- (nrow(a)-1)*sum(apply(a[,1:nsp],2,var))
      for (i in 2:15) wss[i] <- sum(kmeans(a[,1:nsp],centers=i,iter.max = 40)$withinss)
      plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares",main=paste(latlim,lonlim))
      }
    } }
savePlot(paste("cluster_num_v_SS_by_reg_",nsp,"_spp_fc",f,sep=""),type="png")

# identify catch rates by trip
for(rg in 3:4) {
  a <- tpx[tpx$reg==rg & tpx$newfishingcat==1,]
  windows(width=20,height=18); par(mfcol=c(3,2))
  boxplot(a$bet_totpc~a$trip_yr, xlab="Year", ylab="bet_perc",main="Offshore")
  boxplot(a$alb_totpc~a$trip_yr, xlab="Year", ylab="alb_perc")
  boxplot(a$yft_totpc~a$trip_yr, xlab="Year", ylab="yft_perc")
  title(paste("Region",rg,"core"),outer=T,line=-1)
  a <- tpx[tpx$reg==rg & tpx$newfishingcat==2,]
  boxplot(a$bet_totpc~a$trip_yr, xlab="Year", ylab="bet_perc",main="Distant water")
  boxplot(a$alb_totpc~a$trip_yr, xlab="Year", ylab="alb_perc")
  boxplot(a$yft_totpc~a$trip_yr, xlab="Year", ylab="yft_perc")
  title(paste("Region",rg,"core"),outer=T,line=-1)
  savePlot(paste("boxplot_sp_percent_reg_",rg,sep=""),type="png")
}

# for a single sample - run clustering
# 1. define the sample
# 2. set up the input format
# 3. run the cluster
# 4. output, check and save the results
# 5. apply the cluster back to the data


make_clusters <- function(setdat=opx,cldat=tpx,nsp=3,ncl=3,titx="") { 
  a <- na.omit(cldat)
  a <- scale(a[,1:nsp])
  d <- dist(a[,1:nsp], method = "euclidean") # distance matrix
  fit <- hclust(d, method="ward")
  plot(fit, labels = FALSE, hang=-1,  main = titx) # display dendogram  #looks like 3 (or 4)  
  groups <- cutree(fit, k=ncl) # cut tree into ncl clusters
  print(table(groups))
  rect.hclust(fit, k=ncl, border="red")
  clarax <- clara(a[,1:nsp],ncl)             #clustering based upon the percent of spp in total catch of tuna
  kmclust <- kmeans(a[,1:nsp],centers=ncl,iter.max = 40)
  setdat$kmclust <- kmclust$cluster[match(setdat$tripid,cldat$tripid)]
  setdat$claracl <- clarax$clustering[match(setdat$tripid,cldat$tripid)]
  setdat$hclustcl <- groups[match(setdat$tripid,cldat$tripid)]
  return(list(d=d,fit=fit,clarax=clarax,setdat=setdat))
  }

windows(18,18);par(mfrow=c(2,1))
for (r in 3:4) {
  a <- make_clusters(setdat=opx[opx$reg==r & !is.na(opx$tripid),],
      cldat=tpx[tpx$reg==r,c("bet_cpue","alb_cpue","yft_cpue","swo_cpue","tripid","lat","lon")],
      nsp=4,ncl=2,titx=paste("Region",r))
  assign(paste("clusters_4spp_cpue_core_R",r,sep=""),a)
  }
savePlot("cluster_4spp_cpue_core_hclust",type="png")

windows(18,18);par(mfrow=c(2,1))
for (r in 3:4) {
  a <- make_clusters(setdat=opx[opx$reg==r & !is.na(opx$tripid),],
      cldat=tpx[tpx$reg==r,c("bet_totpc","alb_totpc","yft_totpc","swo_totpc","tripid","lat","lon")],
      nsp=4,ncl=2,titx=paste("Region",r))
  assign(paste("clusters_4spp_spcomp_core_R",r,sep=""),a)
  }
savePlot("cluster_4spp_spcomp_core_hclust",type="png")

for(r in 3:4) { 
  a <- get(paste("clusters_4spp_cpue_core_R",r,sep=""))
  print(paste("Region",r))
  print(aggregate(cbind(bet/hooks,alb/hooks,yft/hooks,swo/hooks) ~ kmclust,data=a$setdat,FUN=mean))
  print(aggregate(cbind(bet/hooks,alb/hooks,yft/hooks,swo/hooks) ~ claracl,data=a$setdat,FUN=mean))
  print(aggregate(cbind(bet/hooks,alb/hooks,yft/hooks,swo/hooks) ~ hclustcl,data=a$setdat,FUN=mean))
  }

for(r in 3:4) { 
  a <- get(paste("clusters_4spp_spcomp_core_R",r,sep=""))
  a$setdat$bet_totpc <- with(a$setdat,bet/(bet+alb+yft+swo))
  a$setdat$yft_totpc <- with(a$setdat,yft/(bet+alb+yft+swo))
  a$setdat$alb_totpc <- with(a$setdat,alb/(bet+alb+yft+swo))
  a$setdat$swo_totpc <- with(a$setdat,swo/(bet+alb+yft+swo))
  print(paste("Region",r))
  print(aggregate(cbind(bet_totpc,alb_totpc,yft_totpc,swo_totpc) ~ kmclust,data=a$setdat,FUN=mean))
  print(aggregate(cbind(bet_totpc,alb_totpc,yft_totpc,swo_totpc) ~ claracl,data=a$setdat,FUN=mean))
  print(aggregate(cbind(bet_totpc,alb_totpc,yft_totpc,swo_totpc) ~ hclustcl,data=a$setdat,FUN=mean))
  }

for(r in 3:4) { 
  a <- get(paste("clusters_4spp_spcomp_core_R",r,sep=""))
  if(r==3) core_clust <- a$setdat else core_clust <- rbind(core_clust,a$setdat)
  }
save(core_clust,file="core_clust.RData")

a <- core_clust[,c("claracl","reg")]
core_clust$claracl[a$claracl==2 & a$reg==3] <- 1
core_clust$claracl[a$claracl==1 & a$reg==3] <- 2
a <- core_clust[,c("hclustcl","reg")]
core_clust$hclustcl[a$hclustcl==2 & a$reg==3] <- 1
core_clust$hclustcl[a$hclustcl==1 & a$reg==3] <- 2

core_clust <- core_clust[core_clust$bet+core_clust$yft+core_clust$alb+core_clust$swo > 0,]
write.csv(table(core_clust$hclustcl,core_clust$reg),file="hclust_n_core_reg.csv")
write.csv(with(core_clust,aggregate(cbind(bet,alb,swo,yft)/(bet+alb+swo+yft),by=list(hclustcl,reg),mean)),file="hclust_pc_core_reg.csv")
write.csv(table(core_clust$kmclust,core_clust$reg),file="kmeans_n_core_reg.csv")
write.csv(with(core_clust,aggregate(cbind(bet,alb,swo,yft)/(bet+alb+swo+yft),by=list(kmclust,reg),mean)),file="kmeans_pc_core_reg.csv")
write.csv(table(core_clust$claracl,core_clust$reg),file="clara_n_core_reg.csv")
write.csv(with(core_clust,aggregate(cbind(bet,alb,swo,yft)/(bet+alb+swo+yft),by=list(claracl,reg),mean)),file="clara_pc_core_reg.csv")

# plot sets per cluster per year
windows(height=10,width=8);par(mfrow=c(2,1))
a1 <- aggregate(bet ~ op_yr,data=core_clust[core_clust$hclustcl==1 & core_clust$reg==3,],FUN=length)
a2 <- aggregate(bet ~ op_yr,data=core_clust[core_clust$reg==3,],FUN=length)
plot(a1$op_yr,a1$bet/a2$bet,ylim=c(0,1),ylab="Proportion in YFT cluster",xlab="Year",main="Region 3")
a1 <- aggregate(bet ~ op_yr,data=core_clust[core_clust$hclustcl==1 & core_clust$reg==4,],FUN=length)
a2 <- aggregate(bet ~ op_yr,data=core_clust[core_clust$reg==4,],FUN=length)
plot(a1$op_yr,a1$bet/a2$bet,ylim=c(0,1),ylab="Proportion in YFT cluster",xlab="Year",main="Region 4")
savePlot("nsets_per_cluster_per_yr",type="png")

# plot maps of clusters per decade
a <- with(core_clust[core_clust$hclustcl==1,],tapply(hooks, list(factor(lon,levels=unique(core_clust$lon)),lat,eval(5*floor((op_yr)/5))),sum))
a2 <- with(core_clust,tapply(hooks, list(lon,lat,eval(5*floor((op_yr)/5))),sum))

a <- with(core_clust[core_clust$hclustcl==1,],aggregate(hooks, by=list(lon,lat,eval(5*floor((op_yr)/5))),sum))
a <- merge(a,with(core_clust,expand.grid(Group.1=unique(lon),Group.2=unique(lat),Group.3=seq(1990,2010,5))),all=T)
a2 <- with(core_clust,aggregate(hooks, by=list(lon,lat,eval(5*floor((op_yr)/5))),sum))
a2 <- merge(a2,with(core_clust,expand.grid(Group.1=unique(lon),Group.2=unique(lat),Group.3=seq(1990,2010,5))),all=T)
a[,4] <- a[,4]/a2[,4]
names(a)[1:4] <- c("lon","lat","decade","p")

windows(width=20,height=15);par(mfrow=c(2,2))
for(d in seq(1990,2005,5)) with(a,plot_catchmap(indat=a,vbl=p,latlim=c(-5,10),dcd=d))
savePlot("map_clusters_by_decade",type="png")

# plot nominal CPUE per cluster per year
a <- aggregate(cbind(bet/hooks,yft/hooks) ~ op_yr + reg + hclustcl,data=core_clust,FUN=mean)
names(a)[4:5] <- c("bet_cpue","yft_cpue")
a2 <- aggregate(cbind(bet/hooks,yft/hooks) ~ op_yr + reg,data=core_clust,FUN=mean)
names(a2)[3:4] <- c("bet_cpue","yft_cpue")

windows(10,10);par(mfrow=c(2,2))
for(r in 3:4) {
  x <- a[a$reg==r,]
  x2 <- a2[a2$reg==r,]
  with(x[x$hclustcl==2,],plot(op_yr,bet_cpue,col=2,type="l",ylim=c(0,0.02)))
  with(x[x$hclustcl==1,],lines(op_yr,bet_cpue,col=3))
  with(x2,lines(op_yr,bet_cpue,col=1))
  with(x[x$hclustcl==2,],plot(op_yr,yft_cpue,col=2,type="l",ylim=c(0,0.02)))
  with(x[x$hclustcl==1,],lines(op_yr,yft_cpue,col=3))
  with(x2,lines(op_yr,yft_cpue,col=1))
  }
legend("topright",legend=c("Not clustered","BET cluster","YFT cluster"),lty=1,col=c(1,2,3))
savePlot("Raw_CPUE_by_cluster_core",type="png")

opcore <- op2[op2$lat > -5 & op2$lat < 10 & op2$reg %in% 1:6 & !is.na(op2$tripid),]
tpcore <- aggregate_for_cluster(dat=opcore)
for(r in 3:4) {
  windows(12,12);par(mfrow=c(6,6),mar=c(2,2,2,2))
  for(y in 1978:2010) {
    hist(tpcore[tpcore$reg==r & tpcore$trip_yr==y,]$bet_totpc,main=y,nclass=20,xlim=c(0,1))
    }
  title(paste("Region",r),outer=T,line=1)
  savePlot(paste("PropBET_per_trip_core_R",r),type="png")
  }
for(r in 3:4) {
  windows(12,12);par(mfrow=c(6,6),mar=c(2,2,2,2))
  for(y in 1978:2010) {
    hist(tpcore[tpcore$reg==r & tpcore$trip_yr==y,]$yft_totpc,main=y,nclass=20,xlim=c(0,1))
    }
  title(paste("Region",r),outer=T,line=1)
  savePlot(paste("PropYFT_per_trip_core_R",r),type="png")
  }

