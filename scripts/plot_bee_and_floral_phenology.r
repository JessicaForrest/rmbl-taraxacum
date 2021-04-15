#############################################################
#File created 15 June 2020 by Jessica Forrest
#for plotting RMBL trapnest floral density data,
#using interpolated daily values.
#Goal: visually compare phenology of Taraxacum to that of
#(a) other Cichorieae and (b) other Asteraceae and (c) Aster-specialist bees
#across all years and sites (?). Or possibly just in 2018-2019.
#Much of this code adapted from Lydia Wong's scripts.
#Also, at end of script, add in plots of survival vs lay date to make composite plot (=Fig. 1 of manuscript)
#STILL NEED TO ADD LEGEND FOR PLANTS.
#############################################################

library(reshape2) #for melt and dcast functions

#Load floral data:
setwd("C:/Users/jforr/Dropbox/RMBL/RMBL trapnests")
fls<-read.delim("floral abd data 2013-2019.txt")
summary(fls)
#Use only 2016-19 because we only counted Taraxacum starting in 2016:
fls1619<-subset(fls,Year>=2016)
summary(fls1619)

#Load bee data:
setwd("C:/Users/jforr/Dropbox/RMBL/RMBL trapnests/Trapnests 2019/Charlotte's project")
bees<-read.delim("Taraxacum use 2013-2019 with fitness data CORRECTED.txt")
summary(bees)
bees1619<-subset(bees,Year>=2016)

#Organize floral data:
#Subset all asters:
levels(fls1619$Plant.taxon) 

#All Asteraceae: ----
all.asters <- subset(fls1619, (Plant.taxon == "Erigeron and Aster spp."|
                       Plant.taxon == "Helianthella quinquenervis"|
                       Plant.taxon == "Heliomeris multiflora"|
                       Plant.taxon == "Heterotheca villosa"|
                       Plant.taxon == "Hymenoxys hoopesii"|
                        Plant.taxon == "Other ligulate composites"|
                       Plant.taxon == "Pyrrocoma crocea"|
                       Plant.taxon == "Senecio, Arnica, Solidago spp."|
                         Plant.taxon == "Taraxacum officinale"|
                         Plant.taxon == "Wyethia amplexicaulis"))


#Interpolate daily floral densities for all Asteraceae:----
#create unique factor levels for each year-site-plant
all.asters$YSP <- as.factor(paste(all.asters$Year, all.asters$Site, all.asters$Plant.taxon, sep = "-"))

#empty dataframe to be filled in loop
f_ast1 <- data.frame()

#loop through each year-site-plant
for (i in levels(all.asters$YSP)){
  #subset out current year-site-plant
  focalYSP <- subset(all.asters,YSP == i)
  #create vector containing days of years from when flowering first occured to end of flowering
  day <- c((min(focalYSP$Day.of.year):max(focalYSP$Day.of.year)))
  #interpolation using spline; xout tells r the bounds of the interpolation
  focalSpline <- spline(focalYSP$Day.of.year, focalYSP$Density, xout = day)
  #create vector containing the results of the interpolation (second element in list)
  density <- focalSpline[[2]]
  #creating vectors that repeat the year, site and plant as many times as needed to 
  #build data frame
  focalYear <- focalYSP[1,1]
  focalSite <- focalYSP[1,2]
  focalPlant <- focalYSP[1,6]
  year <- as.vector(rep(focalYear,length(day)))
  site <- as.vector(rep(focalSite,length(day)))
  plant <- as.vector(rep(focalPlant,length(day)))
  
  #bind vectors together
  focalDf <- as.data.frame(cbind(year,site,plant,day,density))
  
  #bind dataframes together for all year-site-plants
  f_ast1 <- rbind(f_ast1, focalDf)
}

#make sure variables are the correct types
f_ast1$day <- as.numeric(as.character(f_ast1$day))
f_ast1$year <- as.numeric(as.character(f_ast1$year))
f_ast1$density <- as.numeric(as.character(f_ast1$density))

#some interpolations ended up negative
#if density was negative, round to 0, otherwise, round to 2 decimal places
f_ast1$density <- ifelse(f_ast1$density < 0, 0, round(f_ast1$density,2))

summary(f_ast1) #Checking.
f_ast1[701:730,]

#Now subset individual aster taxa or groups of taxa and calculate mean daily densities for each subset:----
#non-Cichorieae Asteraceae: ----
asters <- subset(f_ast1, (plant == "Erigeron and Aster spp."|
                            plant == "Helianthella quinquenervis"|
                            plant == "Heliomeris multiflora"|
                            plant == "Heterotheca villosa"|
                            plant == "Hymenoxys hoopesii"|
                            plant == "Pyrrocoma crocea"|
                            plant == "Senecio, Arnica, Solidago spp."|
                            plant == "Wyethia amplexicaulis"))

#Determining daily total aster density (sum across all spp per day at each site):
ast_melt <- melt(asters, id=c("year","site","plant","day"),measure.vars="density")
ast_cast <- dcast(ast_melt, year + site + day ~ variable, sum)

#create unique levels for each site-year  
ast_cast$ys <- as.factor(paste(ast_cast$year,ast_cast$site,sep="-"))

#Determine average floral density for each day of the year across all site-years
melt_ast <- melt(ast_cast, id = c("year", "site", "day", "ys"))
avast <- dcast(melt_ast, day ~ variable, mean) 

#non-Taraxacum Cichorieae: ----
cich <- subset(f_ast1, (plant == "Other ligulate composites"))
#Determining daily total Cichorieae density (not really necessary, but to match formatting of other asters):
cich_melt <- melt(cich, id=c("year","site","plant","day"),measure.vars="density")
cich_cast <- dcast(cich_melt, year + site + day ~ variable, sum)

#create unique levels for each site-year  
cich_cast$ys <- as.factor(paste(cich_cast$year,cich_cast$site,sep="-"))

#Determine average floral density for each day of the year across all site-years
melt_cich <- melt(cich_cast, id = c("year", "site", "day", "ys"))
avcich <- dcast(melt_cich, day ~ variable, mean) 

#Taraxacum: ----
tarax <- subset(f_ast1, (plant == "Taraxacum officinale"))
#Determining daily total Taraxacum density:
tarax_melt <- melt(tarax, id=c("year","site","plant","day"),measure.vars="density")
tarax_cast <- dcast(tarax_melt, year + site + day ~ variable, sum)

#create unique levels for each site-year  
tarax_cast$ys <- as.factor(paste(tarax_cast$year,tarax_cast$site,sep="-"))

#Determine average floral density for each day of the year across all site-years
melt_tarax <- melt(tarax_cast, id = c("year", "site", "day", "ys"))
avtarax <- dcast(melt_tarax, day ~ variable, mean) 

sites<-levels(ast_cast$site)


#Plot floral data, log-transformed:-----
par(mar=c(5,5,0,1))
with(ast_cast,plot(log(density)~day,type = "n",ylim=c(-4,8),xlim=c(152,225),cex.lab=0.8,ylab="",
      xlab="Date",axes=F)) #All the densities >36 are from MM. Use ylim=c(0,12) for log transform; ylim=(-5,5) for log-transform
axis(1,at=c(152,182,213,244),labels=c("Jun 1","Jul 1","Aug 1","Sep 1"),cex.axis=0.8);axis(2,seq(-4,4,2),cex.axis=0.8)
mtext(expression(log(paste("Floral density (capitula per  ", m^2,")"))),side=2,line=2,adj=0.2,cex=0.8)
for (j in c(2016:2019)){
  for (i in c(1:6)){
  with(subset(subset(ast_cast, year==j), site==sites[i]),
              lines(log(density)~day,col=adjustcolor("purple3",alpha=0.5),lwd=1))
  }
}
for (j in c(2018:2019)){
  for (i in c(1:6)){
    with(subset(subset(cich_cast, year==j), site==sites[i]),
         lines(log(density)~day,col=adjustcolor("red3",alpha=0.5),lwd=1))
  }
}
for (j in c(2016:2019)){
  for (i in c(1:6)){
    with(subset(subset(tarax_cast, year==j), site==sites[i]),
         lines(log(density)~day,col=adjustcolor("orange",alpha=0.5),lwd=1))
  }
}
with(avast,lines(log(density)~day,col=adjustcolor("purple3",alpha=0.8),lwd=3))
with(avcich,lines(log(density)~day,col=adjustcolor("red3",alpha=0.8),lwd=3))
with(avtarax,lines(log(density)~day,col=adjustcolor("orange",alpha=0.8),lwd=3))

#Plot bee data:----
#For each bee species, overlay horizontal boxplot of first nesting dates:
with(subset(bees1619, Species=="Osmia coloradensis"),
     boxplot(Day.of.year.first.cell,add=T,horizontal=T,boxwex=1.5,at=6,col=grey(0.7),bty="n",axes=F))
text(155,6.5,"O. coloradensis",font=3,adj=0,cex=0.8)
with(subset(bees1619, Species=="Osmia montana"),
     boxplot(Day.of.year.first.cell,add=T,horizontal=T,boxwex=1.5,at=5,col=grey(0.4),bty="n",axes=F))
text(160,5.2,"O. montana",font=3,adj=0,cex=0.8)
with(subset(bees1619, Species=="Osmia subaustralis"),
     boxplot(Day.of.year.first.cell,add=T,horizontal=T,boxwex=1.5,at=4,col=grey(0.4),bty="n",axes=F))
text(150,4.2,"O. subaustralis",font=3,adj=0,cex=0.8)

#Plot floral data, square-root-transformed:-----
par(mar=c(5,5,0,1))
with(ast_cast,plot(sqrt(density)~day,type = "n",ylim=c(0,7),xlim=c(152,225),cex.lab=0.8,ylab="",
                   xlab="date",axes=F)) 
axis(1,at=c(152,167,182,197,213,228),labels=c("Jun 1","Jun 16", "Jul 1","Jul 16","Aug 1","Aug 16"),cex.axis=0.8);
axis(2,seq(0,5,2),cex.axis=0.8)
mtext(expression(sqrt(paste("Floral density (capitula per  ", m^2,")"))),side=2,line=2,adj=0.1,cex=0.8)
  for (j in c(2016:2019)){
  for (i in c(1:6)){
    with(subset(subset(ast_cast, year==j), site==sites[i]),
         lines(sqrt(density)~day,col=adjustcolor("purple3",alpha=0.5),lwd=1))
  }
}
for (j in c(2018:2019)){
  for (i in c(1:6)){
    with(subset(subset(cich_cast, year==j), site==sites[i]),
         lines(sqrt(density)~day,col=adjustcolor("red3",alpha=0.5),lwd=1))
  }
}
for (j in c(2016:2019)){
  for (i in c(1:6)){
    with(subset(subset(tarax_cast, year==j), site==sites[i]),
         lines(sqrt(density)~day,col=adjustcolor("orange",alpha=0.5),lwd=1))
  }
}
with(avast,lines(sqrt(density)~day,col=adjustcolor("purple3",alpha=0.8),lwd=3))
with(avcich,lines(sqrt(density)~day,col=adjustcolor("red3",alpha=0.8),lwd=3))
with(avtarax,lines(sqrt(density)~day,col=adjustcolor("orange",alpha=0.8),lwd=3))

#Plot bee data:----
#For each bee species, overlay horizontal boxplot of first nesting dates:
with(subset(bees1619, Species=="Osmia coloradensis"),
     boxplot(Day.of.year.first.cell,add=T,horizontal=T,boxwex=0.7,at=6.2,col=grey(0.7),bty="n",axes=F))
text(155,6.4,"O. coloradensis",font=3,adj=0,cex=0.8)
with(subset(bees1619, Species=="Osmia montana"),
     boxplot(Day.of.year.first.cell,add=T,horizontal=T,boxwex=0.7,at=5.6,col=grey(0.4),bty="n",axes=F))
text(162,5.7,"O. montana",font=3,adj=0,cex=0.8)
with(subset(bees1619, Species=="Osmia subaustralis"),
     boxplot(Day.of.year.first.cell,add=T,horizontal=T,boxwex=0.7,at=5,col=grey(0.4),bty="n",axes=F))
text(153,5.1,"O. subaustralis",font=3,adj=0,cex=0.8)

######
#Thoughts on this (from email to Charlotte and Paul, 15 June 2020):
#We could plot an average across sites (and/or years) instead of many individual lines. 
#Bee data represent the distributions of first nesting dates for all individuals (i.e. earliest date of brood cell construction), 
#across all sites and years. 
#(So, if there were more bees at one site in a given year, that site/year is more influential in shaping the bee distribution.) 
#There are plenty of other ways we could show these data (e.g. distributions of median first nesting dates for each site-year), 
#but this seemed relatively easy and relevant. 

######
#PLOTTING FIG. 1 OF MANUSCRIPT:
####
#Merge this plot with overwinter survival graphs (see plot_survival_vs_lay_date.R) 
figure.1<-layout(matrix(c(1,1,1,2,3,4),2,3,byrow=T),widths=c(0.43,0.38,0.19),heights=c(2,2))
#layout.show(figure.1)

#Panel A (phenology):
par(mar=c(3,4.5,0,3),cex=0.8)
with(ast_cast,plot(log(density)~day,type = "n",ylim=c(-4,8),xlim=c(152,225),cex.lab=0.8,ylab="",
                   xlab="Date",axes=F)) #All the densities >36 are from MM. Use ylim=c(0,12) for log transform; ylim=(-5,5) for log-transform
axis(1,at=c(152,182,213,244),labels=c("Jun 1","Jul 1","Aug 1","Sep 1"),cex.axis=0.8);axis(2,seq(-4,4,2),cex.axis=0.8)
mtext(expression(log(paste("capitula per ", m^2))),side=2,line=2.5,adj=0.2,cex=0.8)
for (j in c(2016:2019)){
  for (i in c(1:6)){
    with(subset(subset(ast_cast, year==j), site==sites[i]),
         lines(log(density)~day,col=adjustcolor("purple3",alpha=0.5),lwd=1))
  }
}
for (j in c(2018:2019)){
  for (i in c(1:6)){
    with(subset(subset(cich_cast, year==j), site==sites[i]),
         lines(log(density)~day,col=adjustcolor("red3",alpha=0.5),lwd=1))
  }
}
for (j in c(2016:2019)){
  for (i in c(1:6)){
    with(subset(subset(tarax_cast, year==j), site==sites[i]),
         lines(log(density)~day,col=adjustcolor("orange",alpha=0.5),lwd=1))
  }
}
with(avast,lines(log(density)~day,col=adjustcolor("purple3",alpha=0.8),lwd=3))
with(avcich,lines(log(density)~day,col=adjustcolor("red3",alpha=0.8),lwd=3))
with(avtarax,lines(log(density)~day,col=adjustcolor("orange",alpha=0.8),lwd=3))

#For each bee species, overlay horizontal boxplot of first nesting dates:
with(subset(bees1619, Species=="Osmia coloradensis"),
     boxplot(Day.of.year.first.cell,add=T,horizontal=T,boxwex=1.5,at=6,col=grey(0.7),bty="n",axes=F))
text(155,6.5,"O. coloradensis",font=3,adj=0,cex=0.8)
with(subset(bees1619, Species=="Osmia montana"),
     boxplot(Day.of.year.first.cell,add=T,horizontal=T,boxwex=1.5,at=5,col=grey(0.4),bty="n",axes=F))
text(170,5.2,"O. montana",font=3,adj=1,cex=0.8)
with(subset(bees1619, Species=="Osmia subaustralis"),
     boxplot(Day.of.year.first.cell,add=T,horizontal=T,boxwex=1.5,at=4,col=grey(0.4),bty="n",axes=F))
text(163,4.2,"O. subaustralis",font=3,adj=1,cex=0.8)

#Panel label and plant labels:
text(150,7,"(a)")
text(155,1.5,"Taraxacum",font=3,col="orange")
text(210,-2.5,"other Cichorieae",font=1,col="red3",adj=0)
text(180,2,"other Asteraceae",font=1,col="purple3")

