#Script by JF Oct. 2019
#Goal: to test associations among aster-specialist Osmia phenology, fitness,
#and Taraxacum pollen usage. Data from RMBL, 2013-2019.
#Updated 22 April 2020. 
#Updated 1 July 2020 with corrected data file (using estimated dates each egg was laid, rather than progress recorded in occupant data)
#Updated again 10 July with corrected data file (added 10 previously missing bees). Doesn't qualitatively change anything.
#See bottom for re-analysis excluding VB.

########LOAD PACKAGES AND FILES#####
#load packages#
library(lmerTest)
library(car) #for vif
library(ggplot2)
#library(RColorBrewer)

#set working directory#
setwd("C:/Users/jforrest/Dropbox/RMBL/RMBL trapnests/Trapnests 2019/Charlotte's project")
setwd("C:/Users/jforr/Dropbox/RMBL/RMBL trapnests/Trapnests 2019/Charlotte's project")

#load dataset #
data<-read.delim("Taraxacum use 2013-2019 with fitness data CORRECTED.txt")
summary(data)
data$Proportion.Taraxacum <- data$Sum.of.count.Taraxacum/(data$Sum.of.count.Taraxacum+data$Sum.of.count.Asteraceae+data$Sum.of.count.other)

#####CHECK DATA####
attach(data)
hist(Proportion.Taraxacum) #mostly zeroes
hist(asin(sqrt(Proportion.Taraxacum))) #similar.
hist(Day.of.year.first.cell) #normal
hist(Offspring.produced.Aug1) #right-skewed
hist(sqrt(Offspring.produced.Aug1)) #better; still a bit skewed. Can't use log-transform b/c several zeroes.
hist(Offspring.produced.Aug1^0.25) #very peaked.
hist(Total.offspring) #very right-skewed
hist(log(Total.offspring)) #seems normal (no zeroes so no problem with log transform), but actually not, based on qqplots below.
hist(sqrt(Total.offspring)) #some right-skew.
plot(Proportion.Taraxacum ~ Species)
anova(lm(Proportion.Taraxacum ~ Species))
#Significant differences among species, so including Taraxacum usage and species in the same model is dubious.
fisher.test(table(Taraxacum.Y.N,Species))
detach(data)
#Same conclusion.
#THEREFORE, NEED TO CHECK VIFS FOR MODELS BELOW.

#####Descriptive data by species:####
summary(subset(data,Species=="Osmia montana")) #mean date = 193.3
summary(subset(data,Species=="Osmia subaustralis")) #mean date = 186.2
summary(subset(data,Species=="Osmia coloradensis")) #mean date = 191.3

#How many individuals used Taraxacum, and how many used 100% Taraxacum?
table(data$Taraxacum.Y.N)
table(data$Proportion.Taraxacum,data$Species)

#Plot:
species.plot<-ggplot(data,aes(x=Species,y=Proportion.Taraxacum))+
  scale_color_manual(values=c(adjustcolor(grey(0.7),alpha=0.5),adjustcolor(grey(0.4),alpha=0.5),adjustcolor(grey(0.4),alpha=0.5)), 
                     name="Species",
                     breaks=c("Osmia coloradensis", "Osmia montana","Osmia subaustralis"))+
  geom_point(aes(col=Species,cex=Brood.cells.scored/2),position=position_jitter(width=0.1,height=0)) +
  ylab(expression(paste("Proportion ",italic("Taraxacum")))) +
  xlab("")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(colour="black"))+
  theme(axis.text.x=element_text(face="italic"))+ 
  theme(legend.position="right") +
  guides(color = FALSE)+
  theme(legend.position = "top",
        legend.box = "horizontal")
#To move legend inside the plot:
#species.plot +   theme(legend.position = c(0.6,0.9),legend.direction = "horizontal")
species.plot + scale_size_continuous(name = "Brood cells sampled:",labels=c(1,12),breaks=c(0.5,6))
#With gratuitous violins overlaid:
#species.plot+geom_violin(aes(Species),alpha=0.5)

#####Q1. Do bees that use dandelion pollen begin nesting earlier?#####
#Q1a - treating dandelion usage as a proportion. THIS IS THE PREFERRED APPROACH (MORE INFORMATIVE):+
dand.phen.model<-lmer(Day.of.year.first.cell ~ Proportion.Taraxacum * Species + (1|Site.code) + (1|Year))
summary(dand.phen.model)
anova(dand.phen.model)
vif(dand.phen.model)
plot(dand.phen.model) #heterogeneity of variance
qqnorm(residuals(dand.phen.model))
#Strong effect of Taraxacum usage, no species*T interaction. Simplifying:
dand.phen.model.1<-lmer(Day.of.year.first.cell ~ Proportion.Taraxacum + Species + (1|Site.code) + (1|Year),data=data)
summary(dand.phen.model.1)
anova(dand.phen.model.1) 
vif(dand.phen.model.1)
plot(dand.phen.model.1) #heterogeneity of variance
qqnorm(residuals(dand.phen.model.1))
#Significant effect of species, with subaustralis significantly earlier than coloradensis. Strong effect of Taraxacum usage.
#Highly significant negative effect of Taraxacum usage.
#Re-order to check whether Cephalosmia species differ from one another:
dand.phen.model.1<-lmer(Day.of.year.first.cell ~ Proportion.Taraxacum + relevel(Species,"Osmia montana") + (1|Site.code) + (1|Year))
summary(dand.phen.model.1)
#Yes, O. subaustralis is earlier than both other species, by 7-9 days (after accounting for effect of Taraxacum usage).
#Reversing direction of causation (*should include weights, not done):
dand.phen.model.1<-glmer(Proportion.Taraxacum ~ Day.of.year.first.cell + Species + (1|Site.code) + (1|Year),family = binomial, data=data)
summary(dand.phen.model.1)
anova(dand.phen.model.1) 
#Singular fit.

#Q1b - treating dandelion usage as binary:
#Full model, treating DOY as response, Taraxacum use as binary, and allowing for effect of 
#Taraxacum to vary by species:
dand.phen.model<-lmer(Day.of.year.first.cell ~ Taraxacum.Y.N * Species + (1|Site.code) + (1|Year),data=data)
summary(dand.phen.model)
anova(dand.phen.model)
#Signficant effect of Taraxacum usage, no species*T interaction. Simplifying:
dand.phen.model.1<-lmer(Day.of.year.first.cell ~ Taraxacum.Y.N + Species + (1|Site.code) + (1|Year),data=data)
summary(dand.phen.model.1)
anova(dand.phen.model.1)
vif(dand.phen.model.1)
plot(dand.phen.model.1)
qqnorm(residuals(dand.phen.model.1)) #No problems.
#Significant effect of species, with subaustralis significantly earlier than coloradensis
#Highly significant negative effect of Taraxacum usage.
#Re-order to check whether Cephalosmia species differ from one another:
dand.phen.model.1<-lmer(Day.of.year.first.cell ~ Taraxacum.Y.N + relevel(Species,"Osmia montana") + (1|Site.code) + (1|Year),data=data)
summary(dand.phen.model.1)
#Yes, O. subaustralis is earlier than both other species, by 6-7 days.

####Plotting proportion Taraxacum vs. date (for PPT presentation, not manuscript):####
tarax.vs.date<-ggplot(data,aes(x=Day.of.year.first.cell, y=Proportion.Taraxacum))+
  scale_color_manual(values=rep(adjustcolor("orange",alpha=0.5),3), 
                     name="Species")+
  stat_smooth(method = "glm", method.args = list(family = "binomial"), colour="black",size=1,fill=grey(0.8))+
  geom_point(aes(col=Species),size=2,position=position_jitter(width=0.5,height=0)) +
  #ylab(expression(paste("Proportion ",italic("Taraxacum")))) +
  ylab("Proportion dandelion pollen in nests") +
  xlab("Date of first nesting")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(colour="black"))+
  #theme(legend.text=element_text(face="italic")) +
  scale_x_continuous(breaks = c(167,182,197,213,228),labels = c("Jun 16","Jul 1", "Jul 16","Aug 1","Aug 16"))+
  theme(legend.position="None") 
tarax.vs.date<-tarax.vs.date + scale_y_continuous(breaks = c(0,1),labels = c("0","1"))

tiff("taraxacum.vs.date.laid.tif", units="in", width=5, height=4, res=288,pointsize=12)
tarax.vs.date
dev.off()

####Plotting date vs. Taraxacum binary (Fig. 3A):#### 
cols<-c(adjustcolor("purple3",alpha=0.8),adjustcolor("orange",alpha=0.8))
date.vs.tarax<-ggplot(data,aes(x=Species,y=Day.of.year.first.cell,fill=Taraxacum.Y.N))+
  scale_fill_manual(values=cols, 
                     name=expression(paste(italic("Taraxacum")," in pollen provision:")),
                     breaks=c("N","Y"),
                    labels=c("No","Yes"))+
  scale_x_discrete(breaks=c("Osmia coloradensis", "Osmia montana","Osmia subaustralis"),labels=c("O. coloradensis", "O. montana","O. subaustralis"))  +
  geom_boxplot()+
  ylab("Day of year of first nesting")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(colour="black"))+
  theme(axis.text.x = element_text(face="italic"))  + 
#  theme(legend.position = "top",
#        legend.box = "horizontal")+
  scale_y_continuous(breaks = c(152, 182, 213),labels = c("Jun 1", "Jul 1", "Aug 1"))
date.vs.tarax +   theme(legend.position = "left")

###Plotting date vs. proportion Taraxacum:####
par(mar=c(5,5,1,3))
plot(Day.of.year.first.cell ~ Proportion.Taraxacum, type ="n",
     ylab="Day of year first cell constructed",xlab="Proportion dandelion pollen",data=data)
axis(4,at=c(152,182,213),labels=c("Jun 1","Jul 1","Aug 1"))
plotting.model<-lm(Day.of.year.first.cell~Proportion.Taraxacum,data=data)
#Plot fitted curve & confidence intervals:
df <- data.frame(Proportion.Taraxacum=seq(0,1,0.01))
CI <- predict(plotting.model, newdata=df, interval="confidence")
CI <- as.data.frame(CI) # Coerce the matrix to a dataframe, so we can access the column using the $ operator.
polygon.x <- c(df$Proportion.Taraxacum, rev(df$Proportion.Taraxacum))
polygon.y <- c(CI$lwr, rev(CI$upr))
polygon(x=polygon.x, y=polygon.y, col=adjustcolor("grey",alpha=0.5),border=NA)
lines(x=df$Proportion.Taraxacum, y=CI$fit, lwd=2)
#Adding datapoints by species:
with(subset(data,Species=="Osmia coloradensis"),
     points(Proportion.Taraxacum, Day.of.year.first.cell,
            pch=21,col=adjustcolor("orange",alpha=0.6),bg=adjustcolor("orange",alpha=0.6)))
with(subset(data,Species=="Osmia montana"),
     points(Proportion.Taraxacum, Day.of.year.first.cell,
            pch=21,col=adjustcolor("forestgreen",alpha=0.6),bg=adjustcolor("forestgreen",alpha=0.6)))
with(subset(data,Species=="Osmia subaustralis"),
     points(Proportion.Taraxacum, Day.of.year.first.cell,
            pch=21,col=adjustcolor("navy",alpha=0.6),bg=adjustcolor("navy",alpha=0.6)))

#####Q2. Does nesting earlier predict number of brood cells constructed before Aug. 1?####
#Full model, treating Offspring produced before Aug. 1 as response, 
fitness.phen.model<-lmer(Offspring.produced.Aug1^0.5 ~ Day.of.year.first.cell * Species + (1|Site.code) + (1|Year),data=data)
summary(fitness.phen.model)
anova(fitness.phen.model)
vif(fitness.phen.model)
#Major collinearity problem!
plot(fitness.phen.model)
qqnorm(residuals(fitness.phen.model)) #There is a definite outlier here: Jane (obs. 163)
#Re-running without Jane:
fitness.phen.model<-lmer(Offspring.produced.Aug1^0.5 ~ Day.of.year.first.cell * Species + (1|Site.code) + (1|Year),data=data[-163,])
#No more outlier; results qualitatively unchanged. But collinearity remains, plus singular fit.
#Strong effect of day of year, marginal species*date interaction. Simplifying:
fitness.phen.model.1<-lmer(Offspring.produced.Aug1^0.5 ~ Day.of.year.first.cell + Species + (1|Site.code) + (1|Year),data=data)
summary(fitness.phen.model.1) #almost no variance for site.
anova(fitness.phen.model.1) 
vif(fitness.phen.model.1)
plot(fitness.phen.model.1)
#Removing the interaction solves the collinearity problem, but Jane remains an outlier.
fitness.phen.model.1<-lmer(Offspring.produced.Aug1^0.5 ~ Day.of.year.first.cell + Species + (1|Site.code) + (1|Year),data=data[-163,])
#Results are the same.
summary(fitness.phen.model.1)
anova(fitness.phen.model.1) 
vif(fitness.phen.model.1)
plot(fitness.phen.model.1)
qqnorm(residuals(fitness.phen.model.1)) 
#Predicted number of offspring at mean first nesting date vs. 10 days earlier:
#Jane included:
x.mean<-with(data,mean(Day.of.year.first.cell))
x.vals<-c(x.mean,x.mean-10)
y.vals<-(12.066075-0.051113*x.vals)^2
(y.vals[2]-y.vals[1])/y.vals[1] #This is a measure of effect size (proportional increase in number of offspring from a 9 day advance in first nesting date)
y.vals
##Jane excluded:
x.mean<-with(data[-163,],mean(Day.of.year.first.cell)) 
x.vals<-c(x.mean,x.mean-10)
y.vals<-(11.465459-0.0483*x.vals)^2
(y.vals[2]-y.vals[1])/y.vals[1] #This is a measure of effect size (proportional increase in number of offspring from a 9 day advance in first nesting date)
y.vals

#####Q2b. Does nesting earlier predict number of brood cells overall?####
#Full model, treating total offspring produced as response, 
fitness.phen.model.2<-lmer(sqrt(Total.offspring) ~ Day.of.year.first.cell * Species + (1|Site.code) + (1|Year),data=data)
summary(fitness.phen.model.2)
anova(fitness.phen.model.2)
vif(fitness.phen.model.2)
#Major collinearity problem!
#Strong effect of day of year, no species*date interaction. Simplifying:
fitness.phen.model.3<-lmer(sqrt(Total.offspring) ~ Day.of.year.first.cell + Species + (1|Site.code) + (1|Year),data=data) #with Jane
fitness.phen.model.3<-lmer(sqrt(Total.offspring) ~ Day.of.year.first.cell + Species + (1|Site.code) + (1|Year),data=data[-163,]) #without Jane
summary(fitness.phen.model.3)
anova(fitness.phen.model.3) 
vif(fitness.phen.model.3)
plot(fitness.phen.model.3)
qqnorm(residuals(fitness.phen.model.3)) 
#Removing the interaction solves the collinearity problem. Some non-linearity with log-transform, though.
#Non-linearity goes away with sqrt transform, but then Jane is an outlier. Removing Jane solves that problem.
#Result is the same regardless.
#Jane included:
x.mean<-with(data,mean(Day.of.year.first.cell))
x.vals<-c(x.mean,x.mean-10)
y.vals<-(7.835706-0.026599*x.vals)^2
(y.vals[2]-y.vals[1])/y.vals[1] #This is a measure of effect size (proportional increase in number of offspring from a 9 day advance in first nesting date)
y.vals
#Jane excluded:
x.mean<-with(data[-163,],mean(Day.of.year.first.cell)) 
x.vals<-c(x.mean,x.mean-10)
y.vals<-(7.150305-0.023326*x.vals)^2
(y.vals[2]-y.vals[1])/y.vals[1] #This is a measure of effect size (proportional increase in number of offspring from a 9 day advance in first nesting date)
y.vals

####Plotting brood cells vs day of first nesting (Fig. 3B):####
#ggplot (without Jane), brood cells before Aug. 1 only, colour-coding by Taraxacum usage:
cols<-c(adjustcolor("purple3",alpha=0.8),adjustcolor("orange",alpha=0.8))
brood.cells.vs.date<-ggplot(data[-163,],aes(x=Day.of.year.first.cell,y=Offspring.produced.Aug1))+
  scale_color_manual(values=cols, 
                    name=expression(paste(italic("Taraxacum")," in pollen provision:")),
                    breaks=c("N","Y"),
                    labels=c("No","Yes"))+
#Colour-coding by species instead:
#  scale_color_manual(values=c(adjustcolor("orange",alpha=0.7),adjustcolor("navy",alpha=0.7),adjustcolor("navy",alpha=0.5)), 
#                     name="Species",
#                     breaks=c("Osmia coloradensis", "Osmia montana","Osmia subaustralis"))+
#  stat_smooth(method = "lm", formula = y ~ I(x^2), colour="black",size=1)+
  stat_smooth(method = "lm", colour="black",size=1.5)+
  geom_point(aes(col=Taraxacum.Y.N),position=position_jitter(width=0.5,height=0)) +
  #ylab(expression(sqrt("Brood cells produced before Aug. 1"))) +
  ylab("Brood cells produced before Aug. 1") +
  xlab("Day of year of first nesting")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(colour="black"))
  #theme(legend.text=element_text(face="italic"))
  #theme(legend.position="None") 
brood.cells.vs.date + scale_x_continuous(breaks = c(152, 182, 213),labels = c("Jun 1", "Jul 1", "Aug 1"))

#ggplot (without Jane), total offspring, no colour-coding (for PPT, not MS):
cols<-adjustcolor("orange",alpha=0.8)
brood.cells.vs.date<-ggplot(data[-163,],aes(x=Day.of.year.first.cell,y=Offspring.produced.Aug1))+
  scale_color_manual(values=cols)+
  stat_smooth(method = "lm", colour="black",size=1,fill=grey(0.8))+
  geom_point(aes(col=cols),position=position_jitter(width=0.5,height=0)) +
  #ylab(expression(sqrt("Brood cells produced before Aug. 1"))) +
  #ylab("Brood cells produced") +
  ylab("Brood cells produced before Aug. 1") +
  xlab("Day of year of first nesting")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(colour="black"))+
#theme(legend.text=element_text(face="italic"))
  theme(legend.position="None") 
brood.cells.vs.date<-brood.cells.vs.date + scale_x_continuous(breaks = c(152, 182, 213),labels = c("Jun 1", "Jul 1", "Aug 1"))
tiff("Aug.1.brood.cells.vs.date.tif", units="in", width=5, height=4, res=288,pointsize=12)
brood.cells.vs.date
dev.off()

#####Q3a. Do bees that use dandelion pollen produce more brood cells before Aug. 1?####
#Full model, treating Offspring produced before Aug. 1 (square root transformed) as response, 
#Taraxacum use as proportion, and allowing for effect of Taraxacum to vary by species:
dand.fitness.model<-lmer(Offspring.produced.Aug1^0.5 ~ Proportion.Taraxacum * Species + (1|Site.code) + (1|Year),data=data)
anova(dand.fitness.model)
#Singular fit and no interaction.Simplifying:
dand.fitness.model.1<-lmer(Offspring.produced.Aug1^0.5 ~ Proportion.Taraxacum + Species + (1|Site.code) + (1|Year),data=data)
#Still singular fit.
summary(dand.fitness.model.1)
anova(dand.fitness.model.1) 
#Significant positive effect of Taraxacum usage, no effect of species.
vif(dand.fitness.model.1)
plot(dand.fitness.model.1)
qqnorm(residuals(dand.fitness.model.1)) #Jane is still a bit of an outlier.
#Repeating without Jane:
dand.fitness.model.1<-lmer(Offspring.produced.Aug1^0.5 ~ Proportion.Taraxacum + Species + (1|Site.code) + (1|Year),data=data[-163,])
#Singular fit, but results unchanged.
summary(dand.fitness.model.1)
anova(dand.fitness.model.1) 

####Q3a with Taraxacum binary:####
#Simple analysis to get group means:
t.test(data$Offspring.produced.Aug1~data$Taraxacum.Y.N) #Jane is an outlier.
t.test(Offspring.produced.Aug1~Taraxacum.Y.N,data=data[-163,]) #This is the result without Jane.
#Actual mixed-model analysis:
dand.fitness.model<-lmer(Offspring.produced.Aug1^0.5 ~ Taraxacum.Y.N * Species + (1|Site.code) + (1|Year),data=data)
summary(dand.fitness.model)
anova(dand.fitness.model)
dand.fitness.model.1<-lmer(Offspring.produced.Aug1^0.5 ~ Taraxacum.Y.N + Species + (1|Site.code) + (1|Year),data=data)
dand.fitness.model.1<-lm(Offspring.produced.Aug1^0.5 ~ Taraxacum.Y.N + Species + Site.code + as.factor(Year),data=data)
#Singular fit; no variance for year or site. Treating those as fixed doesn't change things much.
summary(dand.fitness.model.1)
anova(dand.fitness.model.1) 
#Significant positive effect of Taraxacum usage, no effect of species.
vif(dand.fitness.model.1)
plot(dand.fitness.model.1)
qqnorm(residuals(dand.fitness.model.1)) #Jane is still a bit of an outlier.
#Predicted number of offspring with vs. without Taraxacum:
x.vals<-c(0,1)
y.vals<-(1.980422+0.590171*x.vals)^2
(y.vals[2]-y.vals[1])/y.vals[1] #This is a measure of effect size with Jane included. 
y.vals
  #Repeating without Jane:
dand.fitness.model.1<-lmer(Offspring.produced.Aug1^0.5 ~ Taraxacum.Y.N + Species + (1|Site.code) + (1|Year),data=data[-163,])
#Singular fit, but results unchanged.
summary(dand.fitness.model.1)
anova(dand.fitness.model.1) 
vif(dand.fitness.model.1)
plot(dand.fitness.model.1)
qqnorm(residuals(dand.fitness.model.1))
#Predicted number of offspring with vs. without Taraxacum:
x.vals<-c(0,1)
y.vals<-(1.93497+0.52099*x.vals)^2
(y.vals[2]-y.vals[1])/y.vals[1] #This is a measure of effect size without Jane 
y.vals

####Q3b. Do bees that use dandelion pollen produce more brood cells overall?####
#Full model, treating Offspring produced (log transformed) as response, 
#Taraxacum use as proportion, and allowing for effect of Taraxacum to vary by species:
model<-lmer(sqrt(Total.offspring) ~ Proportion.Taraxacum * Species + (1|Site.code) + (1|Year),data=data)
summary(model)
anova(model)
#No effect of Taraxacum usage, no species*T interaction. Simplifying:
model.1<-lmer(sqrt(Total.offspring) ~ Proportion.Taraxacum + Species + (1|Site.code) + (1|Year),data=data)
summary(model.1)
anova(model.1) 
#No effect of Taraxacum usage, no effect of species.
vif(model.1) #no problem.
plot(model.1)
qqnorm(residuals(model.1)) #There is a definite outlier here: Jane (obs. 163)
#Re-running without Jane:
model.1<-lmer(sqrt(Total.offspring) ~ Proportion.Taraxacum + Species + (1|Site.code) + (1|Year),data=data[-163,])
summary(model.1)
anova(model.1) 
#Doesn't change anything.

####Q3b with Taraxacum binary####
#Full model, treating Offspring produced (log transformed) as response, 
#allowing for effect of Taraxacum to vary by species:
model<-lmer(sqrt(Total.offspring) ~ Taraxacum.Y.N * Species + (1|Site.code) + (1|Year),data=data)
summary(model)
anova(model)
#Significant effect of Taraxacum usage, no species*T interaction. Simplifying:
model.1<-lmer(sqrt(Total.offspring) ~ Taraxacum.Y.N + Species + (1|Site.code) + (1|Year),data=data)
#Singular fit, no variance for site.
summary(model.1)
anova(model.1) 
vif(model.1) #no problem.
plot(model.1)
qqnorm(residuals(model.1)) #There is a definite outlier here: Jane (obs. 163)
#Re-running without Jane:
model.1<-lmer(sqrt(Total.offspring) ~ Taraxacum.Y.N + Species + (1|Site.code) + (1|Year),data=data[-163,])
summary(model.1)
anova(model.1) 
#Significant effect of Taraxacum usage.
#Predicted number of offspring with vs. without Taraxacum:
x.vals<-c(0,1)
y.vals<-(2.54991+0.41070*x.vals)^2
(y.vals[2]-y.vals[1])/y.vals[1] #This is a measure of effect size with Jane included
y.vals
#Predicted number of offspring with vs. without Taraxacum:
y.vals<-(2.479834+0.381254*x.vals)^2
(y.vals[2]-y.vals[1])/y.vals[1] #This is a measure of effect size with Jane excluded
y.vals

####Q3a with Taraxacum binary, for individual species (for PPT, not MS) ####
#Full model, treating Offspring produced (log transformed) as response 
col<-subset(data,Species=="Osmia coloradensis")
mon<-subset(data,Species=="Osmia montana")
sub<-subset(data,Species=="Osmia subaustralis")
model<-lmer(sqrt(Offspring.produced.Aug1) ~ Taraxacum.Y.N + (1|Site.code) + (1|Year),data=col)
summary(model)
anova(model)
#O. col.: Marginal effect of Taraxacum usage, and fails to converge, with or without Jane (obs 47). Little variance for site.
#O. mon.: Signifiant effect of Taraxacum usage, singular fit (no variance for site or year)
#O. sub.: No effect of T usage; no problems.
plot(model)
qqnorm(residuals(model))

####Plotting offspring produced <Aug. 1 vs. Taraxacum usage (Fig. 3c):####
#ggplot (Jane excluded):
cols<-c(adjustcolor("purple3",alpha=0.8),adjustcolor("orange",alpha=0.8))
brood.vs.tarax<-ggplot(data[-163,],aes(x=Species,y=Offspring.produced.Aug1,fill=Taraxacum.Y.N))+
  scale_fill_manual(values=cols, 
                    name=expression(paste(italic("Taraxacum")," in pollen provision:")),
                    breaks=c("N","Y"),
                    labels=c("No","Yes"))+
  scale_x_discrete(breaks=c("Osmia coloradensis", "Osmia montana","Osmia subaustralis"),labels=c("O. coloradensis", "O. montana","O. subaustralis"))  +
  geom_boxplot()+
#  xlab(expression(paste(italic("Taraxacum")," in pollen provision"))) +
#  ylab(expression(sqrt("Brood cells produced before Aug. 1"))) +
  ylab("Brood cells produced before Aug. 1") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(colour="black"))+
  theme(axis.text.x=element_text(face="italic"))
#theme(legend.position="None") 
brood.vs.tarax

#Same but for PPT:
cols<-c(adjustcolor("purple3",alpha=0.8),adjustcolor("orange",alpha=0.8))
brood.vs.tarax<-ggplot(data[-163,],aes(x=Species,y=Offspring.produced.Aug1,fill=Taraxacum.Y.N))+
  scale_fill_manual(values=cols, 
                    name=expression(paste(italic("Taraxacum")," in pollen provision:")),
                    breaks=c("N","Y"),
                    labels=c("No","Yes"))+
  scale_x_discrete(breaks=c("Osmia coloradensis", "Osmia montana","Osmia subaustralis"),labels=c("O. coloradensis", "O. montana","O. subaustralis"))  +
  geom_boxplot()+
  ylab("Brood cells produced before Aug. 1") +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(colour="black"))+
  theme(axis.text.x=element_text(face="italic"))+
  theme(legend.position="top") 
brood.vs.tarax

tiff("brood.cells.vs.Tarax.tif", units="in", width=5, height=5, res=288,pointsize=12)
brood.vs.tarax
dev.off()

#ggplot for PPT, not MS. Jane excluded, Taraxacum usage as proportion:
cols<-adjustcolor("orange",alpha=0.8)
brood.vs.tarax<-ggplot(data[-163,],aes(x=Proportion.Taraxacum,y=Offspring.produced.Aug1))+
  scale_color_manual(values=cols)+
  geom_point(aes(col=cols)) +
  xlab(expression(paste("Proportion",italic("Taraxacum")," in pollen provision"))) +
  #  ylab(expression(sqrt("Brood cells produced before Aug. 1"))) +
  ylab("Brood cells produced before Aug. 1") +
  stat_smooth(method = "lm", colour="black",size=1,fill=grey(0.8))+
  #geom_point(aes(col=cols),position=position_jitter(width=0.5,height=0)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(colour="black"))+
  #theme(legend.text=element_text(face="italic"))
  theme(legend.position="None") 
brood.vs.tarax

tiff("Aug.1.brood.cells.vs.prop.Tarax.tif", units="in", width=5, height=4, res=288,pointsize=12)
brood.vs.tarax
dev.off()

##########
#Re-analysis, excluding VB:
data<-subset(data,Site.code!="VB")
summary(data)

#####Q1a - treating dandelion usage as a proportion:####
dand.phen.model<-lmer(Day.of.year.first.cell ~ Proportion.Taraxacum * Species + (1|Site.code) + (1|Year),data=data)
summary(dand.phen.model)
anova(dand.phen.model)
vif(dand.phen.model)
plot(dand.phen.model) #heterogeneity of variance
qqnorm(residuals(dand.phen.model))
#Strong negative effect of Taraxacum usage, no species*T interaction. Simplifying:
dand.phen.model.1<-lmer(Day.of.year.first.cell ~ Proportion.Taraxacum + Species + (1|Site.code) + (1|Year),data=data)
summary(dand.phen.model.1)
anova(dand.phen.model.1) 
vif(dand.phen.model.1)
plot(dand.phen.model.1) #heterogeneity of variance
qqnorm(residuals(dand.phen.model.1))
#Significant effect of species, with subaustralis significantly earlier than coloradensis. Strong effect of Taraxacum usage.
#Highly significant negative effect of Taraxacum usage.
#Re-order to check whether Cephalosmia species differ from one another:
dand.phen.model.1<-lmer(Day.of.year.first.cell ~ Proportion.Taraxacum + relevel(Species,"Osmia montana") + (1|Site.code) + (1|Year),data=data)
summary(dand.phen.model.1)
#Yes, O. subaustralis is earlier than both other species, by 5-9 days (after accounting for effect of Taraxacum usage).

#####Q1b - treating dandelion usage as binary:####
#Full model, treating DOY as response, Taraxacum use as binary, and allowing for effect of 
#Taraxacum to vary by species:
dand.phen.model<-lmer(Day.of.year.first.cell ~ Taraxacum.Y.N * Species + (1|Site.code) + (1|Year),data=data)
summary(dand.phen.model)
anova(dand.phen.model)
#Signficant effect of Taraxacum usage, no species*T interaction. Simplifying:
dand.phen.model.1<-lmer(Day.of.year.first.cell ~ Taraxacum.Y.N + Species + (1|Site.code) + (1|Year),data=data)
summary(dand.phen.model.1)
anova(dand.phen.model.1)
vif(dand.phen.model.1)
plot(dand.phen.model.1)
qqnorm(residuals(dand.phen.model.1)) #No problems.
#Significant effect of species, with subaustralis significantly earlier than coloradensis
#Highly significant negative effect of Taraxacum usage.
#Re-order to check whether Cephalosmia species differ from one another:
dand.phen.model.1<-lmer(Day.of.year.first.cell ~ Taraxacum.Y.N + relevel(Species,"Osmia montana") + (1|Site.code) + (1|Year),data=data)
summary(dand.phen.model.1)
#Yes, O. subaustralis is earlier than both other species, by 6-7 days.

#####Q3a. Do bees that use dandelion pollen produce more brood cells before Aug. 1?####
#Full model, treating Offspring produced before Aug. 1 (square root transformed) as response, 
#Taraxacum use as proportion, and allowing for effect of Taraxacum to vary by species:
dand.fitness.model<-lmer(Offspring.produced.Aug1^0.5 ~ Proportion.Taraxacum * Species + (1|Site.code) + (1|Year),data=data)
anova(dand.fitness.model)
#Singular fit and no interaction.Simplifying:
dand.fitness.model.1<-lmer(Offspring.produced.Aug1^0.5 ~ Proportion.Taraxacum + Species + (1|Site.code) + (1|Year),data=data)
#Still singular fit; no variance for site.
summary(dand.fitness.model.1)
anova(dand.fitness.model.1) 
#Significant positive effect of Taraxacum usage, no effect of species.
vif(dand.fitness.model.1)
plot(dand.fitness.model.1)
qqnorm(residuals(dand.fitness.model.1)) #Jane is still a bit of an outlier.
#Repeating without Jane:
dand.fitness.model.1<-lmer(Offspring.produced.Aug1^0.5 ~ Proportion.Taraxacum + Species + (1|Site.code) + (1|Year),data=data[-145,])
#Singular fit, but results unchanged.
summary(dand.fitness.model.1)
anova(dand.fitness.model.1) 

####Q3a with Taraxacum binary:####
dand.fitness.model<-lmer(Offspring.produced.Aug1^0.5 ~ Taraxacum.Y.N * Species + (1|Site.code) + (1|Year),data=data)
summary(dand.fitness.model)
#Singular fit, no variance for site.
anova(dand.fitness.model)
dand.fitness.model.1<-lmer(Offspring.produced.Aug1^0.5 ~ Taraxacum.Y.N + Species + (1|Site.code) + (1|Year),data=data)
dand.fitness.model.1<-lm(Offspring.produced.Aug1^0.5 ~ Taraxacum.Y.N + Species + Site.code + as.factor(Year),data=data)
#Singular fit; no variance for year or site. Treating those as fixed doesn't change things much.
summary(dand.fitness.model.1)
anova(dand.fitness.model.1) 
#Significant positive effect of Taraxacum usage, no effect of species.
vif(dand.fitness.model.1)
plot(dand.fitness.model.1)
qqnorm(residuals(dand.fitness.model.1)) #Jane is still a bit of an outlier.
#Predicted number of offspring with vs. without Taraxacum:
x.vals<-c(0,1)
y.vals<-(1.9647+0.6433*x.vals)^2
(y.vals[2]-y.vals[1])/y.vals[1] #This is a measure of effect size with Jane included. 
y.vals
#Repeating without Jane:
dand.fitness.model.1<-lmer(Offspring.produced.Aug1^0.5 ~ Taraxacum.Y.N + Species + (1|Site.code) + (1|Year),data=data[-145,])
#Singular fit, but results unchanged.
summary(dand.fitness.model.1)
anova(dand.fitness.model.1) 
vif(dand.fitness.model.1)
plot(dand.fitness.model.1)
qqnorm(residuals(dand.fitness.model.1))
#Predicted number of offspring with vs. without Taraxacum:
x.vals<-c(0,1)
y.vals<-(1.87194+0.65831*x.vals)^2
(y.vals[2]-y.vals[1])/y.vals[1] #This is a measure of effect size without Jane 
y.vals

####Q3b. Do bees that use dandelion pollen produce more brood cells overall?####
#Full model, treating Offspring produced (log transformed) as response, 
#Taraxacum use as proportion, and allowing for effect of Taraxacum to vary by species:
model<-lmer(sqrt(Total.offspring) ~ Proportion.Taraxacum * Species + (1|Site.code) + (1|Year),data=data)
summary(model)
anova(model)
#No effect of Taraxacum usage, no species*T interaction. Simplifying:
model.1<-lmer(sqrt(Total.offspring) ~ Proportion.Taraxacum + Species + (1|Site.code) + (1|Year),data=data)
#Singular fit; little variance for site.
summary(model.1)
anova(model.1) 
#No effect of Taraxacum usage, no effect of species.
vif(model.1) #no problem.
plot(model.1)
qqnorm(residuals(model.1)) #There is a definite outlier here: Jane (obs. 145)
#Re-running without Jane:
model.1<-lmer(sqrt(Total.offspring) ~ Proportion.Taraxacum + Species + (1|Site.code) + (1|Year),data=data[-145,])
summary(model.1)
anova(model.1) 
#Doesn't change anything.

####Q3b with Taraxacum binary####
#Full model, treating Offspring produced (log transformed) as response, 
#Taraxacum use as proportion, and allowing for effect of Taraxacum to vary by species:
model<-lmer(sqrt(Total.offspring) ~ Taraxacum.Y.N * Species + (1|Site.code) + (1|Year),data=data)
summary(model)
anova(model)
#Significant effect of Taraxacum usage, no species*T interaction. Simplifying:
model.1<-lmer(sqrt(Total.offspring) ~ Taraxacum.Y.N + Species + (1|Site.code) + (1|Year),data=data)
#Singular fit, no variance for site.
summary(model.1)
anova(model.1) 
vif(model.1) #no problem.
plot(model.1)
qqnorm(residuals(model.1)) #There is a definite outlier here: Jane (obs. 145)
#Re-running without Jane:
model.1<-lmer(sqrt(Total.offspring) ~ Taraxacum.Y.N + Species + (1|Site.code) + (1|Year),data=data[-145,])
summary(model.1)
anova(model.1) 
#Significant effect of Taraxacum usage.
#Predicted number of offspring with vs. without Taraxacum:
x.vals<-c(0,1)
y.vals<-(2.53697+0.48048*x.vals)^2
(y.vals[2]-y.vals[1])/y.vals[1] #This is a measure of effect size with Jane included
y.vals
#Predicted number of offspring with vs. without Taraxacum:
y.vals<-(2.46675+0.45064*x.vals)^2
(y.vals[2]-y.vals[1])/y.vals[1] #This is a measure of effect size with Jane excluded
y.vals
