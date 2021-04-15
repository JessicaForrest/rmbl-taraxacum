#Script by Megan McAulay 2018, modified by JF Oct. 2019
#See here for brief summary of difference between inclusion of cluster and 
#random terms (the latter = 'frailty' models if allowing random intercepts only):
#https://stat.ethz.ch/pipermail/r-help/2011-June/282065.html
#For more, see: Modeling Survival Data: Extending the Cox Model by T. Therneau and P. Grambsch (2000)
#available electronically from uOttawa library. However, it's old and a bit out of date.
#Plotting component updated 15 June 2020 in response to feedback from Paul CaraDonna.

#load packages
library("survival")
library("survminer")
#library(survMisc)
library(coxme)

#set working directory
setwd("C:/Users/jforr/Dropbox/RMBL/RMBL trapnests/Trapnests 2019/Charlotte's project/Data files for archiving")

####load datasets ####
transfer.data<-read.delim("Taraxacum egg-transfer expt 2018-2019.txt")
summary(transfer.data)
#transfer.data$Days.to.event: Event is cocoon formation (= censored data) or death or last day of observation.
#transfer.data$Survival: Coded as 0 if larva still alive at end of season (includes those in cocoons); 1 if dead.

cephal<-subset(transfer.data,Bee.subgenus=="Cephalosmia")
helic<-subset(transfer.data,Bee.subgenus=="Helicosmia")

all.data<-read.delim("Charlotte's egg-transfer expts merged.txt")
summary(all.data)

####ANALYSIS#####
#Most basic cox ph model (no cluster or random terms, no year or species terms):
survival.test<-coxph(Surv(Days.to.event, Survival)~Pollen.type,data=transfer.data)
summary(survival.test)

#Less basic cox ph model (no cluster or random terms):
survival.test<-coxph(Surv(Days.to.event, Survival)~Pollen.type*Bee.subgenus+as.factor(Year),data=transfer.data)
summary(survival.test)
#Interaction ns, removed:
survival.test.1<-coxph(Surv(Days.to.event, Survival)~Pollen.type+Bee.subgenus+as.factor(Year),data=transfer.data)
summary(survival.test.1)
#Year ns, removed:
survival.test.2<-coxph(Surv(Days.to.event, Survival)~Pollen.type+Bee.subgenus,data=transfer.data)
summary(survival.test.2)
AIC(survival.test,survival.test.1,survival.test.2)

#Ideally we would also account for origin pollen, but it is too confounded with subgenus:
cor.test(as.numeric(Bee.subgenus),as.numeric(Origin.pollen),type=c("spearman"))
#Fortunately, bee subgenus and destination pollen type are independent:
cor.test(as.numeric(Bee.subgenus),as.numeric(Pollen.type),type=c("spearman"))
#Origin pollen type could be analysed in a Helicosmia-only model (see below).

#cox proportional mixed model model with random term for bee nest ID. This can be called a "frailty model"
# - see Austin 2017, A Tutorial on Multilevel Survival Analysis: Methods, Models and Applications.
Survival_coxme<-coxme(Surv(Days.to.event, Survival)~Pollen.type * Bee.subgenus + as.factor(Year) + (1|Bee.source), data=transfer.data)
summary(Survival_coxme)
Survival_coxme.1<-coxme(Surv(Days.to.event, Survival)~Pollen.type + Bee.subgenus + as.factor(Year) + (1|Bee.source), data=transfer.data)
summary(Survival_coxme.1)
#Year ns; simplifying:
Survival_coxme.2<-coxme(Surv(Days.to.event, Survival)~Pollen.type + Bee.subgenus + (1|Bee.source), data=transfer.data)
summary(Survival_coxme.2)
AIC(Survival_coxme.2,Survival_coxme.1, Survival_coxme)

#Simple cox proportional hazards model:
Survival_coxph<-coxph(Surv(Days.to.event, Survival)~Pollen.type + Bee.subgenus + as.factor(Year), data=transfer.data)
summary(Survival_coxph)
Survival_coxph_b<-coxph(Surv(Days.to.event, Survival)~Pollen.type + Bee.subgenus, data=transfer.data)
summary(Survival_coxph_b)

#Cox PH model with cluster term:
Survival_coxph_cluster<-coxph(Surv(Days.to.event, Survival)~Pollen.type + Bee.subgenus + as.factor(Year) + cluster(Bee.source), data=transfer.data)
summary(Survival_coxph_cluster)
Survival_coxph_cluster_b<-coxph(Surv(Days.to.event, Survival)~Pollen.type + Bee.subgenus + cluster(Bee.source), data=transfer.data)
summary(Survival_coxph_cluster_b)

#compare AIC values of coxme and coxph models:
AIC(Survival_coxme,Survival_coxph, Survival_coxph_cluster)
#Mixed-effects model is better, but all basically give the same results.
AIC(Survival_coxme.2,Survival_coxph_b, Survival_coxph_cluster_b)
#Same result.

#If we just look at individual subgenera:
curves<-survfit(Surv(Days.to.event, Survival)~Pollen.type,data=helic)
#summary(curves)
par(mar=c(5,5,1,1))
plot(curves,lty=c(1,2),col=c("slateblue","orange"),lwd=2,conf.int=F,ylab="Proportion surviving",xlab="Days since egg laid",mark.time=T,mark=10)
legend("topright", c("native Asteraceae", "dandelion"), lty = 1:2,adj=0,col=c("slateblue","orange"),bty="n",lwd=2) 
survival.test<-coxme(Surv(Days.to.event, Survival)~Pollen.type*Origin.pollen+ (1|Bee.source),data=helic)
summary(survival.test)
#Interaction ns; removing:
survival.test<-coxme(Surv(Days.to.event, Survival)~Pollen.type+Origin.pollen+ (1|Bee.source),data=helic)
summary(survival.test)
#Origin pollen Taraxacum decreases survival.
#Simplest model, testing only effect of pollen type, to obtain hazard ratio for Pollen type only:
survival.test<-coxme(Surv(Days.to.event, Survival)~Pollen.type+ (1|Bee.source),data=helic)
summary(survival.test)
curves<-survfit(Surv(Days.to.event, Survival)~Pollen.type+Origin.pollen,data=helic)
plot(curves,lty=c(1,1,2,2),col=c("red3","orange"),lwd=2,conf.int=F,ylab="Proportion surviving",xlab="Days since egg laid",mark.Days.to.event=T,mark=10)
legend("bottomleft", c("native Asteraceae", expression(paste(italic("Taraxacum")))), lty = 1:2,
       adj=0,col="black",bty="n",lwd=2) 
text(0,0.2,adj=0,"reared on:")
text(33,0.68,adj=0,"origin Asteraceae",font=1,col="red3")
text(40,0.11,adj=0,expression(paste("origin ",italic("Taraxacum"))),font=3,col="orange")

#Running model with only bees that originated on Asteraceae pollen, since survival is lower
#for those that originate on Taraxacum:
detach(transfer.data)
aster.origin<-subset(transfer.data,Origin.pollen=="Asteraceae")
summary(aster.origin)
attach(aster.origin)
Days.to.event <-Days.to.event
event <- Survival
curves<-survfit(Surv(Days.to.event, Survival)~Pollen.type,data=aster.origin)
par(mar=c(5,5,1,1))
plot(curves,lty=c(1,2),col=c("slateblue","orange"),lwd=2,conf.int=F,ylab="Proportion surviving",xlab="Days since egg laid",mark.time=T,mark=10)
legend("topright", c("native Asteraceae", "dandelion"), lty = 1:2,adj=0,col=c("slateblue","orange"),bty="n",lwd=2) 
Survival_coxme.2<-coxme(Surv(Days.to.event, Survival)~Pollen.type + Bee.subgenus + (1|Bee.source), data=aster.origin)
summary(Survival_coxme.2)
#Bee.subgenus effect ns; simplifying:
Survival_coxme.3<-coxme(Surv(Days.to.event, Survival)~Pollen.type  + (1|Bee.source), data=aster.origin)
summary(Survival_coxme.3)
#Marginal effect of pollen type.
detach(aster.origin)

#What's the situation with "controls", i.e. bees that weren't transfered?
control.data<-subset(all.data,Transfer.type=="Control")
attach(control.data)

Days.to.event <-Days.to.event #Event is cocoon formation (= censored data) or death or last day of observation.
event <- Survival #Coded as 0 if larva still alive at end of season (includes those in cocoons); 1 if dead.
curves<-survfit(Surv(Days.to.event, Survival)~Pollen.type,data=control.data)
summary(curves)
par(mar=c(5,5,1,1))
plot(curves,lty=1,col=c("slateblue","purple","orange"),lwd=2,conf.int=F,ylab="Proportion surviving",xlab="Days since egg laid",mark.time=T,mark=10)
legend("bottomleft", c("native Asteraceae","mix","dandelion"), lty = 1,adj=0,
  col=c("slateblue","purple","orange"),bty="n",lwd=2) 
#sample size is small, but still:
Survival_coxme<-coxme(Surv(Days.to.event, Survival)~Pollen.type  + (1|Bee.source), data=control.data)
#Insufficient information, since no deaths except on Taraxacum pollen.
Survival_coxph<-coxph(Surv(Days.to.event, Survival)~Pollen.type, data=control.data)
#Same result. Anyway, pollen type is confounded with subgenus, so this isn't meaningful.

####Test Assumptions of Cox PH#####

## 1) proportional hazards (PH)

#Schoenfeld residuals are independent of time (assume PH if p values not statistically significant).
#This doesn't work for ME models.
test.ph <-cox.zph(Survival_coxph)
test.ph
#PH assumption met for pollen type but not bee subgenus or year. Try with simpler model:
Survival_coxph_b<-coxph(Surv(Days.to.event, Survival)~Pollen.type + Bee.subgenus, data=transfer.data)
summary(Survival_coxph_b)
test.ph <-cox.zph(Survival_coxph_b)
test.ph
#Subgenus still a problem.
test.ph <-cox.zph(Survival_coxph_cluster)
test.ph
#Not so bad for subgenus, but year is still a problem.
test.ph <-cox.zph(Survival_coxph_cluster_b)
test.ph
#Now subgenus becomes a problem again.

#plot residuals
ggcoxzph(test.ph)
#Hazards seem to increase over time for Helicosmia. This matches survival curves.

## 2) testing influential observations

#dfbeta values (look for major outliers)
ggcoxdiagnostics(Survival_coxph_cluster_b, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())

#deviance residuals (pattern should be fairly symmetric around 0)
ggcoxdiagnostics(Survival_coxph_cluster_b, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())

#Looks OK.

####PLOTTING####

#Preliminary plots----
#Plot survival by pollen type (simple colour plot for PPT):
curves<-survfit(Surv(Days.to.event, Survival)~Pollen.type,data=transfer.data)
tiff("Tarax.survival.tif", units="in", width=5, height=4, res=288,pointsize=12)
par(mar=c(5,5,1,1))
plot(curves,lty=1,lwd=2,conf.int=F,ylab="Proportion surviving",xlab="Days since egg laid",
     col=c(adjustcolor("purple3",alpha=0.7),adjustcolor("orange",alpha=0.7)),mark.time=T,mark=19, cex=1)
legend("topright", c("native Asteraceae", "non-native dandelion"), lty = 1:2,adj=0,col=c(adjustcolor("purple3",alpha=0.7),adjustcolor("orange",alpha=0.7)),bty="n",lwd=2) 
dev.off()

#Plot survival by subgenus:
curves<-survfit(Surv(Days.to.event, Survival)~Bee.subgenus,data=transfer.data)
plot(curves,lty=c(1,2),col=c("slateblue","orange"),lwd=2,conf.int=F,
     ylab="Proportion surviving",xlab="Days since egg laid",mark.time=T,mark=1, cex=0.8)
legend("topright", c("Cephalosmia", "Helicosmia"), lty = 1:2,text.font=3,adj=0,col=c("slateblue","orange"),bty="n",lwd=2) 

#Plot survival by subgenus and pollen type. :
curves<-survfit(Surv(Days.to.event, Survival)~Pollen.type+Bee.subgenus,data=transfer.data)
#summary(curves)
par(mar=c(5,5,1,1))
plot(curves,lty=c(1,1,2,2),col=c("slateblue","orange"),lwd=2,conf.int=F,ylab="Proportion surviving",xlab="Days since egg laid",mark.time=T,mark=10)
legend("bottomleft", c("native Asteraceae", expression(paste(italic("Taraxacum")))), lty = 1:2,
       adj=0,col="black",bty="n",lwd=2) 
text(40,0.69,adj=0,"Cephalosmia",font=3,col="slateblue")
text(38,0.25,adj=0,"Helicosmia",font=3,col="orange")

#Final plot----
#Plot survival by subgenus and pollen type in different panels. This is Fig. 4 of manuscript:
par(mfrow=c(1,2))
#summary(curves)
curves<-survfit(Surv(Days.to.event, Survival)~Pollen.type,data=helic)
par(mar=c(5,5,2,1))
plot(curves,lty=1,lwd=2,conf.int=F,ylab="Proportion surviving",xlab="Days since egg laid",
     col=c(adjustcolor("purple3",alpha=0.7),adjustcolor("orange",alpha=0.7)),mark.time=T,mark=19, cex=1)
text(60,0.95,"Helicosmia",font=3,adj=1,cex=1)
text(60,0.85,"Hazard ratio = 1.5",adj=1,cex=0.9)
mtext("(a)",side=3,adj=0,line=0.5)
curves<-survfit(Surv(Days.to.event, Survival)~Pollen.type,data=cephal)
par(mar=c(5,0,2,6))
plot(curves,lty=1,lwd=2,conf.int=F,ylab=" ",xlab="Days since egg laid",
     col=c(adjustcolor("purple3",alpha = 0.7),adjustcolor("orange",alpha = 0.7)),mark.time=T,mark=19, cex=1,
     axes=F)
axis(1); axis(2, at=seq(0,1,0.2), labels=FALSE); box()
legend("bottomleft", c("native Asteraceae", expression(paste("non-native ",italic("Taraxacum")))), col=c(adjustcolor("purple3",alpha = 0.7),adjustcolor("orange",alpha = 0.7)),
       cex=0.9, lty = 1,adj=0,bty="n",lwd=2,title="Pollen type:",title.adj=0.1)
text(60,0.95,"Cephalosmia",font=3,adj=1,cex=1)
text(60,0.85,"Hazard ratio = 2.0",adj=1,cex=0.9)
mtext("(b)",side=3,adj=0,line=0.5)

####MEGAN'S PLOTTING CODE, USED FOR OECOLOGIA MS####

library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(RColorBrewer)

dev.off()

#autoplot with ggplot2 modifications (with legend)
tiff("SurvivalAnalysis.tiff",width=6,height=6,units='in',res=300)

completefigure <- grid.arrange(
  O.coloradensis_autoplot + theme_bw() + expand_limits(y=0) +
    scale_colour_manual(values=c("#377eb8", "#e41a1c", "#4daf4a")) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"),
      plot.title=element_text(face=4),
      legend.position="none",
      axis.title=element_text(size=10,face="bold"),
      axis.text.y =element_text(size=8, colour="black"),
      axis.text.x =element_text(size=8, colour="black")),
  
  O.iridis_autoplot + theme_bw() + expand_limits(y=0) +
    scale_x_continuous(limits=c(0, 50)) + xlab("Days on provision as a larva") +
    scale_colour_manual(values=c("#377eb8", "#e41a1c", "#4daf4a")) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"),
      plot.title=element_text(face=4),
      legend.position="none",
      axis.title=element_text(size=10,face="bold"),
      axis.text.y =element_text(size=8, colour="black"),
      axis.text.x =element_text(size=8, colour="black")),
  
  O.montana_autoplot + theme_bw() + expand_limits(y=0) +
    scale_x_continuous(limits=c(0, 50)) + xlab("Days on provision as a larva") +
    scale_colour_manual(values=c("#377eb8", "#e41a1c", "#4daf4a")) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title=element_text(face=4),
      legend.position="none", axis.line = element_line(colour = "black"),
      axis.title=element_text(size=10,face="bold"),
      axis.text.y =element_text(size=8, colour="black"),
      axis.text.x =element_text(size=8, colour="black")),
  
  O.subaustralis_autoplot + theme_bw() + expand_limits(y=0) +
    scale_colour_manual(values=c("#377eb8", "#e41a1c", "#4daf4a")) +
    scale_x_continuous(limits=c(0, 50)) + xlab("Days on provision as a larva") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"),
      plot.title=element_text(face=4),
      legend.position="none",
      axis.title=element_text(size=10,face="bold"),
      axis.text.y =element_text(size=8, colour="black"),
      axis.text.x =element_text(size=8, colour="black")),
  
  O.tristella_autoplot + theme_bw() + expand_limits(y=0) +
    scale_colour_manual(values=c("#377eb8", "#e41a1c", "#4daf4a")) +
    scale_x_continuous(limits=c(0, 50)) + xlab("Days on provision as a larva") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      axis.line = element_line(colour = "black"),
      plot.title=element_text(face=4),
      legend.position="none",
      axis.title=element_text(size=10,face="bold"),
      axis.text.y =element_text(size=8, colour="black"),
      axis.text.x =element_text(size=8,colour="black")),
  ncol=2, nrow=3)

dev.off()

