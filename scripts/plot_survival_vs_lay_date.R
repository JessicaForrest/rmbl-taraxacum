#Plot probability of larval overwinter survival for RMBL Osmia, Asteraceae specialists only.
#Panel (a) based on RP survival data (2013-2015); (b) based on survival of 2017 nests at multiple sites.
#See dataset-specific scripts in respective folders for detailed analysis of each.
#Black and white version (for manuscript) in part A; colour version (for PPT, using 2013-15 RP data only) in part B.
#See "plot_bee_and_floral_phenology.r" for integration with panel A of Fig. 1 of manuscript.

library(lme4)
library('tibble')
library('dplyr')
library('readr')

#Load 2013-2015 data:
setwd("C:/Users/jforr/Dropbox/RMBL/RMBL trapnests/Nest-height (Rosy Point)")
setwd("C:/Users/jforrest/Dropbox/RMBL/RMBL trapnests/Nest-height (Rosy Point)")
data.16<-read.delim("RP_nest-height_experiment_2013-2016.txt")
data.16<-subset(data.16,Nest.constructed.after.departure!="y") #laid dates too uncertain.
data.16<-subset(data.16,Bee.subgenus!="Hapsidosmia" & Bee.subgenus!="Melanosmia")

#Load 2017 data:
setwd("C:/Users/jforr/Dropbox/RMBL/RMBL trapnests/Trapnests 2017")
data.17<-read.delim("Osmia 2017 survival to emergence.txt")
data.17<-subset(data.17,Subgenus!="Hapsidosmia" & Subgenus!="NA")

####PART A (B&W):####

#Set up plotting area:
setwd("C:/Users/jforr/Dropbox/RMBL/RMBL trapnests/Trapnests 2019/Charlotte's project")
pdf("survival_vs_date_laid.pdf",width=8,height=3)
layout(matrix(c(1:3),1,3),widths=c(0.43,0.38,0.19))
par(cex=1)

#Plot survival as a function of nest initiation date for 2017, aster specialists only: 
#START HERE TO ADD PANELS TO FIGURE 1:
par(mar=c(4.5,4.5,1,1))
plot(Prop.hosts.survived~DOY.first.egg.laid,data=data.17, type="n",axes=F,ylim=c(-0.05,1.05),
     xlab="Date nest initiated",ylab="Survival through winter")
axis(2,cex.axis=0.8); axis(1,at=c(182,197,213),labels=c("Jul 1", "Jul 16","Aug 1"),cex.axis=0.8)

#Plot model fit:
model.raw<-glmer(Count.of.hosts.survived/Cells.with.host.fate.known~
                   DOY.first.egg.laid + (1|Site/Mother.s.individual.code), 
                 family=binomial, weights=Cells.with.host.fate.known,data=data.17)
summary(model.raw)
x.axis<-seq(169,226,1)
a<-21.14377
b<--0.12174
#Note that these parameters are from a model excluding effect of subgenus
resp<-1/(1+exp(-a-b*x.axis))
lines(x.axis,resp,lwd=2)

#Add datapoints:
with(subset(data.17,Subgenus=="Helicosmia"),
     points(DOY.first.egg.laid,Prop.hosts.survived,pch=21,bg=adjustcolor(grey(0.7), alpha=0.7),col=adjustcolor(grey(0.7),alpha=0.7)))
with(subset(data.17,Subgenus=="Cephalosmia"),
     points(DOY.first.egg.laid,Prop.hosts.survived,pch=21,bg=adjustcolor(grey(0.4),alpha=0.7),col=adjustcolor(grey(0.4),alpha=0.7)))

#Panel label:
text(168,1.05,"(b)",adj=0)

#Plot survival as a function of lay date for 2013-2015, aster specialists only: 
par(mar=c(4.5,2,1,0.5))
plot(Alive.dead.2nd.summer~Est.doy.laid,type="n",axes=F,ylim=c(-0.05,1.05),
     xlab="Date egg laid",ylab="",data=data.16)
axis(2,cex.axis=0.8); axis(1, at=c(167,182,197,213,228),labels=c("Jun 16","Jul 1", "Jul 16","Aug 1","Aug 16"),cex.axis=0.8)

#Points by subgenus:
with(subset(data.16,Bee.subgenus=="Cephalosmia"),
     points(Est.doy.laid,jitter(Alive.dead.2nd.summer,0.15),pch=21,bg=adjustcolor(grey(0.4), alpha=0.7),col=adjustcolor(grey(0.4),alpha=0.7)))
with(subset(data.16,Bee.subgenus=="Helicosmia"),
     points(Est.doy.laid,jitter(Alive.dead.2nd.summer,0.15),pch=21,bg=adjustcolor(grey(0.7), alpha=0.7),col=adjustcolor(grey(0.7),alpha=0.7)))

#Plot model fit:
model.raw<-glmer(Alive.dead.2nd.summer~Est.doy.laid
                 +(1|Mother.ID/Nest.ID),family="binomial",data=data.16)
summary(model.raw)
x.axis<-seq(173,230,1)
a<-11.93658
b<--0.06169
#Note that these parameters are from a model excluding effect of year.
resp<-1/(1+exp(-a-b*x.axis))
lines(x.axis,resp,lwd=2)

text(172,1.05,"(c)",adj=0)

#Legend as separate panel:
par(mar=c(4.5,0,1,1))
plot(Prop.hosts.survived~DOY.first.egg.laid,data=data.17, type="n",axes=F,xlab="",ylab="")
legend("left",legend = c("Helicosmia","Cephalosmia"),bty="n",y.intersp = 1.5,inset = 0.01,
       text.font=3,pch=21,pt.bg=c(adjustcolor(grey(0.7),alpha=0.7),adjustcolor(grey(0.4),alpha=0.7)),
       col=c(adjustcolor(grey(0.7),alpha=0.7),adjustcolor(grey(0.4),alpha=0.7)))

#STOP HERE WHEN ADDING PANELS TO FIGURE 1.

dev.off()

####PART B (colour):####
library(ggplot2)

#Plot survival as a function of lay date for 2013-2015, aster specialists only: 
par(mar=c(4.5,5,1,0.5))

surv.vs.date<-ggplot(data.16,aes(x=Est.doy.laid,y=Alive.dead.2nd.summer))+
        stat_smooth(method = "glm", method.args = list(family = "binomial"), colour="black",size=1,fill=grey(0.8))+
        scale_color_manual(values=c(adjustcolor("orange",alpha=0.5),adjustcolor("orange",alpha=0.5)), 
                           name="Bee.subgenus",
                           breaks=c("Cephalosmia", "Helicosmia"))+
        geom_point(aes(colour=Bee.subgenus),size=2, position=position_jitter(width=0.25,height=0.02)) +
        #geom_rug(position=position_jitter(width=1,height=0))+
        ylab("Survival through winter") +
        xlab("Date egg laid")+
        theme_bw()+
        theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                axis.line = element_line(colour = "black"))+
        theme(axis.text=element_text(colour="black"))+
        scale_x_continuous(breaks = c(167,182,197,213,228),labels = c("Jun 16","Jul 1", "Jul 16","Aug 1","Aug 16"))+
        #theme(legend.text=element_text(face="italic"))
        theme(legend.position="None") 
surv.vs.date<-surv.vs.date + scale_y_continuous(breaks = c(0,1),labels = c("0","1"))

tiff("aster.bees.survival.vs.date.laid.tif", units="in", width=5, height=4, res=288,pointsize=12)
surv.vs.date
dev.off()

