library(flowCore)
library(ggcyto)
library(openCyto)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(extrafont)
library(gridExtra)
library(egg)
library(ggridges)
library(ggsci)
theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6, family = "Arial"), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6, family = "Arial", color= "black"),
                  axis.text.x = element_text(vjust = 0.5, color= "black"),
                  axis.text.y = element_text(vjust = 0.5, color= "black"),
                  axis.title = element_text(size = 6, family = "Arial"), strip.text = element_text(size = 6, family = "Arial", color= "black"),
                  strip.background = element_blank(), legend.title = element_blank()))
mypal<- pal_npg("nrc", alpha = 1)(2)
#non-fluorescent control for background-subtraction
MyFlowSet <- read.flowSet(path="C:/Users/noviello/Documents/Analysis of FlowCytometry for modelling/Flow cytometry/20210604/Experiment_013/allfiles", min.limit=0.01)
chnl <- c("FSC-A", "SSC-A") #are the channels to build the gate on
bou_g <- openCyto:::.boundary(MyFlowSet[[3]], channels = chnl, min = c(0.4e5, 0.20e5), max=c(2e5,1.3e5), filterId = "Boundary Gate")
p <- autoplot(MyFlowSet[[3]], x = chnl[1], y = chnl[2], bins=100)
p +geom_gate(bou_g)
BoundaryFrame<- Subset (MyFlowSet[[3]], bou_g) #is the subset of cells within the gate
chnl <- c("FSC-A", "FSC-H")
sing_g <- openCyto:::.singletGate(BoundaryFrame, channels = chnl)
p <- autoplot(BoundaryFrame, x = chnl[1], y = chnl[2])
p + geom_gate(sing_g)
Singlets_MyFlowSet<-Subset(MyFlowSet, bou_g%&% sing_g)
medianexp<- as.data.frame(fsApply(Singlets_MyFlowSet, each_col, median))
mBFP_neg<-mean(medianexp[c(1:3), 7])
mmCherry_neg<-mean(medianexp[c(1:3), 8])

MyFlowSet <- read.flowSet(path="C:/Users/noviello/Documents/Analysis of FlowCytometry for modelling/Flow cytometry/20210813_timecourse/Experiment_013", min.limit=0.01)
chnl <- c("FSC-A", "SSC-A") #are the channels to build the gate on
bou_g <- openCyto:::.boundary(MyFlowSet[[3]], channels = chnl, min = c(0.4e5, 0.20e5), max=c(2e5,1.3e5), filterId = "Boundary Gate")
p <- autoplot(MyFlowSet[[3]], x = chnl[1], y = chnl[2], bins=100)
p +geom_gate(bou_g)
BoundaryFrame<- Subset (MyFlowSet[[3]], bou_g) #is the subset of cells within the gate
chnl <- c("FSC-A", "FSC-H")
sing_g <- openCyto:::.singletGate(BoundaryFrame, channels = chnl)
p <- autoplot(BoundaryFrame, x = chnl[1], y = chnl[2])
p + geom_gate(sing_g)
Singlets_MyFlowSet<-Subset(MyFlowSet, bou_g%&% sing_g)
# extract normalized median expression----

Transf<-linearTransform(transformationId="defaultLinearTransform", a = 1, b = -mBFP_neg)
Norm_Singlets_MyFlowSet<- transform(Singlets_MyFlowSet, transformList('BV421-A' ,Transf))
Transf<-linearTransform(transformationId="defaultLinearTransform", a = 1, b = -mmCherry_neg)
Norm_Singlets_MyFlowSet<- transform(Norm_Singlets_MyFlowSet, transformList('PE-A' ,Transf))
medianexp<- as.data.frame(fsApply(Norm_Singlets_MyFlowSet, each_col, median))
medianexp<-medianexp[, c(7:8)]


medianexp %>% rownames_to_column()->medianexp
medianexp %>% separate(1, c( NA, NA, "plasmid", "exp", "rep", "time", NA, NA), sep = "_") ->medianexp
medianexp


medianexp$time<-as.numeric(medianexp$time)

medianexp$plasmid<-as.factor(medianexp$plasmid)
medianexp %>% group_by(plasmid) %>% fct_relevel(plasmid, levels=c("SP430", "SP428", "SP430ABA", "SP427", "SP411"))
medianexp %>% dplyr::filter(time!=125)->medianexp
medianexp %>% dplyr::filter (exp == "KD" )->KD
KD %>% dplyr::filter(time==0)->KD0

KD0 %>% group_by(plasmid) %>% summarize( mmcherry=mean(`PE-A`))->t0
merge(x=KD,y=t0,  all.x=T)->KD
KD %>% dplyr::filter(time==150)->KD150
KD150 %>% group_by(plasmid) %>% summarize( mbfp=mean(`BV421-A`))->t150
merge(x=KD,y=t150,  all.x=T)->KD
KD %>% mutate(fc.cherry=`PE-A`/mmcherry, fc.bfp=`BV421-A`/mbfp)->KD
KD %>% dplyr::filter(time!=75)->KD
KD$plasmid<-factor(KD$plasmid, levels= c("SP430", "SP428", "SP430ABA", "SP427", "SP411")) 
KD %>% group_by(plasmid) %>% filter(time>10) %>% 
  summarize(mean.final=mean((`BV421-A`)))->mean.final
KD<-merge(KD, mean.final, all.x = T)
KD %>% group_by(plasmid) %>% filter(time==0) %>% 
  summarize(mean.init=mean((`BV421-A`)))->mean.init
KD<-merge(KD, mean.init, all.x = T)
KD %>% group_by(plasmid) %>% 
  mutate(norm.bfp=(`BV421-A`-mean.init)/(mean.final-mean.init))->KD

KD %>% filter(plasmid=="SP430ABA")->KDSP430ABA
yf=1
y0=0
KDSP430ABA$time->t
KDSP430ABA$norm.bfp ->y
fit<-nls(y ~ yf + (y0 - yf) * exp(-t*(log(2)/t1.2)),
         start = list(t1.2 = 0.8))

coef(fit)[1]
summary(fit)
SP430A.U<-fit
library(nlstools)
overview(SP430A.U)
p.val<-2 * pt(abs(6.12), 23, lower.tail = FALSE)
p<-ggplot(KDSP430ABA,aes(time,norm.bfp))+
  stat_function(fun=function(time){1-exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_KD_SP430ABA_exponential_fitting.pdf", fix, device=cairo_pdf)

p<-ggplot(KDSP430ABA,aes(time,norm.bfp))+
  stat_function(fun=function(time){1-exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(x=c(0,25),y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_segment(x=coef(fit)[1],y=-150,xend=coef(fit)[1] ,yend=0.5, linetype=3, col="black", lwd=0.5)+
  geom_segment(x=-20,y=0.5,xend=coef(fit)[1],yend= 0.5, linetype=3, lwd=0.5, col="black")
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_KD_SP430ABA_exponential_fitting_zoom.pdf", fix, device=cairo_pdf)

KD %>% filter(plasmid=="SP430")->KDSP430
yf=1
y0=0
KDSP430$time->t
KDSP430$norm.bfp ->y
fit<-nls(y ~ yf + (y0 - yf) * exp(-t*(log(2)/t1.2)),
         start = list(t1.2 = 0.8))

coef(fit)[1]
summary(fit)
SP430.U<-fit



p<-ggplot(KDSP430,aes(time,norm.bfp))+
  stat_function(fun=function(time){1-exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_KD_SP430_exponential_fitting.pdf", fix, device=cairo_pdf)

p<-ggplot(KDSP430,aes(time,norm.bfp))+
  stat_function(fun=function(time){1-exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(x=c(0,25),y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_segment(x=coef(fit)[1],y=-150,xend=coef(fit)[1] ,yend=0.5, linetype=3, lwd=0.5, col="black")+
  geom_segment(x=-20,y=0.5,xend=coef(fit)[1],yend= 0.5, linetype=3, lwd=0.5, col="black")
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_KD_SP430_exponential_fitting_zoom.pdf", fix, device=cairo_pdf)

KD %>% filter(plasmid=="SP428")->KDSP428
yf=1
y0=0
KDSP428$time->t
KDSP428$norm.bfp ->y
fit<-nls(y ~ yf + (y0 - yf) * exp(-t*(log(2)/t1.2)),
         start = list(t1.2 = 0.8))

coef(fit)[1]
summary(fit)
SP428.U<-fit



p<-ggplot(KDSP428,aes(time,norm.bfp))+
  stat_function(fun=function(time){1-exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_KD_SP428_exponential_fitting.pdf", fix, device=cairo_pdf)

p<-ggplot(KDSP428,aes(time,norm.bfp))+
  stat_function(fun=function(time){1-exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(x=c(0,25),y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_segment(x=coef(fit)[1],y=-150,xend=coef(fit)[1] ,yend=0.5, linetype=3, lwd=0.5, col="black")+
  geom_segment(x=-20,y=0.5,xend=coef(fit)[1],yend= 0.5, linetype=3, lwd=0.5, col="black")
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_KD_SP428_exponential_fitting_zoom.pdf", fix, device=cairo_pdf)

KD %>% filter(plasmid=="SP427")->KDSP427
yf=1
y0=0
KDSP427$time->t
KDSP427$norm.bfp ->y
fit<-nls(y ~ yf + (y0 - yf) * exp(-t*(log(2)/t1.2)),
         start = list(t1.2 = 0.8))

coef(fit)[1]
summary(fit)
SP427.U<-fit



p<-ggplot(KDSP427,aes(time,norm.bfp))+
  stat_function(fun=function(time){1-exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_KD_SP427_exponential_fitting.pdf", fix, device=cairo_pdf)

p<-ggplot(KDSP427,aes(time,norm.bfp))+
  stat_function(fun=function(time){1-exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(x=c(0,25),y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_segment(x=coef(fit)[1],y=-150,xend=coef(fit)[1] ,yend=0.5, linetype=3, lwd=0.5, col="black")+
  geom_segment(x=-20,y=0.5,xend=coef(fit)[1],yend= 0.5, linetype=3, lwd=0.5, col="black")
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_KD_SP427_exponential_fitting_zoom.pdf", fix, device=cairo_pdf)

KD %>% filter(plasmid=="SP411")->KDSP411
yf=1
y0=0
KDSP411$time->t
KDSP411$norm.bfp ->y
fit<-nls(y ~ yf + (y0 - yf) * exp(-t*(log(2)/t1.2)),
         start = list(t1.2 = 0.8))

coef(fit)[1]
summary(fit)
SP411.U<-fit



p<-ggplot(KDSP411,aes(time,norm.bfp))+
  stat_function(fun=function(time){1-exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_KD_SP411_exponential_fitting.pdf", fix, device=cairo_pdf)

p<-ggplot(KDSP411,aes(time,norm.bfp))+
  stat_function(fun=function(time){1-exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(x=c(0,25),y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_segment(x=coef(fit)[1],y=-150,xend=coef(fit)[1] ,yend=0.5, lty=3, lwd=0.5, col="black")+
  geom_segment(x=-20,y=0.5,xend=coef(fit)[1],yend= 0.5, lty=3, lwd=0.5, col="black")
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_KD_SP411_exponential_fitting_zoom.pdf", fix, device=cairo_pdf)


summary(SP430.U)
summary(SP430A.U)
summary(SP428.U)
summary(SP427.U)
summary(SP411.U)



