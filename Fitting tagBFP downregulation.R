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
mypal<-mypal<- pal_npg("nrc", alpha = 1)(2)
#non-fluorescent control for background extraction

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
medianexp %>% dplyr::filter(exp=="Rev")->REV
REV %>% dplyr::filter(time==0)->REV0
REV0 %>% group_by(plasmid) %>% summarize( mbfp=mean(`BV421-A`))->t0rev
merge(x=REV,y=t0rev,  all.x=T)->REV
REV %>% mutate(fc.bfp=`BV421-A`/mbfp)->REV
REV %>% dplyr::filter(time!=125)->REV
REV %>% dplyr::filter(time!=75)->REV
REV %>% group_by(plasmid) %>% filter(time>10) %>% 
  summarize(mean.final=mean((`BV421-A`)))->mean.final
REV<-merge(REV, mean.final, all.x = T)
REV %>% group_by(plasmid) %>% filter(time==0) %>% 
  summarize(mean.init=mean((`BV421-A`)))->mean.init
REV<-merge(REV, mean.init, all.x = T)
REV %>% group_by(plasmid) %>% 
  mutate(norm.bfp=(`BV421-A`-mean.final)/(mean.init-mean.final))->REV



REV %>% filter(plasmid=="SP430ABA")->REVSP430ABA
yf=0
y0=1
REVSP430ABA$time->t
REVSP430ABA$norm.bfp ->y
fit<-nls(y ~ yf + (y0 - yf) * exp(-t*(log(2)/t1.2)),
         start = list(t1.2 = 0.1))

coef(fit)[1]
summary(fit)
SP430A.D<-fit



p<-ggplot(REVSP430ABA,aes(time,norm.bfp))+
  stat_function(fun=function(time){exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_REV_SP430ABA_exponential_fitting.pdf", fix, device=cairo_pdf)

p<-ggplot(REVSP430ABA,aes(time,norm.bfp))+
  stat_function(fun=function(time){exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(x=c(0,25),y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_segment(x=coef(fit)[1],y=-150,xend=coef(fit)[1] ,yend=0.5, linetype=3, col="black", lwd=0.5)+
  geom_segment(x=-20,y=0.5,xend=coef(fit)[1],yend= 0.5, linetype=3, lwd=0.5, col="black")
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_REV_SP430ABA_exponential_fitting_zoom.pdf", fix, device=cairo_pdf)

REV %>% filter(plasmid=="SP430")->REVSP430
yf=0
y0=1
REVSP430$time->t
REVSP430$norm.bfp ->y
fit<-nls(y ~ yf + (y0 - yf) * exp(-t*(log(2)/t1.2)),
         start = list(t1.2 = 0.1))

coef(fit)[1]
summary(fit)
SP430.D<-fit



p<-ggplot(REVSP430,aes(time,norm.bfp))+
  stat_function(fun=function(time){exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_REV_SP430_exponential_fitting.pdf", fix, device=cairo_pdf)

p<-ggplot(REVSP430,aes(time,norm.bfp))+
  stat_function(fun=function(time){exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(x=c(0,25),y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_segment(x=coef(fit)[1],y=-150,xend=coef(fit)[1] ,yend=0.5, linetype=3, lwd=0.5, col="black")+
  geom_segment(x=-20,y=0.5,xend=coef(fit)[1],yend= 0.5, linetype=3, lwd=0.5, col="black")
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_REV_SP430_exponential_fitting_zoom.pdf", fix, device=cairo_pdf)

REV %>% filter(plasmid=="SP428")->REVSP428
yf=0
y0=1
REVSP428$time->t
REVSP428$norm.bfp ->y
fit<-nls(y ~ yf + (y0 - yf) * exp(-t*(log(2)/t1.2)),
         start = list(t1.2 = 0.1))

coef(fit)[1]
summary(fit)
SP428.D<-fit



p<-ggplot(REVSP428,aes(time,norm.bfp))+
  stat_function(fun=function(time){exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_REV_SP428_exponential_fitting.pdf", fix, device=cairo_pdf)

p<-ggplot(REVSP428,aes(time,norm.bfp))+
  stat_function(fun=function(time){exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(x=c(0,25),y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_segment(x=coef(fit)[1],y=-150,xend=coef(fit)[1] ,yend=0.5, linetype=3, lwd=0.5, col="black")+
  geom_segment(x=-20,y=0.5,xend=coef(fit)[1],yend= 0.5, linetype=3, lwd=0.5, col="black")
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_REV_SP428_exponential_fitting_zoom.pdf", fix, device=cairo_pdf)

REV %>% filter(plasmid=="SP427")->REVSP427
yf=0
y0=1
REVSP427$time->t
REVSP427$norm.bfp ->y
fit<-nls(y ~ yf + (y0 - yf) * exp(-t*(log(2)/t1.2)),
         start = list(t1.2 = 0.1))

coef(fit)[1]
summary(fit)
SP427.D<-fit



p<-ggplot(REVSP427,aes(time,norm.bfp))+
  stat_function(fun=function(time){exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_REV_SP427_exponential_fitting.pdf", fix, device=cairo_pdf)

p<-ggplot(REVSP427,aes(time,norm.bfp))+
  stat_function(fun=function(time){exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(x=c(0,25),y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_segment(x=coef(fit)[1],y=-150,xend=coef(fit)[1] ,yend=0.5, linetype=3, lwd=0.5, col="black")+
  geom_segment(x=-20,y=0.5,xend=coef(fit)[1],yend= 0.5, linetype=3, lwd=0.5, col="black")
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_REV_SP427_exponential_fitting_zoom.pdf", fix, device=cairo_pdf)

REV %>% filter(plasmid=="SP411")->REVSP411
yf=0
y0=1
REVSP411$time->t
REVSP411$norm.bfp ->y
fit<-nls(y ~ yf + (y0 - yf) * exp(-t*(log(2)/t1.2)),
         start = list(t1.2 = 0.1))

coef(fit)[1]
summary(fit)
SP411.D<-fit



p<-ggplot(REVSP411,aes(time,norm.bfp))+
  stat_function(fun=function(time){exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_REV_SP411_exponential_fitting.pdf", fix, device=cairo_pdf)

p<-ggplot(REVSP411,aes(time,norm.bfp))+
  stat_function(fun=function(time){exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(x=c(0,25),y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_segment(x=coef(fit)[1],y=-150,xend=coef(fit)[1] ,yend=0.5, lty=3, lwd=0.5, col="black")+
  geom_segment(x=-20,y=0.5,xend=coef(fit)[1],yend= 0.5, lty=3, lwd=0.5, col="black")
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("C:/Users/noviello/Documents/Paper 2021/Figure 4/Exponential fittings on scaled values/New_REV_SP411_exponential_fitting_zoom.pdf", fix, device=cairo_pdf)


summary(SP430.D)
summary(SP430A.D)
summary(SP428.D)
summary(SP427.D)
summary(SP411.D)
