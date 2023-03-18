#load required packages
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

#set theme for plotting
theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6, color= "black"),
                  axis.text.x = element_text(vjust = 0.5, color= "black"),
                  axis.text.y = element_text(vjust = 0.5, color= "black"),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6, color= "black"),
                  strip.background = element_blank(), legend.title = element_blank()))
mypal<- pal_npg("nrc", alpha = 1)(2)

out_path='plots'

#load non-fluorescent control for background extraction

MyFlowSet <- read.flowSet(path="fcs_files/NFC", min.limit=0.01)
#gating
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


#calculate median
medianexp<- as.data.frame(fsApply(Singlets_MyFlowSet, each_col, median))
#median bfp and and mcherry (mean of medians)
mBFP_neg<-mean(medianexp[c(1:3), 7])
mmCherry_neg<-mean(medianexp[c(1:3), 8])

#load time-course data for fitting
MyFlowSet <- read.flowSet(path="fcs_files/time-course_data", min.limit=0.01)
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

# extract normalized median expression
#1 subtract bfp and mcherry of negative control
Transf<-linearTransform(transformationId="defaultLinearTransform", a = 1, b = -mBFP_neg)
Norm_Singlets_MyFlowSet<- transform(Singlets_MyFlowSet, transformList('BV421-A' ,Transf))
Transf<-linearTransform(transformationId="defaultLinearTransform", a = 1, b = -mmCherry_neg)
Norm_Singlets_MyFlowSet<- transform(Norm_Singlets_MyFlowSet, transformList('PE-A' ,Transf))

#take median and extract
medianexp<- as.data.frame(fsApply(Norm_Singlets_MyFlowSet, each_col, median))
medianexp<-medianexp[, c(7:8)]
medianexp %>% rownames_to_column()->medianexp
medianexp %>% separate(1, c( NA, NA, "plasmid", "exp", "rep", "time", NA, NA), sep = "_") ->medianexp
medianexp$time<-as.numeric(medianexp$time)
medianexp$plasmid<-as.factor(medianexp$plasmid)

#sample "SP430" is dCas9, "SP428" is KRAB-dCas9, "SP430ABA" is Split-KRAB-dCas9, "SP411" is CasRx.
#the samples that we need are those with exp=="Rev"
medianexp %>% dplyr::filter(exp=="Rev")->REV
filter<-dplyr::filter
#Do min-max scaling of BFP level, with max as the mean for time-points after 10h and min as bfp mean at time 0
REV %>% group_by(plasmid) %>% filter(time>10) %>% 
  summarize(mean.final=mean((`BV421-A`)))->mean.final
REV<-merge(REV, mean.final, all.x = T)
REV %>% group_by(plasmid) %>% filter(time==0) %>% 
  summarize(mean.init=mean((`BV421-A`)))->mean.init
REV<-merge(REV, mean.init, all.x = T)
REV %>% group_by(plasmid) %>% 
  mutate(norm.bfp=(`BV421-A`-mean.final)/(mean.init-mean.final))->REV

#fit for KRAB-split-dCas9 for normalized (min-max scaled) data

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
half.time=data.frame(plasmid = 'SP430A', halftime = coef(fit))

p<-ggplot(REVSP430ABA,aes(time,norm.bfp))+
  stat_function(fun=function(time){exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_SP430ABA_fitting.pdf", fix, path=out_path)

#fit for dCas9 for normalized (min-max scaled) data
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
half.time=rbind(half.time,data.frame(plasmid = 'SP430', halftime = coef(fit)))

p<-ggplot(REVSP430,aes(time,norm.bfp))+
  stat_function(fun=function(time){exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_dCas9_fitting.pdf", fix, path=out_path)

#fit for KRAB-dCas9 for normalized (min-max scaled) data

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
half.time=rbind(half.time,data.frame(plasmid = 'SP428', halftime = coef(fit)))

p<-ggplot(REVSP428,aes(time,norm.bfp))+
  stat_function(fun=function(time){exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_KRAB-dCas9_fitting.pdf", fix, path=out_path)

#fit for HDAC4-dCas9 for normalized (min-max scaled) data

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
half.time=rbind(half.time,data.frame(plasmid = 'SP427', halftime = coef(fit)))

p<-ggplot(REVSP427,aes(time,norm.bfp))+
  stat_function(fun=function(time){exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_HDAC4-dCas9_fitting.pdf", fix, path=out_path)

#fit for CasRx for normalized (min-max scaled) data

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
half.time=rbind(half.time,data.frame(plasmid = 'SP411', halftime = coef(fit)))

p<-ggplot(REVSP411,aes(time,norm.bfp))+
  stat_function(fun=function(time){exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_CasRx_fitting.pdf", fix, path=out_path)

#summary of fittings

summary(SP430.D)
summary(SP430A.D)
summary(SP428.D)
summary(SP427.D)
summary(SP411.D)

write_csv(half.time,file='parameters/half_times_downregulation.csv')
