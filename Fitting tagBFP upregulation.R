#load libraries
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
library(nlstools)
#set theme
theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6, color= "black"),
                  axis.text.x = element_text(vjust = 0.5, color= "black"),
                  axis.text.y = element_text(vjust = 0.5, color= "black"),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6, color= "black"),
                  strip.background = element_blank(), legend.title = element_blank()))

out_path='plots'

#load non-fluorescent control for background-subtraction
MyFlowSet <- read.flowSet(path="fcs_files/NFC", min.limit=0.01)
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
#bfp and mcherry of NFC for background subtraction
mBFP_neg<-mean(medianexp[c(1:3), 7])
mmCherry_neg<-mean(medianexp[c(1:3), 8])

#load time-course data
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
#remove background fluorescence and calculate median
Transf<-linearTransform(transformationId="defaultLinearTransform", a = 1, b = -mBFP_neg)
Norm_Singlets_MyFlowSet<- transform(Singlets_MyFlowSet, transformList('BV421-A' ,Transf))
Transf<-linearTransform(transformationId="defaultLinearTransform", a = 1, b = -mmCherry_neg)
Norm_Singlets_MyFlowSet<- transform(Norm_Singlets_MyFlowSet, transformList('PE-A' ,Transf))
medianexp<- as.data.frame(fsApply(Norm_Singlets_MyFlowSet, each_col, median))
medianexp<-medianexp[, c(7:8)]

medianexp %>% rownames_to_column()->medianexp
medianexp %>% separate(1, c( NA, NA, "plasmid", "exp", "rep", "time", NA, NA), sep = "_") ->medianexp
medianexp$time<-as.numeric(medianexp$time)
medianexp$plasmid<-as.factor(medianexp$plasmid)

#the data that we need have exp== "KD"
medianexp %>% dplyr::filter (exp == "KD" )->KD
KD$plasmid<-factor(KD$plasmid, levels= c("SP430", "SP428", "SP430ABA", "SP427", "SP411")) 

#min-max scaling of bfp between time 0 (max) and time >10h (min)
KD %>% group_by(plasmid) %>% filter(time>10) %>% 
  summarize(mean.final=mean((`BV421-A`)))->mean.final
KD<-merge(KD, mean.final, all.x = T)
KD %>% group_by(plasmid) %>% filter(time==0) %>% 
  summarize(mean.init=mean((`BV421-A`)))->mean.init
KD<-merge(KD, mean.init, all.x = T)
KD %>% group_by(plasmid) %>% 
  mutate(norm.bfp=(`BV421-A`-mean.init)/(mean.final-mean.init))->KD



#fitting of scaled values for KRAB-Split-dCas9
KD %>% filter(plasmid=="SP430ABA")->KDSP430ABA
yf=1
y0=0
KDSP430ABA$time->t
KDSP430ABA$norm.bfp ->y
fit<-nls(y ~ yf + (y0 - yf) * exp(-t*(log(2)/t1.2)),
         start = list(t1.2 = 0.8))
SP430A.U<-fit
half.time=data.frame(plasmid = 'SP430A', halftime = coef(fit))

p<-ggplot(KDSP430ABA,aes(time,norm.bfp))+
  stat_function(fun=function(time){1-exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("KD_KRAB-Split-dCas9_fitting.pdf",fix,path=out_path)

#fitting of scaled values for dCas9
KD %>% filter(plasmid=="SP430")->KDSP430
yf=1
y0=0
KDSP430$time->t
KDSP430$norm.bfp ->y
fit<-nls(y ~ yf + (y0 - yf) * exp(-t*(log(2)/t1.2)),
         start = list(t1.2 = 0.8))
SP430.U<-fit
half.time=rbind(half.time,c('SP430',coef(fit)))

p<-ggplot(KDSP430,aes(time,norm.bfp))+
  stat_function(fun=function(time){1-exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("KD_dCas9_fitting.pdf", fix,path=out_path)

#fitting of scaled values for KRAB-dCas9

KD %>% filter(plasmid=="SP428")->KDSP428
yf=1
y0=0
KDSP428$time->t
KDSP428$norm.bfp ->y
fit<-nls(y ~ yf + (y0 - yf) * exp(-t*(log(2)/t1.2)),
         start = list(t1.2 = 0.8))
SP428.U<-fit
half.time=rbind(half.time,c('SP428',coef(fit)))

p<-ggplot(KDSP428,aes(time,norm.bfp))+
  stat_function(fun=function(time){1-exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("KD_KRAB-dCas9_fitting.pdf", fix, path=out_path)

#fitting of scaled values for HDAC4-dCas9

KD %>% filter(plasmid=="SP427")->KDSP427
yf=1
y0=0
KDSP427$time->t
KDSP427$norm.bfp ->y
fit<-nls(y ~ yf + (y0 - yf) * exp(-t*(log(2)/t1.2)),
         start = list(t1.2 = 0.8))
SP427.U<-fit
half.time=rbind(half.time,c('SP427',coef(fit)))

p<-ggplot(KDSP427,aes(time,norm.bfp))+
  stat_function(fun=function(time){1-exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("KD_HDAC4-dCas9_fitting.pdf", fix, path=out_path)

#fitting of scaled values for CasRx
KD %>% filter(plasmid=="SP411")->KDSP411
yf=1
y0=0
KDSP411$time->t
KDSP411$norm.bfp ->y
fit<-nls(y ~ yf + (y0 - yf) * exp(-t*(log(2)/t1.2)),
         start = list(t1.2 = 0.8))

SP411.U<-fit
half.time=rbind(half.time,c('SP411',coef(fit)))

p<-ggplot(KDSP411,aes(time,norm.bfp))+
  stat_function(fun=function(time){1-exp(-time*(log(2)/coef(fit)[[1]]))}, color="black")+
  geom_point(size=.4, alpha=0.7, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(y=c(-0.15,1.2))+
  scale_y_continuous(breaks=(c(0,.25,.5,.75,1)))+
  labs(x= "Time (hours)" , y = "tagBFP (% of final)")+
  scale_color_npg()+
  scale_fill_npg()

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("KD_CasRx_fitting.pdf", fix, path=out_path)

#summary of fittings
summary(SP430.U)
summary(SP430A.U)
summary(SP428.U)
summary(SP427.U)
summary(SP411.U)

write_csv(half.time,file='parameters/half_times_upregulation.csv')

