#load packages
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
library(deSolve)  
#set theme
theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6, family = "Arial"), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6, family = "Arial", color= "black"),
                  axis.text.x = element_text(vjust = 0.5, color= "black"),
                  axis.text.y = element_text(vjust = 0.5, color= "black"),
                  axis.title = element_text(size = 6, family = "Arial"), strip.text = element_text(size = 6, family = "Arial", color= "black"),
                  strip.background = element_blank(), legend.title = element_blank()))
mypal<- pal_npg("nrc", alpha = 1)(2)
#calculate background fluorescence of non-fluorescent control for background subtraction
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
# extract normalized median expression-
#background subtraction
Transf<-linearTransform(transformationId="defaultLinearTransform", a = 1, b = -mBFP_neg)
Norm_Singlets_MyFlowSet<- transform(Singlets_MyFlowSet, transformList('BV421-A' ,Transf))
Transf<-linearTransform(transformationId="defaultLinearTransform", a = 1, b = -mmCherry_neg)
Norm_Singlets_MyFlowSet<- transform(Norm_Singlets_MyFlowSet, transformList('PE-A' ,Transf))
#median
medianexp<- as.data.frame(fsApply(Norm_Singlets_MyFlowSet, each_col, median))
medianexp<-medianexp[, c(7:8)]

medianexp %>% rownames_to_column()->medianexp
medianexp %>% separate(1, c( NA, NA, "plasmid", "exp", "rep", "time", NA, NA), sep = "_") ->medianexp
medianexp$time<-as.numeric(medianexp$time)
medianexp$plasmid<-as.factor(medianexp$plasmid)
medianexp %>% group_by(plasmid) %>% fct_relevel(plasmid, levels=c("SP430", "SP428", "SP430ABA", "SP427", "SP411"))

medianexp %>% dplyr::filter (exp == "KD" )->KD
KD %>% dplyr::filter(time==0)->KD0
#calculate mean mCherry at time 0h (which is the maximal mCherry)
KD0 %>% group_by(plasmid) %>% summarize( mmcherry=mean(`PE-A`))->t0
merge(x=KD,y=t0,  all.x=T)->KD
#calculatemcherry foldchange
KD %>% mutate(fc.cherry=`PE-A`/mmcherry)->KD
#min-max scaling of bfp between time 0 (minimum) and average at time >10h (maximum)
KD %>% group_by(plasmid) %>% filter(time>10) %>% 
  summarize(mean.final=mean((`BV421-A`)))->mean.final
KD<-merge(KD, mean.final, all.x = T)
KD %>% group_by(plasmid) %>% filter(time==0) %>% 
  summarize(mean.init=mean((`BV421-A`)))->mean.init
KD<-merge(KD, mean.init, all.x = T)
KD %>% group_by(plasmid) %>% 
  mutate(norm.bfp=(`BV421-A`-mean.init)/(mean.final-mean.init))->KD


KD$plasmid<-factor(KD$plasmid, levels= c("SP430", "SP428", "SP430ABA", "SP427", "SP411")) 
#ODE-based simulation (package deSolve)
#all the parameters we need have been calculated in the previous steps. 
#mCherry degradation rate is computed in script "ODE-REV.R" and is equal to 0.036
#for KRAB-Split-dCas9 we have 
KD %>% filter(plasmid=="SP430ABA")->KDSP430ABA                   
state<-c(R=0, Y=1/0.036)
ode1<- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-  beta - R* (log(2)/t1.2)
    dY <- (K^n/(K^n + R^n))-alpha*Y
    return (list(c (dR, dY)))
  }
  )
} 
parameters<-c(beta=log(2)/2.7755, t1.2=2.7755,
              K=0.058633 , n=0.806486, alpha=0.036
)
time <- seq(0, 150, by = 0.01)
out430ABA <- ode(y=state,times=time, func=ode1, parms=parameters)
out430ABA[,3]*0.036->out430ABA[,3]
#simulate again, adding delays, in a loop
data.frame()->delays430ABA
for (t in seq(0, 25, by = 0.5)){
  time <- seq(0, 150, by = 0.01)
  delay<-time+t
  out430ABAd <- lsodes(y=state,times=delay, func=ode1, parms=parameters)
  out430ABAd[,3]<-out430ABAd[,3]*0.036
  sim<-cbind(as.data.frame(out430ABAd), t)
  delays430ABA<- rbind(delays430ABA,sim)
}
#take average mcherry at each different time-point of experimental data
KDSP430ABA[, c(4,9)] %>% group_by(time) %>% summarize(m.fc=mean(fc.cherry))->KDSP430ABA.m

#compare simulations with different delays and real data
delays430ABA->del.430ABA
del.430ABA<-merge(del.430ABA, KDSP430ABA.m, all.x = T)
del.430ABA %>% mutate(res=m.fc-Y)->del.430ABA
#calculate minimim absolute error for each delay
del.430ABA %>% group_by(t) %>% 
  summarize(N=sum(!is.na(res), na.rm=T),MAE=(sum(abs(res), na.rm=T))/(N-1)))->sum.res.430ABA
#delay that minimizes the MAE
min(sum.res.430ABA$MAE)
sum.res.430ABA$t[sum.res.430ABA$MAE==min(sum.res.430ABA$MAE)]
#MAE plot for KRAB-Split-dCas9
p<-ggplot(data=sum.res.430ABA, aes(x=t, y=MAE))+
  geom_point(size=.1, alpha=0.4, color="black") +
  geom_point(data=sum.res.430ABA, aes(x=sum.res.430ABA$t[sum.res.430ABA$MAE==min(sum.res.430ABA$MAE)], y=min(MAE)), size=.8, alpha=1, color="#E64B35FF") +
  coord_cartesian(x=c(0,25),y=c(0,0.3))+
  labs(y= "MAE" , x = "Delay (hours)")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=sum.res.430ABA, aes(y=min(MAE)), lty=2)
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("MAE_KD_KRAB-Split-dCas9_mcherry.pdf", fix, device=cairo_pdf)

#simulate again with delay found (6h)
state<-c(R=0, Y=1/0.036)
ode1<- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-  beta - R* (log(2)/t1.2)
    dY <- (K^n/(K^n + R^n))-alpha*Y
    return (list(c (dR, dY)))
  }
  )
} 
parameters<-c(beta=log(2)/2.7755, t1.2=2.7755,
              K=0.058633 , n=0.806486, alpha=0.036
)
time <- seq(0+6, 150+6, by = 0.01)
out430ABAd <- ode(y=state,times=time, func=ode1, parms=parameters)
out430ABAd[,3]*0.036->out430ABAd[,3]

#plot fitting with and without delay
p<-ggplot(KDSP430ABA,aes(time,fc.cherry))+
  geom_point(size=.8, alpha=0.4, color="#E64B35FF") +# adding connecting lines
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out430ABA), aes(time, Y))+
  geom_line(data=data.frame(out430ABAd), aes(time, Y), lty=2)
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("mCherry_KRAB-Split-dCas9.pdf", fix, device=cairo_pdf)

p<-ggplot(KDSP430ABA,aes(time,norm.bfp))+
  geom_point(size=.8, alpha=0.4, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "tagBFP")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out430ABA), aes(time, R))
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("tagBFP_KRAB-Split-dCas9.pdf", fix, device=cairo_pdf)


#Same for HDAC4-dCas9
KD %>% filter(plasmid=="SP427")->KDSP427
state<-c(R=0, Y=1/0.036)
ode1<- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-  beta - R* (log(2)/t1.2)
    dY <- (K^n/(K^n + R^n))-alpha*Y
    return (list(c (dR, dY)))
  }
  )
} 
parameters<-c(beta=log(2)/8.395, t1.2=8.395,
              K=0.417969 , n=3.181102, alpha=0.036
)
time <- seq(0, 150, by = 0.01)
out427 <- ode(y=state,times=time, func=ode1, parms=parameters)
out427[,3]*0.036->out427[,3]

data.frame()->delays427
for (t in seq(0, 25, by = 0.5)){
  time <- seq(0, 150, by = 0.01)
  delay<-time+t
  out427d <- lsodes(y=state,times=delay, func=ode1, parms=parameters)
  out427d[,3]<-out427d[,3]*0.036
  sim<-cbind(as.data.frame(out427d), t)
  delays427<- rbind(delays427,sim)
}
KDSP427[, c(4,9)] %>% group_by(time) %>% summarize(m.fc=mean(fc.cherry))->KDSP427.m
delays427->del.427
del.427<-merge(del.427, KDSP427.m, all.x = T)
del.427 %>% mutate(res=m.fc-Y)->del.427
del.427 %>% group_by(t) %>% 
  summarize(N=sum(!is.na(res), na.rm=T),MAE=sum(abs(res), na.rm=T)/(N-1)))->sum.res.427
min(sum.res.427$MAE)
sum.res.427$t[sum.res.427$MAE==min(sum.res.427$MAE)]

p<-ggplot(data=sum.res.427, aes(x=t, y=MAE))+
  geom_point(size=.1, alpha=0.4, color="black") +
  geom_point(data=sum.res.427, aes(x=sum.res.427$t[sum.res.427$MAE==min(sum.res.427$MAE)], y=min(MAE)), size=.8, alpha=1, color="#E64B35FF") +
  coord_cartesian(x=c(0,25),y=c(0,0.4))+
  labs(y= "MAE" , x = "Delay (hours)")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=sum.res.427, aes(y=min(MAE)), lty=2)

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("MAE_KD_HDAC4-dCas9_mcherry.pdf", fix, device=cairo_pdf)

p<-ggplot(KDSP427,aes(time,fc.cherry))+
  geom_point(size=.8, alpha=0.4, color="#E64B35FF") +# adding connecting lines
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out427), aes(time, Y))
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("mCherry_HDAC4_dCas9.pdf", fix, device=cairo_pdf)

p<-ggplot(KDSP427,aes(time,norm.bfp))+
  geom_point(size=.8, alpha=0.4, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "tagBFP")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out427), aes(time, R))
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("tagBFP_HDAC4-dCas9.pdf", fix, device=cairo_pdf)

#Same for KRAB-dCas9
KD %>% filter(plasmid=="SP428")->KDSP428
state<-c(R=0, Y=1/0.036)
ode1<- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-  beta - R* (log(2)/t1.2)
    dY <- (K^n/(K^n + R^n))-alpha*Y
    return (list(c (dR, dY)))
  }
  )
} 
parameters<-c(beta=log(2)/4.8674, t1.2=4.8674,
              K=0.319691  , n=2.210463, alpha=0.036
)
time <- seq(0, 150, by = 0.01)
out428 <- ode(y=state,times=time, func=ode1, parms=parameters)
out428[,3]*0.036->out428[,3]

data.frame()->delays428
for (t in seq(0, 25, by = 0.5)){
  time <- seq(0, 150, by = 0.01)
  delay<-time+t
  out428d <- lsodes(y=state,times=delay, func=ode1, parms=parameters)
  out428d[,3]<-out428d[,3]*0.036
  sim<-cbind(as.data.frame(out428d), t)
  delays428<- rbind(delays428,sim)
}
delays428->del.428
KDSP428[, c(4,9)] %>% group_by(time) %>% summarize(m.fc=mean(fc.cherry))->KDSP428.m
del.428<-merge(del.428, KDSP428.m, all.x = T)
del.428 %>% mutate(res=m.fc-Y)->del.428
del.428 %>% group_by(t) %>% 
  summarize(N=sum(!is.na(res), na.rm=T),MAE=(sum(abs(res), na.rm=T))/(N-1)))->sum.res.428
min(sum.res.428$MAE)
sum.res.428$t[sum.res.428$MAE==min(sum.res.428$MAE)]

p<-ggplot(data=sum.res.428, aes(x=t, y=MAE))+
  geom_point(size=.1, alpha=0.4, color="black") +
  geom_point(data=sum.res.428, aes(x=sum.res.428$t[sum.res.428$MAE==min(sum.res.428$MAE)], y=min(MAE)), size=.8, alpha=1, color="#E64B35FF") +
  coord_cartesian(x=c(0,25),y=c(0,0.3))+
  labs(y= "MAE" , x = "Delay (hours)")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=sum.res.428, aes(y=min(MAE)), lty=2)

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("MAE_KD_KRAB-dCas9_mcherry.pdf", fix, device=cairo_pdf)

ode1<- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-  beta - R* (log(2)/t1.2)
    dY <- (K^n/(K^n + R^n))-alpha*Y
    return (list(c (dR, dY)))
  }
  )
} 
parameters<-c(beta=log(2)/4.8674, t1.2=4.8674,
              K=0.319691  , n=2.210463, alpha=0.036
)
time <- seq(0+3, 150+3, by = 0.01)
out428d <- ode(y=state,times=time, func=ode1, parms=parameters)
out428d[,3]*0.036->out428d[,3]

p<-ggplot(KDSP428,aes(time,fc.cherry))+
  geom_point(size=.8, alpha=0.4, color="#E64B35FF") +# adding connecting lines
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out428), aes(time, Y))+
  geom_line(data=data.frame(out428d), aes(time, Y), lty=2)
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("mCherry_KRAB-dCas9.pdf", fix, device=cairo_pdf)

p<-ggplot(KDSP428,aes(time,fc.bfp))+
  geom_point(size=.8, alpha=0.4, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "tagBFP")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out428), aes(time, R))

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("tagBFP_KRAB-dCas9.pdf", fix, device=cairo_pdf)

#Same for CasRx
KD %>% filter(plasmid=="SP411")->KDSP411
state<-c(R=0, Y=1/0.036)
ode1<- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-  beta - R* (log(2)/t1.2)
    dY <- (K^n/(K^n + R^n))-alpha*Y
    return (list(c (dR, dY)))
  }
  )
} 
parameters<-c(beta=log(2)/7.168, t1.2=7.168,
              K=0.37203  , n=1.04901, alpha=0.036
)

time <- seq(0, 150, by = 0.01)
out411 <- ode(y=state,times=time, func=ode1, parms=parameters)
out411[,3]*0.036->out411[,3]

data.frame()->delays411
for (t in seq(0, 25, by = 0.5)){
  time <- seq(0, 150, by = 0.01)
  delay<-time+t
  out411d <- lsodes(y=state,times=delay, func=ode1, parms=parameters)
  out411d[,3]<-out411d[,3]*0.036
  sim<-cbind(as.data.frame(out411d), t)
  delays411<- rbind(delays411,sim)
}
delays411 %>% 
  filter()->del.411

KDSP411[, c(4,9)] %>% group_by(time) %>% summarize(m.fc=mean(fc.cherry))->KDSP411.m
del.411<-merge(del.411, KDSP411.m, all.x = T)
del.411 %>% mutate(res=m.fc-Y)->del.411
del.411 %>% group_by(t) %>% 
  summarize(N=sum(!is.na(res), na.rm=T),MAE=(sum(abs(res), na.rm=T))/(N-1)))->sum.res.411
min(sum.res.411$MAE)
sum.res.411$t[sum.res.411$MAE==min(sum.res.411$MAE)]

p<-ggplot(data=sum.res.411, aes(x=t, y=MAE))+
  geom_point(size=.1, alpha=0.4, color="black") +
  geom_point(data=sum.res.411, aes(x=sum.res.411$t[sum.res.411$MAE==min(sum.res.411$MAE)], y=min(MAE)), size=.8, alpha=1, color="#E64B35FF") +
  coord_cartesian(x=c(0,25),y=c(0,0.3))+
  labs(y= "MAE" , x = "Delay (hours)")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=sum.res.411, aes(y=min(MAE)), lty=2)

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("MAE_KD_CasRx_mcherry.pdf", fix, device=cairo_pdf)

p<-ggplot(KDSP411,aes(time,fc.cherry))+
  geom_point(size=.8, alpha=0.4, color="#E64B35FF") +# adding connecting lines
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out411), aes(time, Y))
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("mCherry_CasRx.pdf", fix, device=cairo_pdf)

p<-ggplot(KDSP411,aes(time,norm.bfp))+
  geom_point(size=.8, alpha=0.4, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "tagBFP")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out411), aes(time, R))
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("tagBFP_CasRx.pdf", fix, device=cairo_pdf)

#Same for dCas9
KD %>% filter(plasmid=="SP430")->KDSP430
summarysp430
state<-c(R=0, Y=1/0.036)
ode1<- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-  beta - R* (log(2)/t1.2)
    dY <- (K^n/(K^n + R^n))-alpha*Y
    return (list(c (dR, dY)))
  }
  )
} 
parameters<-c(beta=log(2)/2.3352, t1.2=2.3352,
              K=0.75431  , n=1.00415, alpha=0.036
)

time <- seq(0, 150, by = 0.01)
out430 <- ode(y=state,times=time, func=ode1, parms=parameters)
out430[,3]*0.036->out430[,3]

data.frame()->delays430
for (t in seq(0, 25, by = 0.5)){
  time <- seq(0, 150, by = 0.01)
  delay<-time+t
  out430d <- lsodes(y=state,times=delay, func=ode1, parms=parameters)
  out430d[,3]<-out430d[,3]*0.036
  sim<-cbind(as.data.frame(out430d), t)
  delays430<- rbind(delays430,sim)
}
KDSP430[, c(4,9)] %>% group_by(time) %>% summarize(m.fc=mean(fc.cherry))->KDSP430.m
delays430->del.430
del.430<-merge(del.430, KDSP430.m, all.x = T)
del.430 %>% mutate(res=m.fc-Y)->del.430
del.430 %>% group_by(t) %>% 
  summarize(N=sum(!is.na(res), na.rm=T),MAE=(sum(abs(res), na.rm=T))/(N-1)))->sum.res.430
min(sum.res.430$MAE)
sum.res.430$t[sum.res.430$MAE==min(sum.res.430$MAE)]

p<-ggplot(data=sum.res.430, aes(x=t, y=MAE))+
  geom_point(size=.1, alpha=0.4, color="black") +
  geom_point(data=sum.res.430, aes(x=sum.res.430ABA$t[sum.res.430$MAE==min(sum.res.430$MAE)], y=min(MAE)), size=.8, alpha=1, color="#E64B35FF") +
  coord_cartesian(x=c(0,25),y=c(0,0.3))+
  labs(y= "MAE" , x = "Delay (hours)")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=sum.res.430, aes(y=min(MAE)), lty=2)
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("MAE_KD_dCas9_mcherry.pdf", fix, device=cairo_pdf)

p<-ggplot(KDSP430,aes(time,fc.cherry))+
  geom_point(size=.8, alpha=0.4, color="#E64B35FF") +# adding connecting lines
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out430), aes(time, Y))
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("mCherry_dCas9.pdf", fix, device=cairo_pdf)

p<-ggplot(KDSP430,aes(time,norm.bfp))+
  geom_point(size=.8, alpha=0.4, color="#4DBBD5FF") +# adding connecting lines
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "tagBFP")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out430), aes(time, R))
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("tagBFP_dCas9.pdf", fix, device=cairo_pdf)

