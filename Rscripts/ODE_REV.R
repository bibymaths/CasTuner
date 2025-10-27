#load required libraries
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
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6, color= "black"),
                  axis.text.x = element_text(vjust = 0.5, color= "black"),
                  axis.text.y = element_text(vjust = 0.5, color= "black"),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6, color= "black"),
                  strip.background = element_blank(), legend.title = element_blank()))

out_path='plots'

#load non-fluorescent control for background fluorescence subtraction
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

#load time-course data, remove background fluorescence, take median
MyFlowSet <- read.flowSet(path= "CasTuner-Python/fcs_files/time-course_data", min.limit=0.01)
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

#calculate maximal mCherry (which is mCherry in CasRx experiment at last time point -150h-) for fold change calculation
medianexp %>% dplyr::filter(exp=="Rev")->REV
REV %>% dplyr::filter(time==150)->REV150
filter<-dplyr::filter
REV150 %>% filter(plasmid=="SP411") %>% summarize( mmcherry=mean(`PE-A`))->t0rev2
merge(x=REV,y=t0rev2,  all.x=T)->REV
REV %>% mutate(fc.cherry=`PE-A`/mmcherry)->REV

#mean-max scaling of bfp between time 0 and average at time >10h.
REV %>% group_by(plasmid) %>% filter(time>10) %>% 
  summarize(mean.final=mean((`BV421-A`)))->mean.final
REV<-merge(REV, mean.final, all.x = T)
REV %>% group_by(plasmid) %>% filter(time==0) %>% 
  summarize(mean.init=mean((`BV421-A`)))->mean.init
REV<-merge(REV, mean.init, all.x = T)
REV %>% group_by(plasmid) %>% 
  mutate(norm.bfp=(`BV421-A`-mean.final)/(mean.init-mean.final))->REV

#For this model we assume that for CasRx the repressor goes to 0 very rapidly (indeed the t1.2 is  <0.3h)
#in the case of CasRx, the reporter restores its original level (there is no permanent repression)
#we use CasRx derepression experiment to calculate mCherry half-life

#2 methods work the same
#method 1
REV %>% filter(plasmid=="SP411")->REVSP411
yf=1
y0=mean(REVSP411$fc.cherry[REVSP411$time==0])
REVSP411$time->t
REVSP411$fc.cherry ->y
fit<-nls(y ~ yf + (y0 - yf) * exp(-t*(log(2)/t1.2)),
         start = list(t1.2 = 0.1))
coef(fit)[1]
summary(fit)

mcherry.halflife<-fit

mc<-as.numeric(log(2)/coef(fit)[1])

# 0.0360691
#method 2
fit<-nls(y ~ yf + (y0 - yf) * exp(-t*(a)),
         start = list(a = 0.1))
coef(fit)[1]
summary(fit)
alpha=round(mc,3)
#same result: std.error 0.003 is very small
write_csv(data.frame(alpha),file='parameters/alphamcherry.csv')


# Loading parameters estimated from tBFP time courses and dose-response curves
dplyr::rename->rename
t_down = read.csv('parameters/half_times_downregulation.csv') %>% select(-se)%>% rename(t_down=halftime)
t_up = read.csv('parameters/half_times_upregulation.csv') %>% select(-se) %>% rename(t_up=halftime)
hill_pars = read.csv('parameters/Hill_parameters.csv') 
alpha = read.csv('parameters/alphamcherry.csv')
all_par = t_down %>% left_join(t_up) %>% left_join(hill_pars) %>%merge(alpha)

#Now we simulate the ODE model for CasRx

REV %>% filter(plasmid=="SP411")->REVSP411
sel_par = all_par %>% filter(plasmid=='SP411')
mean(REVSP411$norm.bfp[REVSP411$time==0])->R
mean(REVSP411$fc.cherry[REVSP411$time==0])->Y


state<-c(R=1, Y=Y/sel_par$alpha)
ode1<- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-  - R* (log(2)/t1.2) #in this way this tend to 0 with a dynamic dependent from t1.2 (fitted)
    dY <- (K^n/(K^n +R^n))-alpha*Y
    return (list(c (dR, dY)))
  }
  )
} 
parameters = c(t1.2=sel_par$t_down, K=sel_par$K, n=sel_par$n, alpha=sel_par$alpha)

time <- seq(0, 150, by = 0.01)
out411 <- ode(y=state,times=time, func=ode1, parms=parameters)
out411[,3]*sel_par$alpha->out411[,3]

#Now we simulate the model for dCas9
REV %>% filter(plasmid=="SP430")->REVSP430
sel_par = all_par %>% filter(plasmid=='SP430')
mean(REVSP430$norm.bfp[REVSP430$time==0])->R
mean(REVSP430$fc.cherry[REVSP430$time==0])->Y

state<-c(R=R, Y=Y/sel_par$alpha)
ode1<- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-  - R* (log(2)/t1.2)
    dY <- (K^n)/(K^n+R^n)-alpha*Y
    return (list(c (dR, dY)))
  }
  )
}
parameters = c(t1.2=sel_par$t_down, K=sel_par$K, n=sel_par$n, alpha=sel_par$alpha)


time <- seq(0, 150, by = 0.01)
out430 <- ode(y=state,times=time, func=ode1, parms=parameters)
out430[,3]<-out430[,3]*sel_par$alpha
p<-ggplot(REVSP430,aes(time,fc.cherry))+
  geom_point(size=.6, alpha=0.4, color="#E64B35FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out430), aes(time, Y))
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_mCherry_dCas9_Hill.pdf", fix, path=out_path)

p<-ggplot(REVSP430,aes(time,norm.bfp))+
  geom_point(size=.6, alpha=0.4, color="#4DBBD5FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out430), aes(time, R))
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_tagBFP_dCas9_Hill.pdf", fix, path=out_path)


#Now we simulate the model for KRAB-Split-dCas9
REV %>% filter(plasmid=="SP430ABA")->REVSP430ABA
sel_par = all_par %>% filter(plasmid=='SP430A')
mean(REVSP430ABA$norm.bfp[REVSP430ABA$time==0])->R
mean(REVSP430ABA$fc.cherry[REVSP430ABA$time==0])->Y
state<-c(R=1, Y=Y/sel_par$alpha)
ode1<- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-  - R* (log(2)/t1.2)
    dY <- (K^n)/(K^n+R^n)-alpha*Y
    return (list(c (dR, dY)))
  }
  )
} 

parameters = c(t1.2=sel_par$t_down, K=sel_par$K, n=sel_par$n, alpha=sel_par$alpha)

time <- seq(0, 150, by = 0.01)
out430ABA <- ode(y=state,times=time, func=ode1, parms=parameters)
diagnostics(out430ABA)

#note that "ode" here is inaccurate, so I use another computing method
out430ABA <- lsodes(y=state,times=time, func=ode1, parms=parameters)
out430ABA[,3]<-out430ABA[,3]*sel_par$alpha
p<-ggplot(REVSP430ABA,aes(time,fc.cherry))+
  geom_point(size=.6, alpha=0.4, color="#E64B35FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out430ABA), aes(time, Y))
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_mCherry_KRAB-Split-dCas9_Hill.pdf", fix, path=out_path)

p<-ggplot(REVSP430ABA,aes(time,norm.bfp))+
  geom_point(size=.6, alpha=0.4, color="#4DBBD5FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out430ABA), aes(time, R))
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_tagBFP_KRAB-Split-dCas9_Hill.pdf", fix, path=out_path)

#Now we simulate the model for KRAB-dCas9
REV %>% filter(plasmid=="SP428")->REVSP428
sel_par = all_par %>% filter(plasmid=='SP428')
mean(REVSP428$norm.bfp[REVSP428$time==0])->R
mean(REVSP428$fc.cherry[REVSP428$time==0])->Y
state<-c(R=1, Y=Y/sel_par$alpha)
ode1<- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-  - R* (log(2)/t1.2)
    dY <- (K^n)/(K^n+R^n)-alpha*Y
    return (list(c (dR, dY)))
  }
  )
} 
parameters = c(t1.2=sel_par$t_down, K=sel_par$K, n=sel_par$n, alpha=sel_par$alpha)

time <- seq(0, 150, by = 0.01)
out428 <- ode(y=state,times=time, func=ode1, parms=parameters)

out428[,3]<-out428[,3]*sel_par$alpha
p<-ggplot(REVSP428,aes(time,fc.cherry))+
  geom_point(size=.6, alpha=0.4, color="#E64B35FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out428), aes(time, Y))
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_mCherry_KRAB-dCas9_Hill.pdf", fix, path=out_path)

p<-ggplot(REVSP428,aes(time,norm.bfp))+
  geom_point(size=.6, alpha=0.4, color="#4DBBD5FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out428), aes(time, R))
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_tagBFP_dCas9_Hill.pdf", fix, path=out_path)

#Now we simulate the model for HDAC4-dCas9
REV %>% filter(plasmid=="SP427")->REVSP427
sel_par = all_par %>% filter(plasmid=='SP427')
mean(REVSP427$norm.bfp[REVSP427$time==0])->R
mean(REVSP427$fc.cherry[REVSP427$time==0])->Y
state<-c(R=1, Y=Y/sel_par$alpha)
ode1<- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-  - R* (log(2)/t1.2)
    dY <- (K^n)/(K^n+R^n)-alpha*Y
    return (list(c (dR, dY)))
  }
  )
} 

parameters = c(t1.2=sel_par$t_down, K=sel_par$K, n=sel_par$n, alpha=sel_par$alpha)

time <- seq(0, 150, by = 0.01)
out427 <- ode(y=state,times=time, func=ode1, parms=parameters)

out427[,3]<-out427[,3]*sel_par$alpha
p<-ggplot(REVSP427,aes(time,fc.cherry))+
  geom_point(size=.6, alpha=0.4, color="#E64B35FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out427), aes(time, Y))
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_mCherry_HDAC4-dCas9_Hill.pdf", fix, path=out_path)

p<-ggplot(REVSP427,aes(time,norm.bfp))+
  geom_point(size=.6, alpha=0.4, color="#4DBBD5FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out427), aes(time, R))
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_tagBFP_HDAC4-dCas9_Hill.pdf", fix, path=out_path)

#Simulate the models again but adding delays to find which model (and therefore delay) fits best

#for Krab-Split-dCas9:
REV %>% filter(plasmid=="SP430ABA")->REVSP430ABA
sel_par = all_par %>% filter(plasmid=='SP430A')
mean(REVSP430ABA$norm.bfp[REVSP430ABA$time==0])->R
mean(REVSP430ABA$fc.cherry[REVSP430ABA$time==0])->Y
state<-c(R=R, Y=Y/sel_par$alpha)
ode1<- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-   - R* (log(2)/t1.2)
    dY <- (K^n/(K^n+R^n))-(alpha)*Y
    return (list(c (dR, dY)))
  }
  )
} 

parameters = c(t1.2=sel_par$t_down, K=sel_par$K, n=sel_par$n, alpha=sel_par$alpha)


time <- seq(0, 150, by = 0.01)
out430ABA <- lsodes(y=state,times=time, func=ode1, parms=parameters)
out430ABA[,3]<-out430ABA[,3]*sel_par$alpha
data.frame()->delays430ABA
for (t in seq(0, 25, by = 0.5)){
  time <- seq(0, 150, by = 0.01)
  delay<-time+t
  out430ABAd <- lsodes(y=state,times=delay, func=ode1, parms=parameters)
  out430ABAd[,3]<-out430ABAd[,3]*sel_par$alpha
  sim<-cbind(as.data.frame(out430ABAd), t)
  delays430ABA<- rbind(delays430ABA,sim)
}
delays430ABA %>% 
  filter()->del.430ABA 
REVSP430ABA[, c(4:9)] %>% group_by(time) %>% summarize(m.fc=mean(fc.cherry))->REVSP430ABA.m
del.430ABA<-merge(del.430ABA, REVSP430ABA.m, all.x = T)
del.430ABA %>% mutate(res=m.fc-Y)->del.430ABA

#finding the delay that minimize the standard error of the estimate using minimal absolute error (MAE)
del.430ABA %>% group_by(t) %>% 
  summarize(N=sum(!is.na(res), na.rm=T),MAE=sum(abs(res), na.rm=T)/(N-1))->sum.res.430ABA

sum.res.430ABA$t[sum.res.430ABA$MAE==min(sum.res.430ABA$MAE)]->d430a
#16h
delaysREV<-data.frame(plasmid="SP430A", d_rev=d430a)

p<-ggplot(data=sum.res.430ABA, aes(x=t, y=MAE))+
  geom_point(size=.1, alpha=0.4, color="black") +
  geom_point(data=sum.res.430ABA, aes(x=sum.res.430ABA$t[sum.res.430ABA$MAE==min(sum.res.430ABA$MAE)], y=min(MAE)), size=.8, alpha=1, color="red") +
  coord_cartesian(x=c(0,25),y=c(0,0.3))+
  labs(y= "MAE" , x = "Delay (hours)")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=sum.res.430ABA, aes(y=min(MAE)), lty=2)

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("MAE_REV_KRAB-Split-dCas9_mcherry.pdf", fix, path=out_path)

#Same analysis, for KRAB-dCas9
mean(REVSP428$norm.bfp[REVSP428$time==0])->R
mean(REVSP428$fc.cherry[REVSP428$time==0])->Y
state<-c(R=1, Y=Y/sel_par$alpha)
ode1<- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-  - R* (log(2)/t1.2)
    dY <- (K^n)/(K^n+R^n)-alpha*Y
    return (list(c (dR, dY)))
  }
  )
} 
sel_par = all_par %>% filter(plasmid=='SP428')
parameters = c(t1.2=sel_par$t_down, K=sel_par$K, n=sel_par$n, alpha=sel_par$alpha)

time <- seq(0, 150, by = 0.01)
out428 <- ode(y=state,times=time, func=ode1, parms=parameters)

data.frame()->delays428
for (t in seq(0, 25, by = 0.5)){
  time <- seq(0, 150, by = 0.01)
  delay<-time+t
  out428d <- lsodes(y=state,times=delay, func=ode1, parms=parameters)
  out428d[,3]<-out428d[,3]*sel_par$alpha
  sim<-cbind(as.data.frame(out428d), t)
  delays428<- rbind(delays428,sim)
}
delays428 %>% 
  filter()->del.428

REVSP428[, c(4:9)] %>% group_by(time) %>% summarize(m.fc=mean(fc.cherry))->REVSP428.m
del.428<-merge(del.428, REVSP428.m, all.x = T)
del.428 %>% mutate(res=m.fc-Y)->del.428
del.428 %>% group_by(t) %>% 
  summarize(N=sum(!is.na(res), na.rm=T),MAE=(sum(abs(res), na.rm=T))/(N-1))->sum.res.428
min(sum.res.428$MAE)
sum.res.428$t[sum.res.428$MAE==min(sum.res.428$MAE)]->d428

#18.5h
delaysREV<-rbind(delaysREV,data.frame(plasmid="SP428", d_rev=d428))

p<-ggplot(data=sum.res.428, aes(x=t, y=MAE))+
  geom_point(size=.1, alpha=0.4, color="black") +
  geom_point(data=sum.res.428, aes(x=sum.res.428$t[sum.res.428$MAE==min(sum.res.428$MAE)], y=min(MAE)), size=.8, alpha=1, color="red") +
  coord_cartesian(x=c(0,25),y=c(0,0.3))+
  labs(y= "MAE" , x = "Delay (hours)")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=sum.res.428, aes(y=min(MAE)), lty=2)

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("MAE_REV_KRAB-dCas9_mcherry.pdf", fix, path=out_path)


#Same analysis for HDAC4-dCas9
mean(REVSP427$norm.bfp[REVSP427$time==0])->R
mean(REVSP427$fc.cherry[REVSP427$time==0])->Y
state<-c(R=1, Y=Y/sel_par$alpha)

sel_par = all_par %>% filter(plasmid=='SP427')
parameters = c(t1.2=sel_par$t_down, K=sel_par$K, n=sel_par$n, alpha=sel_par$alpha)

time <- seq(0, 150, by = 0.01)
data.frame()->delays427
for (t in seq(0, 25, by = 0.5)){
  time <- seq(0, 150, by = 0.01)
  delay<-time+t
  out427d <- lsodes(y=state,times=delay, func=ode1, parms=parameters)
  out427d[,3]<-out427d[,3]*sel_par$alpha
  sim<-cbind(as.data.frame(out427d), t)
  delays427<- rbind(delays427,sim)
}
REVSP427[, c(4:9)] %>% group_by(time) %>% summarize(m.fc=mean(fc.cherry))->REVSP427.m
delays427 %>% 
  filter()->del.427
del.427<-merge(del.427, REVSP427.m, all.x = T)
del.427 %>% mutate(res=m.fc-Y)->del.427
del.427 %>% group_by(t) %>% 
  summarize(N=sum(!is.na(res), na.rm=T),MAE=sum(abs(res), na.rm=T)/(N-1))->sum.res.427
min(sum.res.427$MAE)
sum.res.427$t[sum.res.427$MAE==min(sum.res.427$MAE)] ->d427
#6h
delaysREV<-rbind(delaysREV,data.frame(plasmid="SP427", d_rev=d427))

p<-ggplot(data=sum.res.427, aes(x=t, y=MAE))+
  geom_point(size=.1, alpha=0.4, color="black") +
  geom_point(data=sum.res.427, aes(x=sum.res.427$t[sum.res.427$MAE==min(sum.res.427$MAE)], y=min(MAE)), size=.8, alpha=1, color="red") +
  coord_cartesian(x=c(0,25),y=c(0,0.3))+
  labs(y= "MAE" , x = "Delay (hours)")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=sum.res.427, aes(y=min(MAE)), lty=2)

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("MAE_REV_HDAC4-dCas9_mcherry.pdf", fix, path=out_path)

#Same analysis for dCas9
mean(REVSP430$norm.bfp[REVSP430$time==0])->R
mean(REVSP430$fc.cherry[REVSP430$time==0])->Y
state<-c(R=1, Y=Y/sel_par$alpha)

sel_par = all_par %>% filter(plasmid=='SP430')
parameters = c(t1.2=sel_par$t_down, K=sel_par$K, n=sel_par$n, alpha=sel_par$alpha)

time <- seq(0, 150, by = 0.01)
data.frame()->delays430
for (t in seq(0, 25, by = 0.5)){
  time <- seq(0, 150, by = 0.01)
  delay<-time+t
  out430d <- lsodes(y=state,times=delay, func=ode1, parms=parameters)
  out430d[,3]<-out430d[,3]*sel_par$alpha
  sim<-cbind(as.data.frame(out430d), t)
  delays430<- rbind(delays430,sim)
}
REVSP430[, c(4:9)] %>% group_by(time) %>% summarize(m.fc=mean(fc.cherry))->REVSP430.m
delays430 %>% 
  filter()->del.430
del.430<-merge(del.430, REVSP430.m, all.x = T)
del.430 %>% mutate(res=m.fc-Y)->del.430
del.430 %>% group_by(t) %>% 
  summarize(N=sum(!is.na(res), na.rm=T),MAE=(sum(abs(res), na.rm=T))/(N-1))->sum.res.430
min(sum.res.430$MAE)
sum.res.430$t[sum.res.430$MAE==min(sum.res.430$MAE)] ->d430

#0h
delaysREV<-rbind(delaysREV,data.frame(plasmid="SP430", d_rev=d430))

p<-ggplot(data=sum.res.430, aes(x=t, y=MAE))+
  geom_point(size=.1, alpha=0.4, color="black") +
  geom_point(data=sum.res.430, aes(x=sum.res.430$t[sum.res.430$MAE==min(sum.res.430$MAE)], y=min(MAE)), size=.8, alpha=1, color="red") +
  coord_cartesian(x=c(0,25),y=c(0,0.3))+
  labs(y= "MAE" , x = "Delay (hours)")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=sum.res.430, aes(y=min(MAE)), lty=2)

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("MAE_REV_dCas9_mcherry.pdf", fix, path=out_path)

#the delays are
#for dcas9
delay_SP430 = sum.res.430$t[sum.res.430$MAE==min(sum.res.430$MAE)]
#for krab-dcas9
delay_SP428 =sum.res.428$t[sum.res.428$MAE==min(sum.res.428$MAE)]
#for krab-split-dcas9
delay_SP430A =sum.res.430ABA$t[sum.res.430ABA$MAE==min(sum.res.430ABA$MAE)]
#for hdac4-dcas9
delay_SP427 =sum.res.427$t[sum.res.427$MAE==min(sum.res.427$MAE)]


#run again the ODE models, including the delays

REV %>% filter(plasmid=="SP430ABA")->REVSP430ABA
mean(REVSP430ABA$norm.bfp[REVSP430ABA$time==0])->R
mean(REVSP430ABA$fc.cherry[REVSP430ABA$time==0])->Y
state<-c(R=1, Y=Y/sel_par$alpha)
ode1<- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-  - R* (log(2)/t1.2)
    dY <- (K^n)/(K^n+R^n)-alpha*Y
    return (list(c (dR, dY)))
  }
  )
} 

sel_par = all_par %>% filter(plasmid=='SP430A')
parameters = c(t1.2=sel_par$t_down, K=sel_par$K, n=sel_par$n, alpha=sel_par$alpha)

time <- seq(0+delay_SP430A, 150+delay_SP430A, by = 0.01)
out430ABAd <- lsodes(y=state,times=time, func=ode1, parms=parameters)
out430ABAd[,3]<-out430ABAd[,3]*sel_par$alpha


mean(REVSP428$norm.bfp[REVSP428$time==0])->R
mean(REVSP428$fc.cherry[REVSP428$time==0])->Y
state<-c(R=1, Y=Y/sel_par$alpha)
ode1<- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-  - R* (log(2)/t1.2)
    dY <- (K^n)/(K^n+R^n)-alpha*Y
    return (list(c (dR, dY)))
  }
  )
} 

sel_par = all_par %>% filter(plasmid=='SP428')
parameters = c(t1.2=sel_par$t_down, K=sel_par$K, n=sel_par$n, alpha=sel_par$alpha)


time <- seq(0+delay_SP428, 150+delay_SP428, by = 0.01)
out428d <- ode(y=state,times=time, func=ode1, parms=parameters)

out428d[,3]<-out428d[,3]*sel_par$alpha

mean(REVSP427$norm.bfp[REVSP427$time==0])->R
mean(REVSP427$fc.cherry[REVSP427$time==0])->Y
state<-c(R=1, Y=Y/sel_par$alpha)
ode1<- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-  - R* (log(2)/t1.2)
    dY <- (K^n)/(K^n+R^n)-alpha*Y
    return (list(c (dR, dY)))
  }
  )
} 

sel_par = all_par %>% filter(plasmid=='SP427')
parameters = c(t1.2=sel_par$t_down, K=sel_par$K, n=sel_par$n, alpha=sel_par$alpha)


time <- seq(0+delay_SP427, 150+delay_SP427, by = 0.01)
out427d <- ode(y=state,times=time, func=ode1, parms=parameters)
out427d[,3]<-out427d[,3]*sel_par$alpha

p<-ggplot(REVSP427,aes(time,fc.cherry))+
  geom_point(size=.8, alpha=0.4, color="#E64B35FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out427d), aes(time, Y),lty=2)+
  geom_line(data=data.frame(out427), aes(time, Y))
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_mCherry_HDAC4-dCas9_6h_delay.pdf", fix, path=out_path)


p<-ggplot(REVSP428,aes(time,fc.cherry))+
  geom_point(size=.8, alpha=0.4, color="#E64B35FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out428d), aes(time, Y),lty=2)+
  geom_line(data=data.frame(out428), aes(time, Y), lty=1)
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_mCherry_KRAB-dCas9_18.5h_delay.pdf", fix, path=out_path)

p<-ggplot(REVSP430ABA,aes(time,fc.cherry))+
  geom_point(size=.8, alpha=0.4, color="#E64B35FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out430ABAd), aes(time, Y),lty=2)+
  geom_line(data=data.frame(out430ABA), aes(time, Y))

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_mCherry_KRAB-Split-dCas9_16.5h_delay.pdf", fix, path=out_path)

p<-ggplot(REVSP430,aes(time,fc.cherry))+
  geom_point(size=.8, alpha=0.4, color="#E64B35FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out430), aes(time, Y))

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_mCherry_dCas9_delay.pdf", fix, path=out_path)

p<-ggplot(REVSP411,aes(time,fc.cherry))+
  geom_point(size=.8, alpha=0.4, color="#E64B35FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "mCherry")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out411), aes(time, Y))

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_mCherry_CasRx_delay.pdf", fix, path=out_path)

p<-ggplot(REVSP427,aes(time,norm.bfp))+
  geom_point(size=.8, alpha=0.4, color="#4DBBD5FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out427), aes(time, R))
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_tagBFP_HDAC4-dCas9_6h_delay.pdf", fix, path=out_path)


p<-ggplot(REVSP428,aes(time,norm.bfp))+
  geom_point(size=.8, alpha=0.4, color="#4DBBD5FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out428), aes(time, R), lty=1)
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_tagBFP_KRAB-dCas9_18.5h_delay.pdf", fix, path=out_path)



p<-ggplot(REVSP430ABA,aes(time,norm.bfp))+
  geom_point(size=.8, alpha=0.4, color="#4DBBD5FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out430ABA), aes(time, R))

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_tagBFP_KRAB-Split-dCas9_16.5h_delay.pdf", fix, path=out_path)

p<-ggplot(REVSP430,aes(time,norm.bfp))+
  geom_point(size=.8, alpha=0.4, color="#4DBBD5FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out430), aes(time, R))

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_tagBFP_dCas9_delay.pdf", fix, path=out_path)

p<-ggplot(REVSP411,aes(time, norm.bfp))+
  geom_point(size=.8, alpha=0.4, color="#4DBBD5FF") +
  coord_cartesian(x=c(0,150),y=c(0,1.3))+
  labs(x= "Time (hours)" , y = "")+
  scale_color_npg()+
  scale_fill_npg()+
  geom_line(data=data.frame(out411), aes(time, R))

fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("REV_tagBFP_CasRx_delay.pdf", fix, path=out_path)


write_csv(delaysREV,file='parameters/delays_derepression.csv')

