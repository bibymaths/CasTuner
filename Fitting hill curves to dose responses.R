#I am going to fit "Hill functions" to the dTAG-13 dose-response curves of the reporter-repressor:
#load libraries
library(flowCore)
library(openCyto)
library(ggcyto)
library(extrafont)
library(tidyverse)
library(minpack.lm)

#set theme
theme_set(theme_classic() +
            theme(legend.text = element_text(size = 6), panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                  axis.line = element_blank(), axis.text = element_text(size = 6, color= "black"),
                  axis.text.x = element_text(vjust = 0.5, color= "black"),
                  axis.text.y = element_text(vjust = 0.5, color= "black"),
                  axis.title = element_text(size = 6), strip.text = element_text(size = 6, color= "black"),
                  strip.background = element_blank(), legend.title = element_blank()))

out_path='output'

# calculate background fluorescence of non fluorescent control for background subtraction
fMyFlowSet <- read.flowSet(path="fcs_files/NFC", min.limit=0.01)
chnl <- c("FSC-A", "SSC-A") #are the channels to build the gate on
bou_g <- openCyto:::.boundary(fMyFlowSet[[3]], channels = chnl, min = c(0.4e5, 0.20e5), max=c(2e5,1.3e5), filterId = "Boundary Gate")
p <- autoplot(fMyFlowSet[[3]], x = chnl[1], y = chnl[2], bins=100)
p +geom_gate(bou_g)
BoundaryFrame<- Subset (fMyFlowSet[[3]], bou_g) #is the subset of cells within the gate
chnl <- c("FSC-A", "FSC-H")
sing_g <- openCyto:::.singletGate(BoundaryFrame, channels = chnl)
p <- autoplot(BoundaryFrame, x = chnl[1], y = chnl[2])
p + geom_gate(sing_g)
fSinglets_MyFlowSet<-Subset(fMyFlowSet, bou_g%&% sing_g)
fmedianexp<- as.data.frame(fsApply(fSinglets_MyFlowSet, each_col, median))
mBFP_neg<-mean(fmedianexp[1:3,7])
mmCherry_neg<-mean(fmedianexp[1:3,8])

# load dose-response data
MyFlowSet1 <- read.flowSet(path="fcs_files/dose_response_data", min.limit=0.01)
chnl <- c("FSC-A", "SSC-A") #are the channels to build the gate on
bou_g <- openCyto:::.boundary(MyFlowSet1[[2]], channels = chnl, min = c(0.4e5, 0.20e5), max=c(2e5,1.3e5), filterId = "Boundary Gate")
p <- autoplot(MyFlowSet1[[2]], x = chnl[1], y = chnl[2], bins=100)
p +geom_gate(bou_g)
BoundaryFrame<- Subset (MyFlowSet1[[3]], bou_g) #is the subset of cells within the gate
chnl <- c("FSC-A", "FSC-H")
sing_g <- openCyto:::.singletGate(BoundaryFrame, channels = chnl)
p <- autoplot(BoundaryFrame, x = chnl[1], y = chnl[2])
p + geom_gate(sing_g)
Singlets_MyFlowSet1<-Subset(MyFlowSet1, bou_g%&% sing_g)
medianexp<- as.data.frame(fsApply(Singlets_MyFlowSet1, each_col, median))

#subtract background fluorescence and calculate median values
Transf<-linearTransform(transformationId="defaultLinearTransform", a = 1, b = -mBFP_neg)
Norm_Singlets_MyFlowSet1<- transform(Singlets_MyFlowSet1, transformList('BV421-A' ,Transf))
Transf<-linearTransform(transformationId="defaultLinearTransform", a = 1, b = -mmCherry_neg)
Norm_Singlets_MyFlowSet<- transform(Norm_Singlets_MyFlowSet1, transformList('PE-A' ,Transf))
medianexp<- as.data.frame(fsApply(Norm_Singlets_MyFlowSet1, each_col, median))
medianexp<-medianexp[, c(7:8)]

#give right dTAG concentrations to columns
medianexp %>% rownames_to_column()->medianexp
medianexp %>% separate(1, c( "plasmid", "guide", "dTAG", NA), sep = "_") ->medianexp
medianexp %>% separate(3, into = c(NA, "dTAG"), 
                       sep = "(?<=[A-Za-z])(?=[0-9])")->medianexp
medianexp$dTAG[medianexp$dTAG=="1"]<-0
medianexp$dTAG[medianexp$dTAG=="2"]<-0.5
medianexp$dTAG[medianexp$dTAG=="3"]<-1
medianexp$dTAG[medianexp$dTAG=="4"]<-2
medianexp$dTAG[medianexp$dTAG=="5"]<-3
medianexp$dTAG[medianexp$dTAG=="6"]<-5
medianexp$dTAG[medianexp$dTAG=="10"]<-50
medianexp$dTAG[medianexp$dTAG=="8"]<-10
medianexp$dTAG[medianexp$dTAG=="7"]<-8
medianexp$dTAG[medianexp$dTAG=="9"]<-25
medianexp$dTAG[medianexp$dTAG=="11"]<-100
medianexp$dTAG[medianexp$dTAG=="12"]<-500
medianexp$dTAG<-as.numeric(medianexp$dTAG)
medianexp->d4
d4$plasmid<-factor(d4$plasmid, levels=c("430", "428", "430ABA", "427", "411"))

#normalize mCherry to mean mCherry in Non-targeting control
d4 %>% dplyr::filter(guide=="N") %>% group_by(plasmid) %>% mutate(meanNTC= mean(`PE-A`))->meanNTC
d4<-left_join(d4, meanNTC, by=c("plasmid", "dTAG"))
d4 %>%  group_by(plasmid) %>% mutate(fc= `PE-A.x`/ meanNTC)->d4
d4 %>% group_by(plasmid, guide.x) %>% dplyr::filter(dTAG==0) %>%  summarize(max.bfp=mean(`BV421-A.x`)) ->max.bfp
d4<-left_join(d4, max.bfp, by=c("plasmid", "guide.x"))

#normalize bfp to maximal bfp
d4 %>% group_by(plasmid, guide.x) %>% mutate(fc.bfp=`BV421-A.x`/max.bfp)->d4
d4 %>% group_by(plasmid, guide.x) %>% mutate(norm.bfp=(`BV421-A.x`-min(`BV421-A.x`))/(max(`BV421-A.x`)-min(`BV421-A.x`)))->d4

#we need only the cells with targeting guides (guide.x=="G") so we filter for them
filter<-dplyr::filter
d4%>% filter(guide.x=="G")->d4g

#to find the parameters use nlsLM from minpack.lm
#for dCas9, the fit is:
d4g %>%filter(plasmid=="430")->SP430
R=SP430$norm.bfp
y=SP430$fc
data<-data.frame(y,R)
fit = nlsLM(y ~ K^n/(K^n+R^n),
            start = list(K = 0.1, n = 1),
            data = data)
summarysp430<-summary(fit)
parameters=data.frame(plasmid = 'SP430', K =coef(fit)[1], n = coef(fit)[2])

p <-ggplot(SP430,aes(norm.bfp, fc))+
  stat_function(fun=function(norm.bfp){1/(1+(norm.bfp/coef(fit)[1])^coef(fit)[2])}, color="black")+
  geom_point(size=.8, alpha=0.4, color="grey") +# adding connecting lines
  coord_cartesian(y=c(0,1))+
  labs(x= "Normalized repressor level" , y = "Normalized reporter level")
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("Hill_dCas9.pdf", fix, path=out_path)

#for CasRx the fit is
d4g %>% dplyr::filter(plasmid=="411")->SP411
R=SP411$norm.bfp
y=SP411$fc
data<-data.frame(y,R)
fit = nlsLM(y ~ K^n/(K^n+R^n),
            start = list(K = 0.1, n = 1),
            data = data)
summarysp411<-summary(fit)
parameters =rbind(parameters,c('SP411',coef(fit)[1],coef(fit)[2]))

p<-ggplot(SP411,aes(norm.bfp, fc))+
  stat_function(fun=function(norm.bfp){1/(1+(norm.bfp/coef(fit)[1])^coef(fit)[2])}, color="black")+
  geom_point(size=.8, alpha=0.4, color="grey") +# adding connecting lines
  coord_cartesian(y=c(0,1))+
  labs(x= "Normalized repressor level" , y = "Normalized reporter level")
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("Hill_CasRx.pdf", fix, path=out_path)

#for hdac4-dCas9, the fit is
d4g %>%filter(plasmid=="427")->SP427
R=SP427$norm.bfp
y=SP427$fc
data<-data.frame(y,R)
fit = nlsLM( y ~ K^n/(K^n+R^n),
             start = list(K = 0.1, n = 1),
             data = data)
summarysp427<-summary(fit)
parameters =rbind(parameters,c('SP427',coef(fit)[1],coef(fit)[2]))

p<-ggplot(SP427,aes(norm.bfp, fc))+
  stat_function(fun=function(norm.bfp){1/(1+(norm.bfp/coef(fit)[1])^coef(fit)[2])}, color="black")+
  geom_point(size=.8, alpha=0.4, color="grey") +# adding connecting lines
  coord_cartesian(y=c(0,1))+
  labs(x= "Normalized repressor level" , y = "Normalized reporter level")
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("Hill-HDAC4-dCas9.pdf", fix, path=out_path)

#for KRAB-dCas9 the fit is
d4g %>%filter(plasmid=="428")->SP428
R=SP428$norm.bfp
y=SP428$fc
data<-data.frame(y,R)
fit = nlsLM( y ~ K^n/(K^n+R^n),
             start = list(K = 0.1, n = 1),
             data = data)
summarysp428<-summary(fit)
parameters =rbind(parameters,c('SP428',coef(fit)[1],coef(fit)[2]))

p<-ggplot(SP428,aes(norm.bfp, fc))+
  stat_function(fun=function(norm.bfp){1/(1+(norm.bfp/coef(fit)[1])^coef(fit)[2])}, color="black")+
  geom_point(size=.8, alpha=0.4, color="grey") +# adding connecting lines
  coord_cartesian(y=c(0,1))+
  labs(x= "Normalized repressor level" , y = "Normalized reporter level")
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("Hill-KRAB-dCas9.pdf", fix, path=out_path)

#For KRAB-Split-dCas9 the fit is
d4g %>% dplyr::filter(plasmid=="430ABA")->SP430A
R=SP430A$norm.bfp
y=SP430A$fc
data<-data.frame(y,R)
fit = nlsLM( y ~ K^n/(K^n+R^n),
             start = list(K = 0.1, n = 1),
             data = data)
summarysp430A<-summary(fit)
parameters =rbind(parameters,c('SP430',coef(fit)[1],coef(fit)[2]))
p<-ggplot(SP430A,aes(norm.bfp, fc))+
  stat_function(fun=function(norm.bfp){1/(1+(norm.bfp/coef(fit)[1])^coef(fit)[2])}, color="black")+
  geom_point(size=.8, alpha=0.4, color="grey") +# adding connecting lines
  coord_cartesian(y=c(0,1))+
  labs(x= "Normalized repressor level" , y = "Normalized reporter level")
fix <- set_panel_size(p, width = unit(1.5*1.618, "cm"), height = unit(1.5, "cm"))
ggsave("Hill-KRAB-Split-dCas9.pdf", fix, path=out_path)

#summary of fittings
summarysp411
summarysp427
summarysp428
summarysp430A
summarysp430

write_csv(parameters,file='output/Hill_parameters.csv')
