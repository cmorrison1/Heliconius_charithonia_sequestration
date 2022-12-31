########################################################################
########################################################################
########################################################################
#### --------------------------------------------------------------- ###
### ---------------- C. Morrison, C. Nguyen & L. Gilbert ----------- ###
### --- Heliconius charithonia cyanogenic glycoside concentation --- ###
### ------------------------- ANALYSES ----------------------------- ###
### ---------------------------------------------------------------- ###
########################################################################
########################################################################
########################################################################


# Colin Richard Morrison 
# PhD Candidate
# The University of Texas at Austin 
# Department of Integrative Biology 
# Graduate Program In Ecology, Evolution and Behavior
# crmorrison@utexas.edu


getwd()
setwd('~/Documents/CRM_manuscripts_essays/Morrison_etal_charithonia sequestration_2022/data_analysis')


library(multcomp)
library(plyr)
library(ggplot2)
library(ggpubr)
library(ggplot2)
library(nlme)
library(cowplot)
library(plotrix)
#library(glmm)
library(MuMIn)


### --- Make a summary data table with mean, sd, and se plots
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}



### --- Count of adult butterfly sample sizes by host for sample sizes
adult<-read.csv("char.adult.csv") # 
adultN<-adult[!is.na(adult$ug.vol.correct),]

aff=subset(adultN,pair=="P.affinis") # 23
bif=subset(adultN,pair=="P.biflora") # 11
lut=subset(adultN,pair=="P.lutea") # 14
sub=subset(adultN,pair=="P.suberosa") # 12
nrow(sub)


### -------------------------------------------------------------------- ### 
### --- Did larval mass predict cyanide toxicity on different hosts? --- ###
### -------------------------------------------------------------------- ### 
# linear model to check for normality
Lm1=lm(ug.vol.correct~mass.mg.,data=cats)
plot(hist(resid(Lm1))) # fairly normal distribution 

# check for homoscedasticity 
boxplot(ug.vol.correct~host, data = cats) # no heterogeneity of variance 

cats<-cats[!is.na(cats$ug.vol.correct),]
gof<-glm(ug.vol.correct~mass.mg.,data=cats)
anova(gof,test="F")
#         Df Deviance Resid. Df Resid. Dev     F    Pr(>F)    
# NULL                        73    1557420                    
# mass.mg.  1   661276        72     896143 53.13 3.229e-10 ***
r.squaredGLMM(gof) 
#            R2m       R2c
# [1,] 0.4212311 0.4212311

# --- P. lutea  
lut=subset(cats,host=="lut")
l.sum<-glm(ug.vol.correct~mass.mg.,data=lut)
summary(l.sum) 
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.09392   15.46035  -0.006    0.995    
# mass.mg.     0.62109    0.09405   6.604 8.38e-06 ***
r.squaredGLMM(l.sum) 
#            R2m       R2c
# [1,] 0.7315877 0.7315877

# --- P. biflora 
bif=subset(cats,host=="bi")
b.sum<-glm(ug.vol.correct~mass.mg.,data=bif)
summary(b.sum) 
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.92054   22.12901   0.177 0.860750    
# mass.mg.     0.40902    0.09302   4.397 0.000165 ***
r.squaredGLMM(b.sum)
#            R2m       R2c
# [1,] 0.4172841 0.4172841

# --- P. affinis 
aff=subset(cats,host=="aff")
a.sum<-glm(ug.vol.correct~mass.mg.,data=aff)
summary(a.sum) 
#            Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -18.4977    22.9308  -0.807    0.427    
# mass.mg.      1.4806     0.1207  12.271 1.49e-12 ***
r.squaredGLMM(a.sum)
#            R2m       R2c
# [1,] 0.8432086 0.8432086


### ------------------------------------------------------- ### 
### --- Was mass significantly predicted by life stage? --- ###
### ------------------------------------------------------- ### 
# correlation between size and development stage
sizeCor<-glm(mass.mg.~life.stage,data=cats)
summary(sizeCor) 
anova(sizeCor,test="F")
#            Df Deviance Resid. Df Resid. Dev      F    Pr(>F)    
# NULL                          73    1008963                     
# life.stage  1   607056        72     401906 108.75 4.867e-16 ***
r.squaredGLMM(sizeCor) 
#            R2m       R2c
# [1,] 0.5983534 0.5983534


### ------------------------------------------------------------------ ### 
### --- What was the critical life stage where toxicity increased? --- ###
### ------------------------------------------------------------------ ###
# Across instars, host plant controlled for 
cats$life.stage<-as.factor(cats$life.stage)
inst.sum<-glm(ug.vol.correct~life.stage,data=cats)
summary(inst.sum) 
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    31.11      30.28   1.027    0.307    
# life.stage4    48.18      41.99   1.147    0.255    
# life.stage5   257.22      39.31   6.543 4.87e-09 ***
summary(glht(inst.sum, linfct=mcp(life.stage="Tukey")))
#                Estimate Std. Error z value Pr(>|z|) 
# 4 - 3 == 0    45.62      38.21   1.194 0.455997    
# 5 - 3 == 0   181.09      36.25   4.996  < 1e-04 ***
# 5 - 4 == 0   135.48      33.83   4.005 0.000189 ***


### --------------------------------------------------------- ### 
### --- Did larvae from different hosts acquire more CN?  --- ###
### --------------------------------------------------------- ###
# make master database with all hosts
# biflora=read.csv("bi.csv")
# affinis=read.csv("aff.csv")
# lutea=read.csv("lut2.csv")
cats<-rbind.data.frame(lutea[,1:15],affinis[,1:15],biflora[,1:15])
cats$host<-as.factor(cats$host)
cats$log.ug.correct<-log(cats$ug.vol.correct)
cats<-cats[!is.na(cats$ug.vol.correct),]

# compare across each instar
CATstats <- data_summary(cats, varname="ug.vol.correct", 
                         groupnames=c("host","life.stage"))
CATstats
#         host life.stage ug.vol.correct        sd        se
# 1   P. lutea          3       12.41220  22.55222  7.973415
# 2   P. lutea          4       68.50507  26.36306  8.787687
# 3   P. lutea          5      363.30867 280.27727 84.506777
# 4 P. affinis          3       36.08371  25.15559  8.385198
# 5 P. affinis          4       78.51092  36.39747 12.868449
# 6 P. affinis          5      391.55721 175.02527 50.525444
# 7 P. biflora          3       46.07096  37.59345 14.208987
# 8 P. biflora          4       90.74812  82.39701 27.465669
# 9 P. biflora          5      116.35226  63.42440 18.309046

### - 3rd instar 
cats3=subset(cats,life.stage=="3")
i3<-glm(ug.vol.correct~host,data=cats3)
summary(glht(i3, linfct=mcp(host="Tukey")))
#                Estimate Std. Error z value Pr(>|z|)  
# bi - aff == 0     9.987     14.382   0.694   0.7666  
# lut - aff == 0  -23.672     13.867  -1.707   0.2024  
# lut - bi == 0   -33.659     14.770  -2.279   0.0587 .

### - 4th instar 
cats4=subset(cats,life.stage=="4")
i4<-glm(ug.vol.correct~host,data=cats4)
summary(glht(i4, linfct=mcp(host="Tukey")))
#                Estimate Std. Error z value Pr(>|z|)  
# bi - aff == 0     12.24      26.64   0.459    0.890
# lut - aff == 0   -10.01      26.64  -0.376    0.925
# lut - bi == 0    -22.24      25.85  -0.861    0.665

### - 5th instar
cats5=subset(cats,life.stage=="5")
i5<-glm(ug.vol.correct~host,data=cats5)
plot(ug.vol.correct~host,data=cats5)
summary(glht(i5, linfct=mcp(host="Tukey")))
#                Estimate Std. Error z value Pr(>|z|)  
# bi - aff == 0   -275.20      77.95  -3.530  0.00116 **
# lut - aff == 0   -28.25      79.71  -0.354  0.93310   
# lut - bi == 0    246.96      79.71   3.098  0.00543 **


### ------------------------------------------------------------------------ ###
### --- CN concentration different among hosts at different lifestages ? --- ###
### ------------------------------------------------------------------------ ###
### raw adult butterfly CN concentration data
adult<-read.csv("char.adult.csv")
hc2<-adult[,c('number','species','pair','ug','ug.vol.correct','ug.mg')] # select columns I need
hc2[hc2 == "H.charithonia"] <- "adult"
names(hc2)<-c("number","life.stage","host","ug","ug.vol.correct","ug.mg")
head(hc2)

### raw larval CN concentration data
cats<-rbind.data.frame(lutea[,1:15],affinis[,1:15],biflora[,1:15])
cats2 <-cats[,c('number','host','life.stage','ug','ug.vol.correct','ug.mg')] # select columns I need
cats2 <- cats2[, c(2, 1, 3, 4, 5)]
cats2[cats2 == "aff"] <- "P.affinis"
cats2[cats2 == "bi"] <- "P.biflora"
cats2[cats2 == "lut"] <- "P.lutea"
head(cats2)

### combine these DFs
comparison3<-rbind.data.frame(hc2,cats2)
comparison3$host<-as.factor(comparison3$host)
tail(comparison3)

### - 3rd instar 
cats3=subset(comparison3,life.stage=="3")
cats3<-cats3[!is.na(cats3$ug.mg)]
i3<-glm(ug.mg~host,data=cats3)
#plot(hist(residuals(i3)),breaks=4) #normal-shaped error distribution 
#shapiro.test(residuals(i3)) # W = 0.88821, p-value = 0.02993
summary(glht(i3, linfct=mcp(host="Tukey")))
#                Estimate Std. Error z value Pr(>|z|)  
# P. affinis - P. lutea == 0     0.6270     0.3792   1.654    0.223
# P. biflora - P. lutea == 0    -0.1517     0.4038  -0.376    0.925
# P. biflora - P. affinis == 0  -0.7787     0.3932  -1.980    0.117

### - 4th instar 
cats4=subset(comparison3,life.stage=="4")
i4<-glm(ug.mg~host,data=cats4)
#plot(hist(residuals(i4)),breaks=6) #fairly normal-shaped error distribution 
#shapiro.test(residuals(i4)) # W = 0.91019, p-value = 0.0356
summary(glht(i4, linfct=mcp(host="Tukey")))
#                Estimate Std. Error z value Pr(>|z|)  
# P. affinis - P. lutea == 0     0.4394     0.1871   2.348   0.0494 *  
# P. biflora - P. lutea == 0    -0.7828     0.1815  -4.313   <0.001 ***
# P. biflora - P. affinis == 0  -1.2222     0.1871  -6.533   <0.001 ***

### - 5th instar
cats5=subset(comparison3,life.stage=="5")
i5<-glm(ug.mg~host,data=cats5)
#plot(hist(residuals(i5)),breaks=6) #fairly normal-shaped error distribution 
#shapiro.test(residuals(i5)) # W = 0.8363, p-value = 0.0002622
summary(glht(i5, linfct=mcp(host="Tukey")))
#                Estimate Std. Error z value Pr(>|z|)  
# P. affinis - P. lutea == 0    0.01672    0.25668   0.065    0.998    
# P. biflora - P. lutea == 0   -1.27805    0.25668  -4.979 1.56e-06 ***
# P. biflora - P. affinis == 0 -1.29477    0.25104  -5.158  < 1e-06 ***

### - adults
adults=subset(comparison3,life.stage=="adult")
hc2$host<-as.factor(hc2$host)
As<-glm(ug.mg~host,data=hc2)
summary(glht(As, linfct=mcp(host="Tukey")))
#                Estimate Std. Error z value Pr(>|z|) 
# P.biflora - P.affinis == 0  -3.40352    1.19035  -2.859   0.0219 *
# P.lutea - P.affinis == 0     0.06124    1.10070   0.056   0.9999  
# P.suberosa - P.affinis == 0  0.03495    1.15631   0.030   1.0000  
# P.lutea - P.biflora == 0     3.46476    1.30829   2.648   0.0400 *
# P.suberosa - P.biflora == 0  3.43847    1.35541   2.537   0.0540 .
# P.suberosa - P.lutea == 0   -0.02630    1.27740  -0.021   1.0000  


### -------------------------------------------- ###
### --- Did larval mass differ among hosts?  --- ###
### -------------------------------------------- ###
### - Summary data table 
MASSstats <- data_summary(cats, varname="mass.mg.", 
                          groupnames=c("host","life.stage"))
MASSstats
#   host life.stage  mass.mg.        sd        se
# 1   aff          3  36.15420  6.981656  2.207794
# 2   aff          4  56.20412 14.864907  5.255538
# 3   aff          5 282.98783 65.254400 18.837323
# 4    bi          3  76.59891 19.248399  5.803611
# 5    bi          4 197.68792 99.009529 27.460303
# 6    bi          5 296.82558 33.189458  9.580971
# 7   lut          3  41.62075  8.457248  2.990089
# 8   lut          4  77.83256 17.867311  5.955770
# 9   lut          5 247.10564 78.887836 23.785577

cats5=subset(cats,life.stage=="5")
m5<-glm(mass.mg.~host,data=cats5)
# plot(mass.mg.~host,data=cats5)
summary(glht(m5, linfct=mcp(host="Tukey")))
#                Estimate Std. Error z value Pr(>|z|)  
# Passiflora affinis - Passiflora lutea == 0      35.88      28.10   1.277   0.4082  
# Passiflora biflora - Passiflora lutea == 0      64.77      28.10   2.305   0.0552 .
# Passiflora biflora - Passiflora affinis == 0    28.89      27.48   1.051   0.5446  


### ---------------------------------------------------------- ### 
### --- Did adults from some hosts acquire more total CN?  --- ###
### ---------------------------------------------------------- ### 
adult<-read.csv("char.adult.csv")
head(adult)
hc<-adult[,c('number','species','pair','sex','ug','ug.vol.correct','ug.mg',"extraction.date")]
hc$pair<-as.factor(hc$pair)

# remove NAs
hc<-hc[!is.na(hc$ug.vol.correct),]

# check for normality
Lm1=lm(ug.vol.correct~pair,data=hc)
plot(hist(resid(Lm1),breaks = 6)) #  normal distribution 
shapiro.test(resid(Lm1)) # W = 0.95937, p-value = 0.04381

# check for homoscedasticity 
boxplot(ug.vol.correct~pair,data=hc) # Some heterogeneity of variance 

palatable<-lme(ug.vol.correct~pair,weights=varIdent(form=~1|as.factor(pair)),
               random = ~1|sex, method = "REML",data=hc)
summary(palatable)
#                    Value Std.Error DF   t-value p-value
# (Intercept)     419.4452  56.07642 55  7.479885  0.0000
# pairP.biflora  -287.8077  61.03991 55 -4.715075  0.0000
# pairP.lutea      29.2179 114.11951 55  0.256029  0.7989
# pairP.suberosa  -39.3966  72.45681 55 -0.543726  0.5888
r.squaredGLMM(palatable) 
#            R2m       R2c
# 0.6705583 0.6705583
anova(palatable)
# (Intercept)     1    55 131.23927  <.0001
# pair            3    55  14.49937  <.0001

# adult CN descriptive statistics 
ADULTstats <- data_summary(hc, varname="ug.vol.correct", 
                           groupnames=c("pair","species"))
ADULTstats
#        pair       species ug.vol.correct        sd       se
# 1  P.affinis H.charithonia       419.4452 268.93315 56.07644
# 2  P.biflora H.charithonia       131.6374  79.96475 24.11028
# 3    P.lutea H.charithonia       448.6631 371.88945 99.39164
# 4 P.suberosa H.charithonia       380.0485 158.94995 45.88490


######################################
### - LMM of adult concentration - ###
######################################
adult<-read.csv("char.adult.csv")
hc<-adult[,c('number','species','sex','extraction.date','pair','ug','ug.vol.correct','ug.mg')]
adults<-hc[!is.na(hc$ug.mg),]

CONpalatable<-lme(ug.mg~pair,weights=varIdent(form=~1|as.factor(pair)),
                  random = ~1|sex, method = "REML",data=adults)

# check distribution
shapiro.test(resid(CONpalatable)) # W = 0.94884, p-value = 0.01374 --> close to normally distributed
plot(hist(resid(CONpalatable),breaks = 6))

summary(CONpalatable)
#                    Value Std.Error DF   t-value p-value
# (Intercept)     4.787976 0.7872465 55  6.081928  0.0000
# pairP.biflora  -3.346723 0.7964875 55 -4.201852  0.0001
# pairP.lutea     0.124158 1.3446682 55  0.092334  0.9268
# pairP.suberosa  0.002514 1.0232504 55  0.002457  0.9980
r.squaredGLMM(CONpalatable) 
#            R2m       R2c
# 0.7600878 0.7872198
anova(CONpalatable)
#             numDF denDF   F-value p-value
# (Intercept)     1    55 60.32792  <.0001
# pair            3    55 14.47492  <.0001


################################################################ 
### - Correlation between adult and host CN concentration? - ###
################################################################ 
# Host plant data
pass<-read.csv("char.host.csv")
head(pass)
pass2<-pass[,c('species','pair','ug.mg')]
pass2

# Adult data
adult<-read.csv("char.adult.csv")
head(adult)
adult2<-adult[,c('species','pair','ug.mg')]
adult2

# Combine host and adult databases
cn<-rbind.data.frame(pass2,adult2)

### --- Host-Adult concentration descriptive statistics 
concens <- data_summary(cn, varname="ug.mg", 
                        groupnames=c("species","pair"))
concens
#         species       pair      ug.mg         sd          se
# 1 H.charithonia  P.affinis 4.83406502 3.69438698 0.770332936
# 2 H.charithonia  P.biflora 1.43054584 0.70992859 0.214051523
# 3 H.charithonia    P.lutea 4.89530844 4.15248160 1.109797389
# 4 H.charithonia P.suberosa 4.86901195 2.35442912 0.679665145
# 5     P.affinis  P.affinis 0.03268708 0.02747129 0.007930279
# 6     P.biflora  P.biflora 0.22399147 0.23990091 0.042408890
# 7       P.lutea    P.lutea 0.27816993 0.51978591 0.134208146
# 8    P.suberosa P.suberosa 0.02638875 0.02826091 0.008158221

### --- how much higher was adults  that their hosts?
#         species     ug.mg 
# 1 H.charithonia   4.83406502 
#   P.affinis       0.3268708
4.83406502/0.3268708 # = 14.78892 X
# 2 H.charithonia   1.43054584 
#   P.biflora       0.22399147
1.43054584 /  0.22399147 # = 6.386609 X
# 3 H.charithonia   4.89530844 
#   P.lutea         0.27816993
4.89530844/0.27816993 # = 17.59827 X

### --- Correlation across species?
#cor.test(concens[1:4,]$ug.mg, concens[5:8,]$ug.mg)
# t = 0.50439, df = 2, p-value = 0.6641
# cor = 0.3359322 


########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################