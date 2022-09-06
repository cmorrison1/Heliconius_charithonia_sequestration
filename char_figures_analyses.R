########################################################################
########################################################################
########################################################################
#### --------------------------------------------------------------- ###
### ---------------- C. Nguyen, C. Morrison & L. Gilbert ----------- ###
### --- Heliconius charithonia cyanogenic glycoside concentation --- ###
### --------------------- Figures and analyses --------------------- ###
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
setwd('~/Desktop/charithonia_sequestration_Nguyen_Morrison')


library(multcomp)
library(plyr)
library(ggplot2)
library(ggpubr)
library(ggplot2)
library(nlme)
#library(glmm)
#library(MuMIn)


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

##############################################################
##############################################################
############################ PLOTS ########################### 
##############################################################
##############################################################
--------------------------------------------------------------
--------------------------------------------------------------
################################################################ 
### ---///  FIGURE 2: Larval and adult sequestration  ///--- ###
################################################################ 


### ------------------------------------------- ###
### --- Panel A: larvae reared on each host --- ###
### ------------------------------------------- ###
# make master database with all hosts
biflora=read.csv("bi.csv")
affinis=read.csv("aff.csv")
lutea=read.csv("lut.csv")
cats<-rbind.data.frame(lutea[,1:15],affinis[,1:15],biflora[,1:15])
head(cats)

# --- three-paneled scatter plot of mass vs CGs 
# set graphics parmeters
cats$host = factor(cats$host,
                   levels= c("lut","aff","bi"),
                   labels = c("Passiflora lutea","Passiflora affinis",
                              "Passiflora biflora"))
# plot it 
develop<-ggscatter(cats, x = "mass.mg.", y = "ug.vol.correct",
                   add = "reg.line",                     
                   conf.int = TRUE,                     
                   add.params = list(color = "blue",
                                     fill = "lightgray")) +
  facet_grid(. ~ host) +
  panel_border() +
  xlab('Larval mass (mg)') +
  ylab(expression(paste("Cyanogenic glycosides (", mu, "g)"))) +
  theme(panel.grid.major.x = element_line(colour = "grey90"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(face = "italic",size=14),
        strip.background = element_rect(
          color="black", fill="white", size=1.0, linetype="solid"),panel.spacing = unit(0.5, "lines"),
        axis.text.x = element_text(size=14, colour = "grey30"),
        axis.text.y = element_text(size=14, colour = "grey30"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16,vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16,vjust=3),
        plot.margin = margin(t = 20, r = 75, b = -10, l = 0, unit = "pt")
  )

# print this in the dev.new window for a PDF
dev.new(
  title = "CG-mass",
  width = 10,
  height = 6,
  noRStudioGD = TRUE
)
develop

### ---------------------------------------- ###
### --- Panel B: compare larvae & adults --- ###
### ---------------------------------------- ###
# amount of adult CG per host association
adult<-read.csv("char.adult.csv")
hc<-adult[,c('species','pair','ug','ug.vol.correct','ug.mg')]
# summary data table
HCstats <- data_summary(hc, varname="ug.vol.correct", 
                        groupnames=c("pair","species"))
names(HCstats)<-c("host","life.stage","ug.total","sd","se")
HCstats[HCstats == "H.charithonia"] <- "adult"
HCstats

# amount of larval CG per host association
biflora=read.csv("bi.csv")
affinis=read.csv("aff.csv")
lutea=read.csv("lut.csv")
suberosa=read.csv("sub.csv")
cats2<-rbind.data.frame(lutea[,1:15],affinis[,1:15],biflora[,1:15],suberosa[,1:15])
# summary data table
CATstats <- data_summary(cats2, varname="ug.vol.correct", 
                         groupnames=c("host","life.stage"))
names(CATstats)<-c("host","life.stage","ug.total","sd","se")
CATstats[CATstats == "aff"] <- "P.affinis"
CATstats[CATstats == "bi"] <- "P.biflora"
CATstats[CATstats == "lut"] <- "P.lutea"
CATstats[CATstats == "sub"] <- "P.suberosa"
CATstats

# combine these DFs
comparison<-rbind.data.frame(CATstats,HCstats)
comparison

# set graphics color and order parameters 
comparison$host = factor(comparison$host,
                         levels= c("P.lutea","P.affinis",
                                   "P.suberosa","P.biflora"),
                         labels = c("P. lutea","P. affinis",
                                    "P. suberosa","P. biflora"))

surv.colors<-c("white","grey90","grey60","grey30")

# plot it 
allCGs<-ggplot(comparison, aes(x=host, y=ug.total,fill=life.stage)) + 
  geom_bar(width=0.75,
           stat="identity",
           colour = "black",
           position=position_dodge(),
           lwd=0.75) +
  geom_errorbar(aes(ymin=ug.total-se, ymax=ug.total+se),
                width=0.4,
                position=position_dodge(0.75),
                col='black') +
  ggtitle("") +
  xlab("Host plant") + 
  ylab(expression(paste("Cyanogenic glycosides (", mu, "g)"))) +
  labs(fill = "Life stage") +
  geom_hline(yintercept = 0,lwd=0.5,col="grey10") +
  scale_fill_manual(values = surv.colors) + 
  theme(panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_line(colour = "grey90"),
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(vjust=0,hjust=0.5,size = 16, face = "bold", colour = "black"),
        legend.title = element_text(size = 16, face = "bold", colour = "black"),
        legend.text = element_text(size = 14),
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(hjust=1,size=14, colour = "grey30"),
        axis.text.x = element_text(angle = 30,hjust=1,size=14, 
                                   face='italic', colour = "grey30"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 0),size=16,vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20),size=16,vjust=3),
        plot.margin = margin(t = 0, r = 5, b = 20, l = 30, unit = "pt")
  )
# print this in the dev.new window for a PDF
dev.new(
  title = "CG-stages",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
allCGs


# --- Combine into a 2-panel figure 
Fig2=ggarrange(develop,allCGs,
               labels = c("(a)", "(b)"),
               ncol = 1, nrow = 2,
               font.label = list(size = 25))
# dev plot
dev.new(
  title = "Figure2",
  width = 8.5,
  height = 11,
  noRStudioGD = TRUE
)
Fig2

### ------------------------------- ### 
### --- EXTRA exploratory plots --- ###
### ------------------------------- ###

# --- Figure Paired butterfly-host [CG]  
# teneral adult CG data
#adult<-read.csv("char.adult.csv")
#names(adult)
#hc<-adult[,c('species','pair','ug.mg')]
### host plant data CG data
#host<-read.csv("char.host.csv")
#names(host)
#pass<-host[,c('species','pair','ug.mg')]
### combine the host and adult DFs
#cn<-rbind.data.frame(hc,pass)
#cn$ug.mg<-as.numeric(as.character(cn$ug.mg))
#cn<-na.omit(cn)

# summarize stats for bargraph
#stats <- data_summary(cn, varname="ug.mg", 
#                      groupnames=c("species","pair"))

# order the factor levels by increasing CG in host plants
#stats$pair = factor(stats$pair,
#                    levels= c("P.lutea","P.suberosa",
#                              "P.affinis","P.biflora"),
#                    labels = c("P. lutea","P. suberosa",
#                               "P. affinis","P. biflora"))
# plot it
#butter<-ggplot(stats,aes(y=ug.mg, x=pair,fill=species)) + 
#  geom_bar(position="dodge", 
#           stat="identity",
#           col="black",
#           width=0.7,
#           lwd=0.75) +
#  geom_errorbar(aes(ymin=ug.mg-se, ymax=ug.mg+se),
#                width=0.4,
#                position=position_dodge(0.7),
#                col='black') +
#  ylab(expression(paste("Cyanogenic glycosides (", mu, "g/mg)"))) +
#  xlab("Host plant") +
#  scale_fill_manual(values=c("white", "darkgrey",
#                             "darkgrey", "darkgrey",
#                             "darkgrey"),
#                    name="",
#                    labels=c(expression(italic("H. charitonia")), 
#                             "new leaves")) +
#  geom_hline(yintercept = 0,lwd=0.5,col="grey10") +
#  theme(panel.grid.major = element_blank(), 
#        panel.grid.minor = element_line(color="black"),
#        panel.background = element_blank(),
#        plot.title = element_text(vjust=0,hjust=0.5,size = 16, face = "bold", colour = "black"),
#        legend.title = element_text(size = 16, face = "bold", colour = "black"),
#        legend.text = element_text(size = 16, colour = "black"),
#        axis.line = element_line(colour = "black"),
#        axis.text.y = element_text(hjust=1,size=16),
#        axis.text.x = element_text(angle = 50,hjust=1,size=14, face='italic'),
#        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16,vjust=-1),
#        axis.title.x = element_text(margin = margin(t = 20),size=16),
#        plot.margin = margin(t = 10, r = 0, b = 0, l = 0, unit = "pt"))
# butter


----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#################################################################
#################################################################
############################ ANALYSES ########################### 
#################################################################
#################################################################

### ----------------------------------------------------- ### 
### --- Did adults from some hosts acquire more CGs?  --- ###
### ----------------------------------------------------- ### 
adult<-read.csv("char.adult.csv")
hc<-adult[,c('number','species','pair','ug','ug.vol.correct','ug.mg')]
hc$pair<-as.factor(hc$pair)

ADULTstats <- data_summary(hc, varname="ug.vol.correct", 
                                   groupnames=c("pair","species"))
ADULTstats
#        pair       species ug.vol.correct        sd       se
# 1  P.affinis H.charithonia       419.4452 268.93315 56.07644
# 2  P.biflora H.charithonia       131.6374  79.96475 24.11028
# 3    P.lutea H.charithonia       448.6631 371.88945 99.39164
# 4 P.suberosa H.charithonia       380.0485 158.94995 45.88490

hc<-hc[!is.na(hc$ug.vol.correct),]
# check for normality
Lm1=lm(ug.vol.correct~pair,data=hc)
plot(hist(resid(Lm1))) #  normal distribution 
# check for homoscedasticity 
boxplot(ug.vol.correct~pair,data=hc) # Some heterogeneity of variance 

palatable<-lme(ug.vol.correct~pair,weights=varIdent(form=~1|as.factor(pair)),
               random = ~1|number, method = "ML",data=hc)
r.squaredGLMM(palatable) 
#            R2m       R2c
# [1,] 0.6912613 0.7234091
anova(palatable)
#             numDF denDF   F-value p-value
# pair            3    56  14.60833  <.0001

summary(glht(palatable, linfct=mcp(pair="Tukey")))
#                             Estimate Std. Error z value Pr(>|z|)    
# P.biflora - P.affinis == 0   -287.81      59.47  -4.840  < 0.001 ***
# P.lutea - P.affinis == 0       29.22     110.37   0.265  0.99295    
# P.suberosa - P.affinis == 0   -39.40      70.27  -0.561  0.93926    
# P.lutea - P.biflora == 0      317.03      98.50   3.219  0.00653 ** 
# P.suberosa - P.biflora == 0   248.41      49.58   5.010  < 0.001 ***
# P.suberosa - P.lutea == 0     -68.61     105.37  -0.651  0.90883 

### adult sample sizes 
adult<-read.csv("char.adult.csv") # 
adult<-adult[!is.na(adult$ug.vol.correct),]

aff=subset(adult,pair=="P.affinis") # 23
nrow(aff)
bif=subset(adult,pair=="P.biflora") # 11
lut=subset(adult,pair=="P.lutea") # 14
sub=subset(adult,pair=="P.suberosa") # 12


### ------------------------------------------------------------ ### 
### --- Correlation between adult and host CG concentration? --- ###
### ------------------------------------------------------------ ### 
# Host-Adult concentrations
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
# 8    P.suberosa P.suberosa 0.27816993 0.02826091 0.008158221

cor.test(concens[1:4,]$ug.mg, concens[5:8,]$ug.mg)
# t = 0.50439, df = 2, p-value = 0.6641
# cor = 0.3359322 

### how much higher were adults that their hosts?
#         species     ug.mg 
# 1 H.charithonia   4.83406502 
#   P.affinis       0.3268708
4.83406502/0.3268708 # = 14.78892
# 2 H.charithonia   1.43054584 
#   P.biflora       0.22399147
1.43054584 /  0.22399147 # = 6.386609
# 3 H.charithonia   4.89530844 
#   P.lutea         0.27816993
4.89530844/0.27816993 # = 17.59827
# 4 H.charithonia   4.86901195 
#   P.suberosa      0.27816993 
4.86901195 / 0.27816993 # = 17.50373


### ---------------------------------------------------------- ### 
### --- Did larvae from different hosts acquire more CGs?  --- ###
### ---------------------------------------------------------- ###
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

# P. lutea  
lut=subset(cats,host=="lut")
l.sum<-lm(ug.vol.correct~mass.mg.,data=lut)
summary(l.sum) # p-value: 2.024e-07; Adjusted R-squared: 0.6392

# P. biflora 
bif=subset(cats,host=="bi")
b.sum<-lm(ug.vol.correct~mass.mg.,data=bif)
summary(b.sum) # p-value: 0.000165; Adjusted R-squared: 0.4044

# P. affinis 
aff=subset(cats,host=="aff")
a.sum<-lm(ug.vol.correct~mass.mg.,data=aff)
summary(a.sum) # p-value: 1.486e-12; Adjusted R-squared: 0.8423



############################################## 
### Models mass x CGs at different instars ###
############################################## 
# linear model to check for normality
Lm1=lm(log.ug.correct~mass.mg.,data=cats)
plot(hist(resid(Lm1))) # fairly normal distribution 
# check for homoscedasticity 
boxplot(log.ug.correct~host, data = cats) # Some heterogeneity of variance 
#
cats<-cats[!is.na(cats$ug.vol.correct),]
gof<-lm(ug.vol.correct~mass.mg.*host,data=cats)
anova(gof)
#           Df  Sum Sq Mean Sq F value    Pr(>F)    
# mass.mg.       1 1257428 1257428 129.178 < 2.2e-16 ***
# host           2  570341  285170  29.296 3.033e-10 ***
# mass.mg.:host  2  346393  173196  17.793 4.175e-07 ***
r.squaredGLMM(gof) 
#            R2m       R2c
# [1,] 0.7267011 0.7267011


### 3rd instar 
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
summary(glht(i5, linfct=mcp(host="Tukey")))
#                Estimate Std. Error z value Pr(>|z|)  
# bi - aff == 0   -275.20      77.95  -3.530  0.00116 **
# lut - aff == 0   -28.25      79.71  -0.354  0.93310   
# lut - bi == 0    246.96      79.71   3.098  0.00543 **



### ----------------------------- ### 
### ---  exploratory analyses --- ###
### ----------------------------- ### 
# --- Did larval weights vary on different hosts?
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

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################