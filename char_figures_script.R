########################################################################
########################################################################
########################################################################
#### --------------------------------------------------------------- ###
### ---------------- C. Morrison, C. Nguyen & L. Gilbert ----------- ###
### --- Heliconius charithonia cyanogenic glycoside concentation --- ###
### --------------------------- FIGURES ---------------------------- ###
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


############################################################################## 
### ---///  FIGURE 2: Larval and adult cyanide content comparisons  ///--- ###
############################################################################## 
### ------------------------------------------- ###
### --- Panel A: larvae reared on each host --- ###
### ------------------------------------------- ###
# make master database with all hosts
lutea=read.csv("lut.csv")
affinis=read.csv("aff.csv")
biflora=read.csv("bi.csv")
cats<-rbind.data.frame(lutea[,1:15],affinis[,1:15],biflora[,1:15])
head(cats)

# --- three-paneled scatter plot of mass vs CNs 
# set graphics parmeters
cats$host = factor(cats$host,
                   levels= c("lut","aff","bi"),
                   labels = c("Passiflora lutea","Passiflora affinis",
                              "Passiflora biflora"))
cats$life.stage = factor(cats$life.stage)

# fill colors
instar.colors<-c("white","grey75","grey50")

# plot it 
develop<-ggscatter(cats, x = "mass.mg.", y = "ug.vol.correct",
                   color = "black",
                   fill = "life.stage",
                   shape = 21,
                   stroke=1,
                   size=2.5,
                   add = "reg.line",                     
                   conf.int = TRUE,                     
                   add.params = list(color = "blue",
                                     fill = "grey75")) +
  facet_grid(. ~ host) +
  panel_border() +
  xlab('Larval mass (mg)') +
  ylab(expression(paste("CN content (", mu, "g)"))) +
  labs(fill = "Larval instar") +
  scale_fill_manual(values = instar.colors) + 
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(panel.grid.major.x = element_line(colour = "grey90"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size = 16, colour = "black",face = "bold"),
        legend.text = element_text(size = 14),
        legend.position="right",
        strip.text = element_text(face = "italic",size=14),
        strip.background = element_rect(
          color="black", fill="white", size=1.0, linetype="solid"),panel.spacing = unit(0.5, "lines"),
        axis.text.x = element_text(size=14, colour = "grey30"),
        axis.text.y = element_text(size=14, colour = "grey30"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16,vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16,vjust=3),
        plot.margin = margin(t = 20, r = 5, b = -10, l = 0, unit = "pt")
  )

# print this in the dev.new window for a PDF
dev.new(
  title = "CN-mass",
  width = 10,
  height = 6,
  noRStudioGD = TRUE
)

develop

### --------------------------------------------------- ###
### --- Panel B: compare CN content larvae & adults --- ###
### --------------------------------------------------- ###
### amount of adult CN per host association
adult<-read.csv("char.adult.csv")
hc<-adult[,c('species','pair','ug','ug.vol.correct','ug.mg')]

# summary data table
HCstats <- data_summary(hc, varname="ug.vol.correct", 
                        groupnames=c("pair","species"))
names(HCstats)<-c("host","life.stage","ug.total","sd","se")
HCstats[HCstats == "H.charithonia"] <- "adult"
HCstats = HCstats[(HCstats$host !="P.suberosa"),] # remove P. suberosa data
HCstats

### amount of larval CN per host association
biflora=read.csv("bi.csv")
affinis=read.csv("aff.csv")
lutea=read.csv("lut.csv")
cats<-rbind.data.frame(lutea[,1:15],affinis[,1:15],biflora[,1:15])

# summary data table
CATstats <- data_summary(cats, varname="ug.vol.correct", 
                         groupnames=c("host","life.stage"))
names(CATstats)<-c("host","life.stage","ug.total","sd","se")
CATstats[CATstats == "aff"] <- "P.affinis"
CATstats[CATstats == "bi"] <- "P.biflora"
CATstats[CATstats == "lut"] <- "P.lutea"
#CATstats[CATstats == "sub"] <- "P.suberosa"
CATstats

# combine these DFs
comparison<-rbind.data.frame(CATstats,HCstats)
comparison

# set graphics color and order parameters 
comparison$host = factor(comparison$host,
                         levels= c("P.lutea","P.affinis","P.biflora"),
                         labels = c("P. lutea","P. affinis","P. biflora"))

# set colors 
CN.colors<-c("white","grey75","grey50","grey30")

# plot it 
allCNs<-ggplot(comparison, aes(x=host, y=ug.total,fill=life.stage)) + 
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
  ylab(expression(paste("CN content (", mu, "g)"))) +
  labs(fill = "Life stage") +
  geom_hline(yintercept = 0,lwd=0.5,col="grey10") +
  scale_fill_manual(values = CN.colors) + 
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
  title = "CN-stages",
  width = 10,
  height = 10,
  noRStudioGD = TRUE
)
allCNs

### --------------------------------------------------------- ###
### --- Panel C: compare CN concentration larvae & adults --- ###
### --------------------------------------------------------- ###
adult<-read.csv("char.adult.csv")
hc<-adult[,c('species','pair','ug','ug.vol.correct','ug.mg')]
head(hc)

# summary data table of concentrations 
HCconcen <- data_summary(hc, varname="ug.mg", 
                         groupnames=c("pair","species"))
names(HCconcen)<-c("host","life.stage","ug.mg","sd","se")
HCconcen[HCconcen == "H.charithonia"] <- "adult"
HCconcen = HCconcen[(HCconcen$host !="P.suberosa"),] # remove P. suberosa data
HCconcen

# using summary stats on amount of larval CN per host association
cats<-rbind.data.frame(lutea[,1:15],affinis[,1:15],biflora[,1:15])
head(cats)

# summary data table
CATconcen <- data_summary(cats, varname="ug.mg", 
                          groupnames=c("host","life.stage"))
names(CATconcen)<-c("host","life.stage","ug.mg","sd","se")
CATconcen[CATconcen == "aff"] <- "P.affinis"
CATconcen[CATconcen == "bi"] <- "P.biflora"
CATconcen[CATconcen == "lut"] <- "P.lutea"
CATconcen

# combine these DFs
comparison2<-rbind.data.frame(CATconcen,HCconcen)
comparison2

# set graphics color and order parameters 
comparison2$host = factor(comparison2$host,
                          levels= c("P.lutea","P.affinis","P.biflora"),
                          labels = c("P. lutea","P. affinis","P. biflora"))

# set colors 
con.colors<-c("white","grey75","grey50","grey30")

# plot it 
CNconcenBAR<-ggplot(comparison2, aes(x=host, y=ug.mg,fill=life.stage)) + 
  geom_bar(width=0.75,
           stat="identity",
           colour = "black",
           position=position_dodge(),
           lwd=0.75) +
  geom_errorbar(aes(ymin=ug.mg-se, ymax=ug.mg+se),
                width=0.4,
                position=position_dodge(0.75),
                col='black') +
  ggtitle("") +
  xlab("Host plant") + 
  ylab(expression(paste("CN concentration (", mu, "g/mg)"))) +
  labs(fill = "Life stage") +
  geom_hline(yintercept = 0,lwd=0.5,col="grey10") +
  scale_fill_manual(values = con.colors) + 
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
  title = "CN-stages",
  width = 8,
  height = 6,
  noRStudioGD = TRUE
)
CNconcenBAR


# --- Combine into a 2-panel figure 
Fig2=ggarrange(develop,allCNs,CNconcenBAR,
               labels = c("(a)", "(b)", "(c)"),
               ncol = 1, nrow = 3,
               font.label = list(size = 25))
# dev plot
dev.new(
  title = "Figure2",
  width = 9.5,
  height = 15,
  noRStudioGD = TRUE
)
Fig2


############################################################################################# 
### ---/// SUPPLEMENTAL FIGURE 1: LME frequency distributions & shapiro-wilk test  ///--- ###
#############################################################################################
adult<-read.csv("char.adult.csv")
hc<-adult[,c('species','pair','ug','ug.vol.correct','ug.mg')]

#par(mfrow=c(2,1))
### --- (a) CN content LME --- ###
# - make graphable object 
palatable=gls(ug.vol.correct~pair+sex,data=hc, 
              weights=varIdent(form=~1|as.factor(pair)))
shapiro.test(resid(palatable)) # W = 0.96609, p-value = 0.09372 --> normally distributed

# - plot of normal distribution
dev.new(
  title = "(a)",
  width = 6,
  height = 6,
  noRStudioGD = TRUE
)

plot(hist(resid(palatable),breaks = 6),
     main="(a) GLS: Total CN (ug) ~ host + sex",
     xlab="residuals",ylab="frequency",
     xlim=c(-600,700),ylim=c(0,20))
text(x=545,y=16,label=expression(bold(paste("Shapiro-Wilk"))))
text(x=545,y=15,label=expression(paste(italic('P'), " = 0.09")))
text(x=545,y=14,label="W = 0.97")

### --- (b) CN concentration LME --- ###
adults<-adult[!is.na(adult$ug.mg),]

# - make graphable object 
CONpalatable<-gls(ug.mg~pair+sex,data=adults,
                  weights=varIdent(form=~1|as.factor(pair)))
shapiro.test(resid(CONpalatable)) # W = 0.95513, p-value = 0.0273 --> approaching normal distribution

# - plot of normal distribution
dev.new(
  title = "(b)",
  width = 6,
  height = 6,
  noRStudioGD = TRUE
)

plot(hist(resid(CONpalatable),breaks = 6),
     main="(b) GLS: CN (ug/mg) ~ host + sex",
     xlab="residuals",ylab="frequency",ylim=c(0,20))
text(x=8,y=16,label=expression(bold("Shapiro-Wilk")))
text(x=8,y=15,label=expression(paste(italic('P'), " = 0.03")))
text(x=8,y=14,label="W = 0.96")



############################################################################################### 
### ---/// SUPPLEMENNTAL FIGURE 2: Adult (+ P. suberosa) CN content & concentration  ///--- ###
###############################################################################################
adult<-read.csv("char.adult.csv") # 
hc<-adult[,c('species','pair','ug','ug.vol.correct','ug.mg')]

### --- Panel A: Adcult total cyanide amounts --- ###
# summary data table of CN content 
HCcontent <- data_summary(hc, varname="ug.vol.correct", 
                         groupnames=c("pair","species"))
names(HCcontent)<-c("host","life.stage","ug.mg","sd","se")
HCcontent

# set graphics parameters 
HCcontent$host = factor(HCcontent$host,
                       levels= c("P.lutea","P.affinis","P.suberosa","P.biflora"),
                       labels = c("P. lutea","P. affinis","P. suberosa","P. biflora"))
# plot it 
CNcontentBAR<-ggplot(HCcontent, aes(x=host, y=ug.mg)) + 
  geom_bar(width=0.55,
           stat="identity",
           colour = "black",
           position=position_dodge(),
           lwd=0.75) +
  geom_errorbar(aes(ymin=ug.mg-se, ymax=ug.mg+se),
                width=0.3,
                position=position_dodge(0.75),
                col='black') +
  ggtitle("") +
  xlab("") + 
  scale_y_continuous(name=expression(paste("Cyanide content (", mu, "g)")),
                     breaks=seq(0,500,100)) +
  geom_hline(yintercept = 0,lwd=0.5,col="grey10") +
  scale_fill_manual(values = "grey30") + 
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
        plot.margin = margin(t = 10, r = 5, b = 0, l = 0, unit = "pt")
  )
CNcontentBAR


### --- Panel B: Adult concentrations --- ###
# summary data table of CN concentrations 
HCconcen <- data_summary(hc, varname="ug.mg", 
                         groupnames=c("pair","species"))
names(HCconcen)<-c("host","life.stage","ug.mg","sd","se")
HCconcen

# set graphics color and order parameters 
HCconcen$host = factor(HCconcen$host,
                       levels= c("P.lutea","P.affinis","P.suberosa","P.biflora"),
                       labels = c("P. lutea","P. affinis","P. suberosa","P. biflora"))

# plot it 
CNconcenBAR<-ggplot(HCconcen, aes(x=host, y=ug.mg)) + 
  geom_bar(width=0.55,
           stat="identity",
           colour = "black",
           position=position_dodge(),
           lwd=0.75) +
  geom_errorbar(aes(ymin=ug.mg-se, ymax=ug.mg+se),
                width=0.3,
                position=position_dodge(0.75),
                col='black') +
  ggtitle("") +
  xlab("Larval host plant") + 
  scale_y_continuous(name=expression(paste("Cyanide concentration (", mu, "g/mg)")),
                     breaks=seq(0,6,1)) +
  geom_hline(yintercept = 0,lwd=0.5,col="grey10") +
  scale_fill_manual(values = "grey30") + 
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
        plot.margin = margin(t = 10, r = 5, b = 0, l = 0, unit = "pt")
  )
CNconcenBAR

# --- Combine into a 2-panel figure 
FigSS3=ggarrange(CNcontentBAR,CNconcenBAR,
               labels = c("(a)", "(b)"),
               ncol = 1, nrow = 2,
               font.label = list(size = 25))
# dev plot
dev.new(
  title = "FigureSS3",
  width = 6,
  height = 11,
  noRStudioGD = TRUE
)
FigSS3


######################################################################### 
### ---///  Graphic Image: Larval cyanide content comparisons  ///--- ###
######################################################################### 
lutea=read.csv("lut.csv")
affinis=read.csv("aff.csv")
biflora=read.csv("bi.csv")
cats<-rbind.data.frame(lutea[,1:15],affinis[,1:15],biflora[,1:15])
head(cats)

# set graphics parameters
cats$host = factor(cats$host,
                   levels= c("lut","aff","bi"),
                   labels = c("Passiflora lutea","Passiflora affinis",
                              "Passiflora biflora"))
cats$life.stage = factor(cats$life.stage)

# fill colors
instar.colors<-c("white","grey75","grey50")

# plot it 
graphic<-ggscatter(cats, x = "mass.mg.", y = "ug.vol.correct",
                   color = "black",
                   fill = "life.stage",
                   shape = 21,
                   stroke=1,
                   size=2.5,
                   add = "reg.line",                     
                   conf.int = TRUE,                     
                   add.params = list(color = "blue",
                                     fill = "grey75")) +
  facet_grid(. ~ host) +
  panel_border() +
  xlab('Caterpillar mass (mg)') +
  ylab(expression(paste("Cyanides (", mu, "g)"))) +
  guides(fill = FALSE) +
  labs(fill = "Life stage") +
  scale_fill_manual(values = instar.colors) + 
  guides(colour = guide_legend(override.aes = list(size=3))) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position="right",
        strip.text = element_text(face = "italic",size=14),
        strip.background = element_rect(
          color="black", fill="white", size=1.0, linetype="solid"),panel.spacing = unit(0.5, "lines"),
        axis.text.x = element_text(size=14, colour = "grey30"),
        axis.text.y = element_text(size=14, colour = "grey30"),
        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16,vjust=-1),
        axis.title.x = element_text(margin = margin(t = 20,r = 20,b = 20, l = 20),size=16,vjust=3),
        plot.margin = margin(t = 20, r = 5, b = -10, l = 0, unit = "pt")
  )

# print this in the dev.new window for a PDF
dev.new(
  title = "graphic",
  width = 10,
  height = 6,
  noRStudioGD = TRUE
)

graphic


### ------------------------------- ### 
### --- EXTRA exploratory plots --- ###
### ------------------------------- ###
# --- Figure BOXPLOT of life stages CN concentrations on each host
### using raw adult butterfly CN concentration data
#adult<-read.csv("char.adult.csv")
#hc2<-adult[,c('species','pair','ug','ug.vol.correct','ug.mg')] # select columns I need
#hc2 = hc2[(hc2$pair !="P.suberosa"),] # remove P. suberosa data
#hc2[hc2 == "H.charithonia"] <- "adult"
#names(hc2)<-c("life.stage","host","ug","ug.vol.correct","ug.mg")
#head(hc2)
### using raw larval CN concentration data
#biflora=read.csv("bi.csv")
#affinis=read.csv("aff.csv")
#lutea=read.csv("lut.csv")
# cats<-rbind.data.frame(lutea[,1:15],affinis[,1:15],biflora[,1:15])
#cats2 <-cats[,c('host','life.stage','ug','ug.vol.correct','ug.mg')] # select columns I need
#cats2 <- cats2[, c(2, 1, 3, 4, 5)]
#cats2[cats2 == "aff"] <- "P.affinis"
#cats2[cats2 == "bi"] <- "P.biflora"
#cats2[cats2 == "lut"] <- "P.lutea"
#head(cats2)

# combine these DFs
#comparison3<-rbind.data.frame(hc2,cats2)
#comparison3

# set graphics color and order parameters 
#comparison3$host = factor(comparison3$host,
#                          levels= c("P.lutea","P.affinis","P.biflora"),
#                          labels = c("P. lutea","P. affinis","P. biflora"))
# dot fill colors 
#con.colors<-c("white","grey75","grey50","grey25")
# plot it 
#CNconcenBOX<-ggplot(comparison3, 
#                    aes(x=host, y=ug.mg,fill=life.stage)) + 
#  geom_boxplot(position=position_dodge(width=0.75),
#               outlier.shape = NA) +
#  geom_point(position=position_jitterdodge(),
#             pch=21,size=1.25) +
#  ggtitle("") +
#  xlab("Host plant") + 
#  ylab(expression(paste("Cyanide concentration (", mu, "g/mg)"))) +
#  labs(fill = "Life stage") +
#  geom_hline(yintercept = 0,lwd=0.5,col="grey10") +
#  scale_fill_manual(values = con.colors) + 
#  theme(panel.grid.major.y = element_line(colour = "grey90"),
#        panel.grid.major.x = element_blank(),
#        panel.grid.minor.y = element_line(colour = "grey90"),
#        panel.grid.minor.x = element_blank(),
#        panel.background = element_blank(),
#        plot.title = element_text(vjust=0,hjust=0.5,size = 16, face = "bold", colour = "black"),
#        legend.title = element_text(size = 16, face = "bold", colour = "black"),
#        legend.text = element_text(size = 14),
#        axis.line = element_line(colour = "black"),
#        axis.text.y = element_text(hjust=1,size=14, colour = "grey30"),
#        axis.text.x = element_text(angle = 30,hjust=1,size=14, 
#                                   face='italic', colour = "grey30"),
#        axis.title.y = element_text(margin = margin(t = 20,r = 20,b = 20, l = 0),size=16,vjust=-1),
#        axis.title.x = element_text(margin = margin(t = 20),size=16,vjust=3),
#        plot.margin = margin(t = 0, r = 5, b = 20, l = 30, unit = "pt")
#  )
#CNconcenBOX


# --- Figure Paired butterfly-host [CN]  
# teneral adult CN data
#adult<-read.csv("char.adult.csv")
#names(adult)
#hc<-adult[,c('species','pair','ug.mg')]
### host plant data CN data
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

# order the factor levels by increasing CN in host plants
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


########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################