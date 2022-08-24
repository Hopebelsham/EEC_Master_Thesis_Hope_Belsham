### ACCOUNTING FOR STORAGE TIME ####

rm(list=ls())
setwd("/Users/hopebelsham/Desktop/Telomere Project/Data/")
getwd()
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(sqldf)

###########################  Telomere Data  ###################################

# # Reviewing the Telomere Data
# # Checking unique data points in data set
# # Make id column
# data$UIDkey <- paste0(data$BroodRef, data$SocialID, data$BroodYear, data$RearingBrood, 
#                       data$Parent, data$TSRatio, data$FemaleFeedCount, data$MaleFeedCount)
# # print(nrow(data))
# data_unique <- data %>% distinct(UIDkey, .keep_all = TRUE) %>% select
# (-c(Sex, OffspringNo, TValue, SValue, UIDkey))
# 
# # print(nrow(data_unique))
# 
# data_unique$TelomereAttRate <- data_unique

###############################################################################

# Reading in the data
td <- read.csv("td.csv")
# td.no <- sqldf("select distinct(BirdID) from td")
td$CaptureDate <- as.Date(td$CaptureDate, "%d/%m/%Y")
tplot1 <- plot(td$TSRatio)
tplot1 <- abline(tplot1)

#### Accounting for storage time #####
ModelTS <- lm(TSRatio ~ BloodAge + DNAAge, data = td)
td$TSRatioRes <- summary(ModelTS)$residuals
plot(td$TSRatioRes)

# Creating a reduced Telomere data frame, setting variables
td_calcd_2 <- data.frame(BirdID = as.numeric(c()), EndDate = as.Date(c()), 
                       AttrRate = as.numeric(c()), AttrTimeDays = as.numeric(c()),
                       AttrRateYearly = as.numeric(c()), DeathAgeY = as.numeric(c()))

# Calculating attrition rates with new residuals

for (id in unique(td$BirdID)){
  # for (id in x){
  # Subset individual bird data
  df <- subset(td, td$BirdID == id)
  
  # Select useful rows
  df <- df %>% select(TSRatioRes, CaptureDate) %>% arrange(CaptureDate)
  
  # Offset start/end dates
  df1 <- head(df, nrow(df)-1)
  df2 <- tail(df, nrow(df)-1)
  colnames(df1) <- c("StartTSRatio", "StartDate")
  colnames(df2) <- c("EndTSRatio", "EndDate")
  
  # Now pair them up
  df_calc <- cbind(df1, df2)
  df_calc$BirdID = id
  
  # Calc Attrition Rate and time over which this was calc'd
  df_calc$AttrRate <- df_calc$StartTSRatio - df_calc$EndTSRatio
  df_calc$AttrTimeDays <- as.double.difftime(df_calc$EndDate - df_calc$StartDate)
  df_calc$AttrRateYearly <- (df_calc$AttrRate / df_calc$AttrTimeDays) * 365
  
  # Remove extra columns and put back into final df
  df_calc <- df_calc %>% select(BirdID, EndDate, AttrRate, AttrTimeDays, AttrRateYearly)
  td_calcd_2 <- rbind(td_calcd_2, df_calc)
  
}

# Removing all Inf's from data set, then subsetting to remove all values under 50
# time after brooding - needs to be relatively long to see an effect 
td_calcd_2 <- subset(td_calcd_2, is.finite(td_calcd_2$AttrRateYearly))
td_calcd_2 <- subset(td_calcd_2, td_calcd_2$AttrTimeDays > 50)

# Check to see lengthening telomeres vs shortening
print(nrow(subset(td_calcd_2, td_calcd_2$AttrRateYearly < 0)))  # Lengthening 306
print(nrow(subset(td_calcd_2, td_calcd_2$AttrRateYearly > 0)))  # Shortening 311

# Checking the means of the dates for both S and L to see if they are reliable points
mean(subset(td_calcd_2, td_calcd_2$AttrRateYearly < 0)$AttrTimeDays) # 360.5
mean(subset(td_calcd_2, td_calcd_2$AttrRateYearly > 0)$AttrTimeDays) # 408.6

# Exporting telomere attrition rates as a CSV 
# write.csv(td_calcd_2,"/Users/hopebelsham/Desktop/Telomere Project/Data//Telomere_Attrition_Rates.csv", row.names = FALSE)

# Plotting TAR data
p <- ggplot(td_calcd_2, aes(x=AttrTimeDays, y=AttrRateYearly)) + geom_point(alpha=0.125)
print(p)

# sqldf("select distinct(BirdID) from td_calcd_2")
#N=415


########################  Telomere Elongation Stats  ################################


# Analysis to distinguish actual lengthening of telomeres and error in measuring 
td_ana_RES <- td %>% select(BirdID, CaptureDate,TSRatioRes)
td_ana_RES$CaptureDate <- as.Date(td_ana_RES$CaptureDate, "%d/%m/%Y")
td$BirdID <- as.numeric(td$BirdID)

# Creating a data frame with Birds that have 3 telomere measurements
count <- td_ana_RES %>% group_by(BirdID) %>% summarise(count = n())
td_ana_F <- merge(td_ana_RES, count, by="BirdID")
td_MSStats_F <- subset(td_ana_F, td_ana_F$count >=3)

sqldf("select distinct(BirdID) from td_MSStats_F")


# write.csv(td_MSStats_2,"/Users/hopebelsham/Desktop/Telomere Project/Data//ElongationStatsRes.csv", row.names = FALSE)
# Made csv appropriate for TEA

data=read.csv("ESRES.csv")
data <- data %>% select(V1, V2, V3)
data <- na.omit(data)
data$V1 <- as.numeric(data$V1)
data$V2 <- as.numeric(data$V2)
data$V3 <- as.numeric(data$V3)
# N = 227


#Estimating error using individual regressions, Equation 4 in main manuscript
matrixresi<-matrix(0,dim(data)[1],1)  # a matrix to put the residuals of the individual regressions in in
x=1:dim(data)[2] #number of columns is the number of timepoints
j=1
while(j<(dim(data)[1]+1)) #loop the individual regressions for the amount of individuals in the dataset
{
  fit<-lm(t(data[j,])~x)  #individual linear regression
  matrixresi[j,]<-sum((residuals(fit))^2)/(length(x)-2) #residual sum of squares
  j<-j+1
}
sigma1<-mean(c(matrixresi)) #average residual sum of squares across the individuals


#Estimating error under the assumption that TL cannot increase over time, Equation 5 in main manuscript
#first we calculate the increases TL increases over time (between first(1) and last(3) timepoint)
deltaTL=data[,3]-data[,1]

#Next we determine which individuals increase in TL between the first and last timepoint
indexTLincreases=which(deltaTL>0)

#We create a new variable including only the data of individual increases
TLincreases=deltaTL[indexTLincreases]
#sigma2
sigma2=0.5*sum(TLincreases^2)/(length(TLincreases))


#compare both estimates of error variance (sigma1 and sigma2), Equation 6 in main manuscript
vratio=sigma2/sigma1
pvalue<-pf(vratio,length(TLincreases)-1,dim(data)[1]-1,lower.tail=F)
print(pvalue) 
#[1] 0.056

mean(TLincreases)

#designating a set of individuals who show TL increases with a set confidence interval, e.g. 97.5%, equation 7.
upperlimit<-((dim(data)[1]-1)*sigma1)/qchisq(0.025,(dim(data)[1]-1)) 
#Note, change 0.025 in the qchisq to change the confidence at which individual increases are determined. 
#This upperlimit is the upper confidence of the variance determined by the individual
#regressions. This variance is the variance of the upper confidence limit of the 
#underlying normally distributed error function. Using this we can look up individual
#TL increases that are at the boundary (with 95% confidence) of this normal distribution 
#(with standard deviation equal to the upper confidence of sqrt(sigma1)

upperTL=qnorm(0.05,0,sd=sqrt(upperlimit),lower.tail=F)

outsideconfindex=which((0.5*TLincreases)>upperTL) #because TL increases are a result
#from the addition of two equal error distributions divide by 2 (i.e. *0.5).

print(indexTLincreases[outsideconfindex]) #initial row of individual in the dataset
#that shows TL increase beyond the set confidence interval


############################## Parental Care Data #########################################

dvdi <- read.csv("dvdi.csv")
dvdr <- read.csv("dvdr.csv")
pc <- merge(dvdi, dvdr, data.frame, by.x="DVDRef", by.y="DVDRef")
pc.1 <- pc %>% select(BroodRef, OffspringNo, EffectiveLength,FemaleFeedCount, MaleFeedCount, Age)
pc.1 <- pc.1 %>% filter(Age %in% c("7", "11"))
str(pc.1)

####  Female Effective Feed Count
pc.1 <- mutate(pc.1, FemaleEffectiveFeedTime =  pc.1$FemaleFeedCount/pc.1$EffectiveLength)

####  Male Effective Feed Count
pc.1 <- mutate(pc.1, MaleEffectiveFeedTime =  pc.1$MaleFeedCount/pc.1$EffectiveLength)

#### Summing the Male and Female Effective Feed Times per Brood 
pc.2 <- pc.1 %>% group_by(BroodRef) %>% 
  summarise(FemaleEffectiveFeedTime = sum(FemaleEffectiveFeedTime, na.rm=TRUE), 
            MaleEffectiveFeedTime = sum(MaleEffectiveFeedTime, na.rm = TRUE))


# Investigating parental care data for both male and female
# min(data$MaleFeedCount) = 1
# max(data$MaleFeedCount) = 78
# min(data$FemaleFeedCount) = 3
# max(data$FemaleFeedCount) = 58


##############################Brood Data ####################################

# Brood Data
bd <- read.csv("bd.csv")
bd.1 <- bd %>% select (BroodRef, SocialDadID, SocialDadCertain, SocialMumID,
                       SocialMumCertain, BroodYear, Eggs, Hatchlings, EggDate,
                       FledgeNperRearingNest, RearingBrood, Cohort)
bd.1$EggDate <- as.Date(bd.1$EggDate, "%d-%b-%y")


# Merging the data sets
# Merging Parental Care Data and Brood Data by Brood Reference

# MERGE 
nd <- merge(pc.2, bd.1, data.frame, by.x="BroodRef", by.y="RearingBrood")
str(nd)

# Filtering nd by social parents known TRUE to ensure validity of data 
nd.2 <- nd %>% filter(SocialMumCertain == TRUE & SocialDadCertain == TRUE)


nd.2$MumAttrRateYearly <- NA
nd.2$DadAttrRateYearly <- NA


# x <- c(241)
# for (rownum in x){
for (rownum in 1:nrow(nd.2)){
  # FOR MUM
  
  # Subset td_calcd by BirdID
  df <- subset(td_calcd_2, BirdID == nd$SocialMumID[rownum])
  
  # Calc eggdate + 13 (as this covers the time a bird spend before hatching into nestlings)
  # Find closest telomere measurement after 13 days
  df$DateDiff <- df$EndDate - (nd$EggDate[rownum]  + 13)
  df <- df %>% arrange(DateDiff)
  
  selected <- subset(df, df$DateDiff >= 0)[1,]
  
  if (is.na(selected$BirdID)){
    # IF NOT NEEDED, COMMENT OUT
    selected <- subset(df, df$DateDiff < 0)
    selected$DateDiff <- abs(selected$DateDiff)
    selected <- selected %>% arrange(DateDiff)
    selected <- selected[1,]
  }
  
  nd.2$MumAttrRateYearly[rownum] <- selected$AttrRateYearly
  
  # If none, use closest before measurement
  # Do this for mumid and dadid and put into respective columns
  
  df <- subset(td_calcd_2, BirdID == nd$SocialDadID[rownum])
  
  # Calc eggdate + 13 (as this covers the time a bird spends as a nestling)
  # Find closest telomere measurement after 13 days
  df$DateDiff <- df$EndDate - (nd$EggDate[rownum] + 13)
  df <- df %>% arrange(DateDiff)
  
  selected <- subset(df, df$DateDiff >= 0)[1,]
  
  if (is.na(selected$BirdID)){
    # IF NOT NEEDED, COMMENT OUT
    selected <- subset(df, df$DateDiff < 0)
    selected$DateDiff <- abs(selected$DateDiff)
    selected <- selected %>% arrange(DateDiff)
    selected <- selected[1,]
  }
  
  nd.2$DadAttrRateYearly[rownum] <- selected$AttrRateYearly
}


# Removing NAs in both M & D column keeping as many data points as possible 
nd.2_non_na <- nd.2%>% filter_at(vars(MumAttrRateYearly,DadAttrRateYearly),any_vars(!is.na(.)))

# Removing BroodRed.y 
nd.2_non_na <- subset(nd.2_non_na, select = -BroodRef.y)


###################   Reproductive Success Measures   #########################

# Investigating data for variables relating to reproductive success measure 

# TCS = Total Clutch Size (maximum number of eggs found in a nest/brood)
# NS = Nest Success (the probability of at least on egg hatching in a nest)
# EG = Egg Survival (the proportion of eggs surviving to hatch time in successful
#                    nests and is calculated as ES = CSH / TCS) CSH = Clutch size 
#                    at hatch
# HS = Hatching success (proportion of eggs that hatch in a successful nest and is 
#                        calculated as HS = FLN / CSH) FLN = Number of fledglings 
#                        raised)


###### Total Clutch Size ######
nd.2_non_na$TotalClutchSize <- nd.2_non_na$Eggs

###### Hatching Success ######
nd.2_non_na$HatchingSuccess <- nd.2_non_na$FledgeNperRearingNest / nd.2_non_na$Hatchlings

FinalDataRes <- nd.2_non_na

# sqldf("select distinct(SocialDadID) from FinalDataRes") #= 139
# sqldf("select distinct(SocialMumID) from FinalDataRes") #= 147

######################## Removing Outlier's from data ###########################

FinalDataRes <- FinalDataRes%>% filter(MumAttrRateYearly >=-4)
FinalDataRes <- FinalDataRes %>% filter(DadAttrRateYearly >=-1.7)

sqldf("select distinct(SocialDadID) from FinalDataRes")
sqldf("select distinct(SocialMumID) from FinalDataRes")

sqldf("select (MumAttrRateYearly) from FinalDataRes")
sqldf("select (DadAttrRateYearly) from FinalDataRes")

sqldf("select distinct(BroodRef) from FinalDataRes")
pcn <- sqldf("select distinct(BroodRef) from pc")
tdn <- sqldf("select distinct(BirdID) from td")


############################ Linear Modeling ###################################

# install.packages("DHARMa")
# install.packages("gridExtra")
# install.packages("MuMIn")
library(MuMIn)
library(DHARMa)
library(ggplot2)
library(dplyr)
library(broom)
library(ggpubr)
library(lme4)
library(lmerTest)
library(scales)
library(gridExtra)


# Plots

hist(FinalDataRes$EggSurvival)
hist(FinalDataRes$HatchingSuccess)
hist(FinalDataRes$TotalClutchSize)
hist(FinalDataRes$FledgeNperRearingNest)
hist(FinalDataRes$FemaleEffectiveFeedTime)
hist(FinalDataRes$MaleEffectiveFeedTime)
hist(FinalDataRes$Hatchlings)
hist(FinalDataRes$DadAttrRateYearly)
hist(FinalDataRes$MumAttrRateYearly)
ggqqplot(FinalDataRes$DadAttrRateYearly)
ggqqplot(FinalDataRes$MumAttrRateYearly)
ggdensity(FinalDataRes$DadAttrRateYearly) 
ggdensity(FinalDataRes$MumAttrRateYearly)

###############################################################################
grid.arrange(p1,p2,p3,p4)
grid.arrange(p5,p6,p7,p8)
grid.arrange(p9,p10,p11,p12)
###############################################################################

# Total Clutch Size - S for Females
FinalDataTCSF.lmer <- lmer(MumAttrRateYearly ~ TotalClutchSize + (1|Cohort), data = FinalDataRes)
summary(FinalDataTCSF.lmer)
FinalDataTCSM.lmer <- lmer(DadAttrRateYearly ~ TotalClutchSize + (1|Cohort), data = FinalDataRes)
summary(FinalDataTCSM.lmer)

r.squaredGLMM(FinalDataTCSF.lmer)
r.squaredGLMM(FinalDataTCSM.lmer)


# Checking fits
Residuals_FinalDataTCSF.lmer  <- simulateResiduals(FinalDataTCSF.lmer)
plot(Residuals_FinalDataTCSF.lmer)
Residuals_FinalDataTCSM.lmer <- simulateResiduals(FinalDataTCSM.lmer)
plot(Residuals_FinalDataTCSM.lmer)

# Plot
 p3 <- ggplot(FinalDataRes, aes(TotalClutchSize, MumAttrRateYearly)) +
  geom_jitter(width=0.1) +
  scale_x_continuous(breaks=seq(from=1, to=7, by=1))+
  scale_y_continuous(breaks=seq(from=-1.25, to=1.25, by=0.25))+
  xlab("Total Clutch Size (No. of eggs per brood)")+ ylab("Female Attriton Rates per Year") +
  geom_text(x=2, y=-0.75, label="P value = 0.052", colour="darkblue", size=8, (aes(fontface=3)))+
  geom_smooth(method = "lm", se = FALSE)+
  theme(axis.title.x=element_text(size=25), axis.title.y =element_text(size=25))+
  theme(panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",
          size = 2, linetype = "solid"), panel.grid.major = element_line
          (size = 0.5, linetype = 'solid', colour = "white"), 
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
          colour = "white"))

p3



p4 <- ggplot(FinalDataRes, aes(TotalClutchSize, DadAttrRateYearly)) +
  geom_point() +
  xlim(1, 7) +
  ylim(-1, 1.2) +
  labs(x="Total Clutch Size", y="Male Attriton Rates per Year",
       title="Total clutch size against male attrtion rates, no outliers") +
  geom_text(x=2, y=0.75, label="P value = 0.852", colour="darkblue")+
  geom_smooth(method = "lm", se = FALSE)  +
  theme(plot.title = element_text((family="Times New Roman"), hjust = 1, size = 16, face = "bold.italic"),
        axis.title.x = element_text((family="Times New Roman"),size = 15),
        axis.title.y = element_text((family="Times New Roman"),size = 15)) +
  theme(panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",
                                        size = 2, linetype = "solid"), panel.grid.major = element_line
        (size = 0.5, linetype = 'solid', colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white"))
p4
###############################################################################

# Fledged young - NS
FinalDataFYF.lmer <- lmer(MumAttrRateYearly ~ FledgeNperRearingNest + (1|Cohort), data = FinalDataRes)
summary(FinalDataFYF.lmer)
FinalDataFYM.lmer <- lmer(DadAttrRateYearly ~ FledgeNperRearingNest + (1|Cohort), data = FinalDataRes)
summary(FinalDataFYM.lmer)

r.squaredGLMM(FinalDataFYF.lmer)
r.squaredGLMM(FinalDataFYM.lmer)

# Checking fits
Residuals_FinalDataFYF.lmer  <- simulateResiduals(FinalDataFYF.lmer)
plot(Residuals_FinalDataFYF.lmer)
Residuals_FinalDataFYM.lmer <- simulateResiduals(FinalDataFYM.lmer)
plot(Residuals_FinalDataFYM.lmer)

# Plot
p5 <- ggplot(FinalDataRes, aes(FledgeNperRearingNest, MumAttrRateYearly)) +
  geom_point() +
  xlim(1, 5) +
  ylim(-1, 1.2) +
  labs(x="Female Attriton Rates per Year", y="Female Attriton Rates per Year",
       title="Fledged young per brood against female attrtion rates, no outliers") +
  geom_text(x=4, y=1, label="P value = 0.218", colour="darkblue")+
  geom_smooth(method = "lm", se = FALSE)  +
  theme(plot.title = element_text((family="Times New Roman"), hjust = 1, size = 16, face = "bold.italic"),
        axis.title.x = element_text((family="Times New Roman"),size = 15),
        axis.title.y = element_text((family="Times New Roman"),size = 15)) +
  theme(panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",
                                        size = 2, linetype = "solid"), panel.grid.major = element_line
        (size = 0.5, linetype = 'solid', colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white"))
p5

p6 <- ggplot(FinalDataRes, aes(FledgeNperRearingNest, DadAttrRateYearly)) +
  geom_point() +
  xlim(1, 5) +
  ylim(-1, 1.2) +
  labs(x="Male Attriton Rates per Year", y="Male Attriton Rates per Year",
       title="Fledged young per brood male attrtion rates, no outliers") +
  geom_text(x=4, y=1, label="P value = 0.359", colour="darkblue")+
  geom_smooth(method = "lm", se = FALSE)  +
  theme(plot.title = element_text((family="Times New Roman"), hjust = 1, size = 16, face = "bold.italic"),
        axis.title.x = element_text((family="Times New Roman"),size = 15),
        axis.title.y = element_text((family="Times New Roman"),size = 15)) +
  theme(panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",
                                        size = 2, linetype = "solid"), panel.grid.major = element_line
        (size = 0.5, linetype = 'solid', colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white"))
p6

###############################################################################

# Effective Feed Time - NS
FinalDataEFTF.lmer <- lmer(MumAttrRateYearly ~ FemaleEffectiveFeedTime + (1|Cohort), data = FinalDataRes)
summary(FinalDataEFTF.lmer)
FinalDataEFTM.lmer <- lmer(DadAttrRateYearly ~ MaleEffectiveFeedTime + (1|Cohort), data = FinalDataRes)
summary(FinalDataEFTM.lmer)

r.squaredGLMM(FinalDataEFTF.lmer)
r.squaredGLMM(FinalDataEFTM.lmer)

# Checking fits
Residuals_FinalDataEFTF.lmer  <- simulateResiduals(FinalDataEFTF.lmer)
plot(Residuals_FinalDataEFTF.lmer)
Residuals_FinalDataEFTM.lmer <- simulateResiduals(FinalDataEFTM.lmer)
plot(Residuals_FinalDataEFTM.lmer)

# Plot
p7 <- ggplot(FinalDataRes, aes(FemaleEffectiveFeedTime, MumAttrRateYearly)) +
  geom_point() +
  xlim(0, 1) +
  ylim(-1, 1.2) +
  labs(x="Effective Feed Time", y="Female Attriton Rates per Year",
       title="Female effective feed time against female attrtion rates, no outliers") +
  geom_text(x=0.80, y=1, label="P value = 0.838", colour="darkblue") +
  geom_smooth(method = "lm", se = FALSE)  +
  theme(plot.title = element_text((family="Times New Roman"), hjust = 1, size = 16, face = "bold.italic"),
        axis.title.x = element_text((family="Times New Roman"),size = 15),
        axis.title.y = element_text((family="Times New Roman"),size = 15)) +
  theme(panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",
                                        size = 2, linetype = "solid"), panel.grid.major = element_line
        (size = 0.5, linetype = 'solid', colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white"))
p7

p8 <- ggplot(FinalDataRes, aes(MaleEffectiveFeedTime, DadAttrRateYearly)) +
  geom_point() +
  xlim(0, 0.8) +
  ylim(-1, 1.2) +
  labs(x="Effective Feed Time", y="Male Attriton Rates per Year",
       title="Male effective feed timemale attrtion rates, no outliers") +
  geom_text(x=0.60, y=1, label="P value = 0.743", colour="darkblue") +
  geom_smooth(method = "lm", se = FALSE)  +
  theme(plot.title = element_text((family="Times New Roman"), hjust = 1, size = 16, face = "bold.italic"),
        axis.title.x = element_text((family="Times New Roman"),size = 15),
        axis.title.y = element_text((family="Times New Roman"),size = 15)) +
  theme(panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",
                                        size = 2, linetype = "solid"), panel.grid.major = element_line
        (size = 0.5, linetype = 'solid', colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white"))
p8
###############################################################################

# Hatching Success - NS
FinalDataHSF.lmer <- lmer(MumAttrRateYearly ~ HatchingSuccess + (1|Cohort), data = FinalDataRes)
summary(FinalDataHSF.lmer)
FinalDataHSM.lmer <- lmer(DadAttrRateYearly ~ HatchingSuccess + (1|Cohort), data = FinalDataRes)
summary(FinalDataHSM.lmer)

r.squaredGLMM(FinalDataHSF.lmer)
r.squaredGLMM(FinalDataHSM.lmer)

# Checking fits
Residuals_FFinalDataHSF.lmer  <- simulateResiduals(FinalDataHSF.lmer)
plot(Residuals_FFinalDataHSF.lmer)
Residuals_FinalDataHSM.lmer <- simulateResiduals(FinalDataHSM.lmer)
plot(Residuals_FinalDataHSM.lmer)

# Plot
p9 <- ggplot(FinalDataRes, aes(HatchingSuccess, MumAttrRateYearly)) +
  geom_point() +
  xlim(0.1, 1) +
  ylim(-1.2, 1.2) +
  labs(x="Hatching Success", y="Female Attriton Rates per Year",
       title="Hatching Success against female attrtion rates, no outliers") +
  geom_text(x=0.90, y=1, label="P value = 0.661", colour="darkblue") +
  geom_smooth(method = "lm", se = FALSE)  +
  theme(plot.title = element_text((family="Times New Roman"), hjust = 1, size = 16, face = "bold.italic"),
        axis.title.x = element_text((family="Times New Roman"),size = 15),
        axis.title.y = element_text((family="Times New Roman"),size = 15)) +
  theme(panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",
                                        size = 2, linetype = "solid"), panel.grid.major = element_line
        (size = 0.5, linetype = 'solid', colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white"))
p9

p10 <- ggplot(FinalDataRes, aes(HatchingSuccess, DadAttrRateYearly)) +
  geom_point() +
  xlim(0.1, 1) +
  ylim(-1.2, 1.2) +
  labs(x="Hatching Success", y="Male Attriton Rates per Year",
       title="Hatching Success against male attrtion rates, no outliers") +
  geom_text(x=0.90, y=1, label="P value = 0.637", colour="darkblue") +
  geom_smooth(method = "lm", se = FALSE)  +
  theme(plot.title = element_text((family="Times New Roman"), hjust = 1, size = 16, face = "bold.italic"),
        axis.title.x = element_text((family="Times New Roman"),size = 15),
        axis.title.y = element_text((family="Times New Roman"),size = 15)) +
  theme(panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",
                                        size = 2, linetype = "solid"), panel.grid.major = element_line
        (size = 0.5, linetype = 'solid', colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white"))
p10
###############################################################################

# Hatchlings 
FinalDataHF.lmer <- lmer(MumAttrRateYearly ~ Hatchlings + (1|Cohort), data = FinalDataRes)
summary(FinalDataHF.lmer)
FinalDataHM.lmer <- lmer(DadAttrRateYearly ~ Hatchlings + (1|Cohort), data = FinalDataRes)
summary(FinalDataHM.lmer)

r.squaredGLMM(FinalDataHF.lmer)
r.squaredGLMM(FinalDataHM.lmer)

# Checking fits
Residuals_FinalDataHF.lmer  <- simulateResiduals(FinalDataHF.lmer)
plot(Residuals_FinalDataHF.lmer)
Residuals_FinalDataHM.lmer <- simulateResiduals(FinalDataHM.lmer)
plot(Residuals_FinalDataHM.lmer)


# Plot
p11 <- ggplot(FinalDataRes, aes(Hatchlings, MumAttrRateYearly)) +
  geom_point() +
  xlim(1, 6) +
  ylim(-1.2, 1.2) +
  labs(x="Hatchlings", y="Female Attriton Rates per Year",
       title="Number of hatchlings against female attrtion rates, no outliers") +
  geom_text(x=5, y=1, label="P value = 0.168", colour="darkblue") +
  geom_smooth(method = "lm", se = FALSE)  +
  theme(plot.title = element_text((family="Times New Roman"), hjust = 1, size = 16, face = "bold.italic"),
        axis.title.x = element_text((family="Times New Roman"),size = 15),
        axis.title.y = element_text((family="Times New Roman"),size = 15)) +
  theme(panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",
                                        size = 2, linetype = "solid"), panel.grid.major = element_line
        (size = 0.5, linetype = 'solid', colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white"))
p11

p12 <- ggplot(FinalDataRes, aes(Hatchlings, DadAttrRateYearly)) +
  geom_point() +
  xlim(1, 6) +
  ylim(-1.2, 1.2) +
  labs(x="Hatchlings", y="Male Attriton Rates per Year",
       title="Number of hatchlings against male attrtion rates, no outliers") +
  geom_text(x=5, y=1, label="P value = 0.112", colour="darkblue") +
  geom_smooth(method = "lm", se = FALSE)  +
  theme(plot.title = element_text((family="Times New Roman"), hjust = 1, size = 16, face = "bold.italic"),
        axis.title.x = element_text((family="Times New Roman"),size = 15),
        axis.title.y = element_text((family="Times New Roman"),size = 15)) +
  theme(panel.background = element_rect(fill = "#BFD5E3", colour = "#6D9EC1",
                                        size = 2, linetype = "solid"), panel.grid.major = element_line
        (size = 0.5, linetype = 'solid', colour = "white"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "white"))
p12
###############################################################################

