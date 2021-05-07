00011111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111100
00#            DESCRIPTION:                                                                            ##
00#                                                                                                    ##
00#             This R script includes processed data, and generates graphs associated with            ##
00#             the manuscript "How tree species, tree size and topographical location influenced      ##
00#             tree transpiration in northern boreal forests during the historic 2018 drought         ##                                                                               
00#                                                                                                    ##
00#             NOTE: This data and associated scripts are an extract from a larger script and         ##
00#                   database, and consequently, the files might contains additional data that        ##
00#                   was not included in the manuscript. Additionally, some data are already fully    ##
00#                   processed, and the R scripts provided are for reference purposes only            ##
00#                                                                                                    ##
00#             Manuscript DOI: https://doi.org/10.1111/gcb.15601                                      ##
00#                                                                                                    ##
00#                                                                                                    ##
00#            SAP FLOW DETAILS:                                                                       ##
00#                 - Method: Heat dissipation. Granier 1985                                           ##
00#                 - Corrections: Wounding drift & Radial profiles                                    ##
00#                                                                                                    ##
00#            WHOLE-TREE CONDUCTANCE METHOD:                                                          ##
00#                 - Phillips & Oren, 1998                                                            ##
00#                                                                                                    ##
00#                                        KRYCKLAN Sites                                              ##
00#                                                                                                    ##
00#      Created by: Jose GL, March 1, 2019                                                            ##
00#       Last update/edit:                                                                            ##
00#                  Mar 15, 2021                                                                      ##
00011111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111100
00011111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111100

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#
#  Files required
#
000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

# LIST OF FILES REQUIRED   ----------------------------------------------------------------

"stk.sf.met.csv"
"daily.data.csv"
"all.nodes.sapflow.csv"
"wound.drift.fd.csv"
"all.nodes.raw.dat"
"longterm_drought.csv"
"all.eco.met.dat"
"allometric_data.dat"
"tree.diam.info.csv"
"par.estimates.csv"
"histo_data.csv"

# END OF SECTION ----

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#
#  LIBRARIES - PACKAGES
#
000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

# INSTALL - LIBRARIES & PACKAGES ----------------------------------------------------------------

library(robustbase)   
library(scPDSI)
library(lubridate)
library(data.table)
library(dplyr)
library(SPEI)
library(rbin)
library(minpack.lm)
library(biwt)
library(lme4)

# SET - WORKING DIRECTORY ----------------------------------------------------------------

# Get working directory and files
getwd() 
setwd("Data/") 


# END OF SECTION ----

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#
#  ENVIRONMENTAL DATA
#
000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

# OPEN - METEOROLOGICAL AND ECOSYSTEMS DATA  ----------------------------------------------------------------

all.eco.met = read.table(file="Data/all.eco.met", header=T, sep = ",")

#Add date to database

all.eco.met["Date"] <- (Date= as_datetime(seq(as.POSIXct("2016/01/01 00:00:00", tz="UTC"), as.POSIXct("2020/12/30 23:30:00", tz="UTC"), by = "30 min")))

rm(Date)
#max(grep("11:30 PM$", eco_met_678$TIMESTAMP_30))
#tail(eco_met_678$TIMESTAMP_30, 1)

# ESTIMATE - REFERENCE ET & VPD  ----------------------------------------------------------------


#Estimate saturation vapor pressure & VPD [both in kPa]

SVP = 0.6108 * exp((17.27*all.eco.met$Ta_1_1_1)/(all.eco.met$Ta_1_1_1+237.3))
VPD = SVP * (1-all.eco.met$RH_1_1_1/100)

#Estimate slope component 

D = (4098* SVP)/(all.eco.met$Ta_1_1_1 + 237.1)^2
plot(D)

# Estimate atmospheric pressure for a given elevation
elev = 275
P = 101.3*((293-0.0065*elev)/293)^5.26

#Estimate psychrometric constant
#Latent heat vaporizaton = 2.45 MJ K-1
#Specific heat at constant pressure = 1.013 10-3 MK K-1 C-1
#Ratio of molecular weight of water vapour/dry air = 0.622

g =  (0.001013*P)/(0.622*2.45)

# Estimate average of all soil heat flux sensors
all.eco.met["G_All"] <-rowMeans(all.eco.met[grep("^G_1|^G_2|^G_3|^G_4", names(all.eco.met))],  na.rm = TRUE)

# Estimate potential evapotranspiration 
# We have two measurements per hour

ET0 = (0.408 * D * (((all.eco.met$NetRad_1_2_1*0.0036)-(all.eco.met$G_All*0.0036)) + g * (18.75 / (all.eco.met$Ta_1_1_1 + 273)) * 0 * VPD)) / (D + g * (1 + 0.24 * 0))

# Air density 
pa = P /(287.058 *all.eco.met$Ta_1_1_1 + 273.15)

###///////// Section to adjust ET0 to a forest (paused at the moment)

rs = 100/0.5*2.95
ra = (log(20-20*(2/3)/20*0.123)*log(20-20*(2/3)/20*0.123))/( 0.41^2*0.5)

#plot(pa, ylim=c(0,10))
plot(all.eco.met$Ta_1_1_1)

ET0.f = D * (((all.eco.met$NetRad_1_2_1*0.0036)-(all.eco.met$G_All*0.0036)) +  (((1.225*0.001013*VPD/ra)))/
               D+g*(1+(rs/ra)))
rm(pa)
###///////// 

# NetRad is in W/m^2 and we need MJ m^-2 h^-1
# 1 watt = 1 Joule/s
# There are 3600 s in an hour
# There are 1 000 000 Joules in a MJ 
# e.g. 1000 watt/m^2 * 3600 s == 3 600 000 J m^-2 h^-1 == 3.6 MJ m^-2 h^-1
# Or I just multiply by 0.0036 and get the same number
# 1000 * 0.0036 = 3.6
# We have two measurements per hour
# To add up all hours in a day, each needs to be divided over two
# Or simply divide the multiplier by two:
# 1000 watt/m^2 * 0.0018 = 1.8 MJ m^-2 (per measurement)
# "Per measurement" units can be easily added up 


all.eco.met["ETo_mm"]<- ifelse(ET0 < 0, 0, ET0/2) 

all.eco.met["ET0.f"]  <- ET0.f/2

all.eco.met["VPD_kPa"]<- VPD

all.eco.met["PP_ET"] <- all.eco.met$P_1_1_1 - all.eco.met$ETo_mm


all.eco.met["DOY_cont"] <- rep(seq(1,NROW(all.eco.met)/48, 1), times=1, each=48)

par(mar=c(4,4, 2, 2))

plot(aggregate(cbind(ET0.f)~DOY_cont, data=all.eco.met, FUN=sum), type="l", lty=1, col="red")
lines(aggregate(cbind(ETo_mm)~DOY_cont, data=all.eco.met, FUN=sum), col="blue")

#ETo_mm <-subset(all.eco.met, select=c("Date", "ETo_mm"))
names(all.eco.met)

#write.table(ETo_mm, "ETo_mm_mar_15_2021.csv", sep=",", row.names=F)

rm(ETo_mm)
rm(D, ET0, SVP, VPD, elev, P, g)
rm(ET0.f, oa, ra, rs, temp, mod.col)

# Testing different VPD formulas (all should match)
function(testrandom_vpd){
  
  
  SVP = 0.6108 * exp(17.27*all.eco.met$Ta_1_1_1/(all.eco.met$Ta_1_1_1+237.3))
  VPD = SVP * (1-all.eco.met$RH_1_1_1/100)
  
  all.eco.met["VPD_kPa_2"]<- VPD
  
  SVP = 610.78 * exp(all.eco.met$Ta_1_1_1 /(all.eco.met$Ta_1_1_1+237.3)*17.2694 ) 
  
  VPD = SVP *(1-all.eco.met$RH_1_1_1/100)/1000
  
  all.eco.met["VPD_kPa_3"]<- VPD
  
  par(mar = c(3, 3, 0, 0))
  
  junecol <- ifelse(all.eco.met$Month==6, "red", "grey")
  
  plot(all.eco.met$Date, all.eco.met$VPD_kPa)
  
  lines(subset(all.eco.met, all.eco.met$Year==2017, select=c("Date", "VPD_kPa_2")), col="red")
  
  lines(all.eco.met$Date, SVP, col="blue")
  
  
  vpd_agg <- aggregate(cbind(VPD_kPa, VPD_kPa_2, VPD_kPa_3)~Hour, 
                       data=all.eco.met,FUN=mean, na.rm=T)
  
  plot(vpd_agg$Hour, vpd_agg$VPD_kPa, col="green")
  
  lines(vpd_agg$Hour, vpd_agg$VPD_kPa_2, col="red")
  lines(vpd_agg$Hour, vpd_agg$VPD_kPa_3, col="blue")
  
  
}


# ADD - HOUR, DOY, DOY_cont, MONTH AND YEAR ----------------------------------------------------------------

NROW(all.eco.met)

names(all.eco.met)

all.eco.met["DOY_cont"] <- rep(seq(1,NROW(all.eco.met)/48, 1), times=1, each=48)
all.eco.met["DOEY"] <- lubridate::day(all.eco.met$Date)

all.eco.met["Month"] <- lubridate::month(all.eco.met$Date)

all.eco.met["Hour"] <- lubridate::hour(all.eco.met$Date)
all.eco.met["Year"] <- lubridate::year(all.eco.met$Date)



# ADD - BINNED VPD & PAR TO - "all.eco.met"   ----------------------------------------------------------------
# Binned are used to explore the data, at reduced computing power, and help identify trends

labels = seq(0,4, 0.02)

all.eco.met["VPD_0.2"] <-  labels[cut(all.eco.met$VPD_kPa, breaks = seq(0,4.02, 0.02),labels = labels,right = TRUE)]

labels = seq(0,2100, 50)
all.eco.met["PAR_50"] <-  labels[cut(all.eco.met$PPFD_IN_1_2_1, breaks = seq(0,2150, 50),labels = labels,right = TRUE)]

rm(labels)


# ESTIMATE - MEAN VOLUMETRIC WATER CONTENT & SOIL SATURATION ----------------------------------------------------------------

names(all.eco.met)

all.eco.met["VWC_5"] <- (rowMeans(all.eco.met[c(29,34,39,44)],  na.rm = TRUE)/100)
plot(all.eco.met$Date, all.eco.met$VWC_5)
all.eco.met$VWC_5[all.eco.met$VWC_5<0.15] <- NA

all.eco.met["VWC_10"] <- (rowMeans(all.eco.met[c(30,35,40,45)],  na.rm = TRUE)/100)

plot(all.eco.met$Date, all.eco.met$VWC_10)
all.eco.met$VWC_10[all.eco.met$VWC_10<0.15] <- NA

all.eco.met["VWC_15"] <- (rowMeans(all.eco.met[c(31,36,41,46)],  na.rm = TRUE)/100)

par(mar=c(4, 4, 0, 0))
plot(all.eco.met$Date, all.eco.met$VWC_15)
all.eco.met$VWC_15[all.eco.met$VWC_15<0.19] <- NA

all.eco.met["VWC_30"] <- (rowMeans(all.eco.met[c(32,37,42,47)],  na.rm = TRUE)/100)

par(mar=c(4, 4, 0, 0))
plot(all.eco.met$Date, all.eco.met$VWC_30)
all.eco.met$VWC_30[all.eco.met$VWC_30<0.14] <- NA

all.eco.met["VWC_50"] <- (rowMeans(all.eco.met[c(33,38)],  na.rm = TRUE)/100)

plot(all.eco.met$Date, all.eco.met$VWC_50)
all.eco.met$VWC_50[all.eco.met$VWC_50<0.1] <- NA


all.eco.met["VWC0.15"] <- rowMeans(all.eco.met[grep("VWC_5$|VWC_10$|VWC_15", names(all.eco.met))],  na.rm = TRUE)

plot(all.eco.met$Date, all.eco.met$VWC0.15)
all.eco.met$VWC0.15[all.eco.met$VWC0.15<0.18] <- NA

all.eco.met["VWC0.30"] <- rowMeans(all.eco.met[grep("VWC_5$|VWC_10$|VWC_15|VWC_30", names(all.eco.met))],  na.rm = TRUE)

plot(all.eco.met$Date, all.eco.met$VWC0.30)
all.eco.met$VWC0.15[all.eco.met$VWC0.15<0.18] <- NA

# Estimate soil saturation from 0-15 cm 

soilsat <- (all.eco.met$VWC0.15-0.19)/(0.31-0.19)
all.eco.met["S0.15"] <- ifelse(soilsat<0, 0,
                               ifelse(soilsat>1, 1, soilsat))

rm(soilsat)

plot(all.eco.met$Date, all.eco.met$S0.15)


# Estimate missing soil saturation using LOESS
# Method works well if only a small fraction of data is missing

all.eco.met["ID"] <- seq(1:NROW(all.eco.met))
x <- all.eco.met$ID
y <- all.eco.met$S0.15

lo <- loess(y ~ x, na.rm=T, family="gaussian", control=loess.control(surface="direct", statistics="approximate"), span=0.01)

all.eco.met["S0.15_loess"] <- ifelse(predict(lo,x,na.rm=T)<0,0, predict(lo,x,na.rm=T))

all.eco.met["S0.15.nomiss"] <- ifelse(is.na(all.eco.met$S0.15), all.eco.met$S0.15_loess, all.eco.met$S0.15)

rm(x, y, lo)

# Estimate missing VPD using LOESS

x <- all.eco.met$ID
y <- all.eco.met$VPD_kPa

lo <- loess(y ~ x, na.rm=T, family="gaussian", control=loess.control(surface="direct", statistics="approximate"), span=0.001)

all.eco.met["VPD_kPa_loess"] <- ifelse(predict(lo,x,na.rm=T)<0, 0, predict(lo,x,na.rm=T))

all.eco.met["S0.15.nomiss"] <- ifelse(is.na(all.eco.met$VPD_kPa), all.eco.met$VPD_kPa_loess, all.eco.met$S0VPD_kPa.15)

plot(all.eco.met$Date, all.eco.met$VPD_kPa_loess,type="l", col="red", ylim=c(0, 3))

lines(all.eco.met$Date, all.eco.met$VPD_kPa, col="blue")

rm(x, y, lo)



par(mar=c(4, 4, 2, 2 ))
plot(all.eco.met$Date, all.eco.met$S0.15_loess, type="l", col="red", lwd=3)
lines(all.eco.met$Date, all.eco.met$S0.15, col="blue", lwd=3)



# ADD - CONTINUOS WEEK ----------------------------------------------------------------

names(all.eco.met)

all.eco.met["Week"] <- ifelse(all.eco.met$Year==2016,lubridate::week(all.eco.met$Date),
                              ifelse((all.eco.met$Year==2017),lubridate::week(all.eco.met$Date)+53,
                                     ifelse((all.eco.met$Year==2018),lubridate::week(all.eco.met$Date)+106,
                                            ifelse((all.eco.met$Year==2019),lubridate::week(all.eco.met$Date)+159,
                                                   lubridate::week(all.eco.met$Date)+212))))

#plot(all.eco.met$Hour, all.eco.met$ETo_mm)
#plot(all.eco.met$Week, all.eco.met$P_1_1_1)
#plot(all.eco.met$Week, all.eco.met$ETo_mm)


# ESTIMATE - OVERCAST DAYS----------------------------------------------------------------
# [Paused, Feb 19, 2021. Might not be needed]

# Aggregates by day and week to extract max ETo and max PAR

names(all.eco.met)

temp_agg <- aggregate(cbind(Week/48, VPD_kPa, ETo_mm, PPFD_IN_1_2_1)~DOY_cont, 
                      data=all.eco.met,FUN=sum, na.rm=T)

temp_agg_week <- aggregate(cbind(V1, VPD_kPa, ETo_mm)~V1, 
                           data=temp_agg,FUN=max, na.rm=T)

temp_week <- subset(all.eco.met, select=c("Week", "DOY_cont"))

colnames(temp_agg)[c(2)]<-c("Week")
colnames(temp_agg_week)[c(2)]<-c("Week")

temp_agg["ETo_max_mm"] <- merge(temp_agg, temp_agg_week, by="Week", all.x=T, all.y=, sort=T )[8]
temp_agg["PPFD_IN_1_2_1_max"] <- merge(temp_agg, temp_agg_week, by="Week", all.x=T, all.y=, sort=T )[5]

names(temp_agg)

x <- temp_agg$DOY_cont
y <- temp_agg$ETo_max_mm

lo <- loess(y~x, family = c("gaussian"), control = loess.control(surface = c("direct"),
                                                                 statistics = c("approximate")),span=0.01)

temp_agg["ETo_loess"] <- summary(lo)$fitted

par(mar=c(4, 4, 2, 2 ))
plot(temp_agg$DOY_cont, temp_agg[,4], type="l", col="red", lwd=3)
lines(temp_agg$DOY_cont, temp_agg[,8], col="blue", lwd=3)

legend = c("Max ETo", "ETo")

legend("topleft",  legend = legend, col = c("blue", "red"), lty=1, lwd=3 
       ,bty="n",ncol=1, cex=1)
rm(legend)
temp_agg["Overcast"] <- ifelse(temp_agg$ETo_mm<temp_agg$ETo_loess *0.9, "Overcast", "Clear")
temp_agg["Overcast.range"] <- temp_agg$ETo_mm/temp_agg$ETo_loess
plot(temp_agg$DOY_cont, temp_agg$Overcast.range)
tmp_col = ifelse(temp_agg$Overcast=="Clear", "Blue", "red")

plot(temp_agg[,1], temp_agg[,4], type="p", col=tmp_col)

temp_week <-  merge(temp_week, temp_agg, by="DOY_cont", all.x=T, all.y=, sort=T)

all.eco.met["ETo_max_loess"]     <- temp_week$ETo_loess
all.eco.met["Overcast"]          <- temp_week$Overcast
all.eco.met["PPFD_IN_1_2_1_max"] <- temp_week$PPFD_IN_1_2_1_max
all.eco.met["Overcast.range"] <- temp_week$Overcast.range

rm(temp_agg, temp_agg_week, temp_week, x, y, surface, temp, lo, tmp_col)


# ESTIMATE: Daily means ----------------------------------------------------------------

names(all.eco.met)


daily.data <- aggregate(cbind(Date)~DOY_cont, 
                        data=all.eco.met, FUN = min, na.action = na.pass)

daily.data$Date <- date(as.POSIXct(daily.data$Date, origin="1970/01/1", tz="UTC"))

all.eco.met["Date_DOY"] <- date(all.eco.met$Date)

daily.means <- aggregate(cbind(Ta_1_1_1, RH_1_1_1, PPFD_IN_1_2_1, VPD_kPa, 
                               VWC_5, VWC_10, VWC_15, VWC_30, VWC_50, VWC0.15, S0.15)~Date_DOY, 
                         data=all.eco.met, FUN = mean, na.action = na.pass)

#daily.means["Date"] <- date(as.POSIXct(daily.means$Date, origin="1970/01/1", tz="UTC"))

daily.data$Date
daily.means$Date

names(daily.data)
names(daily.means)

# Change name of variable for "merge"
colnames(daily.means)[1]<-c("Date")

#daily.data <- merge(daily.data, daily.means, by="Date", all.x=T, all.y=F)

rm(daily.means)

NROW(daily.data)
NROW(daily.means)

summary_Fd.Q["Date_DOY"] <- date(summary_Fd.Q$Date)

daily.sums <- aggregate(cbind(ETo_mm, All_Q, PPFD_IN_1_2_1.sm = PPFD_IN_1_2_1, VPD_kPa.sm =VPD_kPa, P_1_1_1,
                              Pines_Q, Pine_15DIAM_Q,Pine_25DIAM_Q,Pine_.25DIAM_Q,
                              Spruce_Q, Spruce_15DIAM_Q, Spruce_25DIAM_Q, Spruce_.25DIAM_Q, 
                              Birch_Q, Birch_15DIAM_Q, Birch_25DIAM_Q)~Date_DOY, data=summary_Fd.Q, FUN = sum, na.action = na.omit)


# Change name of variable for "merge"
colnames(daily.sums)[1]<-c("Date")

#daily.sums["Date"] <- date(as.POSIXct((daily.sums$Date/48)-42300, origin="1970/01/1", tz="UTC"))

# Add daylength-corrected VPD
names(summary_Fd.Q)
summary_Fd.Q_daylength <-  subset(summary_Fd.Q, summary_Fd.Q$PPFD_IN_1_2_1>20, select=c("Date","PPFD_IN_1_2_1", "VPD_kPa", "DOY_cont",  "Month", "Year" ),na.action=na.omit)

summary_Fd.Q_daylength["Date"] <- date(as.POSIXct(summary_Fd.Q_daylength$Date, origin="1970/01/1", tz="UTC"))

temp <- aggregate(cbind(DOY_cont, PPFD_IN_1_2_1, VPD_kPa)~Date, data=summary_Fd.Q_daylength, FUN=function(x) c(mean = mean(x), lgt = length(x)), na.action=na.omit)

# Estimate daylength-corrected VPD: Dz
temp["Dz"] <- temp$VPD_kPa[,1]*(temp$VPD_kPa[,2]/48)

# Estimate daylength-corrected PAR: PARz

temp["PARz"] <- temp$PPFD_IN_1_2_1[,1]*(temp$PPFD_IN_1_2_1[,2]/48)

as.numeric(temp$Date)

as.numeric(daily.data$Date)

as.numeric(daily.means$Date)
#The part below is no longer needed

#temp["Date"] <- date(as.POSIXct(temp[["Date"]][,1], origin="1970/01/1", tz="UTC"))

#temp["DOY_cont"]  <- temp[["DOY_cont"]][,1]

#daily.data["DOY"] <- lubridate::yday(daily.data$Date)
#daily.data["Year"] <- lubridate::year(daily.data$Date)

temp.date <- date(seq(as.POSIXct("2016/01/01", tz="UTC"),as.POSIXct("2019/12/31", tz="UTC"), by = "day"))

daily.data <- setNames(as.data.frame(matrix(ncol = 1, nrow = NROW(temp.date))),c("Date"))

daily.data["Date"] <- temp.date
rm(temp.date)

daily.data <- merge(daily.data, daily.sums, by="Date", all.x=T, all.y=F) 

#colnames(daily.data)[3]<-c("DOY_cont")

#daily.data <- merge(daily.data, daily.means, by="Date", all.x=T, all.y=F) 


daily.data["Dz"] <- merge(daily.data, temp, by="Date", all.x=T, all.y=F)$Dz

daily.data["PARz"] <- merge(daily.data, temp, by="Date", all.x=T, all.y=F)$PARz

rm(daily.sums, daily.means)

#I think the stuff below is not needed
names(daily.data)

# Plot test to see if things work 

plot(daily.data$Dz, daily.data$Pines_Q)

rm(temp, temp.2, temp.3, temp.4)

plot(daily.data$Date.x,daily.data$S0.15)



# ADD: Long-term drought data to daily.data  ----------------------------------------------------------------

names(daily.data)
daily.data["Year"] <- lubridate::year(daily.data$Date)
daily.data["Month"] <- lubridate::month(daily.data$Date)
daily.data["DOY"] <- lubridate::yday(daily.data$Date)

daily.data["Week"] <- ifelse(daily.data$Year==2016,lubridate::week(daily.data$Date),
                             ifelse((daily.data$Year==2017),lubridate::week(daily.data$Date)+53,
                                    ifelse((daily.data$Year==2018),lubridate::week(daily.data$Date)+106,
                                           lubridate::week(daily.data$Date)+159)))


longterm_drought$Date
temp.drought <- subset(longterm_drought, longterm_drought$Year>2015)

temp.drought["Week"] <- ifelse(temp.drought$Year==2016,lubridate::week(temp.drought$Date),
                               ifelse((temp.drought$Year==2017),lubridate::week(temp.drought$Date)+53,
                                      ifelse((temp.drought$Year==2018),lubridate::week(temp.drought$Date)+106,
                                             lubridate::week(temp.drought$Date)+159)))


as.Date(temp.drought$Date)
daily.data["spei"] <- merge(daily.data, temp.drought, by="Week", all.x= T, all.y=T)["spei"]

rm(temp.drought)

plot(daily.data$Date, daily.data$spei)

write.table(daily.data, "Data/daily.data.csv", sep=",", row.names=F)
# END OF SECTION ----

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#
#  LONG-TERM ENVIRONMENTAL DATA
#
000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

# OPEN - LONG_TERM DATA  ----------------------------------------------------------------

# Open files
histo_data = read.table(file="Data/histo_data.csv",header=T,sep=",")

names(histo_data)

histo_data["Date"] <- as_datetime(seq(as.POSIXct("1985/01/1 00:00:00", tz="UTC"), as.POSIXct("2020/12/30 23:30:00", tz="UTC"), by = "day"))


# Add Year, Month, Week and DOY 
histo_data["Year"] <- lubridate::year(histo_data$Date)
histo_data["Month"] <- lubridate::month(histo_data$Date)
histo_data["Week"] <- lubridate::week(histo_data$Date)
histo_data["DOY"] <- lubridate::yday(histo_data$Date)

# ADD - ETo TO LONG_TERM DATA  ----------------------------------------------------------------

names(all.eco.met)
temp.eto <- aggregate(cbind(Date, ETo_mm, Ta_1_1_1, RH_1_1_1)~DOY_cont, data=all.eco.met, FUN=sum, na.action = na.omit)
temp.eto["Date"] <- aggregate(cbind(Date, ETo_mm)~DOY_cont, data=all.eco.met, FUN=min, na.action = na.omit)["Date"]
head(temp.eto$Date, 5)

temp.eto$Date <- as.POSIXct(temp.eto$Date, origin="1970/01/1", tz="UTC")

histo_data <- merge(histo_data, temp.eto, by="Date", all.x=T)
#par(mar=c(1,1,1,1))
plot(histo_data$ETo_mm)

tail(histo_data$ETo_mm, 50)

rm(temp.eto)

# ADD - MISSING PP TO histo_data  ----------------------------------------------------------------

#histo_data$precip[is.na(histo_data$precip)] <- 0

histo_data$precip[is.na(histo_data$precip)&histo_data$Year>2015] 

histo_data$Date  <- lubridate::date(histo_data$Date)

histo_data["precip.16.20"] <- NA

daily.data$Date
histo_data$precip.16.19  <- merge(histo_data, daily.data, by="Date", all.x=T, all.y = F)["P_1_1_1.no.mis"]

tail(histo_data$precip.16.19, 300)

str(histo_data)

histo_data["precip.nomiss"] <-  NA

histo_data$precip.nomiss <- unlist(ifelse(histo_data$precip<0,histo_data$precip.16.19, histo_data$precip))



# COUNT - DAYS WITHOUT PP  ----------------------------------------------------------------

# All pp under 0.1 = 0

histo_data["precip.adj"] <- NA


histo_data$precip.adj <- ifelse(histo_data$precip<0.1&is.na(histo_data$precip), 0, histo_data$precip)


histo_data$precip[is.na(histo_data$precip)]
names(histo_data)

head(histo_data, 30)


### Identify dates where no precipitatio occurs continuously 

histo_data["pp.cont"] <- NA
i = 1
j = NA
y=NA
for(i in 1:NROW(histo_data)){
  
  j <- histo_data$precip.nomiss[i]

  if(j>0){y=0}
  
  if(j==0){y=y+1
  }else{y=y}
  
  print(y)
  histo_data$pp.cont[i] <- y
  
  i=i+1
  
}

rm(i, j, y)



histo_data["pp.cont.lngt"] <- NA

i = 1
j = NA
y=NA
for(i in 1:NROW(histo_data)){
  
  #j  <- histo_data$pp.cont[i-1]
  xx <- histo_data$pp.cont[i]
  
  y <- ifelse(xx==0, histo_data$pp.cont[i-1], NA)
  
  histo_data$pp.cont.lngt[i-1] <- y
  print(y)
  # dt.temp$pp.fix[i] <- y
  
  i=i+1
  
}

hist.temp <- subset(histo_data, histo_data$pp.cont.lngt>0&histo_data$Month>3&histo_data$Month<11, select=(pp.cont.lngt))


length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==1])/length(hist.temp$pp.cont.lngt)

histogram <- setNames(as.data.frame(matrix(ncol = 3, nrow = 30)),c("ID", "n", "percent" ))

histogram$ID <- seq(1,30)

histogram$n[1] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==1])
histogram$n[2] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==2])
histogram$n[3] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==3])
histogram$n[4] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==4])
histogram$n[5] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==5])
histogram$n[6] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==6])
histogram$n[7] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==7])
histogram$n[8] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==8])
histogram$n[9] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==9])
histogram$n[10] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==10])
histogram$n[11] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==11])
histogram$n[12] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==22])
histogram$n[13] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==23])
histogram$n[14] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==24])
histogram$n[15] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==25])
histogram$n[16] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==26])
histogram$n[17] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==27])
histogram$n[18] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==28])
histogram$n[19] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==29])
histogram$n[20] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==20])
histogram$n[27] <- length(hist.temp$pp.cont.lngt[hist.temp$pp.cont.lngt==27])

histogram$percent <-  histogram$n/length(hist.temp$pp.cont.lngt)

is.na(hist.temp$pp.cont.lngt)


plot(NULL, xlim = c(1,12), ylim=c(0,0.5),
     ylab="Percent",
     xlab="Length of period (days)")

arrows(histogram$ID, 0, histogram$ID, histogram$percent, length=0, angle=50, code=2, lwd= 15, col="blue")

text(histogram$ID, histogram$percent+0.05, round(histogram$percent, 2),col="deeppink", cex=1)
text(histogram$ID, histogram$percent+0.02, paste("n=", histogram$n,sep=""),col="black", cex=0.8)




# End of this script
names(histo_data)

# This also works great 
#no.pp.count <- rle(histo_data$precip.adj)$lengths[rle(histo_data$precip.adj)$values==0]

hist(no.pp.count)
rm(no.pp.count)
rm(i, j, y, xx, hist.temp)
rm(no.pp.count)
rm(histogram)

#write.table(no.pp.count, "no.pp.count.csv", sep=",", row.names = F)





# SAVE - HISTORIC DATA ----------



# OPEN - HISTORIC DATA ----------


# PREDICT - LONG_TERM ETo  ----------------------------------------------------------------
# Predictions are made with VPD and AirT. VPD is used primarily, and AirT is used when VPD data is missing

summary(histo_data$VPD_data.range)
x = unlist(subset(histo_data, histo_data$VPD_data.range=="Good" , select=c("DATE", "VPD"))[2])
y = unlist(subset(histo_data, histo_data$VPD_data.range=="Good" , select=c("DATE", "ETo_mm"))[2])

# Send model results to variable

plot(x, y)

ETo_model <-suppressWarnings( nlsLM(y ~ b0*(1 - (b1 * exp(-b2 * x ))),
                                    data=histo_data, start=list(b0=1, b1=1, b2=1)))

# Store predicted coefficients for legend
b0_ETo <- summary(ETo_model)$coefficients[1, 1]
b1_ETo <- summary(ETo_model)$coefficients[2, 1]
b2_ETo <- summary(ETo_model)$coefficients[3, 1]

# Created column with predicted ETo
pred.eto <- b0_ETo*(1 - b1_ETo * exp(-b2_ETo *histo_data$VPD)) 
histo_data["ETo_mm_pred.VPD"] <- ifelse(pred.eto<0, 0, pred.eto)
#summary_DOY["ETo_mm_pred"] <- 2.32924126330696 / (1 + 0.337335339562585 * exp( -0.121242767960978 * summary_DOY$Ta_1_1_1 )) ^ (1 / 0.0807243045601035)

plot(histo_data$ETo_mm, ylim=c(0, 6))
lines(histo_data$ETo_mm_pred.VPD, col="red")

rm(x, y,ETo_model, b0_ETo, b1_ETo, b2_ETo, pred.eto)

# Predict ETo using AirT

summary(histo_data$VPD_data.range)

x = unlist(subset(histo_data, histo_data$VPD_data.range=="Good"&histo_data$AirT>0 , select=c("DATE", "AirT"))[2])
y = unlist(subset(histo_data, histo_data$VPD_data.range=="Good"&histo_data$AirT>0 , select=c("DATE", "ETo_mm"))[2])

# Send model results to variable

plot(x, y)

ETo_model <-suppressWarnings( nlsLM(y ~ b0*(1 - (b1 * exp(-b2 * x ))),
                                    data=histo_data, start=list(b0=-2, b1=0.99, b2=0.02)))

# Store predicted coefficients for legend
b0_ETo <- summary(ETo_model)$coefficients[1, 1]
b1_ETo <- summary(ETo_model)$coefficients[2, 1]
b2_ETo <- summary(ETo_model)$coefficients[3, 1]


# Created column with predicted ETo
pred.eto <- b0_ETo*(1 - b1_ETo * exp(-b2_ETo *histo_data$AirT)) 
histo_data["ETo_mm_pred_AirT"] <- ifelse(pred.eto<0, 0, pred.eto)
#summary_DOY["ETo_mm_pred"] <- 2.32924126330696 / (1 + 0.337335339562585 * exp( -0.121242767960978 * summary_DOY$Ta_1_1_1 )) ^ (1 / 0.0807243045601035)
rm(pred.eto)
plot(histo_data$ETo_mm, ylim=c(0, 6))
lines(histo_data$ETo_mm_pred_AirT, col="red")

rm(x, y,ETo_model, b0_ETo, b1_ETo, b2_ETo)

# Estimate a continuous ETo

histo_data["ETo_mm.combined"] <- ifelse(is.na(histo_data$ETo_mm_pred.VPD), histo_data$ETo_mm_pred_AirT, histo_data$ETo_mm_pred.VPD)

plot(histo_data$ETo_mm, ylim=c(0, 6))
lines(histo_data$ETo_mm.combined, col="red")


rm(b0_ETo, b1_ETo, b2_ETo, x, y, ETo_model)


# END OF SECTION ----

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#
#  LONG-TERM DROUGHT
#
000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

# ADD - LONG_TERM DROUGHT INDEX  ----------------------------------------------------------------

# Sum by week
histo_data["ID"] <- NA
histo_data$ID <- 1:NROW(histo_data)

histo_data$ID <- as.numeric(histo_data$ID)

plot(histo_data$Date, histo_data$ETo_mm.combined)
head(histo_data$Date, 5)
names(histo_data)

str(histo_data)
longterm_drought <- aggregate(cbind(precip.nomiss, ETo_mm.combined, Date, ID, DATE)~Week+Year, data=histo_data,FUN=sum, na.rm=T, na.action=NULL)

temp.date <- aggregate(cbind(Date)~Week+Year, data=histo_data,FUN=min, na.rm=F, na.action=NULL)


temp.date["Date"]  <-  histo_data[!duplicated(histo_data[c("Year", "Week")]),]["Date"]

longterm_drought$Date <- temp.date$Date

longterm_drought$Date <- as.Date(longterm_drought$Date)

rm(temp.date)

longterm_drought.month <- aggregate(cbind(precip.nomiss, ETo_mm.combined, Date, ID)~Month+Year, data=histo_data,FUN=sum, na.rm=T, na.action=NULL)
date.month <- aggregate(cbind(precip.nomiss, ETo_mm.combined, Date, ID)~Month+Year, data=histo_data,FUN=min, na.rm=T, na.action=NULL)

longterm_drought.month$ID <- aggregate(cbind(ID)~Month+Year, data=histo_data,FUN=min, na.rm=T, na.action=NULL)["ID"]

date.month$Date <- as.POSIXct(date.month$Date, origin="1970/01/1", tz="UTC")
date.month$Date <- as.Date(date.month$Date)

longterm_drought.month$Date <-date.month$Date
rm(date.month)

longterm_drought["PP_ET"] <- longterm_drought$precip-longterm_drought$ETo_mm.combined
longterm_drought.month["PP_ET.month"] <- longterm_drought.month$precip-longterm_drought.month$ETo_mm.combined

rm(temp_date)
# Estimate SPEI week
spei_out <- spei(longterm_drought$PP_ET, 1, kernel = list(type = 'rectangular', shift = 0),
                 distribution = 'log-Logistic', fit="ub-pwm", na.rm = FALSE,
                 ref.start=NULL, ref.end=NULL, x=FALSE, params=NULL)

plot(spei_out)

# Save to dataframe
longterm_drought["spei"] <- as.data.frame(spei_out[[2]])


# Estimate SPEI month
spei_out <- spei(longterm_drought.month$PP_ET, 1, kernel = list(type = 'rectangular', shift = 0),
                 distribution = 'log-Logistic', fit="ub-pwm", na.rm = FALSE,
                 ref.start=NULL, ref.end=NULL, x=FALSE, params=NULL)


#par(mar = c(10, 10, 10, 10))
#par(mar=c(1,1,1,1))
plot(spei_out)
par("mar")

# Save to dataframe
longterm_drought.month["spei"] <- as.data.frame(spei_out[[2]])

longterm_drought["spei.month"] <- merge(longterm_drought, longterm_drought.month, by="ID", all.x = T)["spei.y"]

#lubridate::date(longterm_drought$Date)

plot(longterm_drought$ID, longterm_drought$spei, type="l", col="blue")
plot(longterm_drought.month$spei, col="red")

rm(spei_out, longterm_drought.month)

# Decided not to estimate the other INDICES

#scPDSI_out <- pdsi(longterm_drought$pp, longterm_drought$ET0_mod.rich_gwth.year, sc=TRUE, AWC=50 )
#longterm_drought["scPDSI"] <- as.data.frame(scPDSI_out[[2]])[1:1749,]


# Replace week 1 and 53 with 0, to close the polygon
# This is done for esthetic purposes. Do not delete if the plot is not a polygon closed every year
longterm_drought$spei <- ifelse(longterm_drought$Week==1, 0, longterm_drought$spei)
longterm_drought$spei <- ifelse(longterm_drought$Week==53, 0, longterm_drought$spei)

longterm_drought.month$spei <- ifelse(longterm_drought.month$Month==1, 0, longterm_drought.month$spei)
longterm_drought.month$spei <- ifelse(longterm_drought.month$Month==12, 0, longterm_drought.month$spei)


# SAVE - LONG_TERM DROUGHT INDEX  ----------------------------------------------------------------
names(longterm_drought)

write.table(longterm_drought, "Data/longterm_drought.csv", sep=",", row.names=FALSE)


# OPEN - LONG_TERM DROUGHT INDEX  ----------------------------------------------------------------

longterm_drought <- read.table(file="Data/longterm_drought.csv",header=T,sep=",")
names(longterm_drought)

# END OF SECTION ----

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#
#  PARAMETERS 
#
000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

# This includes a database with parameters fitted for a various models/equations
# They are compiled here in a single file, to simplify the data analysis
# OPEN - PARAMETER DATABASE ----------------------------------------------------------------

par.estimates <- read.table(file="Data/par.estimates.csv",header=T,sep=",")

# END OF SECTION ----
000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#
#  TREE DATA
#
000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

# OPEN - DATA FOR FIT ALLOMETRIC EQUATIONS [PINE, SPRUCE, BIRCH]  ----------------------------------------------------------------

# Load databases

dt.allomet <- read.table(file="Data/allometric_data.dat",header=T,sep=",", na.strings = "NAN" )

dt.allomet <- subset(dt.allomet, dt.allomet$Outlier=="No")

# Ignore birch data

#dt.allomet_birch = read.table(file="/Users/data/krycklan_birch.csv",header=T,sep=",", na.strings = "NAN" )



colnames(dt.allomet)

# CREATE DATAFRAMES  ----------------------------------------------------------------

dt.allomet.df <-data.frame(ID=dt.allomet$ID, species=dt.allomet$species, status=dt.allomet$status,diamclass=dt.allomet$diamclass_cm,
                           sapdepth = dt.allomet$R_cm, corelength=dt.allomet$corelength_cm, complete_r=dt.allomet$Complete_R,
                           dbh_cm=dt.allomet$DBH_cm)

#dt.allomet_RA.df <- as.data.frame(dt.allomet_RA)

#rm(dt.allomet_pine)

colnames(dt.allomet.df)


# ESTIMATE SAPWOOD AREA FOR ALL SPECIES  ----------------------------------------------------------------

# Estimate sapwood area 

dt.allomet.df["ba_m"] <- (0.7854 * (dt.allomet.df$dbh_cm)^2)/10000

plot(dt.allomet.df$dbh_cm ,dt.allomet.df$sapdepth)

# Force-correct these values. Likely outliers
dt.allomet.df$bark.depth.cm[21] <-3
dt.allomet.df$bark.depth.cm[10] <-1

dt.allomet.df["As.core_m"] <- ((0.7854*dt.allomet.df$corelength^2) - (0.7854*(dt.allomet.df$corelength-(dt.allomet.df$sapdepth*2))^2))/10000

dt.allomet.df["Bark.depth_cm"] <- dt.allomet.df$dbh_cm - dt.allomet.df$corelength

dt.allomet.df["Ab.dbh_m"] <- (0.7854*dt.allomet.df$dbh_cm^2) - (0.7854*((dt.allomet.df$dbh_cm)-(dt.allomet.df$Bark.depth_cm*2))^2)

dt.allomet.df["As.dbh_m"] <- (dt.allomet.df$Ab.dbh_m - (0.7854*(dt.allomet.df$dbh_cm-dt.allomet.df$Bark.depth_cm*2-dt.allomet.df$sapdepth*2)^2))/10000

# Plot to test results
plot(dt.allomet.df$dbh_cm,dt.allomet.df$As.core_m, col=ifelse(dt.allomet.df$species=="Pine", "red", "blue"))

# Plot to test results
plot(dt.allomet.df$dbh_cm, dt.allomet.df$As.dbh_m, col=c(ifelse(dt.allomet.df$species=="Pine", "red", "blue")))

# Next, remove sapwood area estimate, that looks way off
dt.allomet.df$As.dbh_m[21] <- NA

# Below, just to test if the formulas equal cero, and they do
#dt.allomet.df["test_total"]  <- ((0.7854*dt.allomet.df$ave_dbh^2)/10000) - ((0.7854*(dt.allomet.df$ave_dbh-(dt.allomet.df$bark.depth.cm*2))^2)/10000) - ((0.7854*(dt.allomet.df$ave_dbh-(dt.allomet.df$bark.depth.cm*2)-(dt.allomet.df$ave_sapdepth*2))^2)/10000) - dt.allomet.df$As.dbh_m

# Estimate sapwood area for Birch 
# Birch not included in this version. Ignore

dt.allomet_birch["diam_cm"] <- dt.allomet_birch$Circumference_cm/3.1416

A = dt.allomet_birch$Sapwood_cm * 2
B = dt.allomet_birch$Bark_cm * 2
Diam = dt.allomet_birch$diam_cm


dt.allomet_birch["saparea_m"] <-   (3.1415 * (Diam / 2) ^ 2 - 3.1415 * ((Diam - A - B) / 2) ^ 2)/10000

rm(A, B, Diam)

# OPEN - TREE DATA & ESTIMATE As & Ds ----------------------------------------------------------------
# These are the actual trees used for sap flow measurements

# Import diameter for each tree

tree.info = read.table(file="Data/tree.diam.info.csv",header=T,sep=",")

# Dataframe with unique name IDS (needed later on)
temp.names <- as.data.frame(paste(paste("Fd", tree.info$Unique.name.ID[grep("Pine.*_02CM_.*|Spruce.*_02CM_.*", tree.info$Unique.name.ID)], sep = "_"), "radial", sep="."))

tree.dat.temp <- setNames(as.data.frame(matrix(ncol = 3, nrow = NROW(temp.names))),c("TreeID", "DBH", "Ds"))

tree.dat.temp["TreeID"] <- temp.names

rm(temp.names)

tree.dat.temp["DBH"] <- as.data.frame(tree.info$Diameter[grep("Pine.*_02CM_.*|Spruce.*_02CM_.*", tree.info$Unique.name.ID)])

# Add species

tree.dat.temp["Species"] <- NA

tree.dat.temp$Species[grep("Pine", tree.dat.temp$TreeID)] <- "Pine"

tree.dat.temp$Species[grep("Spruce", tree.dat.temp$TreeID)] <- "Spruce"

# Add sapwood area 

tree.dat.temp["As"] <- NA

tree.dat.temp$As[grep("Pine", tree.dat.temp$TreeID)] <- exp(par.estimates[7,2] + par.estimates[7,3] * log(tree.dat.temp$DBH[grep("Pine", tree.dat.temp$TreeID)]))*10000

tree.dat.temp$As[grep("Spruce", tree.dat.temp$TreeID)] <- exp(par.estimates[8,2] + par.estimates[8,3] * log(tree.dat.temp$DBH[grep("Spruce", tree.dat.temp$TreeID)]))*10000

plot(tree.dat.temp$DBH, tree.dat.temp$As)


# Add sapwood depth 

tree.dat.temp["Ds"] <- NA

as_pine   <- exp(par.estimates[5,2] + par.estimates[5,3] * log(tree.dat.temp$DBH))
as_spruce <- exp(par.estimates[6,2] + par.estimates[6,3] * log(tree.dat.temp$DBH))


tree.dat.temp$Ds[tree.dat.temp$Species=="Pine"] <- exp(par.estimates[5,2] + par.estimates[5,3] * log(tree.dat.temp$DBH[tree.dat.temp$Species=="Pine"]))

tree.dat.temp$Ds[tree.dat.temp$Species=="Spruce"] <- exp(par.estimates[6,2] + par.estimates[6,3] * log(tree.dat.temp$DBH[tree.dat.temp$Species=="Spruce"]))


rm(as_pine, as_spruce)

plot(tree.dat.temp$DBH, tree.dat.temp$Ds, col=ifelse(tree.dat.temp$Species=="Pine", "Red", "blue"))

plot(tree.dat.temp$DBH, tree.dat.temp$As, col=ifelse(tree.dat.temp$Species=="Pine", "Red", "blue"))

# END OF SECTION ----

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#
#  SAP FLOW DATA
#
000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

# IMPORT - RAW SAP FLOW DATA ----------------------------------------------------------------

# Read data table

all_nodes_raw = read.table(file="AllNodes_Umea_RAW_alligned_March_01_2021.dat",header=T,sep=",", na.strings = "NAN" )

# Turn it into a data frame

all_nodes_raw.DF <- as.data.frame(all_nodes_raw)

# REMOVE - BAD SENSORS ----------------------------------------------------------------


# Remove columns that have bad data 

env_initial <- new.env()

env_initial$badsensors <- c("NodeONE_28Spruce20_02CM_Random_7.15DIAM_",
                            "NodeTWO_14Spruce19_46CM_Random_15.25DIAM_",
                            "NodeTWO_20Spruce3_46CM_Random_15.25DIAM_",
                            "NodeTWO_28Pine13_02CM_Random_7.15DIAM_",
                            "NodeTWO_50Spruce6_24CM_Random_15.25DIAM_",
                            "NodeTHREE_15Pine7_02CM_Random_15.25DIAM_", 
                            "NodeTHREE_39Pine20_68CM_Random_25.30DIAM_",
                            "NodeTHREE_43Pine5_02CM_Random_25.30DIAM_",
                            "NodeONE_28_Status",
                            "NodeTWO_14_Status",
                            "NodeTWO_20_Status",
                            "NodeTWO_28_Status",
                            "NodeTWO_50_Status",
                            "NodeTHREE_15_Status",
                            "NodeTHREE_39_Status",
                            "NodeTHREE_43_Status")

# Sensor 26 is good, I just made a mistake and label it all bad data
# April 30, 2019. Sensor 26 removed from "Bad sensors"
# June 16, 2019. Sensor 28 added to bad sensors; too much noise
# June 16, 2019. Sensor Node Three 43 added to bad sensors; too much noise
# August 2, 2019. Decide what to do with NodeTwo_35, seems noisy
# August 2, 2019. NodeTwo coded only until July/7
# August 3, 2019, Sensor NodeThree_10 has noise. Process, and decide if the noise stays
# August 3, 2019. Sensor NodeThree_6. Heater broken after May 28, 2019. Fix?
# August 6, 2019. Sensor NodeTwo_50. Added to bad sensors; too too little good data, and too much noisy data
# August 6, 2019. Sensor NodeOne_26. Added to bad sensors; too little good data

# May 31, 2020. Sensor NodeONE_26 was reviewed. 
# May 31, 2020. The rest sensors were confirmed to have unnacceptable amounts of noise 


all_nodes_raw.DF[ ,c(env_initial$badsensors)]  <- list(NULL)

rm(all_nodes_raw)

names(all_nodes_raw.DF)


# FILTER - REMOVE NOISE DUE TO TEMPERATURE GRADIENTS  ----------------------------------------------------------------

# Perform initial filtering based on temperature-induced noise
# Noise was previously marked in another step
# "Good" data is will be marked with the number 6 [or any other number]

# Location for first and last columns with raw mV

min_mV = min(grep("Pine|Spruce|Birch", names(all_nodes_raw.DF)))
max_mV = max(grep("Pine|Spruce|Birch", names(all_nodes_raw.DF)))

min_tmax = min(grep("Status", names(all_nodes_raw.DF)))
max_tmax = max(grep("Status", names(all_nodes_raw.DF)))

# This loop tags noisy data with a random number (6)

for (i in min_tmax:max_tmax) {
  isna <- is.na(all_nodes_raw.DF[i])
  all_nodes_raw.DF[i][isna] <- 6 # Random number, just to ID good data
  i=i+1
}

all_nodes_raw.DF["ID"] <- seq(1, NROW(all_nodes_raw.DF))


# Add Timestamp to raw data base

all_nodes_raw.DF["Date"] <- (Date= as_datetime(seq(as.POSIXct("2016/06/17 00:00:00", tz="UTC"), as.POSIXct("2020/12/30 23:30:00", tz="UTC"), by = "30 min")))

#head(all_nodes_raw.DF$TIMESTAMP_30min, 5)

rm(min_mV,max_mV,min_tmax, max_tmax, i, isna, Date)

# The section below is to split the entire data base. Currently not being used

# Make a list of data that contains all Birch sensors 

#env_initial$list.birch <- grep("Birch", names(all_nodes_raw.DF))

# Subset all data for birch, sensors and status columns 

#all_nodes_raw.birch.DF <-  subset(all_nodes_raw.DF, select=c(1:106, env_initial$list.birch, 325:334, 437:444))

#write.table(all_nodes_raw.birch.DF, "all_nodes_raw.birch.DF.csv", sep=",", row.names=FALSE)

# Next, remove everything associated with Birch

#all_nodes_raw.DF[c(env_initial$list.birch, 325:334, 437:444)] = NULL

colnames(all_nodes_raw.DF)





# IN CASE I NEED TO DELETE TMAX ESTIMATES   ----------------------------------------------------------------
# Not the actual parameters, needs to be edited: June 30, 2019

#all_nodes_raw.DF_aligned <- all_nodes_raw.DF_aligned[, -c(min(grep("^TMax", names(all_nodes_raw.DF_aligned))):max(grep("^TMax", names(all_nodes_raw.DF_aligned))))] 

#dt_irrig <- dt_irrig[, -c(132:312)] 

#names(all_nodes_raw.DF_aligned)

# ADD - MET & ECO DATA TO RAW mV DATABASE -> all_nodes_raw.DF  ----------------------------------------------------------------

names(all_nodes_raw.DF)
names(all.eco.met)

all_nodes_raw.DF["Ta_1_1_1"] <- merge(all_nodes_raw.DF, all.eco.met, by="Date", all.y=F)["Ta_1_1_1.y"]
all_nodes_raw.DF["RH_1_1_1"] <- merge(all_nodes_raw.DF, all.eco.met, by="Date", all.y=F)["RH_1_1_1.y"]
all_nodes_raw.DF["VPD_kPa"] <- merge(all_nodes_raw.DF, all.eco.met, by="Date", all.y=F)["VPD_kPa.y"]
all_nodes_raw.DF["ETo_mm"] <- merge(all_nodes_raw.DF, all.eco.met, by="Date", all.y=F)["ETo_mm.y"]

plot(aggregate(cbind(ETo_mm)~DOY, data=all_nodes_raw.DF, FUN=sum))


# DATA MANAGEMENT PRIOR TO T.MAX DETERMINATION  ----------------------------------------------------------------


# Identify again the columns with Tmax status

min_tmax = min(grep("Status", names(all_nodes_raw.DF)))
max_tmax = max(grep("Status", names(all_nodes_raw.DF)))

# Identify columns with raw mV data 

min_mV = min(grep("Pine|Spruce|Birch", names(all_nodes_raw.DF)))
max_mV = max(grep("Pine|Spruce|Birch", names(all_nodes_raw.DF)))

# Extract names of columns with data and names of columns with data quality information
env_initial$listofnames <- colnames(as.data.frame(all_nodes_raw.DF[min_mV:max_mV]))
env_initial$listofstatus <- colnames(as.data.frame(all_nodes_raw.DF[min_tmax:max_tmax]))

# Create artificial data frame (will be used later for comparison purposes)
dataframe <- (Date= as_datetime(seq(as.POSIXct("2016/06/17 00:00:00", tz="UTC"), as.POSIXct("2020/12/30 23:30:00", tz="UTC"), by = "30 min")))
NROW(all_nodes_raw.DF)
NROW(dataframe)
#rm(dtNODE_ONE_DF_aligned)

# Empty data frame that will be populated in the "For Loop" below
# Needs names of colums for later "merge"

number_of_cols <- length(colnames(all_nodes_raw.DF[min_mV:max_mV]))+7

all_nodes_raw.DF_aligned <- setNames(as.data.frame(matrix(ncol = number_of_cols, nrow = NROW(dataframe))),c("Date", "VPD", "VPD_0.05", "PAR", "PAR_50", "RH", "AirT", colnames(all_nodes_raw.DF[min_mV:max_mV])))

rm(number_of_cols)
# Before aligning these two, I MUST make sure that both datastamps are the same, else, all will be wrong

# Insert original Timestamp in artificial dataframe
all_nodes_raw.DF_aligned$Date <- dataframe

# Insert ID variable to original and artificial dataframes 
all_nodes_raw.DF["ID"] <- seq(1:NROW(all_nodes_raw.DF))
all_nodes_raw.DF_aligned["ID"] <- seq(1:NROW(all_nodes_raw.DF_aligned))

#Insert PAR and VPD to aligned dataframe 
all_nodes_raw.DF_aligned$VPD <- all.eco.met$VPD_kPa[8065:NROW(all.eco.met)]
all_nodes_raw.DF_aligned$PAR <- all.eco.met$PPFD_IN_1_2_1[8065:NROW(all.eco.met)]
#all_nodes_raw.DF_aligned$VPD_0.05 <- all.eco.met$VPD_0.2[8065:NROW(all.eco.met)]
#all_nodes_raw.DF_aligned$PAR_50 <- all.eco.met$PAR_50[8065:NROW(all.eco.met)]
all_nodes_raw.DF_aligned$AirT <- all.eco.met$Ta_1_1_1[8065:NROW(all.eco.met)]


# Insert DOY and Year for original and artificial data frame, using the actual date [will be used later for quality control]
all_nodes_raw.DF_aligned["DOEY"] <- lubridate::yday(dataframe)
all_nodes_raw.DF_aligned["Year"] <- lubridate::year(dataframe)

all_nodes_raw.DF["DOEY"] <- lubridate::yday(dataframe)
all_nodes_raw.DF["Year"] <- lubridate::year(dataframe)

# Originally, I was going to split by year, but then decided agains it; it creates more problems for the script [IDEA DISCARDED]

# Remove incomplete days from the data frame, both original and fake, so that we can add time of day (TOD) and DOY

#all_nodes_raw.DF <- all_nodes_raw.DF[c(1:51600),]

#all_nodes_raw.DF_aligned <- all_nodes_raw.DF_aligned[c(1:51600),]

# Rewrite the Artificial dataframe to march the date on the last column:

#dataframe <- (Date= as_datetime(seq(as.POSIXct("2016/06/17 00:00:00", tz="UTC"), as.POSIXct("2019/06/12 23:30:00", tz="UTC"), by = "30 min")))

# Insert TOD (time of day)
all_nodes_raw.DF["TOD"] <- rep(seq(0,1410, 30), times=NROW(all_nodes_raw.DF)/48, each=1)
names(all_nodes_raw.DF)

all_nodes_raw.DF_aligned["TOD"] <- rep(seq(0,1410, 30), times=NROW(all_nodes_raw.DF_aligned)/48, each=1)


# Insert an artificial DOY (DOY_cont), that will be used to aggregate all the database without having to worry about similar DOY across years

all_nodes_raw.DF_aligned["DOY"] <- ifelse(all_nodes_raw.DF_aligned$Year==2016,lubridate::yday(all_nodes_raw.DF_aligned$Date),
                                          ifelse((all_nodes_raw.DF_aligned$Year==2017),lubridate::yday(all_nodes_raw.DF_aligned$Date)+366,
                                                 ifelse((all_nodes_raw.DF_aligned$Year==2018),lubridate::yday(all_nodes_raw.DF_aligned$Date)+731,
                                                        ifelse((all_nodes_raw.DF_aligned$Year==2019),lubridate::yday(all_nodes_raw.DF_aligned$Date)+1096,
                                                               lubridate::yday(all_nodes_raw.DF_aligned$Date)+1461))))

all_nodes_raw.DF["DOY"]         <- ifelse(all_nodes_raw.DF$Year==2016,lubridate::yday(all_nodes_raw.DF$Date),
                                          ifelse((all_nodes_raw.DF$Year==2017),lubridate::yday(all_nodes_raw.DF$Date)+366,
                                                 ifelse((all_nodes_raw.DF$Year==2018),lubridate::yday(all_nodes_raw.DF$Date)+731,
                                                        ifelse((all_nodes_raw.DF$Year==2019),lubridate::yday(all_nodes_raw.DF$Date)+1096,
                                                               lubridate::yday(all_nodes_raw.DF$Date)+1461))))

NROW(all_nodes_raw.DF_aligned)
NROW(all_nodes_raw.DF)

#colnames(all_nodes_raw.DF_aligned)
#colnames(all_nodes_raw.DF)

#plot(all_nodes_raw.DF$ID, all_nodes_raw.DF$NodeONE_4Spruce14_02CM_Random_7.15DIAM_, ylim=c(0,2),xlim=c(1,4000))

sequence_DOY <- rep(seq(1,NROW(all_nodes_raw.DF_aligned)/48, 1), times=1, each=48)

# Temporary data frame, to merge Tmax 
temp_dataframe <- data.frame(DOY=sequence_DOY)

# Add ID to previous file

temp_dataframe["ID"] <- seq(1,NROW(sequence_DOY))

#colnames(temp_dataframe)

# Names for the Tmax values

Tmaxset <- paste("Tmax", colnames(all_nodes_raw.DF[min_mV:max_mV]), sep = "_")

#rm(Date, dataframe, min_mV, max_mV, min_tmax, max_tmax, sequence_DOY)




# END OF SECTION ----

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
# ESTIMATE - SLOPE OF RAW DATA [SLOPE PAUSED - FEB 25]  ----------------------------------------------------------------

# NOTE - Script extremely slow, run only if absolutelly necessary

# Empty dataframe for slope estimates

min.slope <- min(grep("^Node", names(all_nodes_raw.DF_aligned)))
max.slope <- max(grep("^NodeTHREE", names(all_nodes_raw.DF_aligned)))


slopeset <- paste("slope", colnames(all_nodes_raw.DF_aligned[min.slope:max.slope]), sep = "_")

#rm(min.slope, max.slope, slopeset)

# Create a temporary data frame to add list of K names to previous data frame 
# Uses row number from a previous step (dataframe)

tempdataframe <- setNames(as.data.frame(matrix(ncol = length(slopeset)+1, nrow = NROW(all_nodes_raw.DF_aligned))),c("Date", slopeset))

tempdataframe["Date"] <-all_nodes_raw.DF_aligned$Date

#all_nodes_raw.DF_K <- cbind(all_nodes_raw.DF_aligned, tempdataframe)


all_nodes_raw.DF[id.tree][i+1,]

all_nodes_raw.DF$ID[i]

names(all_nodes_raw.DF)

list.IDS = colnames(all_nodes_raw.DF_aligned[min.slope:max.slope])


l=1
m=1


idslope=slopeset[l]
id.tree=list.IDS[m]

# Slope with roll slide

as.vector(temp)

rm(x, y)

slope <- function(z) coefficients(lm(y ~ x, as.data.frame(z)))["x"]   

#rollapplyr(zoo(temp), 6, slope, by.column = FALSE)

colnames(filtered_data) <- c("ID", "DOY", "ETo_mm", paste(env_initial$listofnames[i]), "slope")


i=1

for (i in 1:NCOL(tempdataframe)){
  
  idslope=slopeset[i]
  id.tree=list.IDS[i]
  
  temp <- subset(all_nodes_raw.DF,all_nodes_raw.DF[id.tree]>0,  select=c("ID",id.tree, "Date"))
  
  colnames(temp) <- c("x","y", "Date")
  
  temp["slope"]  <- NA
  
  temp$slope[4:NROW(temp)] <-  rollapplyr(zoo(temp[1:2]), 4, slope, by.column = FALSE)
  
  tempdataframe[idslope] <-  merge(all_nodes_raw.DF, temp, by="Date", all.x=T)["slope"]
  
  
  i=i+1
  
}


#plot(tempdataframe$Date, tempdataframe$slope_NodeTWO_3Pine20_02CM_Random_25.30DIAM_)

rm(i, idslope, id.tree, temp)

rm(l, list.IDS, m, max_mV, max_tmax, max.slope, min_mV, min_tmax, min.slope, sequence_DOY, slopeset, Tmaxset, slope, dataframe, Date)


# SLOPE - MERGE ESTIMATED SLOPE WITH [SLOPE PAUSED - FEB 25]   ----------------------------------------------------------------

# Merge
all_nodes_raw.slope.DF <-  cbind(all_nodes_raw.DF, tempdataframe)

names(tempdataframe)

names(all_nodes_raw.DF)

names(all_nodes_raw.slope.DF)


# Save slope files

write.table(all_nodes_raw.slope.DF, "all_nodes_raw.slope.DF.Feb192021.csv", sep=",", row.names = F)


plot(tempdataframe$slope_NodeONE_1Spruce19_02CM_West_15.25DIAM_)


plot(all_nodes_raw.slope.DF$NodeTHREE_40Spruce9_02CM_East_15.25DIAM_)



# END OF SECTION ----

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

# TMAX - EXTRACT RELEVANT VALUES  ----------------------------------------------------------------
# Since no slope, I changed the name of the DF from all_nodes_raw.slope.DF, to all_nodes_raw.DF_aligned

names(all_nodes_raw.DF)

names(all_nodes_raw.DF_aligned)

# Initialize FOR loop

# i = the first column with raw data
# j = the first column wit status
# ii = the first name of column

i=1 # Position of first name of columns to be selected
j=min(grep("Status", names(all_nodes_raw.DF))) # Position of first status
jj=min(grep("slope", names(all_nodes_raw.DF)))


length(grep("Status", names(all_nodes_raw.DF)))
length(grep("NodeOne|NodeTWO|NodeTHREE", names(all_nodes_raw.DF)))

mincol <-min(grep("^NodeONE", names(all_nodes_raw.DF)))
maxcol <-max(grep("^NodeTHREE_58Birch", names(all_nodes_raw.DF)))

all_nodes_raw.DF["Hour"] <- lubridate::hour(all_nodes_raw.DF$Date)
all_nodes_raw.DF["hour.night"] <- ifelse(all_nodes_raw.DF$Hour<7, 22.7, 7.22)

for(i in 1:NROW(env_initial$listofnames)){
  
  # Currently not filtering based on slope [March 02, 2021]
  
  filtered_data <- subset(all_nodes_raw.DF, all_nodes_raw.DF$ETo_mm<0.001&all_nodes_raw.DF[j]==6
                          # &all_nodes_raw.DF[jj]>-0.001&all_nodes_raw.DF[jj]<0.001 
                          &all_nodes_raw.DF$hour.night==22.7, 
                          select=c("ID","DOY", env_initial$listofnames[i]))
  
  #### FIX THE LOESS FUNCTION!!!
  
  DOY <- filtered_data$DOY
  ID <- filtered_data$ID
  S <- filtered_data[env_initial$listofnames[i]]
  S <- unlist(S)
  
  
  # Aggregate data when it is filtered & send to a temp file before aligning to data base
  # I have tried various methods, and using the max is the best option. 
  filtered_data_Tmax <- aggregate(cbind(ID, S)~DOY, data=filtered_data, FUN=max, na.rm=F)
  
  
  all_nodes_raw.DF_aligned[env_initial$listofnames[i]] <- merge(temp_dataframe, filtered_data_Tmax, by="ID", all.x=T, all.y=F, sort=T )[4]
  
  #plot(subset(all_nodes_raw.DF,select=c("ID",env_initial$listofnames[i])), ylim=c(0,1))
  #points(ID, S, type="l", col="red")
  
  i=i+1
  j=j+1
  jj=jj+1
  
}

rm(i, j, jj, mincol, maxcol, filtered_data, filtered_data_Tmax, DOY, ID, S)

head(all_nodes_raw.DF_aligned, 20)

# PLOT - EXTRACTED VALUES FROM PREVIOUS STEP  ----------------------------------------------------------------

#names(all_nodes_raw.DF_aligned)

plot(all_nodes_raw.DF_aligned[,1], all_nodes_raw.DF_aligned[,10],
     ylim=c(0,2))

# If all is Ok, remove all temp.variables



warnings()

# TMAX - ESTIMATE LOESS  ----------------------------------------------------------------

names(Tmaxset)

#all_nodes_raw.DF_aligned <-all_nodes_raw.DF_aligned[1:179]

names(all_nodes_raw.DF_aligned)


Tmax_temp <- setNames(as.data.frame(matrix(ncol = length(Tmaxset)+1, nrow = NROW(all_nodes_raw.DF_aligned))),c(Tmaxset))
Tmax_temp["ID"] <- seq(1, NROW(Tmax_temp))


#names(Tmax_temp)
names(all_nodes_raw.DF_aligned)

i = 8
j = 1
for(i in 8:175){
  
  x = all_nodes_raw.DF_aligned[,"ID"]
  y = all_nodes_raw.DF_aligned[,i]
  
  #plot(x, y)
  # Fit LOESS and save values to "tmax"
  #lo <- loess(y ~ x, na.rm=T, family="symmetric", drop.square = FALSE,method="loess", control=loess.control(surface="direct", statistics="approximate"), span=0.02)
  #out = predict(lo,x,na.rm=T)
  
  f <- approxfun(x, y,method="linear", f=0)
  
  out<-f(x)
  
  Tmax_temp[Tmaxset[j]] <- out
  
  #plot(y)
  #lines(out, col="pink")
  
  j=j+1
  i=i+1
}

#rm(x, y, lo, out, i, j, Date, Tmaxset)
rm(out, x, y, f)

# TMAX - MERGE ESTIMATED LOESS WITH: all_nodes_raw.DF_aligned  ----------------------------------------------------------------

names(all_nodes_raw.DF_aligned)

names(Tmax_temp)

#all_nodes_raw.DF_aligned <- all_nodes_raw.DF_aligned[,1:179]
# Merge
all_nodes_raw.DF_aligned <-  cbind(all_nodes_raw.DF_aligned[1:180], Tmax_temp[1:168])

rm(Tmax_temp)
# PLOT - TEST PREVIOUS STEP  ----------------------------------------------------------------
names(all_nodes_raw.DF_aligned)

plot(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned[,8], col="grey", 
     ylim=c(0,1), 
     xlim=c(1,76000))

lines(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned[,8+173], col="blue")

points(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned[,58], col="grey")

lines(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned[,58+173], col="blue")

points(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned[,169], col="grey")

lines(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned[,169+173], col="blue")

points(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned[,100], col="grey")

lines(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned[,100+173], col="blue")

list(names(all_nodes_raw.DF_aligned))

#all_nodes_raw.DF_aligned[178:1000] = NULL


# TMAX - INSERT GOOD DATA TO: all_nodes_raw.DF_aligned  ----------------------------------------------------------------

# Initialize FOR loop

# j = the first column with status

i=1 # Position of first name of columns to be selected
j=min(grep("Status", names(all_nodes_raw.DF)))

mincol <-min(grep("Pine|Spruce", names(all_nodes_raw.DF)))
maxcol <-max(grep("Pine|Birch|Spruce", names(all_nodes_raw.DF)))
for(a in mincol:maxcol){
  
  # Filter data marked as [good data]
  filtered_data <- subset(all_nodes_raw.DF, all_nodes_raw.DF[j]==6, 
                          select=c("ID","DOY", env_initial$listofnames[i]))
  
  # Merge data (to align) and send to data base
  
  all_nodes_raw.DF_aligned[env_initial$listofnames[i]] <- merge(temp_dataframe, filtered_data, by="ID", all.x=T, all.y=F, sort=T )[4]
  
  j=j+1
  i=i+1
  
}

rm(mincol, maxcol, a, i, j, Tmaxset, filtered_data, Tmax_temp, temp_dataframe, sequence_DOY)
rm(a, filtered_data, filtered_data_Tmax, i, j, DOY, ID, S)

rm(max_mV, max_tmax, min_mV, min_tmax, sequence_DOY, maxcol, mincol, dataframe, temp_dataframe)

rm(Date, jj)

# PLOT - TEST PREVIOUS STEP  ----------------------------------------------------------------

plot(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned[,56], col="grey", 
     ylim=c(0,1), 
     xlim=c(50000,66600))

lines(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned[,56+173], col="blue")

points(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned[,58], col="grey")

lines(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned[,58+173], col="blue")

points(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned[,169], col="grey")

lines(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned[,169+173], col="blue")

points(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned[,8], col="grey")

lines(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned[,8+173], col="blue")
#names(all_nodes_raw.DF_aligned)
#all_nodes_raw.DF_aligned[178:1000] = NULL

# PLOT - SECOND TEST - DIFFERENT SENSORS  ----------------------------------------------------------------

names(all_nodes_raw.DF_aligned)

plot(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned$NodeONE_1Spruce19_02CM_West_15.25DIAM_, 
     xlim = c(1,65000),
     ylim = c(0.1,2),
     col="grey", lty=1, las=1, tcl=+0.2)

lines(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned$Tmax_NodeONE_1Spruce19_02CM_West_15.25DIAM_, type="l", col="blue",pch=16,lty=1, lwd=3)

lines(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned$NodeTWO_1Spruce18_68CM_Random_15.25DIAM_, type="l", col="grey",pch=16,lty=1, lwd=3)

lines(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned$Tmax_NodeTWO_1Spruce18_68CM_Random_15.25DIAM_, type="l", col="blue",pch=16,lty=1, lwd=3)


lines(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned$NodeTHREE_51Birch1_02CM_Random, type="l", col="grey",pch=16,lty=1, lwd=3)

lines(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned$Tmax_NodeTHREE_51Birch1_02CM_Random, type="l", col="blue",pch=16,lty=1, lwd=3)


lines(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned$NodeTHREE_20Pine7_02CM_East_15.25DIAM_, type="l", col="grey",pch=16,lty=1, lwd=3)

lines(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned$Tmax_NodeTHREE_20Pine7_02CM_East_15.25DIAM_, type="l", col="blue",pch=16,lty=1, lwd=3)

lines(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned$NodeTHREE_9Pine2_24CM_Random_25.30DIAM_, type="l", col="grey",pch=16,lty=1, lwd=3)

lines(all_nodes_raw.DF_aligned$ID, all_nodes_raw.DF_aligned$Tmax_NodeTHREE_9Pine2_24CM_Random_25.30DIAM_, type="l", col="blue",pch=16,lty=1, lwd=3)


#write.table(all_nodes_raw.DF_aligned, "all_nodes_raw.DF_aligned.May222020.csv", sep=",", row.names=FALSE)
# END OF SECTION ----

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

# ESTIMATE K  ----------------------------------------------------------------

colnames(all_nodes_raw.DF_aligned)

# Set list of names for K

min_K = min(grep("^Node", names(all_nodes_raw.DF_aligned)))
max_K = max(grep("^NodeTHREE", names(all_nodes_raw.DF_aligned)))


Kset <- paste("K", colnames(all_nodes_raw.DF_aligned[min_K:max_K]), sep = "_")


# Create a temporary data frame to add list of K names to previous data frame 
# Uses row number from a previous step (dataframe)

tempdataframe <- setNames(as.data.frame(matrix(ncol = length(Kset)+1, nrow = NROW(all_nodes_raw.DF_aligned))),c("Date", Kset))
names(tempdataframe)
all_nodes_raw.DF_K <- cbind(all_nodes_raw.DF_aligned, tempdataframe)

colnames(all_nodes_raw.DF_K)


# Initialize For loop
# i = first column with raw mV
# j = first column with Tmax
# ii = first column with K

i=min(grep("^Node", names(all_nodes_raw.DF_K)))
j=min(grep("Tmax", names(all_nodes_raw.DF_K)))
ii=min(grep("^K_", names(all_nodes_raw.DF_K)))

# FOR loop starts 

for(i in min(grep("^Node", names(all_nodes_raw.DF_K))):max(grep("^Node", names(all_nodes_raw.DF_K))) ){
  
  K <- (all_nodes_raw.DF_K[j]-all_nodes_raw.DF_K[i]) / all_nodes_raw.DF_K[i]
  
  all_nodes_raw.DF_K[ii] <- K
  
  i=i+1
  j=j+1
  ii=ii+1
  
}

# End For loop

rm(tempdataframe,K, i, ii, j, Kset, max_K, min_K)


head(all_nodes_raw.DF_K$Date, 5)

colnames(all_nodes_raw.DF_K)


# PLOT TO SEE IF K LOOP WORKED  ----------------------------------------------------------------



fdcol <- palette(rainbow(25))  

plot(all_nodes_raw.DF_K$ID, all_nodes_raw.DF_K$K_NodeONE_1Spruce19_02CM_West_15.25DIAM_,type="l",
     ylim = c(-0.1,0.5),
     xlim = c(1,74000),
     col="grey", lty=1, las=1, tcl=+0.2)

lines(all_nodes_raw.DF_K$ID, all_nodes_raw.DF_K$K_NodeONE_9Spruce10_02CM_Random_15.25DIAM_, type="l", col=fdcol[1],pch=16,lty=1, lwd=2)

lines(all_nodes_raw.DF_K$ID, all_nodes_raw.DF_K$K_NodeONE_43Pine13_02CM_Random_15.25DIAM_, type="l", col=fdcol[2],pch=16,lty=1, lwd=2)
lines(all_nodes_raw.DF_K$ID, all_nodes_raw.DF_K$K_NodeTWO_26Pine12_46CM_Random_25.30DIAM_, type="l", col=fdcol[3],pch=16,lty=1, lwd=2)
lines(all_nodes_raw.DF_K$ID, all_nodes_raw.DF_K$K_NodeTWO_28Pine13_02CM_Random_7.15DIAM_, type="l", col=fdcol[4],pch=16,lty=1, lwd=2)
lines(all_nodes_raw.DF_K$ID, all_nodes_raw.DF_K$K_NodeTWO_32Pine9_24CM_Random_7.15DIAM_, type="l", col=fdcol[5],pch=16,lty=1, lwd=2)
lines(all_nodes_raw.DF_K$ID, all_nodes_raw.DF_K$K_NodeONE_51Birch1_02CM_Random, type="l", col=fdcol[6],pch=16,lty=1, lwd=2)
lines(all_nodes_raw.DF_K$ID, all_nodes_raw.DF_K$K_NodeONE_52Birch2_02CM_Random, type="l", col=fdcol[7],pch=16,lty=1, lwd=2)


rm(fdcol)
#write.table(all_nodes_raw.DF_K, "all_nodes_raw.DF_K_May22_2020.csv", sep=",", row.names=F)



# END OF SECTION ----

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

# ESTIMATE SAP FLUX DENSITY  ----------------------------------------------------------------

# Set list of names for flux density

min_Fd = min(grep("^Node", names(all_nodes_raw.DF_K)))
max_Fd = max(grep("^Node", names(all_nodes_raw.DF_K)))


Fdset <- paste("Fd", colnames(all_nodes_raw.DF_K[min_Fd:max_Fd]), sep = "_")


# Create a temporary data frame to add list of K names to previous data frame 
# Uses row number from a previous step (dataframe)

length(Fdset)
names(all_nodes_raw.DF_K)

Fd_temp <- setNames(as.data.frame(matrix(ncol = length(Fdset), nrow = NROW(all_nodes_raw.DF_K))),c(Fdset))

all.nodes.raw.tmax.k.fd <- cbind(all_nodes_raw.DF_K, Fd_temp)

# Parameters 
rm(Fd_temp)


#a = 82.14
#b = 1.105

#(118*0.6^b)*3600 / 10000
#(42.48*0.6^b) 

#a = 42.84
#a = 119
#b = 1.231

# Updated parameters from Peter Richards et al. 2018

a = 85.25
b = 1.55

# Initialize For loop

min_KQ = min(grep("^K_", names(all.nodes.raw.tmax.k.fd)))
max_KQ = max(grep("^K_", names(all.nodes.raw.tmax.k.fd)))

min_FdQ = min(grep("^Fd_", names(all.nodes.raw.tmax.k.fd)))

i=min_KQ
ii=min_FdQ
colnames(all.nodes.raw.tmax.k.fd)
# FOR loop starts 

for(i in min_KQ:max_KQ ){
  
  species <- ifelse("Pine" %in% colnames(all.nodes.raw.tmax.k.fd[i]), "Pine", "Spruce")
  
  if(species=="Pine"){
    a=85.25
    b=1.55
  }else{
    a=85.25 
    b=1.55}
  
  Fd_temp <- (a * all.nodes.raw.tmax.k.fd[i]^b)#*3600 / 10000
  
  all.nodes.raw.tmax.k.fd[ii] <- Fd_temp
  
  i=i+1
  ii=ii+1
  
}

# End For loop

all.nodes.raw.tmax.k.fd["ID"] <- seq(1, NROW(all.nodes.raw.tmax.k.fd))

rm(Fdset, a, b, i, ii, min_Fd, max_Fd, max_KQ, min_KQ, min_FdQ)
rm(species, Fd_temp)


# PLOT TO TEST IF LOOP WORKED  ----------------------------------------------------------------


# Plot to see if flux density loop worked
fdcol <- palette(rainbow(25))  
plot(all.nodes.raw.tmax.k.fd$ID, all.nodes.raw.tmax.k.fd$Fd_NodeONE_1Spruce19_02CM_West_15.25DIAM_,
     ylim = c(0,35),
     xlim = c(1,76000),
     col=fdcol[1], lty=1, las=1, tcl=+0.2, type="l")

names(all.nodes.raw.tmax.k.fd)

lines(all.nodes.raw.tmax.k.fd$ID, all.nodes.raw.tmax.k.fd$Fd_NodeONE_43Pine13_02CM_Random_15.25DIAM_, type="l", col=fdcol[2], pch=16,lty=1, lwd=2)

lines(all.nodes.raw.tmax.k.fd$ID, all.nodes.raw.tmax.k.fd$Fd_NodeTWO_1Spruce18_68CM_Random_15.25DIAM_, type="l", col=fdcol[3],pch=16,lty=1, lwd=2)

lines(all.nodes.raw.tmax.k.fd$ID, all.nodes.raw.tmax.k.fd$Fd_NodeTWO_50Spruce6_24CM_Random_15.25DIAM_, type="l", col=fdcol[4],pch=16,lty=1, lwd=2)

lines(all.nodes.raw.tmax.k.fd$ID, all.nodes.raw.tmax.k.fd$Fd_NodeTHREE_9Pine2_24CM_Random_25.30DIAM_, type="l", col=fdcol[5],pch=16,lty=1, lwd=2)

lines(all.nodes.raw.tmax.k.fd$ID, all.nodes.raw.tmax.k.fd$Fd_NodeONE_51Birch1_02CM_25DIAM_, type="l", col=fdcol[6],pch=16,lty=1, lwd=2)
lines(all.nodes.raw.tmax.k.fd$ID, all.nodes.raw.tmax.k.fd$Fd_NodeONE_52Birch2_02CM_15DIAM_, type="l", col=fdcol[7],pch=16,lty=1, lwd=2)

lines(all.nodes.raw.tmax.k.fd$ID, all.nodes.raw.tmax.k.fd$Fd_NodeTHREE_51Birch1_02CM_Random, type="l", col=fdcol[8],pch=16,lty=1, lwd=2)
lines(all.nodes.raw.tmax.k.fd$ID, all.nodes.raw.tmax.k.fd$Fd_NodeTHREE_55Birch5_02CM_Random, type="l", col=fdcol[9],pch=16,lty=1, lwd=2)

lines(all.nodes.raw.tmax.k.fd$ID, all.nodes.raw.tmax.k.fd$Fd_NodeTHREE_56Birch6_02CM_Random, type="l", col=fdcol[10],pch=16,lty=1, lwd=2)
lines(all.nodes.raw.tmax.k.fd$ID, all.nodes.raw.tmax.k.fd$Fd_NodeTHREE_57Birch7_02CM_Random, type="l", col=fdcol[11],pch=16,lty=1, lwd=2)

rm(fdcol)

names(all.nodes.raw.tmax.k.fd)



# DELETE: ROWS FROM FAKE DATA  ----------------------------------------------------------------

names(all.nodes.raw.tmax.k.fd)

# Birch Node ONE

grep("NodeONE.*Birch", names(all_nodes_raw.DF_Fd))

all.nodes.raw.tmax.k.fd[1:50400,(grep("NodeONE.*Birch", names(all.nodes.raw.tmax.k.fd)))]=NA

# Pine&Sruce Node TWO
grep("NodeTWO.*(51|52|53|54|55|56|57|58)", names(all.nodes.raw.tmax.k.fd))

all.nodes.raw.tmax.k.fd[1:52992,(grep("NodeTWO.*(51|52|53|54|55|56|57|58)", names(all.nodes.raw.tmax.k.fd)))]=NA

#Birch Node THREE

grep("NodeTHREE.*Birch", names(all_nodes_raw.DF_Fd))

all.nodes.raw.tmax.k.fd[1:49968,(grep("NodeTHREE.*Birch", names(all.nodes.raw.tmax.k.fd)))]=NA

# Delete all clear non-Fd columns

grep("Date|^Fd_", names(all.nodes.raw.tmax.k.fd))

#write.table(all.nodes.raw.tmax.k.fd, "all.nodes.raw.tmax.k.fd_Feb212020.csv", sep=" ", row.names = F)


# EXTRACT FD COLUMNS FROM - all.nodes.raw.tmax.k.fd  ----------------------------------------------------------------


all.nodes.raw.Fd <- all.nodes.raw.tmax.k.fd[,(grep("Date|^Fd_", names(all.nodes.raw.tmax.k.fd)))]


names(all.nodes.raw.Fd) 

NROW(all.nodes.raw.Fd)
# SAVE - FD DATA BASE  ----------------------------------------------------------------

#colnames(all.nodes.raw.Fd)

#Export database for further cleaning
write.table(all.nodes.raw.Fd, "all.nodes.raw.Fd.csv", sep=",",row.names = FALSE)
# write.table(all.nodes.raw.Fd, "all.nodes.raw.Fd_Feb_26_2021.csv", sep=",",row.names = FALSE)

# OPEN - RAW FD DATABASE ----------------------------------------------------------------

all.nodes.raw.Fd = read.table(file="/Users/JA/Documents/PostDoc/sapflow_data/all.nodes.raw.Fd_Mar_01_2021.fun.max.csv",header=T,sep=",")

all.nodes.raw.Fd["Date"] <- seq(as.POSIXct("2016/06/17 00:00:00", tz="UTC"), as.POSIXct("2020/12/30 23:30:00", tz="UTC"), by = "30 min")

head(all.nodes.raw.Fd,10)

# ADD - DATE-YEAR-DOY to FD Database ----------------------------------------------------------------


all.nodes.raw.Fd["Date"] <- seq(as.POSIXct("2016/06/17 00:00:00", tz="UTC"), as.POSIXct("2020/12/30 23:30:00", tz="UTC"), by = "30 min")
all.nodes.raw.Fd["Year"] <- lubridate::year(all.nodes.raw.Fd$Date)
all.nodes.raw.Fd["DOY"] <- ifelse(all.nodes.raw.Fd$Year==2016,lubridate::yday(all.nodes.raw.Fd$Date),
                                  ifelse((all.nodes.raw.Fd$Year==2017),lubridate::yday(all.nodes.raw.Fd$Date)+366,
                                         ifelse((all.nodes.raw.Fd$Year==2018),lubridate::yday(all.nodes.raw.Fd$Date)+731,
                                                ifelse((all.nodes.raw.Fd$Year==2019),lubridate::yday(all.nodes.raw.Fd$Date)+1096,
                                                       lubridate::yday(all.nodes.raw.Fd$Date)+1461))))


# END OF SECTION ----
000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#
#  WOUNDING DRIFT CORRECTION
#
000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

# Correction of reductions in sap flux density, resulting from wounding
# CORRECT - WOUNDING DRIFFT  ----------------------------------------------------------------

#temp.names <- paste(paste("Fd", tree.info$Unique.name.ID[grep("Pine.*_02CM_.*|Spruce.*_02CM_.*", tree.info$Unique.name.ID)], sep = "_"), "radial", sep=".")

names(all.nodes.raw.Fd)

slope.set <- as.data.frame(c(names(all.nodes.raw.Fd[3:107]),names(all.nodes.raw.Fd[116:170])))

slope.set <- setNames(as.data.frame(matrix(ncol = 2, nrow = 160)),c("ID", "Species"))

slope.set["ID"] <- as.data.frame(c(names(all.nodes.raw.Fd[3:107]),names(all.nodes.raw.Fd[116:170])))

#slope.set["species"] <- NA

slope.set$Species[grep("Pine", slope.set$ID)]

slope.set["Species"] <- NA 

slope.set[grep("Pine",slope.set$ID), "Species"] <- "Pine"
slope.set[grep("Spruce",slope.set$ID), "Species"] <- "Spruce"
slope.set[grep("Birch",slope.set$ID), "Species"] <- "Birch"

names.slope <- c(names(all.nodes.raw.Fd[3:107]),names(all.nodes.raw.Fd[116:170]))
wound.drift.fd <- setNames(as.data.frame(matrix(ncol = 160+1, nrow = NROW(all.nodes.raw.Fd))),c(names.slope))

rm(names.slope)

# Estimate moving slow from 1, to max slope (from linear equations)

slope.pine <- seq(1, par.estimates[9,2], length.out = NROW(all.nodes.raw.Fd))
slope.spruce <- seq(1, par.estimates[10,2], length.out = NROW(all.nodes.raw.Fd))
slope.birch <- rep(1, NROW(all.nodes.raw.Fd))

names(all.nodes.raw.Fd)
names(wound.drift.fd)

NCOL(all.nodes.raw.Fd)
NCOL(wound.drift.fd)
i = 1
#ii=2


for(i in 1:length(slope.set$ID)){
  
  tree.id <- as.character(slope.set$ID[i])
  fd <- all.nodes.raw.Fd[tree.id]
  species <- slope.set$Species[i]
  
  
  adj.slope <- if(species == "Pine"){slope.pine
  }else if(species == "Spruce"){slope.spruce
  }else if(species == "Birch"){slope.birch}
  
  wound.drift.fd[tree.id] <- 0+adj.slope*fd
  
  #str(tree.id)
  
  i=i+1
}


head(wound.drift.fd, 20)

wound.drift.fd["Date"] <- seq(as.POSIXct("2016/06/17 00:00:00", tz="UTC"), as.POSIXct("2020/12/30 23:30:00", tz="UTC"), by = "30 min")


#plot(wound.drift.fd[,1], wound.drift.fd[,2])

rm(i, fd, tree.id, species, adj.slope, slope.pine, slope.spruce, slope.birch)
rm(slope.set)


#rm(intercept, p.pine, adj.factor, F.pine, rsqr.pine, temp.del, x, y, xlim_fd, ylim_fd)



# SAVE - WOUND.DRIFT DATABASE  ----------------------------------------------------------------

#rm(all.nodes.raw.Fd.rad)
all.nodes.raw.Fd.rad <-  merge(all.nodes.raw.Fd, wound.drift.fd, by="Date", all.x=T, all.y=T, sort=T )

write.table(wound.drift.fd, "Data/wound.drift.fd.csv", sep=",", row.names = FALSE)

# END OF SECTION ----

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000


000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#
#  RADIAL PROFILE PARAMETERS
#
000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

# In this section I fitted one radial profile every 30-minute interval
# A model is used to estimate the parameters for a normalized radial profile (i.e., find the shape of the profile)
# Parameters are used on trees that do not have sensors at various depths to predict expected radial profile
# Radial profiles estimated for Large, Medium and Small trees (both species)
# PREPARE: DAta for radial profiles  ----------------------------------------------------------------


# Database to populate with data and parameters

# Empty database to populate with parameters

dt.ratio.pars <- as.data.frame(seq(as.POSIXct("2016/06/17 00:00:00", tz="UTC"), as.POSIXct("2020/12/30 23:30:00", tz="UTC"), by = "30 min"))

colnames(dt.ratio.pars)[1] <- c("Date")

0000000000000000000000000000000000

# Parameters for PINE - LARGE 

0000000000000000000000000000000000

rm(temp.dt)

## LARGE PINE


depth1cm.ag <- rowMeans(wound.drift.fd[grep("Pine.*_02CM_.*_25.30DIAM_", names(wound.drift.fd))], na.rm=T)
depth3cm.ag <- rowMeans(wound.drift.fd[grep("Pine.*_24CM_.*_25.30DIAM_", names(wound.drift.fd))], na.rm=T)
depth5cm.ag <- rowMeans(wound.drift.fd[grep("Pine.*_46CM_.*_25.30DIAM_", names(wound.drift.fd))], na.rm=T)
depth7cm.ag <- rowMeans(wound.drift.fd[grep("Pine.*_68CM_.*_25.30DIAM_", names(wound.drift.fd))], na.rm=T)

nrow(wound.drift.fd)

dt.ratio.pars["p.1.lg"] <- depth1cm.ag
dt.ratio.pars["p.3.lg"] <- depth3cm.ag
dt.ratio.pars["p.5.lg"] <- depth5cm.ag
dt.ratio.pars["p.7.lg"] <- depth7cm.ag

dt.ratio.pars["p.ra.1.lg"] <- ifelse((dt.ratio.pars$p.1.lg/dt.ratio.pars$p.1.lg)>1, NA, dt.ratio.pars$p.1.lg/dt.ratio.pars$p.1.lg)
dt.ratio.pars["p.ra.3.lg"] <- ifelse((dt.ratio.pars$p.3.lg/dt.ratio.pars$p.1.lg)>1,  NA, dt.ratio.pars$p.3.lg/dt.ratio.pars$p.1.lg)
dt.ratio.pars["p.ra.5.lg"] <- ifelse((dt.ratio.pars$p.5.lg/dt.ratio.pars$p.1.lg)>1, NA, dt.ratio.pars$p.5.lg/dt.ratio.pars$p.1.lg)
dt.ratio.pars["p.ra.7.lg"] <- ifelse((dt.ratio.pars$p.7.lg/dt.ratio.pars$p.1.lg)>1, NA, dt.ratio.pars$p.7.lg/dt.ratio.pars$p.1.lg)

# Delete all files created

rm(depth1cm, depth3cm, depth5cm, depth7cm)
rm(depth1cm.ag, depth3cm.ag, depth5cm.ag, depth7cm.ag)

## MEDIUM PINE

depth1cm.ag <- rowMeans(wound.drift.fd[grep("Pine.*_02CM_.*_15.25DIAM_", names(wound.drift.fd))], na.rm=T)
depth3cm.ag <- rowMeans(wound.drift.fd[grep("Pine.*_24CM_.*_15.25DIAM_", names(wound.drift.fd))], na.rm=T)
depth5cm.ag <- rowMeans(wound.drift.fd[grep("Pine.*_46CM_.*_15.25DIAM_", names(wound.drift.fd))], na.rm=T)
depth7cm.ag <- rowMeans(wound.drift.fd[grep("Pine.*_68CM_.*_15.25DIAM_", names(wound.drift.fd))], na.rm=T)

dt.ratio.pars["p.1.md"] <- depth1cm.ag
dt.ratio.pars["p.3.md"] <- depth3cm.ag
dt.ratio.pars["p.5.md"] <- depth5cm.ag
dt.ratio.pars["p.7.md"] <- depth7cm.ag

dt.ratio.pars["p.ra.1.md"] <- ifelse((dt.ratio.pars$p.1.md/dt.ratio.pars$p.1.md)>1, NA, dt.ratio.pars$p.1.md/dt.ratio.pars$p.1.md)
dt.ratio.pars["p.ra.3.md"] <- ifelse((dt.ratio.pars$p.3.md/dt.ratio.pars$p.1.md)>1,  NA, dt.ratio.pars$p.3.md/dt.ratio.pars$p.1.md)
dt.ratio.pars["p.ra.5.md"] <- ifelse((dt.ratio.pars$p.5.md/dt.ratio.pars$p.1.md)>1, NA, dt.ratio.pars$p.5.md/dt.ratio.pars$p.1.md)
dt.ratio.pars["p.ra.7.md"] <- ifelse((dt.ratio.pars$p.7.md/dt.ratio.pars$p.1.md)>1, NA, dt.ratio.pars$p.7.md/dt.ratio.pars$p.1.md)

# Delete all files created

rm(depth1cm, depth3cm, depth5cm, depth7cm)
rm(depth1cm.ag, depth3cm.ag, depth5cm.ag, depth7cm.ag)

## SMALL PINE

depth1cm.ag <- rowMeans(wound.drift.fd[grep("Pine.*_02CM_.*_7.15DIAM_", names(wound.drift.fd))], na.rm=T)
depth3cm.ag <- rowMeans(wound.drift.fd[grep("Pine.*_24CM_.*_7.15DIAM_", names(wound.drift.fd))], na.rm=T)


dt.ratio.pars["p.1.sm"] <- depth1cm.ag
dt.ratio.pars["p.3.sm"] <- depth3cm.ag

dt.ratio.pars["p.ra.1.sm"] <- ifelse((dt.ratio.pars$p.1.sm/dt.ratio.pars$p.1.sm)>1, NA, dt.ratio.pars$p.1.sm/dt.ratio.pars$p.1.sm)
dt.ratio.pars["p.ra.3.sm"] <- ifelse((dt.ratio.pars$p.3.sm/dt.ratio.pars$p.1.sm)>1,  NA, dt.ratio.pars$p.3.sm/dt.ratio.pars$p.1.sm)

# Delete all files created

rm(depth1cm, depth3cm, depth5cm, depth7cm)
rm(depth1cm.ag, depth3cm.ag, depth5cm.ag, depth7cm.ag)

par(mar=c(3, 3, 0, 0))

plot(dt.ratio.pars$Date[1:75000], dt.ratio.pars$p.ra.1.lg[1:75000], ylim=c(0,1))
lines(dt.ratio.pars$Date, dt.ratio.pars$p.ra.3.lg, col="red")
lines(dt.ratio.pars$Date, dt.ratio.pars$p.ra.5.lg, col="blue")
lines(dt.ratio.pars$Date, dt.ratio.pars$p.ra.7.lg, col="green")


plot(dt.ratio.pars$Date[1:75000], dt.ratio.pars$p.1.lg[1:75000], ylim=c(0,20))
lines(dt.ratio.pars$Date, dt.ratio.pars$p.3.lg, col="red")
lines(dt.ratio.pars$Date, dt.ratio.pars$p.5.lg, col="blue")
lines(dt.ratio.pars$Date, dt.ratio.pars$p.7.lg, col="green")

# Delete all files created

rm(depth1cm, depth3cm, depth5cm, depth7cm)
rm(depth1cm.ag, depth3cm.ag, depth5cm.ag, depth7cm.ag)

## LARGE SPRUCE


depth1cm.ag <- rowMeans(wound.drift.fd[grep("Spruce.*_02CM_.*_25.30DIAM_", names(wound.drift.fd))], na.rm=T)
depth3cm.ag <- rowMeans(wound.drift.fd[grep("Spruce.*_24CM_.*_25.30DIAM_", names(wound.drift.fd))], na.rm=T)
depth5cm.ag <- rowMeans(wound.drift.fd[grep("Spruce.*_46CM_.*_25.30DIAM_", names(wound.drift.fd))], na.rm=T)
depth7cm.ag <- rowMeans(wound.drift.fd[grep("Spruce.*_68CM_.*_25.30DIAM_", names(wound.drift.fd))], na.rm=T)


dt.ratio.pars["s.1.lg"] <- depth1cm.ag
dt.ratio.pars["s.3.lg"] <- depth3cm.ag
dt.ratio.pars["s.5.lg"] <- depth5cm.ag
dt.ratio.pars["s.7.lg"] <- depth7cm.ag

dt.ratio.pars["s.ra.1.lg"] <- ifelse((dt.ratio.pars$s.1.lg/dt.ratio.pars$s.1.lg)>1, NA, dt.ratio.pars$s.1.lg/dt.ratio.pars$s.1.lg)
dt.ratio.pars["s.ra.3.lg"] <- ifelse((dt.ratio.pars$s.3.lg/dt.ratio.pars$s.1.lg)>1,  NA, dt.ratio.pars$s.3.lg/dt.ratio.pars$s.1.lg)
dt.ratio.pars["s.ra.5.lg"] <- ifelse((dt.ratio.pars$s.5.lg/dt.ratio.pars$s.1.lg)>1, NA, dt.ratio.pars$s.5.lg/dt.ratio.pars$s.1.lg)
dt.ratio.pars["s.ra.7.lg"] <- ifelse((dt.ratio.pars$s.7.lg/dt.ratio.pars$s.1.lg)>1, NA, dt.ratio.pars$s.7.lg/dt.ratio.pars$s.1.lg)

# Delete all files created

rm(depth1cm, depth3cm, depth5cm, depth7cm)
rm(depth1cm.ag, depth3cm.ag, depth5cm.ag, depth7cm.ag)

## MEDIUM SPRUCE


depth1cm.ag <- rowMeans(wound.drift.fd[grep("Spruce.*_02CM_.*_15.25DIAM_", names(wound.drift.fd))], na.rm=T)
depth3cm.ag <- rowMeans(wound.drift.fd[grep("Spruce.*_24CM_.*_15.25DIAM_", names(wound.drift.fd))], na.rm=T)
depth5cm.ag <- rowMeans(wound.drift.fd[grep("Spruce.*_46CM_.*_15.25DIAM_", names(wound.drift.fd))], na.rm=T)
depth7cm.ag <- rowMeans(wound.drift.fd[grep("Spruce.*_68CM_.*_15.25DIAM_", names(wound.drift.fd))], na.rm=T)

dt.ratio.pars["s.1.md"] <- depth1cm.ag
dt.ratio.pars["s.3.md"] <- depth3cm.ag
dt.ratio.pars["s.5.md"] <- depth5cm.ag
dt.ratio.pars["s.7.md"] <- depth7cm.ag

dt.ratio.pars["s.ra.1.md"] <- ifelse((dt.ratio.pars$s.1.md/dt.ratio.pars$s.1.md)>1, NA, dt.ratio.pars$s.1.md/dt.ratio.pars$s.1.md)
dt.ratio.pars["s.ra.3.md"] <- ifelse((dt.ratio.pars$s.3.md/dt.ratio.pars$s.1.md)>1,  NA, dt.ratio.pars$s.3.md/dt.ratio.pars$s.1.md)
dt.ratio.pars["s.ra.5.md"] <- ifelse((dt.ratio.pars$s.5.md/dt.ratio.pars$s.1.md)>1, NA, dt.ratio.pars$s.5.md/dt.ratio.pars$s.1.md)
dt.ratio.pars["s.ra.7.md"] <- ifelse((dt.ratio.pars$s.7.md/dt.ratio.pars$s.1.md)>1, NA, dt.ratio.pars$s.7.md/dt.ratio.pars$s.1.md)

# Delete all files created

rm(depth1cm, depth3cm, depth5cm, depth7cm)
rm(depth1cm.ag, depth3cm.ag, depth5cm.ag, depth7cm.ag)

## SMALL SPRUCE

depth1cm.ag <- rowMeans(wound.drift.fd[grep("Spruce.*_02CM_.*_7.15DIAM_", names(wound.drift.fd))], na.rm=T)
depth3cm.ag <- rowMeans(wound.drift.fd[grep("Spruce.*_24CM_.*_7.15DIAM_", names(wound.drift.fd))], na.rm=T)

dt.ratio.pars["s.1.sm"] <- depth1cm.ag
dt.ratio.pars["s.3.sm"] <- depth3cm.ag

dt.ratio.pars["s.ra.1.sm"] <- ifelse((dt.ratio.pars$s.1.sm/dt.ratio.pars$s.1.sm)>1, NA, dt.ratio.pars$s.1.sm/dt.ratio.pars$s.1.sm)
dt.ratio.pars["s.ra.3.sm"] <- ifelse((dt.ratio.pars$s.3.sm/dt.ratio.pars$s.1.sm)>1,  NA, dt.ratio.pars$s.3.sm/dt.ratio.pars$s.1.sm)

# Delete all files created

rm(depth1cm, depth3cm, depth5cm, depth7cm)
rm(depth1cm.ag, depth3cm.ag, depth5cm.ag, depth7cm.ag)



par(mar=c(3, 3, 0, 0))

plot(dt.ratio.pars$Date[1:5000], dt.ratio.pars$s.ra.1.lg[1:5000], ylim=c(0,2))
lines(dt.ratio.pars$Date, dt.ratio.pars$s.ra.3.lg, col="red")
lines(dt.ratio.pars$Date, dt.ratio.pars$s.ra.5.lg, col="blue")
lines(dt.ratio.pars$Date, dt.ratio.pars$s.ra.7.lg, col="green")


plot(dt.ratio.pars$Date[1:5000], dt.ratio.pars$s.1.lg[1:5000], ylim=c(0,20))
lines(dt.ratio.pars$Date, dt.ratio.pars$s.3.lg, col="red")
lines(dt.ratio.pars$Date, dt.ratio.pars$s.5.lg, col="blue")
lines(dt.ratio.pars$Date, dt.ratio.pars$s.7.lg, col="green")

# Delete all files created

rm(depth1cm, depth3cm, depth5cm, depth7cm)
rm(depth1cm.ag, depth3cm.ag, depth5cm.ag, depth7cm.ag)

# Fit parameters for each measurement 
# Test the graph with the small subset (else it might collapse)

# ESTIMATE: Radial profile parameters. Every 30 minutes  ----------------------------------------------------------------

# Try to fit parametres for a small subset

head(dt.ratio.pars, 20)

temp.dt.ratio.pars <- dt.ratio.pars

temp.dt.ratio.pars["pine.a"] <- NA
temp.dt.ratio.pars["pine.b"] <- NA
temp.dt.ratio.pars["pine.c"] <- NA
temp.dt.ratio.pars["pine.d"] <- NA
temp.dt.ratio.pars["pine.e"] <- NA

temp.dt.ratio.pars["spruce.a"] <- NA
temp.dt.ratio.pars["spruce.b"] <- NA
temp.dt.ratio.pars["spruce.c"] <- NA
temp.dt.ratio.pars["spruce.d"] <- NA
temp.dt.ratio.pars["spruce.e"] <- NA

colfunc <- colorRampPalette(c("blueviolet", "cyan", "blueviolet", "red"))


# Main loop for pine - parameters
#plot(temp.dt$Depth, temp.dt$Js, xlim=c(0,1))
i=1
NROW(temp.dt.ratio.pars)
for (i in 1:NROW(temp.dt.ratio.pars)){
  #for (i in 1:200){
  pars.temp <- setNames(as.data.frame(matrix(ncol = 4, nrow = 13)),c("Depth","Js", "Ratio", "Jsnorm"))
  
  pars.temp[1,] <- c(0.2/8, temp.dt.ratio.pars$p.1.lg[i], temp.dt.ratio.pars$p.ra.1.lg[i], NA) # Extra
  pars.temp[2,] <- c(1/8, temp.dt.ratio.pars$p.1.lg[i], temp.dt.ratio.pars$p.ra.1.lg[i], NA)
  pars.temp[3,] <- c(3/8, temp.dt.ratio.pars$p.3.lg[i], temp.dt.ratio.pars$p.ra.3.lg[i], NA)
  pars.temp[4,] <- c(5/8, temp.dt.ratio.pars$p.5.lg[i], temp.dt.ratio.pars$p.ra.5.lg[i], NA)
  pars.temp[5,] <- c(7/8, temp.dt.ratio.pars$p.7.lg[i], temp.dt.ratio.pars$p.ra.7.lg[i], NA)
  pars.temp[6,] <- c(1/7, temp.dt.ratio.pars$p.1.md[i], temp.dt.ratio.pars$p.ra.1.md[i], NA)
  pars.temp[7,] <- c(0.2/7, temp.dt.ratio.pars$p.1.md[i], temp.dt.ratio.pars$p.ra.1.md[i], NA) # Extra
  pars.temp[8,] <- c(3/7, temp.dt.ratio.pars$p.3.md[i], temp.dt.ratio.pars$p.ra.3.md[i], NA)
  pars.temp[9,] <- c(5/7, temp.dt.ratio.pars$p.5.md[i], temp.dt.ratio.pars$p.ra.5.md[i], NA)
  pars.temp[10,] <- c(7/7, temp.dt.ratio.pars$p.7.md[i], temp.dt.ratio.pars$p.ra.7.md[i], NA)
  pars.temp[11,] <- c(1/3, temp.dt.ratio.pars$p.1.sm[i], temp.dt.ratio.pars$p.ra.1.sm[i], NA)
  pars.temp[12,] <- c(0.99,0,0, 0)
  pars.temp[13,] <- c(1,0,0, 0)
  
  pars.temp$Jsnorm <- pars.temp$Js/max(pars.temp$Js, na.rm=T)
  
  if(is.na(pars.temp$Ratio[1])){
    next
  }
  
  ## Fit equation
  
  x = pars.temp$Depth
  y = pars.temp$Ratio
  
  b0 <- summary(model)$parameters[1,1]
  b1 <- summary(model)$parameters[2,1]
  
  try({  
    
    model <- nlsLM(y ~ b0*exp(-0.5*((x-0.125)/b1)^2), 
                   control = nls.lm.control(maxiter= 200),data=pars.temp, 
                   start = list(b0=b0, b1=b1), na.action=NULL)
    
  })
  
  temp.dt.ratio.pars$pine.a[i] <- summary(model)$parameters[1,1]
  temp.dt.ratio.pars$pine.b[i] <- summary(model)$parameters[2,1]
  
  i=i+1
  
}


rm(b0, b1, b2, b3, b4, hour, i, js, month, x, y, colfunc)
rm(model, temp.dt, temp.rd, temp.rd.agg, djs.temp)

plot(temp.dt.ratio.pars$Date[1:1000], temp.dt.ratio.pars$pine.a[1:1000], ylim=c(-1,2), col="red")
lines(temp.dt.ratio.pars$Date[1:1000], temp.dt.ratio.pars$pine.b[1:1000], col="blue")



# Main loop for spruce - parameters


temp.dt.ratio.pars$spruce.a <- NA
temp.dt.ratio.pars$spruce.a <- NA

i=1
NROW(temp.dt.ratio.pars)
for (i in 1:NROW(temp.dt.ratio.pars)){
  #  for (i in 1:900){
  pars.temp <- setNames(as.data.frame(matrix(ncol = 4, nrow = 13)),c("Depth","Js", "Ratio", "Jsnorm"))
  
  pars.temp[1,] <- c(0.2/8, temp.dt.ratio.pars$s.1.lg[i], temp.dt.ratio.pars$s.ra.1.lg[i], NA) # Extra
  pars.temp[2,] <- c(1/8, temp.dt.ratio.pars$s.1.lg[i], temp.dt.ratio.pars$s.ra.1.lg[i], NA)
  pars.temp[3,] <- c(3/8, temp.dt.ratio.pars$s.3.lg[i], temp.dt.ratio.pars$s.ra.3.lg[i], NA)
  pars.temp[4,] <- c(5/8, temp.dt.ratio.pars$s.5.lg[i], temp.dt.ratio.pars$s.ra.5.lg[i], NA)
  pars.temp[5,] <- c(7/8, temp.dt.ratio.pars$s.7.lg[i], temp.dt.ratio.pars$s.ra.7.lg[i], NA)
  pars.temp[6,] <- c(1/7, temp.dt.ratio.pars$s.1.md[i], temp.dt.ratio.pars$s.ra.1.md[i], NA)
  pars.temp[7,] <- c(0.2/7, temp.dt.ratio.pars$s.1.md[i], temp.dt.ratio.pars$s.ra.1.md[i], NA) # Extra
  pars.temp[8,] <- c(3/7, temp.dt.ratio.pars$s.3.md[i], temp.dt.ratio.pars$s.ra.3.md[i], NA)
  pars.temp[9,] <- c(5/7, temp.dt.ratio.pars$s.5.md[i], temp.dt.ratio.pars$s.ra.5.md[i], NA)
  pars.temp[10,] <- c(7/7, temp.dt.ratio.pars$s.7.md[i], temp.dt.ratio.pars$s.ra.7.md[i], NA)
  pars.temp[11,] <- c(1/3, temp.dt.ratio.pars$s.1.sm[i], temp.dt.ratio.pars$s.ra.1.sm[i], NA)
  pars.temp[12,] <- c(0.99,0,0, 0)
  pars.temp[13,] <- c(1,0,0, 0)
  
  pars.temp$Jsnorm <- pars.temp$Js/max(pars.temp$Js, na.rm=T)
  
  if(is.na(pars.temp$Ratio[1])){
    next
  }
  
  ## Fit equation
  
  x = pars.temp$Depth
  y = pars.temp$Ratio
  
  b0 <- 0.4
  b1 <- 0.7
  
  try({  
    
    model <- nlsLM(y ~ b0*exp(-0.5*((x-0.125)/b1)^2), 
                   control = nls.lm.control(maxiter= 200),data=pars.temp, 
                   start = list(b0=0.4, b1=0.7), na.action=NULL)
    
  })
  
  temp.dt.ratio.pars$spruce.a[i] <- summary(model)$parameters[1,1]
  temp.dt.ratio.pars$spruce.b[i] <- summary(model)$parameters[2,1]
  
  i=i+1
  
}


rm(b0, b1, b2, b3, b4, hour, i, js, month, x, y, colfunc)
rm(model, temp.dt, temp.rd, temp.rd.agg, djs.temp)

plot(temp.dt.ratio.pars$Date, temp.dt.ratio.pars$spruce.a, ylim=c(-1,2), col="red")
lines(temp.dt.ratio.pars$Date, temp.dt.ratio.pars$spruce.b, col="blue")

###

i=1
NROW(temp.dt.ratio.pars)
for (i in 1:NROW(temp.dt.ratio.pars)){
  #for (i in 1:500){
  
  pars.temp <- setNames(as.data.frame(matrix(ncol = 4, nrow = 13)),c("Depth","Js", "Ratio", "Jsnorm"))
  
  pars.temp[1,] <- c(0.2/8, temp.dt.ratio.pars$s.1.lg[i], temp.dt.ratio.pars$s.ra.1.lg[i], NA) # Extra
  pars.temp[2,] <- c(1/8, temp.dt.ratio.pars$s.1.lg[i], temp.dt.ratio.pars$s.ra.1.lg[i], NA)
  pars.temp[3,] <- c(3/8, temp.dt.ratio.pars$s.3.lg[i], temp.dt.ratio.pars$s.ra.3.lg[i], NA)
  pars.temp[4,] <- c(5/8, temp.dt.ratio.pars$s.5.lg[i], temp.dt.ratio.pars$s.ra.5.lg[i], NA)
  pars.temp[5,] <- c(7/8, temp.dt.ratio.pars$s.7.lg[i], temp.dt.ratio.pars$s.ra.7.lg[i], NA)
  pars.temp[6,] <- c(1/7, temp.dt.ratio.pars$s.1.md[i], temp.dt.ratio.pars$s.ra.1.md[i], NA)
  pars.temp[7,] <- c(0.2/7, temp.dt.ratio.pars$s.1.md[i], temp.dt.ratio.pars$s.ra.1.md[i], NA) # Extra
  pars.temp[8,] <- c(3/7, temp.dt.ratio.pars$s.3.md[i], temp.dt.ratio.pars$s.ra.3.md[i], NA)
  pars.temp[9,] <- c(5/7, temp.dt.ratio.pars$s.5.md[i], temp.dt.ratio.pars$s.ra.5.md[i], NA)
  pars.temp[10,] <- c(7/7, temp.dt.ratio.pars$s.7.md[i], temp.dt.ratio.pars$s.ra.7.md[i], NA)
  pars.temp[11,] <- c(1/3, temp.dt.ratio.pars$s.1.sm[i], temp.dt.ratio.pars$s.ra.1.sm[i], NA)
  pars.temp[12,] <- c(0.99,0,0, 0)
  pars.temp[13,] <- c(1,0,0, 0)
  
  pars.temp$Jsnorm <- pars.temp$Js/max(pars.temp$Js, na.rm=T)
  
  points(pars.temp$Depth, pars.temp$Js, col=colfunc(1000)[i])
  
  if(is.na(pars.temp$Ratio[1])){
    next
  }
  
  #pars.temp <- pars.temp[!is.na(pars.temp$Ratio),]
  
  
  ## Fit equation
  
  x = pars.temp$Depth
  y = pars.temp$Ratio
  
  b0 <- summary(model)$parameters[1,1]
  b1 <- summary(model)$parameters[2,1]
  
  try({
    
    model <- nlsLM(y ~ b0*exp(-0.5*((x-0.125)/b1)^2), 
                   control = nls.lm.control(maxiter= 200),data=pars.temp, 
                   start = list(b0=b0, b1=b1), na.action=NULL)
    #model <- nlsLM(y ~ b0*(1-exp(-b1*x))*(b0*exp(-b2*x)*x+(b3-b4*x)), 
    #              control = nls.lm.control(maxiter= 500),data=pars.temp, 
    #               start = list(b0=20, b1=1.5, b2=20, b3=0.1, b4=0.1), na.action=NULL)
    
    #start = list(b0=20, b1=1.5, b2=20, b3=0.1, b4=0.1), na.action=NULL)
    
    temp.dt.ratio.pars$spruce.a[i] <- summary(model)$parameters[1,1]
    temp.dt.ratio.pars$spruce.b[i] <- summary(model)$parameters[2,1]
    #temp.dt.ratio.pars$spruce.c[i] <- summary(model)$parameters[3,1]
    #temp.dt.ratio.pars$spruce.d[i] <- summary(model)$parameters[4,1]
    #temp.dt.ratio.pars$spruce.e[i] <- summary(model)$parameters[5,1]
    
    
    
  })
  
  
  #temp.dt <- setNames(as.data.frame(matrix(ncol = 4, nrow = 100)),c("Depth", "Js", "pred", "Jsa"))
  
  #temp.dt["Depth"] <- seq(0,1, length.out = 100)
  
  #temp.dt["Js"] <- rep(dt.ratio.pars$s.1.lg[i], NROW(temp.dt))
  
  b0 <- temp.dt.ratio.pars$spruce.a[i]
  b1 <- temp.dt.ratio.pars$spruce.b[i]
  #b2 <- temp.dt.ratio.pars$spruce.c[i]
  #b3 <- temp.dt.ratio.pars$spruce.d[i]
  #b4 <- temp.dt.ratio.pars$spruce.e[i]
  #hour <- temp.dt.ratio.pars$Hour[i]
  
  #temp.dt["pred"] <- b0*exp(-0.5*((temp.dt$Depth-0.125)/b1)^2)
  
  
  #temp.dt["pred"] <- b0*(1-exp(-b1*temp.dt$Depth))*(b0*exp(-b2*temp.dt$Depth)*temp.dt$Depth+(b3-b4*temp.dt$Depth))
  
  #temp.dt["Jsa"] <- temp.dt$pred*temp.dt$Js
  
  b0 <- NA
  b1 <- NA
  b2 <- NA
  b3 <- NA
  b4 <- NA
  
  #lines(temp.dt$Depth, temp.dt$Jsa, col=colfunc(1000)[i], lwd=1)
  #lines(temp.dt$Depth, temp.dt$pred, col=colfunc(1000)[i], lwd=1)
  
  
  i=i+1
  
}

plot(temp.dt.ratio.pars$Date, temp.dt.ratio.pars$spruce.a, ylim=c(-1,2), col="red")
lines(temp.dt.ratio.pars$Date, temp.dt.ratio.pars$spruce.b, col="blue")

rm(b0, b1, b2, b3, b4, hour, i, js, month, x, y, colfunc)
rm(model, temp.dt, temp.rd, temp.rd.agg, djs.temp)

plot(NULL, ylim=c(0.8,1.1))

lines(temp.dt.ratio.pars$Date, temp.dt.ratio.pars$spruce.a)

warnings()

temp.dt.ratio.pars["Month"] <- lubridate::month(temp.dt.ratio.pars$Date)
temp.dt.ratio.pars["Hour"] <- lubridate::hour(temp.dt.ratio.pars$Date)

# END OF SECTION ----

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#
#  SAP FLOW FROM JS AND RADIAL PROFILES
#
000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

# This section estimates sap flow using radial profiles for each species
# ADJUST - RADIAL PROFILE - MAIN LOOP [RADIAL PROFILE EQUATIONS]  ----------------------------------------------------------------

# Extract names of actual database

sap.set <- paste(names(wound.drift.fd[grep("Pine.*_02CM_.*|Spruce.*_02CM_.*", names(wound.drift.fd))]),"Q" , sep = ".")
list.fd.trees <- names(wound.drift.fd[grep("Pine.*_02CM_.*|Spruce.*_02CM_.*", names(wound.drift.fd))])

# Create empty data frame for new data

sapflow.temp <- setNames(as.data.frame(matrix(ncol = length(sap.set)+1, nrow = NROW(wound.drift.fd))),c( "Date", sap.set))

sapflow.temp$Date <- wound.drift.fd$Date


tree.dat.temp["TreeID2"] <- paste(substr(tree.dat.temp$TreeID,1,nchar(as.character(tree.dat.temp$TreeID))-7),"Q" , sep = ".")

str(tree.dat.temp$TreeID)

# Duplicate to compare effects of mean Js and integrating radial profiles
js.mean.temp <- sapflow.temp

names(wound.drift.fd)[2]
names(sap.set)

temp.dt <- setNames(as.data.frame(matrix(ncol = 10, nrow = 50)),c("Depth", "Ds", "Js", "Ratio", "Jsa", "As", "Asi", "Asi.perc", "Qi", "Qi.perc"))
temp.dt["Depth"] <- seq(0,1, length.out = 50)



i = 1
i=58
NCOL(wound.drift.fd)
for(i in 1:NCOL(wound.drift.fd)){
  
  tree.id <- list.fd.trees[i]
  fd = as.data.frame(wound.drift.fd[, tree.id])
  dbh <- tree.dat.temp$DBH[grep(sap.set[i], tree.dat.temp$TreeID2)]
  species <- tree.dat.temp$Species[grep(sap.set[i], tree.dat.temp$TreeID2)]
  
  Ds = tree.dat.temp$Ds[grep(sap.set[i], tree.dat.temp$TreeID2)]
  As = tree.dat.temp$As[grep(sap.set[i], tree.dat.temp$TreeID2)]
  
  #
  core.length <- ifelse(species=="Pine", -0.1824+0.8824*dbh, 0.7895+0.8342*dbh)
  temp.dt["Ds"] <- seq(0,Ds, length.out = 50)
  temp.dt["As"] <- (0.7854*(core.length-(temp.dt$Ds*2))^2)
  temp.dt["Asi"] <- ifelse(temp.dt$Ds==0,(0.7854*(core.length)^2)-temp.dt$As[1], dplyr::lag(temp.dt$As)-temp.dt$As)
  #temp.dt["Asi.perc"] <- temp.dt$Asi/sum(temp.dt$Asi)
  #plot(temp.dt$Ds[2:50], temp.dt$Asi[2:50])
  
  ii = 1
  #NROW(fd)
  for(ii in 1:NROW(fd)){
    #for(ii in 1:100){

    #####
    if(species=="Pine"){
      b0 = temp.dt.ratio.pars$pine.a[ii]
      b1 = temp.dt.ratio.pars$pine.b[ii]
    }else if(species=="Spruce"){
      b0 = temp.dt.ratio.pars$spruce.a[ii]
      b1 = temp.dt.ratio.pars$spruce.b[ii]
    }
    
    temp.dt["Js"] <- rep(fd[ii,], NROW(temp.dt))
    temp.dt["Ratio"] <- b0*exp(-0.5*((temp.dt$Depth-0.125)/b1)^2)
    temp.dt["Jsa"] <- temp.dt$Js*temp.dt$Ratio
    temp.dt["Qi"] <- (temp.dt$Jsa*temp.dt$Asi)/1000
    
    #lines(temp.dt$Depth, temp.dt$Jsa, col="blue", lwd=0.5)
    #points(1/8, fd[ii,], pch=19, cex=0.5, col="red")
    
    sapflow.temp[ii,tree.id] <- sum(temp.dt$Qi)
    js.mean.temp[ii,tree.id] <- mean(temp.dt$Jsa)
    
    ii=ii+1
  }
  i=i+1
}

names(sapflow.temp)
plot(sapflow.temp$Date, sapflow.temp[,tree.id])
lines(sapflow.temp$Date, sapflow.temp[,3], col="blue")

# Remove all variables related to fitting radial profiles

rm(As, fd, b0, b1, b2, b3, b4, core.length, dbh, Ds, i, ii, rad.set, species, tree.id, temp.dt)
rm(sap.set)
rm(b0_pine, b0_spruce, b1_pine, b1_spruce)
rm(col.del, list.fd.trees)

names(sapflow.temp)

#write.table(sapflow.temp, "sapflow.temp_march.08.2021.csv", sep=" ,", row.names=F)


000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#
#  AGGREGATE SAP FLOW
#
000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000


# AGGREGATE SAP FLOW BY DAY  ----------------------------------------------------------------


sapflow.temp["Year"] <-lubridate::year(sapflow.temp$Date)
sapflow.temp["Month"] <-lubridate::month(sapflow.temp$Date)
sapflow.temp["Week"] <-lubridate::week(sapflow.temp$Date)
sapflow.temp["DOEY"] <-lubridate::yday(sapflow.temp$Date)
sapflow.temp["Hour"] <-lubridate::hour(sapflow.temp$Date)
sapflow.temp["Date.d"] <-lubridate::date(sapflow.temp$Date)

#sapflow.temp["Drought.severe"] <- ifelse(sapflow.temp$Year==2018&sapflow.temp$DOY>179&sapflow.temp$DOY<204, "Drought",  
#                                        ifelse(sapflow.temp$Year==2018&sapflow.temp$DOY>204&sapflow.temp$DOY<230,  "Normal", "NA"))



sapflow.temp["DOY"] <-  ifelse(sapflow.temp$Year==2016,lubridate::yday(sapflow.temp$Date),
                               ifelse((sapflow.temp$Year==2017),lubridate::yday(sapflow.temp$Date)+366,
                                      ifelse((sapflow.temp$Year==2018),lubridate::yday(sapflow.temp$Date)+731,
                                             ifelse((sapflow.temp$Year==2019),lubridate::yday(sapflow.temp$Date)+1096,
                                                    lubridate::yday(sapflow.temp$Date)+1461))))


names(sapflow.temp)

#rm(all.nodes.sapflow)

as.data.frame(sap.set)

all.nodes.sapflow <- aggregate(cbind(Date, 
                                     Fd_NodeONE_2Spruce19_02CM_East_15.25DIAM_,
                                     Fd_NodeONE_3Spruce19_02CM_North_15.25DIAM_,
                                     Fd_NodeONE_4Spruce14_02CM_Random_7.15DIAM_,
                                     Fd_NodeONE_5Spruce19_02CM_South_15.25DIAM_,
                                     Fd_NodeONE_6Pine9_02CM_Random_7.15DIAM_,
                                     Fd_NodeONE_7Pine15_02CM_Random_7.15DIAM_,
                                     Fd_NodeONE_9Spruce10_02CM_Random_15.25DIAM_,
                                     Fd_NodeONE_14Pine12_02CM_Random_15.25DIAM_,
                                     Fd_NodeONE_18Spruce2_02CM_Random_7.15DIAM_,
                                     Fd_NodeONE_19Pine17_02CM_Random_25.30DIAM_,
                                     Fd_NodeONE_20Pine3_02CM_Random_15.25DIAM_,
                                     Fd_NodeONE_23Pine1_02CM_Random_25.30DIAM_,
                                     Fd_NodeONE_24Spruce5_02CM_Random_7.15DIAM_,
                                     Fd_NodeONE_26Pine6_02CM_Random_7.15DIAM_,
                                     Fd_NodeONE_29Pine7_02CM_Random_7.15DIAM_,
                                     Fd_NodeONE_30Spruce11_02CM_Random_25.30DIAM_,
                                     Fd_NodeONE_31Spruce4_02CM_Random_15.25DIAM_,
                                     Fd_NodeONE_35Pine16_02CM_East_15.25DIAM_,
                                     Fd_NodeONE_36Pine16_02CM_South_15.25DIAM_,
                                     Fd_NodeONE_37Pine16_02CM_North_15.25DIAM_,
                                     Fd_NodeONE_38Pine16_02CM_West_15.25DIAM_,
                                     Fd_NodeONE_42Spruce8_02CM_Random_15.25DIAM_,
                                     Fd_NodeONE_43Pine13_02CM_Random_15.25DIAM_,
                                     Fd_NodeONE_50Spruce18_02CM_Random_25.30DIAM_,
                                     Fd_NodeTWO_3Pine20_02CM_Random_25.30DIAM_,
                                     Fd_NodeTWO_8Spruce18_02CM_Random_15.25DIAM_,
                                     Fd_NodeTWO_9Spruce17_02CM_Random_7.15DIAM_,
                                     Fd_NodeTWO_13Spruce19_02CM_Random_15.25DIAM_,
                                     Fd_NodeTWO_15Pine16_02CM_Random_7.15DIAM_,
                                     Fd_NodeTWO_16Spruce3_02CM_Random_15.25DIAM_,
                                     Fd_NodeTWO_21Pine1_02CM_Random_15.25DIAM_,
                                     Fd_NodeTWO_22Spruce2_02CM_Random_7.15DIAM_,
                                     Fd_NodeTWO_23Pine4_02CM_Random_25.30DIAM_,
                                     Fd_NodeTWO_30Pine11_02CM_Random_15.25DIAM_,
                                     Fd_NodeTWO_33Pine9_02CM_Random_7.15DIAM_,
                                     Fd_NodeTWO_34Pine12_02CM_Random_25.30DIAM_,
                                     Fd_NodeTWO_35Spruce5_02CM_Random_7.15DIAM_,
                                     Fd_NodeTWO_36Pine7_02CM_South_25.30DIAM_,
                                     Fd_NodeTWO_37Spruce6_02CM_Random_15.25DIAM_,
                                     Fd_NodeTWO_38Pine7_02CM_East_25.30DIAM_,
                                     Fd_NodeTWO_39Pine7_02CM_North_25.30DIAM_,
                                     Fd_NodeTWO_40Pine7_02CM_West_25.30DIAM_,
                                     Fd_NodeTWO_41Spruce10_02CM_Random_7.15DIAM_,
                                     Fd_NodeTWO_42Pine14_02CM_Random_15.25DIAM_,
                                     Fd_NodeTWO_43Spruce15_02CM_Random_15.25DIAM_,
                                     Fd_NodeTWO_46Spruce8_02CM_South_15.25DIAM_,
                                     Fd_NodeTWO_47Spruce8_02CM_East_15.25DIAM_,
                                     Fd_NodeTWO_48Spruce8_02CM_West_15.25DIAM_,
                                     Fd_NodeTWO_49Spruce8_02CM_North_15.25DIAM_,
                                     Fd_NodeTHREE_2Spruce9_02CM_West_15.25DIAM_,
                                     Fd_NodeTHREE_3Spruce9_02CM_South_15.25DIAM_,
                                     Fd_NodeTHREE_6Spruce13_02CM_Random_25.30DIAM_,
                                     Fd_NodeTHREE_7Pine19_02CM_Random_25.30DIAM_,
                                     Fd_NodeTHREE_10Pine10_02CM_Random_15.25DIAM_,
                                     Fd_NodeTHREE_11Pine7_02CM_West_15.25DIAM_,
                                     Fd_NodeTHREE_12Spruce9_02CM_North_15.25DIAM_,
                                     Fd_NodeTHREE_13Spruce8_02CM_Random_25.30DIAM_,
                                     Fd_NodeTHREE_16Spruce11_02CM_Random_15.25DIAM_,
                                     Fd_NodeTHREE_17Pine2_02CM_Random_25.30DIAM_,
                                     Fd_NodeTHREE_20Pine7_02CM_East_15.25DIAM_,
                                     Fd_NodeTHREE_21Pine7_02CM_South_15.25DIAM_,
                                     Fd_NodeTHREE_25Pine4_02CM_Random_15.25DIAM_,
                                     Fd_NodeTHREE_26Pine3_02CM_Random_25.30DIAM_,
                                     Fd_NodeTHREE_29Pine6_02CM_Random_25.30DIAM_,
                                     Fd_NodeTHREE_32Spruce12_02CM_Random_25.30DIAM_,
                                     Fd_NodeTHREE_35Spruce14_02CM_Random_7.15DIAM_,
                                     Fd_NodeTHREE_36Spruce15_02CM_Random_7.15DIAM_,
                                     Fd_NodeTHREE_37Spruce17_02CM_Random_15.25DIAM_,
                                     Fd_NodeTHREE_40Spruce9_02CM_East_15.25DIAM_,
                                     Fd_NodeTHREE_41Pine1_02CM_Random_15.25DIAM_,
                                     Fd_NodeTHREE_44Pine20_02CM_Random_25.30DIAM_,
                                     Fd_NodeTHREE_48Spruce18_02CM_Random_25.30DIAM_,
                                     Fd_NodeTHREE_49Spruce16_02CM_Random_7.15DIAM_)~Date.d, sapflow.temp, FUN=sum, na.rm=T, na.action=NULL)


# Add DATE, YEAR, MONTH

all.nodes.sapflow["Date"] <- aggregate(cbind(Date)~Date.d, sapflow.temp, FUN=min, na.rm=T, na.action=NULL)[2]

all.nodes.sapflow$Date <- as.POSIXct(all.nodes.sapflow$Date, origin="1970/01/1", tz="UTC")
all.nodes.sapflow$Date <- as.Date(all.nodes.sapflow$Date)
names(all.nodes.sapflow)

# SAVE SAP FLOW DATABASE  ----------------------------------------------------------------

#write.table(all.nodes.sapflow, "Data/all.nodes.sapflow.csv", sep=",", row.names = FALSE)

names(all.notes.sapflow)
# STACK - SAP FLOW BY DAY  ----------------------------------------------------------------


# List of names to stack

list <- names(all.nodes.sapflow)[3:75]


# Create an empty data.frame 

stack.sapflow <- setNames(as.data.frame(matrix(ncol = 4, nrow = length(list)*NROW(all.nodes.sapflow) )),c("ID", "Sapflow", "Date", "DBH"))

temp.date <- seq(as.POSIXct("2016/06/17", tz="UTC"), as.POSIXct("2020/12/30 23:30:00", tz="UTC"), by = "day")

stack.sapflow["Date"] <- rep(temp.date, NROW(list))

rm(temp.date)

# Add ID name to empty data.frame


names(stack.sapflow)

stack.sapflow["ID"] <- NA

ii=1
j= 1
jj= NROW(all.nodes.sapflow) 
number = NROW(all.nodes.sapflow) 
for (i in 1:NROW(list)){
  
  name<-list[ii]
  
  stack.sapflow[j:jj,"ID"] <- rep(name, number)
  
  ii=ii+1
  j=j+number
  jj=jj+number
  
}

tail(stack.sapflow, 30)

rm(i, ii, j, jj)

# Add sap flow to each stack

ii=1
j= 1
jj= NROW(all.nodes.sapflow) 
number = NROW(all.nodes.sapflow) 
for (i in 1:NROW(list)){
  
  name<-list[ii]
  
  stack.sapflow[j:jj,"Sapflow"] <- all.nodes.sapflow[name]
  
  ii=ii+1
  j=j+number
  jj=jj+number
  
}

#Clean
rm(i, ii, j, jj,name, number)
rm(list)


# Add DBH

list <- names(all.nodes.sapflow)[3:75]

i=1
number = NROW(all.nodes.sapflow) 
for (i in 1:NROW(list)){
  
  name<-list[i]
  dbh <- tree.dat.temp$DBH[grep(name, tree.dat.temp$TreeID2)]
  
  
  stack.sapflow$DBH[grep(name,stack.sapflow$ID) ] <- rep(dbh, number)
  
  i=i+1
  
}

#Clean
rm(i,name, number, list)



col.del <-ifelse(stack.sapflow$ID=="Pine", "red", "blue")
plot(stack.sapflow$DBH, stack.sapflow$Sapflow, col=col.del, ylim=c(0,200))

plot(stack.sapflow$Sapflow[grep("Fd_NodeONE_23Pine1_02CM_Random_25.30DIAM_.radial.Q",stack.sapflow$ID)]) 
plot(stack.sapflow$Sapflow[grep("Fd_NodeONE_30Spruce11_02CM_Random_25.30DIAM_.radial.Q",stack.sapflow$ID)]) 


rm(col.del)




# REMOVE - SENSOR WITH UNUSUAL DATA [not running currently] ----------------------------------------------------------------

# Summary: I decided to wait until I fit the linear model, before deciding which data to delete

# DATA labeled as "suspicious"

# Remove high values of suspicious sensor

# Sensor from pine tree with excessively high Fd values in 2016

stack.sapflow$Sapflow[grep("Fd_NodeONE_23Pine1_02CM_Random_25.30DIAM_.radial.Q",stack.sapflow$ID)][1:199]  <- NA

# Sensor from sprue tree, that had signs of beetle damage

stack.sapflow$Sapflow[grep("Fd_NodeONE_30Spruce11_02CM_Random_25.30DIAM_.radial.Q",stack.sapflow$ID)][564:1292] <- NA

plot(stack.sapflow$Sapflow[grep("Fd_NodeONE_23Pine1_02CM_Random_25.30DIAM_.radial.Q",stack.sapflow$ID)])





# ADD - YEAR / MONTH / SPECIES / ETC   ----------------------------------------------------------------



# Add Month to database 

stack.sapflow["Month"] <- lubridate::month(stack.sapflow$Date)
stack.sapflow["Year"] <- lubridate::year(stack.sapflow$Date)
stack.sapflow["Week"] <- lubridate::week(stack.sapflow$Date)
stack.sapflow["DOY"] <- lubridate::yday(stack.sapflow$Date)

stack.sapflow["Species"] <- NA

stack.sapflow$Species[grep("Pine", stack.sapflow$ID)] <- "Pine"
stack.sapflow$Species[grep("Spruce", stack.sapflow$ID)] <- "Spruce"

plot(daily.data$Date, daily.data$spei)

names(daily.data)

# ADD - MAIN ENVIRONMENTAL VARIABLES TO daily.data   ----------------------------------------------------------------


#names(daily.data) <- c("Date", "S0.15", "Dz", "PARz", "Week", "Month", "Year", "Ta_1_1_1", "RH_1_1_1", "VPD_kPa", "S0.15.nomiss", "Ta_1_1_1.x", "RH_1_1_1.x", "VPD_kPa.x", "S0.15.nomiss.x", "P_1_1_1")


daily.data.temp <- subset(daily.data, daily.data$Date>"2016-01-16", select=c("Date","Ta_1_1_1", "RH_1_1_1", "VPD_kPa","S0.15", "S0.15.nomiss", "Dz","Dz.nomiss", "PARz","PARz.nomiss",  "spei", "P_1_1_1.no.mis" ))

plot(daily.data$Date, daily.data$spei)

daily.data.temp["Date"] <- as.Date(daily.data.temp$Date)

stack.sapflow$Date <- as.Date(stack.sapflow$Date)

tail(stack.sapflow$Date, 10)

tail(daily.data.temp$Date, 10)


stk.sf.met <-merge(daily.data.temp, stack.sapflow, by="Date", sort=T)


stk.sf.met$Sapflow[stk.sf.met$Sapflow==0] <- NA 


# Divide sap flow over VPD

stk.sf.met["Q.Dz"] <- stk.sf.met$Sapflow/stk.sf.met$Dz


stk.sf.met["drought.range"] <- ifelse(stk.sf.met$spei<(-1.5), "ASevere",
                                      ifelse(stk.sf.met$spei<(-0.5),"BMild", 
                                             ifelse(stk.sf.met$spei<0.5, "CNormal", "DWet")))




# ADD - QDz, DBH.range, ba, etc to  stk.sf.met.agg ----------------------------------------------------------------


names(stk.sf.met)
stk.sf.met["DBH.range"] <- ifelse(stk.sf.met$DBH<=15, "Small",
                                  ifelse(stk.sf.met$DBH>25, "Large", "Medium"))


stk.sf.met["ba"] <- ((0.7854*stk.sf.met$DBH^2)/10000)
stk.sf.met["Q.ba"] <- stk.sf.met$Sapflow/((0.7854*stk.sf.met$DBH^2)/10000)
stk.sf.met["As"] <- NA


stk.sf.met["Possition"] <- NA

stk.sf.met$Possition[grep("NodeONE", stk.sf.met$ID)] <- "Backslope"
stk.sf.met$Possition[grep("NodeTWO", stk.sf.met$ID)] <- "Shoulderslope"
stk.sf.met$Possition[grep("NodeTHREE", stk.sf.met$ID)] <- "Footslope"

stk.sf.met$As[grep("Pine", stk.sf.met$ID)] <- exp(par.estimates[7,2] + par.estimates[7,3] * log(stk.sf.met$DBH[grep("Pine", stk.sf.met$ID)]))
stk.sf.met$As[grep("Spruce", stk.sf.met$ID)] <- exp(par.estimates[8,2] + par.estimates[8,3] * log(stk.sf.met$DBH[grep("Spruce", stk.sf.met$ID)]))
plot(stk.sf.met$ba, stk.sf.met$As)

stk.sf.met["Q.As"] <- NA

stk.sf.met["Q.As"] <- stk.sf.met$Sapflow/stk.sf.met$As

stk.sf.met["Q.dzba"] <- NA

stk.sf.met["Q.dzba"] <- stk.sf.met$Q.Dz/stk.sf.met$ba

#plot(stk.sf.met$Date, stk.sf.met$Q.ba)
plot(stk.sf.met$Date, stk.sf.met$Q.Dz)

stk.sf.met["Repeat"] <- NA

stk.sf.met$Repeat[grep("Random|South", stk.sf.met$ID)] <- "Unique" 


# Add species when/if needed

stk.sf.met["Species"] <- NA

stk.sf.met$Species[grep("Pine", stk.sf.met$ID)] <- "Pine"
stk.sf.met$Species[grep("Spruce", stk.sf.met$ID)] <- "Spruce"


write.table(stk.sf.met, "Data/stk.sf.met.csv", sep=",", row.names = F)

# END OF SECTION ----


000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
#
#  KEY FIGURES 
#
000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000


000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000


# FIGURE::THREE:: MEAN DAILY SAP FLOW  ----------------------------------------------------------------

##### 

umecol <- c("firebrick1", "firebrick",
            "dodgerblue1","dodgerblue3", 
            "darkorchid1","darkorchid3",
            "cyan2","cyan3",
            "blue1", "blue3",
            "chartreuse2", "chartreuse4")

umecol_adj <- adjustcolor(c("firebrick1", "firebrick",
                            "dodgerblue1","dodgerblue3", 
                            "darkorchid1","darkorchid3",
                            "cyan2","cyan3",
                            "blue1", "blue3",
                            "chartreuse2", "chartreuse4"), alpha.f = 0.1)

#y.d <- c(2017, 2018, 2019)

##### 
split.screen(rbind(c(0.1, 0.75, 0.1, 0.98),
                   c(0.76, 0.98, 0.1, 0.98)))

screen(1)

par(mar = c(0.0, 0, 0, 0))

plot(NULL,  
     type="n", 
     xaxs = "r",
     xaxt="n",
     xlab = "",
     ylab = "",
     ylim = c(-1.5,20),
     xlim= c(1,89),
     col = "blue", pch=16 )

axis(side = 1, at = c(5,14,23, 37, 45, 54, 66, 75, 84),  labels = FALSE)
axis(side=1, at = c(5,14,23, 37, 45, 54, 66, 75, 84), tck=1, lwd = 0, 
     label= c("May", "Jul", "Sept", "May", "Jul", "Sept", "May", "Jul", "Sept"),line = 0, cex.axis=1)


mtext(c("2017", "2018", "2019"), side = 3, at=c(15, 47, 79) ,line = -1.5, cex = 1.2, col="deeppink")

#mtext("All trees", side = 3, at= 0.2, line = -2, cex = 1, col="black")

#abline(0,1, col="grey", lty=3, lwd=3)

#mtext(bquote("PAR"[z] ~'('*mu*'mol' *' '* 'm'^-2 *' '*'s'^-1*')'), side = 1, at=2, outer = FALSE,  cex = 1.1, line =2.5, col = "black")

mtext(expression('Mean sap flow (Q, Liters'*' '* 'd'^-1 * ')'), side = 2, outer = FALSE, cex = 1.1, line = 2, col = "black")

temp.mod <- subset(stk.sf.met, stk.sf.met$Year>2016&
                     stk.sf.met$ID!="Fd_NodeTHREE_7Pine19_02CM_Random_25.30DIAM_.Q"&
                     stk.sf.met$Sapflow>0&
                     stk.sf.met$P_1_1_1.no.mis<1&
                     stk.sf.met$Month>3& stk.sf.met$Month<11, 
                   select=c("Year", "Week", "Date", "PARz", "Dz", "Sapflow", "Month", "spei"))

temp.mod$Date <- lubridate::date(temp.mod$Date)

temp.mod.agg.date <- aggregate(cbind(Sapflow, PARz, Dz , Month, spei)~Date, data=temp.mod, FUN=mean, na.rm=T, na.action=NULL)

temp.mod.agg.date["Week"] <- lubridate::week(temp.mod.agg.date$Date)
temp.mod.agg.date["Year"] <- lubridate::year(temp.mod.agg.date$Date)


temp.mod.agg <- aggregate(cbind(Sapflow, Date, spei)~Week+Year, data=temp.mod.agg.date, FUN=mean, na.rm=T, na.action=NULL)


temp.mod.agg["stderr"] <- aggregate(cbind(Sapflow)~Week+Year, data=temp.mod, FUN=mad)["Sapflow"]

temp.mod.agg["Week"] <- seq(1:NROW(temp.mod.agg))


temp <- subset(temp.mod.agg, temp.mod.agg$Year=="2017")
#temp["Week"] <- seq(1,NROW(temp))

polygon(c(temp$Week, rev(temp$Week)),c((temp[,3]+temp[,6]), rev(temp[,3]-temp[,6])), border=FALSE, col=umecol_adj[3])

lines(temp$Week, temp$Sapflow, col=umecol[9],lty=1, lwd=3)

lines(temp$Week, temp$Sapflow+temp[,6], col=umecol[9],lty=3)
lines(temp$Week, temp$Sapflow-temp[,6], col=umecol[9],lty=3)

# SPEI polygon
col.spei <- ifelse(temp$spei<(-1.5), "red", 
                   ifelse(temp$spei>(-1.5)&temp$spei<(-0.5), "orange",
                          ifelse(temp$spei>(-0.5)&temp$spei<0.5, "darkgreen", "blue")))

points(temp$Week, rep(-0.5, NROW(temp)), pch=15, col = col.spei, cex=1.5)

# // Add SPEI legend

legend(15,-1,legend = c("SPEI:","Severe", "Mild", "Normal", "Wet"), pch=15,  ncol=5,col=c("White","Red", "Orange", "Darkgreen", "Blue"), cex=1, bty="n")

rm(col.spei)

temp <- subset(temp.mod.agg, temp.mod.agg$Year=="2018")

polygon(c(temp$Week, rev(temp$Week)),c((temp[,3]+temp[,6]), rev(temp[,3]-temp[,6])), border=FALSE, col=umecol_adj[3])

lines(temp$Week, temp$Sapflow, col=umecol[9],lty=1, lwd=3)

lines(temp$Week, temp$Sapflow+temp[,6], col=umecol[9],lty=3)
lines(temp$Week, temp$Sapflow-temp[,6], col=umecol[9],lty=3)

# SPEI polygon
col.spei <- ifelse(temp$spei<(-1.5), "red", 
                   ifelse(temp$spei>(-1.5)&temp$spei<(-0.5), "orange",
                          ifelse(temp$spei>(-0.5)&temp$spei<0.5, "darkgreen", "blue")))

points(temp$Week, rep(-0.5, NROW(temp)), pch=15, col = col.spei, cex=1.5)

rm(col.spei)

temp <- subset(temp.mod.agg, temp.mod.agg$Year=="2019")

polygon(c(temp$Week, rev(temp$Week)),c((temp[,3]+temp[,6]), rev(temp[,3]-temp[,6])), border=FALSE, col=umecol_adj[3])


lines(temp$Week, temp$Sapflow, col=umecol[9],lty=1, lwd=3)

lines(temp$Week, temp$Sapflow+temp[,6], col=umecol[9],lty=3)
lines(temp$Week, temp$Sapflow-temp[,6], col=umecol[9],lty=3)

# SPEI polygon
col.spei <- ifelse(temp$spei<(-1.5), "red", 
                   ifelse(temp$spei>(-1.5)&temp$spei<(-0.5), "orange",
                          ifelse(temp$spei>(-0.5)&temp$spei<0.5, "darkgreen", "blue")))

points(temp$Week, rep(-0.5, NROW(temp)), pch=15, col = col.spei, cex=1.5)

rm(col.spei)

rm(temp.mod)


screen(2)

par(mar = c(0.0, 0, 0, 0))

plot(NULL,  
     type="n", 
     xaxs = "r",
     xaxt="n",
     yaxt="n",
     xlab = "",
     ylab = "",
     ylim = c(-1.5,20),
     xlim = c(0,3), col = "blue", pch=16 )

temp.mod  <- subset(stk.sf.met, stk.sf.met$Year>2016&
                      stk.sf.met$ID!="Fd_NodeTHREE_7Pine19_02CM_Random_25.30DIAM_.Q"&
                      stk.sf.met$Sapflow>0&
                      stk.sf.met$P_1_1_1.no.mis<1&
                      stk.sf.met$Month>3& stk.sf.met$Month<11, 
                    select=c("Year", "Week", "Date", "PARz", "Dz", "Sapflow", "Month", "spei"))

temp.mod.agg.date <- aggregate(cbind(Sapflow, PARz, Dz , Month, spei)~Date, data=temp.mod, FUN=mean, na.rm=T, na.action=NULL)
temp.mod.agg.date["Year"] <- lubridate::year(temp.mod.agg.date$Date)

# Extract statistics 

sf <- temp.mod.agg.date$Sapflow
yr <- temp.mod.agg.date$Year

sf.aov <- aov(yr~sf, data=temp.mod.agg.date)

sf.aov <- summary(sf.aov)

p.value <- ifelse(sf.aov[[1]]$`Pr(>F)`[1]<0.0001, 0.0001, sf.aov[[1]]$`Pr(>F)`[1])

f.value <- sf.aov[[1]]$`F value`[1]

# // Add legend & insert predicted parameters for each site

expr <- vector("expression", 2)
expr[[1]] <- bquote(bold("ANOVA"))
expr[[2]] <- bquote(~"F="~ .(round(f.value,2)))
expr[[3]] <- bquote(  ~ "Pvalue="~"<" ~.(format(round(p.value, 5), scientific=F)))

legend(0,20,legend = expr, pch=NA,  ncol=1,col="white", cex=0.8, bty="n")

rm(sf, yr, sf.aov, expr, p.value, f.value)
rm(expr)

b2017 <- data.frame(subset(temp.mod.agg.date, temp.mod.agg.date$Year==2017, select=c("Sapflow")))
b2018 <- data.frame(subset(temp.mod.agg.date, temp.mod.agg.date$Year==2018, select=c("Sapflow")))
b2019 <- data.frame(subset(temp.mod.agg.date, temp.mod.agg.date$Year==2019, select=c("Sapflow")))

par(new = TRUE)

par(mar = c(0, 0, 0, 0))

boxplot(c(b2017, b2018, b2019), col=c(umecol[9]),bty="n",frame=F,
        xaxt="n", 
        yaxt="n",
        ylim = c(-1.5,20))

axis(side = 1, at = c(1,2,3),  labels = FALSE)
axis(side=1, at = c(1,2,3), tck=1, lwd = 0, 
     label= c("2017","2018", "2019"),line = 0, cex.axis=1)


# Just for reference, do not use
#axis(side=2, at = pretty(range(-0.05,15)))

# Add means to box plots

expr <- vector("expression", 2)

expr[[1]] <- bquote(.(round(mean(na.omit(b2017[,])),2)))
expr[[2]] <- bquote(.(round(mean(na.omit(b2018[,])),2)))
expr[[3]] <- bquote(.(round(mean(na.omit(b2019[,])),2)))

loc.x <- c(bquote(.(round(mean(na.omit(b2017[,])),2))), 
           bquote(.(round(mean(na.omit(b2018[,])),2))),
           bquote(.(round(mean(na.omit(b2019[,])),2))))

text(c(1, 2, 3), loc.x, expr,col="white", cex=1.2)

rm(b2017, b2018, b2019,  expr, loc.x)

rm(pchplot, legend)

rm(temp.mod, temp.mod.agg)

close.screen(all.screens = TRUE)
# END OF SECTION ----

# STATISTICAL ANAYSES - FIGURE 3
##### 
temp.barplot <- subset(stk.sf.met, stk.sf.met$Year>2016&
                         stk.sf.met$ID!="Fd_NodeTHREE_7Pine19_02CM_Random_25.30DIAM_.Q"&
                         
                         stk.sf.met$Month>4&stk.sf.met$Month<10&
                         stk.sf.met$Sapflow>0&
                         stk.sf.met$Dz.nomiss>0.1&stk.sf.met$P_1_1_1.no.mis<1,
                       select=c(Date, Species, Dz, Q.ba, spei,Sapflow,Q.Dz, drought.range))
temp.barplot["Year"] <- lubridate::year(temp.barplot$Date)

temp.barplot.ag.date <- aggregate(cbind(Q.Dz, Dz, Sapflow, Year)~Date, data=temp.barplot, FUN=mean, na.rm=T, na.action=NULL)

tapply(temp.barplot.ag.date$Q.Dz, temp.barplot.ag.date$Year, FUN=mean, na.rm=T) 

par(mar = c(2, 5, 0.5, 0.5))

plot(temp.barplot.ag.date$Date, temp.barplot.ag.date$Q.Dz)
# Extract statistics 

sf <- temp.barplot.ag.date$Q.Dz
yr <- temp.barplot.ag.date$Year

aov <- aov(yr~sf, data=temp.barplot.ag.date)

aov <- summary(aov)

rm(temp.barplot,temp.barplot.ag.date)
rm(sf, yr, aov)

# END OF SECTION ----

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

# FIGURE::FOUR:: MEAN SAP FLOW BY TREE SIZE & DROUGHT RANGE ----------------------------------------------------------------

##### 

umecol <- c("firebrick1", "firebrick",
            "dodgerblue1","dodgerblue3", 
            "darkorchid1","darkorchid3",
            "cyan2","cyan3",
            "blue1", "blue3",
            "chartreuse2", "chartreuse4")

st.err <-function(x){sd(x)/sqrt(length(x))}

#ylim.qdz <- c(0,250)
y.d <- c(2018)

#plot(stk.sf.met$Date[stk.sf.met$ID=="Fd_NodeTHREE_7Pine19_02CM_Random_25.30DIAM_.Q"], stk.sf.met$Sapflow[stk.sf.met$ID=="Fd_NodeTHREE_7Pine19_02CM_Random_25.30DIAM_.Q"])
ylim.q <- c(0,9)

size.id <- c("Large", "Medium", "Small")

col.qdz <- c(umecol[1], umecol[3])

par(mfrow =c(1,3)) 

i=1
for (i in 1:3) {
  
  temp.barplot <- subset(stk.sf.met, stk.sf.met$Year%in%y.d&
                           stk.sf.met$Repeat=="Unique"&
                           stk.sf.met$Month>4&stk.sf.met$Month<10&
                           stk.sf.met$Sapflow>0&
                           stk.sf.met$Dz>0.1&stk.sf.met$P_1_1_1.no.mis<1&
                           stk.sf.met$DBH.range==size.id[i], 
                         select=c(Date, Species, Dz, spei,Sapflow,Q.Dz, drought.range))
  temp.barplot$Date <- lubridate::date(temp.barplot$Date)
  
  temp.barplot <- aggregate(cbind(Dz, Sapflow,Q.Dz)~Date+Species+drought.range, data=temp.barplot, FUN=mean, na.rm=T, na.action=NULL)
  
  temp.barplot["ID"] <- paste(temp.barplot$drought.range, temp.barplot$Species, sep=" ")
  
  
  temp.barplot.mn<- as.vector(aggregate(cbind(Sapflow, Q.Dz)~ID, data=temp.barplot, FUN=mean, na.rm=T, na.action=NULL))
  
  
  ylim.qdz <- c(0,max(temp.barplot.mn$Q.Dz)*1.3)
  
  par(mar = c(4, 5, 0.5, 0.5))
  
  plot(NULL,  
       type="n", 
       xaxs = "r",
       xaxt = "n",
       xlab = " ",
       ylab = " ",
       xlim = c(0,6),
       ylim = ylim.qdz, col = "blue", pch=16 )
  
  #axis(side = 2, at = c(0, 5, 10, 15, 20, 25, 30),  labels = FALSE)
  #axis(side=2, at = c(0, 5, 10, 15, 20, 25, 30), tck=0.1, lwd = 0, label= c(0, 5, 10, 15, 20, 25, 30),line = 0, cex.axis=1)
  
  mtext(size.id[i], side = 3,  line = -2, cex = 1, col="deeppink")
  
  mtext(expression('Q'[Dz] ~'(Liters'*' '*'day'^-1 *' '*'kPa'^-1 * ')'), side=2, outer=F, cex=0.8, line=2, col="black")
  mtext(expression('Drought range'), side=1, outer=F, cex=0.8, line=3, col="black")
  
  
  if(i==1){
    
    legend("topleft",  legend = c("Pine", "Spruce"), fill = c(umecol[1], umecol[3]),
           density = c(200, 200, 200), angle = c( 0, 0,0), bty="n", ncol=1, cex=1)
    
  }
  
  
  par(new=TRUE)
  
  barplot( Q.Dz~ ID, temp.barplot.mn,
           #      xlim=c(0,6),
           xaxt= "n",
           yaxt= "n",
           xlab= " ",
           ylab= " ",
           ylim = ylim.qdz,
           space=c(0.1,0.1,0.5,0.1,0.5,0.1,0.5,0.1), 
           col=col.qdz)
  
  temp.barplot.sterr <- as.vector(aggregate(cbind(Sapflow, Q.Dz)~ID, data=temp.barplot, FUN=st.err))[2:3]
  
  arrows(c(0.6, 1.7, 3.2, 4.3,5.8, 6.9, 8.4, 9.5), temp.barplot.mn$Q.Dz-temp.barplot.sterr$Q.Dz, 
         c(0.6, 1.7, 3.2, 4.3,5.8, 6.9, 8.4, 9.5), temp.barplot.mn$Q.Dz+temp.barplot.sterr$Q.Dz, length=0.05, angle=90, code=3,lwd= 2, col="grey")
  
  axis(side=1, at = c(1.3,3.7,6.1, 8.5), tck=0.1, lwd = 0, label= c("Severe", "Mild", "Normal", "Wet"),line = 0, cex.axis=1)
  
  
  # Add average change
  expr <- vector("expression", 2)
  
  expr[[1]] <- bquote(.(round((temp.barplot.mn[1,3]-temp.barplot.mn[5,3])/temp.barplot.mn[5,3]*100,0)))
  expr[[2]] <- bquote(.(round((temp.barplot.mn[2,3]-temp.barplot.mn[6,3])/temp.barplot.mn[6,3]*100,0)))
  
  
  loc.x <- c(temp.barplot.mn[1,3]+temp.barplot.mn[1,3]*0.15,temp.barplot.mn[2,3]+temp.barplot.mn[2,3]*0.15)
  
  text(c(0.6, 1.7), loc.x, expr,col="blue", cex=1)
  
  # Arrows for percent change
  
  arrows(c(5.5, 6.7), temp.barplot.mn$Q.Dz[5:6], 
         c(5.5, 6.7), temp.barplot.mn$Q.Dz[5:6]+temp.barplot.mn$Q.Dz[5:6]*0.15, length=0, angle=0, code=2,lwd= 2, lty=1,col=col.qdz)
  
  arrows(c(5.5, 6.7), temp.barplot.mn$Q.Dz[5:6]+temp.barplot.mn$Q.Dz[5:6]*0.15, 
         c(0.7, 1.9), temp.barplot.mn$Q.Dz[5:6]+temp.barplot.mn$Q.Dz[5:6]*0.15, length=0, angle=0, code=2,lwd= 2, lty=1,col=col.qdz)
  
  arrows(c(0.7, 1.9), temp.barplot.mn$Q.Dz[5:6]+temp.barplot.mn$Q.Dz[5:6]*0.15, 
         c(0.7, 1.9), loc.x+loc.x*0.2, length=0.1, angle=45, code=2,lwd= 2, lty=1,col=col.qdz)
  
  
  ### Extract statistics 
  
  sf <- as.numeric(temp.barplot$Q.Dz)
  ID <- temp.barplot$ID
  
  sf.aov <- aov(sf~ID)
  
  sf.aov <- summary(sf.aov)
  
  p.value <- ifelse(sf.aov[[1]]$`Pr(>F)`[1]<0.001, 0.0001, sf.aov[[1]]$`Pr(>F)`[1])
  
  f.value <- sf.aov[[1]]$`F value`[1]
  
  # // Add legend & insert predicted parameters for each site
  
  expr <- vector("expression", 2)
  expr[[1]] <- bquote(italic("ANOVA"))
  expr[[2]] <- bquote(~"F="~ .(round(f.value,2)))
  expr[[3]] <- bquote(  ~ "Pvalue="~"<" ~.(format(p.value, scientific=F)))
  
  #legend(3.7,max(temp.barplot.mn$Q.Dz)*1.2,legend = expr, pch=NA,  ncol=1,col="white", cex=1, bty="n")
  
  rm(sf, ID, sf.aov, p.value, f.value)
  
  rm(expr, loc.x)
  
  rm(temp.barplot, temp.barplot.mn)
  
  rm(temp.barplot.sterr)
  
  
  i=1+1
  
}

# END OF SECTION ----


# STATISTICAL ANAYSES - FIGURE 4
##### 
y.d <- c(2018) 
i=1
for (i in 1:3) {
  
  temp.barplot <- subset(stk.sf.met, stk.sf.met$Year%in%y.d&
                           stk.sf.met$Repeat=="Unique"&
                           stk.sf.met$Month>4&stk.sf.met$Month<10&
                           stk.sf.met$Sapflow>0&
                           stk.sf.met$Dz>0.1&stk.sf.met$P_1_1_1.no.mis<1&
                           stk.sf.met$DBH.range==size.id[i], 
                         select=c(Date, Species, Dz, spei,Sapflow,Q.Dz, drought.range))
  temp.barplot$Date <- lubridate::date(temp.barplot$Date)
  
  temp.barplot <- aggregate(cbind(Dz, Sapflow,Q.Dz)~Date+Species+drought.range, data=temp.barplot, FUN=mean, na.rm=T, na.action=NULL)
  
  temp.barplot["ID"] <- paste(temp.barplot$drought.range, temp.barplot$Species, sep=" ")
  
  
  temp.barplot.mn<- as.vector(aggregate(cbind(Sapflow, Q.Dz)~drought.range, data=temp.barplot, FUN=mean, na.rm=T, na.action=NULL))
  
  expr <- vector("expression", 1)
  
  expr[[1]] <- bquote(.(round((temp.barplot.mn[1,3]-temp.barplot.mn[3,3])/temp.barplot.mn[3,3]*100,0)))
  expr[[2]] <- bquote(.(round((temp.barplot.mn[2,3]-temp.barplot.mn[6,3])/temp.barplot.mn[6,3]*100,0)))
  
  print(expr)
  
  
  i=i+1
}


### Mean sap flow by year 

temp.mod <- subset(stk.sf.met, stk.sf.met$Year>2016&
                     stk.sf.met$Sapflow>0&
                     stk.sf.met$Repeat=="Unique"&
                     stk.sf.met$P_1_1_1.no.mis<1&
                     stk.sf.met$Dz>0.1&
                     stk.sf.met$Month>3&stk.sf.met$Month<11, 
                   select=c("Year", "Week", "Date", "PARz", "Dz", "Sapflow", "Month"))

b2017 <- data.frame(subset(temp.mod, temp.mod$Year==2017&temp.mod$Dz>0.1, select=c("Sapflow")))
b2018 <- data.frame(subset(temp.mod, temp.mod$Year==2018&temp.mod$Dz>0.1, select=c("Sapflow")))
b2019 <- data.frame(subset(temp.mod, temp.mod$Year==2019&temp.mod$Dz>0.1, select=c("Sapflow")))

mean(b2017, na.rm=T)

mean(na.omit(b2017[,]))
mean(na.omit(b2018[,]))
mean(na.omit(b2019[,]))

mean(temp.mod$Sapflow[temp.mod$Year==2017], na.rm=T)
mean(temp.mod$Sapflow[temp.mod$Year==2018], na.rm=T)
mean(temp.mod$Sapflow[temp.mod$Year==2019], na.rm=T)

rm(temp.mod, b2017, b2018, b2019)

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

# END OF SECTION ----

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

# FIGURE::FIVE:: BOREAL DROUGHT - MEAN SAP FLOW BY TOPOGRAPHIC POSSITION & DROUGHT RANGE ----------------------------------------------------------------

##### 

umecol <- c("firebrick1", "firebrick",
            "dodgerblue1","dodgerblue3", 
            "darkorchid1","darkorchid3",
            "cyan2","cyan3",
            "blue1", "blue3",
            "chartreuse2", "chartreuse4")

st.err <-function(x){sd(x)/sqrt(length(x))}

#ylim.qdz <- c(0,250)
y.d <- c(2018)

ylim.q <- c(0,9)

possition.id <- c("Shoulderslope", "Backslope", "Footslope")

col.qdz <- c(umecol[1], umecol[3])

par(mfrow =c(1,3)) 
i=1
for (i in 1:3) {
  
  temp.barplot <- subset(stk.sf.met, stk.sf.met$Year%in%y.d&
                           stk.sf.met$Repeat=="Unique"&
                           stk.sf.met$Month>3&stk.sf.met$Month<11&
                           stk.sf.met$Sapflow>0&
                           stk.sf.met$Dz>0.1&stk.sf.met$P_1_1_1.no.mis<1&
                           stk.sf.met$Possition==possition.id[i], 
                         select=c(Date, Species, Dz, Ta_1_1_1, spei,Sapflow,Q.Dz, drought.range))
  temp.barplot$Date <- lubridate::date(temp.barplot$Date)
  
  #temp.barplot["Q.Dz"] <- temp.barplot$Q.Dz*(115.8+0.423*temp.barplot$Ta_1_1_1)
  
  temp.barplot <- aggregate(cbind(Dz, Sapflow,Q.Dz)~Date+Species+drought.range, data=temp.barplot, FUN=mean, na.rm=T, na.action=NULL)
  
  temp.barplot["ID"] <- paste(temp.barplot$drought.range, temp.barplot$Species, sep=" ")
  
  temp.barplot.mn<- as.vector(aggregate(cbind(Sapflow, Q.Dz)~ID, data=temp.barplot, FUN=mean, na.rm=T, na.action=NULL))
  
  ylim.qdz <- c(0,max(temp.barplot.mn$Q.Dz)*1.3)
  
  par(mar = c(4, 5, 0.5, 0.5))
  
  plot(NULL,  
       type="n", 
       xaxs = "r",
       xaxt = "n",
       xlab = " ",
       ylab = " ",
       xlim = c(0,6),
       ylim = ylim.qdz, col = "blue", pch=16 )
  
  mtext(possition.id[i], side = 3,  line = -2, cex = 1, col="deeppink")
  
  mtext(expression('Q'[Dz] ~'(Liters'*' '*'day'^-1 *' '*'kPa'^-1 * ')'), side=2, outer=F, cex=0.8, line=2, col="black")
  
  if(i==1){
    
    legend("topleft",  legend = c("Pine", "Spruce"), fill = c(umecol[1], umecol[3]),
           density = c(200, 200, 200), angle = c( 0, 0,0), bty="n", ncol=1, cex=1)
    
  }
  
  par(new=TRUE)
  
  barplot( Q.Dz~ ID, temp.barplot.mn,
           #      xlim=c(0,6),
           xaxt= "n",
           yaxt= "n",
           xlab= " ",
           ylab= " ",
           ylim = ylim.qdz,
           space=c(0.1,0.1,0.5,0.1,0.5,0.1,0.5,0.1), 
           col=col.qdz)
  
  temp.barplot.sterr <- as.vector(aggregate(cbind(Sapflow, Q.Dz)~ID, data=temp.barplot, FUN=st.err))[2:3]
  
  arrows(c(0.6, 1.7, 3.2, 4.3,5.8, 6.9, 8.4, 9.5), temp.barplot.mn$Q.Dz-temp.barplot.sterr$Q.Dz, 
         c(0.6, 1.7, 3.2, 4.3,5.8, 6.9, 8.4, 9.5), temp.barplot.mn$Q.Dz+temp.barplot.sterr$Q.Dz, length=0.05, angle=90, code=3,lwd= 2, col="grey")
  
  axis(side=1, at = c(1.3,3.7,6.1, 8.5), tck=0.1, lwd = 0, label= c("Severe", "Mild", "Normal", "Wet"),line = 0, cex.axis=1)
  
  
  # Add average change
  expr <- vector("expression", 2)
  
  expr[[1]] <- bquote(.(round((temp.barplot.mn[1,3]-temp.barplot.mn[5,3])/temp.barplot.mn[5,3]*100,0)))
  expr[[2]] <- bquote(.(round((temp.barplot.mn[2,3]-temp.barplot.mn[6,3])/temp.barplot.mn[6,3]*100,0)))
  
  
  loc.x <- c(temp.barplot.mn[1,3]+temp.barplot.mn[1,3]*0.15,temp.barplot.mn[2,3]+temp.barplot.mn[2,3]*0.15)
  
  text(c(0.6, 1.7), loc.x, expr,col="blue", cex=1)
  
  # Arrows for percent change
  
  arrows(c(5.5, 6.7), temp.barplot.mn$Q.Dz[5:6], 
         c(5.5, 6.7), temp.barplot.mn$Q.Dz[5:6]+temp.barplot.mn$Q.Dz[5:6]*0.15, length=0, angle=0, code=2,lwd= 2, lty=1,col=col.qdz)
  
  arrows(c(5.5, 6.7), temp.barplot.mn$Q.Dz[5:6]+temp.barplot.mn$Q.Dz[5:6]*0.15, 
         c(0.7, 1.9), temp.barplot.mn$Q.Dz[5:6]+temp.barplot.mn$Q.Dz[5:6]*0.15, length=0, angle=0, code=2,lwd= 2, lty=1,col=col.qdz)
  
  arrows(c(0.7, 1.9), temp.barplot.mn$Q.Dz[5:6]+temp.barplot.mn$Q.Dz[5:6]*0.15, 
         c(0.7, 1.9), loc.x+loc.x*0.2, length=0.1, angle=45, code=2,lwd= 2, lty=1,col=col.qdz)
  
  
  ### Extract statistics 
  
  sf <- as.numeric(temp.barplot$Q.Dz)
  ID <- temp.barplot$ID
  
  #
  
  #wilcox.test(temp.barplot$Q.Dz~temp.barplot$ID, data=temp.barplot)
  
  #
  
  sf.aov <- aov(sf~ID)
  
  sf.aov <- summary(sf.aov)
  
  p.value <- ifelse(sf.aov[[1]]$`Pr(>F)`[1]<0.001, 0.0001, sf.aov[[1]]$`Pr(>F)`[1])
  
  f.value <- sf.aov[[1]]$`F value`[1]
  
  # // Add legend & insert predicted parameters for each site
  
  expr <- vector("expression", 2)
  expr[[1]] <- bquote(italic("ANOVA"))
  expr[[2]] <- bquote(~"F="~ .(round(f.value,2)))
  expr[[3]] <- bquote(  ~ "Pvalue="~"<" ~.(format(p.value, scientific=F)))
  
  #legend(3.7,max(temp.barplot.mn$Q.Dz)*1.2,legend = expr, pch=NA,  ncol=1,col="white", cex=1, bty="n")
  
  rm(sf, ID, sf.aov, p.value, f.value)
  
  rm(expr, loc.x)
  
  rm(temp.barplot, temp.barplot.mn)
  
  rm(temp.barplot.sterr)
  
  
  i=1+1
  
}

##### 

rm(i, possition.id, size.id, co.qdz, labels, ylim.q, ylim.qdz, y.d)
rm(col.qdz, legend, umecol, umecol_adj, st.err)
# END OF SECTION ----

# STATISTICAL ANAYSES - FIGURE 5

##### 

possition.id <- c("Shoulderslope","Backslope",  "Footslope")


i=1
for (i in 1:3) {
  
  temp.barplot <- subset(stk.sf.met, stk.sf.met$Year%in%y.d&
                           stk.sf.met$Repeat=="Unique"&
                           stk.sf.met$Month>4&stk.sf.met$Month<10&
                           stk.sf.met$Sapflow>0&
                           stk.sf.met$Dz>0.1&stk.sf.met$P_1_1_1.no.mis<1&
                           stk.sf.met$Possition==possition.id[i], 
                         select=c(Date, Species, Dz, spei,Sapflow,Q.Dz, drought.range))
  temp.barplot$Date <- lubridate::date(temp.barplot$Date)
  
  temp.barplot <- aggregate(cbind(Dz, Sapflow,Q.Dz)~Date+Species+drought.range, data=temp.barplot, FUN=mean, na.rm=T, na.action=NULL)
  
  temp.barplot["ID"] <- paste(temp.barplot$drought.range, temp.barplot$Species, sep=" ")
  
  
  temp.barplot.mn<- as.vector(aggregate(cbind(Sapflow, Q.Dz)~drought.range, data=temp.barplot, FUN=mean, na.rm=T, na.action=NULL))
  
  print((temp.barplot.mn[1,3]-temp.barplot.mn[3,3])/temp.barplot.mn[3,3])
  
  print((temp.barplot.mn[1,3]-temp.barplot.mn[5,3])/temp.barplot.mn[5,3])
  print((temp.barplot.mn[2,3]-temp.barplot.mn[6,3])/temp.barplot.mn[6,3])
  
}

# END OF SECTION ----

000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000

