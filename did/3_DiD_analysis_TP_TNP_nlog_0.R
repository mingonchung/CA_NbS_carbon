#######
## CA C DiD analysis
## (1) treatment effects between PAs and non-PAs
## 10/14/2023

rm(list=ls())

## number
num <- commandArgs(trailingOnly = TRUE)
num <- as.numeric(num)

## library
library(data.table)
library(did)
library(etwfe)
library(dplyr)
library(ggplot2)

## setwd
setwd("/projects/mich9173/CA_carbon/output/DiD_input")
input <- fread("3_DiD_data_TP_TNP_NN.csv.gz")

## is.na = 0
input$year.mngint[is.na(input$year.mngint)] <- 0

## C log transformation
# CFlux_GPP
input$lCFlux_GPP <- log1p(input$CFlux_GPP)
# CONUS
input$lCONUS <- log1p(input$CONUS)
# CStocks_Live
input$lCStocks_Live <- log1p(input$CStocks_Live)
# LEMMA
input$lLEMMA <- log1p(input$LEMMA)

## subset
ii <- num
c.list <- c("CFlux_GPP", "CONUS", "CStocks_Live", "LEMMA")
lc.list <- c("lCFlux_GPP", "lCONUS", "lCStocks_Live", "lLEMMA")
#c.col.num <- c(4,5,6,7)

### select log or non-log #!!!!!!!!!
## subset by c.list
col.list <- c("V1","layer","year",c.list[ii],"US_L3CODE", "year.mngint", "land.mngint", "treat", "subclass", "ppt_00_02", "tmean_00_02", "DEM_slope", "DEM_elevation", "DEM_aspect", "popden_00_02", "citytrvtime_11")

input$land.mngint[is.na(input$land.mngint)] <- 0
input.sub <- input %>% select(all_of(col.list))
input.sub <- input.sub[complete.cases(input.sub),]

if(ii == 1 | ii == 3) {
  input.sub <- subset(input.sub, year > 1989) #1989; 1994
  
}

rm(list=c("input"))
gc()
set.seed(160617)

for(jj in c(1,2,3)) {
  
  input.sub2 <- subset(input.sub, land.mngint == jj)

  print(paste0(c.list[[ii]]," DiD started"))
  print(paste0("treatment severity ", jj))
  print("input data processed")
  
###################
### (1) DiD analyses w/ did package
# return a class MP object
# estimates of group-time average treatment effects for all groups in all time periods
# estimate group-time average treatment effects without covariates
c.attgt.0 <- att_gt(yname = c.list[[ii]],
                   gname = "year.mngint",
                   idname = "V1",
                   tname = "year",
                   xformla = ~ 1, # without covariates
                   data = subset(input.sub2, treat==0),
                   control_group = "notyettreated")

  c.attgt.1 <- att_gt(yname = c.list[[ii]],
                      gname = "year.mngint",
                      idname = "V1",
                      tname = "year",
                      xformla = ~ 1, # without covariates
                      data = subset(input.sub2, treat==1),
                      control_group = "notyettreated")

print("DiD analysis done")

########
## (1-1) Overall effects of participating in the treatment
# "Overall ATT": we estimate that increasing the minimum wage decreased teen employment by 3.1% and the effect is marginally statistically significant
# type="simple": tends to overweight the effect of early-treated groups simply because we observe more of them during post-treatment periods
# type="dynamic": aggregate group-time effects into an event study plot; averaged into average treatment effects at different lengths of exposure to the treatment
# type="group": aggregate group-time average treatment effects into group-specific average treatment effects; The Overall ATT averages - this parameter is a leading choice as an overall summary effect of participating in the treatment
# type="calendar": aggregations across different time periods; y-axis - the average effect of participating in the treatment in a particular time period for all groups that participated in the treatment in that time period
group_effects.0 <- aggte(c.attgt.0, type = "group", na.rm=TRUE) # simple; dynamic; group; calendar;
summary(group_effects.0)

group_effects.1 <- aggte(c.attgt.1, type = "group", na.rm=TRUE) # simple; dynamic; group; calendar;
summary(group_effects.0)

print("DiD group summary done")

########
## (1-2) Event Studies
# aggregate the group-time average treatment effects
# aggregate the group-time average treatment effects into a small number of parameters
c.dyn.0 <- aggte(c.attgt.0, type = "dynamic", na.rm=TRUE)
c.dyn.1 <- aggte(c.attgt.1, type = "dynamic", na.rm=TRUE)

# "event time": when they first participate in the treatment
# event time=0: the "on impact" effect
# event time=-1: the effect in the period before a unite becomes treated
# checking that this is equal to 0 is potentiall useful as a pre-test
c.dyn.sum.0 <- summary(c.dyn.0)
c.dyn.sum.1 <- summary(c.dyn.1)

# one fails to reject parallel trends in pre-treatment periodsand it looks like somewhat negative effects of the minimum wage on youth employment.
jpeg(paste0("./fig/3_",jj,"_did_dynamic_",c.list[ii],"_0_nlog.jpeg"), width=7, height=4, units="in", res=300)
plot.did.0 <- ggdid(c.dyn.0)
print(plot.did.0)
dev.off()

jpeg(paste0("./fig/3_",jj,"_did_dynamic_",c.list[ii],"_1_nlog.jpeg"), width=7, height=4, units="in", res=300)
plot.did.1 <- ggdid(c.dyn.1)
print(plot.did.1)
dev.off()

print("DiD dynamic plot done")

############
## save output csv
# did all effects
did.df.all.0 <- broom::tidy(group_effects.0)
write.csv(did.df.all.0, paste0("./csv/3_",jj,"_did_all_TP_TNP_",num,"_0_nlog.csv"), row.names=F)

did.df.all.1 <- broom::tidy(group_effects.1)
write.csv(did.df.all.1, paste0("./csv/3_",jj,"_did_all_TP_TNP_",num,"_1_nlog.csv"), row.names=F)

# did evenet effects
did.df.event.0 <- broom::tidy(c.dyn.0)
write.csv(did.df.event.0, paste0("./csv/3_",jj,"_did_event_TP_TNP_",num,"_0_nlog.csv"), row.names=F)

did.df.event.1 <- broom::tidy(c.dyn.1)
write.csv(did.df.event.1, paste0("./csv/3_",jj,"_did_event_TP_TNP_",num,"_1_nlog.csv"), row.names=F)

print("DiD csv saved")

}

