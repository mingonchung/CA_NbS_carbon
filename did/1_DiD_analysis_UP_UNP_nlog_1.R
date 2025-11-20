#######
## CA C DiD analysis
## (1) treatment effects between PAs and non-PAs
## 10/14/2023

rm(list=ls())

## number
num <- commandArgs(trailingOnly = TRUE)
num <- as.numeric(num)

## library
# # install.packages("devtools")
# devtools::install_github("bcallaway11/did") # module load cmake

library(data.table)
library(did)
library(etwfe)
library(dplyr)
library(ggplot2)

## setwd
setwd("/projects/mich9173/CA_carbon/output/DiD_input")
input <- fread("1_DiD_data_UP_UNP_NN.csv.gz")

## is.na = 0
input$year.cpad[is.na(input$year.cpad)] <- 0

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
# , "ppt_00_02", "tmean_00_02", "DEM_slope", "DEM_elevation", "DEM_aspect", "popden_00_02", "citytrvtime_11"
col.list <- c("V1","layer","year",c.list[ii],"US_L3CODE", "year.cpad", "treat", "subclass")
input.sub <- input %>% select(all_of(col.list))
input.sub <- input.sub[complete.cases(input.sub),]

if(ii == 1 | ii == 3) {
  input.sub <- subset(input.sub, year > 1989) #1989; 1994
}
# else if (ii == 2) { # CONUS, 2017
#   input.overyr <- subset(input.sub, year.cpad < 2018 & treat==1)
#   input.overyr.sb <- unique(input.overyr$subclass)
#   input.sub <- subset(input.sub, subclass %in% input.overyr.sb)
#   
# } else if (ii == 4) { # LEMMA, 2016
#   input.overyr <- subset(input.sub, year.cpad < 2017 & treat==1)
#   input.overyr.sb <- unique(input.overyr$subclass)
#   input.sub <- subset(input.sub, subclass %in% input.overyr.sb)
# }

## select pre1985
#input.overyr <- subset(input.sub, year.cpad > 1985 & treat==1)
#input.overyr.sb <- unique(input.overyr$subclass)
#input.sub <- subset(input.sub, subclass %in% input.overyr.sb)

## update established year
estb.yr <- fread("All_cpad_yr.csv.gz")
estb.yr.sub <- subset(estb.yr, year.cpad > 0)

input.sub <- merge(input.sub, estb.yr.sub[,c(-1,-4)], by=c("V1","layer"), all.x=T)

input.overyr <- subset(input.sub, year.cpad.x > 0 & year.cpad.y > 0)
input.overyr.sb <- unique(input.overyr$subclass)
input.sub <- subset(input.sub, subclass %in% input.overyr.sb)

## finalized
names(input.sub)[9] <- "year.cpad"
input.sub$year.cpad[is.na(input.sub$year.cpad)] <- 0
print(nrow(input.sub))

set.seed(160617)

# cpad.sum <- input.sub %>% 
#             group_by(year.cpad) %>% 
#             summarize(n=n()/32)
# as.data.frame(cpad.sum)

## subset
input.overyr <- subset(input.sub, treat==1)
input.overyr.sb <- unique(input.overyr$subclass)
print(length(input.overyr.sb))  

if(length(input.overyr.sb) > 50000) {
input.overyr.sb2 <- sample(input.overyr.sb, 50000, replace = FALSE)
input.sub <- subset(input.sub[,-6], subclass %in% input.overyr.sb2)
}
print(nrow(input.sub))

print(paste0(c.list[[ii]]," DiD started"))
print("input data processed")

###################
### (1) DiD analyses w/ did package
# return a class MP object
# estimates of group-time average treatment effects for all groups in all time periods
# estimate group-time average treatment effects without covariates
c.attgt <- att_gt(yname = c.list[[ii]],
                  gname = "year.cpad",
                  idname = "V1",
                  tname = "year",
                  xformla = ~ 1, # without covariates
                  data = input.sub,
                  control_group = "nevertreated") #nevertreated

# summarize the results
# "P-value for pre-test of parallel trends assumption" a Wlad pre-test the paralle trends assumption # P > 0.05: the parallel trends assumption would not be rejected at conventional significance levels
summary(c.attgt)

print("DiD analysis done")

# plot the results
# red plots: pre-treatment group-time average treatment effects
# blue plots: post-treatment group-time average treatment effects
# 95% simultaneous confidence intervals
# set ylim so that all plots have the same scale along y-axis
# There does not appear to be much evidence against the parallel trends assumption. One fails to rejectusing the Wald test reported in summary ; likewise the uniform confi dence bands cover 0 in all pre-treatment periods.
# There is some evidence of negative effects of the minimum wage on employment. Two group-timeaverage treatment effects are negative and statistically different from 0. These results also suggest thatit may be helpful to aggregate the group-time average treatment effects.
#jpeg("1_did_plot.jpeg", width=7, height=4, units="in", res=300)
#ggdid(c.attgt) #, ylim = c(-.3, .3)
#dev.off()

########
## (1-1) Overall effects of participating in the treatment
# "Overall ATT": we estimate that increasing the minimum wage decreased teen employment by 3.1% and the effect is marginally statistically significant
# type="simple": tends to overweight the effect of early-treated groups simply because we observe more of them during post-treatment periods
# type="dynamic": aggregate group-time effects into an event study plot; averaged into average treatment effects at different lengths of exposure to the treatment
# type="group": aggregate group-time average treatment effects into group-specific average treatment effects; The Overall ATT averages - this parameter is a leading choice as an overall summary effect of participating in the treatment
# type="calendar": aggregations across different time periods; y-axis - the average effect of participating in the treatment in a particular time period for all groups that participated in the treatment in that time period
group_effects <- aggte(c.attgt, type = "group", na.rm=TRUE) # simple; dynamic; group; calendar;
summary(group_effects)

print("DiD group summary done")

########
## (1-2) Event Studies
# aggregate the group-time average treatment effects
# aggregate the group-time average treatment effects into a small number of parameters
c.dyn <- aggte(c.attgt, type = "dynamic", na.rm=TRUE)

# "event time": when they first participate in the treatment
# event time=0: the "on impact" effect
# event time=-1: the effect in the period before a unite becomes treated
# checking that this is equal to 0 is potentiall useful as a pre-test
c.dyn.sum <- summary(c.dyn)
c.dyn.sum

# one fails to reject parallel trends in pre-treatment periodsand it looks like somewhat negative effects of the minimum wage on youth employment.
jpeg(paste0("./fig/1_did_dynamic_",c.list[ii],"_nlog.jpeg"), width=7, height=4, units="in", res=300)
ggdid(c.dyn)
dev.off()

print("DiD dynamic plot done")


############
## save output csv
# did all effects
did.df.all <- broom::tidy(group_effects)
write.csv(did.df.all,paste0("./csv/1_did_all_UP_UNP_",num,"_nlog.csv"), row.names=F)
# did event effects
did.df.event <- broom::tidy(c.dyn)
write.csv(did.df.event,paste0("./csv/1_did_event_UP_UNP_",num,"_nlog.csv"), row.names=F)

rm(list=c("c.attgt", "group_effects", "c.dyn"))
print("did csv saved")

###################
### (2) DiD analyses w/ etwfe package
# Estimate the model
# vcov: cluster standard errors at the individual unit level
if(ii==1) {
  c.etwfe <-
    etwfe(
      fml  = CFlux_GPP ~ 1, # outcome ~ controls 
      tvar = year,        # time variable # calendar
      gvar = year.cpad, # group variable # group
      data = input.sub,       # dataset
      vcov = ~ subclass,  # vcov adjustment (here: clustered)
      cgroup = c("never")) #never
} else if(ii==2) {
  c.etwfe <-
    etwfe(
      fml  = CONUS ~ 1, # outcome ~ controls
      tvar = year,        # time variable
      gvar = year.cpad, # group variable
      data = input.sub,       # dataset
      vcov = ~ subclass,  # vcov adjustment (here: clustered)
      cgroup = c("never")) 
} else if(ii==3) {
  c.etwfe <-
    etwfe(
      fml  = CStocks_Live ~ 1, # outcome ~ controls
      tvar = year,        # time variable
      gvar = year.cpad, # group variable
      data = input.sub,       # dataset
      vcov = ~ subclass,  # vcov adjustment (here: clustered)
      cgroup = c("never")) 
} else if(ii==4) {
  c.etwfe <-
    etwfe(
      fml  = LEMMA ~ 1, # outcome ~ controls
      tvar = year,        # time variable
      gvar = year.cpad, # group variable
      data = input.sub,       # dataset
      vcov = ~ subclass,  # vcov adjustment (here: clustered)
      cgroup = c("never")) 
}

summary(c.etwfe)
print("etwfe done")

########
## (2-2) emfx
# produce a standard margnialeffects object

# # recover the average treatment effect on the treated (ATT)
# # an increase in the minimum wage leads to an approximate 5 percent decrease in teen employment
# c.dtwfe.all <- emfx(c.etwfe, collapse=FALSE, vcov=TRUE) #vcoc=FALSE # simple; group; calendar; event
# c.dtwfe.all
# print("emfx all done")

# type = "event": recover dynamic treatment effects a la an event study
# suggest that the teen disemployment effect of a minimum wage hike is fairly modest at first (3%), but increase over the next few years (>10%)
c.dtwfe.es <- emfx(c.etwfe, type = "calendar", collapse=FALSE, vcov=TRUE) #, vcoc=FALSE
print(c.dtwfe.es, style="data.frame")

print("emfx calendar done")

############
## save output csv
# # a standard margnialeffects
# c.dtwfe.df.all <- broom::tidy(c.dtwfe.all)
# write.csv(c.dtwfe.df.all,paste0("./csv/1_etwfe_all_UP_UNP_",num,"_nlog.csv"), row.names=F)
# event marginal effects
c.dtwfe.df.es <- broom::tidy(c.dtwfe.es)
write.csv(c.dtwfe.df.es,paste0("./csv/1_etwfe_calendar_UP_UNP_",num,"_nlog.csv"), row.names=F)

rm(list=c("c.etwfe","c.dtwfe.df.es")) # , "c.dtwfe.df.all"
print("emfx all & calendar saved")
