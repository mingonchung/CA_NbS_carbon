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
library(fixest)

## setwd
setwd("/projects/mich9173/CA_carbon/output/DiD_input")
input <- fread("2_DiD_data_FP_FNP_NN.csv.gz")

## is.na = 0
input$year.fire[is.na(input$year.fire)] <- 0

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
#  "ppt_00_02", "tmean_00_02", "DEM_slope", "DEM_elevation", "DEM_aspect", "popden_00_02", "citytrvtime_11"
col.list <- c("V1","layer","year",c.list[ii],"US_L3CODE", "year.fire", "wildfire", "treat", "subclass")

input$wildfire[is.na(input$wildfire)] <- 0
input.sub <- input %>% select(all_of(col.list))
input.sub <- input.sub[complete.cases(input.sub),]

if(ii == 1 | ii == 3) {
  input.sub <- subset(input.sub, year > 1989) #1989; 1994
  
} 

gc()
set.seed(160617)

for(jj in c(2,3,4)) {
  
  input.sub2 <- subset(input.sub, wildfire %in% jj)
  input.sub2 <- input.sub2 %>% mutate_if(is.integer,as.numeric)

   input.overyr <- subset(input.sub2, treat==1)
   input.overyr.sb <- unique(input.overyr$subclass)
   print(length(input.overyr.sb))  
   
   input.overyr.sb2 <- sample(input.overyr.sb, 50000, replace = FALSE)
   input.sub2 <- subset(input.sub2, subclass %in% input.overyr.sb2)
  
   ## select treat 1 or 0
   input.sub2 <- subset(input.sub2, treat==0)
   
   # input.sub2 <- input.sub2[complete.cases(input.sub2),]
    print(nrow(input.sub2)) 
  
  print(paste0(c.list[[ii]]," DiD started"))
  print(paste0("wildfire severity ", jj))
  print("input data processed")

###################
### (2) DiD analyses w/ etwfe package
# Estimate the model
# vcov: cluster standard errors at the individual unit level
if(ii==1) {
c.etwfe <-
  etwfe(
    fml  = CFlux_GPP ~ 1, # outcome ~ controls
    tvar = year,        # time variable # calendar
    gvar = year.fire, # group variable # group
    data = input.sub2,       # dataset
    vcov = ~ subclass,  # vcov adjustment (here: clustered)
    cgroup = c("notyet"))
    #xvar = treat) 
} else if(ii==2) {
  c.etwfe <-
    etwfe(
      fml  = CONUS ~ 1, 
      tvar = year,        # time variable
      gvar = year.fire, # group variable
      data = input.sub2,       # dataset
      vcov = ~ subclass,  # vcov adjustment (here: clustered)
      cgroup = c("notyet"))
      #xvar = treat) 
} else if(ii==3) {
  c.etwfe <-
    etwfe(
      fml  = CStocks_Live ~ 1,
      tvar = year,        # time variable
      gvar = year.fire, # group variable
      data = input.sub2,       # dataset
      vcov = ~ subclass,  # vcov adjustment (here: clustered)
      cgroup = c("notyet"))
      #xvar = treat) 
} else if(ii==4) {
  c.etwfe <-
    etwfe(
      fml  = LEMMA ~ 1,
      tvar = year,        # time variable
      gvar = year.fire, # group variable
      data = input.sub2,       # dataset
      vcov = ~ subclass,  # vcov adjustment (here: clustered)
      cgroup = c("notyet"))
      #xvar = treat) 
}

summary(c.etwfe)
print("etwfe done")

########
## (2-2) emfx
# produce a standard margnialeffects object

# recover the average treatment effect on the treated (ATT)
# an increase in the minimum wage leads to an approximate 5 percent decrease in teen employment
c.dtwfe.all <- emfx(c.etwfe, collapse=FALSE, vcov=TRUE) #vcoc=FALSE # simple; group; calendar; event
c.dtwfe.all
print("emfx all done")

# type = "event": recover dynamic treatment effects a la an event study
# suggest that the teen disemployment effect of a minimum wage hike is fairly modest at first (3%), but increase over the next few years (>10%)
c.dtwfe.es <- emfx(c.etwfe, type = "event", collapse=FALSE, vcov=TRUE) #, vcoc=FALSE
print(c.dtwfe.es, style="data.frame")

print("emfx event done")

############
## save output csv
# a standard margnialeffects
c.dtwfe.df.all <- broom::tidy(c.dtwfe.all)
write.csv(c.dtwfe.df.all,paste0("./csv/2_",jj,"_etwfe_all_FP_FNP_",num,"_0_nlog.csv"), row.names=F)
# event marginal effects
c.dtwfe.df.es <- broom::tidy(c.dtwfe.es)
write.csv(c.dtwfe.df.es,paste0("./csv/2_",jj,"_etwfe_event_FP_FNP_",num,"_0_nlog.csv"), row.names=F)

rm(list=c("c.etwfe","c.dtwfe.df.es", "c.dtwfe.df.all"))
print("emfx all & event saved")

gc()

}