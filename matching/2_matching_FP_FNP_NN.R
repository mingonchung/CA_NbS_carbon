#######
## CA C matching
## (1) Comparison between PAs and non-PAs
## 10/14/2023

rm(list=ls())

## number
num <- commandArgs(trailingOnly = TRUE)
num <- as.numeric(num)

## library
library(data.table)
#library(etwfe)
library(stringr)
library(MatchIt)
library(cobalt)
library(ggplot2)
#library(optmatch)

## setwd
setwd("/projects/mich9173/CA_carbon/input/sampling")

# sampling data
all.smpl <- read.csv("All_sampling_sum.csv", header=T, stringsAsFactors = FALSE)

###################
### (1) Input Data
setwd("/projects/mich9173/CA_carbon/input/all_unit_transfer/Dist_sel")
code.list <- list.files(recursive=T, pattern="*_code.csv.gz")

setwd("/projects/mich9173/CA_carbon/input/all_unit_transfer/EnvSoc/")
csv.list <- c("Decid", "Everg", "Herba", "Mixed", "Shrub")
envsoc.list <- list.files(recursive=T, pattern="*.csv.gz")

# layer: 41 (Decid) 42 (Everg) 43 (Mixed)
layer.list <- c(41,42,72,43,71)

### for loop
for(i in c(1,2,4)) {
  
  ## read EnvSoc
  setwd("/projects/mich9173/CA_carbon/input/all_unit_transfer/EnvSoc/")
  envsoc <- fread(envsoc.list[[i]])
  
  ## read code
  setwd("/projects/mich9173/CA_carbon/input/all_unit_transfer/Dist_sel")
  code <- fread(code.list[[i]])

  ## merge
  # add treated year; year.fire
  # add fire severity; wildfire
  envsoc.code <- merge(envsoc[,-1], code[,c(1,2,9,10,24)], by=c("V1", "layer") )
  
  ## subset
  # treat: "F2P1", "F2P2", "F3P1", "F3P2", "F4P1", "F4P2"
  envsoc.code.treat <- subset(envsoc.code, code %in% c("F2P1", "F2P2", "F3P1", "F3P2", "F4P1", "F4P2"))
  
  ## treat sampling

  all.smpl.sub <- subset(all.smpl, layer == layer.list[i])
  all.smpl.V1 <- all.smpl.sub$V1
  
  envsoc.code.treat2 <- subset(envsoc.code.treat, V1 %in% all.smpl.V1)
  envsoc.code.treat2$treat <- 1
  print(paste0(csv.list[i]," - ",nrow(envsoc.code.treat2)))
  
  # control: "F2NP", "F3NP", "F4NP"
  envsoc.code.control <- subset(envsoc.code, code %in% c("F2NP", "F3NP", "F4NP"))
  envsoc.code.control$treat <- 0

  ## rbind all
  envsoc.code.all <- rbind(envsoc.code.treat2, envsoc.code.control)
  
  if(i==1) {
  envsoc.input <- envsoc.code.all} else {
  envsoc.input <- rbind(envsoc.input, envsoc.code.all) 
  }
  
  print(i)
  
}

## subset ecoregion
# delete 7,14,80,81,85 (13?)
#c(1,4,5,6,7,8,9,13,14,78,80,81,85)
eco.sub <- c(1,4,5,8,9,78)
envsoc.input <- subset(envsoc.input, US_L3CODE %in% eco.sub)

## complete.case
envsoc.input$huc8 <- str_sub(envsoc.input$huc12, 1,8)

envsoc.input2 <- envsoc.input[complete.cases(envsoc.input[ , c(2:7,9:11,16:18)]),]

print(paste0(nrow(envsoc.input2)," / ", nrow(envsoc.input)))
print("matching input processed")

###################
## (2) matching
set.seed(160617)

#!!!!!!!!!!!!!!!!
method.name <- "nearest" # nearest; genetic; optimal

match.output <- matchit(treat ~ ppt_00_02 + tmean_00_02 + DEM_slope + DEM_elevation + DEM_aspect + popden_00_02 + citytrvtime_11,
                        data = envsoc.input2,
                        method = method.name, # nearest; genetic; optimal
                        exact = ~ layer + US_L3CODE + huc8 + wildfire,
                        ratio = 1, # one-to-one matching
                        distance = "bart", #glm; gam; gbm; lasso; ridge; elasticnet; rpart; cbps; bart; randomforest; nnet
                        #pop.size = 200, # Genetic population size; bigger is better
                        caliper = 0.25, #c(citytrvtime_11 = 0.25, popden_00_02 = 0.25),
                        mahvars = ~ ppt_00_02 + tmean_00_02 + DEM_slope + DEM_elevation + DEM_aspect + popden_00_02 + citytrvtime_11,
                        replace = FALSE) # matching without replacement        
  print("matching done")
  
  ## save matched data
  match.output.data <- match.data(match.output)
  head(match.output.data)
  
  # setwd
  setwd(paste0("/projects/mich9173/CA_carbon/output/Matching/",method.name))
  write.csv(match.output.data, paste0("2_Matching_data_FP_FNP_",method.name,".csv"), row.names=F)
  
  print("matching data saved")
  
  ###################
  ## (3) Assessing Balance
  ### love plot
  # for both full matching and nearest nieghbor matching simultaneously
  #"m" - mean differences
  #"v" - variance ratios
  #"ks" - Kolmogorov-Smirnov statistics
  #pdf("2_love_plot.pdf", width=4, height=4)
  jpeg(paste0("2_love_plot.jpeg"), width=4, height=4, units="in", res=300)
  
  #love.plot(match.output, binary = "std")
  cobalt::love.plot(treat ~ ppt_00_02 + tmean_00_02 + DEM_slope + DEM_elevation + DEM_aspect + popden_00_02 + citytrvtime_11, data = envsoc.input2, weights=get.w(match.output), 
            stats="m", abs=TRUE, 
            drop.distance = TRUE, 
            thresholds = c(m = .1),
            var.order = "unadjusted", 
            binary = "std",
            limits = list(m = c(0, 1)), #v = c(.3, 6)
            shapes = c("circle", "circle filled"),
            colors = c("blue", "red"),
            sample.names = c("Original", paste0(method.name," Matching")),
            position = "none",
            labels=TRUE,
            title = NULL,
            themes = list(m = theme(legend.position = c(.75, 0.15), 
                                    legend.title = element_blank()))) # legend.key.size = unit(.02, "npc")
  
  dev.off()
  
  print("love plot saved")
  
  ### density plot
  # black line = treated; gray = control; Perfectly overlapping lines indicate good balance 
  # Perfectly overlapping lines indicate good balance
  #pdf("2_density_plot.pdf", width=4, height=4)
  jpeg(paste0("2_density_plot.jpeg"), width=4, height=4, units="in", res=300)
  plot(match.output, type="density", interactive = FALSE, which.xs = c("tmean_00_02", "DEM_elevation", "citytrvtime_11"))
  dev.off()
  print("density plot saved")
  
  ### eQQ plot
  # when values fall on the 45 degree line, the groups are balanced
  #pdf("2_eQQ_plot.pdf", width=4, height=4)
  jpeg(paste0("2_eQQ_plot.jpeg"), width=4, height=4, units="in", res=300)
  plot(match.output, type = "qq", which.xs = c("tmean_00_02", "DEM_elevation", "citytrvtime_11")) 
  dev.off()
  print("eQQ plot saved")
  
  ### eCDF plot
  # Perfectly overlapping lines indicate good balance
  #pdf("2_eCDF_plot.pdf", width=4, height=4)
  jpeg(paste0("2_eCDF_plot.jpeg"), width=4, height=4, units="in", res=300)
  plot(match.output, type = "ecdf", which.xs = c("tmean_00_02", "DEM_elevation", "citytrvtime_11"))
  dev.off()
  print("eCDF plot saved")
  
  ## Standardized Mean Difference (SMD): <.1 and .05 for prognostically important covariates
  ## Variance Ratio (VR): close to 1 indicate good balance, recommendation for VR between .5 and 2.
  ## Empirical Cumulative Density Function statistics (eCDFs): values closer to zero indicating better balance, no specific recommendations
  
  ### summary.matchit()
  # interactions: controls whether balance statistics for all squares and pairwise interactions of covariatesare to be displayed in addition to the covariates.
  # addvariables: allows for balance to be assessed on variables other than those inside the matchit object.
  # standardize: standardized statistics include the standardized mean difference and eCDF statistics; unstandardized statistics include the raw difference in means and eQQ plot statistics
  # pair.dist: controls whether within-pair distances should be computed and displayed
  # un: control whether balance prior to matching should be displayed
  # improvement: whether the percent balance improvement after matching should be displayed
  match.output$X[c("layer")] <- NULL
  match.output$X[c("US_L3CODE")] <- NULL
  match.output$X[c("huc8")] <- NULL
  
  m.sum <- summary(match.output)  # , un = FALSE
  m.sum
  
  print(paste0("summary F done"))
