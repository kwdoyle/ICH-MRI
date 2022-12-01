### New, more robust, script to make heatmaps of any specified subset of the patient cohort
library(openxlsx)
library(reshape2)
library(ggplot2)
library(dplyr)
library(abind)
library(corrplot)
library(reshape2)


setwd("E:/Experiments/BrainFlowMap Stuff/Results_JanClaassen_ICH_DTI")


txtfiles <- list.files(pattern=".txt")
analysis2 <- grep("Analysis2", txtfiles, value=T)
conscdat <- read.xlsx("../MRI DTI set.xlsx", 1)


# extract "columns"; get unique patient IDs from files
split_txt <- strsplit(analysis2, "_")
ID <- sapply(split_txt, "[", 1)
ID.u <- unique(ID)

### behold; the messiest way to replace all instances of "DTI" with "ICH":
# in order to match the names in the text files.
# this does three replacements and uses the output of the prior replacement
# as the input for the next one. it first removes everything after any instance of an "A"
# (because some of the names are like DTI 8A/8B, etc.)
# and then it replaces the DTIs with a space between the number,
# and then replaces the DTIs without a space.
# THEN it removes any blank spaces, which would cause an identical match from the mainTable IDs to appear different
conscdat$DTI.Number <- gsub(" ", "", gsub("DTI", "ICH", gsub("DTI ", "ICH", gsub("A.*", "", conscdat$DTI.Number))))

# convert conscious columns to factor
conscdat[,c("unconsious.on.admission", "unconsious.at.time.of.MRI", "unconscious.on.discharge", "Recovery.of.consciousness")] <- lapply(
  conscdat[,c("unconsious.on.admission", "unconsious.at.time.of.MRI", "unconscious.on.discharge", "Recovery.of.consciousness")],
  as.factor
)



##### make separate arrays for unconscious at dch or not
fa_unconsc <- array(data=NA, dim=c(28,28,0))
fa_consc <- array(data=NA, dim=c(28,28,0))
md_unconsc <- array(data=NA, dim=c(28,28,0))
md_consc <- array(data=NA, dim=c(28,28,0))
ad_unconsc <- array(data=NA, dim=c(28,28,0))
ad_consc <- array(data=NA, dim=c(28,28,0))
rd_unconsc <- array(data=NA, dim=c(28,28,0))
rd_consc <- array(data=NA, dim=c(28,28,0))
fiber_unconsc <- array(data=NA, dim=c(28,28,0))
fiber_consc <- array(data=NA, dim=c(28,28,0))



# set parameters of interest
params <- c("fibercount")



# collect patient cohort to analyze
conscious_time <- "unconscious.on.discharge"
conscious_time <- "unconsious.at.time.of.MRI"
checkUnconscMRI <- 1
"L.side.location"


### All patients
# patients unconscious at discharge (n = 4)
unconscDchID <- filter(conscdat, DTI.Number %in% ID.u & unconscious.on.discharge==1)$DTI.Number
# patients conscious at discharge (n = 10)
conscDchID <- filter(conscdat, DTI.Number %in% ID.u & unconscious.on.discharge==0)$DTI.Number
# patients unconscious at MRI (n = 7)
unconscMRIID <- filter(conscdat, DTI.Number %in% ID.u & unconsious.at.time.of.MRI==1)$DTI.Number
# patients conscious at MRI (n = 7)
conscMRIID <- filter(conscdat, DTI.Number %in% ID.u & unconsious.at.time.of.MRI==0)$DTI.Number


### Only patients unconscious at MRI
# patients unconscious at discharge (n = 4)
unconscMRIunconscDchID <- filter(conscdat, DTI.Number %in% ID.u & unconsious.at.time.of.MRI==1 & unconscious.on.discharge==1)
# patients conscious at discharge (n = 3)
unconscMRIconscDchID <- filter(conscdat, DTI.Number %in% ID.u & unconsious.at.time.of.MRI==1 & unconscious.on.discharge==0)


####### THIS IS NOT HOW THE L AND R SIDEDNESS WILL WORK
# ### Unconscious at MRI patients for recovery of consciousness and for L or R location side
# ## MIGHT NOT WANT TO LOOK AT ONLY UNCONSCIOUS AT MRI PATIENTS FOR THIS, AS FINAL GROUP NUMS WILL BE SMALL
# # patients unconscious at discharge, L side
# unconscDch.LsideID <- filter(conscdat, DTI.Number %in% ID.u & unconsious.at.time.of.MRI==1 & unconscious.on.discharge==1 & L.side.location==1)
# # patients unconscious at discharge, R side
# # (there is one patient with a '?' for the location variable in this)
# unconscDch.RsideID <- filter(conscdat, DTI.Number %in% ID.u & unconsious.at.time.of.MRI==1 & unconscious.on.discharge==1 & L.side.location==0)
# 
# # patients conscious at discharge, L side
# conscDch.LsideID <- filter(conscdat, DTI.Number %in% ID.u & unconsious.at.time.of.MRI==1 & unconscious.on.discharge==0 & L.side.location==1)
# # patients conscious at discharge, R side
# # yeah, for this one, none of the patients have it on the right side.
# conscDch.RsideID <- filter(conscdat, DTI.Number %in% ID.u & unconsious.at.time.of.MRI==1 & unconscious.on.discharge==0 & L.side.location==0)
# 



# get connectivity matrix files for these groups of patients
unconscDchfiles <- grep(paste(unconscDchID, collapse="|"), analysis2, value=T)


# ...maybe this method is actually more confusing.







# CreateHeatmaps <- function(analysis, IDlist, params) {
#   
#   if 
#   
# }






















