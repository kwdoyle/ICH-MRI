library(dplyr)
library(reshape2)
library(openxlsx)

# load data, which falls on different sheets between the two files
dat1 <- read.xlsx("E:/Experiments/ICH_MRI/daily_scores1.xlsx", sheet = 3)
dat2 <- read.xlsx("E:/Experiments/ICH_MRI/daily_scores2.xlsx", sheet = 2)
# make main column names the same
names(dat2)[1:4] <- names(dat1)[1:4]

# melt the command score days and values
dat1.m <- melt(dat1, id.vars=c("MRN", "MRI.Date", "ICU.Admission.Date", "ICU.Discharge.Date"),
               variable.name="Command_Score_Day", value.name="Command_Score")

dat2.m <- melt(dat2, id.vars=c("MRN", "MRI.Date", "ICU.Admission.Date", "ICU.Discharge.Date"),
               variable.name="Command_Score_Day", value.name="Command_Score")

# combine the two datasets
dat <- rbind(dat1.m, dat2.m)

# set dates and turn scores to numeric
dat$MRI.Date <- convertToDateTime(dat$MRI.Date)
dat$ICU.Admission.Date <- convertToDateTime(dat$ICU.Admission.Date)
dat$ICU.Discharge.Date <- convertToDateTime(dat$ICU.Discharge.Date)
dat$Command_Score <- as.numeric(dat$Command_Score)
dat$Command_Score_Day <- as.numeric(as.character(dat$Command_Score_Day))

# can remove rows where Command_Score is missing
rm_ix <- which(is.na(dat$Command_Score))
dat <- dat[-rm_ix,]


# add in column for day-of-ICU-discharge
dat <- dat %>%
  group_by(MRN) %>%
  mutate(Day_of_ICU_Dch = as.numeric(as.Date(ICU.Discharge.Date) - as.Date(ICU.Admission.Date))+1,
         Day_of_MRI = as.numeric(as.Date(MRI.Date) - as.Date(ICU.Admission.Date))+1,
         Day_of_Hosp_Dch = as.numeric(max(Command_Score_Day)))



## Want to get:
# Best CS ICU
# Best CS hospital (i.e., entire stay)
# CS on day 14
# CS on day of ICU discharge
# CS day of hospital discharge (i.e., the last day for each patient)


##### What do we do for patients who didn't have any CS scores done while in the ICU?
##### i.e., were admitted and discharged from the ICU on the same day, and then had
##### CS scores afterwards.

##### A: Keep them as NaS--Ben is going to check for their scores.


##### Some patients had their ICU discharge a month+ away from the admission date,
##### and there are only scores up to 30 days.
##### And some patients are missing a discharge date
##### OR they had two dates entered into the discharge date column in the excel file.

cmd_scores <- dat %>%
  group_by(MRN) %>%
  mutate(Best_CS_ICU = max(Command_Score[Command_Score_Day <= Day_of_ICU_Dch | Day_of_ICU_Dch == 0]),
         ## Need to change this to only look at scores after the MRI day
         Best_CS_Hosp = max(Command_Score[Command_Score_Day >= Day_of_MRI], na.rm=T),
         CS_Day_14 = ifelse(is.na(mean(Command_Score[Command_Score_Day == 14])), 
                            yes=ifelse(is.na(mean(Command_Score[Command_Score_Day == 13])),  # if want a range around day 14, can change this to get the max command score where command score day is between a range
                                       yes=ifelse(is.na(mean(Command_Score[Command_Score_Day == 15])), yes=NA, 
                                                  no=mean(Command_Score[Command_Score_Day == 15]) ), 
                                       no=mean(Command_Score[Command_Score_Day == 13]) ), 
                            no=mean(Command_Score[Command_Score_Day == 14]) ),
         
         CS_ICU_Dch = ifelse(is.na(mean(Command_Score[Command_Score_Day == Day_of_ICU_Dch])), 
                             yes=ifelse(is.na(mean(Command_Score[Command_Score_Day == Day_of_ICU_Dch-1])), 
                                        yes=ifelse(is.na(mean(Command_Score[Command_Score_Day == Day_of_ICU_Dch+1])), 
                                                   yes=mean(Command_Score[Command_Score_Day == max(Command_Score_Day)]), 
                                                   no=mean(Command_Score[Command_Score_Day == Day_of_ICU_Dch+1]) ), 
                                        no=mean(Command_Score[Command_Score_Day == Day_of_ICU_Dch-1]) ), 
                             no=mean(Command_Score[Command_Score_Day == Day_of_ICU_Dch]) ),
         
         CS_Hosp_Dch = ifelse(is.na(mean(Command_Score[Command_Score_Day == max(Command_Score_Day, na.rm=T)])), 
                              yes=ifelse(is.na(mean(Command_Score[Command_Score_Day == max(Command_Score_Day, na.rm=T)-1])), 
                                         yes=ifelse(is.na(mean(Command_Score[Command_Score_Day == max(Command_Score_Day, na.rm=T)+1])), yes=NA, 
                                                    no=mean(Command_Score[Command_Score_Day == max(Command_Score_Day, na.rm=T)+1]) ), 
                                         no=mean(Command_Score[Command_Score_Day == max(Command_Score_Day, na.rm=T)-1]) ), 
                              no=mean(Command_Score[Command_Score_Day == max(Command_Score_Day, na.rm=T)]) ))




# save uncleaned table to analyze the days too
write.xlsx(cmd_scores, "E:/Experiments/ICH_MRI/Command_scores_and_days.xlsx", rowNames=F)


# clean the table
cmd_scores.c <- cmd_scores %>%
  select(MRN, Best_CS_ICU, Best_CS_Hosp, CS_Day_14, CS_ICU_Dch, CS_Hosp_Dch) %>%
  distinct()

# save table to use in ICH_MRI_v6 markdown script to prepare data
write.xlsx(cmd_scores.c, "E:/Experiments/ICH_MRI/Command_scores_from_daily.xlsx", rowNames=F)
