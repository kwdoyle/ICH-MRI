library(sqldf)
library(reshape2)
library(dplyr)
library(openxlsx)


readDB <- function(filepath) {
  library(openxlsx)
  
  out <- list()
  wb <- loadWorkbook(filepath)
  nms <- names(wb)
  for (i in 1:length(nms)) {
    out[[i]] <- readWorkbook(wb, sheet=nms[i])
  }
  names(out) <- nms
  
  return(out)
}



filepath <- "E:/SHOP & ICHOP forms/ICHOP_1_17_18_database_full.xlsx"

# read in data
ichmri <- read.xlsx("E:/Experiments/ICH_MRI/MRI_Merged_Data_158patients.xlsx")
ichop <- readDB(filepath)

# set MRI dates
ichmri$MRI_date_Loc<-as.Date(ichmri$MRI_date_Loc,origin = "1899-12-30")

# remove "."s from column names and change to "_"
for (i in 1:length(names(ichop))) {
  nms <- names(ichop[[names(ichop)[i]]])
  
  names(ichop[[names(ichop)[i]]]) <- gsub("\\.", "_", nms)
  
}



# join 3M, 6M, and 12M outcome tables to the MRIs in ichmri
fu_3M <- ichop$`3 month FU`
fu_6M <- ichop$`6 month FU`
fu_12M <- ichop$`12 month FU`

outcomes <- sqldf(
  "select t1.MRN, 
  t2.Lawton_PSMS as Lawton_PSMS_3M, 
  t2.Lawton_IADLs as Lawton_IADLs_3M, 
  t2.TICS as TICS_3M, 
  t2.Barthel as Barthel_3M, 
  t2.mRS as mRS_3M, 
  t2.GOS as GOS_3M,
  t3.Lawton_PSMS as Lawton_PSMS_6M, 
  t3.Lawton_IADLs as Lawton_IADLs_6M, 
  t3.TICS as TICS_6M, 
  t3.Barthel as Barthel_6M, 
  t3.mRS as mRS_6M, 
  t3.GOS as GOS_6M,
  t4.Lawton_PSMS as Lawton_PSMS_12M, 
  t4.Lawton_IADLs as Lawton_IADLs_12M, 
  t4.TICS as TICS_12M, 
  t4.Barthel as Barthel_12M, 
  t4.mRS as mRS_12M, 
  t4.GOS as GOS_12M
  
  from ichmri as t1
  
  left join fu_3M as t2
  on t1.MRN = t2.MRN
  
  left join fu_6M as t3
  on t1.MRN = t3.MRN
  
  left join fu_12M as t4
  on t1.MRN = t4.MRN"
)



# turn all to numeric
outcomes[,1:ncol(outcomes)] <- lapply(outcomes[,1:ncol(outcomes)], as.numeric)


# melt outcomes
outcomes.m <- melt(outcomes, id="MRN")
# turn variable to character
outcomes.m$variable <- as.character(outcomes.m$variable)

# make new column for follow-up month
outcomes.m$fu_month <- NA
for (i in 1:nrow(outcomes.m)) {
  # get measure
  st <- outcomes.m$variable[i]
  # pull out the month measure
  # (this pulls everything before the last underscore and replaces it with nothing,
  # effectively just getting the month number after that last underscore)
  mth <- gsub("([^_]*)_", "", st)
  # only take the number
  mth <- as.numeric(gsub("M", "", mth))
  #mth <- substr(mth, 1, 1)
  # put number in new month column
  outcomes.m$fu_month[i] <- mth
  # remove month from name
  # (this pulls everything after the last underscore and replaces it with nothing,
  # yielding just the measure name)
  outcomes.m$variable[i] <- gsub("_([^_]*)$", "", st) #gsub(".{3}$", "", st) this bit just pulled the last 3 values, but it didn't work for double digit numbers (since it just got one of the digits)
}


# remove patients who have an mRS of 6 at any fu month
# outcomes.alive <- outcomes.m %>%
#   group_by(MRN) %>%
#   #filter((! any((variable=='mRS' & value==6 | is.na(value))) )) #| any(is.na(value)) )
#   filter(any( (variable == 'mRS' & value != 6) | (variable != 'mRS' & is.na(value)) )) # is.na(value) )  ))


# do some dplyr stuff to check who's missing which outcome score
outcomeSums <- outcomes.m %>%
  group_by(MRN, variable) %>%
  summarize(check_sum = sum(value, na.rm=T))

# find patients who ever had an mRS of 6, because they shouldn't be included in this check
deadPatients <- outcomes.m %>%
  filter(variable=='mRS' & value==6) %>%
  select(MRN) %>%
  distinct()

#outcomeSums.alive <- outcomeSums[which(!(outcomeSums$MRN %in% deadPatients$MRN)), ]


# pull patients who are missing outcome measures
# but not in the list of patients who had an mRS of 6
missingOut <- outcomeSums %>%
  filter(! MRN %in% deadPatients$MRN) %>%
  group_by(MRN) %>%
  filter(check_sum == 0)
#filter(! (variable=="mRS" & check_sum==6))


## 49 patients missing outcomes
length(unique(missingOut$MRN))


# look back at original table to make sure
#View(filter(outcomes.m, MRN %in% missingOut$MRN))






###### Other attempts

mRSGOS <- outcomes.m %>%
  filter((variable=='mRS' | variable=='GOS') & is.na(value)==T)

#table(mRSGOS$variable, mRSGOS$fu_month)

no.mRS.GOS.MRN <- data.frame(MRN=unique(mRSGOS$MRN))

##### So we're only concerned with patients missing mRS or GOS??

## join the MRI dates to these MRNs
MRNDates <- select(ichmri, MRN, MRI_date_Loc)

missing.mRSGOS.MRN.Dates <- merge(no.mRS.GOS.MRN, MRNDates, by="MRN")



# filter, by group, where mRS !=6 (but this removes when mRS == NA too)


test <- data.frame(one=c(1,1,1,2,2,2), two=c("a", "b", "c", "a", "b", "c"), three=c(9, "MISSING", 8, 9, 6, 3))
test2 <- data.frame(one=c(1,1,1,2,2,2,3,3,3), two=c("a", "b", "c", "a", "b", "c", "a", "b", "c"), three=c(9, NA, 8, 9, 6, 3, NA, NA, 3))

# filter where b 

test %>%
  group_by(one) %>%
  #filter(!(two == "b" & three == 6 ))
  filter(any( !(two == "b" & three == 6 ) ) | any( (two == "b" & is.na(three))  ))




test %>%
  group_by(one) %>%
  filter(!any(two=="b" & three==6))











#### Filter entire outcome.m but replace NAs with a string: "MISSING" ####
outcomes.m2 <- outcomes.m



outcomes.m2[is.na(outcomes.m2$value),"value"] <- "MISSING"

checknoout2 <- outcomes.m2 %>%
  group_by(MRN) %>%
  filter(!any(variable=="mRS" & value==6)) %>%
  group_by(MRN, variable) %>%
  filter(all(value=="MISSING"))
#filter(value=="MISSING")



length(unique(checknoout2$MRN))


miss_val <- NA_real_

outcomes$GOSE_3M <- miss_val  # 10 means "missing".
outcomes$GOSE_6M <- miss_val  # 10 means "missing".
outcomes$GOSE_12M <- miss_val  # 10 means "missing".

#### Try to calculate GOSE ####
## The TRUEs on the left of the equation are to assign the value on the right side
## when none of the above 'formulas' are met.
gose <- outcomes %>%
  group_by(MRN) %>%
  mutate(
    GOSE_3M = case_when(GOS_3M == 5 ~ 1,
                        GOS_3M == 4 ~ 2,
                        GOS_3M == 3 ~ 3,
                        GOS_3M == 2 ~ 5,
                        GOS_3M == 1 ~ 7,
                        TRUE ~ GOSE_3M),
    
    GOSE_6M = case_when(GOS_6M == 5 ~ 1,
                        GOS_6M == 4 ~ 2,
                        GOS_6M == 3 ~ 3,
                        GOS_6M == 2 ~ 5,
                        GOS_6M == 1 ~ 7,
                        TRUE ~ GOSE_6M),
    
    GOSE_12M = case_when(GOS_12M == 5 ~ 1,
                         GOS_12M == 4 ~ 2,
                         GOS_12M == 3 ~ 3,
                         GOS_12M == 2 ~ 5,
                         GOS_12M == 1 ~ 7,
                         TRUE ~ GOSE_12M)
    
  ) %>%
  mutate(
    GOSE_3M = case_when((GOSE_3M == 3 | is.na(GOSE_3M)) & mRS_3M == 3 ~ 4,
                        (GOSE_3M == 7 | is.na(GOSE_3M)) & mRS_3M == 0 ~ 8,
                        is.na(GOSE_3M) & mRS_3M == 1 ~ 7,
                        is.na(GOSE_3M) & mRS_3M == 2 ~ 5,
                        is.na(GOSE_3M) & mRS_3M %in% 4:5 ~ 3,
                        is.na(GOSE_3M) & mRS_3M == 6 ~ 1,
                        TRUE ~ GOSE_3M),
    
    GOSE_6M = case_when((GOSE_6M == 3 | is.na(GOSE_6M)) & mRS_6M == 3 ~ 4,
                        (GOSE_6M == 7 | is.na(GOSE_6M)) & mRS_6M == 0 ~ 8,
                        is.na(GOSE_6M) & mRS_6M == 1 ~ 7,
                        is.na(GOSE_6M) & mRS_6M == 2 ~ 5,
                        is.na(GOSE_6M) & mRS_6M %in% 4:5 ~ 3,
                        is.na(GOSE_6M) & mRS_6M == 6 ~ 1,
                        TRUE ~ GOSE_6M),
    
    GOSE_12M = case_when((GOSE_12M == 3 | is.na(GOSE_12M)) & mRS_12M == 3 ~ 4,
                         (GOSE_12M == 7 | is.na(GOSE_12M)) & mRS_12M == 0 ~ 8,
                         is.na(GOSE_12M) & mRS_12M == 1 ~ 7,
                         is.na(GOSE_12M) & mRS_12M == 2 ~ 5,
                         is.na(GOSE_12M) & mRS_12M %in% 4:5 ~ 3,
                         is.na(GOSE_12M) & mRS_12M == 6 ~ 1,
                         TRUE ~ GOSE_12M)
  ) %>%
  mutate(
    GOSE_3M = case_when(is.na(GOSE_3M) & Barthel_3M == 100 ~ 7,
                        TRUE ~ GOSE_3M),
    
    GOSE_6M = case_when(is.na(GOSE_6M) & Barthel_3M == 100 ~ 7,
                        TRUE ~ GOSE_6M),
    
    GOSE_12M = case_when(is.na(GOSE_12M) & Barthel_3M == 100 ~ 7,
                         TRUE ~ GOSE_12M)
  ) %>%
  mutate(
    GOSE_3M = case_when(GOSE_3M == 3 & Lawton_PSMS_3M <= 16 ~ 4, #%in% 8:16 ~ 4,
                        GOSE_3M == 5 & Lawton_PSMS_3M <= 9 ~ 6,  #%in% 6:9 ~ 6,
                        is.na(GOSE_3M) & Lawton_PSMS_3M <= 30 ~ 3,
                        is.na(GOSE_3M) & Lawton_PSMS_3M <= 16 ~ 4,
                        is.na(GOSE_3M) & Lawton_PSMS_3M <= 12 ~ 5,
                        is.na(GOSE_3M) & Lawton_PSMS_3M <= 9 ~ 6,
                        is.na(GOSE_3M) & Lawton_PSMS_3M <= 6 ~ 7,
                        TRUE ~ GOSE_3M),
    
    GOSE_6M = case_when(GOSE_6M == 3 & Lawton_PSMS_6M <= 16 ~ 4, #%in% 8:16 ~ 4,
                        GOSE_6M == 5 & Lawton_PSMS_6M <= 9 ~ 6,  #%in% 6:9 ~ 6,
                        is.na(GOSE_6M) & Lawton_PSMS_6M <= 30 ~ 3,
                        is.na(GOSE_6M) & Lawton_PSMS_6M <= 16 ~ 4,
                        is.na(GOSE_6M) & Lawton_PSMS_6M <= 12 ~ 5,
                        is.na(GOSE_6M) & Lawton_PSMS_6M <= 9 ~ 6,
                        is.na(GOSE_6M) & Lawton_PSMS_6M <= 6 ~ 7,
                        TRUE ~ GOSE_6M),
    
    GOSE_12M = case_when(GOSE_12M == 3 & Lawton_PSMS_12M <= 16 ~ 4, #%in% 8:16 ~ 4,
                         GOSE_12M == 5 & Lawton_PSMS_12M <= 9 ~ 6,  #%in% 6:9 ~ 6,
                         is.na(GOSE_12M) & Lawton_PSMS_12M <= 30 ~ 3,
                         is.na(GOSE_12M) & Lawton_PSMS_12M <= 16 ~ 4,
                         is.na(GOSE_12M) & Lawton_PSMS_12M <= 12 ~ 5,
                         is.na(GOSE_12M) & Lawton_PSMS_12M <= 9 ~ 6,
                         is.na(GOSE_12M) & Lawton_PSMS_12M <= 6 ~ 7,
                         TRUE ~ GOSE_12M)
  ) %>%
  mutate(
    GOSE_3M = case_when((GOSE_3M == 5 | is.na(GOSE_3M)) & Lawton_IADLs_3M <= 15 ~ 6,
                        is.na(GOSE_3M) & Lawton_IADLs_3M <= 20 ~ 5,
                        TRUE ~ GOSE_3M),
    
    GOSE_6M = case_when((GOSE_6M == 5 | is.na(GOSE_6M)) & Lawton_IADLs_3M <= 15 ~ 6,
                        is.na(GOSE_6M) & Lawton_IADLs_3M <= 20 ~ 5,
                        TRUE ~ GOSE_6M),
    
    GOSE_12M = case_when((GOSE_12M == 5 | is.na(GOSE_12M)) & Lawton_IADLs_3M <= 15 ~ 6,
                         is.na(GOSE_12M) & Lawton_IADLs_3M <= 20 ~ 5,
                         TRUE ~ GOSE_12M)
  )


# h <- gose %>%
#   mutate(
#     AAA = case_when((GOSE_3M == 5 | is.na(GOSE_3M)) & Lawton_IADLs_3M <= 15 ~ 6,
#                                          is.na(GOSE_3M) & Lawton_IADLs_3M <= 20 ~ 5)
#   )





# patients who are missing GOSE at any month
newtest <- gose[is.na(gose$GOSE_3M) | is.na(gose$GOSE_6M) | is.na(gose$GOSE_12M), c("MRN", "GOSE_3M", "GOSE_6M", "GOSE_12M")]






