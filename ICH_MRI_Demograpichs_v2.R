library(openxlsx)
library(dplyr)


# read the ichop database
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


# Beware ye who enter here: this is the worst function ever. I am sorry.
FindHemmorageVolumes <- function(MRNs, MRI=TRUE) {
  # for each table, merge the MRNs related to the ICHOP IDs in,
  # then find the volumes
  tmp <- merge(ichopdb$`Admission data`[,c("Subject.ID", "MRN")], ichopdb$`CT Scan (Admission)`, by="Subject.ID", all.y=T)
  # set the initial MRNs for either those with MRI or those without.
  if (MRI==TRUE) {
    tmp <- filter(tmp, MRN %in% MRNs)
  } else {
    tmp <- filter(tmp, !MRN %in% MRNs)
  }
  
  vol <- select(tmp, MRN, Hematoma.volume)
  # set initial main table as result from the asmission CT data
  main_table <- vol
  
  # now begin the nightmareish tunnel of checking every other table.
  missing <- which(is.na(vol$Hematoma.volume))
  if (any(missing)) {
    # look in the next table
    miss_MRN <- vol$MRN[missing]
    
    # do the same as above
    tmp <- merge(ichopdb$`Admission data`[,c("Subject.ID", "MRN")], ichopdb$`CT Scan (2nd in 24 hr)`, by="Subject.ID", all.y=T)
    tmp <- filter(tmp, MRN %in% miss_MRN)
    vol <- select(tmp, MRN, Hematoma.volume)
    
    # put the found ones in here
    found <- which(!is.na(vol$Hematoma.volume))
    for (ix in found) {
      mrn <- vol$MRN[ix]
      mt_ix <- which(main_table$MRN==mrn)
      # replace missing value
      main_table$Hematoma.volume[mt_ix] <- vol$Hematoma.volume[ix]
    }
    
    missing <- which(is.na(vol$Hematoma.volume))
    if (any(missing)) {
      miss_MRN <- vol$MRN[missing]
      
      tmp <- merge(ichopdb$`Admission data`[,c("Subject.ID", "MRN")], ichopdb$`CT Scan (3rd in 24hr)`, by="Subject.ID", all.y=T)
      tmp <- filter(tmp, MRN %in% miss_MRN)
      vol <- select(tmp, MRN, Hematoma.volume)
      
      # stuff for nonmissing
      found <- which(!is.na(vol$Hematoma.volume))
      for (ix in found) {
        mrn <- vol$MRN[ix]
        mt_ix <- which(main_table$MRN==mrn)
        # replace missing value
        main_table$Hematoma.volume[mt_ix] <- vol$Hematoma.volume[ix]
      }
      
      missing <- which(is.na(vol$Hematoma.volume))
      if (any(missing)) {
        miss_MRN <- vol$MRN[missing]
        
        tmp <- merge(ichopdb$`Admission data`[,c("Subject.ID", "MRN")], ichopdb$`CT Scan (3-6 day)`, by="Subject.ID", all.y=T)
        tmp <- filter(tmp, MRN %in% miss_MRN)
        vol <- select(tmp, MRN, Hematoma.volume)
        
        # stuff for found
        found <- which(!is.na(vol$Hematoma.volume))
        for (ix in found) {
          mrn <- vol$MRN[ix]
          mt_ix <- which(main_table$MRN==mrn)
          # replace missing value
          main_table$Hematoma.volume[mt_ix] <- vol$Hematoma.volume[ix]
        }
        
        missing <- which(is.na(vol$Hematoma.volume))
        if (any(missing)) {
          miss_MRN <- vol$MRN[missing]
          
          tmp <- merge(ichopdb$`Admission data`[,c("Subject.ID", "MRN")], ichopdb$`CT Scan (OSH)`, by="Subject.ID", all.y=T)
          tmp <- filter(tmp, MRN %in% miss_MRN)
          vol <- select(tmp, MRN, Hematoma.volume)
          
          found <- which(!is.na(vol$Hematoma.volume))
          for (ix in found) {
            mrn <- vol$MRN[ix]
            mt_ix <- which(main_table$MRN==mrn)
            # replace missing value
            main_table$Hematoma.volume[mt_ix] <- vol$Hematoma.volume[ix]
          }
          
          missing <- which(is.na(vol$Hematoma.volume))
          if (any(missing)) {
            print(paste("No hematoma volumes for", vol$MRN[missing]))
          }
          
          return(main_table)
          
        } else {
          return(main_table)
        }
        
        
      } else {
        return(main_table)
      }
      
    } else {
      return(main_table)
    }
    
    
  } else {
    return(main_table)
  }
  
}


# get table of lobar and deep yes/no for patients
GetLobarDeep <- function(MRNs) {
  tmp <- ichopdb$Etiology %>%
    left_join(select(ichopdb$`Admission data`, Subject.ID, MRN), by="Subject.ID") %>%
    filter(MRN %in% MRNs) %>%
    mutate(lobar = ifelse(grepl("1", Location2), yes=1, 
                          no=ifelse(grepl("6", Location2), yes=1, no=0)),
           deep = ifelse(grepl("2", Location2), yes=1,
                         no=ifelse(grepl("3", Location2), yes=1, no=0))) %>%
    select(MRN, lobar, deep)
  
  return(tmp)
}


GetInfratentForNoMRI <- function(MRNs) {
  tmp <- ichopdb$Etiology %>%
    left_join(select(ichopdb$`Admission data`, Subject.ID, MRN), by="Subject.ID") %>%
    filter(MRN %in% MRNs) %>%
    # non-strict search
    mutate(Has_Infratent_from_ICHOPdat = ifelse(grepl("4", Location2) | grepl("5", Location2), yes=1, no=0)) %>%
    
    # strict search
    #mutate(Has_Infratent_from_ICHOPdat = ifelse(Location2=="4" | Location2=="5" | Location2=="4, 5", yes=1, no=0)) %>%
    select(MRN, Has_Infratent_from_ICHOPdat)
  
  return(tmp)
}


## MRNs needs to be: ichmri[,1, drop=F]
FindIVHScores <- function(MRNs, MRI=TRUE) {
  
  if (class(MRNs) != "data.frame") {
    stop("Need to provide MRNs as a data frame and not just a vector, i.e., use drop=FALSE when subsetting")
  }
  
  if (MRI==TRUE) {
    test_tab <- merge(MRNs, select(ichopdb$`Admission data`, MRN, Subject.ID))
  } else {
    test_tab <- select(ichopdb$`Admission data`, MRN, Subject.ID)
    test_tab <- filter(test_tab, !MRN %in% MRNs$MRN)
  }
  
  
  check1 <- merge(test_tab, select(ichopdb$`CT Scan (Admission)`, Subject.ID, IVH.Score.total.volume), by="Subject.ID", all.x=T)
  check2 <- merge(test_tab, select(ichopdb$`CT Scan (2nd in 24 hr)`, Subject.ID, IVH.Score.total.volume), by="Subject.ID", all.x=T)
  check3 <- merge(test_tab, select(ichopdb$`CT Scan (3rd in 24hr)`, Subject.ID, IVH.Score.total.volume), by="Subject.ID", all.x=T)
  check4 <- merge(test_tab, select(ichopdb$`CT Scan (3-6 day)`, Subject.ID, IVHS), by="Subject.ID", all.x=T)
  check5 <- merge(test_tab, select(ichopdb$`CT Scan (OSH)`, Subject.ID, IVH.Score.total.volume), by="Subject.ID", all.x=T)
  
  all_t <- cbind(check1, IVH2nd=check2[,3], IVH3rd=check3[,3], IVH3to6=check4[,3], IVHOSH=check5[,3])
  names(all_t)[3] <- "IVHAdm"
  
  all_t_full <- all_t %>%
    mutate(All_Scores = ifelse(is.na(IVHAdm), 
                               yes=ifelse(is.na(IVH2nd), 
                                          yes=ifelse(is.na(IVH3rd), 
                                                     yes=ifelse(is.na(IVH3to6), 
                                                                yes=ifelse(is.na(IVHOSH), yes=NA, no=IVHOSH), 
                                                                no=IVH3to6), 
                                                     no=IVH3rd), 
                                          no=IVH2nd), 
                               no=IVHAdm))
  
  
  
  return(all_t_full)
  
}


# functions to use in summarise_all to get the 25th and 75th quantiles seperately,
# and as a single value
get_25th <- function(x, na.rm=T) {
  out <- quantile(x, na.rm=T)[2]
  return(out)
}

get_75th <- function(x, na.rm=T) {
  out <- quantile(x, na.rm=T)[4]
  return(out)
}


FindNumberNA <- function(x) {
  out <- length(which(is.na(x)))
  return(out)
}


runUnivariateRegress <- function(table, depen_vars, indep_vars) {
  
  output_table <- data.frame()
  
  dep <- paste(depen_vars)
  response <- paste("cbind(",dep,") ~ ", sep="")
  
  for (i in 1:length(indep_vars)) {
    
    # see if cbinding a single variable is the same as supplying just the variable. I think it is
    
    effect <- paste(indep_vars[i])
    
    formula <- paste(response, effect, sep="")
    form <- as.formula(formula)
    
    fit <- try(glm(form, family=binomial, data=table), silent=F)
    if (any(class(fit) == "try-error")) {
      message(paste("variable", indep_vars[i]," has error"))
      return(NA)
    }
    smry <- summary(fit)
    
    params <- try(matrix(smry$coefficients[2,], nrow=1, ncol=length(smry$coefficients[2,]), byrow=T), silent=T)
    # add the confidence intervals, but first check if model was able to be fit for this independent variable.
    if (class(params) != "try-error") {
      params <- matrix(c(params, c(confint.default(fit)[2,])), nrow=1)
      colnames(params) <- c(colnames(smry$coefficients), paste("CI", colnames(confint(fit))))
      rownames(params) <- rownames(smry$coefficients)[2]
      
      output_table <- rbind(output_table, as.data.frame(params))
    } else {
      params <- matrix(c(rep(NA, length(output_table))), nrow=1)
      # I think this will break if the very first row to be added has NAs/is a try-error, since then there'll be no colnames/rownames to choose from.
      colnames(params) <- colnames(output_table)
      rownames(params) <- effect
      output_table <- rbind(output_table, as.data.frame(params))
    }
    
    
  }
  
  # mark p < 0.05 w/ star
  output_table$`Pr(>|z|)` <-  round(output_table$`Pr(>|z|)`, digits=3)
  output_table$`Pr(>|z|)` <-  ifelse(output_table$`Pr(>|z|)` < 0.05, yes=paste(output_table$`Pr(>|z|)`, "*"), no=output_table$`Pr(>|z|)`)
  
  
  return(output_table)
  
  
  
}






# list of 158 patients who had MRI
ichmri <- read.xlsx("~/Experiments/ICH_MRI/MRI_Merged_Data_158patients.xlsx")
# ichop database
ichopdb <- readDB("/mnt/drives/nwhome/Data/SHOP & ICHOP forms/ICHOP_1_17_18_database_full.xlsx")
# CS data
new_CS <- read.xlsx("~/Experiments/ICH_MRI/Command_scores_and_days.xlsx")
new_CS$MRN <- as.numeric(new_CS$MRN)
new_CS[,c("MRI.Date", "ICU.Admission.Date", "ICU.Discharge.Date")] <- lapply(new_CS[,c("MRI.Date", "ICU.Admission.Date", "ICU.Discharge.Date")],
                                                                             convertToDate)

# outcome scores
outcomes <- read.xlsx("~/Experiments/ICH_MRI/ICH_MRI_Patient_Outcomes.xlsx")
# create column of "all GOSE" by looking at all month follow ups and
# taking the most recent GOSE, if it exists                             ))
outcomes <- outcomes %>%
  mutate(GOSE_All = if_else(!is.na(GOSE_12M), GOSE_12M, 
                            if_else(!is.na(GOSE_6M), GOSE_6M,
                                    if_else(!is.na(GOSE_3M), GOSE_3M, GOSE_3M))))

# IVH info
IVH_info <- read.csv("~/Experiments/ICH_MRI/Has_IVH_info.csv")


# get ICH onset dates
# get dch dates for the 158 patients from ichopdb
onset_dates <- ichopdb$`Admission data` %>%
  filter(MRN %in% unique(ichmri$MRN)) %>%
  select(MRN, onset.date) %>%
  distinct()

# set dates, obnoxiously because some read in differently
# get indices
dumb_ix <- which(!is.na(as.numeric(onset_dates$onset.date)))
#onset_dates2$onset.date[dumb_ix] <- gsub("-", "/", as.character(convertToDate(onset_dates2$onset.date[dumb_ix])))
# change those dates
actual_dumb_dates <- convertToDate(onset_dates$onset.date[dumb_ix])
# change the other dates
onset_dates$onset.date <- as.Date(onset_dates$onset.date, format="%m/%d/%Y")
# assign the now NAs with the dumb dates
onset_dates$onset.date[dumb_ix] <- actual_dumb_dates
# merge it into the new_CS table
new_CS <- merge(new_CS, onset_dates, by="MRN", all.x=T)
# make time-to columns
new_CS$Time_ICH_to_MRI <- as.numeric(difftime(new_CS$MRI.Date, new_CS$onset.date, units="days"))
new_CS$Time_MRI_to_ICU_dch <- as.numeric(difftime(new_CS$ICU.Discharge.Date, new_CS$MRI.Date, units="days"))

# set dates
ichmri$MRI_date_Loc <- convertToDate(ichmri$MRI_date_Loc)
ichmri$mri_date <- convertToDate(ichmri$mri_date)
ichmri$Discharge_date <- convertToDate(ichmri$Discharge_date)
# apparently there's one missing value in MRI_date_Loc, so use mri_date instead.

# normalize brain volumes
min_val <- which(ichmri$Brain_vol == min(ichmri$Brain_vol))
mBV <- mean(ichmri$Brain_vol[-min_val])
# replace min value
ichmri$Brain_vol[min_val] <- mBV

# normalize hemmorage vol
ichmri$Hg_vol_Norm <- (ichmri$Hg_vol * ichmri$Brain_vol) / mBV
# normalize edema vol
ichmri$Ed_vol_Norm <- (ichmri$Ed_vol * ichmri$Brain_vol) / mBV
# normalize MLS
ichmri$MLS_Norm <- (ichmri$`MLS.[mm]` * ichmri$Brain_vol) / mBV


# admission data mrns as numeric
ichopdb$`Admission data`$MRN <- as.numeric(ichopdb$`Admission data`$MRN)
# admission dates as dates
ichopdb$`Admission data`$CUMC.admission.date <- convertToDate(ichopdb$`Admission data`$CUMC.admission.date)
ichopdb$Discharge$Discharge.Date <- convertToDate(ichopdb$Discharge$Discharge.Date)
ichopdb$`Admission data`$NICU.discharge.date <- convertToDate(ichopdb$`Admission data`$NICU.discharge.date)
ichopdb$`Admission data`$NICU.admission.date <- convertToDate(ichopdb$`Admission data`$NICU.admission.date)


# set missing ed vols to 0 like in the analysis script
ichmri$Ed_vol[which(is.na(ichmri$Ed_vol))] <- 0
# change volumns to mL
ichmri$Hg_vol_mL <- ichmri$Hg_vol / 1000

## etiologies
etio_table <- merge(ichopdb$`Admission data`[,c("Subject.ID", "MRN")], ichopdb$Etiology, by="Subject.ID", all.y=T)

# EVD and clot evac
evd_table <- merge(select(ichopdb$`Admission data`, Subject.ID, MRN), 
                   select(ichopdb$Procedures, Subject.ID, EVD.placement, Surgical.evacuation), by="Subject.ID", all.x=T)


# days in hospital/ICU
new_CS2 <- filter(new_CS, MRN %in% ichmri$MRN)

# patients with no MRI hospital stay length
Stay_Length_noMRI <- ichopdb$`Admission data` %>%
  left_join(ichopdb$Discharge, by = "Subject.ID") %>%
  filter(!MRN %in% ichmri$MRN) %>%
  mutate(Hosp_Stay_Len_noMRI = as.numeric(difftime(Discharge.Date, CUMC.admission.date, units="days")),
         ICU_Stay_Len_noMRI = as.numeric(difftime(NICU.discharge.date, NICU.admission.date, units="days"))) %>%
  select(MRN, Hosp_Stay_Len_noMRI, ICU_Stay_Len_noMRI)






# merge everything into one table
table1 <- ichopdb$`Admission data`[,c("MRN"), drop=F]
table1$MRN <- as.numeric(table1$MRN)
# has MRI
table1$MRI <- ifelse(table1$MRN %in% ichmri$MRN, 1, 0)

# load in CT volumns
CT_vol <- FindHemmorageVolumes(unique(table1$MRN))
CT_vol$Hematoma.volume <- as.numeric(CT_vol$Hematoma.volume)

# lobar deep numbers
lobardeep_table <- GetLobarDeep(unique(table1$MRN))
# has infratent, calculated using ICHOP data (used to calculate for patients without MRI)
noMRIpats <- unique(filter(table1, MRI==0)$MRN)
has_infratent_table_noMRI <- GetInfratentForNoMRI(noMRIpats)
# number ivh via ivh scores
IVH_out <- FindIVHScores(unique(table1[,1,drop=F]))

# data from admission table
table1 <- merge(table1, select(ichopdb$`Admission data`, MRN, Age, Sex, Ethnicity, `GCS.(CUMC)`, FUNC, ICH.Score), my="MRN", all.x=T)
# etiologies
table1 <- merge(table1, select(etio_table, MRN, Etiology), by="MRN", all.x=T)
## 3 month outcomes
table1 <- merge(table1, select(ichopdb$`3 month FU`, MRN, mRS, GOS), by="MRN", all.x=T)
# ct volumes
table1 <- merge(table1, CT_vol, by="MRN", all.x=T)
# infratent and volumes from ichmri
table1 <- merge(table1, select(ichmri, MRN, infra_tent, Hg_vol_Norm, Ed_vol_Norm, MLS_Norm), by="MRN", all.x=T)
# lobar/deep
table1 <- merge(table1, lobardeep_table, by="MRN", all.x=T)
# has infratent from ICHOP for no MRI patients
table1 <- merge(table1, has_infratent_table_noMRI, by="MRN", all.x=T)
# EVD
table1 <- merge(table1, select(evd_table, MRN, EVD.placement, Surgical.evacuation), by="MRN", all.x=T)
# IVH scores
table1 <- merge(table1, select(IVH_out, MRN, All_Scores), by="MRN", all.x=T)
names(table1)[names(table1)=="All_Scores"] <- "IVH_Score"
# CS levels
table1 <- merge(table1, select(new_CS2, MRN, Day_of_ICU_Dch, Day_of_Hosp_Dch, Time_ICH_to_MRI, Time_MRI_to_ICU_dch, onset.date), by="MRN", all.x=T)
# hosp/icu stay lengths for no MRI patients
table1 <- merge(table1, Stay_Length_noMRI, by="MRN", all.x=T)
# GOSE scores
table1 <- merge(table1, select(outcomes, MRN, GOSE_All), by="MRN", all.x=T)

# remove duplicate rows
table1 <- unique(table1)
# edit names
names(table1)[names(table1)=="GCS.(CUMC)"] <- "GCS"

# set to numeric
table1[,c("Age", "GCS", "FUNC", "ICH.Score", "GOS", "Hematoma.volume")] <- lapply(table1[,c("Age", "GCS", "FUNC", "ICH.Score", "GOS", "Hematoma.volume")],
                                                               as.numeric)


# recode etiology, (gender, race), ICH volumes into groups, infra_tent into 0/1, IVH_score into just 0/1 for if have or not 0 vs a number,
# and a new column for if dead or not at 3 months
table1$Sex <- dplyr::recode(table1$Sex, '2'="Female", '1'="Male")
#table1$Ethnicity <- recode(table1$Ethnicity, '0'="White", '1'="Black", '2'="Asian", '3'="Hispanic", '4'="Other", '5'="Unknown")
table1$Ethnicity <- recode(table1$Ethnicity, '0'="White", '1'="Non-white", '2'="Non-white", '3'="Non-white", '4'="Non-white", '5'="Non-white")
table1$New_Etiology <- ifelse(grepl(1, table1$Etiology), yes="HTN", 
                              no=ifelse(grepl(2, table1$Etiology), yes="Amyloid", 
                                        no=ifelse(grepl(7, table1$Etiology) | grepl(9, table1$Etiology), yes="Coagulopathy/Anticoag",
                                                  no="Other")))

# set any '11's to "Other" (because they were included in the grepl(1) check)
table1$New_Etiology <- ifelse(grepl(11, table1$Etiology), yes="Other", no=table1$New_Etiology)
# set NAs back to NAs
table1$New_Etiology <- ifelse(is.na(table1$Etiology), yes=NA, no=table1$New_Etiology)
# ...if someone has both 1 and 7 (or 1 and 9, or 1 and 2), then they are just marked as "HTN"....
# never mind, I get around this by just summing the percentages where ever a 1 occurs and an 11 doesn't occur
# and the same for the other categories.

table1$Volume_Group <- ifelse(table1$Hematoma.volume < 30, yes="<30",
                              no=ifelse(table1$Hematoma.volume >=30 & table1$Hematoma.volume <= 60, yes=">=30_<=60",
                                        no=ifelse(table1$Hematoma.volume > 60, yes=">60", no=NA)))

table1$Has_Infratent <- ifelse(table1$infra_tent != "None", yes=1, no=0)
table1$Has_Infratent <- ifelse(is.na(table1$infra_tent), yes=NA, no=table1$Has_Infratent)

# table1$Has_IVH <- ifelse(table1$IVH_Score != 0, yes=1, no=0)
# table1$Has_IVH <- ifelse(is.na(table1$IVH_Score), yes=NA, no=table1$Has_IVH)
## Now using IVH info from MRI data
table1 <- merge(table1, IVH_info, by="MRN", all.x=T)


table1$Dead_3M <- ifelse(table1$mRS == 6, yes=1, no=0)
table1$Dead_3M <- ifelse(is.na(table1$mRS), yes=NA, no=table1$Dead_3M)


## Convert GOSE to GOS to mitigate the 51 missing regular GOS scores we have
table1$GOS_new <- ifelse(table1$GOSE_All==8 | table1$GOSE_All==7, yes=5,
                         no=ifelse(table1$GOSE_All==6 | table1$GOSE_All==5, yes=4,
                                   no=ifelse(table1$GOSE_All==4 | table1$GOSE_All==3, yes=3,
                                             no=ifelse(table1$GOSE_All==2, yes=2,
                                                       no=ifelse(table1$GOSE_All==1, yes=1, no=NA)))))



## dichotomize age for above and below median
med_Age <- median(table1$Age, na.rm=T)
table1$Age_Dichot <- ifelse(table1$Age < med_Age, yes=0, no=1)



#### make combined 'has infratent' column using the MRI data and the ICHOP no MRI data so I can run a fisher test on it ####
table1$Has_Infratent_Combined <- ifelse(is.na(table1$Has_Infratent), yes=table1$Has_Infratent_from_ICHOPdat, no=table1$Has_Infratent)
#### make combined 'hosp stay' and 'ICU stay' columns using the MRI data and the ICHOP no MRI data ####
table1$Hosp_Stay_Combined <- ifelse(is.na(table1$Day_of_Hosp_Dch), yes=table1$Hosp_Stay_Len_noMRI, no=table1$Day_of_Hosp_Dch)
table1$ICU_Stay_Combined <- ifelse(is.na(table1$Day_of_ICU_Dch), yes=table1$ICU_Stay_Len_noMRI, no=table1$Day_of_ICU_Dch)







#### medians and count/freqs ####

# numeric vars
numeric_vars <- c("Age", "GCS", "ICH.Score", "FUNC", "Day_of_ICU_Dch", "Day_of_Hosp_Dch", "GOS_new", "mRS", "GOSE_All",
                  "Hosp_Stay_Combined", "ICU_Stay_Combined")
# categorical vars
cat_vars <- c("Sex", "Ethnicity", "Etiology", "lobar", "deep", "Has_Infratent", "Volume_Group", "Has_IVH",
              "EVD.placement", "Surgical.evacuation", "Dead_3M",
              "Has_Infratent_from_ICHOPdat", "Has_Infratent_Combined")



## medians [IQR] (get medians and IQRs for the numeric variables)
med_out <- table1 %>%
  select_(.dots=c("MRI", numeric_vars)) %>%
  group_by(MRI) %>%
    summarise_all(funs(median, get_25th, get_75th, .args=list(na.rm=T))) # use the two functions I made here to get the 25th and 75th quantiles

# reorder columns to keep the IQR and median for each parameter together,
# but keep "MRI" as the first column
med_out <- med_out[,c("MRI", sort(names(med_out))[-which(sort(names(med_out))=="MRI")])]





#lapply(table1[,numeric_vars], FindNumberNA)

mri_p <- table1 %>%
  select_(.dots=c("MRI", numeric_vars, "Dead_3M")) %>% 
  filter(MRI==1)

lapply(mri_p, FindNumberNA)



nomri_p <- table1 %>%
  select_(.dots=c("MRI", numeric_vars, "Dead_3M")) %>% 
  filter(MRI==0)


lapply(nomri_p, FindNumberNA)



## n (%) [makes counts and frequencies of the values in each categorical variable]
var_tabs <- list()
for (i in 1:length(cat_vars)) {
  
  assign( cat_vars[i],
          
          table1 %>%
            group_by_("MRI", cat_vars[i]) %>%
            filter_(
              paste('!is.na(', cat_vars[i], ')', sep="")
            ) %>%
            summarize(tot = n()) %>%
            group_by(MRI) %>%
            mutate(tot_sum = sum(tot)) %>%
            mutate(percentage = (tot / tot_sum) * 100)
  )
  
  var_tabs[[i]] <- eval(as.name(cat_vars[i]))
  names(var_tabs)[i] <- cat_vars[i]
}



# do stuff for etiologies seperately, because the data is stored more ..compliclatedly
#ix <- which(  grepl(1, var_tabs$Etiology$Etiology)  & !grepl(11, var_tabs$Etiology$Etiology) )

# This looks through the Etiology counts/percents table and sums up percentages
# for patients with each etiology
sum_percents <- var_tabs$Etiology %>%
  group_by(MRI) %>%
  # the percents won't add to 100 because a patient can be counted in more than one group
  # if they have two or more of these etiologies
  summarise(HTNsum = sum(percentage[grepl(1, Etiology)  & !grepl(11, Etiology)], na.rm=T), # looks for those with "1" and not "11" (HTN)
            AmylSum = sum(percentage[grepl(2, Etiology)], na.rm=T),                        # looks for those with "2" (Amyloid)
            CoagSum = sum(percentage[grepl(7, Etiology) | grepl(9, Etiology)], na.rm=T),   # looks for those with "7" or "9" (all coagulopathy)
            OtherSum = sum(percentage[!(grepl(1, Etiology)  & !grepl(11, Etiology)) &      # everyone without the above groups, and also isn't NA
                                        !(grepl(2, Etiology)) &
                                        !(grepl(7, Etiology) | grepl(9, Etiology)) &
                                        !is.na(Etiology)], na.rm=T),
            
            HTNtot = sum(tot[grepl(1, Etiology)  & !grepl(11, Etiology)], na.rm=T),
            Amyltot = sum(tot[grepl(2, Etiology)], na.rm=T),
            Coagtot = sum(tot[grepl(7, Etiology) | grepl(9, Etiology)], na.rm=T),
            Othertot = sum(tot[!(grepl(1, Etiology)  & !grepl(11, Etiology)) &      # everyone without the above groups, and also isn't NA
                                 !(grepl(2, Etiology)) &
                                 !(grepl(7, Etiology) | grepl(9, Etiology)) &
                                 !is.na(Etiology)], na.rm=T))



## manual etiology fisher
# I think these counts were slighly wrong
# tb <- as.table(matrix(c(77,255,58,64,16,90,5,13), nrow=4, ncol=2, byrow = T))
# fisher.test(tb)



# trying to make table of number of patients in each group
# by multiplying the percentages in each row by the total number of patients in that group
etiology_patients <- matrix(nrow=2, ncol=5)
etiology_patients[,1] <- as.matrix(sum_percents[,1])

# ..or just do it on the percentages.
new_tb <- as.table(as.matrix(sum_percents[,c(2,3,4,5)]))
eti_res <- fisher.test(new_tb) #, workspace = 1e9)







#### Extract FUNC and GCS for the 158 MRI patients
FUNC_GCS_ICH_MRI <- table1 %>%
  filter(MRI==1) %>%
  select(MRN, FUNC, GCS, ICH.Score)

# save it
#write.xlsx(FUNC_GCS_MRI, "E:/Experiments/ICH_MRI/FUNC_GCS_ICH_Score_158_MRI_patients.xlsx", rowNames=F)






#### stats ####

# continuous vars (it's just Age. NOT ANYMORE IT ISN'T)
numeric_test_vars <- c("Age", "ICU_Stay_Combined", "Hosp_Stay_Combined")
# categorical vars
cat_test_vars <- c("Sex", "Ethnicity", "lobar", "deep", "Volume_Group", "Has_IVH",
              "EVD.placement", "Surgical.evacuation", "Dead_3M", "GOS", "mRS", "ICH.Score", "GCS", "FUNC",
              "Has_Infratent_Combined") #"Age_Dichot", 


#### Don't even use this anymore since dichotomizing age
## wilcox test
wilcox_out <- list()
# loop over all continuous variables
for (i in 1:length(numeric_test_vars)) {
  # and create a formula to pass on to the wilcox.test function
  form <- as.formula(paste(numeric_test_vars[i], "~", "MRI"))
  res <- wilcox.test(form, data=table1)#, conf.int=T)
  # save the output in the list with the same name as the current variable
  wilcox_out[[numeric_test_vars[i]]] <- res
  
}







### If I ever need to re=run this again, I'm going to have to replace the
### original GOS scores for the 158 patients by GOS_new



## Fisher test
fisher_out <- list()
# loop over categorical variables
for (i in 1:length(cat_test_vars)) {
  # make a different formula to pass to the xtabs function.
  # the format of this formula is a little different than the one for the wilcox test.
  # this is just because xtabs is finicky.
  form <- as.formula(paste("~", cat_test_vars[i], "+", "MRI"))
  # this creates a contingency table
  tab <- xtabs(form, data=table1)
  # which is then passed to the fisher.test function
  if (all(dim(tab) == c(2, 2))) {
    res <- fisher.test(tab)
  } else {
    res <- fisher.test(tab, hybrid=T)
  }
  #res <- fisher.test(tab)#, workspace=1e8)
  # edit the "data.name" for each result to specify what the input variables are
  # (this gets "removed" when feeding the function just a table (named "tab" in this case))
  # ...don't really even need this at all.
  res$data.name <- paste("MRI", "~", cat_test_vars[i])
  # save output in this list with same name as the current variable again
  fisher_out[[cat_test_vars[i]]] <- res
  
}


## Bonferroni correction on all the fisher test p values
# extract all p values
fisher_p <- c(Age=wilcox_out$Age$p.value, unlist(lapply(fisher_out, "[[", "p.value")), Etiology=eti_res$p.value)
p.adjust(fisher_p, method="bonferroni")



#### logistic regression ####

# make a new table with just the 158 patients and the relevant parameters
logRegressDat <- select(ichmri, MRN, MRI_Cs2, follow2)
logRegressDat <- merge(table1, logRegressDat, by="MRN", all.x=T)
logRegressDat <- logRegressDat %>%
                      filter(MRI==1)

logRegressDat$Hg_vol_mL <- logRegressDat$Hg_vol_Norm / 10000  # in "units" of 10 mL
logRegressDat$Ed_vol_mL <- logRegressDat$Ed_vol_Norm / 10000


new_numeric <- c("ICH.Score", "FUNC", "GCS", "Time_ICH_to_MRI", "Time_MRI_to_ICU_dch", "Hg_vol_mL", "Ed_vol_mL", "MLS_Norm")
new_categ <- c("Has_IVH", "Has_Infratent", "deep", "lobar")

## medians and iqr
# at MRI
med_out2 <- logRegressDat %>%
  select_(.dots=c("MRI_Cs2", new_numeric)) %>%
  group_by(MRI_Cs2) %>%
  summarise_all(funs(median, get_25th, get_75th, .args=list(na.rm=T)))

med_out2 <- med_out2[,c("MRI_Cs2", sort(names(med_out2))[-which(sort(names(med_out2))=="MRI_Cs2")])]

# at Dch
med_out2.d <- logRegressDat %>%
  select_(.dots=c("follow2", new_numeric)) %>%
  group_by(follow2) %>%
  summarise_all(funs(median, get_25th, get_75th, .args=list(na.rm=T)))

med_out2.d <- med_out2.d[,c("follow2", sort(names(med_out2.d))[-which(sort(names(med_out2.d))=="follow2")])]


## counts and freqs
var_tabs2 <- list()
for (i in 1:length(new_categ)) {
  
  assign( new_categ[i],
          
          logRegressDat %>%
            group_by_("follow2", new_categ[i]) %>%
            summarize(tot = n()) %>%
            group_by(follow2) %>%
            mutate(tot_sum = sum(tot)) %>%
            mutate(percentage = (tot / tot_sum) * 100)
  )
  
  var_tabs2[[i]] <- eval(as.name(new_categ[i]))
  names(var_tabs2)[i] <- new_categ[i]
}




# first model for assessment at MRI
glmfit <- glm(MRI_Cs2 ~ ICH.Score + FUNC + GCS + Time_ICH_to_MRI + Time_MRI_to_ICU_dch + lobar + deep +
                Has_Infratent + Hg_vol_mL + Ed_vol_mL + MLS_Norm + Has_IVH, data=logRegressDat)


summary(glmfit)
# odds ratios
exp(coef(glmfit))
exp(confint(glmfit))



# second model for assessment at discharge
glmfit2 <- glm(follow2 ~ ICH.Score + FUNC + GCS + Time_ICH_to_MRI + Time_MRI_to_ICU_dch + lobar + deep +
                 Has_Infratent + Hg_vol_mL + Ed_vol_mL + MLS_Norm + Has_IVH, data=logRegressDat)

summary(glmfit2)
exp(coef(glmfit2))
exp(confint(glmfit2))






#### Do univariate analyses ####
all_params <- c(new_numeric, new_categ)

univar_MRI <- runUnivariateRegress(logRegressDat, "MRI_Cs2", all_params)
univar_Dch <- runUnivariateRegress(logRegressDat, "follow2", all_params)

# convert estimate and CIs into Odds ratios
univar_MRI[,c(1,5,6)] <- lapply(univar_MRI[,c(1,5,6)], exp)
names(univar_MRI)[1] <- "OR"
univar_Dch[,c(1,5,6)] <- lapply(univar_Dch[,c(1,5,6)], exp)
names(univar_Dch)[1] <- "OR"

# subset the tables
univar_OR_MRI <- univar_MRI[,c("OR", "CI 2.5 %", "CI 97.5 %")]
univar_OR_Dch <- univar_Dch[,c("OR", "CI 2.5 %", "CI 97.5 %")]

#write.xlsx(univar_OR_MRI, "E:/Experiments/ICH_MRI/Univariate_logistic_regresss_MRI_v2.xlsx", rowNames=T)
#write.xlsx(univar_OR_Dch, "E:/Experiments/ICH_MRI/Univariate_logistic_regresss_Dch_v2.xlsx", rowNames=T)







#### Check n(%) using logRegressDat ####
var_tabs2 <- list()
for (i in 1:length(cat_vars)) {
  
  assign( cat_vars[i],
          
          logRegressDat %>%
            group_by_("MRI_Cs2", cat_vars[i]) %>%
            filter_(
              paste('!is.na(', cat_vars[i], ')', sep="")
            ) %>%
            summarize(tot = n()) %>%
            group_by(MRI_Cs2) %>%
            mutate(tot_sum = sum(tot)) %>%
            mutate(percentage = (tot / tot_sum) * 100)
  )
  
  var_tabs2[[i]] <- eval(as.name(cat_vars[i]))
  names(var_tabs2)[i] <- cat_vars[i]
}



