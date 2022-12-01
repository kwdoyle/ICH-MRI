library(openxlsx)
library(dplyr)


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

# calculates number of IVH cases there are by going through each patient's
# CT scans and seeing if they have a "total IVH volume"
NumberIVH <- function(MRNs, MRI=TRUE) {
  if (MRI==TRUE) {
    IDs.MRI <- filter(ichopdb$`Admission data`, MRN %in% MRNs)$Subject.ID
  } else {
    IDs.MRI <- filter(ichopdb$`Admission data`, !MRN %in% MRNs)$Subject.ID
  }
  
  
  scores <- NULL
  for (i in 1:length(IDs.MRI)) {
    tmp <- filter(ichop$`CT Scan (3-6 day)`, Subject_ID==IDs.MRI[i])
    if (nrow(tmp) != 0) {
      ivhs <- tmp$IVHS
    } else {
      ivhs <- NA
    }
    
    if (is.na(ivhs)) {
      tmp <- filter(ichop$`CT Scan (3rd in 24hr)`, Subject_ID==IDs.MRI[i])
      if (nrow(tmp) != 0) {
        ivhs <- tmp$IVH_Score_total_volume
      } else {
        ivhs <- NA
      }
      
      if (is.na(ivhs)) {
        tmp <- filter(ichopdb$`CT Scan (2nd in 24 hr)`, Subject.ID==IDs.MRI[i])
        if (nrow(tmp) != 0) {
          ivhs <- tmp$IVH.Score.total.volume
        } else {
          ivhs <- NA
        }
        
        if (is.na(ivhs)) {
          tmp <- filter(ichopdb$`CT Scan (Admission)`, Subject.ID==IDs.MRI[i])
          if (nrow(tmp) != 0) {
            ivhs <- tmp$IVH.Score.total.volume 
          } else {
            ivhs <- NA
          }
          
          if (is.na(ivhs)) {
            tmp <- filter(ichopdb$`CT Scan (OSH)`, Subject.ID==IDs.MRI[i])
            if (nrow(tmp) != 0) {
              ivhs <- tmp$IVH.Score.total.volume 
            } else {
              ivhs <- NA
            }
            
            if (is.na(ivhs)) {
              print(paste("No IVHS for ICHOP ID", IDs.MRI[i]))
            }
          } else {
            scores <- c(scores, ivhs)
          }
        } else {
          scores <- c(scores, ivhs)
        }
      } else {
        scores <- c(scores, ivhs)
      }
    } else {
      scores <- c(scores, ivhs)
    }
  }
  
  
  # how many scores aren't 0
  n_IVH <- length(which(scores != 0))
  pcnt <- (n_IVH / length(scores)) * 100
  # how many are missing
  miss <- length(IDs.MRI) - length(scores)
  
  return(list(n_IVH=n_IVH, percentage=pcnt, n_missing=miss))
  
}



NumberIVH2 <- function(MRNs, MRI=TRUE) {
  if (MRI==TRUE) {
    IDs.MRI <- filter(ichopdb$`Admission data`, MRN %in% MRNs)$Subject.ID
  } else {
    IDs.MRI <- filter(ichopdb$`Admission data`, !MRN %in% MRNs)$Subject.ID
  }
  
  
  scores <- matrix(nrow=length(IDs.MRI), ncol=2)
  for (i in 1:length(IDs.MRI)) {
    curr_mrn <- filter(ichopdb$`Admission data`, Subject.ID == IDs.MRI[i])$MRN
    scores[i,1] <- curr_mrn
    
    tmp <- filter(ichop$`CT Scan (3-6 day)`, Subject_ID==IDs.MRI[i])
    if ()
    if (nrow(tmp) != 0) {
      scores[i,2] <- tmp$IVHS
    } else {
      tmp <- filter(ichop$`CT Scan (3rd in 24hr)`, Subject_ID==IDs.MRI[i])
      if (nrow(tmp) != 0) {
        scores[i,2] <- tmp$IVH_Score_total_volume
      } else {
        tmp <- filter(ichopdb$`CT Scan (2nd in 24 hr)`, Subject.ID==IDs.MRI[i])
        if (nrow(tmp) != 0) {
          scores[i,2] <- tmp$IVH.Score.total.volume
        } else {
          tmp <- filter(ichopdb$`CT Scan (Admission)`, Subject.ID==IDs.MRI[i])
          if (nrow(tmp) != 0) {
            scores[i,2] <- tmp$IVH.Score.total.volume 
          } else {
            tmp <- filter(ichopdb$`CT Scan (OSH)`, Subject.ID==IDs.MRI[i])
            if (nrow(tmp) != 0) {
              scores[i,2] <- tmp$IVH.Score.total.volume 
            } else {
              print(paste("No IVHS for ICHOP ID", IDs.MRI[i]))
            }
          }
        }
      }
    }
    

  }
  
  
  # # how many scores aren't 0
  # n_IVH <- length(which(scores != 0))
  # pcnt <- (n_IVH / length(scores)) * 100
  # # how many are missing
  # miss <- length(IDs.MRI) - length(scores)
  scores <- as.data.frame(scores)
  colnames(scores) <- c("MRN", "score")
  
  return(scores)
  
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





# find Hemmorage volumes for patients
# ...I don't think I need this?
FindHemmorageVolumesBAD <- function(MRNs, MRI=TRUE) {
  if (MRI==TRUE) {
    IDs.MRI <- filter(ichopdb$`Admission data`, MRN %in% MRNs)$Subject.ID
  } else {
    IDs.MRI <- filter(ichopdb$`Admission data`, !MRN %in% MRNs)$Subject.ID
  }
  
  
  scores <- data.frame()
  for (i in 1:length(IDs.MRI)) {
    tmp <- filter(ichopdb$`CT Scan (3-6 day)`, Subject.ID==IDs.MRI[i])
    if (nrow(tmp) != 0) {
      ivhs <- c(IDs.MRI[i], tmp$Hematoma.volume)
    } else {
      ivhs <- NA
    }
    
    if (is.na(ivhs)) {
      tmp <- filter(ichopdb$`CT Scan (3rd in 24hr)`, Subject.ID==IDs.MRI[i])
      if (nrow(tmp) != 0) {
        ivhs <- tmp$Hematoma.volume
      } else {
        ivhs <- NA
      }
      
      if (is.na(ivhs)) {
        tmp <- filter(ichopdb$`CT Scan (2nd in 24 hr)`, Subject.ID==IDs.MRI[i])
        if (nrow(tmp) != 0) {
          ivhs <- tmp$Hematoma.volume
        } else {
          ivhs <- NA
        }
        
        if (is.na(ivhs)) {
          tmp <- filter(ichopdb$`CT Scan (Admission)`, Subject.ID==IDs.MRI[i])
          if (nrow(tmp) != 0) {
            ivhs <- tmp$Hematoma.volume 
          } else {
            ivhs <- NA
          }
          
          if (is.na(ivhs)) {
            tmp <- filter(ichopdb$`CT Scan (OSH)`, Subject.ID==IDs.MRI[i])
            if (nrow(tmp) != 0) {
              ivhs <- tmp$Hematoma.volume 
            } else {
              ivhs <- NA
            }
            
            if (is.na(ivhs)) {
              print(paste("No hematoma volumes for ICHOP ID", IDs.MRI[i]))
            }
          } else {
            scores <- c(scores, ivhs)
          }
        } else {
          scores <- c(scores, ivhs)
        }
      } else {
        scores <- c(scores, ivhs)
      }
    } else {
      scores <- c(scores, ivhs)
    }
  }
  
  
  # how many scores aren't 0
  n_IVH <- length(which(scores != 0))
  pcnt <- (n_IVH / length(scores)) * 100
  # how many are missing
  miss <- length(IDs.MRI) - length(scores)
  
  return(list(n_IVH=n_IVH, percentage=pcnt, n_missing=miss))
  
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


# calculate scores from the ICHOP database (not actually useful)
ICHOPscoreCalc <- function(MRNs, param) {
  filtab <- filter(ichopdb$`Admission data`, MRN %in% MRNs)
  
  md <- median(filtab[,param], na.rm=T)
  qtl <- quantile(filtab[,param], na.rm=T)
  
  return(c(median=md, quantile=qtl))
}

# calculate numbers of lobar and deep
lobardeepNumbers <- function(MRNs, MRI=TRUE) {
  if (MRI==TRUE) {
    tab <- ichopdb$Etiology %>%
      left_join(select(ichopdb$`Admission data`, Subject.ID, MRN), by="Subject.ID") %>%
      filter(MRN %in% MRNs) %>%
      mutate(lobar = ifelse(grepl("1", Location2), yes=1, 
                            no=ifelse(grepl("6", Location2), yes=1, no=0)),
             deep = ifelse(grepl("2", Location2), yes=1,
                           no=ifelse(grepl("3", Location2), yes=1, no=0))) %>%
      select(MRN, Location2, lobar, deep)
  } else {
    tab <- ichopdb$Etiology %>%
      left_join(select(ichopdb$`Admission data`, Subject.ID, MRN), by="Subject.ID") %>%
      filter(!MRN %in% MRNs) %>%
      mutate(lobar = ifelse(grepl("1", Location2), yes=1, 
                            no=ifelse(grepl("6", Location2), yes=1, no=0)),
             deep = ifelse(grepl("2", Location2), yes=1,
                           no=ifelse(grepl("3", Location2), yes=1, no=0))) %>%
      select(MRN, Location2, lobar, deep)
  }
  n_lobar <- table(tab$lobar)[2]
  pcnt_lobar <- n_lobar / sum(table(tab$lobar)) * 100
  
  n_deep <- table(tab$deep)[2]
  pcnt_deep <- n_deep / sum(table(tab$deep)) * 100
  
  out <- list(n_lobar=n_lobar, pcnt_lobar=pcnt_lobar, n_deep=n_deep, pcnt_deep=pcnt_deep)
  return(out)
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








# list of 158 patients who had MRI
ichmri <- read.xlsx("E:/Experiments/ICH_MRI/MRI_Merged_Data_158patients.xlsx")
# don't need to subset unique MRNs, since each patient has their own row
length(unique(ichmri$MRN))==nrow(ichmri)

# ichop database
ichopdb <- readDB("E:/SHOP & ICHOP forms/ICHOP_1_17_18_database_full.xlsx")



# set dates
ichmri$MRI_date_Loc <- convertToDate(ichmri$MRI_date_Loc)
ichmri$mri_date <- convertToDate(ichmri$mri_date)
ichmri$Discharge_date <- convertToDate(ichmri$Discharge_date)
# apparently there's one missing value in MRI_date_Loc, so use mri_date instead.


# set missing ed vols to 0 like in the analysis script
ichmri$Ed_vol[which(is.na(ichmri$Ed_vol))] <- 0




#### GCS on Admission ####

# ok... this is stupid... there were 160 total patients that the SQL query was looking for
# (two of them not being in the list of 158 patients, obvsiously)
# BUT, there wasn't any data for 2 OTHER patients, yeilding 158 patients in this table,
# causing THESE 158 patients to be different than the normal 158 patients.
# so need to filter this for the patients in the study.
# (will be less than 158 in this table then.)
# GCSdat <- read.xlsx("E:/Experiments/ICH_MRI/Demograpics/GCS_ON_ADMISSION_DATE.xlsx")
# GCSdat$MRI_DTM <- convertToDateTime(GCSdat$MRI_DTM)
# GCSdat$ADMIT_DTM <- convertToDateTime(GCSdat$ADMIT_DTM)
# GCSdat$GCS_RECORDDTM <- convertToDateTime(GCSdat$GCS_RECORDDTM)
# # "null" is a value for the GCS score. *sigh*...
# GCSdat$OBSVALUE <- as.numeric(GCSdat$OBSVALUE)
# # filter for the 158 patients
# GCSdat<- GCSdat %>%
#   filter(PATIENT_ID %in% ichmri$MRN)
# 
# # filter for everyone else
# GCSdat_NoMRI<- GCSdat %>%
#   filter(!PATIENT_ID %in% ichmri$MRN)


# try making big conglomerate table
table1 <- ichopdb$`Admission data`[,c("MRN"), drop=F]
table1$MRI <- ifelse(table1$MRN %in% ichmri$MRN, 1, 0)




# oh hot damn, the GCS scores are in the ichop database too
# and nice, they match with the stuff I pulled
# MRI
summary(as.numeric(filter(ichopdb$`Admission data`, MRN %in% ichmri$MRN)$`GCS.(CUMC)`))
# no MRI
summary(as.numeric(filter(ichopdb$`Admission data`, !MRN %in% ichmri$MRN)$`GCS.(CUMC)`))

wilcox.test(as.numeric(filter(ichopdb$`Admission data`, MRN %in% ichmri$MRN)$`GCS.(CUMC)`),
            as.numeric(filter(ichopdb$`Admission data`, !MRN %in% ichmri$MRN)$`GCS.(CUMC)`))

table1 <- merge(table1, select(ichopdb$`Admission data`, MRN, `GCS.(CUMC)`), by="MRN", all.x=T)





## Need to select just the OBSVALUE from the minimum GCS_RECORDDTM for each patient
# well great. if a value was forgotten to be filled in for the min RECORDEDDTM,
# then this chooses the NA for that patient.
# GCSdat.clean <- GCSdat %>%
#   group_by(PATIENT_ID, MRI_DTM, ADMIT_DTM) %>%
#   filter(GCS_RECORDDTM == min(GCS_RECORDDTM, na.rm=T))


# try ..imputing the missing data first?
# for each patient, fill in any missing values with the average GCS score
# (if a patient has no scores, i.e., is just a single row, then it will just put back NA)
GCSdat.imput <- GCSdat %>%
  group_by(PATIENT_ID) %>%
  mutate(OBSVALUE = ifelse(is.na(OBSVALUE), mean(OBSVALUE, na.rm=T), OBSVALUE))


# now get first value from day of admission closest to the MRI assessment
GCSdat.clean <- GCSdat.imput %>%
  group_by(PATIENT_ID, MRI_DTM, ADMIT_DTM) %>%
  filter(GCS_RECORDDTM == min(GCS_RECORDDTM, na.rm=T) | is.na(GCS_RECORDDTM)) %>%
  distinct()
  # this supposedly would keep NAs, but it doesn't.
  # can forcibly keep all NAs though with | is.na(GCS_RECORDEDDTM)
  #filter(GCS_RECORDDTM %in% min(GCS_RECORDDTM, na.rm=T))

median(GCSdat.clean$OBSVALUE, na.rm=T)
quantile(GCSdat.clean$OBSVALUE, na.rm=T)




#### Age, Sex, Ethnicity ####
# ASE <- read.xlsx("E:/Experiments/ICH_MRI/Demograpics/AGE_SEX_ETHNICITY.xlsx")
# # make column of indicator for who had MRI or not
# ## ,,,why 'tf did I make this column?? Just filter patients in here for same MRNs as in ichmri$MRN.
# #ASE$Had_MRI <- ifelse(ASE$PATIENT_ID %in% GCSdat$PATIENT_ID, 1, 0)
# ASE.MRI <- ASE %>%
#   filter(PATIENT_ID %in% ichmri$MRN)
# 
# # # select the 158 patients who had MRI
# # # NOPE this was wrong, because the original 158 patients in GCSdat wasn't the same 158 from ichmri
# # ASE.MRI <- ASE %>%
# #   filter(PATIENT_ID %in% unique(GCSdat$PATIENT_ID))
# 
# ## summary stats
# # mean age
# # mean(ASE.MRI$AGE)
# # sd(ASE.MRI$AGE)
# summary(ASE.MRI$AGE)

### Oooooooomg. the ages were in the ichop database this whole time.
### I didn't need to get them from the CDW after all.
### But at least the numbers match up from the two data pulls.
summary(as.numeric(filter(ichopdb$`Admission data`, MRN %in% ichmri$MRN)$Age))

## for No MRI
summary(as.numeric(filter(ichopdb$`Admission data`, !MRN %in% ichmri$MRN)$Age))


wilcox.test(as.numeric(filter(ichopdb$`Admission data`, MRN %in% ichmri$MRN)$Age),
            as.numeric(filter(ichopdb$`Admission data`, !MRN %in% ichmri$MRN)$Age))


table1 <- merge(table1, select(ichopdb$`Admission data`, MRN, Age), my="MRN", all.x=T)

# % female
# there's a unique row for each patient, so can just table the counts
# length(unique(ASE.MRI$PATIENT_ID))
# nrow(ASE.MRI)
# tot <- nrow(ASE.MRI)
# # could just index the first one, but this way is more """"robust""""
# n_female <- table(ASE.MRI$SEX_CODE)[names(table(ASE.MRI$SEX_CODE))=="F"]
# (n_female / tot) * 100
# 
# # % white
# n_white <- table(ASE.MRI$RACE_CODE)[names(table(ASE.MRI$RACE_CODE))=="W"]
# (n_white / tot) * 100


# this matches too. I guess 2 refers to female
# MRI
sex_MRI <- table(filter(ichopdb$`Admission data`, MRN %in% ichmri$MRN)$Sex)
sex_MRI[2] / sum(sex_MRI)
# no MRI
sex_noMRI <- table(filter(ichopdb$`Admission data`, !MRN %in% ichmri$MRN)$Sex)
sex_noMRI[2] / sum(sex_noMRI)

# this one doesn't. it's off by 2. I guess I should use this instead though
# MRI
eth_MRI <- table(filter(ichopdb$`Admission data`, MRN %in% ichmri$MRN)$Ethnicity)
eth_MRI[1] / sum(eth_MRI)
# no MRI
eth_noMRI <- table(filter(ichopdb$`Admission data`, !MRN %in% ichmri$MRN)$Ethnicity)
eth_noMRI[1] / sum(eth_noMRI)


table1 <- merge(table1, select(ichopdb$`Admission data`, MRN, Sex, Ethnicity), by="MRN", all.x=T)


# test for male/female different ratios
fisher.test(as.table(matrix(c(sex_MRI[2], sex_MRI[1],  
                              sex_noMRI[2], sex_noMRI[1]), nrow=2, ncol=2, byrow = F)))

# test for ethnicity ratios




### Etiology ####
# dch_table <- merge(ichopdb$`Admission data`[,c("Subject.ID", "MRN")], ichopdb$Discharge, by="Subject.ID", all.y=TRUE)
# dch_MRI_dates3 <- merge(MRI_dates, dch_table[,c("MRN", "Discharge.Date"), drop=FALSE], by="MRN", all.x=TRUE)
etio_table <- merge(ichopdb$`Admission data`[,c("Subject.ID", "MRN")], ichopdb$Etiology, by="Subject.ID", all.y=T)
#etio_table$Etiology <- factor(etio_table$Etiology)

table1 <- merge(table1, select(etio_table, MRN, Etiology), by="MRN", all.x=T)

etio_MRI <- etio_table %>%
  filter(MRN %in% ichmri$MRN) %>%
  select(MRN, Etiology)

etio_NotMRI <- etio_table %>%
  filter(!MRN %in% ichmri$MRN) %>%
  select(MRN, Etiology)

# car::recode(etio_table$Etiology, "'1'='HTN'; '2'='Amyloid'; '4'='AVM'; '5'='Aneurysm'; '7'='Coagulopathy+Anticoag'; '8'='Drug Induced';
#        '9'='Coagulopathy+Anticoag'; '0'='Other'; '11'='Unknown';
#             '1, 8'='HTN'; '1, 4'='HTN'; '1, 7'='HTN'; '0, 1'='HTN'; '1, 9'='HTN'; ") # nope can't recode the doubles becasue there's one that has 1 and 2
# # and we're interested in 1 and 2 seperately

# get counts of each etiology
table(etio_MRI$Etiology)
# as percentages
pcnt_MRI <- sort((table(etio_MRI$Etiology) / sum(table(etio_MRI$Etiology))) * 100, decreasing=T)


# get counts of each etiology
table(etio_NotMRI$Etiology)
# as percentages
pcnt_NotMRI <- sort((table(etio_NotMRI$Etiology) / sum(table(etio_NotMRI$Etiology))) * 100, decreasing=T)



### HTN
## MRI
bad_htn_ix1 <- grep(11, names(pcnt_MRI))
# indices with 1
good_htn_ix1 <- grep(1, names(pcnt_MRI))
# 1s with 11 excluded
use_htn_ix1 <- good_htn_ix1[which(!good_htn_ix1 %in% bad_htn_ix1)]
names(pcnt_MRI)[use_htn_ix1]
# sum all percents that contain "1" to get total percent of patients with "1"
sum(pcnt_MRI[use_htn_ix1])

## No MRI
# indices with 11
bad_htn_ix2 <- grep(11, names(pcnt_NotMRI))
# indices with 1
good_htn_ix2 <- grep(1, names(pcnt_NotMRI))
# 1s with 11 excluded
use_htn_ix2 <- good_htn_ix2[which(!good_htn_ix2 %in% bad_htn_ix2)]
names(pcnt_NotMRI)[use_htn_ix2]
# sum all percents that contain "1" to get total percent of patients with "1"
sum(pcnt_NotMRI[use_htn_ix2])




### Amyloid
## MRI
amyl_ix1 <- grep(2, names(pcnt_MRI))
sum(pcnt_MRI[amyl_ix1])

## No MRI
amyl_ix2 <- grep(2, names(pcnt_NotMRI))
sum(pcnt_NotMRI[amyl_ix2])




### Coagulopathy + Anticoag
## MRI
coag_ix1 <- grep(7, names(pcnt_MRI))
anticoag_ix1 <- grep(9, names(pcnt_MRI))
both_ix1 <- c(coag_ix1, anticoag_ix1)

sum(pcnt_MRI[both_ix1])

# No MRI
coag_ix2 <- grep(7, names(pcnt_NotMRI))
anticoag_ix2 <- grep(9, names(pcnt_NotMRI))
both_ix2 <- c(coag_ix2, anticoag_ix2)

sum(pcnt_NotMRI[both_ix2])




### Other
### This is everything BESIDES the ones above
## MRI
# this winds up just being "Unknown" for the MRI ones
sum(pcnt_MRI[-c(use_htn_ix1, amyl_ix1, both_ix1)])

## No MRI
sum(pcnt_NotMRI[-c(use_htn_ix2, amyl_ix2, both_ix2)])






#### Outcomes ####
# outcms <- read.xlsx("E:/Experiments/ICH_MRI/ICH_MRI_Patient_Outcomes.xlsx")
# 
# ## Median mRS at 3mo fu
# median(outcms$mRS_3M, na.rm=T)
# quantile(outcms$mRS_3M, na.rm=T)


# ## number dead
# ndead <- length(which(outcms$mRS_3M == 6 & !is.na(outcms$mRS_3M)))
# (ndead / tot) * 100

MRI_out <- ichopdb$`3 month FU` %>%
  filter(MRN %in% ichmri$MRN) %>%
  mutate(GOS = as.numeric(GOS)) %>%
  select(MRN, mRS, GOS)


noMRI_out <- ichopdb$`3 month FU` %>%
  filter(!MRN %in% ichmri$MRN) %>%
  mutate(GOS = as.numeric(GOS)) %>%
  select(MRN, mRS, GOS)


summary(MRI_out$mRS)
summary(MRI_out$GOS)

summary(noMRI_out$mRS)
summary(noMRI_out$GOS)


table1 <- merge(table1, select(ichopdb$`3 month FU`, MRN, mRS, GOS), by="MRN", all.x=T)


# ## use ichop db instead
# # MRI
# summary(filter(ichopdb$`3 month FU`, MRN %in% ichmri$MRN)$mRS)
# # no MRI
# summary(filter(ichopdb$`3 month FU`, !MRN %in% ichmri$MRN)$mRS)


## MRI
# number dead
# dead_MRI <- table(filter(ichopdb$`3 month FU`, MRN %in% ichmri$MRN)$mRS)
# dead_MRI[names(dead_MRI)==6]
# # percent
# dead_MRI[names(dead_MRI)==6] / sum(dead_MRI) * 100

length(which(MRI_out$mRS==6))
length(which(MRI_out$mRS==6)) / sum(table(MRI_out$mRS)) * 100

## no MRI
# number dead
# dead_noMRI <- table(filter(ichopdb$`3 month FU`, !MRN %in% ichmri$MRN)$mRS)
# dead_noMRI[names(dead_noMRI)==6]
# # percent
# dead_noMRI[names(dead_noMRI)==6] / sum(dead_noMRI) * 100

length(which(noMRI_out$mRS==6))
length(which(noMRI_out$mRS==6)) / sum(table(noMRI_out$mRS)) * 100






## MRI
# number dead
dead_MRI <- table(filter(ichopdb$`3 month FU`, MRN %in% ichmri$MRN)$GOS)
dead_MRI[names(dead_MRI)==6]
# percent
dead_MRI[names(dead_MRI)==6] / sum(dead_MRI) * 100

## no MRI
# number dead
dead_noMRI <- table(filter(ichopdb$`3 month FU`, !MRN %in% ichmri$MRN)$GOS)
dead_noMRI[names(dead_noMRI)==6]
# percent
dead_noMRI[names(dead_noMRI)==6] / sum(dead_noMRI) * 100







#### Data from ICHOP database (ICH score, FUNC score, etc.) ####


# ichopdb$`Admission data`$MRN <- as.numeric(ichopdb$`Admission data`$MRN)
# ichopdb$`Admission data`$ICH.Score <- as.numeric(ichopdb$`Admission data`$ICH.Score)
# ichopdb$`Admission data`$FUNC <- as.numeric(ichopdb$`Admission data`$FUNC)
# 
# # extract the 158 patients
# ichop.adm.MRI <- filter(ichopdb$`Admission data`, MRN %in% ichmri$MRN)
# 
# ## ICH score
# median(ichop.adm.MRI$ICH.Score, na.rm=T)
# quantile(ichop.adm.MRI$ICH.Score, na.rm=T)


MRI_scores <- ichopdb$`Admission data` %>%
  filter(MRN %in% ichmri$MRN) %>%
  select(MRN, FUNC, ICH.Score)

noMRI_scores <- ichopdb$`Admission data` %>%
  filter(!MRN %in% ichmri$MRN) %>%
  select(MRN, FUNC, ICH.Score)


## FUNC score
# median(ichop.adm.MRI$FUNC, na.rm=T)
# quantile(ichop.adm.MRI$FUNC, na.rm=T)
summary(MRI_scores$FUNC)
summary(noMRI_scores$FUNC)

## ICH score
summary(MRI_scores$ICH.Score)
summary(noMRI_scores$ICH.Score)


table1 <- merge(table1, select(ichopdb$`Admission data`, MRN, FUNC, ICH.Score), by="MRN", all.x=T)



#ICHOPscoreCalc(ichmri$MRN, "ICH.Score")








#### ICH volume ####
# volumes are only normalized in the stats script (ICH_MRI_stats_vX), so these are the "raw" values.
# I also think these are in microliters.
# Need to seperate volumns into groups of: <30mL, 30-60mL, >60mL.

# convert to mL
ichmri$Hg_vol_mL <- ichmri$Hg_vol / 1000
# less30 <- ichmri$Hg_vol_mL[ichmri$Hg_vol_mL < 30]
# btwn30and60 <- ichmri$Hg_vol_mL[ichmri$Hg_vol_mL <= 60 & ichmri$Hg_vol_mL >= 30]
# great60 <- ichmri$Hg_vol_mL[ichmri$Hg_vol_mL > 60]

less30 <- ichmri %>%
  filter(Hg_vol_mL < 30) %>%
  select(MRN, Hg_vol_mL)

length(unique(less30$MRN))
length(unique(less30$MRN)) / nrow(ichmri) * 100


btwn30and60 <- ichmri %>%
  filter(Hg_vol_mL <= 60 & Hg_vol_mL >= 30) %>%
  select(MRN, Hg_vol_mL)

length(unique(btwn30and60$MRN))
length(unique(btwn30and60$MRN)) / nrow(ichmri) * 100

great60 <- ichmri %>%
  filter(Hg_vol_mL > 60) %>%
  select(MRN, Hg_vol_mL)

length(unique(great60$MRN))
length(unique(great60$MRN)) / nrow(ichmri) * 100



#a <- merge(table1, select(ichmri, MRN, Hg_vol_mL), by="MRN", all.x=T)


# median(less30)
# quantile(less30)
# 
# median(btwn30and60)
# quantile(btwn30and60)
# 
# median(great60)
# quantile(great60)




#### CT Volumes instead of ICH volumes (so we can get volumes for everybody and not just those with MRI) ####
MRI_CT_vol <- FindHemmorageVolumes(ichmri$MRN, MRI=TRUE)
noMRI_CT_vol <- FindHemmorageVolumes(ichmri$MRN, MRI=FALSE)

MRI_CT_vol$Hematoma.volume <- as.numeric(MRI_CT_vol$Hematoma.volume)
noMRI_CT_vol$Hematoma.volume <- as.numeric(noMRI_CT_vol$Hematoma.volume)



CT_vol <- rbind(MRI_CT_vol, noMRI_CT_vol)

table1 <- merge(table1, CT_vol, by="MRN", all.x=T)



# MRI
n_30_MRI <- length(which(MRI_CT_vol$Hematoma.volume < 30))
n_30_60_MRI <- length(which(MRI_CT_vol$Hematoma.volume >= 30 & MRI_CT_vol$Hematoma.volume <= 60))
n_60_MRI <- length(which(MRI_CT_vol$Hematoma.volume > 60))
tot_MRI <- length(which(!is.na(MRI_CT_vol$Hematoma.volume)))

print("less than 30 mL - MRI")
n_30_MRI
n_30_MRI / tot_MRI * 100
print("30 to 60 mL - MRI")
n_30_60_MRI
n_30_60_MRI / tot_MRI * 100
print("greater than 60 mL - MRI")
n_60_MRI
n_60_MRI / tot_MRI * 100




# no MRI
n_30_noMRI <- length(which(noMRI_CT_vol$Hematoma.volume < 30))
n_30_60_noMRI <- length(which(noMRI_CT_vol$Hematoma.volume >= 30 & noMRI_CT_vol$Hematoma.volume <= 60))
n_60_noMRI <- length(which(noMRI_CT_vol$Hematoma.volume > 60))
tot_noMRI <- length(which(!is.na(noMRI_CT_vol$Hematoma.volume)))

print("less than 30 mL - noMRI")
n_30_noMRI
n_30_noMRI / tot_noMRI * 100
print("30 to 60 mL - noMRI")
n_30_60_noMRI
n_30_60_noMRI / tot_noMRI * 100
print("greater than 60 mL - noMRI")
n_60_noMRI
n_60_noMRI / tot_noMRI * 100







#### Infratentorial number ####
n_infra <- length(which(ichmri$infra_tent != "None"))
(n_infra / tot) * 100

table1 <- merge(table1, select(ichmri, MRN, infra_tent), by="MRN", all.x=T)



#### lobar/deep number ####
# # MRI
# lobardeep_MRI <- ichopdb$Etiology %>%
#   left_join(select(ichopdb$`Admission data`, Subject.ID, MRN), by="Subject.ID") %>%
#   filter(MRN %in% ichmri$MRN) %>%
#   mutate(lobar = ifelse(grepl("1", Location2), yes=1, 
#                         no=ifelse(grepl("6", Location2), yes=1, no=0)),
#          deep = ifelse(grepl("2", Location2), yes=1,
#                        no=ifelse(grepl("3", Location2), yes=1, no=0))) %>%
#   select(MRN, Location2, lobar, deep)
# 
# 
# n_lobar <- table(lobardeep_MRI$lobar)[2]
# n_lobar
# n_lobar / sum(table(lobardeep_MRI$lobar)) * 100
# 
# n_deep <- table(lobardeep_MRI$deep)[2]
# n_deep
# n_deep / sum(table(lobardeep_MRI$deep)) * 100
# 
# 
# 
# # no MRI
# lobardeep_noMRI <- ichopdb$Etiology %>%
#   left_join(select(ichopdb$`Admission data`, Subject.ID, MRN), by="Subject.ID") %>%
#   filter(!MRN %in% ichmri$MRN) %>%
#   mutate(lobar = ifelse(grepl("1", Location2), yes=1, 
#                         no=ifelse(grepl("6", Location2), yes=1, no=0)),
#          deep = ifelse(grepl("2", Location2), yes=1,
#                        no=ifelse(grepl("3", Location2), yes=1, no=0))) %>%
#   select(MRN, Location2, lobar, deep)




# n_lobar2 <- table(lobardeep_noMRI$lobar)[2]
# n_lobar2
# n_lobar2 / sum(table(lobardeep_noMRI$lobar)) * 100
# 
# n_deep2 <- table(lobardeep_noMRI$deep)[2]
# n_deep2
# n_deep2 / sum(table(lobardeep_noMRI$deep)) * 100


# MRI
lobardeepNumbers(ichmri$MRN)
# no MRI
lobardeepNumbers(ichmri$MRN, MRI=FALSE)

lobardeep_table1 <- GetLobarDeep(unique(table1$MRN))
table1 <- merge(table1, lobardeep_table1, by="MRN", all.x=T)








#### number of EVDs and clot evacuation ####
# MRI
ichop.proc.MRI <- ichopdb$Procedures %>%
  left_join(select(ichopdb$`Admission data`, Subject.ID, MRN), by="Subject.ID") %>%
  filter(MRN %in% ichmri$MRN) %>%
  select(MRN, EVD.placement, Surgical.evacuation) %>%
  distinct()

n_EVD <- length(which(ichop.proc.MRI$EVD.placement))
n_EVD
n_EVD / nrow(ichop.proc.MRI) * 100

n_evac <- length(which(ichop.proc.MRI$Surgical.evacuation))
n_evac
n_evac / nrow(ichop.proc.MRI) * 100


# no MRI
ichop.proc.noMRI <- ichopdb$Procedures %>%
  left_join(select(ichopdb$`Admission data`, Subject.ID, MRN), by="Subject.ID") %>%
  filter(!MRN %in% ichmri$MRN) %>%
  select(MRN, EVD.placement, Surgical.evacuation) %>%
  distinct()

n_EVD2 <- length(which(ichop.proc.noMRI$EVD.placement))
n_EVD2
n_EVD2 / nrow(ichop.proc.noMRI) * 100

n_evac2 <- length(which(ichop.proc.noMRI$Surgical.evacuation))
n_evac2
n_evac2 / nrow(ichop.proc.noMRI) * 100




procedures <- rbind(ichop.proc.MRI, ichop.proc.noMRI)

table1 <- merge(table1, procedures, by="MRN", all.x=T)






#### Number IVH ####
## this is unnecessarily complicated.
## Going to check the "IVHS" for patients in the 3-6 day CT scan and,
## if they're not in that table, then check the 3rd in 24hr table.
## if not in that one, check the 2nd in 24 hr table, and so forth.
# MRI
NumberIVH(ichmri$MRN, MRI=TRUE)
# no MRI
NumberIVH(ichmri$MRN, MRI=FALSE)



IVH_out <- FindIVHScores(unique(table1[,1,drop=F]))

table1 <- merge(table1, select(IVH_out, MRN, All_Scores), by="MRN", all.x=T)





#### Days in hospital/ICU ####
new_CS2 <- filter(new_CS, MRN %in% ichmri$MRN)

table1 <- merge(table1, select(new_CS2, MRN, Day_of_ICU_Dch, Day_of_Hosp_Dch), by="MRN", all.x=T)


#### Put everything to make table1 into a new script ####




#### For table 2 ####















# load data to check against
# ...I realized this is the same data I now open above as "ichmri".
#dataaa <- read.xlsx("E:/Experiments/ICH_MRI/MRI_Merged_Data_158patients.xlsx")
# patients conscious at MRI
c_MRI <- ichmri %>%
  filter(MRI_Cs2 == 1) %>%
  select(MRN)
# unconscious at MRI
uc_MRI <- ichmri %>%
  filter(MRI_Cs2 == 0) %>%
  select(MRN)

# conscious at discharge
c_Dch <- ichmri %>%
  filter(follow2 == 1) %>%
  select(MRN)
# unconscious at dch
uc_Dch <- ichmri %>%
  filter(follow2 == 0) %>%
  select(MRN)





#### GCS scores for the above groups ####
# missing data for 2 patients here
GCS_c_MRI <- GCSdat.clean %>%
  filter(PATIENT_ID %in% c_MRI$MRN)
median(GCS_c_MRI$OBSVALUE, na.rm=T)
quantile(GCS_c_MRI$OBSVALUE, na.rm=T)
summary(GCS_c_MRI$OBSVALUE)

length(unique(GCS_c_MRI$PATIENT_ID))


GCS_uc_MRI <- GCSdat.clean %>%
  filter(PATIENT_ID %in% uc_MRI$MRN)
median(GCS_uc_MRI$OBSVALUE, na.rm=T)
quantile(GCS_uc_MRI$OBSVALUE, na.rm=T)

length(unique(GCS_uc_MRI$PATIENT_ID))


GCS_c_Dch <- GCSdat.clean %>%
  filter(PATIENT_ID %in% c_Dch$MRN)
median(GCS_c_Dch$OBSVALUE, na.rm=T)
quantile(GCS_c_Dch$OBSVALUE, na.rm=T)

length(unique(GCS_c_Dch$PATIENT_ID))


GCS_uc_Dch <- GCSdat.clean %>%
  filter(PATIENT_ID %in% uc_Dch$MRN)
median(GCS_uc_Dch$OBSVALUE, na.rm=T)
quantile(GCS_uc_Dch$OBSVALUE, na.rm=T)

length(unique(GCS_uc_Dch$PATIENT_ID))



#### Infraentorial number ####
## At time of MRI
# Conscious
MRI_c_infra <- ichmri %>%
  filter(MRI_Cs2 == 1 & infra_tent != "None") %>%
  select(MRN) %>%
  distinct() %>%
  nrow()

(MRI_c_infra / nrow(c_MRI)) * 100

# Unconscious
MRI_uc_infra <- ichmri %>%
  filter(MRI_Cs2 == 0 & infra_tent != "None") %>%
  select(MRN) %>%
  distinct() %>%
  nrow()

(MRI_uc_infra / nrow(uc_MRI)) * 100


## At Discharge
# Conscious
Dch_c_infra <- dat %>%
  filter(follow2 == 1 & infra_tent != "None") %>%
  select(MRN) %>%
  distinct() %>%
  nrow()

(Dch_c_infra / nrow(c_Dch)) * 100

# Unconscious
Dch_uc_infra <- dat %>%
  filter(follow2 == 0 & infra_tent != "None") %>%
  select(MRN) %>%
  distinct() %>%
  nrow()

(Dch_uc_infra / nrow(uc_Dch)) * 100











#### IVH ####
NumberIVH(c_MRI$MRN)
NumberIVH(uc_MRI$MRN)

NumberIVH(c_Dch$MRN)
NumberIVH(uc_Dch$MRN)


# could use this for the above too I guess
# aa <- FindIVHScores(ichmri[,1, drop=FALSE], MRI=TRUE)
# length(which(is.na(aa$All_Scores)))




#### ICH score ####
ICHOPscoreCalc(c_MRI$MRN, "ICH.Score")
ICHOPscoreCalc(uc_MRI$MRN, "ICH.Score")

ICHOPscoreCalc(c_Dch$MRN, "ICH.Score")
ICHOPscoreCalc(uc_Dch$MRN, "ICH.Score")


ICHOPscoreCalc(c_MRI$MRN, "FUNC")
ICHOPscoreCalc(uc_MRI$MRN, "FUNC")

ICHOPscoreCalc(c_Dch$MRN, "FUNC")
ICHOPscoreCalc(uc_Dch$MRN, "FUNC")





#### Delay between ICH and MRI/Dch ####
# use dch dates from this table:
# ok here is the table.
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
# convert MRNs to numeric
dat$MRN <- as.numeric(dat$MRN)




# get dch dates for the 158 patients from ichopdb
onset_dates <- ichopdb$`Admission data` %>%
  filter(MRN %in% unique(ichmri$MRN)) %>%
  select(MRN, onset.date) %>%
  distinct()

# why do THREE dates that are formatted the SAME way in Excel read in as numbers
# while EVERYTHING ELSE is a date-character-string???
# ...whatever, I think I can get around it by searching for those that are numbers
# and assigning those differently.
# NOPE I can't reassign them because the classes would be different.

# get indices
dumb_ix <- which(!is.na(as.numeric(onset_dates$onset.date)))
#onset_dates2$onset.date[dumb_ix] <- gsub("-", "/", as.character(convertToDate(onset_dates2$onset.date[dumb_ix])))
# change those dates
actual_dumb_dates <- convertToDate(onset_dates$onset.date[dumb_ix])
# change the other dates
onset_dates$onset.date <- as.Date(onset_dates$onset.date, format="%m/%d/%Y")
# assign the now NAs with the dumb dates
onset_dates$onset.date[dumb_ix] <- actual_dumb_dates

# add in the MRI and dch dates
MRI_dates <- merge(onset_dates, ichmri[,c("MRN", "mri_date")], by="MRN", all.x=T)  # , "Discharge_date"
# add in the discharge dates from the ichop database instead
# ...Nope, because most of these dates are missing.
# dch_table <- merge(ichopdb$`Admission data`[,c("Subject.ID", "MRN")], ichopdb$Discharge, by="Subject.ID", all.y=TRUE)
# dch_MRI_dates3 <- merge(MRI_dates, dch_table[,c("MRN", "Discharge.Date"), drop=FALSE], by="MRN", all.x=TRUE)
# use the dat table instead
dch_MRI_dates <- unique(merge(MRI_dates, dat[,c("MRN", "ICU.Discharge.Date")], by="MRN", all.x=T))
# this broke now.
dch_MRI_dates$ICU.Discharge.Date <- as.Date(as.numeric(dch_MRI_dates$ICU.Discharge.Date))
# of Cooouuurrseeee we're missing values here.
# find them from the OLD table I guess, wherein they might be WRONG ANYWAY.
missing_mrn <- dch_MRI_dates[which(is.na(dch_MRI_dates$ICU.Discharge.Date)), "MRN"]
# ugh, but these are most likely the hospital discharge dates.
# OH WELL I guess they'll just be missing then.
filter(ichmri, MRN %in% missing_mrn)$Discharge_date

### Have to correct this patient, because they had a re-bleed but weren't re-entered into ICHOP,
### so the onset.date obtained is from the original visit.
dch_MRI_dates[which(dch_MRI_dates$MRN=='REDACTED'),"onset.date"] <- as.Date("2014-07-01")

# this patient has the hospital discharge instead of the ICU discharge.
# replace it with the dates from Adu/Caroline table
## Don't need this if using the dat table for ICU discharge dates

###### IT MIGHT BE WORTH JUST TAKING THE ICU DISCHARGE DATES DIRECTLY FROM THIS TABLE INSTEAD #####

# find the time differences
dch_MRI_dates$ICH_MRI_Delay_days <- as.numeric(difftime(dch_MRI_dates$mri_date, dch_MRI_dates$onset.date, units="days"))
dch_MRI_dates$MRI_Dch_Delay_days <- as.numeric(difftime(dch_MRI_dates$ICU.Discharge.Date, dch_MRI_dates$mri_date, units="days"))

### Some of these patients had an MRI after they were discharged.
### ..unless the discharge date refers to ICU discharge and not hospital discharge
dch_MRI_dates[which(dch_MRI_dates$MRI_Dch_Delay_days < 0), ]

# overall
summary(dch_MRI_dates$ICH_MRI_Delay_days)
summary(dch_MRI_dates$MRI_Dch_Delay_days)

# for conscious at MRI
summary(filter(dch_MRI_dates, MRN %in% c_MRI$MRN)$ICH_MRI_Delay_days)
summary(filter(dch_MRI_dates, MRN %in% c_MRI$MRN)$MRI_Dch_Delay_days)
# for unconscous at MRI
summary(filter(dch_MRI_dates, MRN %in% uc_MRI$MRN)$ICH_MRI_Delay_days)
summary(filter(dch_MRI_dates, MRN %in% uc_MRI$MRN)$MRI_Dch_Delay_days)
# for conscious at dch
summary(filter(dch_MRI_dates, MRN %in% c_Dch$MRN)$ICH_MRI_Delay_days)
summary(filter(dch_MRI_dates, MRN %in% c_Dch$MRN)$MRI_Dch_Delay_days)
# for unconscious at dch
summary(filter(dch_MRI_dates, MRN %in% uc_Dch$MRN)$ICH_MRI_Delay_days)
summary(filter(dch_MRI_dates, MRN %in% uc_Dch$MRN)$MRI_Dch_Delay_days)




#### Normalize volumes and MLS to mean brain volume ####
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



#### MLS ####
#summary(ichmri$`MLS.[mm]`)

# ## At MRI
# ichmri %>%
#   group_by(MRI_Cs2) %>%
#   mutate(avgMLS = mean(MLS_Norm, na.rm=T),
#          sdMLS = sd(MLS_Norm, na.rm=T)) %>%
#   select(MRI_Cs2, avgMLS, sdMLS) %>%
#   distinct()
# 
# # at Dch
# ichmri %>%
#   group_by(follow2) %>%
#   mutate(avgMLS = mean(MLS_Norm, na.rm=T),
#          sdMLS = sd(MLS_Norm, na.rm=T)) %>%
#   select(follow2, avgMLS, sdMLS) %>%
#   distinct()

# median+IQR instead
summary(filter(ichmri, MRN %in% c_MRI$MRN)$MLS_Norm)
summary(filter(ichmri, MRN %in% uc_MRI$MRN)$MLS_Norm)

summary(filter(ichmri, MRN %in% c_Dch$MRN)$MLS_Norm)
summary(filter(ichmri, MRN %in% uc_Dch$MRN)$MLS_Norm)


#### ICH volume ####
## At MRI
# ichmri %>%
#   group_by(MRI_Cs2) %>%
#   mutate(avgHg = mean(Hg_vol_Norm, na.rm=T) / 1000,  # divide by 1000 to convert to mL
#          sdHg = sd(Hg_vol_Norm, na.rm=T) / 1000) %>%
#   select(MRI_Cs2, avgHg, sdHg) %>%
#   distinct()
# 
# # at Dch
# ichmri %>%
#   group_by(follow2) %>%
#   mutate(avgHg = mean(Hg_vol_Norm, na.rm=T),
#          sdHg = sd(Hg_vol_Norm, na.rm=T)) %>%
#   select(follow2, avgMLS, sdMLS) %>%
#   distinct()

summary(filter(ichmri, MRN %in% c_MRI$MRN)$Hg_vol_Norm / 1000)
summary(filter(ichmri, MRN %in% uc_MRI$MRN)$Hg_vol_Norm / 1000)

summary(filter(ichmri, MRN %in% c_Dch$MRN)$Hg_vol_Norm / 1000)
summary(filter(ichmri, MRN %in% uc_Dch$MRN)$Hg_vol_Norm / 1000)



#### Edema volume ####
summary(filter(ichmri, MRN %in% c_MRI$MRN)$Ed_vol_Norm / 1000)
summary(filter(ichmri, MRN %in% uc_MRI$MRN)$Ed_vol_Norm / 1000)

summary(filter(ichmri, MRN %in% c_Dch$MRN)$Ed_vol_Norm / 1000)
summary(filter(ichmri, MRN %in% uc_Dch$MRN)$Ed_vol_Norm / 1000)






#### lobar deep numbers ####
lobardeepNumbers(c_MRI$MRN)
lobardeepNumbers(uc_MRI$MRN)

lobardeepNumbers(c_Dch$MRN)
lobardeepNumbers(uc_Dch$MRN)













##### Logistic Regression for table 2 #####


#### All this table-set-up-stuff needs to be moved to the logistic regression section below

# set up new table with all these parameters to use in the logistic regression
logRegressDat <- data.frame(MRN=ichmri$MRN, MRI_Cs2=ichmri$MRI_Cs2, follow2=ichmri$follow2)

#### GCS ####
# merge the GCS data from ICHOP into this new table
logRegressDat <- merge(logRegressDat, select(ichopdb$`Admission data`, MRN, `GCS.(CUMC)`), by="MRN", all.x=T)
logRegressDat$`GCS.(CUMC)` <- as.numeric(logRegressDat$`GCS.(CUMC)`)

#### Infratent ####
# merge the infra_tent values
logRegressDat <- merge(logRegressDat, select(ichmri, MRN, infra_tent), by="MRN", all.x=T)
# turn infra_tent into 1 and 0 for not None and None, respectively
logRegressDat$infra_tent <- ifelse(logRegressDat$infra_tent == "None", 0, 1)


# tmpscore <- NumberIVH2(ichmri$MRN, MRI=T)
# a <- merge(logRegressDat, tmpscore, by="MRN", all.x=T)

#### Number IVH ####
IVHscores <- FindIVHScores(ichmri[,1, drop=FALSE], MRI=TRUE)
logRegressDat <- merge(logRegressDat, select(IVHscores, MRN, All_Scores), by="MRN", all.x=T)
logRegressDat$All_Scores <- ifelse(logRegressDat$All_Scores != 0, 1, 0)
names(logRegressDat)[6] <- "Has_IVH"



# now merge in the ICH and FUNC scores
#### ICH & FUNC Scores ####
logRegressDat <- merge(logRegressDat, select(ichopdb$`Admission data`, MRN, ICH.Score, FUNC), by="MRN", all.x=T)



## I think there's only a few more things to add to this.
## also need to load in the new ichmri data, because the Dch numbers are different


#### Lobar and Deep ####
#lobardeep_tmp <- merge(select(logRegressDat, MRN, select(ichopdb$Etiology, Location2))

LobarDeep_tab <- GetLobarDeep(logRegressDat$MRN)
logRegressDat <- merge(logRegressDat, LobarDeep_tab, by="MRN", all.x=T)



#### ICH & Edema volume ####
logRegressDat <- logRegressDat %>%
  left_join(select(ichmri, MRN, Hg_vol_Norm, Ed_vol_Norm, MLS_Norm), by="MRN") %>%
  mutate(Hg_vol_Norm = Hg_vol_Norm / 1000,
         Ed_vol_Norm = Ed_vol_Norm / 1000)



#### Days between onset-MRI & MRI-Dch ####
new_CS <- read.xlsx("E:/Experiments/ICH_MRI/Command_scores_and_days.xlsx")
new_CS$MRN <- as.numeric(new_CS$MRN)
new_CS[,c("MRI.Date", "ICU.Admission.Date", "ICU.Discharge.Date")] <- lapply(new_CS[,c("MRI.Date", "ICU.Admission.Date", "ICU.Discharge.Date")],
                                                                             convertToDate)


tmp_dates <- merge(select(new_CS, MRN, MRI.Date, ICU.Discharge.Date), onset_dates, by="MRN", all.x=T)


logRegressDat <- merge(logRegressDat, tmp_dates, by="MRN", all.x=T)




logRegressDat$Time_ICH_to_MRI <- as.numeric(difftime(logRegressDat$MRI.Date, logRegressDat$onset.date, units="days"))
logRegressDat$Time_MRI_to_ICU_dch <- as.numeric(difftime(logRegressDat$ICU.Discharge.Date, logRegressDat$MRI.Date, units="days"))



# ...at some point, everything got duplicated.
logRegressDat <- unique(logRegressDat)





#### regression ####


glmfit <- glm(MRI_Cs2 ~ ICH.Score + FUNC + `GCS.(CUMC)` + Time_ICH_to_MRI + Time_MRI_to_ICU_dch + lobar + deep +
                infra_tent + Hg_vol_Norm + Ed_vol_Norm + MLS_Norm + Has_IVH, data=logRegressDat)


summary(glmfit)
# odds ratios
exp(coef(glmfit))



glmfit2 <- glm(follow2 ~ ICH.Score + FUNC + `GCS.(CUMC)` + Time_ICH_to_MRI + Time_MRI_to_ICU_dch + lobar + deep +
                infra_tent + Hg_vol_Norm + Ed_vol_Norm + MLS_Norm + Has_IVH, data=logRegressDat)

summary(glmfit2)
exp(coef(glmfit2))





