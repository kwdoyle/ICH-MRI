# 
# ## Extract sedation data for the patients that match
# #seddat <- seddatAll[seddatAll$MRN %in% ichmriMRN, ]
# seddat <- seddatM8GS[seddatM8GS$MRN %in% ichmriMRN, ]
# 
# # now remove duplicate rows
# ##### !!!!! Need to remove the Order_ID column (column 11) when checking for duplicates too
# dup_rows <- which(duplicated(seddat[,-11]))
# seddat <- seddat[-dup_rows, ]
# # look for duplicate rows between the different Record_Types
# # apparently there aren't any
# any(duplicated(seddat[,-1]))
# 
# ### wait ok it might just be the drug + admin time that are duplicated 
# # can pull just the drug and admininstered time
# drug_and_time <- seddat[,c("MRN", "Medication_Name", "Administered_Dtm")]
# # see which are duplicated in the entire table
# dup_drug_times <- which(duplicated(drug_and_time))
# #View(seddat[dup_drug_times,])
# # then pull the drug+time and its match from the entire table
# matches <- plyr::match_df(drug_and_time, seddat[dup_drug_times,])
# # use its rownames / row indices
# match_rownms <- rownames(matches)
# # to pull the full rows from the entire table
# actual_matches <- seddat[rownames(seddat) %in% match_rownms, ]
# # most of these have the correct calculated mg amount in the Med_Amt column, and have a 0 for the other row
# # which is just the same drug entered differently into the other flowsheet
# actual_matches <- actual_matches[order(actual_matches$Medication_Name, actual_matches$Administered_Dtm), ]
# View(actual_matches)
# # but there ARE some in here with the Med_Amt missing
# any(is.na(actual_matches$Medication_Amt_mg))
# actual_matches_nas <- actual_matches[is.na(actual_matches$Medication_Amt_mg), ]
# # there's only 16 entries like this ..that's good I guess.
# ### And all of these actual_matches_nas rows are in the overall med_nas below
# ## so I guess to deal with the missing medication_amt_mg amounts, I can just use the med_nas subset.
# 
# ## The duplicates in the above table always have a value of 0 entered for the med_amount.
# ## this SHOULDN'T be a problem, then, because we're just summing the values all together.
# 
# 

# 
# 
# 
# # find where Medication_Amt_mg is missing
# med_nas <- seddat[which(is.na(seddat$Medication_Amt_mg)), ]
# # "default" concentrations for each drug
# # Propofol: 10 mg/mL
# # Fentanyl: 10 ug/mL
# # Dexmedetomidine: ????
# # Midazolam: 1 mg/mL
# 
# # so for those that have the dispences amount in mL
# # and have the concentration missing (or just set as 1 with no units),
# # replace the concentration with the above amount,
# # then calculate the mg by multiplying the volume by the new concentration.
# 
# # all of these have the concentration missing / NA anyway
# med_nas_tofix <- med_nas[med_nas$UM %in% c("ml", "mL", "ML", "Ml"), ]
# 
# 
# 
# #med_nas %>% mutate(Concentration = replace(Concentration, which(Concentration)))
# med_nas <- med_nas %>%  # all of these new concentrations are in mg
#   mutate(Concentration = case_when(grepl("propofol", Medication_Name, ignore.case=T) ~ 10,
#                                    grepl("fentanyl", Medication_Name, ignore.case=T) ~ 0.01,
#                                    grepl("midazolam", Medication_Name, ignore.case=T) ~ 1))
# 
# 
# ## .I think the med_nas rows are already accounted for in seddat, under a different flowsheet with the value actually entered.
# # or rather, I think it's always the ICU/IO Flowsheet that has the data. If it's missing there THEN we need to calculate it
# View(filter(seddat, MRN==med_nas$MRN[3] & Medication_Name==med_nas$Medication_Name[3] & Administered_Dtm==med_nas$Administered_Dtm[3]))
# View(filter(med_nas, MRN==med_nas$MRN[3] & Medication_Name==med_nas$Medication_Name[3] & Administered_Dtm==med_nas$Administered_Dtm[3]))
# 
# 
# View(filter(seddat, MRN==med_nas$MRN[1] & Medication_Name==med_nas$Medication_Name[1] & Administered_Dtm==med_nas$Administered_Dtm[1]))
# 
# 



######### So it SEEMS like the nurse-entering-into-multiple-flowsheets
# and the some-medication-amnts-missing issues are two seperate things