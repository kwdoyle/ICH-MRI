# This calculates the percentages for patients who were unconscious at MRI and unconscious at discharge.
# it does NOT have the percentages for their time at MRI and time at discharge, since those would be included in the tables from the above chunk.
# which is FINE, because these categories:
# 1) unconscious at MRI and unconscious at discharge
# 2) unconscious at MRI and conscious at discharge
# 3) conscious at MRI and conscious at discharge
# are meant to be one "column" in those brain region figures anyway.
ucMRIucDC <- dat2.m %>%
  filter( !(variable  %in%  c("Cerebellar tonsillar herniation","MLS [mm]", "Old stroke", "Old ICH",  "Uncal herniation (to which side) ","Transtentorial herniation ",
                              "Uncal herniation (to which side)","Transtentorial herniation","sf_hern",
                              "TH_ant_edema_contro",  "TH_ant_edema_ipsi",  "TH_lat_edema_contro",
                              "TH_lat_edema_ipsi",  "TH_med_edema_contro",  "TH_med_edema_ipsi",
                              "TH_post_edema_contro", "TH_post_edema_ipsi", "TH_ant_ICH_contro",
                              "TH_ant_ICH_ipsi",    "TH_lat_ICH_contro",    "TH_lat_ICH_ipsi",
                              "TH_med_ICH_contro",    "TH_med_ICH_ipsi",    "TH_post_ICH_contro",
                              "TH_post_ICH_ipsi") )) %>%
  # do the ICH/edema filtering step here instead
  filter(grepl("ICH", variable) | variable == "IVH") %>%
  select(MRN, follow2, MRI_Cs2, variable, value) %>%
  distinct() %>%
  # numbers of patients from figure 3 are within this table at this point.
  # now select patients with the specific consciousness levels at each time.
  ### unconscious at MRI & unconscious at discharge
  filter(MRI_Cs2 == 0 & follow2 == 0) %>%
  ## calculate percent of lesions in each area by taking the average of the presence/absence data
  # don't need to group by following @ time anymore, since the time in question now refers to the patients
  # who were unconscious/conscious at mri or discharge, filtered by above
  group_by(variable) %>%
  summarise(percent = mean(value, na.rm=T))




# this calcualtes the percentages for patients who were unconscious at MRI and conscious at discharge
ucMRIcDC <- dat2.m %>%
  filter( !(variable  %in%  c("Cerebellar tonsillar herniation","MLS [mm]", "Old stroke", "Old ICH",  "Uncal herniation (to which side) ","Transtentorial herniation ",
                              "Uncal herniation (to which side)","Transtentorial herniation","sf_hern",
                              "TH_ant_edema_contro",  "TH_ant_edema_ipsi",  "TH_lat_edema_contro",
                              "TH_lat_edema_ipsi",  "TH_med_edema_contro",  "TH_med_edema_ipsi",
                              "TH_post_edema_contro", "TH_post_edema_ipsi", "TH_ant_ICH_contro",
                              "TH_ant_ICH_ipsi",    "TH_lat_ICH_contro",    "TH_lat_ICH_ipsi",
                              "TH_med_ICH_contro",    "TH_med_ICH_ipsi",    "TH_post_ICH_contro",
                              "TH_post_ICH_ipsi") )) %>%
  # do the ICH/edema filtering step here instead
  filter(grepl("ICH", variable) | variable == "IVH") %>%
  select(MRN, follow2, MRI_Cs2, variable, value) %>%
  distinct() %>%
  # numbers of patients from figure 3 are within this table at this point.
  # now select patients with the specific consciousness levels at each time.
  ### unconscious at MRI & unconscious at discharge
  filter(MRI_Cs2 == 0 & follow2 == 1) %>%
  ## calculate percent of lesions in each area by taking the average of the presence/absence data
  # don't need to group by following @ time anymore, since the time in question now refers to the patients
  # who were unconscious/conscious at mri or discharge, filtered by above
  group_by(variable) %>%
  summarise(percent = mean(value, na.rm=T))




# this calculates the percentages for patients who were conscious at MRI and conscious at discharge
cMRIcDC <- dat2.m %>%
  filter( !(variable  %in%  c("Cerebellar tonsillar herniation","MLS [mm]", "Old stroke", "Old ICH",  "Uncal herniation (to which side) ","Transtentorial herniation ",
                              "Uncal herniation (to which side)","Transtentorial herniation","sf_hern",
                              "TH_ant_edema_contro",  "TH_ant_edema_ipsi",  "TH_lat_edema_contro",
                              "TH_lat_edema_ipsi",  "TH_med_edema_contro",  "TH_med_edema_ipsi",
                              "TH_post_edema_contro", "TH_post_edema_ipsi", "TH_ant_ICH_contro",
                              "TH_ant_ICH_ipsi",    "TH_lat_ICH_contro",    "TH_lat_ICH_ipsi",
                              "TH_med_ICH_contro",    "TH_med_ICH_ipsi",    "TH_post_ICH_contro",
                              "TH_post_ICH_ipsi") )) %>%
  # do the ICH/edema filtering step here instead
  filter(grepl("ICH", variable) | variable == "IVH") %>%
  select(MRN, follow2, MRI_Cs2, variable, value) %>%
  distinct() %>%
  # numbers of patients from figure 3 are within this table at this point.
  # now select patients with the specific consciousness levels at each time.
  ### unconscious at MRI & unconscious at discharge
  filter(MRI_Cs2 == 1 & follow2 == 1) %>%
  ## calculate percent of lesions in each area by taking the average of the presence/absence data
  # don't need to group by following @ time anymore, since the time in question now refers to the patients
  # who were unconscious/conscious at mri or discharge, filtered by above
  group_by(variable) %>%
  summarise(percent = mean(value, na.rm=T))



## ...combine all the tables so I can plot them all together?
consc_loc_combined <- ucMRIucDC %>%
  left_join(ucMRIcDC, by = "variable")

# ..this should be the same as if I just created this table in one shot.
# wow, it is. Jeez, it could be this simple: