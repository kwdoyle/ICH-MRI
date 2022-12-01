## Counts of patients with specific lesions for each conscious status
## (to be run along with "ICH_MRI_stats" script)

# check areas:

# col 1:
# BF_ICH_C
# PUT_ICH_contro, PUT_ICH_ipsi
# GP_ICH_contro, GP_ICH_ipsi
# IC_ant_ICH_contro, IC_ant_ICH_ipsi, IC_post_ICH_contro, IC_post_ICH_ipsi

# col 2:
# BF_ICH_C
# FCx_ICH_contro, FCx_ICH_ipsi
# PUT_ICH_contro, PUT_ICH_ipsi
# GP_ICH_contro, GP_ICH_ipsi
# IC_ant_ICH_contro, IC_ant_ICH_ipsi, IC_post_ICH_contro, IC_post_ICH_ipsi
# Hypo_ICH_C


checkdat <- loc_analyse.raw3 %>%
  select(MRN, MRI_Cs2, follow2, 
         BF_ICH_C, 
         FCx_ICH_contro, FCx_ICH_ipsi,
         PUT_ICH_contro, PUT_ICH_ipsi, 
         GP_ICH_contro, GP_ICH_ipsi,
         IC_ant_ICH_contro, IC_ant_ICH_ipsi, IC_post_ICH_contro, IC_post_ICH_ipsi,
         Hypo_ICH_C)

# I think I can melt this and then just run table() on it.
# Probably not necessary to do this, but I think it's the "cleaner" way, as this method is applicable
# to other potential analyses.
checkdat.m <- melt(checkdat, id.vars = c("MRN", "MRI_Cs2", "follow2"), variable.name = "Brain_Region", value.name = "ICH_Lesion")


countdat <- checkdat.m %>%
  group_by(MRN) %>%
  # make column that is "conscious at dch from those thast were unconscious at MRI"
  mutate(CS_Dch_UC_MRI = case_when(
    (MRI_Cs2 == 0 & follow2 == 1) ~ 1,
    (MRI_Cs2 == 0 & follow2 == 0) ~ 0,
    # if patient was conscious at MRI, then set value to NA
    MRI_Cs2 == 1 ~ NA_real_),
         # make a column that is a 1 if there's a 1 for ANY brain region per patient indicated in column 1 of table
         ICH_check1 = case_when(
    any(
      ICH_Lesion[Brain_Region == "BF_ICH_C" 
                 | Brain_Region == "PUT_ICH_contro"
                 | Brain_Region == "PUT_ICH_ipsi"
                 | Brain_Region == "GP_ICH_contro"
                 | Brain_Region == "GP_ICH_ipsi"
                 | Brain_Region == "IC_ant_ICH_contro"
                 | Brain_Region == "IC_ant_ICH_ipsi"
                 | Brain_Region == "IC_post_ICH_contro"
                 | Brain_Region == "IC_post_ICH_ipsi"] == 1
      ) ~ 1,
    TRUE ~ 0
    ),
    # make a column that is a 1 if there's a 1 for ANY brain region per patient indicated in column 2 of table
    ICH_check2 = case_when(
      any(
        ICH_Lesion[Brain_Region == "BF_ICH_C" 
                   | Brain_Region == "FCx_ICH_contro"
                   | Brain_Region == "FCx_ICH_ipsi"
                   | Brain_Region == "PUT_ICH_contro"
                   | Brain_Region == "PUT_ICH_ipsi"
                   | Brain_Region == "GP_ICH_contro"
                   | Brain_Region == "GP_ICH_ipsi"
                   | Brain_Region == "IC_ant_ICH_contro"
                   | Brain_Region == "IC_ant_ICH_ipsi"
                   | Brain_Region == "IC_post_ICH_contro"
                   | Brain_Region == "IC_post_ICH_ipsi"
                   | Brain_Region == "Hypo_ICH_C"] == 1
      ) ~ 1,
      TRUE ~ 0
    ))




# select the relevant columns and get distinct rows to get one row per patient
countdat2 <- countdat %>%
  select(MRN, MRI_Cs2, follow2, CS_Dch_UC_MRI, ICH_check1, ICH_check2) %>%
  distinct()




## Making main table

# first two rows:
table(MRI_Cs2 = countdat2$MRI_Cs2, 
      ICH_col1 = countdat2$ICH_check1)

table(MRI_Cs2 = countdat2$MRI_Cs2, 
      ICH_col2 = countdat2$ICH_check2)


# last two rows:
# Note: there's a total of 53 patients who were unconscious at MRI
table(CS_Dch_UC_MRI = countdat2$CS_Dch_UC_MRI, 
      ICH_col1 = countdat2$ICH_check1)

table(CS_Dch_UC_MRI = countdat2$CS_Dch_UC_MRI, 
      ICH_col2 = countdat2$ICH_check2)
