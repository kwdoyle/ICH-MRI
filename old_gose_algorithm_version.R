## old GOSE algorithm

miss_val <- 10

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
    GOSE_3M = case_when((GOSE_3M == 3 | GOSE_3M == miss_val) & mRS_3M == 3 ~ 4,
                        GOSE_3M == miss_val & mRS_3M == 0 ~ 8,
                        GOSE_3M == miss_val & mRS_3M == 1 ~ 7,
                        GOSE_3M == miss_val & mRS_3M == 2 ~ 5,
                        GOSE_3M == miss_val & mRS_3M %in% 4:5 ~ 3,
                        GOSE_3M == miss_val & mRS_3M == 6 ~ 1,
                        TRUE ~ GOSE_3M),
    
    GOSE_6M = case_when((GOSE_6M == 3 | GOSE_6M == miss_val) & mRS_3M == 3 ~ 4,
                        GOSE_6M == miss_val & mRS_3M == 0 ~ 8,
                        GOSE_6M == miss_val & mRS_3M == 1 ~ 7,
                        GOSE_6M == miss_val & mRS_3M == 2 ~ 5,
                        GOSE_6M == miss_val & mRS_3M %in% 4:5 ~ 3,
                        GOSE_6M == miss_val & mRS_3M == 6 ~ 1,
                        TRUE ~ GOSE_6M),
    
    GOSE_12M = case_when((GOSE_12M == 3 | GOSE_12M == miss_val) & mRS_3M == 3 ~ 4,
                         GOSE_12M == miss_val & mRS_3M == 0 ~ 8,
                         GOSE_12M == miss_val & mRS_3M == 1 ~ 7,
                         GOSE_12M == miss_val & mRS_3M == 2 ~ 5,
                         GOSE_12M == miss_val & mRS_3M %in% 4:5 ~ 3,
                         GOSE_12M == miss_val & mRS_3M == 6 ~ 1,
                         TRUE ~ GOSE_12M)
  ) %>%
  mutate(
    GOSE_3M = case_when(GOSE_3M == miss_val & Barthel_3M == 100 ~ 7,
                        TRUE ~ GOSE_3M),
    
    GOSE_6M = case_when(GOSE_6M == miss_val & Barthel_3M == 100 ~ 7,
                        TRUE ~ GOSE_6M),
    
    GOSE_12M = case_when(GOSE_12M == miss_val & Barthel_3M == 100 ~ 7,
                         TRUE ~ GOSE_12M)
  ) %>%
  mutate(
    GOSE_3M = case_when(GOSE_3M == 3 & Lawton_PSMS_3M <= 16 ~ 4, #%in% 8:16 ~ 4,
                        GOSE_3M == 5 & Lawton_PSMS_3M <= 9 ~ 6,  #%in% 6:9 ~ 6,
                        GOSE_3M == miss_val & Lawton_PSMS_3M <= 30 ~ 3,
                        GOSE_3M == miss_val & Lawton_PSMS_3M <= 16 ~ 4,
                        GOSE_3M == miss_val & Lawton_PSMS_3M <= 12 ~ 5,
                        GOSE_3M == miss_val & Lawton_PSMS_3M <= 9 ~ 6,
                        GOSE_3M == miss_val & Lawton_PSMS_3M <= 6 ~ 7,
                        TRUE ~ GOSE_3M),
    
    GOSE_6M = case_when(GOSE_6M == 3 & Lawton_PSMS_6M <= 16 ~ 4, #%in% 8:16 ~ 4,
                        GOSE_6M == 5 & Lawton_PSMS_6M <= 9 ~ 6,  #%in% 6:9 ~ 6,
                        GOSE_6M == miss_val & Lawton_PSMS_6M <= 30 ~ 3,
                        GOSE_6M == miss_val & Lawton_PSMS_6M <= 16 ~ 4,
                        GOSE_6M == miss_val & Lawton_PSMS_6M <= 12 ~ 5,
                        GOSE_6M == miss_val & Lawton_PSMS_6M <= 9 ~ 6,
                        GOSE_6M == miss_val & Lawton_PSMS_6M <= 6 ~ 7,
                        TRUE ~ GOSE_6M),
    
    GOSE_12M = case_when(GOSE_12M == 3 & Lawton_PSMS_12M <= 16 ~ 4, #%in% 8:16 ~ 4,
                         GOSE_12M == 5 & Lawton_PSMS_12M <= 9 ~ 6,  #%in% 6:9 ~ 6,
                         GOSE_12M == miss_val & Lawton_PSMS_12M <= 30 ~ 3,
                         GOSE_12M == miss_val & Lawton_PSMS_12M <= 16 ~ 4,
                         GOSE_12M == miss_val & Lawton_PSMS_12M <= 12 ~ 5,
                         GOSE_12M == miss_val & Lawton_PSMS_12M <= 9 ~ 6,
                         GOSE_12M == miss_val & Lawton_PSMS_12M <= 6 ~ 7,
                         TRUE ~ GOSE_12M)
  ) %>%
  mutate(
    GOSE_3M = case_when((GOSE_3M == 5 | GOSE_3M == miss_val) & Lawton_IADLs_3M <= 15 ~ 6,
                        GOSE_3M == miss_val & Lawton_IADLs_3M <= 20 ~ 5,
                        TRUE ~ GOSE_3M),
    
    GOSE_6M = case_when((GOSE_6M == 5 | GOSE_6M == miss_val) & Lawton_IADLs_3M <= 15 ~ 6,
                        GOSE_6M == miss_val & Lawton_IADLs_3M <= 20 ~ 5,
                        TRUE ~ GOSE_6M),
    
    GOSE_12M = case_when((GOSE_12M == 5 | GOSE_12M == miss_val) & Lawton_IADLs_3M <= 15 ~ 6,
                         GOSE_12M == miss_val & Lawton_IADLs_3M <= 20 ~ 5,
                         TRUE ~ GOSE_12M)
  )