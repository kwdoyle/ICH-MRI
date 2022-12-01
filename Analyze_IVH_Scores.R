library(dplyr)

## Function to read all sheets of an excel file and puts them all into a list
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


ichop <- readDB("E:/SHOP & ICHOP forms/ICHOP_1_17_18_database_full.xlsx")
maindat <- read.xlsx("E:/Experiments/ICH_MRI/MRI_Merged_Data_158patients.xlsx")

vent_names <- c("IVH.Score.R.lateral", "IVH.Score.L.lateral", "IVH.Score.Third", "IVH.Score.Fourth")
table_names <- c("CT Scan (OSH)", "CT Scan (Admission)", "CT Scan (2nd in 24 hr)", "CT Scan (3rd in 24hr)", "CT Scan (3-6 day)")

# pull these column names along with MRN from all the CT tables and put them into one data frame
# with a column indicating which CT scan they're from

all_CT_IVHScore <- data.frame()
for (nm in table_names) {
  
  tmp <- ichop[[nm]] %>%
    select_(.dots=c("Subject.ID", vent_names)) %>%
    left_join(select(ichop$`Admission data`, Subject.ID, MRN),
              by="Subject.ID") %>%
    mutate(CTscan = nm)
  
  tmp2 <- reshape2::melt(tmp, id.vars=c("Subject.ID", "MRN", "CTscan"), 
                         variable.name="ventrical", value.name="IVH_score")
  
  all_CT_IVHScore <- rbind(all_CT_IVHScore, tmp2)
  
}


all_CT_IVHScore.f <- all_CT_IVHScore %>%
  filter(!is.na(IVH_score))





IVH_w_Dch <- 








