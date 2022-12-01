library(openxlsx)



path <- "E:/Experiments/AuditoryExplorationConscious/Sedation data/as_xlsx"

datfiles <- list.files(path=path, pattern="*.xlsx")


ichmri <- read.xlsx("E:/Experiments/ICH_MRI/MRI_Merged_Data_158patients.xlsx")

seddatAll <- data.frame()

for (i in 1:length(datfiles)) {
  print(paste("loading file", i))
  wb <- loadWorkbook(paste(path, "/", datfiles[i], sep=""))
  nms <- names(wb)
  
  for (j in 1:length(nms)) {
    print(paste("loading sheet", j))
    if (j == 1) {
      sht <- readWorkbook(wb, sheet=nms[j], colNames = TRUE)
    } else {
      sht <- readWorkbook(wb, sheet=nms[j], colNames = FALSE)
      names(sht) <- names(seddatAll)
    }
    
    seddatAll <- rbind(seddatAll, sht)
    
  }
}


seddatAll$MRN <- gsub(" ", "", seddatAll$MRN)
seddatAll$MRN <- as.numeric(seddatAll$MRN)





## Here goes nothing...
ichmriMRN <- unique(ichmri$MRN)
sedMRN <- unique(seddatAll$MRN)

length(which(ichmriMRN %in% sedMRN))



## Extract sedation data for the patients that match
seddat <- seddatAll[seddatAll$MRN %in% ichmriMRN, ]







