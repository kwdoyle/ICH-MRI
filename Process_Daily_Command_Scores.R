library(openxlsx)

dat <- read.xlsx("E:/Experiments/ICH_MRI/MRI Imaging list radiology comprehensive 2.22.17.xlsx",
                           detectDates=FALSE) 
#  detect dates doesn't work--it can actually convert to an incorrect date



dat$mri_date <- convertToDateTime(dat$mri_date)
dat$mri_time <- convertToDateTime(dat$mri_time)
dat$`Command.Score.Note.Date.&.Time` <- convertToDateTime(dat$`Command.Score.Note.Date.&.Time`)
dat$Discharge.note.time <- convertToDateTime(dat$Discharge.note.time)
dat$Discharge_date <- convertToDateTime(dat$Discharge_date)
dat$MRI.eval.note.time <- convertToDateTime(dat$MRI.eval.note.time)
dat$MRI.eval.note.date <- convertToDateTime(dat$MRI.eval.note.date)
#dat$Discharge.note.time <- convertToDateTime(dat$Discharge.note.time)
