library(openxlsx)
library(daff)
library(irr)



CalcKappas <- function(df1, df2, loc_names) {
  out <- list()
  for (nm in loc_names) {
    print(paste("Kappa for", nm))
    dat <- as.data.frame(cbind(df1[,nm], df2[,nm]))
    
    kp <- kappa2(ratings=dat, weight="unweighted")
    out[[nm]] <- kp
    print(kp)
    cat('\n')
  }
  return(out)
}




setwd("E:/Experiments/ICH_MRI")

#### Load in both datasets and cohorse the rows&columns to match ####
datAlex <- read.xlsx("E:/Experiments/ICH_MRI/Ben files/MRI_ICH/ICH-FinalCompilationFinal_v4_With_TH.xlsx")
datDavid <- read.xlsx("E:/Experiments/ICH_MRI/ICH-David-measure-lesions.xlsx")

# subset original data by the same columns as in the data Dr. Roh reviewed
datAlex2 <- datAlex[,which(colnames(datAlex) %in% colnames(datDavid))]

# ..the new data has the summarized TH in it, which the original data does not
names(datDavid)[which(!names(datDavid) %in% names(datAlex2))]

# I guess remove these columns from the new data for now.
datDavid2 <- datDavid[,which(colnames(datDavid) %in% colnames(datAlex2))]

# sort columns to have same order
#new_order <- which(names(datAlex2)==names(datDavid2))
new_order <- names(datDavid2)
datAlex3 <- datAlex2[,new_order]

## Now subset original data by the same patients as in the new data
datAlex4 <- datAlex3[datAlex3$MRN %in% datDavid2$MRN, ]

# set new variable names
Alex_Data <- datAlex4[order(datAlex4$MRN), ]
David_Data <- datDavid2[order(datDavid2$MRN), ]

# check to make sure the patient order matches between the two tables
all.equal(Alex_Data$MRN, David_Data$MRN)



### Save these tables, then combine the L and R in the same row...?
# write.xlsx(David_Data, "David_Data.xlsx", rowNames=F)
# write.xlsx(Alex_Data, "Alex_Data.xlsx", rowNames=F)


# or do it this way?
##### HOW TO COMBINE THE _C LOCATIONS?????
##### RIGHT NOW I'M JUST DUPLICATING THEM.
##### Ok, just remove them and the IVH for this duplication table. have to compute those seperately later
Alex_L_names <- Alex_Data[,grep("MRN|Name|MRI.date|MRI.time|_L", names(Alex_Data))]
Alex_R_names <- Alex_Data[,grep("MRN|Name|MRI.date|MRI.time|_R", names(Alex_Data))]
# rename columns
names(Alex_L_names)[grep("_L", names(Alex_L_names))] <- gsub("_L", "", names(Alex_L_names)[grep("_L", names(Alex_L_names))])
names(Alex_R_names)[grep("_R", names(Alex_R_names))] <- gsub("_R", "", names(Alex_R_names)[grep("_R", names(Alex_R_names))])

New_Alex_Data <- rbind(Alex_L_names, Alex_R_names)



David_L_names <- David_Data[,grep("MRN|Name|MRI.date|MRI.time|_L", names(David_Data))]
David_R_names <- David_Data[,grep("MRN|Name|MRI.date|MRI.time|_R", names(David_Data))]
# rename columns
names(David_L_names)[grep("_L", names(David_L_names))] <- gsub("_L", "", names(David_L_names)[grep("_L", names(David_L_names))])
names(David_R_names)[grep("_R", names(David_R_names))] <- gsub("_R", "", names(David_R_names)[grep("_R", names(David_R_names))])

New_David_Data <- rbind(David_L_names, David_R_names)

# row order still matches
all.equal(New_Alex_Data$MRN, New_David_Data$MRN)



## Tables for the non-L/R ones
Alex_Extra <- Alex_Data[,c("Vermis_edema", "Vermis_ICH", "MB_edema_C", "MB_ICH_C")]
David_Extra <- David_Data[,c("Vermis_edema", "Vermis_ICH", "MB_edema_C", "MB_ICH_C")]

loc_names2 <- c("Vermis_edema", "Vermis_ICH", "MB_edema_C", "MB_ICH_C")



#"AntPons_ICH_L"
# dat <- as.data.frame(cbind(New_Alex_Data[,"AntPons_ICH"], New_David_Data[,"AntPons_ICH"]))
# kappa2(ratings=dat, weight="unweighted")



#### Find differences ####
#all.equal(datAlex4, datDavid2)

diffs <- diff_data(Alex_Data, David_Data)
render_diff(diffs)






#### Calculate Kappa ####
# kappa2(ratings=cbind(Alex_Data[,"OCx_edema_R"], David_Data[,"OCx_edema_R"]), weight="unweighted")
# kappa2(ratings=cbind(Alex_Data[,"Cereb_edema_R"], David_Data[,"Cereb_edema_R"]), weight="unweighted")
# 


### do it for all locations
all(names(David_Data)[5:length(names(David_Data))] == names(Alex_Data)[5:length(names(Alex_Data))])
loc_names <- names(New_David_Data)[5:length(names(New_David_Data))]


#names(New_David_Data)[5:length(names(New_David_Data))]

# put outputs into a list
# kappas <- list()
# for (nm in loc_names) {
#   print(paste("Kappa for", nm))
#   dat <- as.data.frame(cbind(New_Alex_Data[,nm], New_David_Data[,nm]))
#   #dat2 <- dat
#   # convert to factor?
#   # this doesn't do anything to help when all are the same and kappa==NA.
#   #dat2[,1:2] <- lapply(dat[,1:2], factor, levels=c(0,1))#apply(dat, MARGIN=2, FUN=as.factor)
#   kp <- kappa2(ratings=dat, weight="unweighted")
#   kappas[[nm]] <- kp
#   print(kp)
#   cat('\n')
# }




kappas_LR <- CalcKappas(New_Alex_Data, New_David_Data, loc_names=loc_names)
kappas_other <- CalcKappas(Alex_Extra, David_Extra, loc_names=loc_names2)




# extract values from each kappa measurement
kappa_vals <- c()
for (i in 1:length(kappas)) {
  kappa_vals <- c(kappa_vals, kappas[[i]]$value)
}

hist(kappa_vals)
summary(kappa_vals)



## Write kappa list to a text file
# a little hack-y, I feel, but whatever.
sink("Kappa_Info.txt")
print(kappas_LR)
sink()


sink("Kappa_Info2.txt")
print(kappas_other)
sink()
