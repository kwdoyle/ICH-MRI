library(openxlsx)
library(dplyr)
library(caret)
library(glmnet)
library(reshape2)
library(pROC)
library(ROCR)
library(ggplot2)
library(lubridate)
library(cvAUC)



meanNormalization <- function(x) {
  res <- (x - mean(x, na.rm=T)) / (max(x, na.rm=T) - min(x, na.rm=T))
  return(res)
}

calcConfInt <- function(data) {
  m    <- mean(data)
  stdv <- sd(data)
  n    <- length(data)
  
  err  <- qt(0.975, df=n-1) * stdv / sqrt(n)
  CI.l <- m - err
  CI.u <- m + err
  
  return(c(lower=CI.l, upper=CI.u))
}

MakePlot <- function(data, group, param, path=NULL, save=FALSE) {
  
  dat <- select_(data, group, param)
  dat2 <- melt(dat, id=group)
  dat2[,group] <- as.factor(dat2[,group])
  
  plt <- ggplot(data=dat2, aes_string(x=group, y="value", fill="variable")) +
    geom_boxplot(width=0.2, outlier.shape=NA) + ggtitle(param) + theme_minimal() +
    geom_jitter(position = position_jitter(width = 0.05, height = 0.1), size=1)
  
  if (save==TRUE) {
    ggsave(filename=paste(path, "/", param, "_", group, ".pdf", sep=""),
           plot=plt, device="pdf")
  }
  
  return(plt)
  
}

RunModel <- function(data, indep_var_ix, dep_var_vals, opposite_vals=NA, n_runs, lambda_grid=NULL, p=0.80, alpha=0.95, nfolds=10, lambda_use="lambda.min") {
  
  AUCs <- c()
  weights_table <- data.frame()
  probs_table <- as.data.frame(matrix(nrow=nrow(data)))
  rownames(probs_table) <- rownames(data)
  lambda_mins <- c()
  for (i in 1:n_runs) {
    # for some reason, this defaulted to rownames.force=T before???? It NEEDS to be true, 
    # or else the row indices for each patient won't carry over in order to make the confusion matrix
    df <- data.matrix(data[,indep_var_ix], rownames.force=T)
    # add dependent variable to beginning of table
    df <- cbind(as.factor(dep_var_vals), df)
    
    # create testing/training sets
    train_ix <- createDataPartition(as.factor(df[,1]), p=p, list=F)
    train <- df[train_ix,]
    test <- df[-train_ix,]
    
    # run model
    modfit <- cv.glmnet(x=train[,2:ncol(train)], y=train[,1], lambda=lambda_grid,
                        family="binomial", type.measure="auc", alpha=alpha, nfolds=nfolds)
    
    # then use the lambda min in the regular glmnet function?
    # modfit <- glmnet(x=train[,2:ncol(train)], y=train[,1], lambda=modfit$lambda.min,
    #                  family="binomial", alpha=alpha)
    lambda_mins <- c(lambda_mins, modfit$lambda.min)
    
    # get coefficients/weights
    weights <- t(coef(modfit, s=lambda_use)[,1])
    weights_table <- rbind(weights_table, weights)
    
    # predict using the test data
    pred.prob <- predict(modfit, newx=data.matrix(test[,2:ncol(test)]),
                         type="response", s=lambda_use)
    
    # save the probabilities for the patients in the current test dataset
    # using the current row indices of the test table and matching to the
    # row indices in the main table
    probs_table <- cbind(probs_table, X=pred.prob[match(rownames(probs_table), rownames(pred.prob))])
    
    # predict new values
    pred <- prediction(pred.prob, test[,1])
    # get AUC
    auc <- performance(pred, "auc")
    # save aucs to vector
    # 'auc' is an S4 class object; need to extract using '@'
    AUCs <- c(AUCs, auc@y.values[[1]])
  }
  
  # set first column as row indices for probs_table and reset column names
  colnames(probs_table) <- as.character(1:length(colnames(probs_table)))
  probs_table[,1] <- as.numeric(rownames(probs_table))
  colnames(probs_table)[1] <- "patient"
  
  # now melt table and calculate confusion matrix, assigning pred.values < 0.5
  # as unconscious
  medProbs <- melt(probs_table, id="patient") %>%
    group_by(patient) %>%
    summarize(medProb = median(value, na.rm=T))
  
  medProbs$class_vals <- ifelse(medProbs$medProb < 0.5, yes=0, no=1)
  # order by patient row index, that way the dep_var_vals match up with the same patient
  # I think this step is unnecessary now that the row-index column is numeric, but w/e
  medProbs <- medProbs[order(medProbs$patient), ]
  medProbs$actual_vals <- dep_var_vals
  medProbs$opposite_vals <- opposite_vals
  cm <- table(actual=medProbs$actual_vals, predict=medProbs$class_vals)
  if (!is.na(opposite_vals)) { # should probably use an any() function on this check to stop the warning message
    cm2 <- table(opposite_meas=medProbs$opposite_vals, predict=medProbs$class_vals)
  }
  
  
  # return weights table and AUCs together in a list
  if (any(!is.na(opposite_vals))) { # should probably use an any() function on this check to stop the warning message
    out <- list(AUC=AUCs, lambda_min=lambda_mins, weights=weights_table, probs=probs_table, medProb=medProbs, confus_mat=cm, oppo_confus_mat=cm2)
  } else {
    out <- list(AUC=AUCs, lambda_min=lambda_mins, weights=weights_table, probs=probs_table, medProb=medProbs, confus_mat=cm)
  }
  
  return(out)
  
  
}



RunModelNew <- function(data, indep_var_ix, dep_var_vals, n_runs, lambda_grid=NULL, p=0.80, alpha=0.95, nfolds=10, lambda_use="lambda.min") {
  df <- data.matrix(data[,indep_var_ix], rownames.force=T)
  # add dependent variable to beginning of table
  df <- cbind(as.factor(dep_var_vals), df)
  
  # create testing/training sets
  train_ix <- createDataPartition(as.factor(df[,1]), p=p, list=F)
  train <- df[train_ix,]
  test <- df[-train_ix,]
  
  # run model
  modfit <- cv.glmnet(x=train[,2:ncol(train)], y=train[,1], lambda=lambda_grid,
                      family="binomial", type.measure="auc", alpha=alpha, nfolds=nfolds)
  
  # then use the lambda min in the regular glmnet function?
  # modfit <- glmnet(x=train[,2:ncol(train)], y=train[,1], lambda=modfit$lambda.min,
  #                  family="binomial", alpha=alpha)
  lambda_min <- modfit$lambda.min
  
  # get coefficients/weights
  weights <- t(coef(modfit, s=lambda_use)[,1])
  #weights_table <- rbind(weights_table, weights)
  
  # predict using the test data
  preds <- predict(modfit, newx=data.matrix(test[,2:ncol(test)]),
                       type="class", s=lambda_use)
  
  # save the probabilities for the patients in the current test dataset
  # using the current row indices of the test table and matching to the
  # row indices in the main table
  #probs_table <- cbind(probs_table, X=pred.prob[match(rownames(probs_table), rownames(pred.prob))])
  
  # predict new values
  #pred <- prediction(pred.prob, test[,1])
  # get AUC
  AUCs <- cvAUC(as.data.frame(as.numeric(preds)), as.data.frame(as.numeric(test[,1])), folds=100)
  #auc <- performance(pred, "auc")
  # save aucs to vector
  # 'auc' is an S4 class object; need to extract using '@'
  AUCs <- c(AUCs, auc@y.values[[1]])
}



ReformatWeightTable <- function(weights, toReturn="unformatted") {
  # this function takes the weight table output from RunModel
  # and reformats it into a better table.
  
  weight_table <- melt(weights) %>%
    filter(variable != "(Intercept)") %>%
    group_by(variable) %>%
    summarise(median_weight = median(value),
              qt_25 = quantile(value)[2],
              qt_75 = quantile(value)[4]) %>%
    arrange(median_weight, qt_25, qt_75)
  
  if (toReturn == "unformatted") {
    return(weight_table)
    
    
  } else if (toReturn == "formatted") {
    
    ROIs <- unique(weight_table$variable)
    ROIs <- unique(gsub("_edema", "", gsub("_ICH", "", ROIs)))
    newtab <- matrix(nrow=length(ROIs), ncol=2)
    rownames(newtab) <- ROIs
    colnames(newtab) <- c("ICH", "Edema")
    # try to assign values
    for (i in 1:nrow(weight_table)) {
      cur_name <- as.character(weight_table$variable[i])
      if (!cur_name %in% c("Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm", "IVH")) {
        # make name match that from newtab
        name_wo_type <- gsub("_edema", "", gsub("_ICH", "", cur_name))
        # find the row to add the value to
        newtab_ix <- which(rownames(newtab) == name_wo_type)
        # now find out if this was ich or edema
        if (grepl("edema", cur_name)) {
          newtab_col <- "Edema"
        } else if (grepl("ICH", cur_name)) {
          newtab_col <- "ICH"
        }
        
        # now paste the median weight and CI together
        val_to_add <- paste(round(weight_table$median_weight[i], 2), "[", 
                            round(weight_table$qt_25[i], 2), " ", round(weight_table$qt_75[i], 2), "]", sep="")
        
        newtab[newtab_ix, newtab_col] <- val_to_add
        
        #nm_check <- unlist(strsplit(cur_name, "_"))
        
      } else {
        # just put the volume ones in the ICH column.
        newtab_ix <- which(rownames(newtab) == cur_name)
        
        val_to_add <- paste(round(weight_table$median_weight[i], 2), "[", 
                            round(weight_table$qt_25[i], 2), " ", round(weight_table$qt_75[i], 2), "]", sep="")
        
        newtab[newtab_ix, "ICH"] <- val_to_add
        
      }
    }
    
    # something here to reorder the rows. need to exclude volumes if table has no volumes
    if (all(c("Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm") %in% weight_table$variable)) {
      name_order <- c("Caudate_contro", "Caudate_ipsi", "PUT_contro", "PUT_ipsi", "GP_contro", "GP_ipsi",
                      "TH_contro", "TH_ipsi", "BF_C", "Hypo_C", "MB_peduncle_contro", "MB_peduncle_ipsi",
                      "MB_C", "Teg_C", "Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm", "IVH")
    } else {
      name_order <- c("Caudate_contro", "Caudate_ipsi", "PUT_contro", "PUT_ipsi", "GP_contro", "GP_ipsi",
                      "TH_contro", "TH_ipsi", "BF_C", "Hypo_C", "MB_peduncle_contro", "MB_peduncle_ipsi",
                      "MB_C", "Teg_C")
    }
    
    newtab <- newtab[name_order, ]
    
    
    return(newtab)
    
  } else {
    stop("toReturn must be either 'unformatted' or 'formatted'")
  }
  
}





# same function as above, but performs Leave-One-Out cross validation (not used here, actually.. this was for testing earlier;
# this model gives very similar results to the train/test on 80% model as performed in the above function)
RunLOOModel <- function(data, indep_var_ix, dep_var_vals, lambda_grid=NULL, p=0.80, alpha=0.95, nfolds=10, lambda_use="lambda.min") {
  
  AUCs <- c()
  weights_table <- data.frame()
  probs_table <- as.data.frame(matrix(nrow=nrow(data)))
  rownames(probs_table) <- rownames(data)
  for (i in 1:nrow(data)) {
    
    df <- data.matrix(data[,indep_var_ix], rownames.force=T)
    # add dependent variable to beginning of table
    df <- cbind(as.factor(dep_var_vals), df)
    
    # create testing/training sets
    train <- df[-i,]
    # need to transpose the matrix because, for SOME REASON,
    # matrices default to transposing the dimensions when only 1
    # row is indexed out. So I have to RE transpose it back.
    test <- t(df[i,])
    
    # run model
    modfit <- cv.glmnet(x=train[,2:ncol(train)], y=train[,1], lambda=lambda_grid,
                        family="binomial", type.measure="auc", alpha=alpha, nfolds=nfolds)
    
    
    # get coefficients/weights
    weights <- t(coef(modfit, s=lambda_use)[,1])
    weights_table <- rbind(weights_table, weights)
    
    # predict using the test data
    pred.prob <- predict(modfit, newx=t(data.matrix(test[,2:ncol(test)])),
                         type="response", s=lambda_use)
    
    # for this function, can just append probabilities to a vector,
    # since the order of probabilities always refers to the row order of the data table
    probs_table[i,1] <- pred.prob
    
    
  }
  
  # set first column as row indices for probs_table and reset column names
  #colnames(probs_table) <- as.character(1:length(colnames(probs_table)))
  #probs_table[,1] <- as.numeric(rownames(probs_table))
  colnames(probs_table)[1] <- "probability"
  
  
  probs_table$class_vals <- ifelse(probs_table$probability < 0.5, yes=0, no=1)
  # order by patient row index, that was the dep_var_vals match up with the same patient
  # I think this step is unnecessary now that the row-index column is numeric, but w/e
  #medProbs <- medProbs[order(medProbs$patient), ]
  probs_table$actual_vals <- dep_var_vals
  cm <- table(actual=probs_table$actual_vals, predict=probs_table$class_vals)
  cm <- confusionMatrix(cm)
  
  # Calculate ROC AUC
  rc <- roc(actual_vals ~ probability, data=probs_table, auc=T, ci=T, na.rm=T)
  #ci <- ci.thresholds(rc)
  rc <- as.numeric(rc$auc)
  
  
  # return weights table and AUCs together in a list
  #out <- list(AUC=AUCs, weights=weights_table, probs=probs_table, medProb=medProbs, confus_mat=cm)
  out <- list(AUC=rc, probs=probs_table, confus_mat=cm)
  return(out)
  
  
}

bestAlphaLambda <- function(data, indep_var_ix, dep_var_vals, lambda.grid, alpha.grid, method="repeatedcv", p=0.8, number=10, repeats=5) {
  library(caret)
  
  trainCtrl <- trainControl(method=method,
                            number=number,
                            repeats=repeats,
                            p=p)
  srchGrd <- expand.grid(.alpha=alpha.grid, .lambda=lambda.grid)
  
  my.train <- train(x=data.matrix(data[,indep_var_ix]), y=as.factor(dep_var_vals),
                    method="glmnet",
                    tuneGrid=srchGrd,
                    trControl=trainCtrl,
                    standardize=TRUE, maxit=1000000)
  
  return(my.train)
}




makeGOSEplot <- function(medProbs, var1, var2) {
  
  fp <- which(medProbs[,var1] == 1 & medProbs[,var2] == 0)
  tn <- which(medProbs[,var1] == 0 & medProbs[,var2] == 0)
  tp <- which(medProbs[,var1] == 1 & medProbs[,var2] == 1)
  fn <- which(medProbs[,var1] == 0 & medProbs[,var2] == 1)
  
  fpTab <- loc_analyse.raw3[fp, c("MRN", "GOSE_All"), drop=F]
  fpTab$class_grp <- "false_pos"
  tnTab <- loc_analyse.raw3[tn, c("MRN", "GOSE_All"), drop=F]
  tnTab$class_grp <- "true_neg"
  tpTab <- loc_analyse.raw3[tp, c("MRN", "GOSE_All"), drop=F]
  tpTab$class_grp <- "true_pos"
  fnTab <- loc_analyse.raw3[fn, c("MRN", "GOSE_All"), drop=F]
  fnTab$class_grp <- "false_neg"
  
  GOSE_class <- rbind(fpTab, tnTab, tpTab, fnTab)
  GOSE_class$GOSE_All <- factor(GOSE_class$GOSE_All, exclude=NULL)
  
  return(GOSE_class)
  
}


# this uses the above function to generate the GOSE contingency table, the fisher test on that table,
# and the plot of how many GOSE scores are in each category of the input model's "opposite" confusion matrix.
GOSEfisherAndPlot <- function(data, pred_val="class_vals", actual_val="opposite_vals") {
  # first use the function that makes the data and then
  # do the fisher test by creating the count table using the GOSE data
  pltDat <- makeGOSEplot(data[["medProb"]], pred_val, actual_val)
  
  summaryTable <- pltDat %>%
    filter(class_grp %in% c("false_pos", "true_neg")) %>%
    mutate(gose_grp=as.factor(ifelse(as.numeric(as.character(GOSE_All)) >= 4, 1, 0)),
           class_grp=as.factor(class_grp)) %>%
    filter(!is.na(gose_grp)) %>%
    group_by(class_grp, gose_grp)
  
  contingTable <- table(summaryTable$gose_grp, summaryTable$class_grp)
  fish_out <- fisher.test(contingTable)
  
  # now make the plot
  plt <- ggplot(pltDat, aes(x=class_grp, fill=GOSE_All)) + 
    geom_bar(stat="count", position="dodge") +
    scale_y_continuous(breaks = scales::pretty_breaks(n=10)) +
    xlab("Classification Group") +
    ggtitle(paste("GOSE Distribution per ", deparse(substitute(data)), " 'opposite' Confusion Matrix Classification", sep=""), 
            subtitle=paste("true positive = conscious; true negative = unconscious; model = ", deparse(substitute(data)), sep="")) +
    theme_bw() +
    coord_flip()
  
  all_dat <- list(contingTable=contingTable,
                  fisher=fish_out,
                  the_plot=plt)
  
  return(all_dat)
  
}




bootstrapCI <- function(data) {
  
  
  
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


# find accuracy, specificity, etc. for a confusion matrix
CMstats <- function(confus_mat) {
  tot <- sum(confus_mat)
  # percentages/rates of type I and II error
  tp <- confus_mat[2,2]
  tn <- confus_mat[1,1]
  fp <- confus_mat[1,2]
  fn <- confus_mat[2,1]
  act_yes <- sum(confus_mat[1,])
  act_no <- sum(confus_mat[2,])
  
  acc <- (tp + tn) / tot # Accuracy
  miss <- (fp + fn) / tot # Misclassification rate
  typeIerr <- fp / act_yes # Type I error rate (predict yes when actually no)
  spec <- tn / act_yes # Specificity (how often correctly predict no), ie, 1-Type I error
  typeIIerr <- fn / act_no # Type II error rate (predict no when actually yes)
  sens <- tp / act_no # Sensitivity (how often correctly predict yes), ie, 1-Type II error
  null_err <- act_no / tot # null error rate (percent wrong if always predicted the majority class)
  cohen <- acc - null_err # Cohen's Kappa (how well classifier performs compared to chance)
  
  #out <- list()
  outv <- c("Accuracy"=acc, "Misclass. rate"=miss, "Type I error"=typeIerr, "Specificity"=spec,
            "Type II error"=typeIIerr, "Sensitivity"=sens, "Null error"=null_err, "Cohen's kappa"=cohen)
  
  return(outv)
}









#### Load and clean data ####

loc_analyse.raw3 <- read.xlsx("E:/Experiments/ICH_MRI/MRI_Merged_Data_158patients.xlsx")
# rename the new location parameters to have a "_C" for "central" in contrast with the "ipsi/contro"
# just kidding, we renamed these parameters within the original excel file instead.
#names(loc_analyse.raw3)[which(names(loc_analyse.raw3) %in% c("Hypo_edema", "BF_edema", "BF_ICH"))] <- c("Hypo_edema_C", "BF_edema_C", "BF_ICH_C")
# MRI eval note date table
MRI_time_dat <- read.xlsx("E:/Experiments/ICH_MRI/Ben files/MRI_ICH/MRI Imaging list radiology comprehensive-br2.xlsx")


# convert dates
loc_analyse.raw3$Discharge_date<-as.Date(loc_analyse.raw3$Discharge_date,origin = "1899-12-30")
loc_analyse.raw3$DEATH_DATE<-as.Date(loc_analyse.raw3$DEATH_DATE,origin = "1899-12-30")
loc_analyse.raw3$MRI_date_Loc<-as.Date(loc_analyse.raw3$MRI_date_Loc,origin = "1899-12-30")
# obnoxious way to get datetimes from excel number "datetimes"
MRI_time_dat$MRI.eval.note.date<-as.POSIXct(as.numeric(MRI_time_dat$MRI.eval.note.date) * (60*60*24) ,origin = "1899-12-30", tz="GMT")

loc_analyse.raw3 <- merge(loc_analyse.raw3, MRI_time_dat[,c("MRN", "MRI.eval.note.date")], by="MRN")

names(loc_analyse.raw3) <- gsub("\\.", "_", names(loc_analyse.raw3))
names(loc_analyse.raw3) <- gsub("\\(|\\)", "", names(loc_analyse.raw3))

# change variable names
#names(loc_analyse.raw3)[which(names(loc_analyse.raw3)=='Uncal_herniation_(to_which_side)')] <- "Uncal_herniation_to_which_side"
names(loc_analyse.raw3)[which(names(loc_analyse.raw3)=='MLS_[mm]')] <- "MLS"  # change this too and then also change the name in the other logisitc regress

# remove column "Hypo_ICH_C", because apparently we're not using it and alsoe
# there is a patient with an NA value here.
# NOPE NOW WE CAN USE IT
#loc_analyse.raw3 <- loc_analyse.raw3[,-which(names(loc_analyse.raw3) == "Hypo_ICH_C")]
# this one too.
loc_analyse.raw3 <- loc_analyse.raw3[,-which(names(loc_analyse.raw3) == "Uncal_herniation_to_which_side")]

# remove the 'summary' Meso_ICH_2 columns
loc_analyse.raw3 <- loc_analyse.raw3[,-which(names(loc_analyse.raw3) == "Meso_ICH_2_ipsi" | names(loc_analyse.raw3) == "Meso_ICH_2_contro")]



## change the additional 6 other patients with missing Edema volume to 0
loc_analyse.raw3[,"Ed_vol"][which(is.na(loc_analyse.raw3[,"Ed_vol"]))] <- 0





### Normalized hemorrhage volume
# removes outlier w/ small brain b/c of incomplete scan (according to Kay)
# Will replace this value with the mean value
min_val <- which(loc_analyse.raw3$Brain_vol == min(loc_analyse.raw3$Brain_vol))
mBV <- mean(loc_analyse.raw3$Brain_vol[-min_val])
# replace min value
loc_analyse.raw3$Brain_vol[min_val] <- mBV


loc_analyse.raw3$Hg_vol_Norm <- (loc_analyse.raw3$Hg_vol * loc_analyse.raw3$Brain_vol) / mBV


# normalize edema volume
## ...these didn't change that much
loc_analyse.raw3$Ed_vol_Norm <- (loc_analyse.raw3$Ed_vol * loc_analyse.raw3$Brain_vol) / mBV


# normalize MLS
loc_analyse.raw3$MLS_Norm <- (loc_analyse.raw3$MLS * loc_analyse.raw3$Brain_vol) / mBV


## Now get the z-score for these normalized brain volumes to normalize them again
## basically so that the values aren't so big
# mHgvN <- mean(loc_analyse.raw3$Hg_vol_Norm)
# sdHgv <- sd(loc_analyse.raw3$Hg_vol_Norm)
# loc_analyse.raw3$Hg_vol_zscore <- (loc_analyse.raw3$Hg_vol_Norm - mHgvN) / sdHgv

#### Normalize all by mean instead
loc_analyse.raw3[,c("Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm")] <- lapply(loc_analyse.raw3[,c("Hg_vol_Norm", "Ed_vol_Norm", "MLS_Norm")], 
                                                                                     meanNormalization)
# 
# # z-score for median shift too
# mMLS <- mean(loc_analyse.raw3$MLS_Norm)
# sdMLS <- sd(loc_analyse.raw3$MLS_Norm)
# loc_analyse.raw3$MLS_zscore <- (loc_analyse.raw3$MLS_Norm - mMLS) / sdMLS
# 
# 
# # z-score for edema volume
# mEdvN <- mean(loc_analyse.raw3$Ed_vol_Norm)
# sdEdv <- sd(loc_analyse.raw3$Ed_vol_Norm)
# loc_analyse.raw3$Ed_vol_zscore <- (loc_analyse.raw3$Ed_vol_Norm - mEdvN) / sdEdv



#### Check median length of stay for patients ####
new_CS <- read.xlsx("E:/Experiments/ICH_MRI/Command_scores_and_days.xlsx")
new_CS$MRN <- as.numeric(new_CS$MRN)

# get only for the 158 with good MRI and don't include columns for command score day and command score
new_CS <- new_CS %>%
  select(MRN:ICU.Discharge.Date, Day_of_ICU_Dch:CS_Hosp_Dch) %>%
  filter(MRN %in% loc_analyse.raw3$MRN) %>%
  distinct()


summary(new_CS$Day_of_ICU_Dch)
summary(new_CS$Day_of_Hosp_Dch)





#### Load sedation data ####

#path <- "E:/Experiments/AuditoryExplorationConscious/Sedation data/as_xlsx"
# path <- "E:/Experiments/ICH_MRI/Sedation/as_xlsx"
# datfiles <- list.files(path=path, pattern="*.xlsx")
# 
# seddatAll <- data.frame()
# for (i in 1:length(datfiles)) {
#   print(paste("loading file", i))
#   wb <- loadWorkbook(paste(path, "/", datfiles[i], sep=""))
#   nms <- names(wb)
#   
#   for (j in 1:length(nms)) {
#     print(paste("loading sheet", j))
#     if (j == 1) {
#       sht <- readWorkbook(wb, sheet=nms[j], colNames = TRUE, startRow=4)
#     } else {
#       sht <- readWorkbook(wb, sheet=nms[j], colNames = FALSE)
#       names(sht) <- names(seddatAll)
#     }
#     
#     seddatAll <- rbind(seddatAll, sht)
#     
#   }
# }
# 
# # remove spaces in MRNs
# seddatAll$MRN <- gsub(" ", "", seddatAll$MRN)
# seddatAll$MRN <- as.numeric(seddatAll$MRN)
# # convert dates
# seddatAll$Administered.Dtm <- as.POSIXct(as.numeric(seddatAll$Administered.Dtm) * (60*60*24) ,origin = "1899-12-30", tz="GMT")
# # other cleanings
# names(seddatAll) <- gsub("\\.", "_", names(seddatAll))
# names(seddatAll)[17] <- "UM2"
# names(seddatAll)[ncol(seddatAll)] <- "Medication_Amt_mg"
# 
# 
# 
# ### Only take patients from "MILSTEIN 8 GARDEN SOUTH", save this table, then delete the original.
# seddatM8GS <- seddatAll %>%
#   filter(LocUnit == "MILSTEIN 8 GARDEN SOUTH")
# 
# # save as RDS file
# saveRDS(seddatM8GS, file="E:/Experiments/ICH_MRI/Sedation/sedation_M8GS.rds")
# excel actualy does work--might as well stick with it then
# although the excel sheet reads slower AND the dates convert back to stupid numbers..
#write.xlsx(seddatM8GS, "E:/Experiments/ICH_MRI/Sedation/sedation_M8GS.xlsx", rowNames=F)

### Only need the RDS file now.
seddatM8GS <- readRDS("E:/Experiments/ICH_MRI/Sedation/sedation_M8GS.rds")




# unique patients
# ichmriMRN <- unique(loc_analyse.raw3$MRN)
# 
# # empty table to bind useful rows to
# newSeddat <- data.frame()
# 
# for (i in 1:length(ichmriMRN)) {
#   # find all unique timepoints for patient in this loop
#   tmpdat <- filter(seddatM8GS, MRN == ichmriMRN[i])
#   timepts <- unique(tmpdat$Administered_Dtm)
#   ICU_fs <- filter(tmpdat, Record_Type == "20 ICU I/O Flowsheet")
#   for (j in 1:nrow(ICU_fs)) {
#     # look at the ICU_fs data and check if any Medication_Amts are missing
#   }
# }
# 
# 
# 
# 
# 
















## Here goes nothing...
ichmriMRN <- unique(loc_analyse.raw3$MRN)
sedMRN <- unique(seddatM8GS$MRN)

# now there are 107 patients out of 158 in the sedation database
# and these are the same between the full dataset and the M8GS subsetted dataset
length(which(ichmriMRN %in% sedMRN))
length(which(ichmriMRN %in% seddatM8GS$MRN))




#### SEDATION CLEANING METHOD ####

# maybe can:
# 0) First filter out rows where Dispensed_Amt AND Medication_Amt_mg are 0 or NA
# 1) remove duplicates based on MRN, Med_Name, Admin_Date
#      !!! But need to remove the row that has an NA or 0 for the med_amount !!! (only remove med_amount==0 or is NA for the DUPLICATES)
#     this will remove any duplicate entries of the same drug, regardless of whether the med_amount was entered or not
# 2) fill in missing concentrations w/ defaults (find by checking where any remaining med_amounts are still 0 or NA)
#     this will take care of those that are actually missing and don't have the value filled in on the other flowsheet/in a second row
# 3) then, FOR THOSE ROWS THAT HAD THE CONCENTRATIONS FILLED IN, calculate med_amount


#### SO I JUST LEARNED:
# The Controlled Substance Flowsheet always has NAs--don't use it
# so BASICALLY I can just filter out where Record_Type=="30 Controlled Substance F/S"
# The Medication Adm Record has the Dispensed Amt in mg, so the Medication_Amt_mg column is also auto-filled in

# there are still "duplicate" rows with differing Medication_Amts after this, but maybe they're NOT duplicates?
seddat <- seddatM8GS[seddatM8GS$MRN %in% ichmriMRN, ]

# remove absolute duplicate rows
dup_rows <- which(duplicated(seddat))
seddat <- seddat[-dup_rows, ]
# then remove rows that are duplicates of the same drug amount, same flowsheet, same time, same patient
dup_rows2 <- which(duplicated(seddat[,-11]))
seddat <- seddat[-dup_rows2, ]


# 0)
seddat <- seddat %>%
  filter( ! ((Dispensed_Amt==0 | is.na(Dispensed_Amt)) & (Medication_Amt_mg==0 | is.na(Medication_Amt_mg)))    ) %>%
  filter(Record_Type != "30 Controlled Substance F/S")

seddat <- seddat[order(seddat$MRN, seddat$Medication_Name, seddat$Administered_Dtm), ]

# all Med_Amts are NA or 0 in this too, so it's the same as the above
# tmptst <- seddat %>%
#   filter(   Dispensed_Amt==0 | is.na(Dispensed_Amt)  )



### the below might not be necessary anymore, depending on whether or not the duplicates aren't actually duplicates




# 
# # 1)
# ### OK SO NOW THERE ARE DUPLICATE TIME ENTRIES WITH DIFFERING MEDICATION AMOUNTS
# ## after filtering out the controlled substance flowsheet, there's only 12 of these.
# ## I GUESS I just have to assume both numbers are valid?
# # seddat <- seddat[-dup_rows, ]
# MRN_drug_time <- seddat[,c("MRN", "Medication_Name", "Administered_Dtm")]
# # which rows, based on just MRN, Med, and Admin time are duplicated
# MRN_drug_time.dup <- which(duplicated(MRN_drug_time))
# # then pull the MRN+drug+time and its match from the entire table...
# matches <- plyr::match_df(MRN_drug_time, seddat[MRN_drug_time.dup,])
# # ...by using its rownames / row indices
# match_rownms <- rownames(matches)
# # so this has actual entire row duplicates AND not-fully-filled row duplicates (because of the multiple flowsheet entries)
# seddat_dups <- seddat[rownames(seddat) %in% match_rownms, ]
# # order it by MRN, Drug, and Adm_time
# seddat_dups <- seddat_dups[order(seddat_dups$MRN, seddat_dups$Medication_Name, seddat_dups$Administered_Dtm), ]
# 
# 


# Now have to go through each row, find if it's a duplicate of the one before it.
# If it is, then check its Medication_Amt_mg value.
# If it is NA or 0, then check the row before it (the corresponding duplicate)
# and see if ITS Medication_Amt_mg value is also NA or 0







#### Replace missing concentrations with the default amount ####




# # "default" concentrations for each drug
# # Propofol: 10 mg/mL
# # Fentanyl: 10 ug/mL
# # Dexmedetomidine: 4 ug/mL
# # Midazolam: 1 mg/mL




# get indices of where Medication Amt is NA
med_na <- which(is.na(seddat$Medication_Amt_mg))
for (i in 1:length(med_na)) {
  # get the current row
  rw <- seddat[med_na[i], ]
  
  # check if there isn't a med amount but IS a concentration
  if (!is.na(rw$Concentration)) {
    if (rw$Concentration != 1) {
      stop(print(paste("something's f'd with row", med_na[i])))
    }
    
  }
  
  # for propofol, check if dispensed amt units are in mL and then multiply the default concentration of 10 mg/mL to get med amount
  if (grepl("propofol", rw$Medication_Name, ignore.case=T)) {
    if (grepl("ml", seddat[med_na[i], "UM"], ignore.case=T)) {
      print(paste("replacing missing value for Propofol on line", med_na[i]))
      seddat[med_na[i], "Medication_Amt_mg"] <- seddat[med_na[i], "Dispensed_Amt"] * 10 # 10 mg/mL
#      print(seddat[med_na[i], "Dispensed_Amt"] * 10)
    } else {
      stop(print(paste("units are DIFFERENT for propofol row", med_na[i])))
    }
    
  } else if (grepl("fentanyl", rw$Medication_Name, ignore.case=T)) {
    if (grepl("ml", seddat[med_na[i], "UM"], ignore.case=T)) {
      print(paste("replacing missing value for Fentynal on line", med_na[i]))
      seddat[med_na[i], "Medication_Amt_mg"] <- seddat[med_na[i], "Dispensed_Amt"] * 0.01 # 10 ug/mL
#      print(seddat[med_na[i], "Dispensed_Amt"] * 0.01)
    } else {
      stop(print(paste("units are diff for fentalyl row", med_na[i])))
    }
      
    } else if (grepl("midazolam", rw$Medication_Name, ignore.case=T)) {
      if (grepl("ml", seddat[med_na[i], "UM"], ignore.case=T)) {
        print(paste("replacing missing value for Midazolam on line", med_na[i]))
        seddat[med_na[i], "Medication_Amt_mg"] <- seddat[med_na[i], "Dispensed_Amt"] * 1 # 1 mg/mL
#        print(seddat[med_na[i], "Dispensed_Amt"] * 1)
      } else {
        stop(print(paste("units are diff for midaz row", med_na[i])))
      }
    } else if (grepl("dexmedetomidine", rw$Medication_Name, ignore.case=T)) {
      if (grepl("ml", seddat[med_na[i], "UM"], ignore.case=T)) {
        print(paste("replacing missing value for Dexmedetomidine on line", med_na[i]))
        seddat[med_na[i], "Medication_Amt_mg"] <- seddat[med_na[i], "Dispensed_Amt"] * 0.004 # 4 ug/mL
#        print(seddat[med_na[i], "Dispensed_Amt"] * 0.004)
      }
    }  else {
      print(paste("missing info for", rw$Medication_Name, "on line", med_na[i]))
    }
    
  }









# rename. originally due to using sqldf, now just because.
ichmri <- loc_analyse.raw3


# this is 2*the half life (in hours), according to Claassen in Annals of Neurol 2016
Prop_2halflife <- 2
Midaz_2halflife <- 24
Dexmed_2halflife <- 4
Fent_2halflife <- 8
Keta_2halflife <- 6




sedMRI <- ichmri %>%
  select(MRN, MRI_eval_note_date) %>%
  # join all of the sedation data and then filter for the different time ranges for each drug further down
  left_join(select(seddat, MRN, Medication_Name, Administered_Dtm, Medication_Amt_mg), by="MRN") %>%
  group_by(MRN) %>%
  # ..there are sometimes many variations on the same drug. I guess people type it differently?
  
  ## THE SEDATION TABLE DOESN'T HAVE CONCENTRATIONS FOR THE INJECTIONS 
  ########### Therefore, it's not included in this because it's only taking the Medication_Amt column, which is made using those concentrations
  ########### Ok, apparently SOME injections still have this column calculated somehow without a concentration
  ########### OKOk so this is because SOME of the injections have their dispensed amount units in micrograms instead of milliliters.
  
  filter( (Medication_Name=="Propofol +R DRIP" & (Administered_Dtm >= MRI_eval_note_date-(60*60*Prop_2halflife) & Administered_Dtm <= MRI_eval_note_date)) |
          (Medication_Name=="Propofol Inj" & (Administered_Dtm >= MRI_eval_note_date-(60*60*Prop_2halflife) & Administered_Dtm <= MRI_eval_note_date)) |
          (Medication_Name=="Propofol Inj +R+" & (Administered_Dtm >= MRI_eval_note_date-(60*60*Prop_2halflife) & Administered_Dtm <= MRI_eval_note_date)) |
          (Medication_Name=="Midazolam HCL Inj +R+" & (Administered_Dtm >= MRI_eval_note_date-(60*60*Midaz_2halflife) & Administered_Dtm <= MRI_eval_note_date)) |
          (Medication_Name=="Midazolam HCl +R+ DRIP" & (Administered_Dtm >= MRI_eval_note_date-(60*60*Midaz_2halflife) & Administered_Dtm <= MRI_eval_note_date)) |
          (Medication_Name=="Dexmedetomidine Drip 400 microgram/100ml" & (Administered_Dtm >= MRI_eval_note_date-(60*60*Dexmed_2halflife) & Administered_Dtm <= MRI_eval_note_date)) |
          (Medication_Name=="Dexmedetomidine Drip 200 microgram/50ml" & (Administered_Dtm >= MRI_eval_note_date-(60*60*Dexmed_2halflife) & Administered_Dtm <= MRI_eval_note_date)) |
          (Medication_Name=="Dexmedetomidine Drip" & (Administered_Dtm >= MRI_eval_note_date-(60*60*Dexmed_2halflife) & Administered_Dtm <= MRI_eval_note_date)) |
          (Medication_Name=="fentaNYL DRIP" & (Administered_Dtm >= MRI_eval_note_date-(60*60*Fent_2halflife) & Administered_Dtm <= MRI_eval_note_date)) |
          (Medication_Name=="Fentanyl DRIP" & (Administered_Dtm >= MRI_eval_note_date-(60*60*Fent_2halflife) & Administered_Dtm <= MRI_eval_note_date)) |
          (Medication_Name=="fentaNYL DRIP." & (Administered_Dtm >= MRI_eval_note_date-(60*60*Fent_2halflife) & Administered_Dtm <= MRI_eval_note_date)) |
          (Medication_Name=="fentaNYL DRIP 5000 microgram/250ml" & (Administered_Dtm >= MRI_eval_note_date-(60*60*Fent_2halflife) & Administered_Dtm <= MRI_eval_note_date)) |
          (Medication_Name=="Fentanyl Citrate Inj +R+" & (Administered_Dtm >= MRI_eval_note_date-(60*60*Fent_2halflife) & Administered_Dtm <= MRI_eval_note_date)) |
          (Medication_Name=="fentaNYL Citrate Inj +R+" & (Administered_Dtm >= MRI_eval_note_date-(60*60*Fent_2halflife) & Administered_Dtm <= MRI_eval_note_date)) |
          (Medication_Name=="Ketamine DRIP" & (Administered_Dtm >= MRI_eval_note_date-(60*60*Keta_2halflife) & Administered_Dtm <= MRI_eval_note_date)) |
          (Medication_Name=="Ketamine Inj" & (Administered_Dtm >= MRI_eval_note_date-(60*60*Keta_2halflife) & Administered_Dtm <= MRI_eval_note_date))  ) %>%

  #mutate(MRI_hr = hour(MRI_eval_note_date), Adm_hr = hour(Administered_Dtm))
  ## make a new row of the drug-at-time-of-assessment value
  ## by picking the medication amount when the administered time falls within an hour(?) of the mri evaluation time
  #mutate(Administered_Dtm %within% interval(MRI_eval_note_date-(60*60), MRI_eval_note_date))
  
  # %within% is something from the lubridate package.
  # This is checking for the medication amount that falls within an hour range of the time of MRI assessment
  mutate(Assessment_Time_Drug = case_when(Administered_Dtm %within% interval(MRI_eval_note_date-(60*60), MRI_eval_note_date) ~ Medication_Amt_mg,
                                       TRUE ~ 0)) %>%
  # so now make a column of just the dates or group by date(MRI_eval_note_date)?? and then take sum of Medication_Amt_mg and the max of Assessment_Time_Drug
  group_by(MRN, MRI_eval_note_date, Medication_Name) %>%
  summarize(Drug_Sum_mg = sum(Medication_Amt_mg, na.rm=T), Assessment_Time_Drug = max(Assessment_Time_Drug, na.rm=T))



## These are the drugs that are missing concentrations
unique(seddat[which(is.na(seddat$Medication_Amt_mg)), "Medication_Name"])

# apparently only 33 patients had sedation before the MRI assessment
length(unique(sedMRI$MRN))

# remove the row with the obviously wrong large value
ix <- which(sedMRI$Assessment_Time_Drug==max(sedMRI$Assessment_Time_Drug))
sedMRI <- sedMRI[-ix, ]



## Make general drug name column in order to sum over different names of the same drug

########### Right now this is ONLY looking at DRIPS and INJECTIONS THAT WERE IN MICROGRAMS instead of mLs for dispensed amount.
sedMRI <- sedMRI %>%
  mutate(New_Drug_Name = case_when(Medication_Name == "Propofol +R DRIP" ~ "Propofol",
                                   Medication_Name == "fentaNYL DRIP" ~ "Fentanyl",
                                   Medication_Name == "Fentanyl Citrate Inj +R+" ~ "Fentanyl",
                                   Medication_Name == "fentaNYL Citrate Inj +R+" ~ "Fentanyl",
                                   Medication_Name == "Dexmedetomidine Drip 200 microgram/50ml" ~ "Dexmedetomidine",
                                   Medication_Name == "Dexmedetomidine Drip" ~ "Dexmedetomidine",
                                   Medication_Name == "Midazolam HCL Inj +R+" ~ "Midazolam",
                                   Medication_Name == "Midazolam HCl +R+ DRIP" ~ "Midazolam"))





# cast table to have drugs as their own column
sedb4MRI <- dcast(sedMRI, MRN + MRI_eval_note_date ~ New_Drug_Name, value.var="Drug_Sum_mg", fill=0, fun.aggregate = sum)
#test2 <- dcast(test, MRN + MRI_eval_note_date ~ New_Drug_Name, value.var="Drug_Sum_mg", fill=0, fun.aggregate = sum)
sedatMRI <- dcast(sedMRI, MRN + MRI_eval_note_date ~ New_Drug_Name, value.var="Assessment_Time_Drug", fill=0, fun.aggregate = sum)
# ...or maybe not

# join it to the main table
#tstnew <- merge(loc_analyse.raw3, sedb4MRI, by="MRI")

# join to main table (just mri and consciousness columns)
sedConscb4MRI <- merge(ichmri[,c("MRN", "MRI_Cs2")], sedb4MRI[,-2], all.x=T)
sedConscatMRI <- merge(ichmri[,c("MRN", "MRI_Cs2")], sedatMRI[,-2], all.x=T)
# fill in the new NAs with 0s

sedConscb4MRI <- as.data.frame(zoo::na.fill(sedConscb4MRI, 0))
sedConscatMRI <- as.data.frame(zoo::na.fill(sedConscatMRI, 0))
#a <- zoo::na.fill(sedConsc, 0)


# using discharge consciousness
sedConscb4MRI_Dch <- merge(ichmri[,c("MRN", "follow2")], sedb4MRI[,-2], all.x=T)
sedConscatMRI_Dch <- merge(ichmri[,c("MRN", "follow2")], sedatMRI[,-2], all.x=T)
# fill in the new NAs with 0s

sedConscb4MRI_Dch <- as.data.frame(zoo::na.fill(sedConscb4MRI_Dch, 0))
sedConscatMRI_Dch <- as.data.frame(zoo::na.fill(sedConscatMRI_Dch, 0))




#### Merge to main table
b4_merge <- sedConscb4MRI[,-2]
at_merge <- sedConscatMRI[,-2]
names(b4_merge)[-1] <- gsub("$", "_before_MRI", names(b4_merge)[-1])
names(at_merge)[-1] <- gsub("$", "_at_MRI", names(at_merge)[-1])

loc_analyse.raw3 <- merge( merge(loc_analyse.raw3, b4_merge, all.x=T), at_merge, all.x=T  )

#all.equal(loc_analyse.raw3, test[,-c(116:123)])

#### Do this regression with JUST the patients who had sedation too ####
#### I might be diluting it by using all patients


## see which drug, if any, has an effect
drugfit_b4 <- glm(MRI_Cs2 ~ Dexmedetomidine * Fentanyl * Midazolam * Propofol, data=sedConscb4MRI, family = "binomial")
summary(drugfit_b4)
# apparently they don't.
# rather, they don't play an affect in our dataset on predicting consciousness.

drugfit_at <- glm(MRI_Cs2 ~ Dexmedetomidine * Fentanyl * Midazolam * Propofol, data=sedConscatMRI, family = "binomial")
summary(drugfit_at)
# fentanyl is marginally significant when looking at just the amounts given at most an hour before the assessment.
# the fact that the coefficients of the others are positive probablly means that many conscious patients had a value of 0 for the drug.

# yes, this gives the same results as above
# drugfit_at_justdoubles <- glm(MRI_Cs2 ~ Dexmedetomidine + Fentanyl + Midazolam + Propofol + 
#                                 Dexmedetomidine:Fentanyl + Dexmedetomidine:Midazolam + Fentanyl:Midazolam +
#                                 Dexmedetomidine:Propofol + Fentanyl:Propofol + Midazolam:Propofol, data=sedConscatMRI, family = "binomial")
# summary(drugfit_at_justdoubles)
# see if this matches
drugfit_at2 <- glm(MRI_Cs2 ~ Dexmedetomidine_at_MRI * Fentanyl_at_MRI * Midazolam_at_MRI * Propofol_at_MRI, data=loc_analyse.raw3, family="binomial")
summary(drugfit_at2)




drugfit_b4_Dch <- glm(follow2 ~ Dexmedetomidine * Fentanyl * Midazolam * Propofol, data=sedConscb4MRI_Dch, family = "binomial")
summary(drugfit_b4_Dch)
# apparently they don't.
# rather, they don't play an affect in our dataset on predicting consciousness.

drugfit_at_Dch <- glm(follow2 ~ Dexmedetomidine * Fentanyl * Midazolam * Propofol, data=sedConscatMRI_Dch, family = "binomial")
summary(drugfit_at_Dch)




#### Sedation Results ####
## drugs before MRI
# Call:
#   glm(formula = MRI_Cs2 ~ Dexmedetomidine * Fentanyl * Midazolam * 
#         Propofol, family = "binomial", data = sedConscb4MRI)
# 
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -1.6510  -1.4250   0.8626   0.8626   1.8128  
# 
# Coefficients: (9 not defined because of singularities)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                  7.969e-01  1.872e-01   4.256 2.08e-05 ***
# Dexmedetomidine                              2.415e+01  3.979e+01   0.607   0.5440    
# Fentanyl                                    -8.900e+00  4.813e+00  -1.849   0.0644 .  
# Midazolam                                    1.831e+01  1.377e+03   0.013   0.9894    
# Propofol                                    -2.476e-03  2.419e-03  -1.024   0.3060    
# Dexmedetomidine:Fentanyl                            NA         NA      NA       NA    
# Dexmedetomidine:Midazolam                           NA         NA      NA       NA    
# Fentanyl:Midazolam                          -1.036e+02  7.280e+03  -0.014   0.9886    
# Dexmedetomidine:Propofol                            NA         NA      NA       NA    
# Fentanyl:Propofol                            6.061e-02  9.782e-02   0.620   0.5355    
# Midazolam:Propofol                                  NA         NA      NA       NA    
# Dexmedetomidine:Fentanyl:Midazolam                  NA         NA      NA       NA    
# Dexmedetomidine:Fentanyl:Propofol                   NA         NA      NA       NA    
# Dexmedetomidine:Midazolam:Propofol                  NA         NA      NA       NA    
# Fentanyl:Midazolam:Propofol                         NA         NA      NA       NA    
# Dexmedetomidine:Fentanyl:Midazolam:Propofol         NA         NA      NA       NA  



## drugs at time of MRI
# Call:
#   glm(formula = MRI_Cs2 ~ Dexmedetomidine * Fentanyl * Midazolam * 
#         Propofol, family = "binomial", data = sedConscatMRI)
# 
# Deviance Residuals: 
#   Min      1Q  Median      3Q     Max  
# -1.687  -1.514   0.858   0.858   1.748  
# 
# Coefficients: (9 not defined because of singularities)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                                  8.097e-01  1.849e-01   4.378  1.2e-05 ***
# Dexmedetomidine                              3.015e+01  6.239e+01   0.483   0.6288    
# Fentanyl                                    -8.372e+01  4.185e+01  -2.000   0.0455 *  
# Midazolam                                   -4.511e+01  3.105e+03  -0.015   0.9884    
# Propofol                                    -4.369e-03  4.672e-03  -0.935   0.3497    
# Dexmedetomidine:Fentanyl                            NA         NA      NA       NA    
# Dexmedetomidine:Midazolam                           NA         NA      NA       NA    
# Fentanyl:Midazolam                           1.281e+03  7.916e+04   0.016   0.9871    
# Dexmedetomidine:Propofol                            NA         NA      NA       NA    
# Fentanyl:Propofol                            2.100e+00  1.663e+02   0.013   0.9899    
# Midazolam:Propofol                                  NA         NA      NA       NA    
# Dexmedetomidine:Fentanyl:Midazolam                  NA         NA      NA       NA    
# Dexmedetomidine:Fentanyl:Propofol                   NA         NA      NA       NA    
# Dexmedetomidine:Midazolam:Propofol                  NA         NA      NA       NA    
# Fentanyl:Midazolam:Propofol                         NA         NA      NA       NA    
# Dexmedetomidine:Fentanyl:Midazolam:Propofol         NA         NA      NA       NA   




#### Do summary stats on sedation data ####

med_summary <- loc_analyse.raw3 %>%
  select(MRI_Cs2, Dexmedetomidine_before_MRI, Dexmedetomidine_at_MRI,
         Fentanyl_before_MRI, Fentanyl_at_MRI,
         Midazolam_before_MRI, Midazolam_at_MRI,
         Propofol_before_MRI, Propofol_at_MRI) %>%
  group_by(MRI_Cs2) %>%
  summarise_all(funs(median, get_25th, get_75th, .args=list(na.rm=T))[])
  #summarise_all(funs(sum, list(na.rm=T)))

med_summary <- med_summary[,c("MRI_Cs2", sort(names(med_summary))[-which(sort(names(med_summary))=="MRI_Cs2")])]



### Need to take mean of only those patients who had the drug
new_drug_table <- loc_analyse.raw3 %>%
  # select just the sedation columns
  select(MRI_Cs2, Dexmedetomidine_before_MRI, Dexmedetomidine_at_MRI,
         Fentanyl_before_MRI, Fentanyl_at_MRI,
         Midazolam_before_MRI, Midazolam_at_MRI,
         Propofol_before_MRI, Propofol_at_MRI) %>%
  # mutate all of them to be binary (i.e., if value exists then 1, else 0)
  mutate_all(funs(binary=ifelse(.!=0, yes=1, no=0)))





sed_vars <- c("Dexmedetomidine_before_MRI", "Dexmedetomidine_at_MRI",
              "Fentanyl_before_MRI", "Fentanyl_at_MRI",
              "Midazolam_before_MRI", "Midazolam_at_MRI",
              "Propofol_before_MRI", "Propofol_at_MRI")


has_drug_summary <- list()
for (i in 1:length(sed_vars)) {
  assign( sed_vars[i],
          
          loc_analyse.raw3 %>%
            select("MRI_Cs2", sed_vars[i]) %>%
            filter_(paste(sed_vars[i], "!= 0")) %>%
            group_by(MRI_Cs2) %>%
            summarise_(mean=paste("mean(",sed_vars[i],")" ),
                       sd=paste("sd(",sed_vars[i],")"))
            
            )
  
  has_drug_summary[[i]] <- eval(as.name(sed_vars[i]))
  names(has_drug_summary)[i] <- sed_vars[i]
}



# loop over all params and do a wilcox test
wilcox_out_sed <- list()
for (i in 1:length(sed_vars)) {
  form <- as.formula(paste(sed_vars[i], "~ MRI_Cs2"))
  use_dat <- loc_analyse.raw3 %>%
    filter_(paste(sed_vars[i], "!= 0"))
  res <- wilcox.test(form, data=use_dat)
  wilcox_out_sed[[sed_vars[i]]] <- res
}
  





## n(%) for sedation
var_tabs_sed <- list()
for (i in 1:length(sed_vars)) {
  
  assign( sed_vars[i],
          
          loc_analyse.raw3 %>%
            # select just the sedation columns
            select(MRI_Cs2, Dexmedetomidine_before_MRI, Dexmedetomidine_at_MRI,
                   Fentanyl_before_MRI, Fentanyl_at_MRI,
                   Midazolam_before_MRI, Midazolam_at_MRI,
                   Propofol_before_MRI, Propofol_at_MRI) %>%
            # mutate all of them to be binary (i.e., if value exists then 1, else 0)
            mutate_all(funs(ifelse(.!=0, yes=1, no=0))) %>%
            # mutate(had_Dex_before = ifelse(Dexmedetomidine_before_MRI != 0, yes=1, no=0),
            #        had_Dex_after) %>%
            
            # then do the counts and frequency calculation for each sedation column
            group_by_("MRI_Cs2", sed_vars[i]) %>%
            summarize(tot = n()) %>%
            group_by(MRI_Cs2) %>%
            mutate(tot_sum = sum(tot)) %>%
            mutate(percentage = (tot / tot_sum) * 100)
  )
  
  var_tabs_sed[[i]] <- eval(as.name(sed_vars[i]))
  names(var_tabs_sed)[i] <- sed_vars[i]
}





### fisher test of proportions of who has sedation
sed_test_vars <- c("Dexmedetomidine_before_MRI", "Dexmedetomidine_at_MRI",
                   "Fentanyl_before_MRI", "Fentanyl_at_MRI",
                   "Midazolam_before_MRI", "Midazolam_at_MRI",
                   "Propofol_before_MRI", "Propofol_at_MRI")



### Add fisher test results to a new list
fish_out_sed <- list()

for (i in 1:length(sed_test_vars)) {
  
  temp_l <-  paste("tot ~ MRI_Cs2 + ",sed_test_vars[i])
  form <-  as.formula(temp_l)
  
  out <- fisher.test(xtabs(form , data=var_tabs_sed[[i]]), workspace = 100000000)
  
  fish_out_sed[[i]] <- out
  names(fish_out_sed)[i] <- names(var_tabs_sed)[i]
}






### Normalize drug values by mean too
loc_analyse.raw3[,c("Dexmedetomidine_before_MRI_MeanNorm","Fentanyl_before_MRI_MeanNorm","Midazolam_before_MRI_MeanNorm","Propofol_before_MRI_MeanNorm",
         "Dexmedetomidine_at_MRI_MeanNorm","Fentanyl_at_MRI_MeanNorm","Midazolam_at_MRI_MeanNorm","Propofol_at_MRI_MeanNorm")] <- sapply(loc_analyse.raw3[,c("Dexmedetomidine_before_MRI","Fentanyl_before_MRI","Midazolam_before_MRI","Propofol_before_MRI",
                                                                                                                                                "Dexmedetomidine_at_MRI","Fentanyl_at_MRI","Midazolam_at_MRI","Propofol_at_MRI")],
                                                                                                                            meanNormalization)







### Check what the average Fentanyl dose for conscious and unconscious patients is at MRI check time










#### Load Lab data ####

### Using a 4 day window only yielded 1 ALT/AST data for 1 more patient...
labs <- read.xlsx("E:/Experiments/ICH_MRI/All_Labs_4d_Window.xlsx")
names(labs) <- c("MRN", "MRI_Test_Date", "Labs", "value", "units", "Primary_Time", "Tdiff")

## need to select the smallest TDIFF for each drug, for each patient, (for each test date -- 
## should be the same if I exclude test date, since each unique MRN corresponds to a unique test date)
# Also select only the relevant lab types

#LabstoUse <- c("Creat", "BUN", "Gluc", "dirBili", "indirBili", "ALT", "AST")
LabstoUse <- c("Creat", "BUN", "Gluc")

labs.final <- labs %>%
  group_by(MRN, MRI_Test_Date, Labs) %>%
  filter(Tdiff == min(Tdiff) & Labs %in% LabstoUse)


# cast the labs table
labs.c <- dcast(labs.final[,1:4], MRN + MRI_Test_Date ~ Labs, value.var = "value", fun.aggregate = mean, na.rm=T)


# join it to the main table
loc_analyse.raw3 <- merge(loc_analyse.raw3, labs.c, by="MRN", all.x=T)


## Impute the 6 missing lab values for across patients by the average of all patients
if (any(is.na(loc_analyse.raw3[,c("BUN", "Creat", "Gluc")]))) {
  loc_analyse.raw3$BUN[which(is.na(loc_analyse.raw3$BUN))] <- mean(loc_analyse.raw3$BUN, na.rm=T)
  loc_analyse.raw3$Creat[which(is.na(loc_analyse.raw3$Creat))] <- mean(loc_analyse.raw3$Creat, na.rm=T)
  loc_analyse.raw3$Gluc[which(is.na(loc_analyse.raw3$Gluc))] <- mean(loc_analyse.raw3$Gluc, na.rm=T)
  
}


# who's missing all labs
# for (i in 1:nrow(loc_analyse.raw3)) {
#   if ( all(is.na(loc_analyse.raw3[i,135:147])) ) {
#     print(loc_analyse.raw3$MRN[i])
#   }
# }
#which(all(is.na(trerst[,135:147])))



# ## Normalize the lab values
# loc_analyse.raw3[,c("ALT_MeanNorm", "AST_MeanNorm", "BUN_MeanNorm", 
#                     "Creat_MeanNorm", "dirBili_MeanNorm", "Gluc_MeanNorm", "indirBili_MeanNorm")] <- lapply(loc_analyse.raw3[,c("ALT", "AST", "BUN", 
#                                                                                                                                 "Creat", "dirBili", "Gluc", "indirBili")], 
#                                                                                      meanNormalization)


## Normalize only the lab values currently being used
loc_analyse.raw3[,c("BUN_MeanNorm", "Creat_MeanNorm", "Gluc_MeanNorm")] <- lapply(loc_analyse.raw3[,c("BUN", "Creat", "Gluc")], meanNormalization)










#### Load outcome data ####

outcomes <- read.xlsx("E:/Experiments/ICH_MRI/ICH_MRI_Patient_Outcomes.xlsx")

# create column of "all GOSE" by looking at all month follow ups and
# taking the most recent GOSE, if it exists
#outcomes$GOSE_All <- if_else(o)
# outcomes %>%
#   mutate(GOSE_All = case_when(!is.na(GOSE_12M) ~ GOSE_12M,
#                               ))
outcomes <- outcomes %>%
  mutate(GOSE_All = if_else(!is.na(GOSE_12M), GOSE_12M, 
                            if_else(!is.na(GOSE_6M), GOSE_6M,
                                    if_else(!is.na(GOSE_3M), GOSE_3M, GOSE_3M))))

# which has the least amount of missing GOSE scores
length(which(is.na(outcomes$GOSE_3M)))
length(which(is.na(outcomes$GOSE_6M)))
length(which(is.na(outcomes$GOSE_12M)))  # it's the 12 month follow-up

# use the 12 month GOSE follow up
loc_analyse.raw3 <- merge(loc_analyse.raw3, outcomes[,c("MRN", "GOSE_All")], by="MRN", all.x=T)


## Find median mRS at 3 months
median(outcomes$mRS_3M, na.rm=T)
quantile(outcomes$mRS_3M, na.rm=T)
# number dead
length(which(outcomes$mRS_3M == 6))
length(which(outcomes$mRS_3M != 6))
length(which(is.na(outcomes$mRS_3M)))




#### FUNC, GCS, and ICH scores ####
FUNC_GCS_ICH <- read.xlsx("E:/Experiments/ICH_MRI/FUNC_GCS_ICH_Score_158_MRI_patients.xlsx")
# join to loc_analyse.raw3
loc_analyse.raw3 <- merge(loc_analyse.raw3, FUNC_GCS_ICH, by="MRN", all.x=T)
# impute missing FUNC scores with the median
loc_analyse.raw3$FUNC[which(is.na(loc_analyse.raw3$FUNC))] <- median(loc_analyse.raw3$FUNC, na.rm=T)
# impute missing ICH scores with the median too
loc_analyse.raw3$ICH.Score[which(is.na(loc_analyse.raw3$ICH.Score))] <- median(loc_analyse.raw3$ICH.Score, na.rm=T)





#### get all parameters ####
ich_ix <- grep("ICH", names(loc_analyse.raw3))
ed_ix <- grep("edema", names(loc_analyse.raw3))
other_ix <- which(names(loc_analyse.raw3) == "Hg_vol_MeanNorm" | names(loc_analyse.raw3) == "MLS_MeanNorm" | 
                    names(loc_analyse.raw3) == "Ed_vol_MeanNorm" | names(loc_analyse.raw3) == "Old_stroke" | # remove old stroke below?
                    names(loc_analyse.raw3) == "IVH")
med_b4_ix <- which(names(loc_analyse.raw3) == "Dexmedetomidine_before_MRI_MeanNorm" | names(loc_analyse.raw3) == "Fentanyl_before_MRI_MeanNorm" |
                  names(loc_analyse.raw3) == "Midazolam_before_MRI_MeanNorm" | names(loc_analyse.raw3) == "Propofol_before_MRI_MeanNorm" )
## Use this one for sed meds
med_at_ix <- which(names(loc_analyse.raw3) == "Dexmedetomidine_at_MRI_MeanNorm" | names(loc_analyse.raw3) == "Fentanyl_at_MRI_MeanNorm" |
                     names(loc_analyse.raw3) == "Midazolam_at_MRI_MeanNorm" | names(loc_analyse.raw3) == "Propofol_at_MRI_MeanNorm" )

## Lab values
# lab_ix <- which(names(loc_analyse.raw3) %in% c("ALT_MeanNorm", "AST_MeanNorm", "BUN_MeanNorm","Creat_MeanNorm",
#                                                "dirBili_MeanNorm", "Gluc_MeanNorm", "indirBili_MeanNorm"))

lab_ix <- which(names(loc_analyse.raw3) %in% c("BUN_MeanNorm","Gluc_MeanNorm"))  # "Creat_MeanNorm", MIGHT be a problem

## all_ix is ONLY including sedation at time of MRI
## JUST KIDDING it has both now
all_meds <- c(med_b4_ix, med_at_ix)

all_ix <- c(ich_ix, ed_ix, other_ix, all_meds, lab_ix) #med_b4_ix

# all_ix without sedation for "All_Param_Dch"
all_ix_dch <- c(ich_ix, ed_ix, other_ix, lab_ix)




# all params without volumes & MLS
other_ix2 <- which(names(loc_analyse.raw3) == "Old_stroke" |
                    names(loc_analyse.raw3) == "IVH")

all_ix2 <- c(ich_ix, ed_ix, other_ix2, med_at_ix, lab_ix)  #med_b4_ix 



## all params minus sedation and labs
all_ix3 <- c(ich_ix, ed_ix, other_ix)



# remove TH subgroups and MesoC summary params
TH_nms2 <- grep("TH", names(loc_analyse.raw3)[all_ix], value=T)
rm_ix2 <- which(TH_nms2 != "TH_ICH_contro" & TH_nms2 != "TH_ICH_ipsi" & 
                  TH_nms2 != "TH_edema_ipsi" & TH_nms2 != "TH_edema_contro")
rm_ix_nms2 <- TH_nms2[rm_ix2]
## This is the good one.
TH_ix_to_remove2 <- which(names(loc_analyse.raw3)[all_ix] %in% rm_ix_nms2)

#### Add names to the rm_ix_nms2 to remove additional parameters
# need to also remove the ant pons and teg ones and replace them with the C ones
new_rm_nms <- c(rm_ix_nms2, "MesoC_ICH_contro", "MesoC_ICH_ipsi", "BS_ICH", "ICH.Score", "Old_ICH", "Old_stroke",
                "Teg_ICH_contro", "Teg_ICH_ipsi", "Teg_edema_contro", "Teg_edema_ipsi",
                "AntPons_ICH_contro", "AntPons_ICH_ipsi", "AntPons_edema_contro", "AntPons_edema_ipsi")
Extra_ix_to_remove <- which(names(loc_analyse.raw3)[all_ix] %in% new_rm_nms)


## Just get meso and volumes
meso_names <- c("FCx_ICH_ipsi", "FCx_ICH_contro", "FCx_edema_ipsi", "FCx_edema_contro",
                "TH_ICH_ipsi", "TH_ICH_contro", "TH_edema_ipsi", "TH_edema_contro",
                "GP_ICH_ipsi", "GP_ICH_contro", "GP_edema_ipsi", "GP_edema_contro",
                "Caudate_ICH_ipsi", "Caudate_ICH_contro", "Caudate_edema_ipsi", "Caudate_edema_contro",
                "PUT_ICH_ipsi", "PUT_ICH_contro", "PUT_edema_ipsi", "PUT_edema_contro",
                "MB_peduncle_ICH_ipsi", "MB_peduncle_ICH_contro", "MB_peduncle_edema_ipsi", "MB_peduncle_edema_contro",
                "MB_ICH_C",
                "Teg_ICH_ipsi", "Teg_ICH_contro", "Teg_edema_ipsi", "Teg_edema_contro",
                "Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm")

meso_ix <- which(names(loc_analyse.raw3) %in% meso_names)


# mesoC without FC
meso_names_vol <- c("TH_ICH_ipsi", "TH_ICH_contro", "TH_edema_ipsi", "TH_edema_contro",
                 "GP_ICH_ipsi", "GP_ICH_contro", "GP_edema_ipsi", "GP_edema_contro",
                 "Caudate_ICH_ipsi", "Caudate_ICH_contro", "Caudate_edema_ipsi", "Caudate_edema_contro",
                 "PUT_ICH_ipsi", "PUT_ICH_contro", "PUT_edema_ipsi", "PUT_edema_contro",
                 "MB_peduncle_ICH_ipsi", "MB_peduncle_ICH_contro", "MB_peduncle_edema_ipsi", "MB_peduncle_edema_contro",
                 "MB_ICH_C",
                 "Teg_ICH_ipsi", "Teg_ICH_contro", "Teg_edema_ipsi", "Teg_edema_contro",
                 "Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm")

meso_ix_vol <- which(names(loc_analyse.raw3) %in% meso_names_vol)





# Brainstem + Diencephalon + vol
BSD_names_vol <- c("Teg_ICH_C", "Teg_edema_C",
                   #"AntPons_ICH_C", "AntPons_edema_C",
                   "BF_ICH_C", "BF_edema_C", "Hypo_ICH_C", "Hypo_edema_C",
                   "MB_ICH_C", "MB_edema_C",
                    "TH_ICH_ipsi", "TH_ICH_contro", "TH_edema_ipsi", "TH_edema_contro",
                    "GP_ICH_ipsi", "GP_ICH_contro", "GP_edema_ipsi", "GP_edema_contro",
                    "Caudate_ICH_ipsi", "Caudate_ICH_contro", "Caudate_edema_ipsi", "Caudate_edema_contro",
                    "PUT_ICH_ipsi", "PUT_ICH_contro", "PUT_edema_ipsi", "PUT_edema_contro",
                    "MB_peduncle_ICH_ipsi", "MB_peduncle_ICH_contro", "MB_peduncle_edema_ipsi", "MB_peduncle_edema_contro",
                    "Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm")


BSD_ix_vol <- which(names(loc_analyse.raw3) %in% BSD_names_vol)




# Brainstem + Diencephalon + vol + IVH
BSD_names_vol_IVH <- c("Teg_ICH_C", "Teg_edema_C",
                   #"AntPons_ICH_C", "AntPons_edema_C",
                   "BF_ICH_C", "BF_edema_C", "Hypo_ICH_C", "Hypo_edema_C",
                   "MB_ICH_C", "MB_edema_C",
                   "TH_ICH_ipsi", "TH_ICH_contro", "TH_edema_ipsi", "TH_edema_contro",
                   "GP_ICH_ipsi", "GP_ICH_contro", "GP_edema_ipsi", "GP_edema_contro",
                   "Caudate_ICH_ipsi", "Caudate_ICH_contro", "Caudate_edema_ipsi", "Caudate_edema_contro",
                   "PUT_ICH_ipsi", "PUT_ICH_contro", "PUT_edema_ipsi", "PUT_edema_contro",
                   "MB_peduncle_ICH_ipsi", "MB_peduncle_ICH_contro", "MB_peduncle_edema_ipsi", "MB_peduncle_edema_contro",
                   "Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm",
                   "IVH")


BSD_ix_vol_IVH <- which(names(loc_analyse.raw3) %in% BSD_names_vol_IVH)




# Brainstem + Diencephalon
BSD_names <- c("Teg_ICH_C", "Teg_edema_C",
                   #"AntPons_ICH_C", "AntPons_edema_C",
                   "BF_ICH_C", "BF_edema_C", "Hypo_ICH_C", "Hypo_edema_C",
                   "MB_ICH_C", "MB_edema_C",
                   "TH_ICH_ipsi", "TH_ICH_contro", "TH_edema_ipsi", "TH_edema_contro",
                   "GP_ICH_ipsi", "GP_ICH_contro", "GP_edema_ipsi", "GP_edema_contro",
                   "Caudate_ICH_ipsi", "Caudate_ICH_contro", "Caudate_edema_ipsi", "Caudate_edema_contro",
                   "PUT_ICH_ipsi", "PUT_ICH_contro", "PUT_edema_ipsi", "PUT_edema_contro",
                   "MB_peduncle_ICH_ipsi", "MB_peduncle_ICH_contro", "MB_peduncle_edema_ipsi", "MB_peduncle_edema_contro")


BSD_ix <- which(names(loc_analyse.raw3) %in% BSD_names)





# mesoC without FC & w/o vols
meso_names3 <- c("TH_ICH_ipsi", "TH_ICH_contro", "TH_edema_ipsi", "TH_edema_contro",
                 "GP_ICH_ipsi", "GP_ICH_contro", "GP_edema_ipsi", "GP_edema_contro",
                 "Caudate_ICH_ipsi", "Caudate_ICH_contro", "Caudate_edema_ipsi", "Caudate_edema_contro",
                 "PUT_ICH_ipsi", "PUT_ICH_contro", "PUT_edema_ipsi", "PUT_edema_contro",
                 "MB_peduncle_ICH_ipsi", "MB_peduncle_ICH_contro", "MB_peduncle_edema_ipsi", "MB_peduncle_edema_contro",
                 "MB_ICH_C",
                 "Teg_ICH_ipsi", "Teg_ICH_contro", "Teg_edema_ipsi", "Teg_edema_contro")

meso_ix3 <- which(names(loc_analyse.raw3) %in% meso_names3)




# indices of the volume and MLS data
vols_ix <- which(names(loc_analyse.raw3) %in% c("Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm"))

# vols with IVH
vols_IVH_ix <- which(names(loc_analyse.raw3) %in% c("Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm", "IVH"))


# indices for FUNC and GCS score
func_ix <- which(names(loc_analyse.raw3) == "FUNC")
GCS_ix <- which(names(loc_analyse.raw3) == "GCS")
ICHscore_ix <- which(names(loc_analyse.raw3) == "ICH.Score")






### New model of Vols + Meso ONLY for patients that have a lesion in the ipsi thalamus region
## So I think I need to use the indices for vol+meso, but use a subsetted version
## of loc_analyse.raw3 that only has patients who have a 1 for the regions listed below:
ipsi_TH_names <- grep("TH", grep("ipsi", names(loc_analyse.raw3), ignore.case=T, value=T), ignore.case=T, value=T)

loc_analyse.ipsi_TH <- loc_analyse.raw3 %>%
  filter_at(vars(one_of(ipsi_TH_names)), any_vars(.==1))

# seems like there's only 62 patients with a lesion in any of the ipsi TH regions
nrow(loc_analyse.ipsi_TH)








# check for specific areas
#View(filter(loc_analyse.raw3, follow2==0 & ))
ipsi_PUT_names <- grep("PUT", grep("ipsi", names(loc_analyse.raw3), ignore.case=T, value=T), ignore.case=T, value=T)
TH_names <- grep("TH", names(loc_analyse.raw3), ignore.case=T, value=T)
TH_names <- TH_names[which(TH_names != "DEATH_DATE")]

tmp <- loc_analyse.raw3 %>%
  filter(follow2==0) %>%
  filter_at(vars(one_of(ipsi_PUT_names)), any_vars(.==1)) %>%
  filter_at(vars(one_of(TH_names)), all_vars(.==0))
  

# the above does do the same as this:
# tmp2 <- loc_analyse.raw3 %>%
#   filter(follow2==0 & (PUT_edema_ipsi==1 | PUT_ICH_ipsi==1))



# unconscious patients without any lesion in TH regions
uc.noTH <- loc_analyse.raw3 %>%
  filter(follow2==0) %>%
  filter_at(vars(one_of(TH_names)), all_vars(.==0))





#### Find best alpha and lambda for each model ####
lambda.grid <- 10^seq(2, -2, length=100)
#alpha.grid <- seq(0, 1, length=10)
alpha.grid <- seq(0, 1, 0.01)
# 
# 
# # alpha seems to vary between seperate runs of this...
# All_Params_MRI_bestparams <- bestAlphaLambda(loc_analyse.raw3, 
#                                              indep_var_ix=all_ix[-Extra_ix_to_remove], 
#                                              dep_var_vals=loc_analyse.raw3$MRI_Cs2, 
#                                              lambda.grid, alpha.grid)
# 
# 
# 
# # ...........this always chooses a different alpha each time it's run.
# All_Params_MRI_bestparams$bestTune












#### The 4 main models
## MRI_Cs2 ~ All_Params
## MRI_Cs2 ~ Mesocircuit
## follow2 ~ All_Params
## follow2 ~ Mesocircuit
## ...and variations on them


#### New Main Models ####
# all @ MRI
# has all areas, volumes+MLS, sedation meds at time of MRI, lab values (BUN, Gluc, Creat)
All_params_MRI <- RunModel(data=loc_analyse.raw3, 
                           indep_var_ix=all_ix[-Extra_ix_to_remove], 
                           dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                           n_runs=500)

All_params_Dch <- RunModel(data=loc_analyse.raw3, 
                           indep_var_ix=all_ix[-Extra_ix_to_remove], 
                           dep_var_vals=loc_analyse.raw3$follow2,
                           n_runs=500)




# All minus sedation and labs
All_noSedorLabs_MRI <- RunModel(data=loc_analyse.raw3, 
                                indep_var_ix=all_ix3[-Extra_ix_to_remove], 
                                dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                                n_runs=500)

All_noSedorLabs_Dch <- RunModel(data=loc_analyse.raw3, 
                                indep_var_ix=all_ix3[-Extra_ix_to_remove], 
                                dep_var_vals=loc_analyse.raw3$follow2,
                                n_runs=500)



## same model, but completely lasso
All_noSedorLabs_Dch_lasso <- RunModel(data=loc_analyse.raw3, 
                                indep_var_ix=all_ix3[-Extra_ix_to_remove], 
                                dep_var_vals=loc_analyse.raw3$follow2,
                                n_runs=500,
                                alpha=1)



# has just MesoC.
Just_Meso_MRI <- RunModel(data=loc_analyse.raw3, 
                          indep_var_ix=meso_ix3, 
                          dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                          n_runs=500)


Just_Meso_Dch <- RunModel(data=loc_analyse.raw3, 
                          indep_var_ix=meso_ix3, 
                          dep_var_vals=loc_analyse.raw3$follow2,
                          n_runs=500)





# MesoC w/ volumes
Meso_Vol_MRI <- RunModel(data=loc_analyse.raw3, 
                          indep_var_ix=meso_ix_vol, 
                          dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                          n_runs=500)



Meso_Vol_Dch <- RunModel(data=loc_analyse.raw3, 
                          indep_var_ix=meso_ix_vol, 
                          dep_var_vals=loc_analyse.raw3$follow2,
                          n_runs=500)



Meso_Vol_Dch2 <- RunModel(data=loc_analyse.raw3, 
                          indep_var_ix=meso_ix_vol, 
                          dep_var_vals=loc_analyse.raw3$follow2,
                          opposite_vals=loc_analyse.raw3$MRI_Cs2,
                          n_runs=500)




#### BSD Models ####



BSD_Vol_Dch <- RunModel(data=loc_analyse.raw3,
                        indep_var_ix=BSD_ix_vol,
                        dep_var_vals=loc_analyse.raw3$follow2,
                        n_runs=500)


BSD_Dch <- RunModel(data=loc_analyse.raw3,
                        indep_var_ix=BSD_ix,
                        dep_var_vals=loc_analyse.raw3$follow2,
                        n_runs=500)











BSD_Vol_MRI <- RunModel(data=loc_analyse.raw3,
                        indep_var_ix=BSD_ix_vol,
                        dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                        n_runs=500)


BSD_MRI <- RunModel(data=loc_analyse.raw3,
                    indep_var_ix=BSD_ix,
                    dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                    n_runs=500)




# final ones, now with IVH and including volumes, that we're using
BSD_Vol_IVH_Dch <- RunModel(data=loc_analyse.raw3,
                            indep_var_ix=BSD_ix_vol_IVH,
                            dep_var_vals=loc_analyse.raw3$follow2,
                            n_runs=500)


BSD_Vol_IVH_MRI <- RunModel(data=loc_analyse.raw3,
                            indep_var_ix=BSD_ix_vol_IVH,
                            dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                            n_runs=500)


# check these also with a forced 50/50 class split (and need to use 5 folds instead of 10)
BSD_Vol_IVH_Dch_p0.5_nfolds5 <- RunModel(data=loc_analyse.raw3,
                                        indep_var_ix=BSD_ix_vol_IVH,
                                        dep_var_vals=loc_analyse.raw3$follow2,
                                        p=0.5,
                                        nfolds=5,
                                        n_runs=500)


BSD_Vol_IVH_MRI_p0.5_nfolds5 <- RunModel(data=loc_analyse.raw3,
                                        indep_var_ix=BSD_ix_vol_IVH,
                                        dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                                        p=0.5,
                                        nfolds=5,
                                        n_runs=500)


# opposite predict
BSD_Vol_IVH_Dch_oppo <- RunModel(data=loc_analyse.raw3,
                            indep_var_ix=BSD_ix_vol_IVH,
                            dep_var_vals=loc_analyse.raw3$follow2,
                            opposite_vals=loc_analyse.raw3$MRI_Cs2,
                            n_runs=500)













Vols_IVH_MRI <- RunModel(data=loc_analyse.raw3,
                          indep_var_ix=vols_IVH_ix,
                          dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                          n_runs=500)


Vols_IVH_Dch <- RunModel(data=loc_analyse.raw3,
                          indep_var_ix=vols_IVH_ix,
                          dep_var_vals=loc_analyse.raw3$follow2,
                          n_runs=500)










## opposite predict measures for figure 5
BSD_Vol_Dch_oppo <- RunModel(data=loc_analyse.raw3,
                        indep_var_ix=BSD_ix_vol,
                        dep_var_vals=loc_analyse.raw3$follow2,
                        opposite_vals=loc_analyse.raw3$MRI_Cs2,
                        n_runs=500)



#### Attempt to run Meso + Volume model for just patients with lesions in the ipsi TH regions ####
# not sure if want to run on @MRI or @Dch, so do both for now

# I was afraid of this. There's too few observations to do the cross-validation (only 62 patients now)
Meso_Vol_ipsi_TH_MRI <- RunModel(data=loc_analyse.ipsi_TH, 
                         indep_var_ix=meso_ix_vol, 
                         dep_var_vals=loc_analyse.ipsi_TH$MRI_Cs2,
                         n_runs=500)


Meso_Vol_ipsi_TH_Dch <- RunModel(data=loc_analyse.ipsi_TH, 
                          indep_var_ix=meso_ix_vol, 
                          dep_var_vals=loc_analyse.ipsi_TH$follow2,
                          #opposite_vals=loc_analyse.raw3$MRI_Cs2,
                          n_runs=500,
                          nfolds=5)




# meso + vols + ICH score
Meso_Vol_ICHscore_Dch <- RunModel(data=loc_analyse.raw3, 
                          indep_var_ix=c(meso_ix_vol, ICHscore_ix), 
                          dep_var_vals=loc_analyse.raw3$follow2,
                          n_runs=500)

# see how well ICH score does on its own
ichscore_rc <- roc(formula=follow2 ~ ICH.Score, data=loc_analyse.raw3, auc=T, ci=T, na.rm=T)
# so this does very well but, ..of course it does. that's the whole point of this score, right?
# I think looking at how well it does on its own is the wrong way to go about asking
# how relevant it is to be using the other MRI parameters..
print(ichscore_rc)
ci.auc(ichscore_rc, method="bootstrap")


Meso_Vol_Dchp0.5 <- RunModel(data=loc_analyse.raw3, 
                          indep_var_ix=meso_ix_vol, 
                          dep_var_vals=loc_analyse.raw3$follow2,
                          opposite_vals=loc_analyse.raw3$MRI_Cs2,
                          p=0.5,
                          nfolds=5,
                          n_runs=500)

# how many of the "false positives" are actually conscious at dch
length(which(Meso_Vol_Dch2$medProb$class_vals==1 & Meso_Vol_Dch2$medProb$opposite_vals==0 & Meso_Vol_Dch2$medProb$actual_vals==1))
length(which(Meso_Vol_Dch2$medProb$class_vals==0 & Meso_Vol_Dch2$medProb$opposite_vals==0 & Meso_Vol_Dch2$medProb$actual_vals==1))
# and for the other cases
length(which(Meso_Vol_Dch2$medProb$class_vals==1 & Meso_Vol_Dch2$medProb$opposite_vals==1 & Meso_Vol_Dch2$medProb$actual_vals==1))
length(which(Meso_Vol_Dch2$medProb$class_vals==0 & Meso_Vol_Dch2$medProb$opposite_vals==1 & Meso_Vol_Dch2$medProb$actual_vals==1))


# check the patients that are unconscious at MRI to see what their MRI_Cs3 scores are
pred_uncon <- which(Meso_Vol_Dch2$medProb$opposite_vals==0 & Meso_Vol_Dch2$medProb$class_vals==0)
pred_con <- which(Meso_Vol_Dch2$medProb$opposite_vals==0 & Meso_Vol_Dch2$medProb$class_vals==1)

pred_uncon_MRI3 <- loc_analyse.raw3[pred_uncon, "MRI_Cs3"]
pred_con_MRI3 <- loc_analyse.raw3[pred_con, "MRI_Cs3"]

pred_MRI3 <- matrix(nrow=length(pred_uncon_MRI3)+length(pred_con_MRI3), ncol=3)
pred_MRI3$pred_val <- c(rep(0, length(pred_uncon_MRI3)), rep(1, length(pred_con_MRI3)))

table(pred_uncon_MRI3)
table(pred_con_MRI3)


mri3tbl <- as.table(matrix(c(8,4,16,25), nrow=2, ncol=2, byrow = T))
fisher.test(mri3tbl)



# Find best CS during hopsital for these patients:
ix_to_check <- which(Meso_Vol_Dch2$medProb$class_vals==1 & Meso_Vol_Dch2$medProb$opposite_vals==0)
mrns_to_check <- ichmri$MRN[ix_to_check]
# get scores for these patients
best_cs_hosp_check_meso_vol_dch <- new_CS %>%
  filter(MRN %in% mrns_to_check) %>%
  select(MRN, Best_CS_Hosp) %>%
  distinct()

length(which(best_cs_hosp_check_meso_vol_dch$Best_CS_Hosp >= 4))


ix_to_check2 <- which(Meso_Vol_Dch2$medProb$class_vals==0 & Meso_Vol_Dch2$medProb$opposite_vals==0)
mrns_to_check2 <- ichmri$MRN[ix_to_check2]
# get scores for these patients
best_cs_hosp_check_meso_vol_dch2 <- new_CS %>%
  filter(MRN %in% mrns_to_check2) %>%
  select(MRN, Best_CS_Hosp) %>%
  distinct()

length(which(best_cs_hosp_check_meso_vol_dch2$Best_CS_Hosp >= 4))


tab_check <- as.table(matrix(c(5,7,11,30), nrow=2, ncol=2))
fisher.test(tab_check)



# has just volumes+MLS
# mean AUC: 0.7508381
Just_Vols_MRI <- RunModel(data=loc_analyse.raw3,
                     indep_var_ix=vols_ix,
                     dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                     n_runs=500)

## Test running volume model with p=0.5
## need to change the number of folds in the cross-validation,
## otherwise get too few observations per group.
# mean AUC: 0.7440529
Just_Vols_MRI_p0.5 <- RunModel(data=loc_analyse.raw3,
                          indep_var_ix=vols_ix,
                          dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                          p=0.5,
                          nfolds=5,
                          n_runs=500)


# volumes, and GCS
Vols_GCS_MRI <- RunModel(data=loc_analyse.raw3,
                          indep_var_ix=c(vols_ix, GCS_ix),
                          dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                          n_runs=500)


## Just see how well GCS does on its own
gcs_rc <- roc(formula=MRI_Cs2 ~ GCS, data=loc_analyse.raw3, auc=T, ci=T, na.rm=T)
# so this also does very well but, ..of course it does. that's the whole point of this score, right?
# I think looking at how well it does on its own is the wrong way to go about asking
# how relevant it is to be using the other MRI parameters..
print(gcs_rc)
# not sure we need this, as the roc function includes a 95% CI, which is pretty similar to this bootstrapped CI
# or maybe we want the bootstrapped CI instead
ci.auc(gcs_rc, method="bootstrap")


# ok this doesn't work. Just run a logistic regression
## can't run this with only one parameter??
# GCS_MRI <- RunModel(data=loc_analyse.raw3,
#                          indep_var_ix=GCS_ix,
#                          dep_var_vals=loc_analyse.raw3$MRI_Cs2,
#                          n_runs=500)

# well ..I mean, this WILL be significant because we KNOW it already is a predictor.
gcs_mri <- glm(MRI_Cs2 ~ GCS, family="binomial", data=loc_analyse.raw3)
summary(gcs_mri)

# p=0.70
# mean AUC: 0.7486828
# Just_Vols_MRI_p0.7 <- RunModel(data=loc_analyse.raw3,
#                                indep_var_ix=vols_ix,
#                                dep_var_vals=loc_analyse.raw3$MRI_Cs2,
#                                p=0.7,
#                                n_runs=1000)



Just_Vols_Dch <- RunModel(data=loc_analyse.raw3,
                          indep_var_ix=vols_ix,
                          dep_var_vals=loc_analyse.raw3$follow2,
                          n_runs=500)



# has just the sedation meds at time of MRI
Just_Sedation_MRI <- RunModel(data=loc_analyse.raw3, 
                         indep_var_ix=med_at_ix, 
                         dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                         n_runs=500)

Just_Sedation_Dch <- RunModel(data=loc_analyse.raw3, 
                              indep_var_ix=med_at_ix, 
                              dep_var_vals=loc_analyse.raw3$follow2,
                              n_runs=500)



# has just the lab values
# Try custom lambda grid with this and see if it works.
# Yes, this fixes it, and the error is described here:
# https://github.com/lmweber/glmnet-error-example/blob/master/glmnet_error_example.R

# The distribution of AUCs for this model has a kurtosis of like 1e9999999e9999 (it is large.)
# and the spread around 0.5 is so small.
# I guess this is due to the custom lambda sequence, but then I'm not sure how to necessarily fix it.
#lambda_grid=seq(3e-3, 31e-2, length=25)  # this one more closely matches the lambda sequence from one of the model fits
## changing the lambda grid doesn't do anything. just use the one defined above in the "find best alpha" part.
#lambda_grid=10^seq(-4,-1,length=100) #lambda_grid=10^seq(1,-2,length=100)  #lambda_grid=10^seq(10,-2,length=100)
Just_Labs_MRI <- RunModel(data=loc_analyse.raw3, 
                          indep_var_ix=lab_ix, 
                          dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                          n_runs=500,
                          lambda_grid = lambda.grid)



Just_Labs_Dch <- RunModel(data=loc_analyse.raw3, 
                          indep_var_ix=lab_ix, 
                          dep_var_vals=loc_analyse.raw3$follow2,
                          n_runs=500,
                          lambda_grid=lambda.grid)


# what if I desperately try to keep rerunning the model with default lambda_grid
# until it doesn't error out.
# for (i in 1:100000000) {
#   try(
#     test <- RunModel(data=loc_analyse.raw3, 
#                               indep_var_ix=lab_ix, 
#                               dep_var_vals=loc_analyse.raw3$MRI_Cs2,
#                               n_runs=1000), silent=T
#   )
#   if (class(test) == "list") {
#     Just_Labs_MRI <- test
#     break
#   }
# }




# save models into a list so I can loop over them
# ...I don't actually ever use this.
Models <- list(All_params_MRI=All_params_MRI, All_params_Dch=All_params_Dch,
               All_noSedorLabs_MRI=All_noSedorLabs_MRI, All_noSedorLabs_Dch=All_noSedorLabs_Dch,
               Just_Meso_MRI=Just_Meso_MRI, Just_Meso_Dch=Just_Meso_Dch,
               Just_Vols_MRI=Just_Vols_MRI, Just_Vols_Dch=Just_Vols_Dch,
               Just_Sedation_MRI=Just_Sedation_MRI, #Just_Sedation_Dch=Just_Sedation_Dch,
               Just_Labs_MRI=Just_Labs_MRI)#, Just_Labs_Dch=Just_Labs_Dch)






#### Get CS levels for patients ####
true_neg_ix <- All_params_MRI$medProb$patient[which(All_params_MRI$medProb$class_vals==0 & All_params_MRI$medProb$actual_vals==0)]
false_neg_ix <- All_params_MRI$medProb$patient[which(All_params_MRI$medProb$class_vals==0 & All_params_MRI$medProb$actual_vals==1)]
true_pos_ix <- All_params_MRI$medProb$patient[which(All_params_MRI$medProb$class_vals==1 & All_params_MRI$medProb$actual_vals==1)]
false_pos_ix <- All_params_MRI$medProb$patient[which(All_params_MRI$medProb$class_vals==1 & All_params_MRI$medProb$actual_vals==0)]

true_neg_mrn <- loc_analyse.raw3$MRN[true_neg_ix]
false_neg_mrn <- loc_analyse.raw3$MRN[false_neg_ix]
true_pos_mrn <- loc_analyse.raw3$MRN[true_pos_ix]
false_pos_mrn <- loc_analyse.raw3$MRN[false_pos_ix]

# get the scores from the new_CS table
true_pos_CS <- filter(new_CS, MRN %in% true_pos_mrn)
# how many are actually conscious
length(which(true_pos_CS$Best_CS_Hosp >= 4))

true_neg_CS <- filter(new_CS, MRN %in% true_neg_mrn)
length(which(true_neg_CS$Best_CS_Hosp >= 4))

false_pos_CS <- filter(new_CS, MRN %in% false_pos_mrn)
length(which(false_pos_CS$Best_CS_Hosp >= 4))

false_neg_CS <- filter(new_CS, MRN %in% false_neg_mrn)
length(which(false_neg_CS$Best_CS_Hosp >= 4))



## Now, looking back at these numbers, I have no idea where they were pulled from.

# fisher test on the numbers
tab <- as.table(matrix(c(8,6,29,10), nrow=2, ncol=2))
fisher.test(tab)


# test on the GOSE numbers
# (I think these numbers come from the GOSE Distribution plots lower in the script)
tab2 <- as.table(matrix(c(1,13,12,27), nrow=2, ncol=2))
fisher.test(tab2)

# using GOSE of >= 3 instead of >= 4
tab3 <- as.table(matrix(c(4,10,20,19), nrow=2, ncol=2))
fisher.test(tab3)


tab4 <- as.table(matrix(c(4,8,22,19), nrow=2, ncol=2))
fisher.test(tab4)





#### LOO Models ####

All_params_MRI_LOO <- RunLOOModel(data=loc_analyse.raw3, 
                           indep_var_ix=all_ix[-Extra_ix_to_remove], 
                           dep_var_vals=loc_analyse.raw3$MRI_Cs2)

All_params_Dch_LOO <- RunLOOModel(data=loc_analyse.raw3, 
                           indep_var_ix=all_ix[-Extra_ix_to_remove], 
                           dep_var_vals=loc_analyse.raw3$follow2)




# All minus sedation and labs
All_noSedorLabs_MRI_LOO <- RunLOOModel(data=loc_analyse.raw3, 
                                indep_var_ix=all_ix3[-Extra_ix_to_remove], 
                                dep_var_vals=loc_analyse.raw3$MRI_Cs2)

All_noSedorLabs_Dch_LOO <- RunLOOModel(data=loc_analyse.raw3, 
                                indep_var_ix=all_ix3[-Extra_ix_to_remove], 
                                dep_var_vals=loc_analyse.raw3$follow2)



# has just MesoC.
Just_Meso_MRI_LOO <- RunLOOModel(data=loc_analyse.raw3, 
                          indep_var_ix=meso_ix3, 
                          dep_var_vals=loc_analyse.raw3$MRI_Cs2)


Just_Meso_Dch_LOO <- RunLOOModel(data=loc_analyse.raw3, 
                          indep_var_ix=meso_ix3, 
                          dep_var_vals=loc_analyse.raw3$follow2)



# has just volumes+MLS
Just_Vols_MRI_LOO <- RunLOOModel(data=loc_analyse.raw3,
                          indep_var_ix=vols_ix,
                          dep_var_vals=loc_analyse.raw3$MRI_Cs2)



Just_Vols_Dch_LOO <- RunLOOModel(data=loc_analyse.raw3,
                          indep_var_ix=vols_ix,
                          dep_var_vals=loc_analyse.raw3$follow2)



# has just the sedation meds at time of MRI
Just_Sedation_MRI_LOO <- RunLOOModel(data=loc_analyse.raw3, 
                              indep_var_ix=med_at_ix, 
                              dep_var_vals=loc_analyse.raw3$MRI_Cs2)

Just_Sedation_Dch_LOO <- RunLOOModel(data=loc_analyse.raw3, 
                              indep_var_ix=med_at_ix, 
                              dep_var_vals=loc_analyse.raw3$follow2)


Just_Labs_MRI_LOO <- RunLOOModel(data=loc_analyse.raw3, 
                                 indep_var_ix=lab_ix, 
                                 dep_var_vals=loc_analyse.raw3$MRI_Cs2)

Just_Labs_Dch_LOO <- RunLOOModel(data=loc_analyse.raw3, 
                                 indep_var_ix=lab_ix, 
                                 dep_var_vals=loc_analyse.raw3$follow2)











#### other models ####



# all @ MRI
# has all areas, volumes+MLS, sedation meds at time of MRI, lab values
# but n_runs = 1000 instead of 500
All_params_MRI <- RunModel(data=loc_analyse.raw3, 
                           indep_var_ix=all_ix[-Extra_ix_to_remove], 
                           dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                           n_runs=1000)

# All_params_MRI_best_alpha <- RunModel(data=loc_analyse.raw3, 
#                            indep_var_ix=all_ix[-Extra_ix_to_remove], 
#                            dep_var_vals=loc_analyse.raw3$MRI_Cs2,
#                            alpha=All_Params_MRI_bestparams$bestTune$alpha,
#                            n_runs=1000)


# all w/o volumes&MLS @ MRI
# has all areas, sedation meds at time of MRI, and lab values
All_params_noVol_MRI <- RunModel(data=loc_analyse.raw3, 
                                 indep_var_ix=all_ix2[-Extra_ix_to_remove], 
                                 dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                                 n_runs=1000)


# Just sedation meds
# has just the sedation meds at time of MRI
Sedation_MRI <- RunModel(data=loc_analyse.raw3, 
                         indep_var_ix=med_at_ix, 
                         dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                         n_runs=1000)


# mesoC without FC @ MRI
# has just mesoC and volumes+MLS
Meso_MRI <- RunModel(data=loc_analyse.raw3, 
                     indep_var_ix=meso_ix2, 
                     dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                     n_runs=1000)





# all @ Dch
# has all areas, volumes+MLS, sedation meds at time of MRI, and lab values
All_params_Dch <- RunModel(data=loc_analyse.raw3, 
                           indep_var_ix=all_ix[-Extra_ix_to_remove], 
                           dep_var_vals=loc_analyse.raw3$follow2,
                           n_runs=1000)





# all w/o volumes&MLS @ Dch
All_params_noVol_Dch <- RunModel(data=loc_analyse.raw3, 
                                 indep_var_ix=all_ix2[-Extra_ix_to_remove], 
                                 dep_var_vals=loc_analyse.raw3$follow2,
                                 n_runs=1000)


# mesoC without FC @ Dch
Meso_Dch <- RunModel(data=loc_analyse.raw3, 
                     indep_var_ix=meso_ix2, 
                     dep_var_vals=loc_analyse.raw3$follow2,
                     n_runs=1000)


# vols @ MRI
Vols_MRI <- RunModel(data=loc_analyse.raw3,
                     indep_var_ix=vols_ix,
                     dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                     n_runs=1000)


# vols @ Dch
Vols_Dch <- RunModel(data=loc_analyse.raw3,
                     indep_var_ix=vols_ix,
                     dep_var_vals=loc_analyse.raw3$follow2,
                     n_runs=1000)


# Just Meso @ MRI
Just_Meso_MRI <- RunModel(data=loc_analyse.raw3, 
                          indep_var_ix=meso_ix3, 
                          dep_var_vals=loc_analyse.raw3$MRI_Cs2,
                          n_runs=1000)


# Just Meso @ Dch
Just_Meso_Dch <- RunModel(data=loc_analyse.raw3, 
                          indep_var_ix=meso_ix3, 
                          dep_var_vals=loc_analyse.raw3$follow2,
                          n_runs=1000)










#### Plots of the 4 models together ####
# histdat <- data.frame(All_params_MRI$AUC,
#                       All_params_Dch$AUC,
#                       Meso_MRI$AUC,
#                       Meso_Dch$AUC,
#                       Vols_MRI$AUC,
#                       Vols_Dch$AUC,
#                       Just_Meso_MRI$AUC,
#                       Just_Meso_Dch$AUC,
#                       All_params_noVol_Dch$AUC,
#                       All_params_noVol_MRI$AUC,
#                       Sedation_MRI$AUC)


### Put AUCs of all models together here
histdat <- data.frame(All_params_MRI$AUC,
                      All_params_Dch$AUC,
                      All_noSedorLabs_MRI$AUC,
                      All_noSedorLabs_Dch$AUC,
                      Just_Meso_MRI$AUC,
                      Just_Meso_Dch$AUC,
                      Just_Vols_MRI$AUC,
                      Just_Vols_Dch$AUC,
                      Meso_Vol_MRI$AUC,
                      Meso_Vol_Dch$AUC)
                      #Just_Sedation_MRI$AUC)#,
                      #Just_Sedation_Dch$AUC,
                      #Just_Labs_MRI$AUC)
                      #Just_Labs_Dch$AUC)


# just the new models
histdat <- data.frame(BSD_Vol_Dch$AUC,
                      BSD_Dch$AUC,
                      BSD_Vol_MRI$AUC,
                      BSD_MRI$AUC,
                      Just_Vols_Dch$AUC,
                      Just_Vols_MRI$AUC)


# models with IVH included. Also include the 'all params' models
histdat <- data.frame(BSD_Vol_Dch$AUC,
                      BSD_Dch$AUC,
                      BSD_Vol_MRI$AUC,
                      BSD_MRI$AUC,
                      Just_Vols_Dch$AUC,
                      Just_Vols_MRI$AUC,
                      BSD_Vol_IVH_Dch$AUC,
                      BSD_Vol_IVH_MRI$AUC,
                      Vols_IVH_MRI$AUC,
                      Vols_IVH_Dch$AUC)




## Ones we're interested in for the most recent paper submission
histdat <- data.frame(All_params_MRI$AUC,
                      All_params_Dch$AUC,
                      BSD_Vol_IVH_Dch$AUC,
                      BSD_Vol_IVH_MRI$AUC)



histdatNOLABS <- data.frame(All_params_MRI$AUC,
                      All_params_Dch$AUC,
                      All_noSedorLabs_MRI$AUC,
                      All_noSedorLabs_Dch$AUC,
                      Just_Meso_MRI$AUC,
                      Just_Meso_Dch$AUC,
                      Just_Vols_MRI$AUC,
                      Just_Vols_Dch$AUC,
                      Just_Sedation_MRI$AUC)







# # just the discharge models
# histdatDch <- data.frame(All_params_Dch$AUC,
#                          Meso_Dch$AUC,
#                          Vols_Dch$AUC,
#                          Just_Meso_Dch$AUC,
#                          All_params_noVol_Dch$AUC)


# remove ",AUC" from variable names and then melt
colnames(histdat) <- gsub(".AUC", "", colnames(histdat))
histdat2 <- melt(histdat)


colnames(histdatNOLABS) <- gsub(".AUC", "", colnames(histdatNOLABS))
histdat2NOLABS <- melt(histdatNOLABS)

# colnames(histdatDch) <- gsub(".AUC", "", colnames(histdatDch))
# histdatDch2 <- melt(histdatDch)


# histograms are kind of messy, with large or small bin sizes
ggplot(data=histdat2, aes(x=value, fill=variable)) +
  geom_histogram(alpha=0.2, position="identity", bins=100) + theme_minimal()


# kernel density plots look much nicer than histograms for these data
ggplot(data=histdat2NOLABS, aes(x=value, fill=variable)) +
  geom_density(alpha=0.3, position="identity") + 
  xlab("AUC") + theme_minimal()

# this one does indeed look weird
ggplot(data=filter(histdat2, variable=="Just_Labs_MRI"), aes(x=value, fill=variable)) +
  geom_density(alpha=0.3, position="identity") + 
  xlab("AUC") + theme_minimal()

ggplot(data=histdatDch2, aes(x=value, fill=variable)) +
  geom_density(alpha=0.3, position="identity") + 
  xlab("AUC") + theme_minimal()


# violin plot
### I guess we COULD put ALL the models into a single violin plot
### since this is basically just a boxplot + density plot combined
### and allows for non-overlap of all the models
### but it's a bit ...less intuitive to look at?
ggplot(data=histdat2, aes(x=reorder(variable, value), y=value, fill=variable)) +
  geom_violin(alpha=0.7, position="identity", draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_hline(yintercept = 0.5) +
  #geom_jitter(height=0, width=0.1, alpha=0.1) +
  ylab("AUC") + xlab("Model") + theme_minimal() +
  coord_flip()



## split the violin plot into 2 plots: one with MRI models and the other with Dch models
histdat.newnms <- histdat2
histdat.newnms$variable <- gsub("noSedorLabs", "no_Sed_or_Labs", histdat.newnms$variable)
# MRI
ggplot(data=filter(histdat.newnms, grepl("_MRI", variable)), aes(x=reorder(variable, value), y=value, fill=variable)) +
  geom_violin(alpha=0.7, position="identity", draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_hline(yintercept = 0.5) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=9)) +
  #geom_jitter(height=0, width=0.1, alpha=0.1) +
  ylab("AUC") + xlab("Model") + theme_minimal() + labs(fill="Model") +
  ggtitle("Model Comparisons", subtitle="Predicting Consciousness at time of MRI") +
  coord_flip()


## Dch
ggplot(data=filter(histdat.newnms, grepl("_Dch", variable)), aes(x=reorder(variable, value), y=value, fill=variable)) +
  geom_violin(alpha=0.7, position="identity", draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_hline(yintercept = 0.5) +
  scale_y_continuous(breaks = scales::pretty_breaks(n=9)) +
  #geom_jitter(height=0, width=0.1, alpha=0.1) +
  ylab("AUC") + xlab("Model") + theme_minimal() + labs(fill="Model") +
  ggtitle("Model Comparisons", subtitle="Predicting Consciousness at discharge") +
  coord_flip()





# # vol
# ggplot(data=filter(histdat.newnms, grepl("_Vol", variable)), aes(x=reorder(variable, value), y=value, fill=variable)) +
#   geom_violin(alpha=0.7, position="identity", draw_quantiles = c(0.25, 0.5, 0.75)) +
#   geom_hline(yintercept = 0.5) +
#   scale_y_continuous(breaks = scales::pretty_breaks(n=9)) +
#   #geom_jitter(height=0, width=0.1, alpha=0.1) +
#   ylab("AUC") + xlab("Model") + theme_minimal() + labs(fill="Model") +
#   ggtitle("Model Comparisons", subtitle="Predicting Consciousness at time of MRI") +
#   coord_flip()
# 
# 
# ## no vol
# ggplot(data=filter(histdat.newnms, !grepl("_Vol", variable)), aes(x=reorder(variable, value), y=value, fill=variable)) +
#   geom_violin(alpha=0.7, position="identity", draw_quantiles = c(0.25, 0.5, 0.75)) +
#   geom_hline(yintercept = 0.5) +
#   scale_y_continuous(breaks = scales::pretty_breaks(n=9)) +
#   #geom_jitter(height=0, width=0.1, alpha=0.1) +
#   ylab("AUC") + xlab("Model") + theme_minimal() + labs(fill="Model") +
#   ggtitle("Model Comparisons", subtitle="Predicting Consciousness at discharge") +
#   coord_flip()









#### GOSE bar plots for patients in confusion matrix ####
All_params_MRI$confus_mat
View(All_params_MRI$medProb)

# interested in mainly these two
fp_pat <- which( All_params_MRI$medProb$class_vals==1 & All_params_MRI$medProb$actual_vals==0 )
tn_pat <- which( All_params_MRI$medProb$class_vals==0 & All_params_MRI$medProb$actual_vals==0 )

# but also plot GOSE for these two as well
tp_pat <- which( All_params_MRI$medProb$class_vals==1 & All_params_MRI$medProb$actual_vals==1 )
fn_pat <- which( All_params_MRI$medProb$class_vals==0 & All_params_MRI$medProb$actual_vals==1 )

# create new subset table with the values for these patients,
# and add an identifier column of each classification grouping
## (I've also just learned that you can keep a single column from a data frame AS a data frame if you use 'drop=FALSE')
fp_tab <- loc_analyse.raw3[fp_pat, "GOSE_All", drop=FALSE]
fp_tab$class_grp <- "false_pos (n=38)"

tn_tab <- loc_analyse.raw3[tn_pat, "GOSE_All", drop=FALSE]
tn_tab$class_grp <- "true_neg (n=15)"

tp_tab <- loc_analyse.raw3[tp_pat, "GOSE_All", drop=FALSE]
tp_tab$class_grp <- "true_pos (n=96)"

fn_tab <- loc_analyse.raw3[fn_pat, "GOSE_All", drop=FALSE]
fn_tab$class_grp <- "false_neg (n=9)"


# rbind them all together
GOSE_class <- rbind(fp_tab, tn_tab, tp_tab, fn_tab)
# convert NAs to "Missing"
#GOSE_class$GOSE_All[which(is.na(GOSE_class$GOSE_All))] <- "Missing"

# or recode as factor instead and don't drop NAs.
# this way, the NAs are just grey
GOSE_class$GOSE_All <- factor(GOSE_class$GOSE_All, exclude=NULL)







# this is the same as GOSE_class now
test <- makeGOSEplot(All_params_MRI$medProb, "class_vals", "actual_vals")

# one for the "weird" confusion matrix -- using MRI values as the actual values
MesoVolGOSEdat <- makeGOSEplot(Meso_Vol_Dch2$medProb, "class_vals", "opposite_vals")

#how many faslepos(?) have gose >=4 / <4 (it's how many predicted conscious using cs @ dch?)


# one for the new BSD model
BSDVolGOSEdat <- makeGOSEplot(BSD_Vol_Dch_oppo$medProb, "class_vals", "opposite_vals")

# and now one for the BSD + vol + IVH model
BSDVolIVHGOSEdat <- makeGOSEplot(BSD_Vol_IVH_Dch_oppo$medProb, "class_vals", "opposite_vals")




#### Make contingency table of predicted Dch status vs. GOSE >=4 / <4 (old way before using function) ####

## true neg, false pos
## find gose >=4 and < 4
# so where class_val == 1 and class_val == 0,
# find how many have gose >= 4 and gose < 4
MesoVolGOSEdatSummary <- MesoVolGOSEdat %>%
  filter(class_grp %in% c("false_pos", "true_neg")) %>%
  mutate(gose_grp=as.factor(ifelse(as.numeric(as.character(GOSE_All)) >= 4, 1, 0)),
         class_grp=as.factor(class_grp)) %>%
  filter(!is.na(gose_grp)) %>%
  group_by(class_grp, gose_grp) #%>%
  # or do I just want to save the table up to the above line and then run table() on it
  # summarise(count = n()) %>%
  # tidyr::complete(gose_grp, fill=list(count=0)) %>%
  # distinct()
  
MesoVolGOSEtable <- table(MesoVolGOSEdatSummary$gose_grp, MesoVolGOSEdatSummary$class_grp)
fisher.test(MesoVolGOSEtable)


MesoVolGOSEcounts <- MesoVolGOSEdat %>%
  filter(class_grp %in% c("false_pos", "true_neg") & !is.na(GOSE_All)) %>%
  group_by(class_grp, GOSE_All) %>%
  summarise(count = n())

# table of just the 53 patients
MesoVolGOSEdatJust53 <- MesoVolGOSEdat %>%
  filter(class_grp %in% c("false_pos", "true_neg"))
  

MesoVolGOSE_reshape <- reshape2::dcast(as.data.frame(MesoVolGOSEcounts), GOSE_All ~ class_grp, value.var="count", fill=0)




# do it for BSD
BSDVolGOSEdatSummary <- BSDVolGOSEdat %>%
  filter(class_grp %in% c("false_pos", "true_neg")) %>%
  mutate(gose_grp=as.factor(ifelse(as.numeric(as.character(GOSE_All)) >= 4, 1, 0)),
         class_grp=as.factor(class_grp)) %>%
  filter(!is.na(gose_grp)) %>%
  group_by(class_grp, gose_grp)


BSDVolGOSEtable <- table(BSDVolGOSEdatSummary$gose_grp, BSDVolGOSEdatSummary$class_grp)
fisher.test(BSDVolGOSEtable)



# and now BSD Vol IVH
BSDVolIVHGOSEdatSummary <- BSDVolIVHGOSEdat %>%
  filter(class_grp %in% c("false_pos", "true_neg")) %>%
  mutate(gose_grp=as.factor(ifelse(as.numeric(as.character(GOSE_All)) >= 4, 1, 0)),
         class_grp=as.factor(class_grp)) %>%
  filter(!is.na(gose_grp)) %>%
  group_by(class_grp, gose_grp)


BSDVolIVHGOSEtable <- table(BSDVolIVHGOSEdatSummary$gose_grp, BSDVolIVHGOSEdatSummary$class_grp)
fisher.test(BSDVolIVHGOSEtable)



#### GOSE plots (old way before using function) ####

ggplot(GOSE_class, aes(x=class_grp, fill=GOSE_All)) + 
  geom_bar(stat="count", position="dodge") +
  scale_y_continuous(breaks = scales::pretty_breaks(n=10)) +
  xlab("Classification Group") +
  ggtitle("GOSE Distribution per Confusion Matrix Classification", 
          subtitle="true positive = conscious; true negative = unconscious; model = All_params_MRI") +
  theme_bw() +
  coord_flip()


ggplot(test, aes(x=class_grp, fill=GOSE_All)) + 
  geom_bar(stat="count", position="dodge") +
  scale_y_continuous(breaks = scales::pretty_breaks(n=10)) +
  xlab("Classification Group") +
  ggtitle("GOSE Distribution per Confusion Matrix Classification", 
          subtitle="true positive = conscious; true negative = unconscious; model = All_params_MRI") +
  theme_bw() +
  coord_flip()




ggplot(MesoVolGOSEdat, aes(x=class_grp, fill=GOSE_All)) + 
  geom_bar(stat="count", position="dodge") +
  scale_y_continuous(breaks = scales::pretty_breaks(n=10)) +
  xlab("Classification Group") +
  ggtitle("GOSE Distribution per Meso_Vol 'opposite' Confusion Matrix Classification", 
          subtitle="true positive = conscious; true negative = unconscious; model = Meso_Vol_Dch2") +
  theme_bw() +
  coord_flip()



ggplot(BSDVolGOSEdat, aes(x=class_grp, fill=GOSE_All)) + 
  geom_bar(stat="count", position="dodge") +
  scale_y_continuous(breaks = scales::pretty_breaks(n=10)) +
  xlab("Classification Group") +
  ggtitle("GOSE Distribution per BSD_Vol 'opposite' Confusion Matrix Classification", 
          subtitle="true positive = conscious; true negative = unconscious; model = BSD_Vol_Dch") +
  theme_bw() +
  coord_flip()



ggplot(BSDVolIVHGOSEdat, aes(x=class_grp, fill=GOSE_All)) + 
  geom_bar(stat="count", position="dodge") +
  scale_y_continuous(breaks = scales::pretty_breaks(n=10)) +
  xlab("Classification Group") +
  ggtitle("GOSE Distribution per BSD_Vol_IVH 'opposite' Confusion Matrix Classification", 
          subtitle="true positive = conscious; true negative = unconscious; model = BSD_Vol_IVH_Dch") +
  theme_bw() +
  coord_flip()




#### Fisher and GOSE plots using function ####
BSDVolIVHGOSE_outInfo <- GOSEfisherAndPlot(BSD_Vol_IVH_Dch_oppo)

BSDVolIVHGOSE_outInfo$contingTable
BSDVolIVHGOSE_outInfo$fisher
BSDVolIVHGOSE_outInfo$the_plot







#### GOSE and FUNC scores for the 53 unconscious patients at MRI ####
GOSE_FUNC_uncon_MRI <- loc_analyse.raw3 %>%
  filter(MRI_Cs2==0) %>%
  #select(MRN, FUNC, GOSE_All) %>%
  mutate(GOSE_grp = ifelse(GOSE_All >= 5, 1, 0), #ifelse(!is.na(GOSE_All), 0, NA)  # don't actually need to account for any NAs
         FUNC_grp = ifelse(FUNC >= 9, 1, 0)) %>% #ifelse(!is.na(FUNC), 0, NA)
  select(MRN, GOSE_All, GOSE_grp, FUNC, FUNC_grp)


GOSE_FUNC_table <- xtabs(~ GOSE_grp + FUNC_grp, data=GOSE_FUNC_uncon_MRI)
print(GOSE_FUNC_table)

fisher.test(GOSE_FUNC_table)

GOSE_FUNC_nondichot <- GOSE_FUNC_uncon_MRI %>%
  mutate(FUNC=ifelse(FUNC >=9, 1, 0)) %>%
  group_by(FUNC, GOSE_All) %>%
  summarise(count = n())

GOSE_FUNC_reshape <- reshape2::dcast(as.data.frame(GOSE_FUNC_nondichot), GOSE_All ~ FUNC, value.var="count", fill=0)



#### COmpare GOSE/FUNC and GOSE/model predict CS @ Dch ####
# can see by how much the odds ratios differ
fisher.test(GOSE_FUNC_table)
fisher.test(MesoVolGOSEtable)


# calculate accuracy, etc. for each confusion matrix:


CMstats(GOSE_FUNC_table)
# these metrics don't really mean anything for this table,
# since it's not comparing the predicted outcome to actual outcome
# -- it's comparing predicted outcome (trained using CS @ dch) to actual outcome it wasn't trained on (CS @ MRI)
CMstats(MesoVolGOSEtable)


#mantelhaen test (no idea how to set up the input data other than manually.)
input_man <- array(c(26,4,10,9,25,13,11,0),
                   dim=c(2,2,2),
                   dimnames=list(GOSE = c("0", "1"),
                                 Consc = c("0", "1"),
                                 tablename = c("FUNC", "Model")))

# this finds that the true common odds ratio is 1.
mantelhaen.test(input_man)






### try a logistici rgreeirsion instead
GOSE_FUNC_uncon_MRI
MesoVolGOSEdatJust53

# combine these tables and dichotomize relevant variables
cm_compare_cb <- merge(GOSE_FUNC_uncon_MRI, MesoVolGOSEdatJust53[,c("MRN","class_grp")], by="MRN")

cm_compare_fit <- glm(GOSE_grp ~ FUNC_grp + class_grp, data=cm_compare_cb, family="binomial")
summary(cm_compare_fit)

# so this shows that the class_groups from the model aren't significant in predicting
# although, the coefficient is much larger for the class_group than the func score (-17 vs 2),
# but the slopes aren't really comparable since one is categorical and the other is numeric (and binary)
# even though the p value for class_group is basically 1.
# the GOSE group (>=4 or <4),
# BUT I think the GOSE group is extrapolated too far from the original output of the model
# to mean anything.
# They are the scores of the patients who the model predicted should be consc or uncosc @ dch, but comparing to their CS @ MRI
# which the model wasn't informed about at all. Like, I get that's the point, to see what those patients prior status was,
# but it still feels weird to be analyzing model output based off of parameters it didn't even use.


#### Confidence Intervals and avg AUC for the 4 models ####
# ...I don't even think I needed the list of all models.
## save this output to another text file
# print("Average AUC & 95% CI for All_params_MRI")
# mean(histdat$All_params_MRI)
# calcConfInt(histdat$All_params_MRI)
# 
# print("Average AUC & 95% CI for All_params_Dch")
# mean(histdat$All_params_Dch)
# calcConfInt(histdat$All_params_Dch)
# 
# print("Average AUC & 95% CI for All_noSedorLabs_MRI")
# mean(histdat$All_noSedorLabs_MRI)
# calcConfInt(histdat$All_noSedorLabs_MRI)
# 
# 
# calcConfInt(histdat$Meso_MRI)
# calcConfInt(histdat$Meso_Dch)

sink("E:/Experiments/ICH_MRI/Mean_AUC_and_CI_BSD_and_All_Params_Dch_without_sedation_10-2-2018.txt")
for (i in 1:ncol(histdat)) {
  print(paste("Average AUC & 95% CI for", names(histdat)[i]))
  print(mean(histdat[,i]))
  print(calcConfInt(histdat[,i]))
  cat('\n')
}
sink()




# for new models with IVH
sink("E:/Experiments/ICH_MRI/Mean_AUC_and_CI_BSD_vol_IVH.txt")
for (i in 1:ncol(newmods)) {
  print(paste("Average AUC & 95% CI for", names(newmods)[i]))
  print(mean(newmods[,i]))
  print(calcConfInt(newmods[,i]))
  cat('\n')
}
sink()










#### ANOVA on the AUCs? ####
# I mean.. they're normally distributed continuous values... so why not..?
aovdat <- melt(histdat)

aovfit <- aov(value ~ variable, data=aovdat)
anova(aovfit)
TukeyHSD(aovfit)










#### Boxplots of volumes for each conscious state ####
path <- "E:/Experiments/ICH_MRI/figures"

# 3 levels
hgplt_mri_Cs3 <- MakePlot(loc_analyse.raw3, "MRI_Cs3", "Hg_vol", path=path, save=F)
edplt_mri_Cs3 <- MakePlot(loc_analyse.raw3, "MRI_Cs3", "Ed_vol", path=path, save=F)
mlsplt_mri_Cs3 <- MakePlot(loc_analyse.raw3, "MRI_Cs3", "MLS", path=path, save=F)

hgplt_dch_Cs3 <- MakePlot(loc_analyse.raw3, "follow3", "Hg_vol", path=path, save=F)
edplt_dch_Cs3 <- MakePlot(loc_analyse.raw3, "follow3", "Ed_vol", path=path, save=F)
mlsplt_dch_Cs3 <- MakePlot(loc_analyse.raw3, "follow3", "MLS", path=path, save=F)


# 2 levels
hgplt_mri_Cs2 <- MakePlot(loc_analyse.raw3, "MRI_Cs2", "Hg_vol", path=path, save=F)
edplt_mri_Cs2 <- MakePlot(loc_analyse.raw3, "MRI_Cs2", "Ed_vol", path=path, save=F)
mlsplt_mri_Cs2 <- MakePlot(loc_analyse.raw3, "MRI_Cs2", "MLS", path=path, save=F)

hgplt_dch_Cs2 <- MakePlot(loc_analyse.raw3, "follow2", "Hg_vol", path=path, save=F)
edplt_dcj_Cs2 <- MakePlot(loc_analyse.raw3, "follow2", "Ed_vol", path=path, save=F)
mlsplt_dch_Cs2 <- MakePlot(loc_analyse.raw3, "follow2", "MLS", path=path, save=F)



# show the plots
hgplt_mri_Cs3 
edplt_mri_Cs3 
mlsplt_mri_Cs3
hgplt_dch_Cs3 
edplt_dch_Cs3 
mlsplt_dch_Cs3
hgplt_mri_Cs2 
edplt_mri_Cs2 
mlsplt_mri_Cs2
hgplt_dch_Cs2 
edplt_dcj_Cs2 
mlsplt_dch_Cs2










#### Boxplots of coefficients #### (maybe make these violin plots?)
# All params - MRI
ggplot(filter(melt(All_params_MRI$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(All_params_MRI$weights)$value), max(melt(All_params_MRI$weights)$value), by=50), 1)) +
  ggtitle("Parameter Coefficients: All Parameters MRI") +
  coord_flip()


# MesoC - MRI
ggplot(filter(melt(Meso_MRI$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(Meso_MRI$weights)$value), max(melt(Meso_MRI$weights)$value), by=1), 1)) +
  ggtitle("Parameter Coefficients: Mesocircuit MRI") +
  coord_flip()


# All params - Dch
ggplot(filter(melt(All_params_Dch$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(All_params_Dch$weights)$value), max(melt(All_params_Dch$weights)$value), by=50), 1)) +
  ggtitle("Parameter Coefficients: All Parameters Dch") +
  coord_flip()


# MesoC - Dch
ggplot(filter(melt(Meso_Dch$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(Meso_Dch$weights)$value), max(melt(Meso_Dch$weights)$value), by=1), 1)) +
  ggtitle("Parameter Coefficients: Mesocircuit Dch") +
  coord_flip()



ggplot(filter(melt(All_params_noVol_Dch$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(All_params_noVol_Dch$weights)$value), max(melt(All_params_noVol_Dch$weights)$value), by=50), 1)) +
  ggtitle("Parameter Coefficients: all params no volume discharge") +
  coord_flip()




ggplot(filter(melt(All_params_MRI_best_alpha$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(All_params_MRI_best_alpha$weights)$value), max(melt(All_params_MRI_best_alpha$weights)$value), by=1), 1)) +
  ggtitle("Parameter Coefficients: all params no volume discharge") +
  coord_flip()




ggplot(filter(melt(Sedation_MRI$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(Sedation_MRI$weights)$value), max(melt(Sedation_MRI$weights)$value), by=1), 1)) +
  ggtitle("Parameter Coefficients: sedation") +
  coord_flip()









#### New Model weights ####
## still looks bad.
ggplot(filter(melt(All_params_MRI$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(All_params_MRI$weights)$value), max(melt(All_params_MRI$weights)$value), by=50), 1)) +
  ggtitle("Parameter Coefficients: All params MRI") +
  coord_flip()



ggplot(filter(melt(All_params_Dch$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(All_params_Dch$weights)$value), max(melt(All_params_Dch$weights)$value), by=50), 1)) +
  ggtitle("Parameter Coefficients: All params Dch") +
  coord_flip()



ggplot(filter(melt(All_noSedorLabs_MRI$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(All_noSedorLabs_MRI$weights)$value), max(melt(All_noSedorLabs_MRI$weights)$value), by=50), 1)) +
  ggtitle("Parameter Coefficients: All params w/o Sedation or Labs MRI") +
  coord_flip()



ggplot(filter(melt(All_noSedorLabs_Dch$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(All_noSedorLabs_Dch$weights)$value), max(melt(All_noSedorLabs_Dch$weights)$value), by=50), 1)) +
  ggtitle("Parameter Coefficients: All params w/o Sedation or Labs Dch") +
  coord_flip()



ggplot(filter(melt(All_noSedorLabs_Dch_lasso$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(All_noSedorLabs_Dch_lasso$weights)$value), max(melt(All_noSedorLabs_Dch_lasso$weights)$value), by=50), 1)) +
  ggtitle("Parameter Coefficients: All params w/o Sedation or Labs Dch - alpha=1") +
  coord_flip()



ggplot(filter(melt(Just_Meso_MRI$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(Just_Meso_MRI$weights)$value), max(melt(Just_Meso_MRI$weights)$value), by=50), 1)) +
  ggtitle("Parameter Coefficients: Just Mesocircuit MRI") +
  coord_flip()


# this one's weights aren't TERRIBLE looking
ggplot(filter(melt(Just_Meso_Dch$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(Just_Meso_Dch$weights)$value), max(melt(Just_Meso_Dch$weights)$value), by=50), 1)) +
  ggtitle("Parameter Coefficients: Just Mesocircuit Dch") +
  coord_flip()


# this look fine, of course.
ggplot(filter(melt(Just_Vols_MRI$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(Just_Vols_MRI$weights)$value), max(melt(Just_Vols_MRI$weights)$value), by=1), 1)) +
  ggtitle("Parameter Coefficients: Just Vols MRI") +
  coord_flip()


# same here; these are fine
ggplot(filter(melt(Just_Vols_Dch$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(Just_Vols_Dch$weights)$value), max(melt(Just_Vols_Dch$weights)$value), by=1), 1)) +
  ggtitle("Parameter Coefficients: Just Vols Dch") +
  coord_flip()


# pretty "bad"; really, the weights are just small. (except for outliers.)
ggplot(filter(melt(Just_Sedation_MRI$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(Just_Sedation_MRI$weights)$value), max(melt(Just_Sedation_MRI$weights)$value), by=1), 1)) +
  ggtitle("Parameter Coefficients: Just Sedation MRI") +
  coord_flip()



ggplot(filter(melt(Just_Sedation_Dch$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(Just_Sedation_Dch$weights)$value), max(melt(Just_Sedation_Dch$weights)$value), by=1), 1)) +
  ggtitle("Parameter Coefficients: Just Sedation Dch") +
  coord_flip()


# not bad, just small
ggplot(filter(melt(Just_Labs_MRI$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(Just_Labs_MRI$weights)$value), max(melt(Just_Labs_MRI$weights)$value), by=1), 1)) +
  ggtitle("Parameter Coefficients: Just Labs MRI") +
  coord_flip()




### For Just Meso @ Dch
ggplot(filter(melt(Just_Meso_Dch$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(Just_Meso_Dch$weights)$value), max(melt(Just_Meso_Dch$weights)$value), by=1), 1)) +
  ggtitle("Parameter Coefficients: Just Meso Dch") +
  coord_flip()




### BSD
# order goes ICH, Edema ipsi  then ICH, Edema contra
param_order <- c("MB_peduncle_ICH_ipsi", "MB_peduncle_edema_ipsi", "MB_peduncle_ICH_contro", "MB_peduncle_edema_contro",
                 "Teg_ICH_C", "Teg_edema_C",
                 "MB_ICH_C", "MB_edema_C",
                 "TH_ICH_ipsi", "TH_edema_ipsi", "TH_ICH_contro", "TH_edema_contro",
                 "Hypo_ICH_C", "Hypo_edema_C",
                 "BF_ICH_C", "BF_edema_C",
                 "Caudate_ICH_ipsi", "Caudate_edema_ipsi", "Caudate_ICH_contro", "Caudate_edema_contro",
                 "PUT_ICH_ipsi", "PUT_edema_ipsi", "PUT_ICH_contro", "PUT_edema_contro",
                 "GP_ICH_ipsi", "GP_edema_ipsi", "GP_ICH_contro", "GP_edema_contro",
                 "Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm",
                 "IVH")

# also after grep-renaming contro to contra
new_param_order <- c("Caudate_ICH_ipsi", "Caudate_edema_ipsi", "Caudate_ICH_contra", "Caudate_edema_contra",
                     "PUT_ICH_ipsi", "PUT_edema_ipsi", "PUT_ICH_contra", "PUT_edema_contra",
                     "GP_ICH_ipsi", "GP_edema_ipsi", "GP_ICH_contra", "GP_edema_contra",
                     "TH_ICH_ipsi", "TH_edema_ipsi", "TH_ICH_contra", "TH_edema_contra",
                     "BF_ICH_C", "BF_edema_C",
                     "Hypo_ICH_C", "Hypo_edema_C",
                     "MB_peduncle_ICH_ipsi", "MB_peduncle_edema_ipsi", "MB_peduncle_ICH_contra", "MB_peduncle_edema_contra",
                     "MB_ICH_C", "MB_edema_C",
                     "Teg_ICH_C", "Teg_edema_C",
                     "Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm",
                     "IVH")

## even newer order for the "All_Params" models
newer_param_order <- c("FCx_ICH_ipsi","FCx_edema_ipsi","FCx_ICH_contra","FCx_edema_contra",
                     "PCx_ICH_ipsi","PCx_edema_ipsi","PCx_ICH_contra","PCx_edema_contra",
                     "TCx_ICH_ipsi","TCx_edema_ipsi","TCx_ICH_contra","TCx_edema_contra",
                     "OCx_ICH_ipsi","OCx_edema_ipsi","OCx_ICH_contra","OCx_edema_contra",
                     "INS_ICH_ipsi","INS_edema_ipsi","INS_ICH_contra","INS_edema_contra",
                     "IC_ant_ICH_ipsi","IC_ant_edema_ipsi","IC_ant_ICH_contra","IC_ant_edema_contra",
                     "IC_post_ICH_ipsi","IC_post_edema_ipsi","IC_post_ICH_contra","IC_post_edema_contra",
                     
                     "Caudate_ICH_ipsi", "Caudate_edema_ipsi", "Caudate_ICH_contra", "Caudate_edema_contra",
                     "PUT_ICH_ipsi", "PUT_edema_ipsi", "PUT_ICH_contra", "PUT_edema_contra",
                     "GP_ICH_ipsi", "GP_edema_ipsi", "GP_ICH_contra", "GP_edema_contra",
                     "TH_ICH_ipsi", "TH_edema_ipsi", "TH_ICH_contra", "TH_edema_contra",
                     "BF_ICH_C", "BF_edema_C",
                     "Hypo_ICH_C", "Hypo_edema_C",
                     "MB_peduncle_ICH_ipsi", "MB_peduncle_edema_ipsi", "MB_peduncle_ICH_contra", "MB_peduncle_edema_contra",
                     "MB_ICH_C", "MB_edema_C",
                     "Teg_ICH_C", "Teg_edema_C",
                     
                     "Cereb_ICH_ipsi", "Cereb_edema_ipsi", "Cereb_ICH_contra", "Cereb_edema_contra",
                     "Vermis_ICH", "Vermis_edema",
                     
                     "Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm",
                     "IVH",
                     
                     "Propofol_at_MRI_MeanNorm","Propofol_before_MRI_MeanNorm",
                     "Midazolam_at_MRI_MeanNorm","Midazolam_before_MRI_MeanNorm",
                     "Dexmedetomidine_at_MRI_MeanNorm","Dexmedetomidine_before_MRI_MeanNorm",
                     "Fentanyl_at_MRI_MeanNorm","Fentanyl_before_MRI_MeanNorm",
                     "Gluc_MeanNorm","BUN_MeanNorm")






newer_param_order_dch <- c("FCx_ICH_ipsi","FCx_edema_ipsi","FCx_ICH_contra","FCx_edema_contra",
                       "PCx_ICH_ipsi","PCx_edema_ipsi","PCx_ICH_contra","PCx_edema_contra",
                       "TCx_ICH_ipsi","TCx_edema_ipsi","TCx_ICH_contra","TCx_edema_contra",
                       "OCx_ICH_ipsi","OCx_edema_ipsi","OCx_ICH_contra","OCx_edema_contra",
                       "INS_ICH_ipsi","INS_edema_ipsi","INS_ICH_contra","INS_edema_contra",
                       "IC_ant_ICH_ipsi","IC_ant_edema_ipsi","IC_ant_ICH_contra","IC_ant_edema_contra",
                       "IC_post_ICH_ipsi","IC_post_edema_ipsi","IC_post_ICH_contra","IC_post_edema_contra",
                       
                       "Caudate_ICH_ipsi", "Caudate_edema_ipsi", "Caudate_ICH_contra", "Caudate_edema_contra",
                       "PUT_ICH_ipsi", "PUT_edema_ipsi", "PUT_ICH_contra", "PUT_edema_contra",
                       "GP_ICH_ipsi", "GP_edema_ipsi", "GP_ICH_contra", "GP_edema_contra",
                       "TH_ICH_ipsi", "TH_edema_ipsi", "TH_ICH_contra", "TH_edema_contra",
                       "BF_ICH_C", "BF_edema_C",
                       "Hypo_ICH_C", "Hypo_edema_C",
                       "MB_peduncle_ICH_ipsi", "MB_peduncle_edema_ipsi", "MB_peduncle_ICH_contra", "MB_peduncle_edema_contra",
                       "MB_ICH_C", "MB_edema_C",
                       "Teg_ICH_C", "Teg_edema_C",
                       
                       "Cereb_ICH_ipsi", "Cereb_edema_ipsi", "Cereb_ICH_contra", "Cereb_edema_contra",
                       "Vermis_ICH", "Vermis_edema",
                       
                       "Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm",
                       "IVH",

                       "Gluc_MeanNorm","BUN_MeanNorm")





# Dch
BSD_Vol_Dch_newnames <- BSD_Vol_Dch$weights
names(BSD_Vol_Dch_newnames) <- gsub("contro", "contra", names(BSD_Vol_Dch_newnames))


ggplot(filter(melt(BSD_Vol_Dch_newnames), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot(outlier.shape=NA) + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  # axis scales used to be: limits=c(-7.5, 1.5), breaks=seq(-7.5, 1.5, by=0.3)
  # ok, never mind, this scale of limits=c(-8, 3), breaks=seq(-10, 10, by=0.2) messes it up for some of the MRI ones.
  scale_y_continuous(limits=c(-10, 10), breaks=seq(-10, 10, by=0.4)) +  #breaks=round(seq(min(melt(BSD_Vol_Dch$weights)$value), max(melt(BSD_Vol_Dch$weights)$value), by=1), 1)) +
  ggtitle("Brainstem + Diencephalon + Volumes + MLS Model", subtitle="Parameter Coefficients at Discharge") +
  ylab("Coefficient") +
  xlab("Parameter / Brain Region") +
  scale_x_discrete(limits=rev(new_param_order)) +
  coord_flip()




# test why some of the boxplot quantiles are seemingly plotted wrong
# ggplot(filter(melt(BSD_Vol_Dch$weights), variable == "MB_peduncle_edema_contro"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
#   geom_boxplot(outlier.shape=NA) + theme_bw() +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
#   scale_y_continuous(limits=c(0, 10), breaks=seq(0, 10, by=0.1)) +  #breaks=round(seq(min(melt(BSD_Vol_Dch$weights)$value), max(melt(BSD_Vol_Dch$weights)$value), by=1), 1)) +
#   ggtitle("Brainstem + Diencephalon + Volumes + MLS Model", subtitle="Parameter Coefficients at Discharge") +
#   ylab("Coefficient") +
#   xlab("Parameter / Brain Region") +
#   #scale_x_discrete(limits=rev(param_order)) +
#   coord_flip()
# 
# boxplot(filter(melt(BSD_Vol_Dch$weights), variable == "MB_peduncle_edema_contro")$value)
# 
# why <- filter(melt(BSD_Vol_Dch$weights), variable == "MB_peduncle_edema_contro")$value
# 
# qplot(x="a", y=why, geom="boxplot")
# ggplot(filter(melt(BSD_Vol_Dch$weights), variable == "MB_peduncle_edema_contro"), aes(x=variable, y=value)) + geom_boxplot() +
#   # it does appear to be due to the limits and breaks. For example, using limits=c(0,10) and breaks=seq(0, 10, by=0.1) shows the correct ~0.8 3rd quartile
#   scale_y_continuous(limits=c(-10, 10), breaks=seq(-10, 10, by=0.2))
# 
# # what could a good scale be?
# ggplot(filter(melt(BSD_Vol_Dch$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
#   geom_boxplot(outlier.shape=NA) + theme_bw() +
#   theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
#   # setting the breaks from -10 to 10 by 0.2 and then just limiting the axis from -8 to 3 *seems* to work
#   scale_y_continuous(limits=c(-8, 3), breaks=seq(-10, 10, by=0.2)) +  #breaks=round(seq(min(melt(BSD_Vol_Dch$weights)$value), max(melt(BSD_Vol_Dch$weights)$value), by=1), 1)) +
#   ggtitle("Brainstem + Diencephalon + Volumes + MLS Model", subtitle="Parameter Coefficients at Discharge") +
#   ylab("Coefficient") +
#   xlab("Parameter / Brain Region") +
#   scale_x_discrete(limits=rev(param_order)) +
#   coord_flip()




ggplot(filter(melt(BSD_Dch$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot(outlier.shape=NA) + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  # axis scales used to be: limits=c(-2.5, 1.5), breaks=seq(-2.5, 1.5, by=0.1)
  scale_y_continuous(limits=c(-10, 10), breaks=seq(-10, 10, by=0.4)) +  #breaks=round(seq(-6.4, 1.6, by=1), 1)) +
  ggtitle("Brainstem + Diencephalon Model", subtitle="Parameter Coefficients at Discharge") +
  ylab("Coefficient") +
  xlab("Parameter / Brain Region") +
  scale_x_discrete(limits=rev(param_order)) +
  coord_flip()



# MRI
ggplot(filter(melt(BSD_Vol_MRI$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot(outlier.shape=NA) + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  # axis scales used to be: limits=c(-7.5, 1.5), breaks=seq(-7.5, 1.5, by=0.3)
  scale_y_continuous(limits=c(-10, 10), breaks=seq(-10, 10, by=0.4)) +  #breaks=round(seq(min(melt(BSD_Vol_Dch$weights)$value), max(melt(BSD_Vol_Dch$weights)$value), by=1), 1)) +
  ggtitle("Brainstem + Diencephalon + Volumes + MLS Model", subtitle="Parameter Coefficients at MRI") +
  ylab("Coefficient") +
  xlab("Parameter / Brain Region") +
  scale_x_discrete(limits=rev(param_order)) +
  coord_flip()


ggplot(filter(melt(BSD_MRI$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot(outlier.shape=NA) + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  # axis scales used to be: limits=c(-2.5, 1.5), breaks=seq(-2.5, 1.5, by=0.1)
  scale_y_continuous(limits=c(-10, 10), breaks=seq(-10, 10, by=0.4)) +  #breaks=round(seq(-6.4, 1.6, by=1), 1)) +
  ggtitle("Brainstem + Diencephalon Model", subtitle="Parameter Coefficients at MRI") +
  ylab("Coefficient") +
  xlab("Parameter / Brain Region") +
  scale_x_discrete(limits=rev(param_order)) +
  coord_flip()



# just volumes
ggplot(filter(melt(Just_Vols_MRI$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) +
  geom_boxplot(outlier.shape=NA) + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(limits=c(-10, 10), breaks=seq(-10, 10, by=0.4)) +  #breaks=round(seq(-6.4, 1.6, by=1), 1)) +
  ggtitle("Volumes", subtitle="Parameter Coefficients at MRI") +
  ylab("Coefficient") +
  xlab("Parameter / Brain Region") +
  coord_flip()


ggplot(filter(melt(Just_Vols_Dch$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) +
  geom_boxplot(outlier.shape=NA) + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(limits=c(-10, 10), breaks=seq(-10, 10, by=0.4)) +  #breaks=round(seq(-6.4, 1.6, by=1), 1)) +
  ggtitle("Volumes", subtitle="Parameter Coefficients at Discharge") +
  ylab("Coefficient") +
  xlab("Parameter / Brain Region") +
  coord_flip()




#### Weight boxplots for the two final modesl of BSD_Vol_IVH ####
# BSD_Vol_IVH_Dch$AUC,
# +                       BSD_Vol_IVH_MRI$AUC,
# +                       Vols_IVH_MRI$AUC,
# +                       Vols_IVH_Dch$AUC)
BSD_Vol_IVH_Dch_newnames <- BSD_Vol_IVH_Dch$weights
names(BSD_Vol_IVH_Dch_newnames) <- gsub("contro", "contra", names(BSD_Vol_IVH_Dch_newnames))



ggplot(filter(melt(BSD_Vol_IVH_Dch_newnames), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot(outlier.shape=NA) + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  # axis scales used to be: limits=c(-2.5, 1.5), breaks=seq(-2.5, 1.5, by=0.1)
  scale_y_continuous(limits=c(-10, 10), breaks=seq(-10, 10, by=0.4)) +  #breaks=round(seq(-6.4, 1.6, by=1), 1)) +
  ggtitle("Brainstem + Diencephalon + Volumes + MLS + IVH Model", subtitle="Parameter Coefficients at Discharge") +
  ylab("Coefficient") +
  xlab("Parameter / Brain Region") +
  scale_x_discrete(limits=rev(new_param_order)) +
  coord_flip()




BSD_Vol_IVH_MRI_newnames <- BSD_Vol_IVH_MRI$weights
names(BSD_Vol_IVH_MRI_newnames) <- gsub("contro", "contra", names(BSD_Vol_IVH_MRI_newnames))


ggplot(filter(melt(BSD_Vol_IVH_MRI_newnames), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot(outlier.shape=NA) + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  # axis scales used to be: limits=c(-2.5, 1.5), breaks=seq(-2.5, 1.5, by=0.1)
  scale_y_continuous(limits=c(-10, 10), breaks=seq(-10, 10, by=0.4)) +  #breaks=round(seq(-6.4, 1.6, by=1), 1)) +
  ggtitle("Brainstem + Diencephalon + Volumes + MLS + IVH Model", subtitle="Parameter Coefficients at MRI") +
  ylab("Coefficient") +
  xlab("Parameter / Brain Region") +
  scale_x_discrete(limits=rev(new_param_order)) +
  coord_flip()



#### Also redo the weight boxplots for all params ####

All_params_Dch_newnames <- All_params_Dch$weights
names(All_params_Dch_newnames) <- gsub("contro", "contra", names(All_params_Dch_newnames))

ggplot(filter(melt(All_params_Dch_newnames), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot(outlier.shape=NA) + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  # axis scales used to be: limits=c(-2.5, 1.5), breaks=seq(-2.5, 1.5, by=0.1)
  scale_y_continuous(limits=c(-10, 10), breaks=seq(-10, 10, by=0.4)) +  #breaks=round(seq(-6.4, 1.6, by=1), 1)) +
  ggtitle("All Parameters Model", subtitle="Parameter Coefficients at Discharge") +
  ylab("Coefficient") +
  xlab("Parameter / Brain Region") +
  scale_x_discrete(limits=rev(newer_param_order)) +
  coord_flip()




All_params_MRI_newnames <- All_params_MRI$weights
names(All_params_MRI_newnames) <- gsub("contro", "contra", names(All_params_MRI_newnames))

ggplot(filter(melt(All_params_MRI_newnames), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_boxplot(outlier.shape=NA) + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  # axis scales used to be: limits=c(-2.5, 1.5), breaks=seq(-2.5, 1.5, by=0.1)
  scale_y_continuous(limits=c(-10, 10), breaks=seq(-10, 10, by=0.4)) +  #breaks=round(seq(-6.4, 1.6, by=1), 1)) +
  ggtitle("All Parameters Model", subtitle="Parameter Coefficients at MRI") +
  ylab("Coefficient") +
  xlab("Parameter / Brain Region") +
  scale_x_discrete(limits=rev(newer_param_order)) +
  coord_flip()















ggplot(filter(melt(Vols_IVH_MRI$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) +
  geom_boxplot(outlier.shape=NA) + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(limits=c(-10, 10), breaks=seq(-10, 10, by=0.4)) +  #breaks=round(seq(-6.4, 1.6, by=1), 1)) +
  ggtitle("Volumes + IVH", subtitle="Parameter Coefficients at MRI") +
  ylab("Coefficient") +
  xlab("Parameter / Brain Region") +
  scale_x_discrete(limits=rev(c("Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm", "IVH"))) +
  coord_flip()


ggplot(filter(melt(Vols_IVH_Dch$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) +
  geom_boxplot(outlier.shape=NA) + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(limits=c(-10, 10), breaks=seq(-10, 10, by=0.4)) +  #breaks=round(seq(-6.4, 1.6, by=1), 1)) +
  ggtitle("Volumes + IVH", subtitle="Parameter Coefficients at Discharge") +
  ylab("Coefficient") +
  xlab("Parameter / Brain Region") +
  scale_x_discrete(limits=rev(c("Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm", "IVH"))) +
  coord_flip()






# as violin plots -- they look awful
ggplot(filter(melt(BSD_Vol_Dch$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_violin() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(breaks=round(seq(min(melt(BSD_Vol_Dch$weights)$value), max(melt(BSD_Vol_Dch$weights)$value), by=1), 1)) +
  ggtitle("Brainstem + Diencephalon + Volumes + MLS Model", subtitle="Parameter Coefficients") +
  ylab("Coefficient") +
  xlab("Parameter / Brain Region") +
  coord_flip()


ggplot(filter(melt(BSD_Dch$weights), variable != "(Intercept)"), aes(x=reorder(variable, -value, FUN=median), y=value)) + 
  geom_violin() + theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.25)) +
  scale_y_continuous(limits=c(-2.5, 1.5), breaks=seq(-2.5, 1.5, by=0.1)) +  #breaks=round(seq(-6.4, 1.6, by=1), 1)) +
  ggtitle("Brainstem + Diencephalon Model", subtitle="Parameter Coefficients") +
  ylab("Coefficient") +
  xlab("Parameter / Brain Region") +
  coord_flip()


#### Can use scale_x_discrete(limits=c(list, names, in, order, here)) to reorder the parameter names



#### Table of weights ####

# BSD_Vol_Dch_Weight_Table <- melt(BSD_Vol_Dch$weights) %>%
#   filter(variable != "(Intercept)") %>%
#   group_by(variable) %>%
#   summarise(median_weight = median(value),
#             qt_25 = quantile(value)[2],
#             qt_75 = quantile(value)[4]) %>%
#   arrange(median_weight, qt_25, qt_75)





#### Reformat/force table into the format Ben wrote ####


BSD_Vol_Dch_Weight_Table <- ReformatWeightTable(BSD_Vol_Dch$weights, toReturn="unformatted")

# with volumes
BSD_Vol_Dch_Weight_Table_formatted <- ReformatWeightTable(BSD_Vol_Dch$weights, toReturn="formatted")
BSD_Vol_MRI_Weight_Table_formatted <- ReformatWeightTable(BSD_Vol_MRI$weights, toReturn="formatted")

# without volumes
BSD_Dch_Weight_Table_formatted <- ReformatWeightTable(BSD_Dch$weights, toReturn="formatted")
BSD_MRI_Weight_Table_formatted <- ReformatWeightTable(BSD_MRI$weights, toReturn="formatted")

# just volumes
Just_Vols_MRI_Weight_Table <- ReformatWeightTable(Just_Vols_MRI$weights, toReturn="unformatted")
Just_Vols_Dch_Weight_Table <- ReformatWeightTable(Just_Vols_Dch$weights, toReturn="unformatted")


# For BSD Vol IVH
BSD_Vol_IVH_Dch_weight_Table_formatted <- ReformatWeightTable(BSD_Vol_IVH_Dch$weights, toReturn="formatted")
BSD_Vol_IVH_MRI_weight_Table_formatted <- ReformatWeightTable(BSD_Vol_IVH_MRI$weights, toReturn="formatted")

# name_order <- c("Caudate_contro", "Caudate_ipsi", "PUT_contro", "PUT_ipsi", "GP_contro", "GP_ipsi",
#                 "TH_contro", "TH_ipsi", "BF_C", "Hypo_C", "MB_peduncle_contro", "MB_peduncle_ipsi",
#                 "MB_C", "Teg_C")


#write.xlsx(BSD_Vol_Dch_Weight_Table, "E:/Experiments/ICH_MRI/BDS Boxplots/BSD_Boxplot_Table.xlsx")

write.xlsx(BSD_Vol_Dch_Weight_Table_formatted, "E:/Experiments/ICH_MRI/Weight Tables/BSD_Vol_Dch_weights.xlsx", rowNames=T)
write.xlsx(BSD_Vol_MRI_Weight_Table_formatted, "E:/Experiments/ICH_MRI/Weight Tables/BSD_Vol_MRI_weights.xlsx", rowNames=T)
write.xlsx(BSD_Dch_Weight_Table_formatted, "E:/Experiments/ICH_MRI/Weight Tables/BSD_Dch_weights.xlsx", rowNames=T)
write.xlsx(BSD_MRI_Weight_Table_formatted, "E:/Experiments/ICH_MRI/Weight Tables/BSD_MRI_weights.xlsx", rowNames=T)

write.xlsx(BSD_Vol_IVH_Dch_weight_Table_formatted, "E:/Experiments/ICH_MRI/Weight Tables/BSD_Vol_IVH_Dch_weights.xlsx", rowNames=T)
write.xlsx(BSD_Vol_IVH_MRI_weight_Table_formatted, "E:/Experiments/ICH_MRI/Weight Tables/BSD_Vol_IVH_MRI_weights.xlsx", rowNames=T)



## This is all done in the function now
# ROIs <- unique(BSD_Vol_Dch_Weight_Table$variable)
# ROIs <- unique(gsub("_edema", "", gsub("_ICH", "", ROIs)))
# newtab <- matrix(nrow=length(ROIs), ncol=2)
# rownames(newtab) <- ROIs
# colnames(newtab) <- c("ICH", "Edema")
# # try to assign values
# for (i in 1:nrow(BSD_Vol_Dch_Weight_Table)) {
#   cur_name <- as.character(BSD_Vol_Dch_Weight_Table$variable[i])
#   if (!cur_name %in% c("Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm")) {
#     # make name match that from newtab
#     name_wo_type <- gsub("_edema", "", gsub("_ICH", "", cur_name))
#     # find the row to add the value to
#     newtab_ix <- which(rownames(newtab) == name_wo_type)
#     # now find out if this was ich or edema
#     if (grepl("edema", cur_name)) {
#       newtab_col <- "Edema"
#     } else if (grepl("ICH", cur_name)) {
#       newtab_col <- "ICH"
#     }
#     
#     # now paste the median weight and CI together
#     val_to_add <- paste(round(BSD_Vol_Dch_Weight_Table$median_weight[i], 2), "[", 
#                         round(BSD_Vol_Dch_Weight_Table$qt_25[i], 2), " ", round(BSD_Vol_Dch_Weight_Table$qt_75[i], 2), "]", sep="")
#     
#     newtab[newtab_ix, newtab_col] <- val_to_add
#     
#     #nm_check <- unlist(strsplit(cur_name, "_"))
#     
#   } else {
#     # just put the volume ones in the ICH column.
#     newtab_ix <- which(rownames(newtab) == cur_name)
#     
#     val_to_add <- paste(round(BSD_Vol_Dch_Weight_Table$median_weight[i], 2), "[", 
#                         round(BSD_Vol_Dch_Weight_Table$qt_25[i], 2), " ", round(BSD_Vol_Dch_Weight_Table$qt_75[i], 2), "]", sep="")
#     
#     newtab[newtab_ix, "ICH"] <- val_to_add
#     
#   }
# }








#### Stacked bar plot of each region for conscious and unconscious patients for TH patients ####
barplot_dat <- melt(loc_analyse.ipsi_TH[c(which(names(loc_analyse.raw3) %in% c("MRN", "follow2")), all_ix3[-Extra_ix_to_remove])], id.vars=c("MRN", "follow2"))#, "follow", "MRI_Cs", "Discharge_date", "DEATH_DATE", "MRI_date_Loc", "mri_date", "infra_tent", 
                                                #"supra_tent", "supra_tent2", "infra_tent2", "follow2", "follow3", "MRI_Cs2", "MRI_Cs3"))

bad_params <- c("Old_stroke", "Old_ICH", "Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm")

barplot_dat2 <- barplot_dat %>%
  filter(! variable %in% bad_params )

# only where value==1 (this will remove any of the continuous measures, but that's fine for now I guess)
barplot_dat3 <- barplot_dat2 %>%
  filter(value==1)

# make a summarized count table first?
barplot_counts <- barplot_dat3 %>%
  group_by(variable, value, follow2) %>%
  summarise(count = n())

# summarized table to take averages (i.e., percents) instead
barplot_percents <- barplot_dat2 %>%
  group_by(variable, follow2) %>%
  summarise(mean = mean(as.numeric(value)))

barplot_counts$follow2 <- as.factor(barplot_counts$follow2)
barplot_percents$follow2 <- as.factor(barplot_percents$follow2)

ggplot(barplot_counts, aes(x=reorder(variable, count), y=count, fill=follow2)) + 
  geom_bar(stat="identity") + #geom_histogram(stat="count")
  coord_flip()


ggplot(barplot_percents, aes(x=reorder(variable, mean), y=mean, fill=follow2)) + 
  geom_bar(stat="identity") +  #, position="dodge"
  ylab("percentage") + xlab("") +
  coord_flip()













