library(openxlsx)
library(tidyverse)
library(caret)
library(glmnet)
library(reshape2)
library(pROC)
library(ROCR)
library(lubridate)
library(cvAUC)


# Normalized hemorrhage volume
# removes outlier w/ small brain b/c of incomplete scan
# Will replace this value with the mean value
normalizeByBrainVol <- function(x, brain_vol_data) {
  # replace outlier smallest brain volumn with the mean
  # Warning: this doesn't perminantly change that outlier value to the mean unlike before.
  # not sure if it matters though, since I don't think we directly use the Brain_vol param for anything
  minval <- which(brain_vol_data == min(brain_vol_data))
  mBV <- mean(brain_vol_data[-minval])
  brain_vol_data[minval] <- mBV
  
  res <- (x * brain_vol_data) / mBV
}


meanNormalization <- function(x) {
  res <- (x - mean(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
}


calcConfInt <- function(data) {
  m    <- mean(data)
  stdv <- sd(data)
  n    <- length(data)
  
  err  <- qt(0.975, df = n-1) * stdv / sqrt(n)
  CI.l <- m - err
  CI.u <- m + err
  
  return(c(lower = CI.l, upper = CI.u))
}


MakePlot <- function(data, group, param, path = NULL, save = FALSE) {
  
  dat <- select_(data, group, param)
  dat2 <- melt(dat, id = group)
  dat2[, group] <- as.factor(dat2[, group])
  
  plt <- ggplot(data = dat2, aes_string(x = group, y = "value", fill = "variable")) +
    geom_boxplot(width = 0.2, outlier.shape = NA) + ggtitle(param) + theme_minimal() +
    geom_jitter(position = position_jitter(width = 0.05, height = 0.1), size = 1)
  
  if (save == TRUE) {
    ggsave(filename = paste(path, "/", param, "_", group, ".pdf", sep = ""),
           plot = plt, device = "pdf")
  }
  return(plt)
}

RunModel <- function(data, indep_var_ix, dep_var_vals, opposite_vals = NA, n_runs, lambda_grid = NULL,
                     p = 0.80, alpha = 0.95, nfolds = 10, lambda_use = "lambda.min") {
  
  AUCs <- c()
  weights_table <- data.frame()
  probs_table <- as.data.frame(matrix(nrow = nrow(data)))
  rownames(probs_table) <- rownames(data)
  lambda_mins <- c()
  for (i in 1:n_runs) {
    # for some reason, this defaulted to rownames.force = T before???? It NEEDS to be true, 
    # or else the row indices for each patient won't carry over in order to make the confusion matrix
    df <- data.matrix(data[, indep_var_ix], rownames.force = T)
    # add dependent variable to beginning of table
    df <- cbind(as.factor(dep_var_vals), df)
    
    # create testing/training sets
    train_ix <- createDataPartition(as.factor(df[, 1]), p = p, list = F)
    train <- df[train_ix, ]
    test <- df[-train_ix, ]
    
    # run model
    modfit <- cv.glmnet(x = train[,2:ncol(train)], y = train[,1], lambda = lambda_grid,
                        family = "binomial", type.measure = "auc", alpha = alpha, nfolds = nfolds)
    
    # # then use the lambda min in the regular glmnet function
    # modfit <- glmnet(x = train[,2:ncol(train)], y = train[,1], lambda = modfit$lambda.min,
    #                  family = binomial", alpha = alpha)
    lambda_mins <- c(lambda_mins, modfit$lambda.min)
    
    # get coefficients/weights
    weights <- t(coef(modfit, s = lambda_use)[, 1])
    weights_table <- rbind(weights_table, weights)
    
    # predict using the test data
    pred.prob <- predict(modfit, newx = data.matrix(test[, 2:ncol(test)]),
                         type = "response", s = lambda_use)
    
    # save the probabilities for the patients in the current test dataset
    # using the current row indices of the test table and matching to the
    # row indices in the main table
    probs_table <- cbind(probs_table, X = pred.prob[match(rownames(probs_table), rownames(pred.prob))])
    
    # predict new values
    pred <- prediction(pred.prob, test[, 1])
    # get AUC
    auc <- performance(pred, "auc")
    # save aucs to vector
    # 'auc' is an S4 class object; need to extract using '@'
    AUCs <- c(AUCs, auc@y.values[[1]])
  }
  
  # set first column as row indices for probs_table and reset column names
  colnames(probs_table) <- as.character(1:length(colnames(probs_table)))
  probs_table[, 1] <- as.numeric(rownames(probs_table))
  colnames(probs_table)[1] <- "patient"
  
  # now melt table and calculate confusion matrix, assigning pred.values < 0.5 as unconscious
  medProbs <- melt(probs_table, id = "patient") %>%
    group_by(patient) %>%
    summarize(medProb = median(value, na.rm = T))
  
  medProbs$class_vals <- ifelse(medProbs$medProb < 0.5, yes = 0, no = 1)
  # order by patient row index, that way the dep_var_vals match up with the same patient
  # I think this step is unnecessary now that the row-index column is numeric, but whatever
  medProbs <- medProbs[order(medProbs$patient), ]
  medProbs$actual_vals <- dep_var_vals
  medProbs$opposite_vals <- opposite_vals
  cm <- table(actual = medProbs$actual_vals, predict = medProbs$class_vals)
  if (!is.na(opposite_vals)) { # should probably use an any() function on this check to stop the warning message
    cm2 <- table(opposite_meas = medProbs$opposite_vals, predict = medProbs$class_vals)
  }
  
  # return weights table and AUCs together in a list
  if (any(!is.na(opposite_vals))) { # should probably use an any() function on this check to stop the warning message
    out <- list(AUC = AUCs, lambda_min = lambda_mins, weights = weights_table, probs = probs_table, medProb = medProbs, 
                confus_mat = cm, oppo_confus_mat = cm2)
  } else {
    out <- list(AUC = AUCs, lambda_min = lambda_mins, weights = weights_table, probs = probs_table, medProb = medProbs, 
                confus_mat = cm)
  }
  
  return(out)
}


ReformatWeightTable <- function(weights, toReturn = "unformatted") {
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
    newtab <- matrix(nrow = length(ROIs), ncol = 2)
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
                            round(weight_table$qt_25[i], 2), " ", round(weight_table$qt_75[i], 2), "]", sep = "")
        
        newtab[newtab_ix, newtab_col] <- val_to_add
        
      } else {
        # just put the volume ones in the ICH column.
        newtab_ix <- which(rownames(newtab) == cur_name)
        
        val_to_add <- paste(round(weight_table$median_weight[i], 2), "[", 
                            round(weight_table$qt_25[i], 2), " ", round(weight_table$qt_75[i], 2), "]", sep = "")
        
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


makeGOSEplot <- function(medProbs, var1, var2) {
  
  fp <- which(medProbs[, var1] == 1 & medProbs[, var2] == 0)
  tn <- which(medProbs[, var1] == 0 & medProbs[, var2] == 0)
  tp <- which(medProbs[, var1] == 1 & medProbs[, var2] == 1)
  fn <- which(medProbs[, var1] == 0 & medProbs[, var2] == 1)
  
  fpTab <- loc_analyse.raw3[fp, c("MRN", "GOSE_All"), drop = F]
  fpTab$class_grp <- "false_pos"
  tnTab <- loc_analyse.raw3[tn, c("MRN", "GOSE_All"), drop = F]
  tnTab$class_grp <- "true_neg"
  tpTab <- loc_analyse.raw3[tp, c("MRN", "GOSE_All"), drop = F]
  tpTab$class_grp <- "true_pos"
  fnTab <- loc_analyse.raw3[fn, c("MRN", "GOSE_All"), drop = F]
  fnTab$class_grp <- "false_neg"
  
  GOSE_class <- rbind(fpTab, tnTab, tpTab, fnTab)
  GOSE_class$GOSE_All <- factor(GOSE_class$GOSE_All, exclude = NULL)
  
  return(GOSE_class)
}


# this uses the above function to generate the GOSE contingency table, the fisher test on that table,
# and the plot of how many GOSE scores are in each category of the input model's "opposite" confusion matrix.
GOSEfisherAndPlot <- function(data, pred_val = "class_vals", actual_val = "opposite_vals") {
  # first use the function that makes the data and then
  # do the fisher test by creating the count table using the GOSE data
  pltDat <- makeGOSEplot(data[["medProb"]], pred_val, actual_val)
  
  summaryTable <- pltDat %>%
    filter(class_grp %in% c("false_pos", "true_neg")) %>%
    mutate(gose_grp = as.factor(ifelse(as.numeric(as.character(GOSE_All)) >= 4, 1, 0)),
           class_grp = as.factor(class_grp)) %>%
    filter(!is.na(gose_grp)) %>%
    group_by(class_grp, gose_grp)
  
  contingTable <- table(summaryTable$gose_grp, summaryTable$class_grp)
  fish_out <- fisher.test(contingTable)
  
  plt <- ggplot(pltDat, aes(x = class_grp, fill = GOSE_All)) + 
    geom_bar(stat = "count", position = "dodge") +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    xlab("Classification Group") +
    ggtitle(paste("GOSE Distribution per ", deparse(substitute(data)), 
                  " 'opposite' Confusion Matrix Classification", sep = ""), 
            subtitle = paste("true positive = conscious; true negative = unconscious; model = ", 
                             deparse(substitute(data)), sep = "")) +
    theme_bw() +
    coord_flip()
  
  all_dat <- list(contingTable = contingTable,
                  fisher = fish_out,
                  the_plot = plt)
  
  return(all_dat)
}


# functions to use in summarise_all to get the 25th and 75th quantiles seperately, and as a single value
get_25th <- function(x, na.rm = T) {
  out <- quantile(x, na.rm = T)[2]
}


get_75th <- function(x, na.rm = T) {
  out <- quantile(x, na.rm = T)[4]
}


# find accuracy, specificity, etc. for a confusion matrix
CMstats <- function(confus_mat) {
  tot <- sum(confus_mat)
  # percentages/rates of type I and II error
  tp <- confus_mat[2, 2]
  tn <- confus_mat[1, 1]
  fp <- confus_mat[1, 2]
  fn <- confus_mat[2, 1]
  act_yes <- sum(confus_mat[1, ])
  act_no <- sum(confus_mat[2, ])
  
  acc <- (tp + tn) / tot  # Accuracy
  miss <- (fp + fn) / tot  # Misclassification rate
  typeIerr <- fp / act_yes  # Type I error rate (predict yes when actually no)
  spec <- tn / act_yes  # Specificity (how often correctly predict no), ie, 1-Type I error
  typeIIerr <- fn / act_no  # Type II error rate (predict no when actually yes)
  sens <- tp / act_no  # Sensitivity (how often correctly predict yes), ie, 1-Type II error
  null_err <- act_no / tot  # null error rate (percent wrong if always predicted the majority class)
  cohen <- acc - null_err  # Cohen's Kappa (how well classifier performs compared to chance)
  
  outv <- c("Accuracy" = acc, "Misclass. rate" = miss, "Type I error" = typeIerr, "Specificity" = spec,
            "Type II error" = typeIIerr, "Sensitivity" = sens, "Null error" = null_err, "Cohen's kappa" = cohen)
  
  return(outv)
}


#---- Load and clean data ---
# main data, generated from ICH_MRI_v8.Rmd
loc_analyse.raw3 <- read.xlsx("../MRI_Merged_Data_158patients.xlsx")
# MRI eval note date table
MRI_time_dat <- read.xlsx("../Ben files/MRI_ICH/MRI Imaging list radiology comprehensive-br2.xlsx")
# convert dates
loc_analyse.raw3[, c("Discharge_date", "DEATH_DATE", "MRI_date_Loc")] <- lapply(
  loc_analyse.raw3[, c("Discharge_date", "DEATH_DATE", "MRI_date_Loc")],
  as.Date,
  origin = "1899-12-30"
)
# obnoxious way to get datetimes from excel number "datetimes"
MRI_time_dat$MRI.eval.note.date <- as.POSIXct(as.numeric(MRI_time_dat$MRI.eval.note.date) * (60*60*24),
                                              origin = "1899-12-30", tz = "GMT")
# join mri time data
loc_analyse.raw3 <- merge(loc_analyse.raw3, 
                          MRI_time_dat[, c("MRN", "MRI.eval.note.date")], 
                          by = "MRN")
# remove symbols from column names
names(loc_analyse.raw3) <- gsub("\\.", "_", names(loc_analyse.raw3))
names(loc_analyse.raw3) <- gsub("\\(|\\)", "", names(loc_analyse.raw3))
# change variable names and remove columns
loc_analyse.raw3 <- loc_analyse.raw3 %>%
  rename("MLS" = "MLS_[mm]") %>%
  select(-Uncal_herniation_to_which_side, -Meso_ICH_2_ipsi, -Meso_ICH_2_contro) %>%
  # change the additional 6 other patients with missing Edema volume to 0
  mutate(Ed_vol = ifelse(is.na(Ed_vol), yes = 0, no = Ed_vol))
  
# Normalize params by brain volume
loc_analyse.raw3[, c("Hg_vol_Norm", "Ed_vol_Norm", "MLS_Norm")] <- lapply(
  loc_analyse.raw3[, c("Hg_vol", "Ed_vol", "MLS")],
  normalizeByBrainVol,
  brain_vol_data = loc_analyse.raw3$Brain_vol
)

# and then normalize all by mean
loc_analyse.raw3[, c("Hg_vol_MeanNorm", "Ed_vol_MeanNorm", "MLS_MeanNorm")] <- lapply(
  loc_analyse.raw3[, c("Hg_vol_Norm", "Ed_vol_Norm", "MLS_Norm")], meanNormalization
)


#---- Check median length of stay for patients ----

# TODO never got around to finish "cleaning" the script LOOOL.
