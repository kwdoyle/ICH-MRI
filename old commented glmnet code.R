
# params to use:
# exclude <- NULL
# for (i in 1:length(names(loc_analyse.raw3))) {
#   if (all(loc_analyse.raw3[,i]==0, na.rm=T)) {
#     exclude <- c(exclude, i)
#   }
# }
# 
# edema_nms <- grep("edema", names(loc_analyse.raw3))
# ich_nms <- grep("ICH", names(loc_analyse.raw3))
# 
# xclude_nms <- names(loc_analyse.raw3)[exclude]
# 
# #names(loc_analyse.raw3)[-c(1,exclude,edema_nms,84:97)]
# colstouse <- names(loc_analyse.raw3)[c(ich_nms)]
# # remove cols with all 0
# colstouse <- colstouse[!colstouse %in% xclude_nms]
# 
# # other cols to use
# othercols <- c(3:81,96:104)


# all params without the the summarized ones we made
#othercols <- c(3:81,96:104)
# not using TH subgroups
#othercols <- c(3:65,80:88)
# not using edema
#othercols <- c(3:37, 52:58)









# change all to numeric


# # use select w/ dplyr instead.
# tst <- loc_analyse.raw3 %>%
#   dplyr::select(IVH, Teg_ICH_contro, )

# remove columns via grep + manual
# edema_cols <- grep("edema", names(loc_analyse.raw3))
# which(names(loc_analyse.raw3) == "sf_hern" | names(loc_analyse.raw3) == "Transtentorial.herniation" |
#       names(loc_analyse.raw3) == "Cerebellar.tonsillar.herniation" | names(loc_analyse.raw3) == "MLS" |
#       names(loc_analyse.raw3) == "Hg_vol" | names(loc_analyse.raw3) == "Ed_vol" | names(loc_analyse.raw3) == "Brain_vol" |
#         )









### New analysis 
# fitnew <- cv.glmnet(x=data.matrix(loc_analyse.raw3[,colstouse]), y=as.factor(loc_analyse.raw3$MRI_Cs2),
#                  family="binomial", type.measure="auc")
# coef(fitnew, s="lambda.min")
# plot(fitnew)


### New analysis 
# fitnew <- cv.glmnet(x=data.matrix(loc_analyse.raw3[,colstouse]), y=as.factor(loc_analyse.raw3$MRI_Cs2),
#                  family="binomial", type.measure="auc")
# coef(fitnew, s="lambda.min")
# plot(fitnew)







# rm_ix <- which(names(loc_analyse.raw3[,colstouse_Hgvol_MLS]) %in% TH_to_rm)
# colstouse_Hgvol_MLS[-rm_ix]








#names(loc_analyse.raw3)[colstouse_Hgvol_MLS][-TH_ix_to_remove][-Meso_ix_to_remove]

## Apparently this accepts family="multinomial" too, so we can do an analysis with the 3-level outcome?
## not sure if we want to define the 'weights' parameter.

# fit1 <- glmnet(x=data.matrix(loc_analyse.raw3[,colstouse]), y=as.factor(loc_analyse.raw3$MRI_Cs2), family="binomial")
# summary(fit1)
# fit1$beta
# fit1$lambda
# plot(fit1, xvar="dev", label=T)
# 
# ## This needs to accept a DATA.MATRIX instead of a normal matrix or as.matrix.
# predict(fit1, newx=data.matrix(loc_analyse.raw3[,colstouse])[1:5,])
# 
# print(fit1)
# 
# 
# 
# # cross-validation fit. This will choose the best model automatically,
# # as opposed to returning all possible model which is done above.
# cvfit <- cv.glmnet(x=data.matrix(loc_analyse.raw3[,colstouse]), y=as.factor(loc_analyse.raw3$MRI_Cs2), family="binomial")
# plot(cvfit)
# # value of lambda that gives the minimum mean corss-validated error
# # these are the coeficients for each "useful" param using this lambda
# coef(cvfit, s="lambda.min")
# 
# # predict with this? predicts the first 5 observations?
# predict(cvfit, newx=data.matrix(loc_analyse.raw3[,colstouse])[1:20,], type="class", s="lambda.min")
# 
# 
# # fit using misclassification error for the cross-validation
# cvfit2 <- cv.glmnet(x=data.matrix(loc_analyse.raw3[,colstouse]), y=as.factor(loc_analyse.raw3$MRI_Cs2),
#                     family="binomial", type.measure="class")  #

# find out which parameters give an error
# for (i in 2:length(othercols)) {
#   cvfit3 <- cv.glmnet(x=data.matrix(loc_analyse.raw3[,othercols[1:i]]), y=as.factor(loc_analyse.raw3$MRI_Cs2),
#                     family="binomial", type.measure="class")
# }
# all seem good now


#### NONE OF THIS MODEL SELECTION MATTERS
#### Can maybe average the results from n-repeated model runs, though, and see if those differ between the different models
# The results are always a little different each time the model is run
# --even if it's the same model being run twice.
# So the small differences between these models are mostly
# just due to chance.

# repeat multiple times to get an average output per model
# conf_matrices <- list()
# for (i in 1:500) {
#   fit <-  cv.glmnet(x=data.matrix(loc_analyse.raw3[,othercols]), y=as.factor(loc_analyse.raw3$MRI_Cs2),
#                     family="binomial", type.measure="class")
#   
#   predvals <- predict(fit, newx=data.matrix(loc_analyse.raw3[,othercols]), type="class", s="lambda.min")
#   actvals <- as.factor(loc_analyse.raw3$MRI_Cs2)
#   cm <- table(actvals, predvals)
#   conf_matrices[[i]] <- cm
# }
# 
# avg_cm <- round(apply(simplify2array(conf_matrices), c(1,2), mean))
# 

# for (i in 1:ncol(df)) {
#   if (any(is.na(df[,i]))==T) {
#     print(i)
#   }
# }











#names(loc_analyse.raw3)[colstouse_Hgvol_MLS][-TH_ix_to_remove][-Meso_ix_to_remove]


# is it the predict function that outputs differently?
# nope.
# for (i in 1:70) {
#   predict(fit, newx=data.matrix(loc_analyse.raw3[,othercols]), type="class", s="lambda.min")
#   actvals <- as.factor(loc_analyse.raw3$MRI_Cs2)
#   cm <- table(actvals, predvals)
#   print(cm)
# }

# 

### The stuff below is basically what I put into the above function
#  cvfit3 <- cv.glmnet(x=data.matrix(loc_analyse.raw3[,othercols]), y=as.factor(loc_analyse.raw3$MRI_Cs2),
#                     family="binomial", type.measure="class")
#  
# 
# 
# 
# 
# plot(cvfit2)
# coef(cvfit2, s="lambda.min")
# coef(cvfit3, s="lambda.min")
# 
# predict(cvfit2, newx=data.matrix(loc_analyse.raw3[,colstouse])[1:20,], type="class", s="lambda.min")
# 
# # make a confusion matrix
# #predvals <- predict(cvfit2, newx=data.matrix(loc_analyse.raw3[,colstouse]), type="class", s="lambda.min")
# predvals <- predict(cvfit3, newx=data.matrix(loc_analyse.raw3[,othercols]), type="class", s="lambda.min")
# actvals <- as.factor(loc_analyse.raw3$MRI_Cs2)
# cm <- table(actvals, predvals)
# tot <- sum(cm)
# # percentages/rates of type I and II error
# # nope, can't use the same denominator for all
# #round(table(actvals, predvals) / sum(table(actvals, predvals)) * 100, digits=1)
# tp <- cm[2,2]
# tn <- cm[1,1]
# fp <- cm[1,2]
# fn <- cm[2,1]
# act_yes <- sum(cm[1,])
# act_no <- sum(cm[2,])
# 
# ## Accuracy
# acc <- (tp + tn) / tot
# ## Misclassification rate
# miss <- (fp + fn) / tot
# ## Type I error rate (predict yes when actually no)
# typeIerr <- fp / act_yes
# # Specificity (how often correctly predict no), ie, 1-Type I error
# spec <- tn / act_yes
# ## Type II error rate (predict no when actually yes)
# typeIIerr <- fn / act_no
# # Sensitivity (how often correctly predict yes), ie, 1-Type II error
# sens <- tp / act_no
# # null error rate (percent wrong if always predicted the majority class)
# null_err <- act_no / tot
# # Cohen's Kappa (how well classifier performs compared to chance)
# cohen <- acc - null_err


### The model has a pretty low specificity--doesn't predict unconscious (0) well,
### but it has a high sensitivity--predicts conscious (1) well.







### stuff from the sedation cleaning
# only propofol
# propdat <- seddat %>%
#   filter(Medication_Name == "Propofol +R DRIP")

# join propofol data
# tstdat <- sqldf::sqldf(
#   "select t1.MRN,
#    t1.MRI_eval_note_date,
#    t2.Medication_Name,
#    t2.Dispensed_Amt,
#    t2.Administered_Dtm
#   
#    from ichmri as t1
#   
#    left join propdat as t2
#    on t1.MRN = t2.MRN
#    --and t2.Administered_Dtm between DATEADD(hh, 6, t1.MRI_eval_note_date) and DATEADD(hh, -6, t1.MRI_eval_note_date)
#    --and t2.Administered_Dtm between t1.MRI_eval_note_date+60*60*200 and t1.MRI_eval_note_date - 60*60*200"
# )
# 
# 
# 
# aaa <- sqldf::sqldf(
#   "select Administered_Dtm,
#   dateadd(hour, 3, Administered_Dtm)
#   
#   from propdat"
# )



## use dplyr to join
# tstdat <- ichmri %>%
#   select(MRN, MRI_eval_note_date) %>%
#   ## Propofol
#   left_join(select(filter(seddat, Medication_Name=="Propofol +R DRIP"), MRN, Medication_Name, Administered_Dtm, Dispensed_Amt), by="MRN") %>%  #Medication_Name=="Propofol +R Drip" & 
#   group_by(MRN) %>%
#   filter(Administered_Dtm >= MRI_eval_note_date-(60*60*2) & Administered_Dtm <= MRI_eval_note_date)

#filter((between(Administered_Dtm, MRI_eval_note_date-(60*60*2), MRI_eval_note_date))) #| #MRI_eval_note_date-(60*60*2) <= Administered_Dtm <= MRI_eval_note_date) |
#(Medication_Name=="Midazolam HCL Inj +R+" & between(Administered_Dtm, MRI_eval_note_date-(60*60*24), MRI_eval_note_date))) #MRI_eval_note_date-(60*60*24) <= Administered_Dtm <= MRI_eval_note_date))






### Turn this into a function:
## Find optimal alpha
# using caret
lambda.grid <- 10^seq(2, -2, length=100)
alpha.grid <- seq(0, 1, length=10)

# set up cross validation method for train function
trainCtrl <- trainControl(method = "repeatedcv",
                          number = 10,
                          repeats = 5)

srchGrd <- expand.grid(.alpha=alpha.grid, .lambda=lambda.grid)


set.seed(2018)
my.train <- train(x=data.matrix(loc_analyse.raw3[,meso_ix]), y=as.factor(loc_analyse.raw3$MRI_Cs2), #all_ix[-Extra_ix_to_remove]]
                  method="glmnet",
                  tuneGrid=srchGrd,
                  trControl=trainCtrl,
                  standardize=TRUE, maxit=1000000)

my.train$bestTune