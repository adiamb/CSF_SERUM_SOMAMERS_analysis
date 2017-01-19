require(data.table);require(readr);require(plyr);require(dplyr);require(samr);require(glmnet);require(limma)
combined = read_csv('~/Dropbox/CSF metabolites/CSFProteinAptamers/combined_FINAL_Jan18.csv')
combined$Age = as.numeric(combined$Age)
combined$BMI = as.numeric(combined$BMI)
combined$`Interval from sleepiness(year)` = as.numeric(combined$`Interval from sleepiness(year)`)
combined$`Interval from cataplexy(year)`=as.numeric(combined$`Interval from cataplexy(year)`)
combined_train = combined[!is.na(combined$Age),]
combined_test = combined[is.na(combined$Age),]


##########glmnet to predict missing ages ####################
combined_train$Gender[combined_train$Gender == "m"] = "M"
combined_train$Gender[combined_train$Gender == "."] = "M"
combined_train$Gender[is.na(combined_train$Gender)] = "M"
combined_train$Race[combined_train$Race!="Caucasian"] = "others"
combined_train$BMI[is.na(combined_train$BMI)] = mean(combined_train$BMI, na.rm =T)
preds = model.matrix(Age~Gender+Race,data = combined_train[-c(1:4,  9, 10)])[,-1]
preds2 = as.matrix(cbind(preds, combined_train[-c(1:6, 9, 7, 10)]))
y = combined_train$Age
fit1 = cv.glmnet(x=preds2, y=y, type.measure = "deviance", nfolds = 10, alpha = 1)
plot(fit1)
fit1$lambda
fit1$lambda.min
fit1$glmnet.fit$beta[which(fit1$glmnet.fit$beta[,97]!=0),97] %>% sort() 









csf = read_csv('~/Dropbox/CSF metabolites/CSFProteinAptamers/CSF_FINAL_Jan18.csv')
serum = read_csv('~/Dropbox/CSF metabolites/CSFProteinAptamers/SERUM_FINAL_Jan18.csv')
jap = read_csv('~/Dropbox/CSF metabolites/CSFProteinAptamers/JAPCSF_FINAL_Jan18.csv')

