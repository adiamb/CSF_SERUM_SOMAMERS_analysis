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
fit1 = cv.glmnet(x=preds2, y=y, type.measure = "deviance", nfolds = 10, alpha = 0.6)
plot(fit1)
fit1$lambda
fit1$lambda.min
lasso_selection=fit1$glmnet.fit$beta[which(fit1$glmnet.fit$beta[,88]!=0),88] %>% sort() %>% as.data.frame()
###loop to plot all the proteins wrt to age##########################
th=theme(axis.text.x = element_text(size = 12), axis.text.y=element_text(size=12), axis.title.x = element_text(size = 16), axis.title.y =element_text(size =16))
for (i in rownames(lasso_selection)){
  plot1= ggplot(combined_train, aes(Age, combined_train[[i]]))+geom_point(size = 3, alpha=0.4)++ggtitle(label = paste(i))+ylab(label = paste(i))+th
  
}




###########################SAMR#################################
y= combined_train$Age
x = t(preds2)
d = list(x=x, y=y, geneid=rownames(x), logged2=T)

samr.obj<-samr(d,  resp.type="Quantitative", assay.type = "array", center.arrays = F, nperms = 1000)
delta.table <- samr.compute.delta.table(samr.obj)
delta.table
delta=6
samr.plot(samr.obj,delta)
siggenes.table<-samr.compute.siggenes.table(samr.obj,delta, d, delta.table)
siggenes.table
siggenes.table$genes.up
siggenes.table$genes.lo
write.csv(siggenes.table$genes.up, file='Samr_upregulatedwithage.csv')
write.csv(siggenes.table$genes.lo, file='Samr_downregulatedwithage.csv')

require(ggplot2)
ggplot(combined_train, aes(Age, combined_train$Apo_E2_CSF))+geom_point()
ggplot(combined_train, aes(Age, `HPLN1_SERUM`))+geom_point()
csf = read_csv('~/Dropbox/CSF metabolites/CSFProteinAptamers/CSF_FINAL_Jan18.csv')
serum = read_csv('~/Dropbox/CSF metabolites/CSFProteinAptamers/SERUM_FINAL_Jan18.csv')
jap = read_csv('~/Dropbox/CSF metabolites/CSFProteinAptamers/JAPCSF_FINAL_Jan18.csv')

