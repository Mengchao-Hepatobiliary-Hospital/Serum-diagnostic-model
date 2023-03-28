#The data in file "DAPdata" is DIA-MS quantitative profiles of 34 HCC-related differentially abundant proteins after pre-processing and median normalization.
load("./DAPdata.Rdata")

#Loading packages and preprocessing data
#install.packages("mlbench")
#install.packages("caret")
library(mlbench)
library(caret)
colnames(DAPdata)[1]<-"rank"
DAPdata$rank<-ifelse(grepl("HCC",DAPdata$rank),1,0)
DAPdata$rank<-as.factor(DAPdata$rank)
# prepare training scheme
control <- trainControl(method="repeatedcv", number=10, repeats=10)
# train the LVQ model
LVQmodel <- train(rank~., data=DAPdata, method="lvq", preProcess="scale", trControl=control)
# estimate variable importance
importance <- varImp(LVQmodel, scale=FALSE)
#save importance
write.csv(importance$importance,"importance.csv")
# summarize importance
print(importance)
# plot importance
plot(importance)
