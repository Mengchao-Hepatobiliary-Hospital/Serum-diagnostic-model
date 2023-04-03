#Data in file"Primary data.Rdata" is used to run the code in file"Primary code.R" 
#Loading data
load("./Primary data.Rdata")


#LVQ model for feature selection in HCC-related differentially abundant proteins
#The data in file "DAPdata" is DIA-MS quantitative profiles of 34 HCC-related differentially abundant proteins after pre-processing and median normalization.
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


#Randomforest model
#The quantitative profiles in file "PRMdata" is exactly the same as that in file "Table S6. Targeted proteomic abundance profiles of 325 samples.csv".
#You can run the code using "RPMdata" directly or after tidying and simplifying "Table S6"
#Loading packages and and preprocessing data
#install.packages("randomForest")
#install.packages("pROC")
library(randomForest)
library(pROC)
Traingset_data<-PRMdata[PRMdata$Data.Set=="Training set",]
Testingset_data<-PRMdata[PRMdata$Data.Set=="Testing set",]
Validationset_data<-PRMdata[PRMdata$Data.Set=="Validation set",]
colnames(Traingset_data)[3]<-"rank"
colnames(Testingset_data)[3]<-"rank"
colnames(Validationset_data)[3]<-"rank"
Traingset_data$rank<-ifelse(grepl("HCC",Traingset_data$rank),1,0)
Testingset_data$rank<-ifelse(grepl("HCC",Testingset_data$rank),1,0)
Validationset_data$rank<-ifelse(grepl("HCC",Validationset_data$rank),1,0)

#Training, testing, verification of P4 model
#P4 model training
P4_randomforest <- randomForest(rank~HABP2+CD163+AFP+PIVKA.II,
                                  data = Traingset_data,
                                  ntree =500,
                                  mtry=2,
                                  importance=TRUE,
                                  proximity=TRUE)
P4_pre_ran1 <- predict(P4_randomforest,newdata=Traingset_data)
P4_obs_p_ran1 = data.frame(P4_pre_ran1,obs=Traingset_data$rank)
#P4 model correction by testing set
P4_pre_ran2 <- predict(P4_randomforest,newdata=Testingset_data)
P4_obs_p_ran2 = data.frame(P4_pre_ran2,obs=Testingset_data$rank)
#Verification of P4 model in validation set
P4_pre_ran3 <- predict(P4_randomforest,newdata=Validationset_data)
P4_obs_p_ran3 = data.frame(P4_pre_ran3,obs=Validationset_data$rank)

#Training, testing, verification of AFP+PIVKA.II model
#AFP+PIVKA.II model training
AP_randomforest <- randomForest(rank~AFP+PIVKA.II,
                                data = Traingset_data,
                                ntree =500,
                                mtry=1,
                                importance=TRUE,
                                proximity=TRUE)
AP_pre_ran1 <- predict(AP_randomforest,newdata=Traingset_data)
AP_obs_p_ran1 = data.frame(AP_pre_ran1,obs=Traingset_data$rank)
#P4 model correction by testing set
AP_pre_ran2 <- predict(AP_randomforest,newdata=Testingset_data)
AP_obs_p_ran2 = data.frame(AP_pre_ran2,obs=Testingset_data$rank)
#Verification of P4 model in validation set
AP_pre_ran3 <- predict(AP_randomforest,newdata=Validationset_data)
AP_obs_p_ran3 = data.frame(AP_pre_ran3,obs=Validationset_data$rank)

#Visualization of ROC curve
#Visualization of ROC curve in training set
auc_AFP1 <- roc(rank~AFP,data=Traingset_data, smooth=T)
plot(auc_AFP1,col="green3",print.auc=F,font=7,cex.axis=1.4,cex.lab=1.4,cex.main=1.6,
     max.auc.polygon=F,print.thres=F,legend=T,main="Training Set")
auc_PIVKAII1 <- roc(rank~PIVKA.II,data=Traingset_data, smooth=T)
lines(auc_PIVKAII1,col="blue")
auc_AP1 <- roc(obs~AP_pre_ran1,data=AP_obs_p_ran1, smooth=F)
lines(auc_AP1,col="darkorchid")
auc_P4_1 <- roc(obs~P4_pre_ran1,data=P4_obs_p_ran1, smooth=F)
lines(auc_P4_1,col="red2")
legend("bottomright",legend = c("P4",
                                "AFP+PIVKA-II","PIVKA-II","AFP"),
       fill = c("red2","darkorchid","blue","green3"),
       col = c("red2","darkorchid","blue","green3"),
       text.font = 7,cex = 1.4)

#Visualization of ROC curve in testing set
auc_AFP2 <- roc(rank~AFP,data=Testingset_data, smooth=T)
plot(auc_AFP2,col="green3",print.auc=F,font=7,cex.axis=1.4,cex.lab=1.4,cex.main=1.6,
     max.auc.polygon=F,print.thres=F,legend=T,main="Testing Set")
auc_PIVKAII2 <- roc(rank~PIVKA.II,data=Testingset_data, smooth=T)
lines(auc_PIVKAII2,col="blue")
auc_AP2 <- roc(obs~AP_pre_ran2,data=AP_obs_p_ran2, smooth=T)
lines(auc_AP2,col="darkorchid")
auc_P4_2 <- roc(obs~P4_pre_ran2,data=P4_obs_p_ran2, smooth=T)
lines(auc_P4_2,col="red2")
legend("bottomright",legend = c("P4",
                                "AFP+PIVKA-II","PIVKA-II","AFP"),
       fill = c("red2","darkorchid","blue","green3"),
       col = c("red2","darkorchid","blue","green3"),
       text.font = 7,cex = 1.4)

#Visualization of ROC curve in validation set
auc_AFP3 <- roc(rank~AFP,data=Validationset_data, smooth=T)
plot(auc_AFP3,col="green3",print.auc=F,font=7,cex.axis=1.4,cex.lab=1.4,cex.main=1.6,
     max.auc.polygon=F,print.thres=F,legend=T,main="Validation Set")
auc_PIVKAII3 <- roc(rank~PIVKA.II,data=Validationset_data, smooth=T)
lines(auc_PIVKAII3,col="blue")
auc_AP3 <- roc(obs~AP_pre_ran3,data=AP_obs_p_ran3, smooth=T)
lines(auc_AP3,col="darkorchid")
auc_P4_3 <- roc(obs~P4_pre_ran3,data=P4_obs_p_ran3, smooth=T)
lines(auc_P4_3,col="red2")
legend("bottomright",legend = c("P4",
                                "AFP+PIVKA-II","PIVKA-II","AFP"),
       fill = c("red2","darkorchid","blue","green3"),
       col = c("red2","darkorchid","blue","green3"),
       text.font = 7,cex = 1.4)

#Validation and visualization of ROC curve in prospective validation set of HCC early diagnosis
colnames(prospective_data)[3]<-"rank"
prospective_data$rank<-ifelse(grepl("HCC",prospective_data$rank),1,0)
pros_randomforest <- randomForest(rank~HABP2+CD163+AFP+PIVKA.II,
                                data = Traingset_data,
                                ntree =500,
                                mtry=2,
                                importance=TRUE,
                                proximity=TRUE)
pros_pre_ran <- predict(pros_randomforest,newdata=prospective_data)
pros_obs_p_ran = data.frame(pros_pre_ran,obs=prospective_data$rank)
auc_pros <- roc(obs~pros_pre_ran,data=pros_obs_p_ran, smooth=T)
plot(auc_pros,col="red2",print.auc=F,font=7,cex.axis=1.4,cex.lab=1.4,cex.main=1.6,
     max.auc.polygon=F,print.thres=F,legend=T,main="Prospective Validation Set")
