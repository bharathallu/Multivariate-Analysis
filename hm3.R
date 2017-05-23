# 2 problem PCA


X_bar <- colMeans(dta)
S<- var(dta)
egn_S <- eigen(S)
# PCA
dta_pca <- princomp(dta)
names(dta_pca)
summary(dta_pca)
dta_pca$loadings
## 53% 0f variance is explained by first PCA
## 27% of variance is explained by second PCA

dim(dta_pca$loadings)
dta_pca$loadings[,1]
ee$vectors[,1]
## the first PCA is equal to the first eigen vector of S
dta_pca$loadings[,2]
egn_S$vectors[,2]
## the second PCA is equal to the second eigen vector of S


PC_S <- as.matrix(dta) %*% egn_S$vectors
PC_var<- var(PC_S)
PC_var[1,1]
egn_S$values[1]
# variance of first PC is equal to the first eigen value of S
PC_var[2,2]
egn_S$values[2]
# variance of second PC is equal to the second eigen value of S

PCA1<- as.matrix(dta)%*%dta_pca$loadings[,1]
PCA2<- as.matrix(dta)%*%dta_pca$loadings[,2]
par(mfrow=c(1,1))
plot(PCA1,PCA2)
# multivariate normality from PCA

which.min(PCA_1) # 63 week  The weighted differnce of most significant variables exxon and shell is min.
which.max(PCA_1) # 56 week . The weighted differnce of most significant variables exxon and shell is max
which.min(PCA_2) # 63 week  The weighted differnce of most significant variables exxon and shell is min.
which.max(PCA_2)


#1
rm(list=ls())
dta_hof<- read.csv(file = "C:/STAT 636/Data/hof_data.csv",header = T,sep = ",")

num_vars <- c("H", "HR", "RBI", "AVG", "SLG", "OBP")
X <- as.matrix(dta_hof[, num_vars])
X_st <- scale(X, center = TRUE, scale = TRUE)
DTA_st <- data.frame(dta_hof$HOF, X_st)
colnames(DTA_st) <- c("HOF", num_vars)

library(MASS)

kappa <- seq(from = 0, to = .5, by = 0.01)

lda_out <- lda(HOF ~ H + HR + RBI + AVG + SLG + OBP, data = DTA_st, CV = TRUE)

qda_out <- qda(HOF ~ H + HR + RBI + AVG + SLG + OBP, data = DTA_st, CV = TRUE)

bacc_lda <- sens_lda <- spec_lda <- ppv_lda <- npv_lda <- numeric(length(kappa))

bacc_qda <- sens_qda <- spec_qda <- ppv_qda <- npv_qda <- numeric(length(kappa))

for(i in 1:length(kappa)) { 
  
  cls_lda <- lda_out$posterior[, 2] > kappa[i]
  
  sens_lda[i] <- mean(cls_lda[DTA_st$HOF=="Y"]==T)
  spec_lda[i] <- mean(cls_lda[DTA_st$HOF=="N"]==F)
  ppv_lda[i] <- mean(DTA_st$HOF[cls_lda==T]=="Y")
  npv_lda[i] <- mean(DTA_st$HOF[cls_lda==F]=="N")
  bacc_lda[i] <- (sens_lda[i]+3*spec_lda[i])/4

  cls_qda <- qda_out$posterior[, 2] > kappa[i]
  
  sens_qda[i] <- mean(cls_qda[DTA_st$HOF=="Y"]==T)
  spec_qda[i] <- mean(cls_qda[DTA_st$HOF=="N"]==F)
  ppv_qda[i] <- mean(DTA_st$HOF[cls_qda==T]=="Y")
  npv_qda[i] <- mean(DTA_st$HOF[cls_qda==F]=="N")
  bacc_qda[i] <- (sens_qda[i]+3*spec_qda[i])/4
  
}


plot(kappa,bacc_lda,type="l",ylim=c(0.8,0.95),lwd=5,col="blue",main = "BACC",ylab = "Bacc")
lines(kappa,bacc_qda,type="l",ylim=c(0,8,0.95),lwd=5,col="green")

kappa[which.max(bacc_lda)]
kappa[which.max(bacc_qda)]

sens_lda[which.max(bacc_lda)] 
spec_lda[which.max(bacc_lda)] 
ppv_lda[which.max(bacc_lda)] 
npv_lda[which.max(bacc_lda)] 
bacc_lda[which.max(bacc_lda)] 

sens_qda[which.max(bacc_qda)] 
spec_qda[which.max(bacc_qda)] 
ppv_qda[which.max(bacc_qda)] 
npv_qda[which.max(bacc_qda)] 
bacc_qda[which.max(bacc_qda)] 


