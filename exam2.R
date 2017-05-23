rm(list=ls())

library(ISLR)
?Auto
names(Auto)
data1 <- Auto[,1:7]
data1_sc <- scale(data1,center = T,scale = T)

### Hierarchical Clustering using complete linkage ####

hc<- hclust(d=dist(data1_sc),method="complete")

# merge height plot ##
plot(hc$height,main="Plot of Merge Heights")

# dendogram plot
par(mfrow=c(1,1))
plot(hc)

## force the data into 3 clusters ###
hc_cluster<- cutree(hc,k=3)

# composition of clusters.
par(mfrow=c(2,4))
with(data.frame(data1_sc),boxplot(mpg~hc_cluster,xlab='cluster',ylab="mpg"))
with(data.frame(data1_sc),boxplot(cylinders~hc_cluster,xlab='cluster',ylab="cylinders"))
with(data.frame(data1_sc),boxplot(displacement~hc_cluster,xlab='cluster',ylab="displacement"))
with(data.frame(data1_sc),boxplot(horsepower~hc_cluster,xlab='cluster',ylab="horsepower"))
with(data.frame(data1_sc),boxplot(weight~hc_cluster,xlab='cluster',ylab="weight"))
with(data.frame(data1_sc),boxplot(acceleration~hc_cluster,xlab='cluster',ylab="acceleration"))
with(data.frame(data1_sc),boxplot(year~hc_cluster,xlab='cluster',ylab="year"))

## PCA ##

data1_PCA <- princomp(data1_sc)
summary(data1_PCA)

####71.58 % var is explained by first PC, 12.36 by second PC##
data1_PCA$loadings[,1:2]

### PC1 highlights the observations which have extreme engine specs. Also the four variables it highlights
# more number of cylinders = more displacement=more weight==more HP ##

new1 <- data1_sc %*% data1_PCA$loadings[,1]
new2 <- data1_sc %*% data1_PCA$loadings[,2]

### scatter plot of pc1, pc2 ###
par(mfrow=c(1,1))
plot(new1,new2,xlab="PC1",ylab="PC2",col="blue")


###2nd question ###
rm(list=ls())
library(ISLR)
names(Default)
head(Default)

### proportion of customers defaulted ###
default_cust <- sum(Default$default =="Yes")/nrow(Default)

### proportion of students ###
stu_cust <- sum(Default$student == "Yes")/nrow(Default)

### proportion of default by student status ###
tab <- table(Default$default,Default$student)
### proportion of student defaulters ###
tab[2,2]/sum(Default$student == "Yes")
### proportion of non-student defaulters ###
tab[2,1]/sum(Default$student == "No")

### Histograms for Balance and income ###
par(mfrow=c(1,2))
hist(Default$income,xlab="Income",ylab="customers",main="Histogram of Income",col="green")
hist(Default$balance,xlab="Balance",ylab="customers",main="Histogram of Balances",col="red")
summary(Default$balance)

summary(Default$income)

### scatter plot of Income vs balance###
par(mfrow=c(1,1))
plot(Default$income,Default$balance,xlab="Income",ylab="Balance",col="red",main="Income vs Balance")

### proportion of 0 balance customers###
sum(Default$balance==0)/nrow(Default)



###Logistic regression ###
rm(list=ls())
library(ISLR)
 names(Default)
logit <- glm(default~ student +balance +income + balance*income, family = binomial, data = Default)
summary(logit)
logit_pre <- predict(logit,Default,type="response",cv=T)

thresh_seq <- c(seq(from = 0, to = 0.15, length = 30),
                seq(from = 0.20, to = 1.00, length = 5))

sens <- c()
spec <- c()

for(i in 0:length(thresh_seq)){

class_def <- ifelse(logit_pre >= thresh_seq[i],"Yes","No")

sens[i] <- mean(class_def[Default$default=="Yes"]=="Yes")
spec[i] <-  mean(class_def[Default$default=="No"]=="No")

}
### ROC curve ###
par(mfrow=c(1,1))
plot(1-spec,sens,type="b",main="ROC curve",col="blue")

### threshold required ###
thresh_seq[sens==max(sens[1-spec <= 0.2])]


class_def1 <- ifelse(logit_pre >= thresh_seq[sens==max(sens[1-spec <= 0.2])],"Yes","No")
table(Default$default,class_def1)

### Results ###

### predicts default for student ###
ifelse(predict(logit, data.frame(student="Yes",balance = 1500, income= 25000),type="response") > thresh_seq[sens==max(sens[1-spec <= 0.2])],"Yes","No")

### predicts not default for this customer###
ifelse(predict(logit, data.frame(student="No",balance = 1000, income= 40000),type="response") > thresh_seq[sens==max(sens[1-spec <= 0.2])],"Yes","No")
