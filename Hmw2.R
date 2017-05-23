# problem 1
x=matrix(c(3,6,4,4,5,7,4,7),nrow=4,byrow = T)
#mle of mean vector is xbar. and variance is S.
X_bar<-colMeans(x) # mean vector
n<-nrow(x) 
S<- var(x)*((n-1)/n) # Sigma^2/n

#problem 2
# a) follows a chi square distribution. with 6 dof Sum of squares of standard normal variables
# b) distribution of sqrt(n)(x-mu) follows normal distribution with (0,sigma) 
# x bar follows normal distribution with (mu, sigma/n)
#c) chi square distribution
# d) it follows scaled F distribution with 6,54 dof.

# problem 3
library(MASS)
dta1<-read.csv("C:/Users/jayabharath/Desktop/STAT 636/Hmw2/used_cars.csv",sep=",",header=T)
qqnorm(x_lt <- dta1$Age); qqline(x_lt)

bc_a<- boxcox(dta1$Age~1)
title(main="Age")
bc_p<- boxcox(dta1$Price~1)
title(main="Price")

lam_a <- bc_a$x[which.max(bc_a$y)] # finding lambda for age
lam_p <- bc_p$x[which.max(bc_p$y)] # finding lambda for price

dta_tr_a <- ((dta1$Age^lam_a)-1)/lam_a # transformation for age
dta_tr_p <- ((dta1$Price^lam_p)-1)/lam_p # transformation for price

par(mfrow=c(2,2))
hist(dta_tr_a,xlab = "",main="Age") # histograms after transformation.
hist(dta_tr_p,xlab = "",main="Price")

qqnorm(dta_tr_a,main="Age") # qq plot for age
qqline(dta_tr_a)

qqnorm(dta_tr_p,main="Price") # qqplot for price
qqline(dta_tr_p)

hist(dta1$Age)
hist(dta1$Price)

# multivariate power transformation
library(car)
summary(powerTransform(cbind(dta1$Age,dta1$Price)~1,))

n<- nrow(dta1)
lambda_seq <- seq(from = -2, to = 2, length = 100)
obj <- matrix(NA, nrow = 100, ncol = 100)

csld <- colSums(log(dta1))
for(i in 1:100) {
  for(j in 1:100) {
    X_l <- dta1
    lambda <- lambda_seq[c(i, j)]
    
    for(k in 1:2) {
      if(lambda[k] != 0) {
        X_l[, k] <- (X_l[, k] ^ lambda[k] - 1) / lambda[k]
      } else {
        X_l[, k] <- log(X_l[, k])
      }
    }
    S <- var(X_l)
    
    obj[i, j] <- -(n / 2) * log(det(S)) + (lambda - 1) %*% csld
  }
}
par(mfrow = c(1, 1))
contour(lambda_seq, lambda_seq, obj, xlab = expression(lambda[1]), 
        ylab = expression(lambda[2]),xlim=c(-2,3),ylim=c(-2,3))
points(1.27, 0.03, pch = 20, cex = 2, col = "red")
text(0.3, 0.8, expression(paste(hat(lambda), "' = [1.27, 0.031]", sep = "")), lwd = 2)

lambda_seq[which(obj==max(obj),arr.ind=TRUE)]

dev.off()

#problem 4
dta2<-read.csv("C:/STAT 636/Hmw2/sweat.csv",sep=",",header=T)
dta2

#qqplot for sweat
qqnorm(dta2$Sweat,main="sweat")
qqline(dta2$Sweat)
hist(dta2$Sweat)
#qqplot for sodium
qqnorm(dta2$Sodium,main="sodium")
qqline(dta2$Sodium)
hist(dta2$Sodium)
#qqplot for potassium
qqnorm(dta2$Potassium,main="potassium")
qqline(dta2$Potassium)
hist(dta2$Potassium)

#pairwise scatter plots
pairs(dta2)

#b 
center<- colMeans(dta2)
sigma_s<- var(dta2)
eig<-eigen(sigma_s)
eis<- eigen(solve(sigma_s))
# axes and half lengths
pri_axis<-eig$vectors[,1];pri_hl<- sqrt(eig$values[1])*sqrt(ncol(dta2)*(nrow(dta2)-1)/(nrow(dta2)*(nrow(dta2)-ncol(dta2)))*qf(.95,ncol(dta2),nrow(dta2)-ncol(dta2)))
sec_axis<-eig$vectors[,2];sec_hl<- sqrt(eig$values[2])*sqrt(ncol(dta2)*(nrow(dta2)-1)/(nrow(dta2)*(nrow(dta2)-ncol(dta2)))*qf(.95,ncol(dta2),nrow(dta2)-ncol(dta2)))
ter_axis<-eig$vectors[,3];ter_hl<- sqrt(eig$values[3])*sqrt(ncol(dta2)*(nrow(dta2)-1)/(nrow(dta2)*(nrow(dta2)-ncol(dta2)))*qf(.95,ncol(dta2),nrow(dta2)-ncol(dta2)))
pri_axis;pri_hl
sec_axis;sec_hl
ter_axis;ter_hl

# partc
a1=c(1,0,0)
a<- a1%*%center
b<- sqrt((t(a1)%*%sigma_s%*%a1)/nrow(dta2))
c<- (((nrow(dta2)-1)*ncol(dta2))/(nrow(dta2)-ncol(dta2)))*qf(.95,ncol(dta2),nrow(dta2)-ncol(dta2))
conf.intrvl_sw<- cat(a-(b*sqrt(c)),a+(b*sqrt(c))) # ci for sweat

a2=c(0,1,0)
a_s<- a2%*%center
b_s<- sqrt((t(a2)%*%sigma_s%*%a2)/nrow(dta2))
c<- (nrow(dta2)-1)*ncol(dta2)/(nrow(dta2)-ncol(dta2))*qf(.95,ncol(dta2),nrow(dta2)-ncol(dta2))
conf.intrvl_sw<- cat(a_s-(b_s*sqrt(c)),a_s+(b_s*sqrt(c))) # ci for sodium

a3=c(0,0,1)
a_k<- a3%*%center
b_k<- sqrt((t(a3)%*%sigma_s%*%a3)/nrow(dta2))
c<- (nrow(dta2)-1)*ncol(dta2)/(nrow(dta2)-ncol(dta2))*qf(.95,ncol(dta2),nrow(dta2)-ncol(dta2))
conf.intrvl_sw<- cat(a_k-(b_k*sqrt(c)),a_k+(b_k*sqrt(c))) # ci for potassium

#partd  bonferroni intervals

# bonferroni interval for sweat
bonf.CI_sw<- cat((center[1]+(qt(.05/6,19)*sqrt(sigma_s[1,1]/nrow(dta2)))),center[1]-(qt(.05/6,19)*sqrt(sigma_s[1,1]/nrow(dta2))))
#bonferroni interval for sodium
bonf.CI_na<- cat((center[2]+(qt(.05/6,19)*sqrt(sigma_s[2,2]/nrow(dta2)))),center[2]-(qt(.05/6,19)*sqrt(sigma_s[2,2]/nrow(dta2))))
#bonferroni interval for potassium
bonf.CI_k<- cat((center[3]+(qt(.05/6,19)*sqrt(sigma_s[3,3]/nrow(dta2)))),center[3]-(qt(.05/6,19)*sqrt(sigma_s[3,3]/nrow(dta2))))


# part E

library(Hotelling)
library(ICSNP)
HotellingsT2(dta2,mu=c(4,45,10))

# part f
mu<- c(4,45,10)
tstat<- sqrt(nrow(dta2)*(t(center-mu)%*%solve(sigma_s)%*%(center-mu)))
tstat
cricval<- sqrt((nrow(dta2)-1)*ncol(dta2)/(nrow(dta2)-ncol(dta2))*qf(.95,ncol(dta2),nrow(dta2)-ncol(dta2)))
cricval
tstat<cricval

# returns True, therefore point lies inside ellipse and is within 95% confidence interval
# This is in line with the result of hypothesis test.

pval<-pf(tstat*(nrow(dta2)-ncol(dta2))/((nrow(dta2)-1)*ncol(dta2)),ncol(dta2),(nrow(dta2)-ncol(dta2)))
1-pval
1-pval>.05 # returns true. Therefore we fail to reject null that mean vector is (4,45,10)





# part G
## A function to compute T2.
n <- 20
B <- 500

T2_f <- function(X, mu_0) {
  ## The covariance matrices under the null and unrestricted scenarios.
  S <- var(X)
  S_0 <- (t(X) - mu_0) %*% t(t(X) - mu_0) / (n - 1)
  
  # Compute T2 if the sample covariance matrices are non-singular.
  T2 <- NA
  if(det(S) > 0 & det(S_0) > 0) {
    Lambda <- (det(S) / det(S_0)) ^ (n / 2)
    T2 <-  Lambda
  }
  
  return(T2)
}

sim_f <- function(mu, mu_0 = c(4,45,10)) {
  ## Simulate a sample from the multivariate t.
  X <- dta2
  
  ## Observed value of T2.
  T2_0 <- T2_f(X, mu_0)
  T2_0_scaled <- (n - p) / ((n - 1) * p) * T2_0
  
  set.seed(101)
  T2_b <- rep(NA, B)
  X_0 <- t(t(X) - colMeans(X) + mu_0)
  for(b in 1:B) {
    ii <- sample(1:n, replace = TRUE)
    T2_b[b] <- T2_f(X_0[ii, ], mu_0)
  }
  T2_b_scaled <- (n - p) / ((n - 1) * p) * T2_b
  p_value_boot <- mean(T2_b_scaled >= T2_0_scaled, na.rm = TRUE)
  
  return(1-p_value_boot)
}





# problem 5
dta3<-read.csv("C:/Users/jayabharath/Desktop/STAT 636/Hmw2/peanut.csv",sep=",",header=T)
names(dta3)
attach(dta3)

Location <- as.character(Location)

Variety <- as.character(Variety)

x_1<- dta3[Location==1,-2]
x_2<- dta3[Location==2,-2]

x1_bar<- colMeans(x_1[,2:4])
x2_bar<- colMeans(x_2[,2:4])
x1var<- var(x_1[,2:4])
x2var<- var(x_2[,2:4])

y1<- as.matrix(dta3[,3:5])
factor<- dta3[,1]
treat<- as.factor(factor)
m1<- manova(y1~treat)
summary(m1,test="Wilks") #  location effect

y2<- as.matrix(dta3[,3:5])
factor1<- dta3[,2]
treat1<- as.factor(factor1)
m2<- manova(y2~treat1)
summary(m2,test="Wilks") # there is effect of variety.

#two way MANOVA
y1
m3<- manova(y1~treat+treat1+treat*treat1)
summary(m3,test="Wilks") # no interaction but there is factor 1 and factor 2 effect.

#5 b

n <- 2
p <- 3
g <- 2 
b <- 3
attach(dta3)

## Summary statistics.
x_bar <- colMeans(dta3[, 3:5])
x_bar_lk <- rbind(colMeans(dta3[Location == 1 & Variety == 5, 3:5]),
                  colMeans(dta3[Location == 2 & Variety == 5, 3:5]),
                  colMeans(dta3[Location == 1 & Variety == 6, 3:5]),
                  colMeans(dta3[Location == 2 & Variety == 6, 3:5]),
                  colMeans(dta3[Location == 1 & Variety == 8, 3:5]),
                  colMeans(dta3[Location == 2 & Variety == 8, 3:5]))



x_bar_l_dot <- rbind(colMeans(dta3[Location == 1, 3:5]), colMeans(dta3[Location == 2, 3:5]))
x_bar_dot_k <- rbind(colMeans(dta3[Variety == 5, 3:5]), colMeans(dta3[Variety == 6, 3:5]),colMeans(dta3[Variety == 8, 3:5]))



## Components for MANOVA.
SSP_cor <- SSP_fac_1 <- SSP_fac_2 <- SSP_int <- SSP_res <- matrix(0, nrow = p, ncol = p)
for(l in 1:g) {
  SSP_fac_1 <- SSP_fac_1 + b * n * t(x_bar_l_dot[l, , drop = FALSE] - x_bar) %*% 
    (x_bar_l_dot[l, , drop = FALSE] - x_bar)
}
for(k in 1:b){
  SSP_fac_2 <- SSP_fac_2 + g * n * t(x_bar_dot_k[k, , drop = FALSE] - x_bar) %*% 
    (x_bar_dot_k[k, , drop = FALSE] - x_bar)
}
for(k in 1:b) {
  for(l in 1:g){
    SSP_int <- SSP_int + n * t(x_bar_lk[(k - 1) * 2 + l, , drop = FALSE] -  x_bar_l_dot[l, , drop = FALSE] - x_bar_dot_k[k, , drop = FALSE] + x_bar) %*% 
      (x_bar_lk[(k - 1) * 2 + l, , drop = FALSE] - x_bar_l_dot[l, , drop = FALSE] - x_bar_dot_k[k, , drop = FALSE] + x_bar)
  }
}
for(l in 1:g) {
  for(k in 1:b) {
    for(r in 1:n){
      SSP_res <- SSP_res + t(as.matrix(dta3[(l - 1) * 3 * n + (k - 1) * n + r, 3:5]) - x_bar_lk[(l - 1) * 3 + k, , drop = FALSE]) %*% 
        (as.matrix(dta3[(l - 1) * 3 * n + (k - 1) * n + r, 3:5]) - x_bar_lk[(l - 1) * 3 + k, , drop = FALSE])
      SSP_cor <- SSP_cor + t(as.matrix(dta3[(l - 1) * 3 * n + (k - 1) * n + r, 3:5]) - x_bar) %*% (as.matrix(dta3[(l - 1) * 3 * n + (k - 1) * n + r, 3:5]) - x_bar)
    }
  }
}



##
## Inference.
##

## No interaction.
Lambda1 <- det(SSP_res) / det(SSP_int + SSP_res)
1 - pf((((g * b * (n - 1) - p + 1) / 2) / ((abs((g - 1) * (b - 1) - p) + 1) / 2)) * 
         (1 - Lambda1) / Lambda1, abs((g - 1) * (b - 1) - p) + 1, g * b * (n - 1) - p + 1)

## There is an effect of location
Lambda2 <- det(SSP_res) / det(SSP_fac_1 + SSP_res)
1 - pf((((g * b * (n - 1) - p + 1) / 2) / ((abs((g - 1) - p) + 1) / 2)) * 
         (1 - Lambda2) / Lambda2, abs((g - 1) - p) + 1, g * b * (n - 1) - p + 1)

## There is an effect of variety.
Lambda3 <- det(SSP_res) / det(SSP_fac_2 + SSP_res)
1 - pf((((g * b * (n - 1) - p + 1) / 2) / ((abs((b - 1) - p) + 1) / 2)) * 
         (1 - Lambda3) / Lambda3, abs((b - 1) - p) + 1, g * b * (n - 1) - p + 1)

summary(manova(y1~treat+treat1+treat*treat1), test = "Wilks")

#6

library(MASS)
hof<-read.csv("C:/Users/jayabharath/Desktop/STAT 636/Hmw2/hof_data.csv",sep=",",header=T)

num_vars <- c("H", "HR", "RBI", "AVG", "SLG", "OBP")
X <- as.matrix(hof[, num_vars]) # extracting those variables
X_st <- scale(X, center = TRUE, scale = TRUE) #standardising
DTA_st <- data.frame(hof$HOF, X_st) 

colnames(DTA_st) <- c("HOF", num_vars)
# part-a

lda_out <- lda(HOF ~ H + HR + RBI + AVG + SLG + OBP, data = DTA_st)
(lda_out$scaling) # coefficients of dicriminants. players in HOF have high H,HR and low RBI,SLG.

ld<- X_st%*%lda_out$scaling
ld[hof$HOF=="Y"]
max(ld[hof$HOF=="N"])
max(ld[hof$HOF=="Y"])
min(ld[hof$HOF=="N"])
min(ld[hof$HOF=="Y"])

# partb

t.test(X_st[,1])

nh<- table(hof$HOF)
attach(hof)

X_bar_Y_H<- mean(X_st[hof$HOF=="Y",1])
X_bar_N_H<- mean(X_st[hof$HOF=="N",1])
s_Y_H<- var(X_st[hof$HOF=="Y",1])
s_N_H<- var(X_st[hof$HOF=="N",1])
s_po<-  sqrt((((nh[2]-1)*s_Y_H^2)+ ((nh[1]-1)*s_N_H^2))/(nh[1]+nh[2]-2))
tsta_H<- (X_bar_Y_H - X_bar_N_H)/(s_po*(sqrt((1/nh[1])+(1/nh[2]))))

X_bar_Y_HR<- mean(X_st[hof$HOF=="Y",2])
X_bar_N_HR<- mean(X_st[hof$HOF=="N",2])
s_Y_HR<- var(X_st[hof$HOF=="Y",2])
s_N_HR<- var(X_st[hof$HOF=="N",2])
s_po<-  sqrt((((nh[2]-1)*s_Y_HR^2)+ ((nh[1]-1)*s_N_HR^2))/(nh[1]+nh[2]-2))
tsta_HR<- (X_bar_Y_HR - X_bar_N_HR)/(s_po*(sqrt((1/nh[1])+(1/nh[2]))))

X_bar_Y_RBI<- mean(X_st[hof$HOF=="Y",3])
X_bar_N_RBI<- mean(X_st[hof$HOF=="N",3])
s_Y_RBI<- var(X_st[hof$HOF=="Y",3])
s_N_RBI<- var(X_st[hof$HOF=="N",3])
s_po<-  sqrt((((nh[2]-1)*s_Y_RBI^2)+ ((nh[1]-1)*s_N_RBI^2))/(nh[1]+nh[2]-2))
tsta_RBI<- (X_bar_Y_RBI - X_bar_N_RBI)/(s_po*sqrt((1/nh[1])+(1/nh[2])))

X_bar_Y_AVG<- mean(X_st[hof$HOF=="Y",4])
X_bar_N_AVG<- mean(X_st[hof$HOF=="N",4])
s_Y_AVG<- var(X_st[hof$HOF=="Y",4])
s_N_AVG<- var(X_st[hof$HOF=="N",4])
s_po<-  sqrt((((nh[2]-1)*s_Y_AVG^2)+ ((nh[1]-1)*s_N_AVG^2))/(nh[1]+nh[2]-2))
tsta_AVG<- (X_bar_Y_AVG - X_bar_N_AVG)/(s_po*sqrt((1/nh[1])+(1/nh[2])))

X_bar_Y_SLG<- mean(X_st[hof$HOF=="Y",5])
X_bar_N_SLG<- mean(X_st[hof$HOF=="N",5])
s_Y_SLG<- var(X_st[hof$HOF=="Y",5])
s_N_SLG<- var(X_st[hof$HOF=="N",5])
s_po<-  sqrt((((nh[2]-1)*s_Y_SLG^2)+ ((nh[1]-1)*s_N_SLG^2))/(nh[1]+nh[2]-2))
tsta_SLG<- (X_bar_Y_SLG - X_bar_N_SLG)/(s_po*sqrt((1/nh[1])+(1/nh[2])))


X_bar_Y_OBP<- mean(X_st[hof$HOF=="Y",6])
X_bar_N_OBP<- mean(X_st[hof$HOF=="N",6])
s_Y_OBP<- var(X_st[hof$HOF=="Y",6])
s_N_OBP<- var(X_st[hof$HOF=="N",6])
s_po<-  sqrt((((nh[2]-1)*s_Y_OBP^2)+ ((nh[1]-1)*s_N_OBP^2))/(nh[1]+nh[2]-2))
tsta_OBP<- (X_bar_Y_OBP - X_bar_N_OBP)/(s_po*sqrt((1/nh[1])+(1/nh[2])))

# part-c
n_k <- table(hof$HOF)
X_bar_Y <- colMeans(X_st[hof$HOF=="Y",])
X_bar_N <- colMeans(X_st[hof$HOF=="N",])

S_Y <- var(X_st[hof$HOF=="Y",])
S_N <- var(X_st[hof$HOF=="N",])

X_bar_Y_H <- colMeans(X_st[hof$HOF=="Y",-1])
X_bar_N_H <- colMeans(X_st[hof$HOF=="N",-1])

S_Y.H <- var(X_st[hof$HOF=="Y",-1])
S_N.H <- var(X_st[hof$HOF=="N",-1])

S_pold <- (((n_k[2]-1)*S_Y)+((n_k[1]-1)*S_N))/(n_k[2]+n_k[1]-2)
S_pold_H<- (((n_k[2]-1)*S_Y.H)+((n_k[1]-1)*S_N.H))/(n_k[2]+n_k[1]-2)

T2_F<- t(X_bar_Y- X_bar_N) %*% solve(S_pold*((1/n_k[2])+(1/n_k[1]))) %*% (X_bar_Y- X_bar_N)
T2_H<- t(X_bar_Y_H - X_bar_N_H) %*% solve(S_pold_H*((1/n_k[2])+(1/n_k[1]))) %*% (X_bar_Y_H - X_bar_N_H)
 
# Partial F statistic

F_stat_H<- ((n_k[1]+n_k[2]-2-6+1)*(T2_F- T2_H))/(n_k[1]+n_k[2]-2+ T2_H)

X_bar_Y_H <- colMeans(X_st[hof$HOF=="Y",-2])
X_bar_N_H <- colMeans(X_st[hof$HOF=="N",-2])

S_Y.H <- var(X_st[hof$HOF=="Y",-2])
S_N.H <- var(X_st[hof$HOF=="N",-2])

S_pold_H<- (((n_k[2]-1)*S_Y.H)+((n_k[1]-1)*S_N.H))/(n_k[2]+n_k[1]-2)

T2_H<- t(X_bar_Y_H - X_bar_N_H) %*% solve(S_pold_H*((1/n_k[2])+(1/n_k[1]))) %*% (X_bar_Y_H - X_bar_N_H)

F_stat_HR<- ((n_k[1]+n_k[2]-2-6+1)*(T2_F- T2_H))/(n_k[1]+n_k[2]-2+ T2_H)

############
X_bar_Y_H <- colMeans(X_st[hof$HOF=="Y",-3])
X_bar_N_H <- colMeans(X_st[hof$HOF=="N",-3])

S_Y.H <- var(X_st[hof$HOF=="Y",-3])
S_N.H <- var(X_st[hof$HOF=="N",-3])

S_pold_H<- (((n_k[2]-1)*S_Y.H)+((n_k[1]-1)*S_N.H))/(n_k[2]+n_k[1]-2)

T2_H<- t(X_bar_Y_H - X_bar_N_H) %*% solve(S_pold_H*((1/n_k[2])+(1/n_k[1]))) %*% (X_bar_Y_H - X_bar_N_H)

F_stat_RBI<- ((n_k[1]+n_k[2]-2-6+1)*(T2_F- T2_H))/(n_k[1]+n_k[2]-2+ T2_H)

############
X_bar_Y_H <- colMeans(X_st[hof$HOF=="Y",-4])
X_bar_N_H <- colMeans(X_st[hof$HOF=="N",-4])

S_Y.H <- var(X_st[hof$HOF=="Y",-4])
S_N.H <- var(X_st[hof$HOF=="N",-4])

S_pold_H<- (((n_k[2]-1)*S_Y.H)+((n_k[1]-1)*S_N.H))/(n_k[2]+n_k[1]-2)

T2_H<- t(X_bar_Y_H - X_bar_N_H) %*% solve(S_pold_H*((1/n_k[2])+(1/n_k[1]))) %*% (X_bar_Y_H - X_bar_N_H)

F_stat_AVG<- ((n_k[1]+n_k[2]-2-6+1)*(T2_F- T2_H))/(n_k[1]+n_k[2]-2+ T2_H)
############

X_bar_Y_H <- colMeans(X_st[hof$HOF=="Y",-5])
X_bar_N_H <- colMeans(X_st[hof$HOF=="N",-5])

S_Y.H <- var(X_st[hof$HOF=="Y",-5])
S_N.H <- var(X_st[hof$HOF=="N",-5])

S_pold_H<- (((n_k[2]-1)*S_Y.H)+((n_k[1]-1)*S_N.H))/(n_k[2]+n_k[1]-2)

T2_H<- t(X_bar_Y_H - X_bar_N_H) %*% solve(S_pold_H*((1/n_k[2])+(1/n_k[1]))) %*% (X_bar_Y_H - X_bar_N_H)

F_stat_SLG<- ((n_k[1]+n_k[2]-2-6+1)*(T2_F- T2_H))/(n_k[1]+n_k[2]-2+ T2_H)
############
X_bar_Y_H <- colMeans(X_st[hof$HOF=="Y",-6])
X_bar_N_H <- colMeans(X_st[hof$HOF=="N",-6])

S_Y.H <- var(X_st[hof$HOF=="Y",-6])
S_N.H <- var(X_st[hof$HOF=="N",-6])

S_pold_H<- (((n_k[2]-1)*S_Y.H)+((n_k[1]-1)*S_N.H))/(n_k[2]+n_k[1]-2)

T2_H<- t(X_bar_Y_H - X_bar_N_H) %*% solve(S_pold_H*((1/n_k[2])+(1/n_k[1]))) %*% (X_bar_Y_H - X_bar_N_H)

F_stat_OBP<- ((n_k[1]+n_k[2]-2-6+1)*(T2_F- T2_H))/(n_k[1]+n_k[2]-2+ T2_H)

# part 6 d
X_st[hof$HOF=="Y",]%*%lda_out$scaling
X_st[hof$HOF=="N",]%*%lda_out$scaling
par(mfrow=c(1,1))
hist(X_st[hof$HOF=="N",]%*%lda_out$scaling,main="HOF",xlab="Discriminant value",col="red",xlim=c(-3,5),probability = T)
hist(X_st[hof$HOF=="Y",]%*%lda_out$scaling, xlab="Discriminant value",col="green",add=T,probability = T)
legend("topright",c("HOF=N","HOF=Y"),col=c("red","green"),bty="n",pch=c(15,15),pt.bg = c("red","green"),cex=1.2)

