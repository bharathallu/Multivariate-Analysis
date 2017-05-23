dta_ex <- read.csv("C:/Users/jayabharath/Desktop/STAT 636/Exam/consumer_pref.csv",header = T,sep=",")
summary(dta_ex)
str(dta_ex)



# normality check
pairs(dta_ex) # Multi vaiate normality is good
qqnorm(dta_ex$Taste,main="taste")
qqline(dta_ex$Taste)

qqnorm(dta_ex$Value,main="value")
qqline(dta_ex$Value)

qqnorm(dta_ex$Flavor,main="flavor")
qqline(dta_ex$Flavor)

qqnorm(dta_ex$Snack,main="snack")
qqline(dta_ex$Snack)

qqnorm(dta_ex$Energy,main="energy")
qqline(dta_ex$Energy)

pairs(dta_ex) # Multi vaiate normality is good

i<- 100
# normality check

t<- matrix(rep(0,100),byrow = T,ncol=1)
for (i in 1:100){
  t[i]<- ((as.matrix(dta_ex[i,])-as.matrix(t(center)))%*%solve(sigma_s)%*%t(as.matrix(dta_ex[i,])-as.matrix(t(center)))<= qchisq(.5,5))
}
mean(t) # approximately 50% of data is less than chi square, therefore multivariate normal.

mean(as.matrix(t(dta_ex-center))%*%solve(sigma_s)%*%(as.matrix(dta_ex-center)) <= qchisq(.5,5))

# hotelling T square
library(Hotelling)
library(ICSNP)
HotellingsT2(dta_ex,mu=c(0,0,0,0,0))

# Eigen values & vectors
center<- colMeans(dta_ex)
sigma_s<- var(dta_ex)
eig<-eigen(sigma_s)
eis<- eigen(solve(sigma_s))
# axes and half lengths
pri_axis<-eig$vectors[,1];pri_hl<- sqrt(eig$values[1])*sqrt(ncol(dta_ex)*(nrow(dta_ex)-1)/(nrow(dta_ex)*(nrow(dta_ex)-ncol(dta_ex)))*qf(.95,ncol(dta_ex),nrow(dta_ex)-ncol(dta_ex)))
sec_axis<-eig$vectors[,2];sec_hl<- sqrt(eig$values[2])*sqrt(ncol(dta_ex)*(nrow(dta_ex)-1)/(nrow(dta_ex)*(nrow(dta_ex)-ncol(dta_ex)))*qf(.95,ncol(dta_ex),nrow(dta_ex)-ncol(dta_ex)))
ter_axis<-eig$vectors[,3];ter_hl<- sqrt(eig$values[3])*sqrt(ncol(dta_ex)*(nrow(dta_ex)-1)/(nrow(dta_ex)*(nrow(dta_ex)-ncol(dta_ex)))*qf(.95,ncol(dta_ex),nrow(dta_ex)-ncol(dta_ex)))
qua_axis<-eig$vectors[,4];qua_hl<- sqrt(eig$values[4])*sqrt(ncol(dta_ex)*(nrow(dta_ex)-1)/(nrow(dta_ex)*(nrow(dta_ex)-ncol(dta_ex)))*qf(.95,ncol(dta_ex),nrow(dta_ex)-ncol(dta_ex)))
pen_axis<-eig$vectors[,5];pen_hl<- sqrt(eig$values[5])*sqrt(ncol(dta_ex)*(nrow(dta_ex)-1)/(nrow(dta_ex)*(nrow(dta_ex)-ncol(dta_ex)))*qf(.95,ncol(dta_ex),nrow(dta_ex)-ncol(dta_ex)))
pri_axis;pri_hl
sec_axis;sec_hl
ter_axis;ter_hl
qua_axis;qua_hl
pen_axis;pen_hl

# ci without large sample inference

a1=c(1,0,0,0,0)
a<- a1%*%center
b<- sqrt((t(a1)%*%sigma_s%*%a1)/nrow(dta_ex))
c<- (((nrow(dta_ex)-1)*ncol(dta_ex))/(nrow(dta_ex)-ncol(dta_ex)))*qf(.95,ncol(dta_ex),nrow(dta_ex)-ncol(dta_ex))
conf.intrvl_sw<- cat(a-(b*sqrt(c)),a+(b*sqrt(c))) # ci for taste

a2=c(0,1,0,0,0)
a_s<- a2%*%center
b_s<- sqrt((t(a2)%*%sigma_s%*%a2)/nrow(dta_ex))
c<- (nrow(dta_ex)-1)*ncol(dta_ex)/(nrow(dta_ex)-ncol(dta_ex))*qf(.95,ncol(dta_ex),nrow(dta_ex)-ncol(dta_ex))
conf.intrvl_sw<- cat(a_s-(b_s*sqrt(c)),a_s+(b_s*sqrt(c))) # ci for value

a3=c(0,0,1,0,0)
a_k<- a3%*%center
b_k<- sqrt((t(a3)%*%sigma_s%*%a3)/nrow(dta_ex))
c<- (nrow(dta_ex)-1)*ncol(dta_ex)/(nrow(dta_ex)-ncol(dta_ex))*qf(.95,ncol(dta_ex),nrow(dta_ex)-ncol(dta_ex))
conf.intrvl_sw<- cat(a_k-(b_k*sqrt(c)),a_k+(b_k*sqrt(c))) # ci for flavor

a1=c(0,0,0,1,0)
a<- a1%*%center
b<- sqrt((t(a1)%*%sigma_s%*%a1)/nrow(dta_ex))
c<- (((nrow(dta_ex)-1)*ncol(dta_ex))/(nrow(dta_ex)-ncol(dta_ex)))*qf(.95,ncol(dta_ex),nrow(dta_ex)-ncol(dta_ex))
conf.intrvl_sw<- cat(a-(b*sqrt(c)),a+(b*sqrt(c))) # ci for snack

a1=c(0,0,0,0,1)
a<- a1%*%center
b<- sqrt((t(a1)%*%sigma_s%*%a1)/nrow(dta_ex))
c<- (((nrow(dta_ex)-1)*ncol(dta_ex))/(nrow(dta_ex)-ncol(dta_ex)))*qf(.95,ncol(dta_ex),nrow(dta_ex)-ncol(dta_ex))
conf.intrvl_sw<- cat(a-(b*sqrt(c)),a+(b*sqrt(c))) # ci for energy

# bonferroni without large sample inference

# bonferroni interval for taste
bonf.CI_sw<- cat((center[1]+(qt(.05/10,99)*sqrt(sigma_s[1,1]/nrow(dta_ex)))),center[1]-(qt(.05/10,99)*sqrt(sigma_s[1,1]/nrow(dta_ex))))
#bonferroni interval for value
bonf.CI_na<- cat((center[2]+(qt(.05/10,99)*sqrt(sigma_s[2,2]/nrow(dta_ex)))),center[2]-(qt(.05/10,99)*sqrt(sigma_s[2,2]/nrow(dta_ex))))
#bonferroni interval for flavor
bonf.CI_k<- cat((center[3]+(qt(.05/10,99)*sqrt(sigma_s[3,3]/nrow(dta_ex)))),center[3]-(qt(.05/10,99)*sqrt(sigma_s[3,3]/nrow(dta_ex))))
#bonferroni interval for snack
bonf.CI_k<- cat((center[4]+(qt(.05/10,99)*sqrt(sigma_s[4,4]/nrow(dta_ex)))),center[3]-(qt(.05/10,99)*sqrt(sigma_s[4,4]/nrow(dta_ex))))
#bonferroni interval for energy
bonf.CI_k<- cat((center[5]+(qt(.05/10,99)*sqrt(sigma_s[5,5]/nrow(dta_ex)))),center[3]-(qt(.05/10,99)*sqrt(sigma_s[5,5]/nrow(dta_ex))))



# bootstrap T2

## A function to compute T2.
X<- dta_ex
n <- 100
B <- 1000
p<- 5
T2_f <- function(X, mu_0) {
  ## The covariance matrices under the null and unrestricted scenarios.
  S <- var(X)
  S_0 <- (t(X) - mu_0) %*% t(t(X) - mu_0) / (n - 1)
  
  # Compute T2 if the sample covariance matrices are non-singular.
  T2 <- NA
  if(det(S) > 0 & det(S_0) > 0) {
    Lambda <- (det(S) / det(S_0)) ^ (n / 2)
    T2 <- (n - 1) * (1 / (Lambda ^ (2 / n)) - 1)
  }
  
  return(T2)
}

sim_f <- function(mu, mu_0 = c(0,0,0,0,0)) {
  ## Simulate a sample from the multivariate t.
  X <- dta_ex
  
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
  
  return(p_value_boot)
}

