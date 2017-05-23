####
#### For Analytics.
####

## Load data.
DTA <- read.csv("hof_data.csv")

## Extract a few offensive statistics (numerical variables).
num_vars <- c("H", "HR", "RBI", "AVG", "SLG", "OBP")
X <- as.matrix(DTA[, num_vars])
X_st <- scale(X, center = TRUE, scale = TRUE)
DTA_st <- data.frame(DTA$HOF, X_st)
colnames(DTA_st) <- c("HOF", num_vars)

p <- ncol(X)

## Summary statistics.
x_bar <- colMeans(X)
S <- var(X)
R <- cor(X)

##
## Principal component analysis.
##

pca <- prcomp(X, center = TRUE, scale = TRUE)
pca
summary(pca)

## Eigenvalues and eigenvectors of R.
egn_R <- eigen(R)
egn_R

## Compute PCs.
PC <- scale(X, center = TRUE, scale = TRUE) %*% egn_R$vectors
apply(PC, 2, var)
egn_R$values

##
## Linear discriminant analysis.
##

library(MASS)

lda_out <- lda(HOF ~ H + HR + RBI + AVG + SLG + OBP, data = DTA_st)
lda_pred <- predict(lda_out, newdata = DTA_st)
LDs <- lda_pred$x
our_linear_combo <- X_st %*% rep(1 / 6, 6)

par(mfrow = c(1, 2))
boxplot(LDs[DTA$HOF == "Y", 1], LDs[DTA$HOF == "N", 1], ylim = c(-3, 8))
boxplot(our_linear_combo[DTA$HOF == "Y"], our_linear_combo[DTA$HOF == "N"], 
  ylim = c(-3, 8))

cor(LDs[, 1], DTA_st$H)
cor(LDs[, 1], DTA_st$HR)
cor(LDs[, 1], DTA_st$RBI)
cor(LDs[, 1], DTA_st$AVG)
cor(LDs[, 1], DTA_st$SLG)
cor(LDs[, 1], DTA_st$OBP)

## Partial F statistic for H.
n_k <- table(DTA$HOF)
X_bar_Y <- colMeans(X_st[DTA$HOF == "Y", ])
X_bar_N <- colMeans(X_st[DTA$HOF == "N", ])
X_bar_Y_red <- X_bar_Y[-1]
X_bar_N_red <- X_bar_N[-1]
S_Y <- var(X_st[DTA$HOF == "Y", ])
S_N <- var(X_st[DTA$HOF == "N", ])
S_Y_red <- var(X_st[DTA$HOF == "Y", -1])
S_N_red <- var(X_st[DTA$HOF == "N", -1])
S_po <- ((n_k[2] - 1) * S_Y + (n_k[1] - 1) * S_N) / (n_k[1] + n_k[2] - 2)
S_po_red <- ((n_k[2] - 1) * S_Y_red + (n_k[1] - 1) * S_N_red) / (n_k[1] + n_k[2] - 2)

T2_full <- t(X_bar_Y - X_bar_N) %*% solve(S_po) %*% (X_bar_Y - X_bar_N) / 
  (1 / n_k[1] + 1 / n_k[2])
T2_red <- t(X_bar_Y_red - X_bar_N_red) %*% solve(S_po_red) %*% 
  (X_bar_Y_red - X_bar_N_red) / (1 / n_k[1] + 1 / n_k[2])

F_stat_H <- (n_k[1] + n_k[2] - 2 - 6 + 1) * ((T2_full - T2_red) / (n_k[1] + n_k[2] - 2 + T2_red))


library(Hotelling)

hot_full <- hotelling.test(H + HR + RBI + AVG + SLG + OBP ~ HOF, data = DTA_st)
hot_red <- hotelling.test(HR + RBI + AVG + SLG + OBP ~ HOF, data = DTA_st)


(n_1 + n_2 - p - 1) / ((n_1 + n_2 - 2) * p) * T2



## Picture of first LD.
with(DTA, boxplot(LDs[HOF == "Y"], LDs[HOF == "N"]))

## Investigate players with high values of LD who are not in HOF.
ii <- (1:nrow(DTA))[DTA$HOF == "N" & LDs > 2]
DTA_st[ii, ]


## Compute linear discriminant.
ld <- X_st %*% lda_out$scaling
with(DTA, boxplot(ld[HOF == "Y"], ld[HOF == "N"]))

oo <- order(abs(ld), decreasing = TRUE)
data.frame(DTA_st$HOF[oo[1:20]], round(ld[oo[1:20]], 2), round(X_st[oo[1:20], ], 2))

## Resubstitution error.
HOF_hat <- predict(lda_out, data = DTA)
with(DTA, table(HOF, HOF_hat$class))

####
#### For DL.
####

## Load data.
DTA <- read.csv("HOF_tr.csv")

n <- nrow(DTA)

##
## PCA on some of the numerical variables.
##

## These are a few offensive statistics. Big numbers signal greater offensive production. 
num_vars <- c("H", "HR", "RBI", "AVG", "SLG", "OBP")
X <- as.matrix(DTA[, num_vars])
head(X)

p <- ncol(X)

## Standardized variables.
X_st <- scale(X, center = TRUE, scale = TRUE)
DTA_st <- data.frame(DTA$HOF, X_st)
colnames(DTA_st) <- c("HOF", num_vars)

## Summary statistics.
x_bar <- colMeans(X)
S <- var(X)
R <- cor(X)

## Eigenvalues and eigenvectors of R.
ee <- eigen(R)
lambda <- ee$values
ee <- ee$vectors

## The first three PCs explain about 94% of variance. The first PC represents a weighted 
## average of all the variables and will distinguish players based on their overall 
## stats. The second PC represents a weighted difference, comparing AVG and OBP to HR, 
## RBI, and SLG. Some players may not have exceptional "power" statistics (such as home 
## runs, runs batted in, and slugging percentage, but they may nevertheless have been 
## very valuable offensive players (e.g. Rod Carew, with high batting average and on base 
## percentage but relatively few home runs or RBIs). The third PC is harder to interpret. 
## We'll talk about it later if we have time.
pca <- prcomp(X_st)
pca
summary(pca)
plot(pca)

## Extract the first three PCs.
Y <- X_st %*% pca$rotation[, 1:3]
head(Y)

## Variance of PCs equal lambdas.
apply(Y, 2, var)
lambda[1:3]

## Correlations with variables...
cor(Y[, 1], X_st[, 1])
ee[1, 1] * sqrt(lambda[1])

ee * kronecker(matrix(1, nrow = p, ncol = 1), sqrt(matrix(lambda, nrow = 1, ncol = p)))

##
## Linear discriminant analysis.
##

library(MASS)

lda_out <- lda(HOF ~ H + HR + RBI + AVG + SLG + OBP, data = DTA_st)

## Compute linear discriminants.
ld <- X_st %*% lda_out$scaling

## Visualize linear discriminants.
with(DTA_st, hist(ld[HOF == "N"], prob = TRUE, xlim = c(-2.5, 6), ylim = c(0, 0.55), 
  xlab = "linear discriminants", main = "", col = "pink"))
with(DTA_st, hist(ld[HOF == "Y"], prob = TRUE, xlim = c(-2.5, 6), ylim = c(0, 0.55), 
  col = "lightgreen", add = TRUE))

## Investigate extreme values of linear discriminants.
ld_oo <- order(abs(ld))
DTA[head(ld_oo), ]
DTA[tail(ld_oo), ]

## Resubstitution error.
HOF_hat <- predict(lda_out, data = DTA_st)
with(DTA_st, table(HOF, HOF_hat$class))



