## These data come from gene expression microarray experiments on small, round blue-cell 
## tumors (SRCTs), described in Khan et al (2001). Nature Medicine, 7: 673-679.  There 
## are four different SRBCT tumor types: neuroblastoma (NB), rhabdomyosarcoma (RMS), non-
## Hodgkin lymphoma (NHL), and the Ewing family of tumors (EWS).  The column names 
## indicate the tumor type, with 1 corresponding to NHL, 2 to EWS, 3 to NB, and 4 to RMS.

## Load data.
data = as.matrix(read.csv("khan.csv"))
id = colnames(data)
for(i in 1:ncol(data))
  id[i] = substr(id[i], 2, 2)

m = nrow(data)
n = ncol(data)
n_k = as.numeric(table(id))

## Fix up the column names
c_names_1 = rep("1", n_k[1])
c_names_2 = rep("2", n_k[2])
c_names_3 = rep("3", n_k[3])
c_names_4 = rep("4", n_k[4])

for(i in 1:n_k[1])
  c_names_1[i] = paste(c_names_1[i], "_", i, sep = "")
for(i in 1:n_k[2])
  c_names_2[i] = paste(c_names_2[i], "_", i, sep = "")
for(i in 1:n_k[3])
  c_names_3[i] = paste(c_names_3[i], "_", i, sep = "")
for(i in 1:n_k[4])
  c_names_4[i] = paste(c_names_4[i], "_", i, sep = "")
c_names = c(c_names_1, c_names_2, c_names_3, c_names_4)
tmp = id
tmp[id == 1] = c_names_1
tmp[id == 2] = c_names_2
tmp[id == 3] = c_names_3
tmp[id == 4] = c_names_4
c_names = tmp
colnames(data) = c_names

## Basic summary statistics
data = as.matrix(data)
X = model.matrix(~ factor(id) - 1)
grp_means = t(t(data %*% X) / n_k)
grp_sds = sqrt(t(t(((data - grp_means %*% t(X)) ^ 2) %*% X) / (n_k - 1)))

####
#### EDA.
####

#### Boxplots and scatterplots 

## Some example gene-specific boxplots.
par(mfrow = c(2, 2))
boxplot(data[1, id == 1], data[1, id == 2], data[1, id == 3], data[1, id == 4])
boxplot(data[2, id == 1], data[2, id == 2], data[2, id == 3], data[2, id == 4])
boxplot(data[3, id == 1], data[3, id == 2], data[3, id == 3], data[3, id == 4])
boxplot(data[4, id == 1], data[4, id == 2], data[4, id == 3], data[4, id == 4])

## Array-specific boxplots by groups.
par(mfrow = c(2, 2))
boxplot(data[, id == 1], ylim = c(-6, 3), main = "Tumor type 1")
boxplot(data[, id == 2], ylim = c(-6, 3), main = "Tumor type 2")
boxplot(data[, id == 3], ylim = c(-6, 3), main = "Tumor type 2")
boxplot(data[, id == 4], ylim = c(-6, 3), main = "Tumor type 2")

## Overall comparisons of groups, in terms of both means and s.d.s.
par(mfrow = c(1, 2))
boxplot(grp_means[, 1], grp_means[, 2], grp_means[, 3], grp_means[, 4], names = 1:4, 
  main = "Means by tumor type")
boxplot(grp_sds[, 1], grp_sds[, 2], grp_sds[, 3], grp_sds[, 4], names = 1:4, 
  main = "S.D.s by tumor type")

## Pairwise scatterplots
par(mfrow = c(1, 2))
plot(data[, 1], data[, 2], pch = 20, col = "grey", xlab = "Array 1", ylab = "Array 2")
abline(0, 1, lty = 2)
lines(lowess(data[, 1], data[, 2]), lwd = 2, col = "blue")
plot((data[, 1] + data[, 2]) / 2, data[, 2] - data[, 1], pch = 20, col = "grey", 
  xlab = "A: (Array 1 + Array 2) / 2", ylab = "M: Array 2 - Array 1")
abline(0, 0, lty = 2)
lines(lowess((data[, 1] + data[, 2]) / 2, data[, 2] - data[, 1]), lwd = 2, col = "blue")

## Mean, variance relationship?
par(mfrow = c(2, 2))
plot(grp_means[, 1], grp_sds[, 1], xlab = "means_1", ylab = "sds_1", pch = 20, 
  col = "grey")
lines(lowess(grp_means[, 1], grp_sds[, 1]), lwd = 2, col = "blue")
plot(grp_means[, 2], grp_sds[, 2], xlab = "means_2", ylab = "sds_2", pch = 20, 
  col = "grey")
lines(lowess(grp_means[, 2], grp_sds[, 2]), lwd = 2, col = "blue")
plot(grp_means[, 3], grp_sds[, 3], xlab = "means_3", ylab = "sds_3", pch = 20, 
  col = "grey")
lines(lowess(grp_means[, 3], grp_sds[, 3]), lwd = 2, col = "blue")
plot(grp_means[, 4], grp_sds[, 4], xlab = "means_4", ylab = "sds_4", pch = 20, 
  col = "grey")
lines(lowess(grp_means[, 4], grp_sds[, 4]), lwd = 2, col = "blue")

#### Clustering.                

## Clustering on samples.
KM = kmeans(t(data), 4)
table(id, KM$cluster)

par(mfrow = c(1, 1))
HC = hclust(dist(t(data)), method = "complete")
plot(HC)
plot(HC$height)

#### SVD.                       

data_c = t(scale(t(data), center = T, scale = F))
ss = svd(data_c)
dd = ss$d ^ 2 / sum(ss$d ^ 2)
round(dd, 2)

oo = order(id)
plot(1:83, ss$v[oo, 1], type = "l", col = "green", ylim = c(-0.4, 0.4), lwd = 2)
lines(1:83, ss$v[oo, 2], col = "lightblue")
lines(1:83, ss$v[oo, 3], col = "pink")
abline(0, 0, lty = 2, col = "black")
for(i in cumsum(n_k)[1:3])
  lines(c(i, i), c(-1, 1))
  
####
#### Differential expression analysis, looking for differences between the four tumor 
#### classes.
####

p_vals = rep(NA, m)
for(i in 1:m) {
  d_i = as.numeric(data[i, ])
  fit = summary(lm(d_i ~ 1 + factor(id)))$fstat
  p_vals[i] = 1 - pf(fit[1], fit[2], fit[3])
}

## Histogram of p-values.
par(mfrow = c(1, 1))
hh = hist(p_vals, prob = T, main = "P-values", xlab = "")
abline(1, 0, lty = 2)

## Estimate FDR and q-values.
pi_0_hat = mean(hh$density[8:10])
library(qvalue)
qq = qvalue(p_vals)
qq$pi0
pi_0_hat

FDR_hat = rep(NA, m)
for(i in 1:m) {
  FDR_hat[i] = (pi_0_hat * m * p_vals[i]) / sum(p_vals <= p_vals[i])
}

q_vals = rep(NA, m)
for(i in 1:m) {
  q_vals[i] = min(FDR_hat[p_vals >= p_vals[i]])
}

## P-values vs. q-values.
plot(p_vals[order(p_vals)], q_vals[order(p_vals)], type = "s", 
  xlab = "P-values", ylab = "Q-values")
abline(0, 1, lty = 2)

## Number of genes selected at each q-value threshold.
S_hat = rep(NA, m)
for(i in 1:m) {
  S_hat[i] = sum(q_vals <= q_vals[i])
}

oo = order(q_vals)
plot(q_vals[oo], S_hat[oo], type = "s", xlab = "q-value", ylab = 
  "number of features selected")

## Can obtain 1000 features with an FDR < 0.01.  
which.min(abs(S_hat - 1000))
S_hat[2123]
q_vals[2123]
  
######
###### Classification.
######

## Set aside roughly 2 / 3 of the samples for training. Pick these randomly, within each 
## tumor type.
set.seed(101)

tr_ii <- NULL
for(i in 1:4) {
  jj <- (1:n)[as.numeric(id) == i]
  n_i <- length(jj)
  n_tr <- floor(n_i * 2 / 3)
  
  tr_ii <- c(tr_ii, sample(jj, n_tr, replace = FALSE))
}

Y_df <- data.frame("GRP" = factor(id), t(data))
Y_train_df <- Y_df[tr_ii, ]
Y_test_df <- Y_df[-tr_ii, ]

####
#### Naive Bayes is 86% accurate.
####

library(e1071)

fit_NB <- naiveBayes(GRP ~ ., data = Y_train_df)
pred_NB <- predict(fit_NB, newdata = Y_test_df)

table(pred_NB, Y_test_df$GRP)
mean(pred_NB == Y_test_df$GRP)

####
#### KNN. Choice of K doesn't seem to matter much, with the accuracies apparently bobbing 
#### randomly about 87% accuracy or so.
####

library(class)

pred_KNN <- rep(NA, 10)
for(k in 1:10) {
  pred_KNN_k <- knn(Y_train_df[, -1], Y_test_df[, -1], Y_train_df$GRP, k = k)
  pred_KNN[k] <- mean(pred_KNN_k == Y_test_df$GRP)
}

plot(1:10, pred_KNN, xlab = "K", ylab = "Accuracy", xaxt = "n", type = "l", lwd = 2)
axis(1, at = 1:10)

####
#### Lasso is nearly 100% accurate, using about 30 features.
####

library(glmnet)

set.seed(101)

cv_lasso <- cv.glmnet(as.matrix(Y_df[, -1]), Y_df$GRP, alpha = 1, family = "multinomial", 
  type.measure = "class")

plot(cv_lasso)

grid <- cv_lasso$lambda
fit_lasso <- glmnet(as.matrix(Y_df[, -1]), Y_df$GRP, alpha = 1, family = "multinomial", 
  type.multinomial = "grouped", lambda = grid)

par(mfrow = c(2, 2))
plot(fit_lasso, xvar = "lambda")

lambda_min <- cv_lasso$lambda.min
cf_lasso <- cbind(coef(fit_lasso)[[1]][-1, match(lambda_min, fit_lasso$lambda)], 
  coef(fit_lasso)[[2]][-1, match(lambda_min, fit_lasso$lambda)], 
  coef(fit_lasso)[[3]][-1, match(lambda_min, fit_lasso$lambda)], 
  coef(fit_lasso)[[4]][-1, match(lambda_min, fit_lasso$lambda)])
  
features_lasso <- (1:m)[rowSums(cf_lasso) != 0]

####
#### Feature selection, then KNN. As lasso suggests, high accuracy can be obtained with 
#### a small fraction of the total number of features. Let's try carrying out formal 
#### feature selection via feature-specific statistical tests. We'll just do K = 5 for 
#### KNN.
####

## Let's do five-fold cross-validation.  Our sample sizes in each group don't divide 
## nicely into five groups, so we'll use these sample sizes for testing in each fold.  
fold_1_n = c(2, 6, 4, 5)
fold_2_n = c(2, 6, 4, 5)
fold_3_n = c(2, 6, 4, 4)
fold_4_n = c(2, 6, 3, 5)
fold_5_n = c(3, 5, 3, 6)

fold_n = rbind(fold_1_n, fold_2_n, fold_3_n, fold_4_n, fold_5_n)

## Now randomly select the above numbers from each class.
class_1_samples = sample((1:n)[id == 1])
class_2_samples = sample((1:n)[id == 2])
class_3_samples = sample((1:n)[id == 3])
class_4_samples = sample((1:n)[id == 4])

fold_1_samples = c(class_1_samples[1:2], class_2_samples[1:6], class_3_samples[1:4], 
  class_4_samples[1:5])
fold_2_samples = c(class_1_samples[3:4], class_2_samples[7:12], class_3_samples[5:8], 
  class_4_samples[6:10])
fold_3_samples = c(class_1_samples[5:6], class_2_samples[13:18], class_3_samples[9:12], 
  class_4_samples[11:14])
fold_4_samples = c(class_1_samples[7:8], class_2_samples[19:24], class_3_samples[13:15], 
  class_4_samples[15:19])
fold_5_samples = c(class_1_samples[9:11], class_2_samples[25:29], class_3_samples[16:18], 
  class_4_samples[20:25])
  
fold_samples = list(fold_1_samples, fold_2_samples, fold_3_samples, fold_4_samples, 
  fold_5_samples)
  
## Since we have 4 classes, let's use an F-statistic, which is a measure of between-class 
## variability to within-class variability.  Here's a function to compute F-statistics 
## for a given set of intensities.
F_stat = function(Y, ID) {
  n_k = as.numeric(table(ID))
  n = sum(n_k)
  k = length(n_k)
  
  X = model.matrix(~ factor(ID) - 1)
  Y_bars = t(t(Y %*% X) / n_k)
  Y_bars_overall = drop(Y %*% rep(1 / n, n))
  
  F_num = drop(t(t((Y_bars - Y_bars_overall) ^ 2) * n_k) %*% rep(1, k))
  F_den = drop((Y - Y_bars %*% t(X)) ^ 2 %*% rep(1, n))
  
  return(F_num / F_den)
}

## Now do CV. For each CV fold, we'll build and assess the accuracy of feature set sizes 
## of 2, 3, ..., 30, 40, 50, ..., 2300. We need to do the feature selection within each 
## CV fold to avoid model selection bias.
f_sizes = c(2:30, seq(40, 2300, by = 10))
acc_f_b = matrix(NA, nrow = length(f_sizes), ncol = 5)
for(b in 1:5) {
  cat(".")

  Y_test = data[, fold_samples[[b]]] 
  Y_train = data[, -fold_samples[[b]]]
  id_test = id[fold_samples[[b]]]
  id_train = id[-fold_samples[[b]]]
  n_test = length(id_test)
  n_train = length(id_train)
  
  ## Feature selection.
  FF = F_stat(Y_train, id_train)
  oo = order(FF, decreasing = T)
  
  for(f in 1:length(f_sizes)) {
    feature_set = oo[1:f_sizes[f]]
    
    y_train = Y_train[feature_set, ]
    y_test = Y_test[feature_set, ]
    
    pred_KNN = knn(t(y_train), t(y_test), id_train, k = 5)
    acc_f_b[f, b] = mean(pred_KNN == id_test)
  }
}

acc_f = rowMeans(acc_f_b)

## Here's the big picture. We obtain basically 100% accuracy with about 30 features. Note 
## that we are actually *worse* off when including more features, beyond a certain point. 
## This is common, as eventually we're just adding noise to the classifier.
par(mfrow = c(1, 1))
plot(f_sizes, acc_f, xlab = "Feature set size", ylab = "CV-based accuracy", lwd = 2, 
  type = "s")

## Here's the same picture, but zoomed in on just the lower subset sizes (2, 3, ..., 30).
plot(f_sizes[1:29], acc_f[1:29], xlab = "Feature set size", ylab = "CV-based accuracy", 
  lwd = 2, type = "s")

## Boxplots for the top 4 features. Two features that characterize class 4 and 2 features 
## that characterize class 2. An optimal subset might include just one of each, if the 
## others provide redundant information.
FF = F_stat(data, id)
oo = order(FF, decreasing = T)

par(mfrow = c(2, 2))
boxplot(data[oo[1], id == 1], data[oo[1], id == 2], data[oo[1], id == 3], 
  data[oo[1], id == 4])
boxplot(data[oo[2], id == 1], data[oo[2], id == 2], data[oo[2], id == 3], 
  data[oo[2], id == 4])
boxplot(data[oo[3], id == 1], data[oo[3], id == 2], data[oo[3], id == 3], 
  data[oo[3], id == 4])
boxplot(data[oo[4], id == 1], data[oo[4], id == 2], data[oo[4], id == 3], 
  data[oo[4], id == 4])
  
features_KNN_fs <- (1:m)[oo[1:30]]

####
#### Classification tree.
####

library(tree)

## Because of the relatively low sample sizes in the different tumor types, the 'tree' 
## function ends up only using 3 of the features, achieving a training accuracy of about 
## 96%.
fit_tree <- tree(GRP ~ ., data = Y_train_df)
summary(fit_tree)
plot(fit_tree)
text(fit_tree)

## The test accuracy is only about 76%, making it inferior to all the methods considered 
## above (in this example, at least). Note, however, the very simple interpretation of 
## the classifier. In terms of defining a diagnostic tool, it would be straightforward to 
## talk in terms like "if this gene has expression below this threshold, then if that 
## gene has expression above that threshold, then you are tumor type X."
pred_tree <- predict(fit_tree, newdata = Y_test_df, type = "class")
table(pred_tree, Y_test_df$GRP)

## While the tree is already small, we could try pruning it further. According to CV, the 
## estimated test error is 87% for the full tree, or for the tree with one fewer terminal 
## node.
set.seed(101)
cv_tree <- cv.tree(fit_tree, FUN = prune.misclass)
cv_tree

plot(cv_tree$size, cv_tree$dev, type = "b")

## Here is the tree pruned to 4 leaves. Basically the same as the tree from above, just 
## with the second split around X1 removed. The test accuracy is actually the same as 
## we saw above, about 76%. Note that CV overestimated this.
pruned_tree <- prune.misclass(fit_tree, best = 4)
plot(pruned_tree)
text(pruned_tree)

pred_pruned_tree <- predict(pruned_tree, newdata = Y_test_df, type = "class")
table(pred_pruned_tree, Y_test_df$GRP)

features_tree <- c(1, 107, 1319)

## Overlap between features selected by tree, KNN + feature selection, lasso.
feature_overlap <- matrix(NA, nrow = 3, ncol = 3)
rownames(feature_overlap) <- colnames(feature_overlap) <- c("Tree", "KNN_fs", "Lasso")
feature_overlap[1, 2] <- mean(features_tree %in% features_KNN_fs)
feature_overlap[1, 3] <- mean(features_tree %in% features_lasso)
feature_overlap[2, 3] <- mean(features_KNN_fs %in% features_lasso)

####
#### Random forest.
####

library(randomForest)

set.seed(101)

## The bagged tree achieves an accuracy of nearly 100%. 
fit_rf <- randomForest(GRP ~ ., data = Y_train_df, mtry = 2, importance = TRUE)
fit_rf

pred_rf <- predict(fit_rf, newdata = Y_test_df)
table(pred_rf, Y_test_df$GRP)

imp <- importance(fit_rf)
imp[order(imp[, "MeanDecreaseGini"], decreasing = TRUE), ][1:20, ]
varImpPlot(fit_rf)

####
#### SVM.
####

library(e1071)

set.seed(101)

## The 'tune' function takes a long time to run, if we want to evaluate cost together 
## with polynomial degree or cost together with radial kernel gamma parameter. After 
## running for a while on my own, it appears that a support vector (linear) classifier is 
## the most accurate of the polynomial and radial kernel classifiers. Accuracy does not 
## appear to change much with 'cost', so I'll just use 'cost' = 1.
fit_svm <- svm(GRP ~ ., data = Y_train_df, kernel = "linear", cost = 1)
pred_svm <- predict(fit_svm, newdata = Y_test_df)
table(pred_svm, Y_test_df$GRP)

####
#### Neural networks. I am not going to cover them. See 'nnet' package.
####





