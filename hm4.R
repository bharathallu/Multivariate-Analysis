###2

dt<- matrix(c(5,4,1,-2,-1,1,3,1),nrow=4,byrow=T)
init <- matrix(c(2,2.5,2,-0.5),nrow=2, byrow=T)
k1<- kmeans(dt, init) 








c<- matrix(c(0,3,2,5,3,0,4,1,2,4,0,7,5,1,7,0), nrow = 4,byrow=T)
d<- as.dist(c)
hc_single<- hclust(d,method="single")
hc_single$merge
hc_single$height
plot(hc_single)

hc_comp<- hclust(d,method="complete")
hc_comp$merge
hc_comp$height
plot(hc_comp)

hc_avg<- hclust(d,method="average")
hc_avg$merge
hc_avg$height
plot(hc_avg)


dta<- read.csv(file="C:/STAT 636/New folder/cereal.csv",header = T,sep=",",row.names=1)
X<- as.matrix(dta)

# Hierarchical clustering using complete linkage 

hc<- hclust(d=dist(X),method="complete")
# dendogram plot
par(mfrow=c(1,1))
plot(hc)

# height plot to determine number of clusters.
plot(hc$height)
hc_cluster<- cutree(hc,k=3)

# composition of clusters.
par(mfrow=c(2,4))
with(dta,boxplot(Calories~hc_cluster,xlab='cluster',ylab="calories"))
with(dta,boxplot(Protein~hc_cluster,xlab='cluster',ylab="Protein"))
with(dta,boxplot(Fat~hc_cluster,xlab='cluster',ylab="Fat"))
with(dta,boxplot(Sodium~hc_cluster,xlab='cluster',ylab="Sodium"))
with(dta,boxplot(Fiber~hc_cluster,xlab='cluster',ylab="Fiber"))
with(dta,boxplot(Carbohydrates~hc_cluster,xlab='cluster',ylab="Carbohydrates"))
with(dta,boxplot(Sugar~hc_cluster,xlab='cluster',ylab="Sugar"))
with(dta,boxplot(Potassium~hc_cluster,xlab='cluster',ylab="Potassium"))

# k means clustering and plot

set.seed(101)
kc<- kmeans(X,3)
plot(X,col= kc$cluster,pch=19)
points(kc$centers,col=1:3,pch=8)
table(kc$cluster)


library(cluster)
library(fpc)

plotcluster(X,kc$cluster,pch=18)

# plotting first two principal components

pca<- princomp(X)
pca$loadings
pca1<- X %*% pca$loadings[,1]
pca2<- X %*% pca$loadings[,2]
par(mfrow=c(1,1))
plot(pca1,pca2,pch=19,col=as.factor(kc$cluster),main="First Two PCAs")
legend("topright",c("cluster1","cluster2","cluster3"),col=c(3,1,2),pch=19)
summary(pca)


dta[kc$cluster==2,]  
dta[kc$cluster==3,] 


# comparision of k means with Hierarchical

table(hc_cluster)
table(kc$cluster)

par(mfrow=c(1,1))
plotcluster(X,kc$cluster,pch=18,main="kmeans")
legend("topright",c("cluster1","cluster2","cluster3"),col=c(2,1,3),pch=19,cex=.6)
#text(X,kc$cluster, labels=row.names(dta),cex =0.3)
plotcluster(X,hc_cluster,pch=18,main="hierarchical")
legend("topright",c("cluster1","cluster2","cluster3"),col=c(3,1,2),pch=19)
mat<- cbind(kc$cluster,hc_cluster)


####4

dta1 <- as.matrix(read.table("C:/STAT 636/New folder/T12-8.DAT", header = FALSE))
pot_type <- c("A", "B", "C", "D")
pot_site <- paste("P", 0:6, sep = "_")
rownames(dta1) <- pot_site
colnames(dta1) <- pot_type

##multi dimensional scaling
D<- dist(as.matrix(dta1))
dta1_2<- cmdscale(d=D,k=2)
plot(dta1_2,ylim=c(-40,30))
text(dta1_2,labels=row.names(dta1),pos=3)

###bi plot

par(mfrow=c(1,1))
dta1_st<- scale(dta1,center = T, scale = T)
a<- princomp(dta1)
a1<- princomp(dta1_st)
biplot(a,expand=2.5,ylim=c(-2,2),xlim=c(-2,0.5))
biplot(a1,expand=2.5,ylim=c(-1.5,2),xlim=c(-2,1.5))
