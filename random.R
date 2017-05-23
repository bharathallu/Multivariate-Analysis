cv <- function(x) sd(x)/mean(x)
cv(1:10)
cv<-function(x) var(x)/sd(x)
cv(1:10)
?lapply
gcd <- function(a,b) {
   if (b == 0) return(a)
  else return(gcd(b, a %% b))
  }
scores<-data.frame()
scores<-edit(scores)
scores
 points <- data.frame(
  label=c("Low", "Mid", "High"),
   lbound=c( 0, 0.67, 1.64),
   ubound=c(0.674, 1.64, 2.33)
  )
 ?read.table
 v <- c(10, 20, 30)
 names(v) <- c("Moe", "Larry", "Curly")
 data.frame(x=1:3, y=c("NY", "CA", "IL"))
 f <- factor(c("Win","Win","Lose","Tie","Win","Lose"))
 comb <- stack(list(fresh=freshmen, soph=sophomores, jrs=juniors))
 ?cat
 iter <- stats::rpois(1, lambda = 10)
 ## print an informative message
 cat("iteration = ", iter <- iter + 1, "\n")
  years <- list(Kennedy=1960, Johnson=1964, Carter=1976, Clinton=1994)
  values <- pnorm(-2:2)
   names <- c("far.left", "left", "mid", "right", "far.right")
   for (i in names(lst)) cat("The", i, "limit is", lst[[i]], "\n")
   mat=matrix(c(1,2,3,1,2,3,1,2,3,4,2,2,2,2,5,6),4,4)
 
   apply(mat,c(1,2),mean)
   