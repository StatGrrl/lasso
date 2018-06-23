# BOOTSTRAP STANDARD ERRORS
boot.se <- function(X, B, FUN, ...) {
  # INPUT:
  # x........vector of random variables, or
  #          matrix/data.frame where each column is a random variable
  # B........number of bootstrap replications
  # FUN......function to calculate statistic of x
  #          with arguments ...
  #
  # OUTPUT:
  # returns bootstrap standard error of statistic 
  # for each random variable in X
  #
  # example: calculate standard error of the mean of vector x
  #          boot.se(x, B=200, mean)
  #          sd(x)/sqrt(length(x))
  #
  #          add arguments to mean function with ...
  #          boot.se(x, B=200, mean, na.rm=TRUE)
  
  FUN <- match.fun(FUN)
  if (is.null(dim(X))) X <- as.matrix(X)
  c <- ncol(X)
  boot.fun <- matrix(0,B,c)
  for (i in 1:B) {
    boot.sample <- apply(X, 2, sample, replace=TRUE)
    boot.fun[i,] <- apply(boot.sample, 2, FUN, ...)
  }
  apply(boot.fun, 2, sd)
}
