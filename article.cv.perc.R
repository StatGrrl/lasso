cv.perc <- function(X, y, K=10, M=100, theta=0.95) {
  # INPUT:
  # X...........predictor matrix
  # y...........response vector
  # K...........number of CV folds
  # M...........number of times to repeat CV
  # theta.......percentile of CV estimates
  
  lambda <- rep(0,M)
  for (i in 1:M)
    lambda[i] <- cv.glmnet(X, y, nfolds=K, standardize=F, intercept=F)$lambda.min
  quantile(lambda,theta)
}