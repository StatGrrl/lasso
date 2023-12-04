pe.path <- function(method, data) {
  # INPUT:
  # method......lasso (default) - requires lars
  #             ridge - requires MASS
  #             forward / backward / exhaustive - requires leaps
  # data........simulation data needed:
  #             list with beta0, beta, X.test, y.test, X.train, y.train
  
  source("estimation.R")
    
  # DATA
  beta0 <- data$beta0
  int <- beta0!=0
  beta <- data$beta
  p <- length(beta)
  n.train <- length(data$y.train[[1]])
  
  # TEST DATA
  X.test <- data$X.test
  y.test <- data$y.test
  y.test.mu <- beta0 + X.test %*% beta
  n.test <- length(y.test)
  
  # OUTPUT
  iter <- length(data$y.train)
  if (method == "lasso") comp <- seq(0, 1, by=0.01)
  if (method == "ridge") comp <- seq(0, 50, by=0.2)
  if (method %in% c("forward","backward")) comp <- 1:p
  n.comp <- length(comp)
  path <- list()
  cov.train <- matrix(0, iter, n.comp)
  df.aprx <- matrix(0, iter, n.comp)
  yhat.test <- rep(list(matrix(0, n.test, iter)), n.comp)

  for (i in 1:iter) {
    
    # TRAINING DATA
    X <- data$X.train[[i]]
    y <- data$y.train[[i]]
    y.mu <- beta0 + X %*% beta
    
    # ESTIMATION
    est <- estimate(X, y, int=int, method=method, comp=comp)
    
    # PREDICTION
    if (int) X.fit <- cbind(1,X) else X.fit <- X
    fit <- apply(est$coef, 1, function(a) X.fit %*% a)
    if (int) X.pred <- cbind(1,X.test) else X.pred <- X.test
    pred <- apply(est$coef, 1, function(a) X.pred %*% a)
    path[[i]] <- data.frame(comp=comp)
    path[[i]]$df <- est$df
    path[[i]]$Nonzero <- rowSums(est$coef!=0)
    path[[i]]$Train <- apply(fit, 2, function(a) sum((y - a)^2)/n.train)
    path[[i]]$Test <- apply(pred, 2, function(a) sum((y.test - a)^2)/n.test)

    for (j in 1:n.comp) {
      cov.train[i,j] <- sum((y - y.mu) * fit[,j])
      df.aprx[i,j] <-  path[[i]]$Nonzero[j]
      yhat.test[[j]][,i] <- pred[,j]
    }
  }
  
  # bias, var, mse, df
  mse <- data.frame(comp, matrix(0,n.comp,3))
  names(mse) <- c("comp", "df", "Sq Bias", "Variance")
  mse[,"df"] <- colMeans(cov.train)/data$sigma^2
  for (i in 1:n.comp) {
    yhat.mean <- apply(yhat.test[[i]], 1, mean)
    yhat.var <- apply(yhat.test[[i]], 1, var)
    mse[i,"Sq Bias"] <- sum((yhat.mean - y.test.mu)^2)/n.test
    mse[i,"Variance"] <- sum(yhat.var)/n.test
  }
  mse$MSE <- mse[,"Sq Bias"] + mse[,"Variance"]
  mse[,"df.aprx"] <- colMeans(df.aprx)
  
  list(path=path, mse=mse)
}
