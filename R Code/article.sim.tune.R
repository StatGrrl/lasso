sim.tune <- function(data) {
  # INPUT:
  # data........simulation data from sim.data
  
  ## DATA
  beta <- data$beta
  p <- length(beta)
  q <- sum(beta != 0)
  A <- which(beta != 0)
  
  # Validation data
  X.valid <- data$X.valid
  y.valid <- data$y.valid
  n.valid <- data$n.valid
  
  iter <- length(data$y.train)

  ## OUTPUT
  select.names <- c("Validation.Set", "5.fold.CV", "10.fold.CV", "LOOCV",
                    "GCV", "Cp", "AIC", "5.fold.CV.1SE", "10.fold.CV.1SE",
                    "Percentile.CV", "Kappa", "PASS", "BIC", "Modified.BIC")
  select.num <- length(select.names)
  select.index <- rep(0, select.num)
  perf.names <- c("n", "lambda", paste0("beta.X",1:p), "df", "MSE", 
                  "Corr.Nonzero", "Inc.Nonzero", "Corr.Zero", "Inc.Zero", "PCS", "PIS")
  perf.num <- length(perf.names)
  perf <- data.frame(matrix(0, nrow=iter, ncol=perf.num))
  names(perf) <- perf.names
  PIP <- rep(0, iter)
  best <- rep(list(perf), select.num)
  names(best) <- select.names
  
  for (i in 1:iter) {
    # Taining data
    X <- data$X.train[[i]]
    y <- data$y.train[[i]]
    n <- data$n.train
    
    # Standardize data
    ymean <- mean(y)
    y <- y-mean(y)  
    xmean <- colMeans(X)
    xnorm <- sqrt(n-1)*apply(X,2,sd)
    X <- scale(X, center=xmean, scale=xnorm)
    
    ## LASSO ESTIMATION
    require(glmnet)
    model <- glmnet(X, y, family="gaussian", alpha=1, standardize=F, intercept=F)
    lambda <- model$lambda
    betahat <- as.matrix(t(coef(model)))[,-1]
    betahat <- scale(betahat, center=FALSE, scale=xnorm)
    path <- data.frame(n=n, lambda=lambda*n, beta=betahat)
    
    ## PERFORMANCE MEASURES
    
    df <- rowSums(betahat!=0) # degrees of freedom
    path$df <- df
    diff <- apply(betahat, 1, function(a) a - beta)
    path$MSE <- apply(diff, 2, function(a) t(a) %*% data$X.sigma %*% a) # mean squared error
    path$Corr.Nonzero <- apply(betahat, 1, function(a) sum(a[A]!=0)) # coef correctly nonzero
    path$Inc.Nonzero <- apply(betahat, 1, function(a) sum(a[-A]!=0)) # coef incorrectly nonzero
    path$Corr.Zero <- apply(betahat, 1, function(a) sum(a[-A]==0)) # coef correctly zero
    path$Inc.Zero <- apply(betahat, 1, function(a) sum(a[A]==0)) # coef incorrectly zero
    path$PCS <- as.integer(path$Inc.Nonzero + path$Inc.Zero == 0) # probability correct subset
    path$PIS <- as.integer(path$Corr.Nonzero==q) # probability including subset
    PIP[i] <- as.integer(any(path$PCS==1)) # probability in path
    
    ## MODEL SELECTION
    
    # Validation set
    pred.v <- apply(betahat, 1, function(a) X.valid %*% a)
    valid <- apply(pred.v, 2, function(a) sum((y.valid - a)^2)/n.valid)
    
    # Cross-validation
    cv5 <- cv.glmnet(X, y, nfolds=5, standardize=F, intercept=F)
    cv10 <- cv.glmnet(X, y, nfolds=10, standardize=F, intercept=F)
    loocv <- cv.glmnet(X, y, nfolds=n, standardize=F, intercept=F)
    
    # Percentile CV
    source("article.cv.perc.R")
    perc <- cv.perc(X, y, K=5)
    
    # Kappa and PASS
    require(pass)
    split <- pass(cbind(X,y), lambda.grid=lambda*n)
    
    # Information Crietria
    X <- data$X.train[[i]]
    y <- data$y.train[[i]]
    pred <- apply(betahat, 1, function(a) X %*% a)
    RSS <- apply(pred, 2, function(a) sum((y - a)^2))
    ols <- lm(y~X-1)
    sigmahat2 <- summary(ols)$sigma^2
    
    gcv <- (1/n)*RSS/((1-df/n)^2)
    cp <- (1/n)*(RSS + 2*df*sigmahat2)
    aic <- (1/n)*(RSS/sigmahat2 + 2*df)
    bic <- (1/n)*(RSS/sigmahat2 + log(n)*df)
    mbic <- log(RSS/n) + df*(log(n)/n)*(sqrt(n)/p)

    select.index[1] <- which.min(valid)
    select.index[2] <- which(lambda==cv5$lambda.min)
    select.index[3] <- which(lambda==cv10$lambda.min)
    select.index[4] <- which(lambda==loocv$lambda.min)
    select.index[5] <- which.min(gcv)
    select.index[6] <- which.min(cp)
    select.index[7] <- which.min(aic)
    select.index[8] <- which(lambda==cv5$lambda.1se)
    select.index[9] <- which(lambda==cv10$lambda.1se)
    select.index[10] <- which(lambda==sort(lambda)[findInterval(perc, sort(lambda))])
    select.index[11] <- which(lambda*n==split$lambda.kappa)
    select.index[12] <- which(lambda*n==split$lambda.pass)
    select.index[13] <- which.min(bic)
    select.index[14] <- which.min(mbic)
        
    for (j in 1:select.num)
      best[[j]][i,] <- path[select.index[j],]
  }
  
  # Average Performance Measures and Median MSE with SE
  bestave <- data.frame(matrix(0, select.num, perf.num+3))
  names(bestave) <- c("Method", perf.names, "Median.MSE","SE.Median.MSE")
  source("article.boot.se.R")
  for (i in 1:select.num) {
    bestave[i,1:perf.num+1] <- colMeans(best[[i]])
    bestave[i,"Median.MSE"] <- median(best[[i]]$MSE)
    bestave[i,"SE.Median.MSE"] <- boot.se(best[[i]]$MSE, B=200, median)
    bestave[i,1] <- select.names[i]
  }
  PIP <- mean(PIP)
  out <- list(best=best, bestave=bestave, PIP=PIP)
  out
}
