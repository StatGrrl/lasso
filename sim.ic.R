sim.ic <- function(method, data) {
  # INPUT:
  # method......lasso - requires lars
  #             ridge - requires MASS
  #             forward / backward / exhaustive - requires leaps
  #             ols
  # data........all simulation data needed:
  #             beta0, beta, rho, X.sigma, sigma, X.test, y.test, 
  #             X.train, y.train, X.valid, y.valid
  #
  # OUTPUT:
  # The true prediction error and estimatesvar selection performance measures
  #
  # path........list of size iter, entire path at each iteration
  # pathave.....averages over path, incl mse and df estimates
  # best........list containing the best model selected at each iteration
  #             by each of the prediction estimates
  # bestave.....average results, median mse and bootstrap standard errors

  source("estimation.R")
  source("resample.R")
    
  # DATA
  beta0 <- data$beta0
  int <- beta0!=0
  beta <- data$beta
  p <- length(beta)
  d <- sum(beta != 0)
  d.index <- which(beta != 0)
  
  # TEST DATA
  X.test <- data$X.test
  y.test <- data$y.test
  y.mu <- beta0 + X.test %*% beta
  n.test <- length(y.test)
  
  # VALIDATION DATA
  X.valid <- data$X.valid
  y.valid <- data$y.valid
  n.valid <- length(y.valid)
  
  # OUTPUT
  iter <- length(data$y.train) # number iterations

  # all fields collected
  beta.names <- paste("b",1:p, sep=".")
  if (int) beta.names <- c("b.0", beta.names)
  beta.prob <- paste0("prob.",beta.names)
  out.names <- c("Method", "IC.Method", "sigma", "rho", "comp", "df", 
                 beta.names, beta.prob,
                 "MSE", "PE", "Train", "Test", 
                 "Nonzero", "Corr.Nonzero", "Inc.Nonzero", 
                 "Corr.Zero", "Inc.Zero", "Corr.Subset", 
                 "Incl.Corr.Subset", "In.Path", 
                 "IC",  "iter")
  var.fixed <- 3 # Method, sigma, rho are fixed
  n.out <- length(out.names)
  beta.pos <- which(out.names %in% beta.names)
  prob.pos <- which(out.names %in% beta.prob)
  
  # model complexity index
  if (method == "ols") comp <- c(p,d)
  if (method == "lasso") {
    comp <- seq(0, 1, by=0.01)
    mode <- "fraction"
  }
  if (method == "ridge") {
    comp <- seq(0, 50, by=0.2)
    mode <- "lambda"
  }
  if (method %in% c("forward","backward")) {
    comp <- 1:p
    mode <- ""
  }
  n.comp <- length(comp)
  
  # collect covariance of y and yhat (training yhat)
  cov.train <- matrix(0, iter, n.comp)
  
  # collect predictions on test set
  yhat <- rep(list(matrix(0, n.test, iter)), n.comp)

  # information criteria and cv
  if (method == "ols")
    ic.names <- c("Full","Oracle") else 
      ic.names <- c("Valid", "5.fold.CV", "10.fold.CV", "LOOCV", "GCV", "Cp", 
                    "BIC", "AIC", "mBIC", "kappa", "CV.Percentile", 
                    "5.fold.CV.1SE", "10.fold.CV.1SE" )
  ic.getmin <- 9 # get minimum position for IC (up till mBIC)
  ic.pe <- 7 # IC estimates of PE (up till BIC)
  n.ic <- length(ic.names)
  
  #collect path
  path <- list()
  
  # position of best model in path
  best.select <- as.data.frame(matrix(0,iter,n.ic))
  names(best.select) <- c(ic.names)
  
  # best models
  tmp <- as.data.frame(matrix(0, nrow=iter, ncol=n.out))
  colnames(tmp) <- out.names
  best <- rep(list(tmp), n.ic)
  names(best) <- ic.names
  rm(tmp)

  for (i in 1:iter) {
    
    # TRAINING DATA
    X <- data$X.train[[i]]
    y <- data$y.train[[i]]
    n.train <- length(y)
    mu <- beta0 + X %*% beta
    
    
    # ESTIMATION
    ols <- estimate(X, y, int=int, method="ols")
    sigma2hat <- ols$rss[[1]]/(n.train - p -int)
    if (method == "ols") {
      oracle <- estimate(X, y, int=int, method="oracle", true=d.index)
      est <- list(coef=rbind(ols$coef, oracle$coef), df=comp+int, 
                  fit=list(cbind(ols$fit[[1]],oracle$fit[[1]])), 
                           rss=list(c(ols$rss[[1]],oracle$rss[[1]])))
    } else {
      est <- estimate(X, y, int=int, method=method, comp=comp, mode=mode)
    }
    path[[i]] <- data.frame(Method=toupper(method), sigma=data$sigma,
                            rho=data$rho, comp=comp, df=est$df, 
                            b=est$coef, prob.b=est$coef!=0)
    names(path[[i]])[beta.pos-1] <- beta.names
    names(path[[i]])[prob.pos-1] <- beta.prob

    if (int) {
      betahat <- est$coef[,-1]
      X.pred <- cbind(1,X.test)
      X.v <- cbind(1,X.valid)
    } else {
      betahat <- est$coef
      X.pred <- X.test
      X.v <- X.valid
    }
    
    # PREDICTION PERFORMANCE
    diff <- apply(betahat, 1, function(a) a - beta)
    pred <- apply(est$coef, 1, function(a) X.pred %*% a)
    path[[i]]$MSE <- apply(diff, 2, function(a) t(a) %*% data$X.sigma %*% a)
    path[[i]]$PE <- path[[i]]$MSE + data$sigma^2
    path[[i]]$Train <- est$rss[[1]]/n.train
    path[[i]]$Test <- apply(pred, 2, function(a) sum((y.test - a)^2)/n.test)
    for (j in 1:n.comp) {
      cov.train[i,j] <- sum((y - mu) * est$fit[[1]][,j])
      yhat[[j]][,i] <- pred[,j]
    }
    
    # VARIABLE SELECTION PERFORMANCE
    path[[i]]$"Nonzero" <- rowSums(est$coef!=0)
    path[[i]]$"Corr.Nonzero" <- apply(est$coef, 1, function(a) sum(a[d.index]!=0))
    path[[i]]$"Inc.Nonzero" <- apply(est$coef, 1, function(a) sum(a[-d.index]!=0))
    path[[i]]$"Corr.Zero" <- apply(est$coef, 1, function(a) sum(a[-d.index]==0))
    path[[i]]$"Inc.Zero" <- apply(est$coef, 1, function(a) sum(a[d.index]==0))
    path[[i]]$"Corr.Subset" <- as.integer(path[[i]]$"Inc.Nonzero" + 
                                            path[[i]]$"Inc.Zero" == 0)
    path[[i]]$"Incl.Corr.Subset" <- as.integer(path[[i]]$"Corr.Nonzero" == d)
    path[[i]]$"In.Path" <- as.integer(any(path[[i]]$"Corr.Subset"==1))
    
    if (method != "ols") {
      # VALIDATION SET, CV AND INFORMATION CRITERIA
      cv5 <- cv.kfold(X, y, int, method, comp, k=5, mode=mode)
      cv10 <- cv.kfold(X, y, int, method, comp, k=10, mode=mode)
      LOOCV <- cv.kfold(X, y, int, method, comp, k=n.train, mode=mode)
      cvperc <- cv.perc(X, y, int, method, comp, k=10, mode=mode, reps=20)
      kappa <- res.kappa(X, y, int, method, comp, mode=mode, reps=20)
      pred.v <- apply(est$coef, 1, function(a) X.v %*% a)
      path[[i]]$Valid <- apply(pred.v, 2, function(a) sum((y.valid - a)^2)/n.valid)
      path[[i]]$"5.fold.CV" <- cv5$cv
      path[[i]]$"10.fold.CV" <- cv10$cv
      path[[i]]$"LOOCV" <- LOOCV$cv
      path[[i]]$GCV <- (1/n.train)*est$rss[[1]]/(1-est$df/n.train)^2
      path[[i]]$Cp <- (est$rss[[1]] + 2*est$df*sigma2hat)/n.train
      path[[i]]$BIC <- est$rss[[1]]/n.train + log(n.train)*est$df*sigma2hat/n.train
      path[[i]]$AIC <- est$rss[[1]]/(n.train*sigma2hat) + 2*est$df/n.train
      path[[i]]$mBIC <- log(est$rss[[1]]/n.train) + 
        est$df*(log(n.train)/n.train)*(sqrt(n.train)/p)
      path[[i]]$"kappa" <- 0
      path[[i]]$"CV.Percentile" <- 0
      path[[i]]$"5.fold.CV.1SE" <- cv5$cv
      path[[i]]$"10.fold.CV 1SE" <- cv10$cv

      best.select[i, "kappa"] <- kappa$index
      best.select[i, "CV.Percentile"] <- cvperc$index
      best.select[i, "5.fold.CV.1SE"] <- cv5$index.1SE
      best.select[i, "10.fold.CV.1SE"] <- cv10$index.1SE
    }

    # MODEL SELECTION
    perf.end <- n.out-var.fixed
    for (j in 1:n.ic) {
      if (method == "ols") best.select[i,j] <- j else
        if (j <= ic.getmin) 
          best.select[i,j] <- which.min(path[[i]][,perf.end+j])
      which <- best.select[i,j]    
      best[[j]][i,"Method"] <- toupper(method)
      best[[j]][i,"IC.Method"] <- ic.names[j]
      best[[j]][i,3:4] <- path[[i]][which,2:3]
      best[[j]][i,5:(perf.end+1)] <- path[[i]][which,4:perf.end]
      if (method != "ols")
        best[[j]][i,"IC"] <- path[[i]][which, perf.end+j]
      best[[j]][i,"iter"] <- i
    }

    col.end <- if (method=="ols") perf.end else 
      ncol(path[[i]])-n.ic+ic.pe
    col.st <- var.fixed+1
    if (i==1) 
      pathave <- path[[i]][,1:col.end] else
        pathave[,-(1:col.st)] <- pathave[,-(1:col.st)] +
        path[[i]][,(col.st+1):col.end]
    path[[i]] <- cbind(path[[i]][,1:col.st], 
                       path[[i]][,c("PE","Train","Test")])
  }

  
  # path averages
  pathave[,-(1:col.st)] <- pathave[,-(1:col.st)]/iter
  pathave$"df.cov" <- colMeans(cov.train)/data$sigma^2
  pathave$"Sq.Bias" <- rep(0,n.comp)
  pathave$"Variance" <- rep(0,n.comp)
  for (i in 1:n.comp) {
    yhat.mean <- apply(yhat[[i]], 1, mean)
    yhat.var <- apply(yhat[[i]], 1, var)
    pathave[i,"Sq.Bias"] <- sum((yhat.mean - y.mu)^2)/n.test
    pathave[i,"Variance"] <- sum(yhat.var)/n.test
  }
  pathave$"MSE.Test" <- pathave[,"Sq.Bias"] + pathave[,"Variance"]

  for (i in 1:length(best)) {
    best.index <- best.select[,i]
    best[[i]] <- 
      cbind(best[[i]], pathave[best.index,
                        c("df.cov","Sq.Bias","Variance","MSE.Test")])
    rownames(best[[i]])=NULL
  }
  
  # average performance measures
  col.n <- n.out+4
  col.no.ave <- 1:(var.fixed+1)
  col.ave <- (1:col.n)[-col.no.ave]
  bestave <- as.data.frame(matrix(0, n.ic, col.n+5))
  names(bestave) <- c(names(best[[1]]), "Median.MSE","SE.Median.MSE",
                     "Sq.Bias.Betas", "Var.Betas", "MSE.Betas")
  for (i in 1:n.ic) {
    bestave[i,col.no.ave] <- best[[i]][1,col.no.ave]
    bestave[i,col.ave] <- colMeans(best[[i]][,col.ave])
    bestave[i,"Median.MSE"] <- median(best[[i]]$MSE)
    bestave[i,"SE.Median.MSE"] <- 
      boot.se(best[[i]]$MSE, B=200, median)
    bestave[i,"Sq.Bias.Betas"] <- 
      sum((bestave[i,beta.pos]-beta)^2)
    bestave[i,"Var.Betas"] <- 
      sum(apply(best[[i]][,beta.pos], 2, var))
    bestave[i,"MSE.Betas"] <- 
      bestave[i,"Sq.Bias.Betas"] + bestave[i,"Var.Betas"]
  }
  bestave <- bestave[,-n.out]

  results <- if (method == "ols") 
    list(bestave=bestave, best=best) else
      list(bestave=bestave, best=best, pathave=pathave, path=path)
  results
}
