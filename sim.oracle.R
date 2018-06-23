sim.oracle <- function(method, data) {
  # INPUT
  # method......lasso / alasso / tlasso / slasso - requires lars
  #             enet / aenet - requires elastic net
  #             rlasso - requires relaxo
  #             scad / mcp - requires ncvreg
  #             ols
  # data........all simulation data needed:
  #             beta0, beta, corr, snr, cond, X.sigma, sigma, X.train, y.train
  #
  # OUTPUT:
  # The true prediction error and variable selection performance measures
  #

  source("estimation.R")
  source("resample.R")
    
  # DATA
  beta0 <- data$beta0
  int <- beta0!=0
  beta <- data$beta
  p <- length(beta)
  D <- which(beta != 0)
  d <- length(D)
  
  # OUTPUT
  iter <- length(data$y.train) # number iterations
  n <- length(data$y.train[[1]])

  # model complexity index
  if (method %in% c("lasso", "alasso","slasso",
                    "tlasso","enet","aenet")) {
    comp <- seq(0, 1, by=0.01)
    mode <- "fraction"
  }
  
  if (method %in% c("scad","mcp")) {
    comp <- exp(seq(1.5,-6,len=100))
    mode <- "lambda"
  }
  
  if (method=="rlasso") {
    comp <- 1:p
    mode <- "lambda"
  }
  
  # all fields collected
  beta.names <- paste("b",1:p, sep=".")
  if (int) beta.names <- c("b.0", beta.names)
  beta.prob <- paste0("prob.",beta.names)
  out.names <- c("Method", "IC.Method", "corr", "n" ,"snr", "cond",
                 "compat", "resteig", "irr",
                 "comp", "tune", beta.names, beta.prob,
                 "Nonzero", "Corr.Nonzero", "Inc.Nonzero","Corr.Subset", 
                 "Incl.Corr.Subset", "In.Path", "MSE", "iter")
  n.out <- length(out.names)

  # best models
  tmp <- as.data.frame(matrix(0, iter, n.out))
  colnames(tmp) <- out.names
  best <- list()
  
  # cv and kappa
  if (method == "ols")
    best$Oracle <- best$Full <- tmp else 
  if (method %in% c("tlasso","aenet"))
    best$cv <- tmp else
      best$kappa <- best$cv <- tmp
  n.ic <- length(best)
  for (i in 1:n.ic) {
    best[[i]][,"Method"] <- toupper(method)
    best[[i]][,"IC.Method"] <- names(best)[i]
    best[[i]][,"corr"] <- data$corr
    best[[i]][,"n"] <- n
    best[[i]][,"compat"] <- data$compat
    best[[i]][,"resteig"] <- data$resteig
  }
  var.fixed <- 6 
  
  for (i in 1:iter) {
    
    # TRAINING DATA
    X <- data$X.train[[i]]
    y <- data$y.train[[i]]
    mu <- beta0 + X %*% beta
    
    
    # ESTIMATION
    if (method == "ols") {
      ols <- estimate(X, y, int=int, method="ols")
      oracle <- estimate(X, y, int=int, method="oracle", true=D)
      est <- list(coef=rbind(ols$coef, oracle$coef), 
                  comp.grid=matrix(c(p,d,NA,NA),2))
      ic.select <- 1:2
    } else {
      est <- estimate(X, y, int=int, method=method, comp=comp, mode=mode)
      if (method %in% c("tlasso","aenet"))
        ic.select <- est$cv$index.min else {
          cv <- cv.kfold(X, y, int, method, comp=est$comp, 
                         tune=est$tune, mode=mode, k=10)
          kappa <- res.kappa(X, y, int, method, comp=est$comp, 
                             tune=est$tune, mode=mode, reps=20)
          ic.select <- c(cv$index.min, kappa$index)
        }
    }

    for (j in 1:n.ic) {
      best[[j]][i, "snr"] <- data$snr[i]
      best[[j]][i, "cond"] <- data$cond[i]
      best[[j]][i, "irr"] <- data$irr[i]
      best[[j]][i, "comp"] <- est$comp.grid[ic.select[j],1]
      best[[j]][i, "tune"] <- est$comp.grid[ic.select[j],2]
      betahat <- est$coef[ic.select[j],]
      best[[j]][i, beta.names] <- betahat
      best[[j]][i, beta.prob] <- betahat!=0
      best[[j]][i, "Nonzero"] <- sum(betahat!=0)
      best[[j]][i, "Corr.Nonzero"] <- sum(betahat[D]!=0)
      best[[j]][i, "Inc.Nonzero"] <- sum(betahat[-D]!=0)
      best[[j]][i, "Corr.Subset"] <- best[[j]][i, "Corr.Nonzero"]==d &
                                    best[[j]][i, "Nonzero"]==d
      best[[j]][i, "Incl.Corr.Subset"] <- best[[j]][i, "Corr.Nonzero"] == d
      best[[j]][i, "In.Path"] <-  
        any(rowSums(est$coef!=0)==d &  rowSums(est$coef[,D]!=0)==d)
      if (int) betahat <- betahat[-1]
      diff <- matrix(betahat - beta, ncol=1)
      best[[j]][i, "MSE"] <- t(diff) %*% data$X.sigma %*% diff
      best[[j]][i,"iter"] <- i
    }    
  }

  # average performance measures
  col.no.ave <- 1:var.fixed
  col.ave <- (var.fixed+1):n.out
  bestave <- as.data.frame(matrix(0, n.ic, n.out+5))
  names(bestave) <- c(out.names, "Median.MSE","SE.Median.MSE",
                     "Sq.Bias.Betas", "Var.Betas", "MSE.Betas")
  for (i in 1:n.ic) {
    bestave[i,col.no.ave] <- best[[i]][1,col.no.ave]
    bestave[i,col.ave] <- colMeans(best[[i]][,col.ave])
    bestave[i,"Median.MSE"] <- median(best[[i]]$MSE)
    bestave[i,"SE.Median.MSE"] <- 
      boot.se(best[[i]]$MSE, B=200, median)
    bestave[i,"Sq.Bias.Betas"] <- 
      sum((bestave[i,beta.names]-beta)^2)
    bestave[i,"Var.Betas"] <- 
      sum(apply(best[[i]][,beta.names], 2, var))
    bestave[i,"MSE.Betas"] <- 
      bestave[i,"Sq.Bias.Betas"] + bestave[i,"Var.Betas"]
  }
  bestave <- bestave[,-n.out]
  list(bestave=bestave, best=best)
}
