estimate <- function(X, y, method, comp=NULL, mode="fraction", 
                     int=T, true=NULL, tune=NULL, k=10, std=TRUE) {
  # INPUT:
  # X...........matrix or data frame of predictors
  # y...........response vector
  # method......subset selection: 
  #             ols / oracle
  #             forward / backward / exhaustive - requires leaps
  #             ridge - requires MASS
  #             lasso / alasso / tlasso / slasso - requires lars
  #             enet / aenet - requires elastic net
  #             rlasso - requires relaxo
  #             scad / mcp - requires ncvreg
  # comp........complexity parameter, scalar or vector
  #             if comp=NULL, defaults are:
  #             subset selection - 1:p
  #             lars, elasticnet & relaxo - lambda at each step
  # mode........mode of complexity parameter
  #             ridge - lambda
  #             lars - step / fraction / norm / lambda
  #             elasticnet - step / fraction / norm / penalty
  #             relaxo - step / lambda (can only specify step)
  #             mcp / scad - lambda
  # int.........logical, include intercept?
  # true........index of true nonzero parameters, only needed for oracle
  # tune........second tuning parameter, may be vector
  #             alasso, slasso - weight calc, default=1
  #             enet, aenet - ridge penalty parameter, default=0.5
  #             tune must be fixed for aenet (scalar)
  #             relaxo - relaxation parameter, default=0.3
  #             scad, mcp - size est control, scad default=3.7, mcp default=3
  # k...........for k-fold cv, only tlasso and aenet which do sequential CV
  # std.........standardize?
  #
  # OUTPUT:
  # Parameter estimates and df, rss, fitted model
  
  
  p <- ncol(X)
  n <- nrow(X)
  X <- as.matrix(X)
  xnames <- colnames(X)
  if (is.null(xnames)) xnames <- paste0("X",1:p)
  comp.call <- comp
  if (is.null(tune)) 
    tune <- switch(method, "alasso" = 1,
                   "slasso" = 1,
                   "enet" = 0.5,
                   "aenet" = 0.5,
                   "rlasso" = 0.3,
                   "scad" = 3.7,
                   "mcp" = 3,
                   NA)
  
  # standardize
    ymean <- mean(y)
    v <- y-mean(y)  
    xmean <- colMeans(X)
  if (std) {
    xnorm <- sqrt(n-1)*apply(X,2,sd)
    Z <- scale(X, center=xmean, scale=xnorm)
  } else Z <- scale(X, center=xmean, scale=F)
  
  # OLS
  if (method %in% c("ols","alasso","slasso")) {
    model <- lm(v~Z-1)
    coef <- coef(model)
    ols <- coef
    se <- summary(model)$coef[,2]
    df <- p+int
  }
  
  # OLS ORACLE
  if (method == "oracle") {
    d <- length(true)
    model <- lsfit(Z[,true], v, intercept=F)
    coef <- rep(0,p)
    coef[true] <- coef(model)
    comp <- d
    mode <- "p"
    df <- d+int
  }

  # SUBSET SELECTION
  if (method %in% c("forward","backward","exhaustive")) {
    require(leaps)
    model <- regsubsets(v~Z, data.frame(Z,y=y), nvmax=p, 
                        method=method, intercept=F)
    if (is.null(comp)) {
      comp <- 1:p
      mode <- "p"
    }
    b <- coef(model,comp)
    which <- summary(model)$which
    len <- length(comp)
    coef <- matrix(0, nrow=len, ncol=p)
    for (j in 1:len) coef[j, which[j,]] <- b[[j]]
    df <- (1:len)+int
  }

  # RIDGE REGRESSION
  if (method == "ridge") {
    require(MASS)
    mode="lambda"
    if (is.null(comp)) comp <- seq(0, 50, by=0.2)
    model <- lm.ridge(v~Z-1, lambda=comp)
    lambda <- model$lambda
    coef <- coef(model)
    len <- length(comp)
    df <- rep(0,len)
    e <- eigen(t(Z)%*%Z, only.values=T)$values
    for (i in 1:len) for (j in 1:p) df[i] <- df[i]+e[j]/(e[j]+comp[i])
  }
  
  if (!method %in% c("ols","oracle","forward","backward","exhaustive","ridge"))
    coef <- NULL
  for (i in 1:length(tune)) {
    
    # WEIGHTS FOR ADAPTIVE LASSO AND SEA-LASSO
    if (method %in% c("alasso","slasso")) {
      w <- 1/abs(ols)^tune[i]
      if (method == "slasso") w <- w*se
      Z <- scale(Z, center=FALSE, scale=w)
    }
    
    # LASSO, ADAPTIVE AND SEA-LASSO
    if (method %in% c("lasso","alasso","slasso","tlasso")) {
      require(lars)
      model <- lars(Z, v, type="lasso", intercept=F, normalize=F)
      lambda <- model$lambda
      if (is.null(comp)) {
        comp <- model$lambda
        mode <- "lambda"
      }
      init <- predict(model, Z, s=comp, mode=mode, type="coefficients")$coef
      if (length(comp)==1) init <- matrix(init,nrow=1)
      if (method %in% c("alasso","slasso"))
        init <- scale(init, center=FALSE, scale=w)
      coef <- rbind(coef, init)
    }
    
    # ELASTIC NET
    if (method %in% c("enet","aenet")) {
      require(elasticnet)
      model <- enet(Z, v, lambda=tune[i], intercept=F, normalize=F)
      lambda <- model$penalty
      full <- length(lambda)
      if (is.null(comp)) {
        comp <- model$penalty[-c(1,2,full)]
        lambda <- comp
        mode <- "penalty"
      }
      init <- predict(model, Z, s=comp, mode=mode, type="coefficients")$coef
      if (is.null(dim(init))) init <- matrix(init,nrow=1)
      coef <- rbind(coef, init)
    }
    
    # MCP AND SCAD
    if (method %in% c("scad","mcp")) {
      require(ncvreg)
      mode <- "lambda"
      if (is.null(comp)) {
        model <- ncvreg(Z, v, family="gaussian", penalty=toupper(method), 
                        gamma=tune[i])
        comp <- model$lambda
      } else
        model <- ncvreg(Z, v, family="gaussian", penalty=toupper(method), 
                        gamma=tune[i], lambda=comp)
      lambda <- model$lambda
      coef <- rbind(coef, t(model$beta)[,-1])
    }
  }

  # TWO-STAGE LASSO AND ADAPTIVE ELASTIC NET
  if (method %in% c("tlasso","aenet")) {
    source("resample.R")
    parent <- substring(method,2)
    if (is.na(tune)) tune <- NULL
    cv <- cv.kfold(Z, v, int, method=parent, comp=comp, k=k, mode=mode, tune=tune)
    comp.init <- cv$comp.min
    w <- 1/abs(init[cv$index.min,])
    sub <- which(w!=Inf)
    if (length(sub)==1) {
      coef <- init
    } else {
      Z <- scale(Z[,sub], center=FALSE, scale=w[sub])
      if (method == "tlasso") {
        model <- lars(Z, v, type="lasso", intercept=F, normalize=F)
        lambda <- model$lambda
      }
      if (method == "aenet") {
        model <- enet(Z, v, lambda=tune, intercept=F, normalize=F)
        lambda <- model$penalty
      }
      coef <- matrix(0, nrow=length(comp), ncol=p)
      coef[,sub] <- predict(model, Z, s=comp, mode=mode, type="coefficients")$coef
      coef <- scale(coef, center=FALSE, scale=w)
      cv <- cv.kfold(Z, v, int, method=parent, comp=comp, k=k, mode=mode, tune=tune, std=F)
    }
    if (is.null(tune)) tune <- NA
  }
  
  # RELAXED LASSO
  if (method == "rlasso") {
    require(relaxo)
    mode <- "lambda"
    if (!any(tune==1)) {
      tune <- c(tune,1)
      rm1 <- T
    } else rm1 <- F
    model <- relaxo(Z, v, phi=tune)
    coef <- model$beta[order(model$phi),]
    lambda <- unique(model$lambda)
    if (is.null(comp)) comp <- 1:p
    if (length(setdiff(comp,lambda))>0) {
      which <- which(model$lambda[order(model$phi)] %in% lambda[comp])
      coef <- coef[which,]
    }
    if (rm1) {
      comp.grid <- expand.grid(comp=comp, tune=tune)
      which <- comp.grid$tune != 1
      coef <- coef[which,]
      tune <- setdiff(tune,1)
    }
  }
  
  # ORIGINAL SCALE
  if (is.null(dim(coef))) coef <- matrix(coef, nrow=1)  
  if (!method %in% c("ols","oracle","forward","backward","exhaustive","ridge"))
    df <- rowSums(coef!=0)
  if (std) 
    coef <- scale(coef, center=FALSE, scale=xnorm)
  if (int) {
    b0 <- apply(coef, 1, function(a) ymean - xmean%*%a)
    coef <- cbind(b0, coef)  
    X <- cbind(1,X)
    xnames <- c("(Int)", xnames)
  }
  colnames(coef) <- xnames

  if (method=="ols") {
    comp <- p
    mode <- "p"
  }
  comp.grid <- expand.grid(comp=comp, tune=tune)
  fit <- list()
  rss <- list()
  len <- length(comp)
  for (i in 1:length(tune)) {
    fit[[i]] <- apply(coef[1:len+(i-1)*len,,drop=F], 1, function(a) X %*% a)
    rss[[i]] <- apply(fit[[i]], 2, function(a) sum((y - a)^2))
  }
  
  if (all(is.na(tune))) tune <- NULL
  out <- list(model=model, coef=coef, df=df, fit=fit, 
       rss=rss, comp=comp, mode=mode, tune=tune, comp.grid=comp.grid)
  if (!method %in% c("ols","oracle","forward","backward","exhaustive"))
    out$lambda <- lambda
  if (method %in% c("aenet","alasso","slasso","tlasso")) out$w <- w
  if (method %in% c("tlasso","aenet")) {
    out$cv <- cv
    out$comp.init <- comp.init
  }
  out
}

lasso.se <- function(Z, v, alpha, sigma, lambda, type="tibs") {
  
  if (type=="tibs") {
    p <- length(alpha)
    d <- which(alpha!=0)
    w <- rep(0,p)
    w[d] <- 1/abs(alpha[d])
    W <- diag(w)*lambda
  }
  
  if (type=="osb") {
    r <- v - Z %*% alpha
    t <- sum(abs(alpha))
    c <- t(Z) %*% r
    lam <- max(abs(c))
    W <- c %*% t(c) /(t * lam)
  }
  
  R <- t(Z)%*%Z
  inv <- solve(R + W)
  covmat <- sigma^2 * inv %*% R %*% inv
  sqrt(diag(covmat))
}