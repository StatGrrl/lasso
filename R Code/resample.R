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

# CV SELECTION AND 1SE RULE
cv.select <- function(cv, cv.se, comp) {
  # INPUT:
  # cv......vector of cv error estimates
  # cv.se...standard errors of cv estimates
  # comp....model complexity parameter
  #         vector or matrix [comp, tune] for validating 2nd parameter
  #
  # OUTPUT:
  # returns position and value at minimum 
  # and within 1 standard error (1 SE rule)
  # if comp is supplied, also returns its value at min and 1SE

  len <- length(cv)
  index.min <- which.min(cv)
  cv.min <- cv[index.min]
  if (!is.null(dim(comp))) {
    comp.min <- comp[index.min,1]
    tune.min <- comp[index.min,2]
    tlen <- length(unique(comp[,2]))
    if (tlen > 1) {
      cv <- subset(cv, comp[,2]==tune.min)
      cv.se <- subset(cv.se, comp[,2]==tune.min)
      comp <- subset(comp, tune==tune.min)
      len <- length(cv)
    }
  } else comp.min <- comp[index.min]
  cv.lo <- cv.min - cv.se[index.min]
  cv.up <- cv.min + cv.se[index.min]
  index.1SE <- min((1:len)[cv >= cv.lo & cv <= cv.up])
  cv.1SE <- cv[index.1SE]
  if (!is.null(dim(comp))) {
    comp.1SE <- comp[index.1SE,1]
    tune.1SE <- comp[index.1SE,2]
  } else comp.1SE <- comp[index.1SE]
  out <- list(index.min=index.min, index.1SE=index.1SE, 
       cv.min=cv.min, cv.1SE=cv.1SE,
       comp.min=comp.min, comp.1SE=comp.1SE)
  if (!is.null(dim(comp))) {
    out$tune.min <-tune.min
    out$tune.1SE <- tune.1SE
    out$sub.cv <- cv
    out$sub.cv.se <- cv.se
    out$sub.comp <- comp[,1]
  }
  out
}

# K-FOLD CV & LOOCV
cv.kfold <- function(X, y, int=TRUE, method, comp, k=10, 
                     mode="fraction", tune=NULL, std=TRUE) {
  # INPUT:
  # X.........training data X
  # y.........training data y
  # int.......intercept term, logical
  # method....depends on function estimate
  # comp......sequence of points indicating the complexity of the model
  # k.........number of folds
  #           5 / 10 (default) / n (LOOCV)
  # mode......mode of complexity parameter
  # tune......tune can be fixed or a vector
  #           for validating over second parameter
  # std.......standardize?
  #
  # OUTPUT:
  # List containing CV error and its standard error,
  # position, value and complexity at the minimum & using 1SE rule
  
  source("estimation.R")
  
  # split data into k folds
  n <- nrow(X)
  if (k==n) folds <- 1:n  else repeat {
    folds <- sample(k, n, replace=TRUE, prob=rep(1/k, k))
    if (all(1:k %in% folds)) break
  }

  # cross-validate each model
  rep <- length(comp)
  if (!is.null(tune)) rep <- rep*length(tune)
  cv <- matrix(0, rep, k)
  for (i in 1:k) {
    
    # training
    y.train <- y[folds!=i]
    X.train <- X[folds!=i,, drop=FALSE]
    est <- estimate(X.train, y.train, int=int, method=method, 
                    comp=comp, mode=mode, tune=tune, std=std)
    
    # testing
    y.test <- y[folds==i]
    X.test <- X[folds==i,, drop=FALSE]
    n.test <- length(y.test)
    if (int) X.pred <- cbind(1,X.test) else X.pred <- X.test
    pred <- apply(est$coef, 1, function(a) X.pred %*% a)
    if (n.test == 1) pred <- matrix(pred, nrow=1)
    cv[,i] <- apply(pred, 2, function(a) sum((y.test - a)^2)/n.test)
  }
  cv.mean <- apply(cv, 1, mean)
  cv.se <- apply(cv, 1, sd)
  comp.grid <- if (is.null(tune)) comp else
    expand.grid(comp=comp, tune=tune)
  select <- cv.select(cv.mean, cv.se, comp.grid)
  out <- list(cv=cv.mean, cv.se=cv.se, comp=comp, mode=mode, 
              index.min=select$index.min, index.1SE=select$index.1SE, 
              cv.min=select$cv.min, cv.1SE=select$cv.1SE, 
              comp.min=select$comp.min, comp.1SE=select$comp.1SE,
              tune=tune)
  if (length(tune) > 1) {
    out$tune.min <- select$tune.min
    out$tune.1SE <- select$tune.1SE
    out$comp.grid <- comp.grid
    out$sub.cv <- select$sub.cv
    out$sub.cv.se <- select$sub.cv.se
    out$sub.comp <- select$sub.comp
  }
  out
}

cv.plot <- function(cv.obj, xlab="", logx=F, topx=NULL, toplab=NULL,
                    title=NULL, alt.x=NULL, alt.min=NULL) {
  # INPUT:
  # cv.obj....object from cv.kfold
  # sub.......subset, quote(logical) 
  # xlab......x axis label, x is complexity index used for CV
  # logx......-log(x), logical
  # topx......secondary x axis
  # toplab....label for secondary x axis
  # title.....plot title
  # alt.x.....alternative x axis same length a comp
  # alt.min...value of alt.x at minimum CV
  #
  # OUTPUT:
  # Plot of CV error, standard error band and position of min & 1SE
  
  # compute points
  if (length(cv.obj$tune) > 1) {
    cv.obj$cv <- cv.obj$sub.cv
    cv.obj$cv.se <- cv.obj$sub.cv.se
    cv.obj$comp <- cv.obj$sub.comp
  }
  cv.upper <- cv.obj$cv+cv.obj$cv.se
  cv.lower <- cv.obj$cv-cv.obj$cv.se
  c <- cv.obj$comp
  c.min <- cv.obj$comp.min
  c.1SE <- cv.obj$comp.1SE
  if (logx) {
    c <- -log(c)
    c.min <- -log(c.min)
    c.1SE <- -log(c.1SE)
  }
  if (!is.null(alt.x)) {
    c <- alt.x
    c.min <- alt.min
    c.1SE <- NULL
  }
  
  # graphics parameters
  par(mar=c(2,3,1,2)+0.1, mgp=c(1.5,0.25,0), cex.axis=0.75, tcl=-0.25)
  if (!is.null(topx)) par(mar=par("mar")+c(0,0,2,0))
  if (!is.null(title)) par(mar=par("mar")+c(0,0,1,0))
  ylim <- c(min(cv.lower),max(cv.upper))
  
  # plot CV error
  plot(c, cv.obj$cv, type="l", ylab="CV Error", xlab="", ylim=ylim)
  mtext(xlab, side=1, line=0.5, at=par("usr")[2], las=1)
  
  # standard error bands
  points(c, cv.upper, pch=95)
  points(c, cv.lower, pch=95)
  for (i in seq_along(c)) segments(x0=c[i], y0=cv.lower[i], y1=cv.upper[i]*0.997)
  
  # position of min & 1SE
  points(c.min, cv.obj$cv.min, col="red", pch=16)
  abline(v=c.min, lty=2)
  text(c.min, ylim[2], labels=" Min", adj=0, cex=0.8)
  if (!is.null(c.1SE)) {
    points(c.1SE, cv.obj$cv.1SE, col="red", pch=16)
    abline(v=c.1SE, lty=2)
    text(c.1SE, ylim[2], labels=" 1SE", adj=0, cex=0.8)
  }
  # secondary axis
  if (!is.null(topx)) {
    axis(side=3, at=c, labels=topx)
    mtext(toplab, side=3, line=0.25, adj=1, at=par("usr")[2] )
  }
  
  # title
  if (!is.null(title)) {
    if (!is.null(topx)) line <- 2 else line <- 0.75
    if (!is.null(title)) title(title, line=line)
  }
}

# PERCENTILE CV
cv.perc <- function(X, y, int=TRUE, method, comp, k=10, 
                    reps=10, perc=0.95, mode="fraction", tune=NULL, std=T) {
  # INPUT:
  # X.........training data X
  # y.........training data y
  # int.......intercept term, logical
  # method....depends on function estimate
  # comp......sequence of points indicating the complexity of the model
  # k.........number of folds
  #           5 / 10 (default) / n (LOOCV)
  # reps......number of cv repititions
  # perc......percentile of complexity parameter
  # mode......mode of complexity parameter, must specify if lambda/penalty
  # tune......secondary tuning parameter, must be fixed
  #           only for estimation, not validated
  # std.......standardize?
  #
  # OUTPUT:
  # List containing index and complexity at the percentile
  
  # repeat cross-validation
  cv.comp <- rep(0,reps)
  for (i in 1:reps) {
    cv <- cv.kfold(X, y, int, method, comp, k, mode, tune, std)
    cv.comp[i] <- cv$comp.min
  }
  # find percentile
  if (mode %in% c("lambda","penalty"))
    best.perc <- quantile(cv.comp, perc) else
      best.perc <- quantile(cv.comp, 1-perc)
  best.index <- findInterval(best.perc, comp)
  list(index=best.index, comp=comp[best.index])
}

# KAPPA SIMILARITY MEASURE
kappa.sim <- function(a,a1,a2) {
  # INPUT:
  # a..........full set
  # a1.........elements in first set
  # a2.........elements in second set
  #
  # OUTPUT:
  # kappa coefficient of agreement between two sets a1 and a2
  
  len.a <- length(a)
  len.a1 <- length(a1)
  len.a2 <- length(a2)
  if ((len.a1==0 & len.a2==0) | (len.a1==len.a & len.a2==len.a))
    kap <- -1 else {
    a1c <- setdiff(a, a1)
    a2c <- setdiff(a, a2)
    obs <- (length(intersect(a1,a2)) + length(intersect(a1c,a2c)))/len.a
    chance <- (length(intersect(a1,a2)) + length(intersect(a1,a2c)))*
      (length(intersect(a1,a2)) + length(intersect(a1c,a2)))/len.a^2 +
      (length(intersect(a1,a2c)) + length(intersect(a1c,a2c)))*
      (length(intersect(a1c,a2)) + length(intersect(a1c,a2c)))/len.a^2
    kap <- (obs-chance)/chance
    }
  kap
}

# RESAMPLING KAPPA
res.kappa <- function(X, y, int=TRUE, method, comp, 
                     reps=10, perc=0.9, mode="fraction", tune=NULL, std=T) {
  # INPUT:
  # X.........training data X
  # y.........training data y
  # int.......intercept term, logical
  # method....depends on function estimate
  # comp......sequence of points indicating the complexity of the model
  # reps......number of repititions
  # perc......percentile of complexity parameter
  # mode......mode of complexity parameter, must specify if lambda/penalty
  # tune......second tuning parameter, may be vector
  # std.......standardize?
  #
  # OUTPUT:
  # List containing index and complexity using kappa coefficient
 
  source("estimation.R")
  
  n <- nrow(X)
  p <- ncol(X)
  len <- length(comp)
  tlen <- length(tune)
  if (tlen > 1) len <- len*tlen
  kap.hat <- matrix(0, len, reps, perc)
  for (i in 1:reps) {
    # split data into 2 folds
    split <- sample(n, floor(n/2))

    # estimate model on each fold
    est1 <- estimate(X[split,], y[split], int=int, method=method, 
                     comp=comp, mode=mode, tune=tune, std=std)
    active1 <- apply(est1$coef, 1, function(x) which(x!=0))
    est2 <- estimate(X[-split,], y[-split], int=int, method=method, 
                     comp=comp, mode=mode, tune=tune, std=std)
    active2 <- apply(est2$coef, 1, function(x) which(x!=0))
    
    # calculate kappa coefficient
    for (j in 1:len)
      kap.hat[j,i] <- kappa.sim(1:p, active1[[j]], active2[[j]])
  }  
  
  # find best kappa
  kap.ave <- apply(kap.hat, 1, mean)
  kap.ratio <- kap.ave/max(kap.ave)
  if (mode %in% c("lambda","penalty"))
    best.index <- min(which(kap.ratio>=perc)) else
      best.index <- max(which(kap.ratio>=perc))
  out <- list(kap=kap.ave, index=best.index)
  if (is.null(tune)) tune <- NA
  if (tlen <= 1) {
    out$comp <- comp[best.index] 
    out$tune <- tune } else {
    out$comp.grid <- expand.grid(comp=comp, tune=tune)
    out$comp <- out$comp.grid[best.index,1]
    out$tune <- out$comp.grid[best.index,2]
  }
  out
}