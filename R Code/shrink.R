func.rss <- function(X, y, alpha) {
  # INPUT:
  # X...........matrix of predictors
  # y...........response vector
  # alpha.......regression parameter matrix with each 
  #             column corresponding to a predictor
  #             vectors are coerced to matrices
  #
  # OUTPUT:
  # values of residual sum of squares
  
  if (is.null(dim(alpha))) alpha <- matrix(alpha, ncol=ncol(x))
  apply(alpha, 1, function(x) sum((y - X %*% x)^2))
}

check.param <- function(method, lambda, tune) {
  # INPUT: 
  # method......name of shrinkage method
  #             "subset", "ridge", "lasso", "bridge", 
  #             "rlasso", "alasso", "flasso",
  #             "nenet","enet", "aenet", "oscar", "pacs", "scad", "mcp", 
  #             "glassp", "cap", "gbridge"
  # lambda......shrinkage parameter, may be vector
  # tune........secondary tuning parameter, may be vector
  #
  # OUTPUT:
  # checks that tuning parameters are within their constraints
  
  if (any(lambda<0)) stop("lambda must be nonnegative")
  if (method %in% c("bridge","alasso","flasso","oscar","pacs") & 
        any(tune<0)) stop("tune must be positive for this method")
  if (method %in% c("rlasso","nenet","enent","aenet","gbridge") & 
        (any(tune<0) | any(tune>1))) 
    stop("tune must be between 0 and 1 for this method")
  if (method=="scad" & any(tune<=2)) 
    stop("tune must be greater than 2 for this method")
  if (method=="mcp" & any(tune<=1)) 
    stop("tune must be greater than 1 for this method")
  if (method=="cap" & any(tune<1)) 
    stop("tune must be greater or equal to 1 for this method")
}

penalty <- function(method, alpha, vars, lambda, tune=NULL, b.init) {
  # INPUT: 
  # method......name of shrinkage method
  #             "subset", "ridge", "lasso", "bridge", "rlasso", "alasso", "flasso",
  #             "nenet","enet", "aenet", "oscar", "pacs", "scad", "mcp",
  #             "glasso", "cap", "gbridge"
  # alpha.......regression parameter matrix with each column 
  #             corresponding to a predictor
  #             vectors are coerced to matrices
  #             for group penalties, alpha is a list containing a
  #             regression parameter matrix for each group, 
  #             vector coerced to matrices
  # vars........number of predictors (to ensure conversion is done correctly)
  #             for group penalties, vars is a vector with group sizes
  # lambda......shrinkage parameter, must be scalar
  # tune........secondary tuning parameter, must be scalar
  #             for cap penalty, tune is a vector 
  #             specifying the Lq norm for each group
  #             if it is given as a scalar, 
  #             the same Lq norm is used for all groups
  # b.init......initial estimate require by some penalties, 
  #             same length as each alpha column
  #
  # OUTPUT:
  # values of the penalty function
  
  check.param(method, lambda, tune)
  if (method %in% c("glasso","cap","gbridge")) {
    len <- length(alpha)
    for (i in 1:len) {
      if (is.null(dim(alpha[[i]]))) 
        alpha[[i]] <- matrix(alpha[[i]],ncol=vars[i])
    }
  } else {
    if (is.null(dim(alpha))) alpha <- matrix(alpha, ncol=vars)
    len <- nrow(alpha)
    alpha.abs <- abs(alpha)
    if (method %in% c("alasso","aenet","pacs")) {
      if (any(b.init==0)) stop("initial estimates cannot be zero")
      b.abs <- abs(b.init)
    }
  }

  if (method == "subset") pen <- lambda^2*rowSums(alpha.abs!= 0)/2
  
  if (method == "ridge") pen <- lambda*rowSums(alpha.abs^2)
  
  if (method == "lasso") pen <- lambda*rowSums(alpha.abs)
  
  if (method == "bridge") {
    if (tune==Inf)
      pen <- lambda*apply(alpha.abs, 1, max) else
        pen <- lambda*rowSums(alpha.abs^tune)
  }
  
  if (method == "rlasso") pen <- lambda*tune*rowSums(alpha.abs)
  
  if (method == "alasso") pen <- lambda*rowSums(alpha.abs/b.abs^tune)
  
  if (method == "flasso") {
    pen <- lambda*rowSums(alpha.abs)
    for (i in 1:len) for (j in 2:vars) 
      pen[i] <- pen[i] + lambda*tune*abs(alpha.abs[i,j]-alpha.abs[i,j-1])
  }
  
  if (method == "nenet") 
    pen <- lambda*(1-tune)*rowSums(alpha.abs) + lambda*tune*rowSums(alpha.abs^2)
  
  if (method == "enet") 
    pen <- (1+tune*lambda)*(lambda*(1-tune)*rowSums(alpha.abs) + 
                              lambda*tune*rowSums(alpha.abs^2))
  
  if (method == "aenet") 
    pen <- (1+tune*lambda)*(lambda*(1-tune)*rowSums(alpha.abs/b.abs^tune) + 
                              lambda*tune*rowSums(alpha.abs^2))

  if (method == "oscar") {
    pen <- lambda*(1-tune)*rowSums(alpha.abs)
    for (i in 1:len) for (j in 1:vars) 
      pen[i] <- pen[i] + lambda*tune*max(alpha.abs[i,j:vars])
  }

  if (method == "pacs") {
    pen <- lambda*rowSums(alpha.abs/b.abs^tune)
    for (i in 1:len) for (j in 1:vars) {
      pen[i] <- pen[i] + 
        lambda*sum(abs(alpha[i,(j+1):vars]-alpha[i,j])/
                     abs(b.init[(j+1):vars]-b.init[j])^tune) +
        lambda*sum(abs(alpha[i,(j+1):vars]+alpha[i,j])/
                     abs(b.init[(j+1):vars]+b.init[j])^tune)
    }
  }

  if (method == "scad") 
    pen <- rowSums(lambda*alpha.abs*(alpha.abs <= lambda) - 
    (alpha.abs^2-2*tune*lambda*alpha.abs+lambda^2)/
      (2*(tune-1))*(alpha.abs>lambda & alpha.abs<=tune*lambda) +
    ((tune+1)*lambda^2/2)*(alpha.abs>tune*lambda))
  
  if (method == "mcp") 
    pen <- rowSums((lambda*alpha.abs-alpha.abs^2/(2*tune))*(alpha.abs<=tune*lambda) +
    (tune*lambda^2/2)*(alpha.abs>tune*lambda))
  
  if (method=="glasso") {
    grp.norms <- matrix(0, nrow(alpha[[1]]), len)
    for (i in 1:len) grp.norms[,i] <- vars[i]*sqrt(rowSums(alpha[[i]]^2))
    pen <- lambda*rowSums(grp.norms)
  }
  
  if (method=="cap") {
    if (length(tune)==1) tune <- rep(tune, len)
    grp.norms <- matrix(0, nrow(alpha[[1]]), len)
    for (i in 1:len) {
      grp.norms[,i] <- if (tune[i]==Inf) 
        apply(alpha[[i]], 1, function(x) max(abs(x))) else 
        rowSums(abs(alpha[[i]])^tune[i])
    }
    pen <- lambda*rowSums(grp.norms)
  }
  
  if (method=="gbridge") {
    grp.norms <- matrix(0, nrow(alpha[[1]]), len)
    for (i in 1:len) 
      grp.norms[,i] <- (vars[i]^(1-tune))*(rowSums(abs(alpha[[i]])))
    pen <- lambda*rowSums(grp.norms^tune)
  }
  pen
}

threshold <- function(method, lambda=2, tune=NULL, ols=seq(-10,10,by=0.1)) {
  # INPUT: 
  # method......name of shrinkage method
  #             "subset", "ridge", "lasso", "bridge", "rlasso", "alasso", 
  #             "nenet", "scad", "mcp"
  # lambda......shrinkage parameter, may be vector
  # tune........secondary tuning parameter, may be vectors
  # ols.........least squares values which will be thresholded
  #
  # OUTPUT:
  # values of the penalty function
  
  check.param(method, lambda, tune)
  ols.abs <- abs(ols)
  ols.sign <- sign(ols)
  ols.len <- length(ols)
  if (is.null(tune)) tune <- 1
  reps <- length(lambda)*length(tune)
  alpha <- matrix(0, ols.len, reps)
  tun.param <- matrix(0, reps, 2)
  names(tun.param) <- c("lambda","tune")
  k <- 1
  for (l in lambda) for (t in tune) {
    if (method=="subset") alpha[,k] <- ols*(ols.abs > l)
    if (method=="ridge") alpha[,k] <- ols/(1 + l)
    if (method=="lasso") alpha[,k] <- ols.sign*(ols.abs - l)*(ols.abs >= l)
    if (method=="bridge") {
      alpha[,k] <- ols + ols.sign*l*t*ols.abs^(t-1)
      alpha[ols==0,] <- 0
    }
    if (method=="rlasso") 
      alpha[,k] <- ols.sign*(ols.abs - t*l)*(ols.abs >= l)
    if (method=="alasso") {
      alpha[,k] <- ols.sign*(ols.abs - l/(ols.abs^t))*(ols.abs > l/(ols.abs^t))
      alpha[ols==0,] <- 0
    }
    if (method=="nenet") 
      alpha[,k] <- ols.sign*(ols.abs - l*(1-t))*(ols.abs >= l*(1-t))/(1+t*l)
    if (method=="scad") 
      alpha[,k] <- ols.sign*(ols.abs - l)*(ols.abs >= l & ols.abs <= 2*l) +
      (((t-1)*ols - ols.sign*t*l)/(t-2))*(ols.abs > 2*l & ols.abs <= t*l) +
      ols*(ols.abs > t*l)
    if (method=="mcp") 
      alpha[,k] <- ols.sign*(ols.abs - l)*(1-1/t)*(ols.abs >= l & ols.abs <= t*l) +
      ols*(ols.abs > t*l)
    tun.param[k,] <- c(l, t)
    k <- k+1
  }
  alpha
}

shrink.labs <- function(method, lambda, tune=NULL, leg=F) {
  # INPUT: 
  # method......name of shrinkage method
  #             "subset", "ridge", "lasso", "bridge", "rlasso", "alasso", "flasso",
  #             "nenet","enet", "aenet", "oscar", "pacs", "scad", "mcp"
  # lambda......shrinkage parameter, may be vector
  # tune........secondary tuning parameter, may be a vector
  # leg.........should legend labels be produced?
  #
  # OUTPUT:
  # labels corresponding to each shrinkage method
  
  method.subs <- switch(method, 
                       "subset"="SS", "ridge"="R",  "lasso"="L",
                       "bridge"="B", "rlasso"="RL", "alasso"="AL",
                       "flasso"="FL", "nenet"="NE", "enet"="E",
                       "aenet"="AE", "oscar"="O", "pacs"="P",
                       "scad"="S", "mcp"="M", "glasso"="GL", 
                       "cap"="C", "gbridge"="GB", "slasso"="SL") 
  alpha.lab <- eval(substitute(paste0("hat(alpha)^", m), list(m=method.subs)))
  alpha.lab <- parse(text=alpha.lab)
  pen.lab <- eval(substitute(paste0("P[", m, "]"), list(m=method.subs)))
  pen.lab <- parse(text=pen.lab)
  
  if (is.null(tune)) tune <- 1
  leg.names <- vector("expression", length(lambda)*length(tune))
  k <- 1
  for (i in lambda) for (j in tune) {
    leg.names[k] <- bquote(paste(lambda == .(i)))[2]
    if (method %in% c("bridge","cap","gbridge")) 
      leg.names[k] <- bquote(paste(lambda == .(i) * "," ~ gamma == .(j)))[2]
    if (method == "rlasso") 
      leg.names[k] <- bquote(paste(lambda == .(i) * "," ~ phi == .(j)))[2]
    if (method %in% c("alasso","aenet","pacs","slasso")) 
      leg.names[k] <- bquote(paste(lambda == .(i) * "," ~ zeta == .(j)))[2]
    if (method %in% c("flasso", "nenet", "enet", "oscar")) 
      leg.names[k] <- bquote(paste(lambda == .(i) * "," ~ psi == .(j)))[2]
    if (method %in% c("scad", "mcp")) 
      leg.names[k] <- bquote(paste(lambda == .(i) * "," ~ xi == .(j)))[2]
    k <- k+1
  }

  method.name <- switch(method, 
                        "subset"="Subset Selection",
                        "ridge"="Ridge Regression",
                        "lasso"="LASSO",
                        "bridge"="Bridge Estimates",
                        "rlasso"="Relaxed LASSO",
                        "alasso"="Adaptive LASSO",
                        "flasso"="Fused LASSO",
                        "nenet"="Naive Elastic Net",
                        "enet"="Elastic Net",
                        "aenet"="Adaptive Elastic Net",
                        "oscar"="OSCAR",
                        "pacs"="PACS",
                        "scad"="SCAD", 
                        "mcp"="MCP", 
                        "glasso"="Group LASSO",
                        "cap"="CAP",
                        "gbridge"="Group Bridge")
  title <- method.name
  if (!leg) 
    title <- parse(text=paste(gsub(" ","~ ", title), 
                              deparse(leg.names[[1]]), sep="~"))
  
  labs <- list(alpha=alpha.lab, penalty=pen.lab, title=title)
  if (leg) labs$legend <- leg.names
  labs
}

plot.shrink <- function(type, method, lambda=2, tune=NULL, 
                        grid=seq(-10,10,by=0.1), b.init=grid, 
                        lab=NULL, title=NULL, lty=1, lwd=1, 
                        col="black", ymax=NULL) {
  # INPUT: 
  # type............type of plot
  #                 "threshold", "penalty"
  # method..........name of shrinkage method
  #                 threshold plot depends on function threshold
  #                 penalty plot depends on function penalty
  # lambda..........shrinkage parameter, overlayed if vector
  # tune............secondary tuning parameter, overlayed if vector
  # grid............vector of values for which to calculate the function
  #                 for theshold this must be ols, for others alpha
  # b.init..........initial estimate require by some penalties
  # lab.............vector c(xlab,ylab). 
  #                 If NULL, will be calculated using shrink.labs
  #                 Use c("","") to omit
  # title...........specified title
  #                 If NULL, will be calculated using shrink.labs
  #                 Use "" to omit
  # lty, lwd, col...graphical parameters
  #                 recycled for length(lambda)*length(tune)
  # ymax............value to limit y axis
  #
  # OUTPUT:
  # plot of either the thresholding function or penalty function
  
  if (is.null(tune)) tune <- 1
  ad <- length(lambda)*length(tune)
  par(mar=c(2,2,1.5,1)+0.1, mgp=c(1.5,0.25,0), 
      cex.main=1.2, cex.axis=0.75, tcl=-0.25)
  if (is.null(lab) | is.null(title) | ad > 1) {
    if (ad > 1) {
      leg <- T
      if (method %in% c("subset", "ridge", "lasso")) {
        par(mar=par("mar")+c(0,0,0,4))
      } else par(mar=par("mar")+c(0,0,0,6))
      pos <- ypos <- 1.04*max(grid) - 0.04*min(grid)
      lty <- rep(lty, length=ad)
      lwd <- rep(lwd, length=ad)
      col <- rep(col, length=ad)
    } else leg <- F
    getlabs <- shrink.labs(method, lambda, tune, leg)
    if (is.null(lab)) {
      ols.lab <- expression(hat(alpha))
      alpha.lab <- getlabs$alpha
      pen.lab <- getlabs$penalty
    } else {
      xlab <- lab[1]
      ylab <- lab[2]
    }
    if (is.null(title)) title <- getlabs$title
  }
  
  if (type=="threshold") {
    alpha <- threshold(method, lambda, tune, ols=grid)
    plot(grid, grid, type="l", xlab="", ylab="", main=title, lty=3)
    if (is.null(lab)) {
      xlab <- ols.lab
      ylab <- alpha.lab
    }
    mtext(xlab, side=1, line=0.5, at=par("usr")[2], las=1)
    mtext(ylab, side=2, line=0.5, at=par("usr")[4], las=1)
    for (i in 1:ad) {
      if (method=="bridge") 
        lines(alpha[,i], grid, lty=lty[i], lwd=lwd[i], col=col[i]) else
        lines(grid, alpha[,i], lty=lty[i], lwd=lwd[i], col=col[i]) 
    }
  }
  
  if (type=="penalty") {
    if (method %in% c("alasso","aenet","pacs") & (any(b.init==0))) {
      inactive <- which(b.init==0)
      grid <- grid[-inactive]
      b.init <- b.init[-inactive]
    }
    pen <- NULL
    for (i in seq_along(lambda)) for (j in seq_along(tune)) {
      pen <- cbind(pen, penalty(method, alpha=grid, vars=1, 
                                lambda[i], tune[j], b.init) )
    }
    if (is.null(ymax))  ymax <- max(pen) 
    plot(grid, rep(0,length(grid)), type="n", xlab="", ylab="", 
         main=title, ylim=c(0,ymax))
    if (is.null(lab)) {
      xlab <- alpha.lab
      ylab <- pen.lab
    }
    mtext(xlab, side=1, line=0.5, at=par("usr")[2], las=1)
    mtext(ylab, side=2, line=0.5, at=1.05*ymax, las=1)      
    for (k in 1:ncol(pen)) lines(grid, pen[,k], lty=lty[k], 
                                 lwd=lwd[k], col=col[k])
    ypos <- 1.04*ymax
  }
  if (leg)
    legend(pos, ypos, legend=getlabs$legend, cex=0.8, 
                  col=col, lwd=lwd, lty=lty, bty="n", xpd=T)  
}

contours.rss <- function(X, y, method, sd.step=10, title=NULL) {
  # INPUT: 
  # X...............matrix of predictors
  # y...............response vector
  # method..........name of shrinkage method
  #                 ridge, lasso
  # sd.step.........used to vary the range of regression parameters
  #                 the range is calculated based on 
  #                 the least squares coefficient
  #                 and sd.steps * its standard deviation
  # title...........specified title
  #
  # OUTPUT:
  # plot of constrained problem, RSS contours touching constraint region
  
  source("estimation.R")
  
  # data
  n <- nrow(X)
  p <- ncol(X)
  
  # standardize
  v <- y - mean(y)
  sxx <- sqrt(n-1) * apply(X, 2, sd)
  Z <- scale(X, center=T, scale=sxx)
  
  # ols
  ols <- lm(v~Z-1)
  ols.coef <- ols$coef
  ols.sd <- coef(summary(ols))[,2]
  ols.rss <- sum(ols$residuals^2)
  ols.l2 <- sum(ols.coef^2)
  ols.l1 <- sum(abs(ols.coef))
  
  s <- seq(0.2, 0.9, by=0.1)
  if (method == "ridge") comp <- seq(0, ols.l2, by=0.2)
  if (method == "lasso") comp <- s
  est <- estimate(Z, v, int=F, method=method, comp=comp, mode="fraction")
  if (method=="ridge") {
    ratio.l2 <- round(rowSums(est$coef^2)/ols.l2,2)
    uniq.index <-  which(!duplicated(ratio.l2))
    uniq.val <-  ratio.l2[uniq.index]
    s.index <- uniq.index[sapply(s, function(x) 
                                 which.min(abs(uniq.val-x)))]
    est$coef <- est$coef[s.index,]
    est$df <- est$df[s.index]
    est$lambda <- est$lambda[s.index]
    pen.levels <- rowSums(est$coef^2)
  }
  if (method=="lasso") pen.levels <- s*ols.l1
  rss.levels <- func.rss(Z, v, est$coef)
  
  alpha.step <- 101
  steps <- seq(-1,1,length=alpha.step)[-c(1,alpha.step)]
  alpha1 <- ols.coef[1] + sd.step*ols.sd[1]*steps
  alpha2 <- ols.coef[2] + sd.step*ols.sd[2]*steps
  minmax <- c(min(alpha1),max(alpha1))
  maxind <- which.max(abs(minmax))
  a1ext <- steps*(minmax[maxind])
  if (maxind==1) alpha1 <- c(alpha1, a1ext[a1ext>minmax[-maxind]])
  if (maxind==2) alpha1 <- c(a1ext[a1ext<minmax[-maxind]], alpha1)
  rss <- outer(alpha1, alpha2, function(a,b) func.rss(Z, v, cbind(a,b)))
  pen <- outer(alpha1, alpha2, function(a,b) 
              penalty(method, cbind(a,b), vars=2, lambda=1))
  
  par(mar=c(3,3,0.5,1)+0.1, mgp=c(1.5,0.25,0), cex.axis=0.75, tcl=-0.25)
  if (!is.null(title)) par(mar=par("mar")+c(0,0,1,0))
  xlab <- expression(alpha[1])
  ylab <- expression(alpha[2])
  contour(x=alpha1, y=alpha2, z=rss, method="simple", drawlabels=F, 
          levels=rss.levels, lwd=2, main=title, xlab=xlab, ylab=ylab)
  contour(x=alpha1, y=alpha2, z=pen, method="simple", drawlabels=F, 
          col="red", add=T, levels=pen.levels, lwd=2)
  for (i in seq_along(pen.levels)) points(x=est$coef[i,1], 
                                          y=est$coef[i,2], pch=17, col="red")
  points(x=ols.coef[1], y=ols.coef[2], pch=16)
  abline(h=0, v=0, lty=3)
  est$coef <- est$coef/sxx
  dimnames(est$coef) <- list(s, names(X))
  est
}

plot.normball <- function(method, dim, tune=NULL) {
  # INPUT: 
  # method......name of shrinkage method
  # dim.........dimension of norm ball, 2 or 3 (3D requires misc3d)
  #             dim must be 3 for group penalties, where
  #             first two dimensions are a group
  # tune........secondary tuning parameter, scalar
  #             may be a vector for cap penalty
  #
  # OUTPUT:
  # plot of norm ball for penalties using Lq norms
  
  if (dim==2) {
    alpha <- seq(-1,1,len=500)
    pen <- outer(alpha, alpha, function(x,y) 
      penalty(method, cbind(x,y), vars=2, lambda=1, tune))
    image(alpha, alpha, pen, breaks=c(0,1), col="red", 
          axes=F, xlab="", ylab="")
  }
  if (dim==3) {
    require(misc3d)
    vars <- 3
    a <- seq(-1,1,len=50)
    alpha <- expand.grid(x = a, y = a, z = a)
    if (method %in% c("glasso", "cap", "gbridge")) {
      vars <- c(2,1)
      alpha <- list(alpha[,1:2], alpha[,3])
    }
    pen <- penalty(method, alpha, vars=vars, lambda=1, tune)
    v <- array(pen, rep(length(a),3))
    con <- computeContour3d(v, max(v), level=1)
    drawScene(makeTriangles(con))
  }
}
