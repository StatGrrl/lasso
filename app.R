# Code for Application in Chapter 7
source("estimation.R")
source("resample.R")
source("plots.R")
library(tables)
booktabs()
library(Cairo)

# DATA

# data set from lars package
# standardized data: has variable names
library(lars)
data(diabetes)
xnames <- colnames(diabetes$x)

# data from website
# original scale
url <- "http://web.stanford.edu/~hastie/Papers/LARS/diabetes.data"
dat <- read.table(url, head=T)
names(dat) <- c(xnames,"y")

y <- dat[,11]
x <- dat[,-11]
p <- ncol(x)
n <- nrow(x)

# condition number of correlation matrix
eigval <- eigen(cor(x), only.values=T)$values
cond <- sqrt(max(eigval)/min(eigval))

# Plot data
pdf(file = "app-diabetes-data.pdf", width=12, height=12)
plot.corr(x,y,col="steelblue")
dev.off()

# standardized data
ymean <- mean(y)
v <- y-mean(y)
xmean <- colMeans(x)
xnorm <- sqrt(n-1)*apply(x,2,sd)
z <- scale(x, center=xmean, scale=xnorm)
# check
sum(z-diabetes$x)

# training and test split
set.seed(131)
test <- sample(n, floor(n/4))

# ESTIMATION

# training data
z.train <- z[-test,]
v.train <- v[-test]

# methods and tuning parameters
method = c("forward", "ridge", "lasso", 
           "alasso", "slasso", "rlasso", 
           "enet", "mcp", "scad", "tlasso", "aenet", "ols")
tune <- list(forward=NULL, ridge=NULL, lasso=NULL,
             alasso=c(0.5,1,2), slasso=c(0.5,1,2), rlasso=seq(0,1,by=0.1),
             enet=c(0.1,0.5,1,5,10), mcp=seq(2.5,3.5,by=0.1), 
             scad=seq(3,4,by=0.1), tlasso=NULL, aenet=NULL, ols=NULL)

# estimate parameters, perform CV and kappa
est <- list()
cv <- list()
kappa <- list()
for (i in 1:9) {
  if (method[i] %in% c("alasso","slasso","enet"))
    comp <- seq(0,1,by=0.01) else
      comp <- NULL
  est[[i]] <- estimate(z.train, v.train, int=F, method[i], tune=tune[[i]], comp=comp) 
  cv[[i]] <- cv.kfold(z.train, v.train, int=F, method=method[i], comp=est[[i]]$comp, 
                      mode=est[[i]]$mode, tune=est[[i]]$tune, k=10)
 if (is.null(tune[[i]])) t <- NULL else t <- cv[[i]]$tune.min
  kappa[[i]] <- res.kappa(z.train, v.train, int=F, method=method[i], comp=est[[i]]$comp, 
                      mode=est[[i]]$mode, tune=t, reps=20)
}
names(est) <- names(cv) <- names(kappa) <- method[1:9]

tune$aenet <- cv$enet$tune.min
for (i in 10:12) {
  est[[i]] <- estimate(z.train, v.train, int=F, method[i], tune=tune[[i]],
                       comp=seq(0,1,by=0.01), mode="fraction") 
  if (i!=12) cv[[i]] <- est[[i]]$cv
}
names(est)[10:12] <- method[10:12]
names(cv)[10:11] <- method[10:11]

# MODEL SELECTION
norm.ols <- sum(abs(est[[12]]$coef))
tmp <-matrix(0,11,2, dimnames=list(toupper(method[-12]),c("tune","comp")))
comp.best <- list(cv=tmp, kappa=tmp[-(10:11),])
tmp <-matrix(0,11,p, dimnames=list(toupper(method[-12]),xnames))
coef.best <- list(cv=tmp, kappa=tmp[-(10:11),])
rm(tmp)
for (i in 1:11) {
  cv.index <- cv[[i]]$index.min
  index <- cv.index
  if (i < 10) {
    if (is.null(tune[[i]])) 
      kappa.index <- kappa[[i]]$index else 
        kappa.index <- which(est[[i]]$comp.grid[,1]==kappa[[i]]$comp & 
                               est[[i]]$comp.grid[,2]==kappa[[i]]$tune)
    index <- c(index, kappa.index)
  }
  for (j in seq_along(index)) {
    # Best models
    comp.best[[j]][i,] <- rev(t(est[[i]]$comp.grid[index[j],]))
    coef.best[[j]][i,] <- t(est[[i]]$coef[index[j],])
  }


  # Plot CV error and coefficient path
  title <- toupper(method[i])
  
  if (method[i]=="forward") xlab="p"
  if (method[i]=="rlasso") xlab="step"
  if (method[i] %in% c("ridge","lasso","mcp","scad")) xlab=expression(lambda)
  if (method[i] %in% c("alasso","slasso","tlasso","enet","aenet")) xlab="s"

  # cv plot
  filename <- paste0("app-cv-", method[i], ".pdf")
  pdf(file = filename, width=9, height=7)
  cv.plot(cv[[i]], xlab=xlab, title=title)
  dev.off()
  
  # coefficient path
  plot.dat <- cbind(index=1:nrow(est[[i]]$comp.grid),
                    est[[i]]$comp.grid,
                    est[[i]]$coef,
                    norm.frac= rowSums(abs(est[[i]]$coef))/norm.ols)
    
  if (is.null(tune[[i]])) 
    sub <- plot.dat else
      sub <- subset(plot.dat, tune==comp.best[[1]][i,1])

  if (method[i] %in% c("forward","rlasso")) 
    xname <- "comp" else {
      xname <- "norm.frac"
      xlab <- "\u2113   Fraction"
    }
  
  cv.val <- sub[sub$index==cv.index, xname]
  vline <- cv.val
  vlty <- 3
  vlab <- "CV"
  if (i < 10) {
    kappa.val <- sub[sub$index==kappa.index, xname]
    vline <- c(vline, kappa.val)
    vlty <- c(vlty, 2)
    #vlab <- c(vlab, "kappa")
    vlab <- c(vlab, expression(kappa))
  } 
    
  ylab <- "Standardized Coefficients"
  filename <- paste0("app-coef-", method[i], ".pdf")
  CairoPDF(file = filename, width=9, height=7)
  plot.lines(sub, xname, ynames=xnames, 
             show.min=rep(F,p), type="l", 
             lty=c(1,1,1,2,2,1,1,1,2,2),
             col=c("black","limegreen","tomato","lightskyblue",
                   "firebrick","magenta","steelblue","forestgreen",
                   "darkviolet","darkorange"),
             lab=c(xlab,ylab), title=title, vline=vline, vlty=vlty, vlab=vlab)
  if (!method[i] %in% c("forward","rlasso")) 
    mtext("1", side=1, cex=0.6, line=1.6, adj=0.46)
  dev.off()
  
}

# Tuning parameters

norm.best <- lapply(coef.best, abs)
norm.best <- lapply(norm.best, rowSums)
norm.frac <- lapply(norm.best, function(x) x/norm.ols)

tune.best <- as.data.frame(comp.best$cv)
names(tune.best)[2] <- "comp.cv"
tune.best$comp.init <- NA
tune.best$comp.init[10] <- est[[10]]$comp.init
tune.best$comp.init[11] <- est[[11]]$comp.init
tune.best$comp.kappa <- NA
tune.best$comp.kappa[1:9] <- comp.best$kappa[,2]
tune.best$frac.cv <- norm.frac$cv
tune.best$frac.kappa <- NA
tune.best$frac.kappa[1:9] <- norm.frac$kappa
tab <- as.tabular(round(tune.best, 2))
latex(tab)

# BEST MODELS

# coefficients
coef.std <- lapply(coef.best, round)
# table cv
tab <- as.tabular(coef.std$cv)
latex(tab)
# table kappa
tab <- as.tabular(coef.std$kappa)
latex(tab)

# Standard errors
ols <- est[[12]]$model
n.train <- nrow(z.train)
sigma <- sqrt(sum(ols$residuals^2)/(n.train-p-1))
lasso1 <- coef.best$cv[3,]
lambda1 <- comp.best$cv[3,1]
lasso2 <- coef.best$kappa[3,]
lambda2 <- comp.best$kappa[3,1]
coef.best$se <- matrix(0,9,p,
                       dimnames=list(c("OLS", "OLS.SE", "OLS.P",
                                       "LASSO1", "LASSO1.TIBS", "LASSO1.OSB", 
                                       "LASSO2", "LASSO2.TIBS", "LASSO2.OSB"),xnames))
coef.best$se[1,] <- est[[12]]$coef
coef.best$se[2,] <- summary(ols)$coef[,2]
coef.best$se[3,] <- round(summary(ols)$coef[,4],4)
coef.best$se[4,] <- lasso1
coef.best$se[5,] <- lasso.se(z, v, lasso1, sigma, lambda1, type="tibs")
coef.best$se[6,] <- lasso.se(z, v, lasso1, sigma, lambda1, type="osb")
coef.best$se[7,] <- lasso2
coef.best$se[8,] <- lasso.se(z, v, lasso2, sigma, lambda2, type="tibs")
coef.best$se[9,] <- lasso.se(z, v, lasso2, sigma, lambda2, type="osb")

# table se
coef.std <- lapply(coef.best, round, digits=3)
tab <- as.tabular(coef.std$se)
latex(tab)

# PREDICTION

# original scale
coef.orig <- coef.best
coef.orig[[3]] <- NULL
coef.orig <- lapply(coef.orig, scale, center=F, scale=xnorm)
for (i in 1:2) {
  int <- apply(coef.orig[[i]], 1, function(a) ymean - xmean%*%a)
  coef.orig[[i]] <- cbind(int=int,coef.orig[[i]])
}

# test data
x.test <- as.matrix(cbind(1, x[test,]))
y.test <- y[test]
n.test <- length(y.test)

# test error
pred.cv <- apply(coef.orig$cv, 1, function(x) sum((y.test-x.test%*%x)^2)/n.test)
pred.kappa <- apply(coef.orig$kappa, 1, function(x) sum((y.test-x.test%*%x)^2)/n.test)
dat.train <- data.frame(y=y[-test],x=x[-test,] )
dat.test <- data.frame(y=y[test],x=x[test,] )
pred.ols <- sum((y.test - 
                   predict(lm(y~., data=dat.train), dat.test, type="response"))^2)/n.test
pred <- as.data.frame(cbind(CV=pred.cv, kappa=pred.kappa))
pred[10:11,2] <- NA

# table test error
pred.tab <- round(pred)
tab <- as.tabular(pred.tab)
latex(tab)

# plot test error
pred.fig=cbind(Method=toupper(method[-12]),pred)
pred.fig$Method <- factor(pred.fig$Method, 
                          levels=c("FORWARD","RIDGE","LASSO","RLASSO","ALASSO",
                                   "SLASSO","TLASSO","ENET","AENET","MCP","SCAD"))
lab <- levels(pred.fig$Method)
lab[1] <- "FWD"
pdf(file = "app-testerror.pdf", width=11, height=7)
plot.bar(pred.fig, c("CV","kappa"), ylab="Test Error", ycol=c("forestgreen","firebrick"), 
         xname="Method", xlab="", border=NA, ylim=c(0,3500), hline=pred.ols, hlty=2,
         labels=lab)
text(x=par("usr")[2]*1.025, y=pred.ols, labels="OLS", cex=0.75, xpd=T)
dev.off()
