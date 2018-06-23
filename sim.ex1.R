# Code for Simulations in Section 6.2

# DATA
library(tables)
source("sim.data.R")
source("sim.ic.R")
source("plots.R")

beta=c(3, 1.5, 0, 0, 2, 0, 0, 0)
p <- length(beta)
n <- 25
sigma <- c(1,3,6)
rho <- c(0, 0.5, 0.9)
vary <- expand.grid(sigma=sigma, rho=rho)
v <- nrow(vary)

data.ex1 <- list()
for (i in 1:v) {
  data.ex1[[i]] <- sim.data(iter=100, n.train=n, n.test=200, n.valid=n, 
                            sigma=vary$sigma[i], corr="AR", 
                            rho=vary$rho[i], beta, beta0=0, seed=2020)
}
snr <- sapply(data.ex1, "[[", 7)
cond <- sapply(data.ex1, "[[", 8)
vary$snr <- colMeans(snr)
vary$cond <- colMeans(cond)

booktabs()
# table snr and condition number
tab <- tabular( Factor(rho, "$\\rho$")~ 
                  Format(digits=3)*
                  (Heading("SNR")*snr+Heading("$\\kappa_2(\\mathbf{\\Sigma})$")*cond)*
                  Heading()*identity*
                  Factor(sigma, "$\\sigma$"), data=vary)
latex(tab)

# SIMULATIONS AND PLOTS
method = c("lasso", "ridge", "forward", "ols")

for (i in seq_along(method)) {
  for (j in 1:nrow(vary)) {
      
    # simulations
    print(paste(method[i], j, Sys.time()), quote=F)
    sim.tmp <- sim.ic(method=method[i],data=data.ex1[[j]])
    if (i+j == 2) {
      best <- do.call("rbind", sim.tmp$best)
      bestave <- sim.tmp$bestave
    } else {
      best <- rbind(best, do.call("rbind", sim.tmp$best))
      bestave <- rbind(bestave, sim.tmp$bestave)
    }
    
    # do plots on path data now and discard data
    if (method[i] != "ols") {
      pathave <- if (i+j == 2) sim.tmp$pathave else 
        rbind(pathave, sim.tmp$pathave)

      title <- substitute(paste(m, ": ", sigma == n, ", ", rho == o), 
                          list(m=toupper(method[i]), 
                               n=vary$sigma[j], o=vary$rho[j]))
      xlab <- switch(method[i], "forward"="p", "lasso"="s", 
                     "ridge"=expression(-log(lambda)))
      if (method[i] == "ridge") logx <- T else logx <- F

      # PREDICTION ERROR
      
      # test error and training error for each sample
      filename <- paste("sim-ex1-", j, "-", method[i], 
                        "-test-train.pdf", sep = "")
      pdf(file = filename, width=9, height=7)
      plot.path(sim.tmp$path, xname="comp", ynames=c("Train","Test"), 
                col=c("steelblue","firebrick"), 
                title=title, logx=logx, topx=sim.tmp$pathave$df.cov, 
                toplab="df", lab=c(xlab,"Prediction Error"))       
      dev.off()
      
      # test error and true prediction error for each sample
      filename <- paste("sim-ex1-", j, "-", method[i], 
                        "-test-pe.pdf", sep = "")
      pdf(file = filename, width=9, height=7)
      plot.path(sim.tmp$path, xname="comp", ynames=c("Test","PE"),
                col=c("firebrick","forestgreen"), 
                title=title, logx=logx, topx=sim.tmp$pathave$df.cov, 
                toplab="df", lab=c(xlab,"Prediction Error"))
      dev.off()
      
      ### set y limits for sigma=3 and rho=0.5
      if (j==5) {
        # test error and training error for each sample
        filename <- paste("sim-ex1-5-train-", method[i], ".pdf", sep = "")
        pdf(file = filename, width=9, height=7)
        plot.path(sim.tmp$path, xname="comp", ynames=c("Train","Test"), 
                  col=c("steelblue","firebrick"), 
                  logx=logx, topx=sim.tmp$pathave$df.cov, 
                  toplab="df", lab=c(xlab,"Prediction Error"), ylim=c(5,25))
        dev.off()
        
        # test error and true prediction error for each sample
        filename <- paste("sim-ex1-5-pe-", method[i], ".pdf", sep = "")
        pdf(file = filename, width=9, height=7)
        plot.path(sim.tmp$path, xname="comp", ynames=c("Test","PE"),
                  col=c("firebrick","forestgreen"), 
                  logx=logx, topx=sim.tmp$pathave$df.cov, 
                  toplab="df", lab=c(xlab,"Prediction Error"), ylim=c(5,25))
        dev.off()
      }
    }
  }
}
rm(sim.tmp, filename, title, xlab, logx)

#save workspace
save.image("sim.ex1.RData")

# save data in Excel
library(RODBC)
xlsFile <- odbcConnectExcel("sim-ex1.xls", readOnly = FALSE)
sqlSave(xlsFile, best)
sqlSave(xlsFile, bestave)
sqlSave(xlsFile, pathave)
odbcCloseAll()
rm(xlsFile)
detach("package:RODBC", unload=TRUE)

# PLOTS
ic <- c("Valid","5.fold.CV", "10.fold.CV", "LOOCV", "GCV", "Cp", "BIC")

for (i in 1:3) {
  for (j in 1:nrow(vary)) {
    
    subrow <- quote(Method==toupper(method[i]) & 
                      sigma==vary$sigma[j] & rho==vary$rho[j])
    title <- substitute(paste(m, ": ", sigma == n, ", ", rho == o), 
                        list(m=toupper(method[i]), 
                             n=vary$sigma[j], o=vary$rho[j]))

    ## AVERAGES OVER PATH
    sub <- subset(pathave, eval(subrow))
    xlab <- switch(method[i], "forward"="p", "lasso"="s", 
                   "ridge"=expression(-log(lambda)))
    if (method[i] == "ridge") logx <- T else logx <- F
    
    # PREDICTION ERROR
    
    # average curves for true PE and all estimates
    filename <- paste("sim-ex1-", j, "-", method[i], "-pe-all.pdf", sep = "")
    pdf(file = filename, width=9, height=7)
    plot.lines(sub, xname="comp", ynames=c("PE","Test", ic,"Train"),
               col=c("forestgreen","firebrick","limegreen","black",
                     "darkmagenta","lightskyblue","darkorange","gray",
                     "goldenrod","steelblue"),
               lty=c(1,1,2,3,4,2,2,4,1,1), lab=c(xlab,"Prediction Error"), 
               logx=logx, topx=c("df.cov","df"), title=title)     
    dev.off()
    
    # PE Split
    filename <- paste("sim-ex1-", j, "-", method[i], "-mse.pdf", sep = "")
    pdf(file = filename, width=9, height=7)
    plot.lines(sub, xname="comp", ynames=c("Sq.Bias", "Variance", "MSE.Test"), 
               show.min=c(F,F,T), type="l", lwd=c(2,2,2), lty=c(4,2,1), 
               col=c("darkgreen","darkorchid","firebrick"),
               logx=logx, topx=c("df.cov","df"), 
               lab=c(xlab,"Mean Squared Error"), title=title)
    dev.off()

    ### set y limits for sigma=3 and rho=0.5
    if (j==5) {
      # average curves for true PE and all estimates
      filename <- paste("sim-ex1-5-ic-", method[i], ".pdf", sep = "")
      pdf(file = filename, width=9, height=7)
      plot.lines(sub, xname="comp", ynames=c("PE","Test", ic,"Train"),
                 col=c("forestgreen","firebrick","limegreen","black",
                       "darkmagenta","lightskyblue","darkorange","gray",
                       "goldenrod","steelblue"),
                 lty=c(1,1,2,3,4,2,2,4,1,1), lab=c(xlab,"Prediction Error"), 
                 logx=logx, topx=c("df.cov","df"), ylim=c(5,25))     
      dev.off()
      
      # PE Split
      filename <- paste("sim-ex1-5-mse-", method[i], ".pdf", sep = "")
      pdf(file = filename, width=9, height=7)
      plot.lines(sub, xname="comp", ynames=c("Sq.Bias", "Variance", "MSE.Test"), 
                 show.min=c(F,F,T), type="l", lwd=c(2,2,2), lty=c(4,2,1), 
                 col=c("darkgreen","darkorchid","firebrick"),
                 logx=logx, topx=c("df.cov","df"), 
                 lab=c(xlab,"Mean Squared Error"), ylim=c(0,10))
      dev.off()      
    }
    
    # COEFFICIENTS
    
    # coefficient averages
    filename <- paste("sim-ex1-", j, "-", method[i], "-beta-ave.pdf", sep = "")
    pdf(file = filename, width=7, height=7)
    plot.lines(sub, xname="comp", ynames=paste0("b.",1:p), 
               show.min=rep(F,p), type="l", lty=c(1,2,1,3,4,1,2,1),
               col=c("forestgreen","steelblue","limegreen","black",
                     "firebrick","lightskyblue","darkmagenta","goldenrod"),
               logx=logx, topx=c("df.cov","df"), 
               lab=c(xlab,"Coefficient Average"), title=title)
    dev.off()
    
    if (method[i] != "ridge") {
      # coefficient probabilities
      filename <- paste("sim-ex1-", j, "-", method[i], "-beta-prob.pdf", sep = "")
      pdf(file = filename, width=7, height=7)
      plot.lines(sub, xname="comp", ynames=paste0("prob.b.",1:p), 
                 show.min=rep(F,p), type="l", lty=c(1,2,1,3,4,1,2,1),
                 col=c("forestgreen","steelblue","limegreen","black",
                       "firebrick","lightskyblue","darkmagenta","goldenrod"),
                 logx=logx, topx=c("df.cov","df"), 
                 lab=c(xlab,"Coefficient Inclusion Probability"), title=title)
      dev.off()
      
      # nonzero coefficients
      filename <- paste("sim-ex1-", j, "-", method[i], "-beta-nonzero.pdf", sep = "")
      pdf(file = filename, width=7, height=7)
      plot.lines(sub, xname="comp", ynames=c("Corr.Nonzero","Inc.Nonzero"),
                 show.min=c(F,F), type="l", lty=c(2,1),
                 col=c("forestgreen","firebrick"),
                 logx=logx, topx=c("df.cov","df"), 
                 lab=c(xlab,"Inclusion Probability"), title=title)
      dev.off()
      
      # degrees of freedom
      filename <- paste("sim-ex1-", j, "-", method[i], "-beta-nonzero.pdf", sep = "")
      pdf(file = filename, width=7, height=7)
      plot.lines(sub, xname="Nonzero", ynames=c("Nonzero","df.cov"),
                 show.min=c(F,F), type="l", lty=c(2,1),
                 col=c("forestgreen","firebrick"),
                 lab=c(xlab,"Degrees of Freedom"), title=title)
      dev.off()
    }
    
    ## BEST MODELS
    sub <- subset(best, eval(subrow))
    
    # Information criteria comparison
    
    # Boxplots MSE
    if (method[i] == "ridge") {
      filename <- paste("sim-ex1-", j, "-", method[i], "-ic-mse.pdf", sep = "")
      pdf(file = filename, width=9, height=7)
      plot.box(sub, ynames="MSE", ylab="Mean Squared Error", 
               x1name="IC.Method", x1lab="", xlab.rot=T, title=title, outline=F)
      dev.off()
    } else {
      filename <- paste("sim-ex1-", j, "-", method[i], "-ic-mse.pdf", sep = "")
      pdf(file = filename, width=9, height=7)
      plot.box(sub, ynames="MSE", ylab="Mean Squared Error", 
               x1name="IC.Method", x1lab="", xlab.rot=T,
               zname="Corr.Subset", zcol="steelblue", title=title, outline=F)
      dev.off()
      
      # Bar plots subsets
      filename <- paste("sim-ex1-", j, "-", method[i], "-ic-select.pdf", sep = "")
      pdf(file = filename, width=9, height=7)
      plot.bar(sub, ynames=c("In.Path","Incl.Corr.Subset","Corr.Subset"), 
               xname="IC.Method", lab.rot=T, title=title, xlab="", 
               ycol=c("steelblue","firebrick","forestgreen"), hline=0.5, border=NA)
      dev.off()
      
      # Barplots nonzero parameters
      filename <- paste("sim-ex1-", j, "-", method[i], "-ic-nonzero.pdf", sep = "")
      pdf(file = filename, width=9, height=7)
      plot.bar(sub, ynames=c("Corr.Nonzero","Inc.Nonzero"), 
               xname="IC.Method", lab.rot=T, title=title, xlab="", 
               ycol=c("forestgreen","firebrick"), beside=F, hline=3, border=NA)
      dev.off()      
    }  
  }
}
rm(filename, title, subrow, sub, xlab, logx)

# MSE Split Multipanel
for (i in 1:3) {
  sub <- subset(pathave, Method==toupper(method[i]))
  if (method[i] == "ridge") logx <- T else logx <- F
  filename <- paste("sim-ex1-", method[i], "-mse.pdf", sep = "")
  pdf(file = filename, width=9, height=9)
  par(mfrow=c(3,3), oma=c(1.5,1.5,1.5,0))
  for (j in 1:nrow(vary)) {
    xaxt <- if (j %in% 7:9) "s" else "n"
    yaxt <- if (j %in% c(1,4,7)) "s" else "n"
    lim <- max(subset(sub,rho==vary$rho[j],select=42:44))
    plot.lines(subset(sub, sigma==vary$sigma[j] & rho==vary$rho[j]), 
               xname="comp", ynames=c("Sq.Bias", "Variance", "MSE.Test"), 
               show.min=c(F,F,T), type="l", lwd=c(2,2,2), lty=c(4,2,1), 
               col=c("darkgreen","darkorchid","firebrick"),
               topx= c("df.cov",""), ylim=c(0,lim), legend=F, 
               xaxt=xaxt, yaxt=yaxt, logx=logx)
  }
  title(toupper(method[i]), outer=T)
  mtext(expression(sigma==1), side=1, outer=T, at=1/6)
  mtext(expression(sigma==3), side=1, outer=T, at=3/6)
  mtext(expression(sigma==6), side=1, outer=T, at=5/6)
  mtext(expression(rho==0.75), side=2, outer=T, at=1/6)
  mtext(expression(rho==0.5), side=2, outer=T, at=3/6)
  mtext(expression(rho==0), side=2, outer=T, at=5/6)
  dev.off()
}
rm(i,j, lim, logx, xaxt, yaxt)

# BEST MODELS 5-fold CV & kappa - compare methods
beta.names = paste0("b.",1:p)
beta.prob <- paste0("prob.",beta.names)

for (j in 1:nrow(vary)) {
  subrow <- quote(IC.Method %in% c("5.fold.CV","kappa","Full","Oracle") & 
                    sigma==vary$sigma[j] & rho==vary$rho[j])
  title <- substitute(paste(sigma == n, ", ", rho == o), 
                      list(n=vary$sigma[j], o=vary$rho[j]))
  
  # Box plots MSE
  sub <- subset(best, eval(subrow))
  sub$Method2 <- as.factor(paste(sub$Method, sub$IC.Method, " "))
  levels(sub$Method2) <- gsub("."," ", levels(sub$Method2), fixed=T)
  filename <- paste("sim-ex1-", j, "-best-mse.pdf", sep = "")
  pdf(file = filename, width=9, height=7)
  plot.box(sub, ynames="MSE", ylab="Mean Squared Error", 
           x1name="Method2", x1lab="", xlab.rot=T,
           zname="Corr.Subset", zcol="steelblue", title=title, outline=F)
  dev.off()
  
  # Barplots subsets
  sub2 <- sub[sub$Method %in% c("FORWARD","LASSO"),]
  sub2$Method2 <- droplevels(sub2$Method2)
  filename <- paste("sim-ex1-", j, "-best-select.pdf", sep = "")
  pdf(file = filename, width=7, height=9)
  plot.bar(sub2, ynames=c("In.Path","Incl.Corr.Subset","Corr.Subset"), 
           xname="Method2", lab.rot=T, title=title, xlab="", 
           ycol=c("steelblue","firebrick","forestgreen"), hline=0.5, hlty=2,
           ylim=c(0,1), border=NA)
  dev.off()
  
  # Barplots nonzero parameters
  filename <- paste("sim-ex1-", j, "-best-nonzero.pdf", sep = "")
  pdf(file = filename, width=7, height=9)
  plot.bar(sub2, ynames=c("Corr.Nonzero","Inc.Nonzero"), 
           xname="Method2", lab.rot=T, title=title, xlab="", 
           ycol=c("forestgreen","firebrick"), beside=F, hline=3, hlty=2,
           ylim=c(0,6), border=NA)
  dev.off()   
  
  # BETAS
  lev <- levels(sub$Method2)
  for (i in 1:p) {
    title <- substitute(paste(sigma == n, ", ", rho == o), 
                        list(n=vary$sigma[j], o=vary$rho[j]))
      
    # Boxplots Betas
    filename <- paste("sim-ex1-", j, "-best-beta-box", i, ".pdf", sep = "")
    pdf(file = filename, width=9, height=7)
    blab  <- eval(substitute(paste0("beta[", m, "]"), list(m=i)))
    blab <- parse(text=blab)
    plot.box(sub, ynames=beta.names[i], ylab=blab, 
             x1name="Method2", x1lab="", xlab.rot=T,
             znames=beta.prob[i], zcol="steelblue", title=title, hline=beta[i])
    dev.off()
    
    # Histograms Betas

    # with density
    filename <- paste("sim-ex1-", j, "-best-beta-histdens-", i, ".pdf", sep = "")
    pdf(file = filename, width=12, height=6)
    par(mfrow=c(2,4), oma=c(0,0,3,0), mar=c(1.5,3,2,1)+0.1, 
        mgp=c(1.5,0.25,0), cex.axis=0.75, tcl=-0.25)
    for (k in lev[c(5,7,1,3, 6,8,2,4)]) {
      b <- sub[sub$Method2==k, beta.names[i]]
      plot.hist(b, ylab=blab, head=k, col="steelblue", norm=F, dens=T, breaks="Scott")
    }
    title(substitute(paste(sigma == n, ", ", rho == o), 
                           list(n=vary$sigma[j], o=vary$rho[j])), 
          outer=T, cex.main=1.5)
    dev.off()
    
    # with normal
    filename <- paste("sim-ex1-", j, "-best-beta-histnorm-", i, ".pdf", sep = "")
    pdf(file = filename, width=12, height=6)
    par(mfrow=c(2,4), oma=c(0,0,3,0), mar=c(1.5,3,2,1)+0.1, 
        mgp=c(1.5,0.25,0), cex.axis=0.75, tcl=-0.25)
    for (k in lev[c(5,7,1,3, 6,8,2,4)]) {
      b <- sub[sub$Method2==k, beta.names[i]]
      plot.hist(b, ylab=blab, head=k, col="steelblue", norm=T, dens=F, breaks="Scott")
    }
    title(substitute(paste(sigma == n, ", ", rho == o), 
                     list(n=vary$sigma[j], o=vary$rho[j])), 
          outer=T, cex.main=1.5)
    dev.off()
    
    # with both density and normal
    filename <- paste("sim-ex1-", j, "-best-beta-hist", i, ".pdf", sep = "")
    pdf(file = filename, width=12, height=6)
    par(mfrow=c(2,4), oma=c(0,0,3,0), mar=c(1.5,3,2,1)+0.1, 
        mgp=c(1.5,0.25,0), cex.axis=0.75, tcl=-0.25)
    for (k in lev[c(5,7,1,3, 6,8,2,4)]) {
      b <- sub[sub$Method2==k, beta.names[i]]
      plot.hist(b, ylab=blab, head=k, col="steelblue", norm=T, dens=T, breaks="Scott")
    }
    title(substitute(paste(sigma == n, ", ", rho == o), 
                     list(n=vary$sigma[j], o=vary$rho[j])), 
          outer=T, cex.main=1.5)
    dev.off()
    
    # QQ-plots betas
    filename <- paste("sim-ex1-", j, "-best-beta-qq-", i, ".pdf", sep = "")
    pdf(file = filename, width=12, height=6)
    par(mfrow=c(2,4), oma=c(0,0,3,0), mar=c(1.5,3,2,1)+0.1, 
        mgp=c(1.5,0.25,0), cex.axis=0.75, tcl=-0.25)
    for (k in lev[c(5,7,1,3, 6,8,2,4)]) {
      b <- sub[sub$Method2==k, beta.names[i]]
      qqnorm(b, ylab=blab, xlab=NULL, main=k, pch=16, col="steelblue")
      qqline(b)
    }
    title(substitute(paste(sigma == n, ", ", rho == o), 
                     list(n=vary$sigma[j], o=vary$rho[j])), 
          outer=T, cex.main=1.5)
    dev.off()
  }   
}


## TABLES
booktabs()

# max model Variance and Sq Bias
null <- which(pathave$Nonzero==0)
full <- which((pathave$Method=="FORWARD" & pathave$comp==p) |
                (pathave$Method=="LASSO" & pathave$comp==1) |
                (pathave$Method=="RIDGE" & pathave$comp==0) )
(tab <- rbind(tabular(Factor(Method)*Factor(rho, "$\\rho$")~
                        Heading()*max*Format(scientific=F)*
                        (Variance+Heading("Squared Bias")*Sq.Bias)*
                        Factor(sigma, "$\\sigma$"),
                      data=pathave[-c(null,full),]),
              tabular(Factor(IC.Method,"Method",c("OLS","ORACLE"))*
                        Factor(rho, "$\\rho$")~
                        Heading()*mean*Format(scientific=F)*
                        (Variance+Heading("Squared Bias")*Sq.Bias)*
                        Factor(sigma, "$\\sigma$"),
                      data=subset(best, Method=="OLS"))))
latex(tab, options=list(justification="n{2}{3}"))

# Model Variance, bias at min MSE
grp <- pathave[,c(1:3,44)]
grp$grp <- paste(grp[,1],grp[,2],grp[,3])
minMSE <- tapply(grp$MSE.Test, grp$grp, min)
minpos <- which(pathave$MSE.Test %in% minMSE)
(tab <- tabular(Factor(rho, "$\\rho$")*Factor(Method)~
                  Format(scientific=F)*Heading()*identity*
                  (Variance+Heading("Squared Bias")*Sq.Bias)*
                  Factor(sigma, "$\\sigma$"),
                data=pathave[minpos,]))
latex(tab,options=list(justification="n{2}{3}"))

# Nonzero variables and df at min MSE
(tab <- tabular(Factor(Method)*Factor(rho, "$\\rho$")~
                  Format(digits=3)*Heading()*identity*
                  (Heading("Number Variables ")*Nonzero+Heading("df")*df.cov)*
                  Factor(sigma, "$\\sigma$"),
                data=pathave[minpos,]))
latex(tab)

# Best models from cv and kappa

# median mse & se
(tab <- tabular(Factor(IC.Method,"",levelnames=list("5-fold CV","kappa"))*
                  Factor(rho, "$\\rho$")*Factor(Method) ~ 
                  Factor(sigma,"",levelnames=c("$\\sigma=1$","$\\sigma=3$","$\\sigma=6$"))
                  *Format(digits=2)*(Heading("Median MSE")*Median.MSE + 
                     Heading()*SE.Median.MSE+Format(digits=1)*
                     Heading("PCS")*Corr.Subset)*
                 Heading()*identity,
                 data=bestave[bestave$IC.Method %in% c("5.fold.CV","kappa"),]))
latex(tab)
(tab <- tabular(Factor(IC.Method,"",levelnames=list("OLS","Oracle"))*
                  Factor(rho, "$\\rho$")~ 
                  Factor(sigma,"",levelnames=c("$\\sigma=1$","$\\sigma=3$","$\\sigma=6$"))
                *Format(digits=2)*(Heading("Median MSE")*Median.MSE + 
                     Heading()*SE.Median.MSE+Heading("PCS")*Corr.Subset)*
                  Heading()*identity,
                data=bestave[bestave$Method== "OLS",]))
latex(tab)

# beta variance and bias
(tab <- tabular(Factor(IC.Method,"")*
                  Factor(rho, "$\\rho$")*Factor(Method) ~ 
                  Factor(sigma,"",levelnames=c("$\\sigma=1$","$\\sigma=3$","$\\sigma=6$"))*
                  (Heading("Variance")*Var.Betas+Heading("Squared Bias")*Sq.Bias.Betas)*
                  Heading()*identity,
                data=bestave[bestave$IC.Method %in% c("5.fold.CV","kappa"),]))
latex(tab,options=list(justification="n{2}{3}"))
(tab <- tabular(Factor(IC.Method,"",levelnames=list("OLS","Oracle"))*
                  Factor(rho, "$\\rho$")~ 
                  Factor(sigma,"",levelnames=c("$\\sigma=1$","$\\sigma=3$","$\\sigma=6$"))*
                  (Heading("Variance")*Var.Betas+Heading("Squared Bias")*Sq.Bias.Betas)*
                  Heading()*identity,
                data=bestave[bestave$Method== "OLS",]))
latex(tab)
