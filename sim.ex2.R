library(tables)
booktabs()
source("sim.data.R")
source("sim.oracle.R")
source("plots.R")
source("color.R")

# DATA

beta=c(3, 1.5, 0, 0, 2, 0, 0, 0)
p <- length(beta)
n <- seq(25,500,by=25)
corr <- c("AR","CS","IR") 
vary <- expand.grid(n=n, corr=corr)
v <- nrow(vary)

data.ex2 <- list()
for (i in 1:v) {
  data.ex2[[i]] <- sim.data(iter=100, n.train=vary$n[i], 
                            corr=as.character(vary$corr[i]),
                            sigma=3, rho=0.5, beta=beta, beta0=0, 
                            seed=2020)
}
snr <- sapply(data.ex2, "[[", 7)
cond <- sapply(data.ex2, "[[", 8)
irr <- sapply(data.ex2, "[[", 11)
vary$snr <- colMeans(snr)
vary$cond <- colMeans(cond)
vary$compat <- sapply(data.ex2, "[[", 9)
vary$resteig <- sapply(data.ex2, "[[", 10)
vary$irr.max <- apply(irr, 2, max)
vary$irr.false <- apply(irr, 2, function(x) sum(x>=1)/100)

# plot average snr
snr.ave <- data.frame(n=n, AR=vary$snr[vary$corr=="AR"],
                     CS=vary$snr[vary$corr=="CS"],
                     IR=vary$snr[vary$corr=="IR"])
pdf(file = "sim-ex2-data-snr.pdf", width=9, height=7)
plot.lines(snr.ave, "n", corr, col=c("steelblue","firebrick","forestgreen"), 
           show.min=rep(F,3), type="b", pch=c(15,16,17), 
           lab=c("n","Average Signal to Noise Ratio"))
dev.off()

# plot average condition number
cond.ave <- data.frame(n=n, AR=vary$cond[vary$corr=="AR"],
                      CS=vary$cond[vary$corr=="CS"],
                      IR=vary$cond[vary$corr=="IR"])
pdf(file = "sim-ex2-data-cond.pdf", width=9, height=7)
plot.lines(cond.ave, "n", corr, col=c("steelblue","firebrick","forestgreen"), 
           show.min=rep(F,3), type="b", pch=c(15,16,17), 
           lab=c("n","Average Condition Number"))
dev.off()

# plot irrepresentable condition
# maximum value
irr.max <- data.frame(n=n, AR=vary$irr.max[vary$corr=="AR"],
                       CS=vary$irr.max[vary$corr=="CS"],
                       IR=vary$irr.max[vary$corr=="IR"])
pdf(file = "sim-ex2-data-irr-line.pdf", width=9, height=7)
plot.lines(irr.max, "n", corr, col=c("steelblue","firebrick","forestgreen"), 
           show.min=rep(F,3), type="b", pch=c(15,16,17), 
           lab=c("n","Maximum Irrepresentable Condition"), hline=1, hlty=2)
dev.off()

# probability that irr cond is false
irr.false <- data.frame(n=n, AR=vary$irr.false[vary$corr=="AR"],
                        CS=vary$irr.false[vary$corr=="CS"],
                        IR=vary$irr.false[vary$corr=="IR"])
# adjacent bars
pdf(file = "sim-ex2-data-irr-bar1.pdf", width=9, height=7)
plot.bar(irr.false, corr, ylab="P(Irrepresentable Condition < 1)", 
         ycol=c("steelblue","firebrick","forestgreen"),
         xname="n", xlab="n", labels=ifelse(n%%100==0, as.character(n), ""),
         border=NA)
dev.off()

# overlay bars
pdf(file = "sim-ex2-data-irr-bar2.pdf", width=9, height=7)
plot.bar(irr.false, rev(corr), ylab="P(Irrepresentable Condition < 1)", 
         ycol=c("forestgreen","firebrick","steelblue"),
         xname="n", xlab="n", labels=ifelse(n%%100==0, as.character(n), ""),
         border=NA, overlay=T)
dev.off()

# table compat and rest eigen conditions
tab <- tabular( Factor(corr, "Correlation Structure")~ 
                  Heading("Condition")*mean*
                  (Heading("Compatibility")*compat +
                  Heading("Restricted Eigenvalue")* resteig), data=vary)
latex(tab)

# SIMULATIONS

method = c("lasso","alasso","slasso","tlasso","enet","aenet","mcp","scad","rlasso", "ols")

for (i in seq_along(method)) {
  for (j in 1:nrow(vary)) {
    
    # simulations
    print(paste(method[i], j, Sys.time()), quote=F)
    sim.tmp <- sim.oracle(method=method[i],data=data.ex2[[j]])
    if (i+j == 2) {
      best <- do.call("rbind", sim.tmp$best)
      bestave <- sim.tmp$bestave
    } else {
      best <- rbind(best, do.call("rbind", sim.tmp$best))
      bestave <- rbind(bestave, sim.tmp$bestave)
    }
  }
}

#save workspace
save.image("sim.ex2.RData")
rm(sim.tmp)

# save data in Excel
lasso.meth <- c("LASSO", "ALASSO", "SLASSO", "TLASSO", "RLASSO")
best1 <- subset(best, Method %in% lasso.meth)
best2 <- subset(best, !Method %in% lasso.meth)
library(RODBC)
xlsFile <- odbcConnectExcel("sim-ex2.xls", readOnly = FALSE)
sqlSave(xlsFile, best1)
sqlSave(xlsFile, best2)
sqlSave(xlsFile, bestave)
odbcCloseAll()
rm(xlsFile, best1, best2)
detach("package:RODBC", unload=TRUE)

# PLOTS

# consistency

METH1 <- toupper(method)[!method %in% c("ols","tlasso","aenet")]
METH2 <- c(METH1, c("TLASSO","AENET"))
col <- c("forestgreen", "limegreen","goldenrod","darkmagenta",
         "steelblue","lightskyblue","darkorange",
         "firebrick","magenta")

for (i in corr) {
  
  sub <- subset(bestave, corr==i & IC.Method=="cv")
  colnames <- METH2
  
  # probability in path
  ip <- data.frame(n=n)
  for (k in colnames) {
    ip <- cbind(ip, sub$"In.Path"[sub$Method==k])
  }
  names(ip) <- c("n", colnames)
  filename <- paste0("sim-ex2-ip-", i, ".pdf")
  pdf(file = filename, width=9, height=7)
  plot.lines(ip, "n", colnames, show.min=rep(F,9), title=i,
             col=col, lty=c(1,1,1,3,5,5,1,1,3), 
             lab=c("n","Probability In Path"))
  dev.off()
  
  for (j in c("cv","kappa")) {
    
    sub <- subset(bestave, corr==i & IC.Method==j)
    title <- substitute(paste(m, ": ", n), list(m=i, n=j))
    if (j=="kappa") colnames <- METH1 else 
      colnames <- METH2

    # probability correct subset
    pcs <- data.frame(n=n)
    for (k in colnames) {
      if (!(j=="kappa" & k %in% c("TLASSO","AENET")))
        pcs <- cbind(pcs, sub$"Corr.Subset"[sub$Method==k])
    }
    names(pcs) <- c("n", colnames)
    filename <- paste0("sim-ex2-pcs-", i, "-", j, ".pdf")
    pdf(file = filename, width=9, height=7)
    plot.lines(pcs, "n", colnames, show.min=rep(F,9), title=title,
               col=col, lty=c(1,1,1,3,5,5,1,1,3), 
               lab=c("n","Probability Selecting Correct Subset"))
    dev.off()

    # probability include correct subset
    ics <- data.frame(n=n)
    for (k in colnames) {
      ics <- cbind(ics, sub$"Incl.Corr.Subset"[sub$Method==k])
    }
    names(ics) <- c("n", colnames)
    filename <- paste0("sim-ex2-ics-", i, "-", j, ".pdf")
    pdf(file = filename, width=9, height=7)
    plot.lines(ics, "n", colnames, show.min=rep(F,9), title=title,
               col=col, lty=c(1,1,1,3,5,5,1,1,3), 
               lab=c("n","Probability Including Correct Subset"))
    dev.off()
  
    # mse estimates
    mse.est <- data.frame(n=n)
    for (k in colnames) {
      mse.est <- cbind(mse.est, sub$"MSE.Betas"[sub$Method==k])
    }
    names(mse.est) <- c("n", colnames)
    filename <- paste0("sim-ex2-mse-est-", i, "-", j, ".pdf")
    pdf(file = filename, width=9, height=7)
    plot.lines(mse.est, "n", colnames, show.min=rep(F,9), title=title,
               col=col, lty=c(1,1,1,3,5,5,1,1,3), 
               lab=c("n","MSE Estimates"))
    dev.off()
    
    # median mse predictions
    mse.pred <- data.frame(n=n)
    for (k in colnames) {
      mse.pred <- cbind(mse.pred, sub$"Median.MSE"[sub$Method==k])
    }
    names(mse.pred) <- c("n", colnames)
    filename <- paste0("sim-ex2-mse-pred-", i, "-", j, ".pdf")
    pdf(file = filename, width=9, height=7)
    plot.lines(mse.pred, "n", colnames, show.min=rep(F,9), title=title,
               col=col, lty=c(1,1,1,3,5,5,1,1,3), 
               lab=c("n","Median MSE Predictions"))
    dev.off()
    
  }
}

# Small sample comparison with sim ex1
best.small <- subset(best, n==25 & corr=="AR")
best.small$Method2 <- as.factor(paste(best.small$Method, best.small$IC.Method, " "))
levels(best.small$Method2) <- gsub("."," ", levels(best.small$Method2), fixed=T)
bestave.small <- subset(bestave, n==25 & corr=="AR")
bestave.small$Method2 <- as.factor(paste(bestave.small$Method, bestave.small$IC.Method, " "))
levels(bestave.small$Method2) <- gsub("."," ", levels(bestave.small$Method2), fixed=T)
title <- expression("AR: "* italic(n)==25)

# PREDICTION

# Box plots MSE
pdf(file = "sim-ex2-small-mse.pdf", width=11, height=7)
plot.box(best.small, ynames="MSE", ylab="Mean Squared Error", 
         x1name="Method2", x1lab="", xlab.rot=T,  
         zname="Corr.Subset", zcol="steelblue", title=title, outline=F)
dev.off()

# table median mse & se
(tab <- tabular(Factor(Method) ~ Factor(IC.Method,"",levelnames=list("10-fold CV","kappa"))*
                  Format(digits=2)*(Heading("Median MSE")*Median.MSE + 
                                      Heading()*SE.Median.MSE+Format(digits=1)*
                                      Heading("PCS")*Corr.Subset)*Heading()*identity,
                data=bestave.small[bestave.small$Method!="OLS",]))
tab[c(1,9),4:6] <- ""
latex(tab)

(tab <- tabular(Factor(Method) ~ Factor(IC.Method,"",levelnames=list("OLS","Oracle"))*
                  Format(digits=2)*(Heading("Median MSE")*Median.MSE + 
                                      Heading()*SE.Median.MSE+Format(digits=1)*
                                      Heading("PCS")*Corr.Subset)*Heading()*identity,
                data=bestave.small[bestave.small$Method=="OLS",]))
latex(tab)

# ESTIMATION

# betas
beta.names = paste0("b.",1:p)
beta.prob <- paste0("prob.",beta.names)
lev <- levels(best.small$Method2)
for (i in 1:p) {
  
  blab  <- eval(substitute(paste0("beta[", m, "]"), list(m=i)))
  blab <- parse(text=blab)

  # Boxplots Betas
  filename <- paste("sim-ex2-small-beta-box", i, ".pdf", sep = "")
  pdf(file = filename, width=9, height=7)
  plot.box(best.small, ynames=beta.names[i], ylab=blab, 
           x1name="Method2", x1lab="", xlab.rot=T,
           znames=beta.prob[i], zcol="steelblue", hline=beta[i])
  dev.off()
  
  # Histograms Betas
  filename <- paste("sim-ex2-small-beta-histnorm-", i, ".pdf", sep = "")
  pdf(file = filename, width=12, height=6)
  par(mfrow=c(2,4), mar=c(1.5,3,2,1)+0.1, 
      mgp=c(1.5,0.25,0), cex.axis=0.75, tcl=-0.25)
  for (k in lev[c(12,2,4,14, 13,3,5,15)]) {
    b <- best.small[best.small$Method2==k, beta.names[i]]
    plot.hist(b, ylab=blab, head=k, col="steelblue", norm=T, dens=F, breaks="Scott")
  }
  dev.off()

  filename <- paste("sim-ex2-small-beta-histnorm-", i, "b.pdf", sep = "")
  pdf(file = filename, width=12, height=6)
  par(mfrow=c(2,4), mar=c(1.5,3,2,1)+0.1, 
      mgp=c(1.5,0.25,0), cex.axis=0.75, tcl=-0.25)
  for (k in lev[c(1,8,16,6, 18,9,17,7)]) {
    b <- best.small[best.small$Method2==k, beta.names[i]]
    plot.hist(b, ylab=blab, head=k, col="steelblue", norm=T, dens=F, breaks="Scott")
  }
  dev.off()
}

# bias and variance of estimates
(tab <- tabular(Factor(Method)~Factor(IC.Method,"",levelnames=list("10-fold CV","kappa"))*
                  (Heading("Variance")*Var.Betas+Heading("Squared Bias")*Sq.Bias.Betas)*
                  Heading()*identity,
                data=bestave.small[bestave.small$Method != "OLS",]))
tab[c(1,9),3:4] <- ""
latex(tab)

(tab <- tabular(Factor(Method)~Factor(IC.Method,"",levelnames=list("OLS","Oracle"))*
                  (Heading("Variance")*Var.Betas+Heading("Squared Bias")*Sq.Bias.Betas)*
                  Heading()*identity,
                data=bestave.small[bestave.small$Method == "OLS",]))
latex(tab)

# SELECTION

# Barplots subsets
sub <- best.small[best.small$Method != "OLS",]
sub$Method2 <- droplevels(sub$Method2)
pdf(file = "sim-ex2-small-select.pdf", width=11, height=7)
plot.bar(sub, ynames=c("In.Path","Incl.Corr.Subset","Corr.Subset"), 
         xname="Method2", lab.rot=T, title=title, xlab="", 
         ycol=c("steelblue","firebrick","forestgreen"), hline=0.5, 
         ylim=c(0,1), border=NA)
dev.off()

# Barplots nonzero parameters
pdf(file = "sim-ex2-small-nonzero.pdf", width=11, height=7)
plot.bar(sub, ynames=c("Corr.Nonzero","Inc.Nonzero"), 
         xname="Method2", lab.rot=T, title=title, xlab="", 
         ycol=c("forestgreen","firebrick"), beside=F, hline=3, 
         ylim=c(0,6), border=NA)
dev.off()   


