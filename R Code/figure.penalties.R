source("shrink.R")

method <- c("subset", "ridge", "lasso")
for (i in 1:length(method)) {
  filename <- paste0("fig-ortho-", method[i], ".pdf")
  pdf(file = filename, width=7, height=7)
  plot.shrink(type="threshold", method[i], lambda=2, lwd=2)
  dev.off()
  filename <- paste0("fig-pen-", method[i], ".pdf")
  pdf(file = filename, width=7, height=7)
  plot.shrink(type="penalty", method[i], lambda=2, lwd=2)
  dev.off()
}

method <- "bridge"
lambda <- c(1,2,4)
tune <- list(c(2,2.5,3,12), c(1, 1.4, 1.6,2), c(1,0.75,0.5,0.25))
col <- c("black","steelblue","maroon","olivedrab")
lty <- c(1,2,5,4)
lwd<- c(2,2,2,2)
ymax <- c(200, 80, 30)
title <- list(expression("Bridge Estimates: Convex with " * 
                           gamma %in% group("[", list(2,infinity), ")")),
              expression("Bridge Estimates: Convex with " * 
                           gamma %in% group("[", list(1,2), "]")),
              expression("Bridge Estimates: Concave with " * 
                           gamma %in% group("(", list(0,1), ")")))
for (i in 1:3) {
  filename <- paste0("fig-ortho-bridge", i, ".pdf")
  pdf(file = filename, width=7, height=7)
  plot.shrink(type="threshold", method, lambda=lambda[i], tune=tune[[i]], 
              title=title[[i]], col=col, lty=lty, lwd=lwd)
  dev.off()
  filename <- paste0("fig-pen-bridge", i, ".pdf")
  pdf(file = filename, width=7, height=7)
  plot.shrink(type="penalty", method, lambda=lambda[i], tune=tune[[i]], 
              title=title[[i]], col=col, lty=lty, lwd=lwd, ymax=ymax[i])
  dev.off()
}

method <- c("rlasso", "alasso")
lambda <- c(4,4)
tune <- list(c(1,0.75,0.5,0.25), c(0,0.5,1,2))
for (i in 1:2) {
  filename <- paste0("fig-ortho-", method[i], ".pdf")
  pdf(file = filename, width=8, height=7)
  plot.shrink(type="threshold", method[i], lambda=lambda[i], tune=tune[[i]], 
              col=col, lty=lty, lwd=lwd)
  dev.off()
  filename <- paste0("fig-pen-", method[i], ".pdf")
  pdf(file = filename, width=8, height=7)
  plot.shrink(type="penalty", method[i], lambda=lambda[i], tune=tune[[i]], 
              col=col, lty=lty, lwd=lwd)
  dev.off()
}

filename <- paste0("fig-ortho-nenet.pdf")
pdf(file = filename, width=8, height=7)
plot.shrink(type="threshold", method="nenet", lambda=2, 
            tune=c(0,0.3,0.6,1), lty=c(1,2,5,1), lwd=c(2,2,2,1),
            col=c("black","steelblue","maroon","black"))
dev.off()

filename <- paste0("fig-pen-nenet.pdf")
pdf(file = filename, width=8, height=7)
plot.shrink(type="penalty", method="nenet", lambda=2, 
            tune=c(0,0.3,0.6,1), lty=c(1,2,5,1), lwd=c(2,2,2,1),
            col=c("black","steelblue","maroon","black"), ymax=80)
dev.off()

filename <- paste0("fig-ortho-scad.pdf")
pdf(file = filename, width=7, height=7)
plot.shrink(type="threshold", method="scad", lambda=2, tune=3.7, lwd=2) 
dev.off()

filename <- paste0("fig-pen-scad.pdf")
pdf(file = filename, width=7, height=7)
plot.shrink(type="penalty", method="scad", lambda=2, tune=3.7, lwd=2) 
dev.off()

filename <- paste0("fig-ortho-mcp.pdf")
pdf(file = filename, width=7, height=7)
plot.shrink(type="threshold", method="mcp", lambda=2, tune=3, lwd=2) 
dev.off()

filename <- paste0("fig-pen-mcp.pdf")
pdf(file = filename, width=7, height=7)
plot.shrink(type="penalty", method="mcp", lambda=2, tune=3, lwd=2) 
dev.off()

