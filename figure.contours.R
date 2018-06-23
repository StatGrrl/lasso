source("sim.data.R")
source("shrink.R")
beta <- c(5, 2)
rho <- c(-0.9, -0.75, -0.5, 0.5, 0.75, 0.9)
method <- c("ridge", "lasso")
sd.steps <- c(18,22,36,26,22,18)
for (i in seq_along(rho)) {
  data <- sim.data(iter=1, n.train=100, sigma=3, rho=rho[i], beta=beta, beta0=0, seed=111)
  X <- data$X.train[[1]]
  y <- data$y.train[[1]]
  #cor(data.frame(X,y))
  for (j in seq_along(method)) {
    method.name <- switch(method[j], "ridge"="Ridge Regression", "lasso"="LASSO")
    title <- substitute(paste(m, ": ", rho==n), list(m=method.name, n=rho[i]))
    filename <- paste("figure-geom-", method[j], i, ".pdf", sep = "")
    pdf(file = filename, width=9, height=7)
    contours.rss(X, y, method[j], sd.step=sd.steps[i], title=title)
    dev.off()
  }
}
