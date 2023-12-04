source("shrink.R")

method <- c("ridge", "lasso")
for (i in 1:2) for (j in 2:3) {
  filename <- paste0("fig-norm-", method[i], "-", j, "D.pdf")
  pdf(file = filename, width=7, height=7)
  plot.normball(method[i], dim=j)
  dev.off()
}

method <- "bridge"
tune <- c(0.5, 1, 1.5, 2, 3, Inf)
for (i in tune) for (j in 2:3) {
  filename <- paste0("fig-norm-", method, "-", i, "-", j, "D.pdf")
  pdf(file = filename, width=7, height=7)
  plot.normball(method, dim=j, i)
  dev.off()
}

method <- c("enet", "oscar", "flasso")
for (i in 1:3) for (j in 2:3) {
  filename <- paste0("fig-norm-", method[i], "-", j, "D.pdf")
  pdf(file = filename, width=7, height=7)
  plot.normball(method[i], dim=j, tune=0.5)
  dev.off()
}

method <- c("glasso", "cap", "gbridge")
tune <- c(0, Inf, 0.5)
for (i in seq_along(method)) {
  filename <- paste0("fig-norm-", method[i], "-", j, "D.pdf")
  pdf(file = filename, width=7, height=7)
  plot.normball(method[i], dim=3, tune=tune[i])
  dev.off()
}
