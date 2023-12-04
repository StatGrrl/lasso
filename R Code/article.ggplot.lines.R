ggplot.lines <- function(data, x, y, xlab, ylab, main, type) {
  require(ggplot2)
  if (type=="pred") {
    names <- c("Validation.Set", "5.fold.CV", "10.fold.CV", "LOOCV",
               "GCV", "Cp", "AIC")
    main <- "Prediction Methods"
  }
  if (type=="select") {
    names <- c("5.fold.CV.1SE", "10.fold.CV.1SE", "Percentile.CV", 
               "Kappa", "PASS", "BIC", "Modified.BIC")
    main <- "Selection Methods"
  }
  subset <- subset(data, Method %in% names)
  
  ggplot(data=subset,  aes_string(x=x, y=y, group="Method", shape="Method")) +
    geom_point(fill="black") + geom_line(aes(linetype=Method)) +  
    labs(list(title=main, x=xlab, y=ylab)) +
    theme_bw() + theme(legend.key = element_blank()) +
    scale_linetype_discrete(breaks=names, labels=gsub("."," ", names, fixed=T)) +
    scale_shape_manual(values=c(21:25,3,4), breaks=names, labels=gsub("."," ", names, fixed=T))
}

load("sim.results.RData")

# Median MSE
pdf(file = "cons-mse-pred.pdf", width=9, height=7)
ggplot.lines(bestave, x="n", y="Median.MSE", "Sample Size", "Median Mean Squared Error", type="pred")
dev.off()

pdf(file = "cons-mse-select.pdf", width=9, height=7)
ggplot.lines(bestave, x="n", y="Median.MSE", "Sample Size", "Median Mean Squared Error", type="select")
dev.off()

# PCS
pdf(file = "cons-pcs-pred.pdf", width=9, height=7)
ggplot.lines(bestave, x="n", y="PCS", "Sample Size", "Probability of Choosing the Correct Subset", type="pred")
dev.off()

pdf(file = "cons-pcs-select.pdf", width=9, height=7)
ggplot.lines(bestave, x="n", y="PCS", "Sample Size", "Probability of Choosing the Correct Subset", type="select")
dev.off()

# PIS
pdf(file = "cons-pis-pred.pdf", width=9, height=7)
ggplot.lines(bestave, x="n", y="PIS", "Sample Size", "Probability of Including the Correct Subset", type="pred")
dev.off()

pdf(file = "cons-pis-select.pdf", width=9, height=7)
ggplot.lines(bestave, x="n", y="PIS", "Sample Size", "Probability of Including the Correct Subset", type="select")
dev.off()
