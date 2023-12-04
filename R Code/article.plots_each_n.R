load("sim.results.RData")
library(ggplot2)
library(reshape2)

names <- c("Validation.Set", "5.fold.CV", "10.fold.CV", "LOOCV",
           "GCV", "Cp", "AIC", "5.fold.CV.1SE", "10.fold.CV.1SE", "Percentile.CV", 
           "Kappa", "PASS", "BIC", "Modified.BIC")
best$Method <- rownames(best)
for (i in names)
  best$Method[grepl(i,rownames(best))] <- i
PCS <- bestave[,c(1,2,18)]
names(PCS)[3] <- "PCSave"
best1 <- merge(best,PCS,by=c("Method","n"))

for (i in 1:v) {
  # Boxplots MSE
  filename <- paste0("n", vary$n[i], "-box-mse.pdf")
  pdf(file=filename, width=9, height=7)
  print(ggplot(data=subset(best1, n==vary$n[i]), aes(x=Method, y=MSE)) + 
          geom_boxplot(aes(fill=PCSave)) +
          scale_x_discrete(name="", limits=names, labels=gsub("."," ", names, fixed=T)) +
          ylab("Mean Squared Error") + ggtitle(substitute(n == a, list(a=vary$n[i]))) + 
          scale_fill_gradient(low="white", high="grey20", name="PCS") +
          theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1)))
  dev.off()
  
  # Barplots number of nonzero variables
  nonzero <- melt(subset(bestave, n==vary$n[i], select=c(1,14,15)),id.var="Method")
  filename <- paste0("n", vary$n[i], "-bar-nonzero.pdf")
  pdf(file=filename, width=9, height=7)
  print(ggplot(data=nonzero, aes(x=Method, y=value, fill=variable)) + 
          geom_bar(stat="identity") + geom_hline(yintercept=3, linetype="dashed") +
          scale_x_discrete(name="", limits=names, labels=gsub("."," ", names, fixed=T)) +
          ylab("Number of Nonzero Coefficients") + ggtitle(substitute(n == a, list(a=vary$n[i]))) + 
          scale_fill_manual(name="True Value", labels=c("Nonzero","Zero"), values=c("grey70","grey80")) +
          theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1)))
  dev.off()
  
  # Barplots probabilities
  prob <- melt(subset(bestave, n==vary$n[i], select=c(1,18,19)),id.var="Method")
  filename <- paste0("n", vary$n[i], "-bar-prob.pdf")
  pdf(file=filename, width=9, height=7)
  print(ggplot(data=prob, aes(x=Method, y=value, fill=variable)) + 
          geom_bar(stat="identity", position="dodge") + geom_hline(yintercept=0.5, linetype="dashed") +
          scale_x_discrete(name="", limits=names, labels=gsub("."," ", names, fixed=T)) +
          ylab("Probability") + ggtitle(substitute(n == a, list(a=vary$n[i]))) + 
          scale_fill_manual(name="Correct Subset", labels=c("Chosen","Included"), values=c("grey70","grey80")) +
          theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1)))
  dev.off()
}
