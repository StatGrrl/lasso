## DATA
source("article.sim.data.R")

beta=c(3, 1.5, 0, 0, 2, 0, 0, 0)
p <- length(beta)
n <- seq(25,500,by=25)
vary <- data.frame(n=n)
v <- nrow(vary)

data.art <- list()
for (i in 1:v) {
  data.art[[i]] <- sim.data(n.train=vary$n[i], n.test=0, beta=beta, seed=2020)
}
irr <- sapply(data.art, "[[", 10)
vary$irr.max <- apply(irr, 2, max)
vary$irr.false <- paste0(as.character(apply(irr, 2, function(x) sum(x>=1))),
                         "% False")
vary$irr.false[vary$irr.max<1] <- ""

# Save workspace
save.image("sim.data.RData")

## SIMULATIONS
source("article.sim.tune.R")

vary$PIP <- 0
for (i in 1:v) {
  print(paste(vary$n[i], Sys.time()), quote=F)
  sim.tmp <- sim.tune(data.art[[i]])
  vary$PIP[i] <- sim.tmp$PIP
  if (i==1) {
    best <- do.call("rbind", sim.tmp$best)
    bestave <- sim.tmp$bestave
  } else {
    best <- rbind(best, do.call("rbind", sim.tmp$best))
    bestave <- rbind(bestave, sim.tmp$bestave)
  }
}

# Save workspace
save.image("sim.results.RData")

# Save data in Excel
library(RODBC)
xlsFile <- odbcConnectExcel("sim-results.xls", readOnly = FALSE)
sqlSave(xlsFile, best)
sqlSave(xlsFile, bestave)
sqlSave(xlsFile, PIP)
odbcCloseAll()
detach("package:RODBC", unload=TRUE)

## PLOTS

# PLOT DATA
library(ggplot2)

# Conditions of data
cond <- rbind(data.frame(n=n, measure="Maximum Irrepresentable Condition", value=vary$irr.max, label=vary$irr.false),
              data.frame(n=n, measure="Probability Correct Subset Lies in LASSO Path", value=vary$PIP, label=""))

pdf(file = "conditions.pdf", width=7, height=7)
ggplot(data=cond,  aes(x=n, y=value)) +
  geom_point() + geom_line() + facet_wrap(~measure, nrow=2, scales="free") +
  geom_text(data=cond, mapping=aes(x=n, y=value, label=label, size=4, hjust=-0.1)) +
  geom_hline(yintercept=1, linetype="dashed") +
  labs(list(x="Sample Size", y="")) +
  theme_bw() + theme(legend.key = element_blank()) +
  theme(strip.background = element_rect(fill="white")) + theme(legend.position = "none")
dev.off()

# Consistency of MSE and PCS
names <- c("Validation.Set", "5.fold.CV", "10.fold.CV", "LOOCV",
           "GCV", "Cp", "AIC", "5.fold.CV.1SE", "10.fold.CV.1SE", "Percentile.CV", 
           "Kappa", "PASS", "BIC", "Modified.BIC")
bestave$Method.Type <- ifelse(bestave$Method %in% names[1:7],"Prediction Method", "Selection Method")
consist <- rbind(data.frame(bestave[,c(22,1,2)],Measure="Median MSE",Value=bestave$Median.MSE), 
                 data.frame(bestave[,c(22,1,2)],Measure="PCS",Value=bestave$PCS))

pdf(file = "consist.pdf", width=15, height=7)
ggplot(data=consist,  aes(x=n, y=Value, group=Method, shape=Method)) +
  geom_point(fill="black") + geom_line(aes(linetype=Method)) + facet_wrap(Measure~Method.Type, scales="free") +
  labs(list(x="Sample Size", y="")) +
  theme_bw() + theme(legend.key = element_blank()) +
  theme(strip.background = element_rect(fill="white")) +
  scale_linetype_discrete(breaks=names, labels=gsub("."," ", names, fixed=T)) +
  scale_shape_manual(values=c(21:25,3,4,21:25,3,4), breaks=names, labels=gsub("."," ", names, fixed=T))
dev.off()

# Split by n<=200, MSE, PCS and PIS, Nonzero Coefficients

best$Method <- rownames(best)
for (i in names)
  best$Method[grepl(i,rownames(best))] <- i
PCS <- bestave[,c(1,2,18)]
names(PCS)[3] <- "PCSave"
best1 <- merge(best,PCS,by=c("Method","n"))

# Boxplots MSE (without outliers)
source("article.geom_boxplot_noOutliers.R")
box <- ggplot(data=subset(best1, n==25|n==200), aes(x=Method, y=MSE, fill=PCSave)) + 
  geom_boxplot_noOutliers() + facet_wrap(~n) +
  scale_x_discrete(name="", limits=names, labels=gsub("."," ", names, fixed=T)) +
  ylab("Mean Squared Error") + 
  scale_fill_gradient(low="white", high="grey20", name="PCS") +
  theme_bw() + theme(axis.text.x=element_blank()) +
  theme(strip.background = element_rect(fill="white"))

# Barplots probabilities
library(reshape2)
prob <- melt(subset(bestave, n==25|n==200, select=c(1,2,18,19)),id.var=c("Method","n"))
bar.prob <- ggplot(data=prob, aes(x=Method, y=value, fill=variable))+ facet_wrap(~n) +
  geom_bar(stat="identity", position="dodge") + geom_hline(yintercept=0.5, linetype="dashed") +
  scale_x_discrete(name="", limits=names, labels=gsub("."," ", names, fixed=T)) +
  ylab("Probability") + 
  scale_fill_manual(name="Correct Subset", labels=c("Chosen (PCS)","Included (PIS)"), values=c("grey70","grey80")) +
  theme_bw() + theme(axis.text.x=element_blank()) +
  theme(strip.background = element_rect(fill="white"))

# Barplots number of nonzero variables
nonzero <- melt(subset(bestave, n==25|n==200, select=c(1,2,14,15)),id.var=c("Method","n"))
bar.nonzero <- ggplot(data=nonzero, aes(x=Method, y=value, fill=variable)) + facet_wrap(~n) +
  geom_bar(stat="identity") + geom_hline(yintercept=3, linetype="dashed") +
  scale_x_discrete(name="", limits=names, labels=gsub("."," ", names, fixed=T)) +
  ylab("Number of Nonzero Coefficients") + 
  scale_fill_manual(name="True Value", labels=c("Nonzero","Zero"), values=c("grey70","grey80")) +
  theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1)) +
  theme(strip.background = element_rect(fill="white"))

library(grid)
pdf(file="results.pdf", width=15, height=10.5)
grid.draw(rbind(ggplotGrob(bar.prob), ggplotGrob(box), ggplotGrob(bar.nonzero), size="first"))
dev.off()

## TABLES

library(tables)
booktabs()

(tab <- as.tabular(vary[,-3]))
latex(tab)

# Consistency of MSE, PCS, PIS

(tab <- tabular(Factor(n)~Factor(Method,"",levelnames=gsub("."," ", names, fixed=T))*
                  Heading()*Median.MSE*Heading()*identity, data=bestave))
latex(tab)

(tab <- tabular(Factor(n)~Factor(Method,"",levelnames=gsub("."," ", names, fixed=T))*
                  Heading()*PCS*Heading()*identity, data=bestave))
latex(tab)

(tab <- tabular(Factor(n)~Factor(Method,"",levelnames=gsub("."," ", names, fixed=T))*
                  Heading()*PIS*Heading()*identity, data=bestave))
latex(tab)

# Median MSE & SE, PCS, PIS

(tab <- tabular(Factor(n)*Factor(Method, levelnames=gsub("."," ", names, fixed=T))~
                  Format(format(digits=1))*
                  (Heading("Median MSE")*Median.MSE + 
                     Heading("Bootstrap SE")*SE.Median.MSE + 
                     PCS + PIS)*Heading()*identity, 
                data=subset(bestave,n==25)))
latex(tab)

save.image("sim.results2.RData")
