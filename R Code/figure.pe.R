source("sim.data.R")
source("pe.path.R")
source("plots.R")

beta <- rep(0,45)
beta[c(1:5,11:15,21:25)] <- c(rep(2,15))
data.pe <- sim.data(n.train=60, n.test=60, sigma=4, 
                  rho=0.25, beta=beta, beta0=0, seed=100)

lasso.pe <- pe.path(method="lasso", data=data.pe)

pdf(file="fig-pe-lasso1.pdf", width=9, height=7)
plot.path(lasso.pe$path, xname="comp", ynames=c("Train", "Test"), 
          topx=lasso.pe$mse$df, toplab="df",
          col=c("deepskyblue","tomato"), lab=c("s","Prediction Error"))
dev.off()

pdf(file="fig-pe-lasso2.pdf", width=9, height=7)
plot.lines(lasso.pe$mse, xname="comp", ynames=c("Sq Bias", "Variance", "MSE"), 
           show.min=c(F,F,T), type="l", lwd=c(2,2,2), lty=c(4,2,1), 
           col=c("darkgreen","darkorchid", "tomato"), 
           lab=c("s","Mean Squared Error"), topx=c("df","df"))
dev.off()

save.image("figure.pe.RData")

plot.lines(lasso.pe$mse, xname="df.aprx", ynames=c("df","df.aprx"),
           show.min=c(F,F), type="l", lty=c(2,1),
           col=c("forestgreen","firebrick"))