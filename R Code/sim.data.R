sim.data <- function(iter=100, n.train=100, n.test=0, n.valid=0, 
                     sigma=1, corr="AR", rho=0, beta, beta0=0, 
                     seed=floor(runif(1,1,10000))) {

  # Generates all X ~ N(0,1), y ~ N(beta0 + X %*% beta, sigma^2)

  # INPUT:
  # iter......number iterations
  #           default=100
  # n.train...size training samples
  #           default=100
  # n.test....size test sample
  #           default=0 - no test data
  # n.valid...size validation sample
  #           default=0 - no validation data
  # sigma.....sqrt of error variance
  #           default=1 - standard normal
  # corr......correlation structure of X, cor(Xi, Xi) = 1
  #           AR - cor(Xi, Xj) = rho^abs(i-j)
  #           CS - cor(Xi, Xj) = rho
  #           IR - cor(Xi, XA) = rho for beta[i]==0 & beta[A]!=0
  # rho.......correlation coefficient
  #           default=0 - orthogonal
  # beta......parameter vector
  # beta0.....intercept term
  # seed......seed for random number generation
  #           default=floor(runif(1,1,10000))


  require(MASS)
  set.seed(seed)
  
  # X STRUCTURE
  p <- length(beta)
  D <- which(beta != 0)
  Dc <- (1:p)[-D]
  d <- length(D)
  X.names <- paste0("X",1:p)
  X.mu <- rep(0,p)
  if (corr=="AR") X.sigma <- toeplitz(rho^(0:(p-1)))
  if (corr=="CS") {
    X.sigma <- matrix(rho, p, p)
    diag(X.sigma) <- 1
  }
  if (corr=="IR") {
    X.sigma <- diag(p)
    X.sigma[D,Dc[1]] <- X.sigma[Dc[1],D] <- rho
  }
    
  # compatibility condition
  quad.eig <- t(beta) %*% X.sigma %*% beta
  compat <- drop(d*quad.eig/(sum(abs(beta[D])))^2) # want > 1
  # restricted eigenvalue condition
  resteig <- drop(quad.eig/sum(beta[D]^2)) # want > 1

  X.train <- list()
  y.train <- list()
  snr <- rep(0,iter)
  cond <- rep(0,iter)
  irr <- rep(0,iter)

  for (i in 1:iter) {
    
    # TRAINING DATA
    X.train[[i]] <- mvrnorm(n.train, X.mu, X.sigma)
    colnames(X.train[[i]]) <- X.names
    y.train[[i]] <- beta0 + X.train[[i]] %*% beta + rnorm(n.train, 0, sigma)

    # signal to noise ratio
    snr[i] <- sqrt(sum((X.train[[i]] %*% beta)^2)) / sigma^2
    # condition number
    eigval <- eigen(cor(X.train[[i]]), only.values=T)$values
    cond[i] <- sqrt(max(eigval)/min(eigval))
    # irrepresentable condition
    irr[i] <- max(t(X.train[[i]][,Dc])%*%X.train[[i]][,D] %*%
                      solve(t(X.train[[i]][,D])%*%X.train[[i]][,D]) %*% 
                      sign(beta[D])) # want < 1
  }

  # TEST DATA
  if (n.test > 0) {
    X.test <- mvrnorm(n.test, X.mu, X.sigma)
    colnames(X.test) <- X.names
    y.test <- beta0 + X.test %*% beta + rnorm(n.test, 0, sigma)
  }
  
  # VALIDATION DATA
  if (n.valid > 0) {
    X.valid <- mvrnorm(n.valid, X.mu, X.sigma)
    colnames(X.valid) <- X.names
    y.valid <- beta0 + X.valid %*% beta + rnorm(n.valid, 0, sigma)
  }

  # OUTPUT:
  data <- list(beta0=beta0, beta=beta, rho=rho, sigma=sigma, X.sigma=X.sigma, 
               corr=corr, snr=snr, cond=cond, compat=compat, resteig=resteig,
               irr=irr, X.train=X.train, y.train=y.train)
  if (n.test > 0) {
    data$X.test <- X.test
    data$y.test <- y.test
  }
  if (n.valid > 0) {
    data$X.valid <- X.valid
    data$y.valid <- y.valid
  }
  data
}
