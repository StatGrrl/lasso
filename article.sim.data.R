sim.data <- function(iter=100, n.train=100, n.test=1000, n.valid=500, 
                     sigma=3, corr="AR", rho=0.5, beta,
                     seed=floor(runif(1,1,10000))) {

  # Generates all X ~ N(0,1), y ~ N(X %*% beta, sigma^2)

  # INPUT:
  # iter......number iterations
  #           default=100
  # n.train...size training samples
  #           default=100
  # n.test....size test sample
  #           default=1000
  # n.valid...size validation sample
  #           default=500
  # sigma.....sqrt of error variance
  #           default=3
  # corr......correlation structure of X, cor(Xi, Xi) = 1
  #           AR - cor(Xi, Xj) = rho^abs(i-j)
  #           CS - cor(Xi, Xj) = rho
  #           IR - cor(Xi, XA) = rho for beta[i]==0 & beta[A]!=0
  #           default=AR
  # rho.......correlation coefficient
  #           default=0.5
  # beta......parameter vector
  # seed......seed for random number generation
  #           default=floor(runif(1,1,10000))


  require(MASS)
  set.seed(seed)
  
  # X STRUCTURE
  p <- length(beta)
  A <- which(beta != 0)
  Ac <- (1:p)[-A]
  q <- length(A)
  X.names <- paste0("X",1:p)
  X.mu <- rep(0,p)
  if (corr=="AR") X.sigma <- toeplitz(rho^(0:(p-1)))
  if (corr=="CS") {
    X.sigma <- matrix(rho, p, p)
    diag(X.sigma) <- 1
  }
  if (corr=="IR") {
    X.sigma <- diag(p)
    X.sigma[A,Ac[1]] <- X.sigma[Ac[1],A] <- rho
  }
    
  # compatibility condition
  quad.eig <- t(beta) %*% X.sigma %*% beta
  compat <- drop(q*quad.eig/(sum(abs(beta[A])))^2) # want > 1
  # restricted eigenvalue condition
  resteig <- drop(quad.eig/sum(beta[A]^2)) # want > 1

  X.train <- list()
  y.train <- list()
  snr <- rep(0,iter)
  cond <- rep(0,iter)
  irr <- rep(0,iter)

  for (i in 1:iter) {
    
    # TRAINING DATA
    X.train[[i]] <- mvrnorm(n.train, X.mu, X.sigma)
    colnames(X.train[[i]]) <- X.names
    y.train[[i]] <- X.train[[i]] %*% beta + rnorm(n.train, 0, sigma)

    # signal to noise ratio
    snr[i] <- sqrt(sum((X.train[[i]] %*% beta)^2)) / sigma^2
    # condition number
    eigval <- eigen(cor(X.train[[i]]), only.values=T)$values
    cond[i] <- sqrt(max(eigval)/min(eigval))
    # irrepresentable condition
    if (p-q > 0)
      irr[i] <- max(t(X.train[[i]][,Ac])%*%X.train[[i]][,A] %*%
                      solve(t(X.train[[i]][,A])%*%X.train[[i]][,A]) %*% 
                      sign(beta[A])) else irr[i] <- 0 # want < 1
  }

  # OUTPUT:
  data <- list(beta=beta, rho=rho, sigma=sigma, X.sigma=X.sigma, 
               corr=corr, snr=snr, cond=cond, compat=compat, resteig=resteig,
               irr=irr, n.train=n.train, X.train=X.train, y.train=y.train)
  # TEST DATA
  if (n.test > 0) {
    X.test <- mvrnorm(n.test, X.mu, X.sigma)
    colnames(X.test) <- X.names
    y.test <- X.test %*% beta + rnorm(n.test, 0, sigma)
    data$n.test <- n.test
    data$X.test <- X.test
    data$y.test <- y.test
  }
  
  # VALIDATION DATA
  if (n.valid > 0) {
    X.valid <- mvrnorm(n.valid, X.mu, X.sigma)
    colnames(X.valid) <- X.names
    y.valid <- X.valid %*% beta + rnorm(n.valid, 0, sigma)
    data$n.valid <- n.valid
    data$X.valid <- X.valid
    data$y.valid <- y.valid
  }
  data
}
