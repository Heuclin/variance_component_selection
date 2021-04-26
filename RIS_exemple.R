library(extraDistr)

n <- 200                      # number of individuals
q <- 20                      # number of variables
T <- 10                       # number of time points

se2 <- 2                     # residual variance 
rho <- 0.7                    # autoregressive parameter on the residual correlation       structure 
nb_noeuds <- round(T/3)       # knot number for B-spline basis

#### Simulation
sim <- list()

# intercept over time
sim$mu <-  1 + sin(pi*1:T/20) 

# plot the intercept onver time
plot(sim$mu, t='l'); abline(0, 0)


# among the 50 variables, we fixe the first four variables with non zero effect over time
sim$beta <- matrix(0, T, q)
sim$beta[, 1] <- 4 + cumsum(rnorm(T, 0, 2))
sim$beta[, 2] <- 2*cos(pi * (1:T-25)/15)+1:T/50
sim$beta[, 3] <- c(rep(2, round(T/3)), rep(0, round(T/3)), rep(1, T-2*round(T/3)))
sim$beta[, 4] <- 2*60 / (25 + (1:T - T/2)^2)

# plot variable effects over time
par(mfrow = c(4,1), mar=c(2, 4,0.5, 1))
for(j in 1:4) {plot(sim$beta[, j], ylab = paste0("beta_", j)); lines(sim$beta[, j]); abline(0, 0, lty = 2)}

sim$rho <-  rho
sim$se2 <- se2

# SNPs simulation
X <- scale(matrix(sample(0:2, n*q, replace = T), n, q))

# construction of the autoregressive residual correlation structure
Gamma <- matrix(sim$rho, T, T); for(i in 1:T) for(j in 1:T) Gamma[i, j] <- sim$rho^abs(i-j) 

# simulation of the response variable over time
Y <- matrix(NA, n, T)
for(i in 1:n){
  Y[i, ] <- sim$mu + sim$beta %*% as.numeric(X[i, ]) + t(mvtnorm::rmvnorm(1, rep(0, T), sim$se2*Gamma))   
}

# plot response variable over time 
par(mfrow= c(1, 1), mar = c(4, 4, 4, 1))
matplot(t(Y), t='l', main = "matplot of Y", xlab = "time", ylab = "")



y <- c(Y)
y <- y - mean(y)
x <- pracma::repmat(X, T, 1)
ind <- as.factor(kronecker(1:T, rep(1, n)))


source("ris_mcmc.R")
Rcpp::sourceCpp('fct.cpp')


niter <- 5000
burnin <- 1000
thinin <- 10


chain_hs <- ris_model(
  y=y, X=x, niv=ind, prior_fixed = "HS", prior_rand = "HS", 
  niter=niter, burnin = burnin, thinin = thinin,
  correlation_rand = TRUE,  # unknown RIS correlation matrix
  individual_effect = TRUE, # add an individual random effect
  AR_residual = TRUE, # supposed an AR(1) structur on residuals of one individual
  method = "folded", 
  K=30, # shrinkage parameter on angle parameters (for the RIS correlation matrix estimation)
  adaptive = FALSE, # adaptive parameter K (beta version)
  id_bloc_R =  q+1,  # Size of blocs of the RIS correlation matrix R (q+1 for one full matrix)
  epsilon_K=1,  # parameter for the MH step of K
  epsilon_angle=0.1, # parameter for the MH step of angles
  epsilon_rho=0.01 # parameter for the MH step of rho (AR(1) residual structur)
)

# fixed effects 
plot(apply(chain_hs$beta, 2, median))
abline(0, 0)

# random effects
plot(apply(chain_hs$sdu, 2, median)[-1])



# RIS correlation matrix estimation
tmp <- apply(chain_hs$R0, 2, mode)

# R0 estimated at each iteration
R_est <- diag(q+1)
R_est[lower.tri(R_est)] <- tmp
R_est <- as.matrix(Matrix::forceSymmetric(R_est, uplo = "L"))

library(corrplot)

par(mfrow=c(1, 1))
corrplot(R_est_HS, title = "", mar = c(1, 1, 1, 1), 
         tl.cex=0.4, tl.col = 1)







# effects
beta_hs <- apply(chain_hs$beta, 2, mean)
blup_hs <- matrix(NA, q+1, T)
tmp <- rep(1:(q+1), T)
i <- 1
for(i in 1:(q+1)){
    tmp2 <- chain_hs$u[, tmp==i] * chain_hs$sdu[, i]
    blup_hs[i, ] <- apply(tmp2, 2, mean)
}

par(mfrow=c(3, 7), mar=c(2, 2, 2, 1))
# Intercept
plot(colMeans(Y), t='l', main="Intercept", col = 2, lty=1, lwd = 1); abline(0, 0, lty = 3)
lines(mean(Y)+blup_hs[1, ], main="Intercept", t='l', lwd=2, lty=1)
# Effects
for(i in 1:q){
  plot(blup_hs[i+1, ] + beta_hs[i], main = i, t='l', lty=1, lwd=2, ylim = c(-5, 5)); 
  abline(0, 0, lty=3)
  lines(sim$beta[, i], col=2)
}







