

# Bayesian variance component selection in animal model or random intercept and slope (RIS) model.

B. Heuclin, F. Mortier, C. Trottier, S. Tisné and M. Denis

04/26/2021

![](logo.png)

We propose a folded version of the horseshoe prior 
to shrink toward zero standard deviation associated to non-relevant random effect. This work is related to the paper: 
"Continuous shrinkage priors for fixed and random effects selection in linear mixed models: application to genetic mapping", B. Heuclin, F. Mortier, C. Trottier, S. Tisné and M. Denis (submitted).


`RIS_mcmc.R` file contains a R code of the MCMC sampler algorithm for the random intercept and slope (RIS) model using a horseshoe prior for selection of fixed effects and either folded horseshoe, Cauchy or spike-and-slab priors for the identification of relevant random effects through the selection of their associated standard-deviation parameters. The main function is `ris_model`.

`fct.cpp` file contains a cpp function to improve the speed of the RIS MCMC sampler algorithm.

`RIS_exemple.R` file is a R exemple on simulated dataset to run the MCMC sampler of the RIS model (same exemple as below).


`AM_mcmc.R` file contains a R code of the MCMC sampler algorithm for the animal model using either folded horseshoe, Cauchy or spike-and-slab priors for the identification of relevant random effects through the selection of their associated standard-deviation parameters. The main function is `animal_model`.

`AM_exemple.R` file is a R exemple on simulated dataset to run the MCMC sampler of the animal model (same exemple as below).

`simulated_IBD_matrices.Rdata` file is a Rdata object contaning 50 simulated IBD matrices of size 140x140.

## RIS model


### Data simulation 
```{r}
library(extraDistr)

n <- 200                      # number of individuals
q <- 50                       # number of variables
T <- 20                       # number of time points

se2 <- 4                      # residual variance 
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

```


### Fit the varying coefficient model with spike-and-slab prior for variable selection with the main function `VCM_fct` 
```{r}
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

```


### Parameter estimations

#### Fixed effect dans standard deviations of random effects
```{r}
# fixed effects 
plot(apply(chain_hs$beta, 2, median))
abline(0, 0)

# random effects
plot(apply(chain_hs$sdu, 2, median)[-1])
```

#### RIS correlation matrix
```{r}
tmp <- apply(chain_hs$R0, 2, mode)

# R0 estimated at each iteration
R_est <- diag(q+1)
R_est[lower.tri(R_est)] <- tmp
R_est <- as.matrix(Matrix::forceSymmetric(R_est, uplo = "L"))

par(mfrow=c(1, 1))
corrplot::corrplot(R_est_HS, title = "", mar = c(1, 1, 1, 1), 
         tl.cex=0.4, tl.col = 1)
```


#### plot of the evolution over time of regression effects 
```{r}
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
```







---
## Animal model

#### Data simulation 



Load 50 simulated IBD matrices of size 140x140
```{r}
load("simulated_IBD_matrices.Rdata")

n = nrow(A[[1]]); n
q = length(A); q
```

#### Calculate the SVD decomposition of the IBD matrices
```{r}
svdA <- list()
for(k in 1:q){
  svdA[[k]] <- list()
  svd.tmp <- svd(A[[k]])
  svdA[[k]]$u <- svd.tmp$u
  svdA[[k]]$d <- svd.tmp$d
}
```

Simulated responce variables
```{r}
sim <- list()
sim$se2 <- 1
sim$mu <- 3
sim$sdu <- rep(0, q)
sim$sdu[c(1, 5, 10, 15)] <- 2

sim$U <- matrix(NA, n, q)
for(i in 1:q){
  sim$U[, i] <- tcrossprod( svdA[[i]]$u, rmvn(1, rep(0, n), diag(sqrt(svdA[[i]]$d)), isChol = T))
}

sim$Y <- sim$mu + sim$U %*% sim$sdu + c(rmvn(n, 0, sqrt(sim$se2)))
```



#### Fit the animal model with the folded horseshoe prior for the selection of the random effect associated to each IBD matrices using the function `animal_model` in the `AM_mcmc.R Rscript`.

```{r}
library(extraDistr)
source("AM_mcmc.R")

niter <- 15000
burn <- 5000
thin <- 10

chain_hs <- animal_model(y = sim$Y, svdA = svdA, prior="HS", 
                         niter = niter, burnin = burn, thinin = thin, 
                         method = "folded")
```

#### Estimations:
```{r}
sdu_hat_hs = apply(chain_hs$sdu, 2, median)

plot(sdu_hat_hs, ylim=c(0, 3))
points(sim$sdu, col=4)
legend("topright", c("folded horseshoe", "Simulated"), pch=1, col=c(1, 4))


mean(chain_hs$se2); sim$se2
mean(chain_hs$mu); sim$mu
```

#### It is also possible to use folded Cauchy prior or folded spike-and-slab prior for standard deviation selection
```{r}
# folded Cauchy
chain_c <- animal_model(y = sim$Y, svdA = svdA, prior="C", 
                        niter = niter, burnin = burn, thinin = 1, 
                        method = "folded")

sdu_hat_c = apply(chain_c$sdu, 2, median)
plot(sdu_hat_c, ylim=c(0, 3))
points(sim$sdu, col=4)
legend("topright", c("folded Cauchy", "Simulated"), pch=1, col=c(1, 4))


# folded spike-and-slab
chain_ss <- animal_model(y = sim$Y, svdA = svdA, prior="SS", 
                         niter = niter, burnin = burn, thinin = 1,
                         method = "folded")

sdu_hat_ss = apply(chain_ss$sdu, 2, median, na.rm=TRUE)
plot(sdu_hat_ss, ylim=c(0, 3))
points(sim$sdu, col=4)
legend("topright", c("folded S&S", "Simulated"), pch=1, col=c(1, 4))


par(mar = c(2, 2, 1, 1), mfrow = c(1, 1))
plot(colMeans(1*chain_ss$g, na.rm = TRUE), ylim = c(0, 1))
```
#### It is also possible to use truncated distribution instead of folded by setting parameter `method` to "truncated".

```{r}
chain_hs <- animal_model(y = sim$Y, svdA = svdA, prior="HS", niter = niter, burnin = burn, thinin = thin, method = "truncated")

sdu_hat_hs = apply(chain_hs$sdu, 2, median)
plot(sdu_hat, ylim=c(-1, 3))
points(sim$sdu, col=4)
legend("topright", c("folded horseshoe", "Simulated"), pch=1, col=c(1, 4))
```
