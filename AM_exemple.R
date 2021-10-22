# ______________________________________________________________________________________________
# Prior Horseshoe pour selection de composantes de la variance dans un modèle animal
# 13/01/2020
# ______________________________________________________________________________________________

rm(list=ls())
#source("Effet_aleatoire_program.R")
library(mvnfast)
library(invgamma)
library(tmvtnorm)
library(mvtnorm)


# cores = 10
# doParallel::registerDoParallel(cores = cores)


# load("data_aek_loba/data_group_B.Rdata")
# 
# q <- length(marker$Position_cM)
# q
# 
# diff <- marker$Position_cM[2:q] - marker$Position_cM[1:(q-1)]
# 
# treshold <- 45 #(cM)
# i=1
# id.rm <- NULL
# for(i in 1:16){
#   print(i)
#   idx.tmp <- which(marker$Chromosome == i)
#   start <- range(idx.tmp)[1]
#   end <- range(idx.tmp)[2]
#   
#   j=1
#   for(j in (idx.tmp)[-length(idx.tmp)]){
#     if(!j %in% id.rm){
#       diff <- abs(marker$Position_cM[j]-marker$Position_cM[(j+1):end])
#       id.rm <- c(id.rm, which(diff < treshold) + j)
#     }
#   }
# }
# id.rm <- unique(id.rm)
# length(id.rm)
# # marker[id.rm, ]
# 
# q-length(id.rm)
# 
# id.keep <- which(! 1:q %in% id.rm)
# # marker[id.keep, ]
# svdMat <- svdM[id.keep]
# 
# 
# n <- 140
# q <- 50
# 
# A <- list()
# for(k in 1:q){
#   A[[k]] <- svdMat[[k]]$u %*% diag(svdMat[[k]]$d) %*% t(svdMat[[k]]$u)
#   A[[k]] <- A[[k]][1:140, 1:140]
# }
# 
# save(A, file="simulated_IBD_matrices.Rdata")

load("simulated_IBD_matrices.Rdata")

n <- nrow(A[[1]])
q <- length(A); q

svdA <- list()
for(k in 1:q){
  svdA[[k]] <- list()
  svd.tmp <- svd(A[[k]])
  svdA[[k]]$u <- svd.tmp$u
  svdA[[k]]$d <- svd.tmp$d
}

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


# MCMC --------------------------------------------------------------------

source("AM_mcmc.R")

niter <- 5000
burn <- 1000
thin <- 10

chain_hs <- animal_model(y = sim$Y, svdA = svdA, prior="HS", niter = niter, burnin = burn, thinin = thin, method = "folded")

sdu_hat_hs = apply(chain_hs$sdu, 2, median)
plot(sdu_hat_hs, ylim=c(-1, 3))
points(sim$sdu, col=4)
legend("topright", c("folded horseshoe", "Simulated"), pch=1, col=c(1, 4))


mean(chain_hs$se2); sim$se2
mean(chain_hs$mu); sim$mu



# Cauchy and SS -----------------------------------------------------------

chain_c <- animal_model(y = sim$Y, svdA = svdA, prior="C", niter = niter, burnin = burn, thinin = 1, method = "folded")

sdu_hat_c = apply(chain_c$sdu, 2, median)
plot(sdu_hat_c, ylim=c(-1, 4))
points(sim$sdu, col=4)
legend("topright", c("folded Cauchy", "Simulated"), pch=1, col=c(1, 4))


chain_ss <- animal_model(y = sim$Y, svdA = svdA, prior="SS", niter = niter, burnin = burn, thinin = 1, method = "folded")

sdu_hat_ss = apply(chain_ss$sdu, 2, median, na.rm=TRUE)
plot(sdu_hat_ss, ylim=c(0, 3))
points(sim$sdu, col=4)
legend("topright", c("folded S&S", "Simulated"), pch=1, col=c(1, 4))


par(mar = c(2, 2, 1, 1), mfrow = c(1, 1))
plot(colMeans(1*chain_ss$g, na.rm = TRUE), ylim = c(0, 1))


# Truncated ---------------------------------------------------------------

chain_hs <- animal_model(y = sim$Y, svdA = svdA, prior="HS", niter = niter, burnin = burn, thinin = thin, method = "truncated")

sdu_hat_hs = apply(chain_hs$sdu, 2, median)
plot(sdu_hat, ylim=c(-1, 3))
points(sim$sdu, col=4)
legend("topright", c("folded horseshoe", "Simulated"), pch=1, col=c(1, 4))




