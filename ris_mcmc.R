

require(truncnorm)
require(mvnfast)
library(Matrix)





# Heuclin Benjamin
#
#   MAJ du 10/11/2020
#     correction du folded (abs)
# 
#   MAJ du 29/10/2020
#     modification pour estimer un R0 par block
#     ajout du paramètre R0_restriction
# 
#  MAJ 28/10/2020
#   ajout du paramètre adaptive = TRUE pour la mise à jour de R0 avec une etape adatatatif Metropolis  
# 
#  MAJ 27/10/2020
#   correction d'un bug, oublie de décommenter la mise à jour de gamma_0_inv
#   modification de la fonction rtnorm  (donnée par Fred)
#
#  MAJ 26/10/202
#     possibilité d'échantilloner s suivant une loi normale tronquée
#     joint or marginal truncated MH step
# 
#  MAJ le 23/10/2020
#   optimisation de la l'echantillonage de R0, suppression du double calcul de l'inverse de R0
#   optimisation du kronecker product entre Gamma_0_inv et diag(), passage des deux matrices en sparse
#
#  MAJ le 22/10/2020




dtnorm <- function(x, mu=0., sigma=1., A=-Inf, B=Inf){
  # Author: Dimitris Rizopoulos
  stopifnot(A<B)
  dnorm(x, mu, sigma)/(pnorm(B, mu, sigma)-pnorm(A, mu, sigma))
}

rtnorm <- function(n, mu=0., sigma=1., A=-Inf, B=Inf){
  # Author: Dimitris Rizopoulos
  stopifnot(A<B)
  mu+sigma*qnorm(runif(n, pnorm(A, mu, sigma), pnorm(B, mu, sigma)))
}



#' Calculate the empirical mode of a distribution from a sample vector
#'
#' @param x vector, samples from a distribution 
#'
#' @return Return the empirical mode
#' @export
#'
#' @examples 
#' a <- rnorm(1000)
#' mode(a)
mode <- function(x){
  return(density(x)$x[which.max(density(x)$y)])
}




#' Caculate the logarithm of the sum of two values from their logarithms
#'
#' @param log_a numeric value
#' @param log_b numeric value
#'
#' @return Return log(exp(log_a) + expt(log_b))
#' @export
#'
#' @examples
#' a <- 2
#' b <- 1
#' log(a+b)
#' 
#' log_a <- log(a)
#' log_b <- log(b)
#' sumlogs(log_a, log_b)
sumlogs <- function(log_a, log_b){
  M <- max(log_a, log_b)
  return(M+log(sum(exp(c(log_a, log_b)-M))))
}








#' Calculate the Cholesky factor of a correlation matrix from hyperspherical parameterization 
#'
#' @param theta_mat lower triangular matrix of the angles associated to the hyperspherical parameterization of the Cholesky factor.
#' 
#' @return Return the Cholesky lower trianguar matrix of a correlation matrix from hyperspherical parameterization 
#' @seealso \code{\link[randcorr]{randcorr}, \link{sample_thetaPourhamandi}}
#' @keywords 
#' @references \url{https://www.sciencedirect.com/science/article/abs/pii/S0167715215002011}
#' @export
#'
#' @examples 
#' theta <- sample_thetaPourhamandi(10)
#' corrPourhamandi(theta)
corrPourhamandi <- function (theta_mat) 
{
  q <- ncol(theta_mat)
  # q  <- (1+sqrt(1+4*2*length(theta)))/2
  L = matrix(0, q, q)
  # theta_mat <- matrix(0, q, q)
  # theta_mat[lower.tri(theta_mat)] <- theta
  L[1, 1] <- 1
  for (i in 2:q){
    # L[i, 2:i] = cumprod(sin(theta_mat[i, 1:i - 1]))
    L[i, 1:i] = c(1, cumprod(sin(theta_mat[i, 1:(i-1)])))
  }
  cos_mat = cos(theta_mat)
  cos_mat[upper.tri(cos_mat)] = 0
  B = L * cos_mat # elements products:
  return(B)
  
  # R = tcrossprod(B)
  # return(list(R=R, B=B))
}






#' Uptade the Cholesky factor of a correlation matrix from hyperspherical parameterization after changing one angle 
#' @title  Uptade the Cholesky factor
#'
#' @param B oldest Cholesky lower trianguar matrix of a correlation matrix from hyperspherical parameterization 
#' @param theta_mat new lower triangular matrix of the angles associated to the hyperspherical parameterization of the Cholesky factor.
#' @param i row indice of the new value in theta_mat
#' @param j column indice of the new value in theta_mat
#'
#' @return Return the new Cholesky lower trianguar matrix of a correlation matrix from hyperspherical parameterization 
#' @seealso \link{sample_thetaPourhamandi}}, \code{\link[randcorr]{randcorr}
#' @export
#'
#' @examples
#' theta <- sample_thetaPourhamandi(10)
#' B <- corrPourhamandi(theta)
#' theta[3, 1] <- pi/2
#' B_new <- corrPourhamandi_new(B, theta, 3, 1)
#' 
corrPourhamandi_new <- function (B, theta_mat, i, j) 
{
  B[i, j:i] <- c(1, cumprod(sin(theta_mat[i, 1:(i-1)])))[j:i]*
    c(cos(theta_mat[i, j:(i-1)]), 1)
  return(B)
}







#' Generate a q x q random angles matrix associated to hyperspherical parameterization of the Cholesky factor of a correlation matrix
#' @title Generate a q x q random angles matrix
#'
#' \code{\link[randcorr]{randcorr.sample.sink}}
#'
#' @param q dimension of correlation matrix
#'
#' @return A random q x q angles matrix
#' @details Angle[i, j] is sampled from a distribution with probability density function sin^k(theta) (0 < theta < pi) with k= q-j using the efficient sampling algorithm \code{\link[randcorr]{randcorr.sample.sink}}
#' 
#' @references 
#' 
#'  [1] Makalic, E. & Schmidt, D. F.
#'  An efficient algorithm for sampling from sin^k(x) for generating random correlation matrices
#'  arXiv:1809.05212, 2018 \url{https://arxiv.org/abs/1809.05212}
#' 
#'  [2] Mohsen Pourahmadi and Xiao Wang,
#'  Distribution of random correlation matrices: Hyperspherical parameterization of the Cholesky factor,
#'  Statistics & Probability Letters, Volume 106, Pages 5-12, 2015.
#'  
#'  [3] Enes Makalic and Daniel F. Schmidt
#'  An efficient algorithm for sampling from sin^k(x) for generating random correlation matrices,
#'  arXiv:1809.05212, 2018.
#' @seealso \link{sample_thetaPourhamandi}}, \code{\link[randcorr]{randcorr.sample.sink}}
#' @export
#'
#' @examples
sample_thetaPourhamandi <- function (q) 
{
  theta = matrix(0, q, q)
  # browser()
  for (j in 1:(q - 1)) theta[(j + 1):q, j] = randcorr::randcorr.sample.sink((q - j))
  # for (j in 1:(q - 1)) theta[(j + 1):q, j] = randcorr::randcorr.sample.sink(1)
  return(theta)
}





#' Empirical estimation of the jointe posterior probabilities of inclusion of variables using the spike and slab prior.
#'
#' @param gamma matrix of posterior samples of indicator inclusion associated to the spike and slab prior. Samples are in row and variables are in column.
#' @param names names of the variables
#' @param n numver of the highest probabilities which should be return. The default value is 3.
#'
#' @return Return the n higher jointe posterior probabilities of inclusion
#' @export
#'
#' @examples
#' names <- paste0("x.", 1:15)
#' gamma <- matrix(0, 1000, 15)
#' set.seed(1)
#' my_prob <- c(1, 0.9, 0.8, rep(0, 6), 0.5, 0.3, rep(0, 4))
#' for( i in 1:15) gamma[, i] <- sample(0:1, 1000, prob=c(1 - my_prob[i], my_prob[i]), replace = TRUE)
#' proba_jointe(gamma, names, 5)
#' 
proba_jointe <- function(gamma, names, n = 3){
  require(stringr)
  idx <- 1:ncol(gamma)
  proba.g.jointe <- sort(table(apply(1*gamma, 1, function(l) paste(l, collapse=""))) / nrow(gamma), decreasing = T ) [1:n]
  for(i in 1:n){
    st <- names(proba.g.jointe)[i]
    names(proba.g.jointe)[i] <- paste(names[str_locate_all(st, "1")[[1]][, 1]], collapse = ", ")
  }
  return(proba.g.jointe)
}




#__________________________________________________________________________________________________________
#__________________________________________________________________________________________________________

lcst_Norm <- function(K, q){
  commun <- (2*K+(1:(q-1)))*0.5
  lnum <- lgamma(commun+1)
  lden <- lgamma(commun+0.5) #lgamma(commun+0.5)
  delta <- sum((lnum-lden-0.5*log(pi))*(1:(q-1)))
  return(delta)
}




# l_p_R0_K_fct <- function(q, theta_vec, K){
#   # browser()
#   res <- 2*K*sum(log(sin(theta_vec))) + lcst_Norm(K, q)
#   # for(j in 1:(q-1)){
#   #   # print(j * (log(gamma( (2*K+1)/2 + 1 )) -0.5*log(pi) - log(gamma((2*K+j+1)/2)) ) )
#   #   res <- res + j * (lgamma( (2*K+j)/2 + 1 ) -0.5*log(pi) - lgamma((2*K+j+1)/2))
#   # }
#   return(res)
# }





#__________________________________________________________________________________________________________
#__________________________________________________________________________________________________________



#' MCMC algorithm 
#' 
#' New version (22/10/2020) correction de bugs dans l'échantillonage de R0 
#' + optimisation de l'échantillonage des effets aleatoires u
#'
#' @param y n-vector of the centered response variable
#' @param x n x q covariable matrix 
#' @param ind n-vector indicating the level of each observation
#' @param prior_fixed character indicating which prior must be used to achieve selection of fixed effects. "SS" denote the Spike and Slab prior, "HS" denote the Horseshoe prior. The default value is "SS".
#' @param prior_rand character indicating which prior must be used to achieve selection of variance componantes. "SS" denote the Spike and Slab prior, "HS" denote the Horseshoe prior, "cauchy" denote the cauchy prior, "DG" denote the double-gamma prior, "L" denote the Laplace prior and "NIG" denote the Normal Inverse-Gamma prior. The default value is "SS".
#' @param correlation logical value indicating whether to learn the correlation matrix of the random effects. The default value is TRUE.
#' @param individual_effect logical value indicating whether to add a random individual effect. The default value is FALSE.
#' @param niter positive integer, indicating the number of MCMC iterations to perform, including the burn-in. Has to be larger than  nburn + 1. The default value is 10000.
#' @param burnin non-negative integer, indicating the number of iterations discarded as burn-in. Has to be smaller than to niter - 1. The default value is round(niter / 2).
#' @param thinin positive integer, indicating the degree of thinning to be performed. Every nthin draw is kept and returned. The default value is 1, implying that every draw is kept.
#' @param epsilon positive, real number. Determines the standard deviation of the proposal distribution for the Metropolis Hastings step for theta. The default value is 0.3.
#' @param epsilon_a_s positive, real number. Only available for prior_rand = "DG". Determines the standard deviation of the proposal distribution for the Metropolis Hastings step for a_xi. The default value is 1.
#'
#' @details For details concerning the algorithm please refer to the paper by ...
#' @return The value returned is a list object containing samples for parameters:
#' 
#' \itemize{
#' \item  alpha, only if individual_effect==TRUE
#' \item  sd_alpha, only if individual_effect==TRUE
#' \item  beta
#' \item  g_b, only if prior_fixed=="SS"
#' \item  pi_beta, only if prior_fixed=="SS"
#' \item  sdu
#' \item  sdu_tild
#' \item  g_s, only if prior_fixed=="SS"
#' \item  pi_sdu, only if prior_fixed=="SS"
#' \item  se2
#' \item  a_s, only if prior_rand=="DG"
#' \item  kappa2, only if prior_rand=="DG"
#' \item  omega2_beta
#' \item  omega2_sdu
#' \item  theta
#' \item  u
#' \item  tau2
#' }
#'  The list contains also:
#'  \itemize{
#' \item  acc_rate_theta: acceptence rate of the MH step for sample theta, only if correlation == TRUE.
#' \item  acc_rate_a_s:  acceptence rate of the MH step for sample a_s, only if prior_rand=="DG".
#' }
#' 
#' @note     To cite this package please reference: 
#' 
#' @export
#'
#' @examples
ris_model <- function(y, X, niv, prior_fixed = "HS", prior_rand="HS", 
                     correlation_rand = TRUE, 
                     individual_effect = TRUE,
                     AR_residual = TRUE,
                     niter, burnin, thinin, 
                     method = "folded",
                     K=10, 
                     adaptive = FALSE, 
                     id_bloc_R = ncol(X)+1,
                     epsilon_K=1,
                     epsilon_angle=0.1, 
                     epsilon_rho=0.01
                     )
{
  if(!prior_fixed  %in% c("SS", "HS")) stop("prior_fixed choice must be 'SS'or 'HS'." ) 
  if(!prior_rand  %in% c("SS", "HS", "DG", "cauchy", "NIG", "L")) stop("prior_rand choice must be 'SS', 'HS', 'DG', 'cauchy', 'L' or 'NIG'." )
  if(burnin >= niter ) stop("'burnin' must be lower than 'niter'")
  
  print("Variance components selection in RIS models")
  # if(prior_rand=="HS") print("Horseshoe prior for variance component in RIS models")
  # if(prior_rand=="cauchy") print("Cauchy prior for variance component in RIS models")
  # if(prior_rand=="NIG") print("Gaussian prior for variance component in RIS models")
  
  # Settings _________________________________________________________________
  levels(niv) <- 1:length(unique(levels(niv)))
  n_niv       <- length(levels(niv))
  no          <- length(y)
  ni          <- table(niv)
  m <- ni[1]
  q_fix       <- ncol(X)
  q_rand      <- q_fix + 1
  # K=0
  index <- unlist(lapply(1:(q_fix), function(l) 2*K+rep(q_rand - l, q_rand - l)))
  
  
  
  #sort by increasing factor
  y       <- y[order(niv)]
  X       <- X[order(niv), ]
  niv <- sort(niv)
  J       <- t(as(niv,  Class = "sparseMatrix"))
  L       <- matrix(0, (q_rand)^2, q_rand)
  pos     <- cbind(which(diag(q_rand)>0), 1:(q_rand))
  L[pos]  <- 1
  
  if(individual_effect){
    table_ind <- table(niv)
    zz <- matrix(0, length(y), max(table_ind))
    for(k in 1:n_niv) zz[niv==k, 1:table_ind[k]] <- diag(table_ind[k])
  } 
  
  # Outputs _________________________________________________________________
  chain             <- list()
  if(individual_effect) chain$alpha       <- matrix(0, floor(niter - burnin)/thinin, ncol(zz))
  if(individual_effect) chain$sd_alpha    <- rep(NA, floor(niter - burnin)/thinin)
  # chain$mu          <- rep(NA, floor(niter - burnin)/thinin)
  chain$beta        <- matrix(0, floor(niter - burnin)/thinin, q_fix)
  if(prior_fixed=="SS") chain$g_b         <- matrix(0, floor(niter - burnin)/thinin, q_fix)
  if(prior_fixed=="SS") chain$pi_beta     <- rep(NA, floor(niter - burnin)/thinin)
  chain$sdu         <- matrix(0, floor(niter - burnin)/thinin, q_rand)
  chain$sdu_tild    <- matrix(0, floor(niter - burnin)/thinin, q_rand)
  chain$s    <- matrix(0, floor(niter - burnin)/thinin, q_rand)
  chain$g_s         <- matrix(0, floor(niter - burnin)/thinin, q_rand)
  if(prior_rand=="SS") chain$pi_sdu         <- rep(NA, floor(niter - burnin)/thinin)
  chain$tau2 <- matrix(NA, floor(niter - burnin)/thinin, 2) # by convention tau[, 1]=tau_beta et tau[, 2]=tau_sdu
  chain$omega2_beta <- matrix(1, floor(niter - burnin)/thinin, q_fix) 
  chain$omega2_sdu  <- matrix(1, floor(niter - burnin)/thinin, q_rand)
  chain$se2         <- rep(NA, floor(niter - burnin)/thinin)
  chain$sd_t        <-  matrix(NA, floor(niter - burnin)/thinin, n_niv)
  chain$rho         <- rep(NA, floor(niter - burnin)/thinin)
  chain$u           <- matrix(NA, floor(niter - burnin)/thinin, n_niv*q_rand)
  if(prior_rand=="DG") chain$a_s        <- rep(NA, floor(niter - burnin)/thinin)
  if(prior_rand=="DG") chain$kappa2     <- rep(NA, floor(niter - burnin)/thinin)
  # chain$k         <- rep(NA, floor(niter - burnin)/thinin)
  
  # Initialisation __________________________________________________________
  pi_beta <- pi_sdu <- 0.5
  mu <- 0 # mean(y)
  if(individual_effect) alpha <- rep(0, ncol(zz))
  if(individual_effect) sd_alpha <- 10000
  U                 <-  matrix(rnorm(q_rand*n_niv), q_rand, n_niv)
  u                 <- c(U)
  # u <- init$u
  # U               <- matrix(u, q_rand, n_niv)
  beta <- chain$beta[1, ] <- rnorm(q_fix)
  sdu_tild <- s <- rep(1, q_rand) #rnorm(q_rand)
  if(method == "folded" | method == "truncated") sdu_tild <- abs(sdu_tild)
  sdu <- abs(sdu_tild)
  se2 <- chain$se2[1] <- 1 # rgamma(1, 2, 2)
  
  
  # tmp <- 0*diag(q_rand)
  # tmp[lower.tri(tmp)] <- 
  # tmp 
  # lower.tri(tmp)
  # 
  # df_id_theta <- data.frame(id = 1: (q_rand*(q_rand-1)/2), id_row = rep(1:q_rand, q_rand), id_col =  rep(1:q_rand, each = q_rand))
  # df_id_theta <- df_id_theta[lower.tri(tmp), ]
  # df_id_theta
  # 
  # 
  
  #---- init RO----
  if(is.null(id_bloc_R)) id_bloc_R <- q_rand # 1:(q_rand*(q_rand-1)/2)
  chain$K         <- matrix(NA, floor(niter - burnin)/thinin, length(id_bloc_R))
  chain$R0          <- matrix(0, floor(niter - burnin)/thinin, 0.5*q_rand*(q_rand-1))
  K <- rep(K, length(id_bloc_R)) # sample(1:100, length(id_bloc_R))
  
  if(correlation_rand){
    # id_tot <- 1:(q_rand*(q_rand-1)/2)
    # id_R0_1 <- R0_restriction
    # id_R0_0 <- which(! id_tot %in% id_R0_1 )
    
    id_bloc_vec <- rep(1:length(id_bloc_R), times = id_bloc_R)
    
    epsilon_K <- rep(3, length(id_bloc_R))
    epsilon_angle <- rep(epsilon_angle, length(id_bloc_R))
    df_id_theta <- iR0_list <- id_lower.tri <- B_list <- R0_list <- theta <- theta_mat <- chain$theta <- theta_hist <- theta_hist_epsilon <- K_hist <- K_hist_epsilon <- list()
    # chain$K         <- matrix(NA, floor(niter - burnin)/thinin, length(id_bloc_R))
    # chain$R0          <- matrix(0, floor(niter - burnin)/thinin, 0.5*q_rand*(q_rand-1))
    
    i=2
    for(i in 1:length(id_bloc_R)){
      chain$theta[[i]]       <- matrix(0, floor(niter - burnin)/thinin, 0.5*id_bloc_R[i]*(id_bloc_R[i]-1))
      theta_hist[[i]]        <- matrix(0, niter, 0.5*id_bloc_R[i]*(id_bloc_R[i]-1))
      theta_hist_epsilon[[i]] <- matrix(epsilon_angle[i], niter, 0.5*id_bloc_R[i]*(id_bloc_R[i]-1))
      K_hist[[i]]        <- matrix(0, niter, 1)
      K_hist_epsilon[[i]] <- matrix(epsilon_K[i], niter, 1)
      
      theta_mat[[i]] <- matrix(0, id_bloc_R[i], id_bloc_R[i]); 
      theta_mat[[i]][lower.tri(theta_mat[[i]])] <- truncnorm::rtruncnorm(length(theta_mat[[i]][lower.tri(theta_mat[[i]])]), mean=pi/2, sd=0.5, a=0, b=pi) 
      chain$theta[[i]][1, ] <- theta_hist[[i]][1, ] <- theta[[i]] <- theta_mat[[i]][lower.tri(theta_mat[[i]])]
      if(id_bloc_R[i] == 1){
        B_list[[i]] <- 1
        R0_list[[i]] <- 1
        iR0_list[[i]] <- 1
        id_lower.tri[[i]] <- NULL
        df_id_theta[[i]] <- NULL
      }else{
        B_list[[i]] <- corrPourhamandi(theta_mat[[i]])
        R0_list[[i]] <- tcrossprod(corrPourhamandi(theta_mat[[i]]))
        iR0_list[[i]] <- solve(R0_list[[i]])
        id_lower.tri[[i]] <- lower.tri(theta_mat[[i]])
        df_id_theta[[i]] <- data.frame(id = 1: (id_bloc_R[i]*(id_bloc_R[i]-1)/2), id_row = rep(1:id_bloc_R[i], id_bloc_R[i])[id_lower.tri[[i]]], id_col =  rep(1:id_bloc_R[i], each = id_bloc_R[i])[id_lower.tri[[i]]])
      }
      # indices for angle of pourhamadi correlation matrix
    }
    
    B <- as.matrix(Matrix::bdiag(B_list))
    R0 <- Matrix::bdiag(R0_list)
    iR0 <- Matrix::bdiag(iR0_list)
  }else{
    theta_mat <- matrix(0, q_rand, q_rand); theta_mat[lower.tri(theta_mat)] <- pi/2
    theta <- theta_mat[lower.tri(theta_mat)]
    B <- R0 <- iR0 <- diag(q_rand)
    id_lower.tri <- lower.tri(theta_mat)
    
  }
  
  
  tau2 <- chain$tau2[1, ] <- c(1, 1)
  if(prior_fixed=="SS"){  
    g_b <- chain$g_b[1, ] <- sample(c(0, 1), q_fix, replace = TRUE)
    beta[g_b==0] <- 0
    chain$beta[1, ] <- beta
  }else{ # Horseshoe
    g_b <- rep(1, q_fix)
    tau2 <- chain$tau2[1, ] <- c(1e-5, 1)
  }
  xi2 <- rep(1, 2)
  omega2_beta <- chain$omega2_beta[1, ] <- rep(1, q_fix)
  nu2_beta <- rep(1, q_fix)
  omega2_sdu <- chain$omega2_sdu[1, ] <- rep(1, q_rand)
  nu2_sdu <- rep(1, q_rand)
  if(prior_rand == "cauchy"){
    nu2_sdu <- rep(2, q_rand)
  }
  if(prior_rand=="SS"){
    g_s <- chain$g_s[1, ] <- sample(c(0, 1), q_rand, replace = TRUE)
    sdu[g_s == 0] <- 0
  }else{
    g_s <- chain$g_s[1, ] <- rep(1, q_rand)
  }
  if(prior_rand == "DG"){
    epsilon_a_s=1
    d1 <- 0.001
    d2 <- 0.001
    kappa_2 <- 20 #rgamma(1, shape = d1, rate = d1)
    
    b_w <- 10
    a_s <- 0.1 # rexp(1, b_w)
    l_p_omega_sdu_i <- sum(dgamma(omega2_sdu, a_s, a_s*kappa_2/2, log = TRUE))
    l_p_a_s_i <- dexp(a_s, b_w, log=TRUE)
  }
  
  #___________________________________________________________
  
  diag_niv <- Diagonal(n=n_niv)
  diag_q_rand <- Diagonal(q_rand)
  diag_q_rand_2 <- diag(q_rand)
  tZ <- KhatriRao(t(J),  t(cbind(1, X)))
  Z <- t(tZ)
  zu <- crossprod(tZ, crossprod(kron_cpp(U, diag_q_rand_2), L))
  zus <- zu %*% sdu_tild
  xzu <- cbind(X, zu)
  iR0 <- chol2inv(t(B))
  counter_tot <- counter_1 <- counter_a_s <- counter_rho <- counter_1_K <- counter_tot_K <- 0
  z_list <- list(); for(k in 1:n_niv) z_list[[k]] <- cbind(1, X[niv==k, ])
  if(individual_effect){
    zz_alpha <- zz %*% alpha
    zztzz <- crossprod(zz)
    svd_zztzz <- svd(zztzz)
  }else{ 
    zz_alpha <- rep(0, length(y))
  }
  par <- c(beta, sdu_tild)
  
  ii <- 1
  
  
  if(AR_residual) rho <- runif(1) else rho <- 0 
  Gamma_0 <- matrix(rho, n_niv, n_niv); for(i in 1:n_niv) for(j in 1:n_niv) Gamma_0[i, j] <- rho^abs(i-j) # matrice intervenant dans la variance résiduelle
  diag_ind <- Diagonal(n=ni[1])
  
  Gamma_0_inv <- update_Gamma_inv(diag_niv, rho)
  Gamma.inv <- kron_cpp_s_s(Gamma_0_inv, diag_ind)
  sd_t <- rep(1, n_niv)
  x <- cbind(1, X[niv == 1, ])
  
  # MCMC ___________________________________________________________
  
  print("0 %")
  for(it in 2:niter){
    # print(it)
    # if(it %% (niter/100) == 0) print(paste((it %/% (niter/100))*1, "%"))
    if(it %% (niter/10) == 0) print(paste((it %/% (niter/10))*10, "%"))
    # if(it == 100) browser()
    
    # --------------------- update mu
    # mu <- rnorm(1, mean(y-zz_alpha - xzu %*% par), sd = sqrt(se2/no))
    
    # --------------------- alpha update  (individual effect) ----
    if(individual_effect){
      
      # if(it == 20) browser()
      
      one_Gamma_0_inv_one <- t(rep(1, n_niv)) %*% Gamma_0_inv %*% rep(1, n_niv)
      is_alpha <- 1/ c(one_Gamma_0_inv_one[1, 1] / se2 + 1 / sd_alpha^2)
      iS_alpha <- is_alpha * diag(m)
      mu_alpha <- is_alpha * crossprod(zz, Gamma.inv %*% (y - X %*% beta - zu %*% sdu_tild))/se2
      # mu_alpha <- is_alpha / se2 * t(zz) %*% Gamma.inv %*% (y - X %*% beta - zu %*% sdu_tild)
      # mu_alpha <- is_alpha / se2 * kronecker(matrix(1, 1, n_niv) %*% Gamma_0_inv, diag(m)) %*% (y - X %*% beta - zu %*% sdu_tild)
      # mu_alpha <- is_alpha / se2 * kronecker( colSums(Gamma_0_inv), diag(m)) %*% (y - X %*% beta - zu %*% sdu_tild)
      alpha <- c(mvnfast::rmvn(1, mu_alpha, iS_alpha))
      
      # iS_alpha2 <- diag(1/diag(crossprod(zz, Gamma.inv%*% zz)/se2 + diag(ncol(zz))/sd_alpha^2))
      # mu_alpha2 <- crossprod(iS_alpha2, crossprod(zz, Gamma.inv %*% (y - X %*% beta - zu %*% sdu_tild))/se2)
      # # alpha <- c(mvnfast::rmvn(1, mu_alpha2, iS_alpha2))
      
      zz_alpha <- zz %*% alpha
      sd_alpha <-  sqrt(1/rgamma(1, 1/2+ncol(zz)/2, rate=1/2 + crossprod(alpha)/2))
    }
    
    # browser()
    # --------------------- beta update -----
    if(prior_fixed=="SS"){
      for(k in 1:q_fix){
        y_tild <- y - zz_alpha - X[, -k] %*% beta[-k] - zus
        Sigma_bk <- 1/(crossprod(X[, k], Gamma.inv %*% X[, k])/se2 + 1/(se2 * tau2[1]))
        # Sigma_bk <- 1/(crossprod(x[, k]) * one_Gamma_0_inv_one /se2 + 1/(se2 * tau2[1]))
        # log_R <- as.numeric(log(pi_beta)- log(1-pi_beta) - 0.5 *log(crossprod(X[, k], Gamma.inv %*% X[, k]) + tau2[1]) -0.5 * log(tau2[1]) + (0.5 * t(y_tild) %*% Gamma.inv %*% X[, k] %*% Sigma_bk %*% t(X[, k]) %*% Gamma.inv %*% y_tild / (se2)^2 ))
        log_R <- as.numeric(log(pi_beta)- log(1-pi_beta) - 0.5 *log(crossprod(X[, k], Gamma.inv %*% X[, k]) + tau2[1]) -0.5 * log(tau2[1]) + (0.5 * t(y_tild) %*% Gamma.inv %*% X[, k] %*% Sigma_bk %*% t(X[, k]) %*% Gamma.inv %*% y_tild / (se2)^2 ))
        g_b[k] <- 1*(log(runif(1)) < log_R - sumlogs(log_R, 0))
        
        if(g_b[k]==0) {
          beta[k] <- 0
        }else{
          beta[k] <- rnorm(1, as.numeric(Sigma_bk * crossprod(X[, k], Gamma.inv %*% y_tild))/se2, sqrt(as.numeric(Sigma_bk)))
        }
      }
      pi_beta <- rbeta(1, 1+sum(g_b), 1+q_fix-sum(g_b));
    }else{ # horseshoe
      if(prior_rand == "L"){ tmp_V <- c(omega2_beta) }else{ tmp_V <- c(tau2[1]*omega2_beta) }
      # iV <- Diagonal(x = 1/tmp_V)
      # iS_part_2 <- solve(crossprod(X, Gamma.inv %*% X) + iV)
      # mu_part_2 <- crossprod(iS_part_2, crossprod(X, Gamma.inv %*% (y - zz_alpha - zus)))
      # beta <- as.vector(mvnfast::rmvn(1, mu_part_2, se2*iS_part_2))
      
      SVD <- svd(crossprod(X, Gamma.inv %*% X))
      d <- sqrt(se2) * diag(1/sqrt(SVD$d + 1/tmp_V))
      mu_part <- crossprod(SVD$u %*% d, crossprod(X, Gamma.inv %*% (y - zz_alpha - zus))) /se2
      beta <- as.vector(SVD$u %*% d %*% (rnorm(q_fix) + mu_part) )
      
    }
    xbeta <- X %*% beta
    
    
    
    # --------------------- sdu update  -----
    if(method == "folded"){
      U_sign <- diag(sign(s)) %*% U  #U * sign(s)
      my_zu <- crossprod(tZ, crossprod(kron_cpp(U_sign, diag_q_rand_2), L))
    }else{
      my_zu <- zu
    }
    
    if(prior_rand=="SS"){
      for(k in 1:q_rand){
        y_tild <- y - zz_alpha - xbeta - zu[, -k] %*% sdu[-k]
        Sigma_sk <- 1/(crossprod(zu[, k], Gamma.inv %*% zu[, k])/se2 + 1/(se2 * tau2[2]))
        log_R <- as.numeric(log(pi_sdu)- log(1-pi_sdu) - 0.5 *log(crossprod(zu[, k], Gamma.inv %*% zu[, k]) + tau2[2]) -0.5 * log(tau2[2]) + (0.5 * t(y_tild) %*% Gamma.inv %*% zu[, k] %*% Sigma_sk %*% t(zu[, k]) %*% Gamma.inv %*% y_tild / (se2)^2 ))
        
        g_s[k] <- 1*(log(runif(1)) < (log_R - sumlogs(log_R, 0)))
        if(g_s[k]==0) {
          s[k] <- sdu[k] <- sdu_tild[k] <- 0
        }else{
          if(method == "truncated"){
            # sdu[k] <- sdu_tild[k] <- truncnorm::rtruncnorm(1, as.numeric(Sigma_sk *crossprod(zu[, k], Gamma.inv %*% y_tild))/se2, sqrt(as.numeric(Sigma_sk)), a=0)
            sdu_tild[k] <- rtnorm(1, as.numeric(Sigma_sk *crossprod(my_zu[, k], Gamma.inv %*% y_tild))/se2, sqrt(as.numeric(Sigma_sk)), A=0)
            if(is.na(sdu_tild[k]) | is.infinite(sdu_tild[k])){
              sdu_tild[k] <- 0
            }
            sdu[k] <- sdu_tild[k]
          }
          if(method== "folded"){
            s[k] <- rnorm(1, as.numeric(Sigma_sk *crossprod(my_zu[, k], Gamma.inv %*% y_tild))/se2, sqrt(as.numeric(Sigma_sk)))
            sdu_tild[k] <- sdu[k] <- abs(s[k])
          }
          if(method == "none"){
            sdu_tild[k] <- rnorm(1, as.numeric(Sigma_sk *crossprod(my_zu[, k], Gamma.inv %*% y_tild))/se2, sqrt(as.numeric(Sigma_sk)))
            sdu[k] <- abs(sdu_tild[k])
          }
        }
      }
      # print(round(s, 2))
      # pi_sdu <- rbeta(1, 1+sum(g_s), 1+q_rand-sum(g_s));
      # pi_sdu <- 0.1 #rbeta(1, 1+sum(g_s), 1+q_rand-sum(g_s));
      # print(pi_sdu)
      # print(sum(g_s > 0.5))
    }else{ # horseshoe or Cauchy
      if(method == "none"){
        if(prior_rand == "L"){ tmp_V <- omega2_sdu }else{ tmp_V <- tau2[2]*omega2_sdu }
        iV <- Diagonal(x = 1/tmp_V)
        svd_tmp <- svd(crossprod(zu, Gamma.inv %*% zu))
        id_tmp <- which(svd_tmp$d > 1e-5)
        svd_tmp$u <- svd_tmp$u[, id_tmp]
        svd_tmp$d <- svd_tmp$d[id_tmp]
        iS_part <- solve(svd_tmp$u %*% diag(svd_tmp$d) %*% t(svd_tmp$u) + iV)
        
        # iS_part <- solve(crossprod(zu, Gamma.inv %*% zu) + iV)
        mu_part <- crossprod(iS_part, crossprod(zu, Gamma.inv %*%(y - zz_alpha - xbeta)))
        # if(it == 10) browser()
        
        sdu_tild <- as.vector(mvnfast::rmvn(1, mu_part, se2*iS_part))
        # sdu_tild <- as.vector(tmvtnorm::rtmvnorm(1, as.numeric(mu_part), se2*iS_part, lower=rep(0, q_rand)))
        if(method == "folded"){
          s <- sdu_tild
          sdu <- sdu_tild <- abs(sdu_tild) 
        }else{
          sdu <- abs(sdu_tild)
        }
      }
      if(method == "folded"){
        # if(it == 10) browser()
        
        # U_sign <- U * sign(s)
        # my_zu <- crossprod(tZ, crossprod(kron_cpp(U_sign, diag_q_rand_2), L))
        
        if(prior_rand == "L"){ tmp_V <- omega2_sdu }else{ tmp_V <- tau2[2]*omega2_sdu }
        iV <- Diagonal(x = 1/tmp_V)
        svd_tmp <- svd(crossprod(my_zu, Gamma.inv %*% my_zu))
        id_tmp <- which(svd_tmp$d > 1e-5)
        svd_tmp$u <- svd_tmp$u[, id_tmp]
        svd_tmp$d <- svd_tmp$d[id_tmp]
        iS_part <- solve(svd_tmp$u %*% diag(svd_tmp$d) %*% t(svd_tmp$u) + iV)
        
        # iS_part <- solve(crossprod(zu, Gamma.inv %*% zu) + iV)
        mu_part <- crossprod(iS_part, crossprod(my_zu, Gamma.inv %*%(y - zz_alpha - xbeta)))
        # if(it == 10) browser()
        
        s <- as.vector(mvnfast::rmvn(1, mu_part, se2*iS_part))
        sdu <- sdu_tild <- abs(s) 
      }
      # if(method == "folded"){
      #   if(it == 10) browser()
      #   
      #   dim(U)
      #   sign(s)
      #   
      #   U[, 1:3]
      #   U_sign <- U * sign(s)
      #   U_sign[, 1:3]
      #   my_zu <- crossprod(tZ, crossprod(kron_cpp(U, diag_q_rand_2), L))
      #   
      #   if(prior_rand == "L"){ tmp_V <- omega2_sdu }else{ tmp_V <- tau2[2]*omega2_sdu }
      #   iV <- Diagonal(x = 1/tmp_V)
      #   svd_tmp <- svd(crossprod(zu, Gamma.inv %*% zu))
      #   id_tmp <- which(svd_tmp$d > 1e-5)
      #   svd_tmp$u <- svd_tmp$u[, id_tmp]
      #   svd_tmp$d <- svd_tmp$d[id_tmp]
      #   iS_part <- solve(svd_tmp$u %*% diag(svd_tmp$d) %*% t(svd_tmp$u) + iV)
      #   
      #   # iS_part <- solve(crossprod(zu, Gamma.inv %*% zu) + iV)
      #   mu_part <- crossprod(iS_part, crossprod(zu, Gamma.inv %*%(y - zz_alpha - xbeta)))
      #   # if(it == 10) browser()
      #   
      #   sdu_tild <- as.vector(mvnfast::rmvn(1, mu_part, se2*iS_part))
      #   
      #   s <- sdu_tild
      #   sdu <- sdu_tild <- abs(sdu_tild) 
      # }
      if(method == "truncated"){
        zuGI_FM <- crossprod(zu, Gamma.inv)
        zuGIzu_FM <- diag(zuGI_FM %*% zu)
        tmp_y <- y - zz_alpha - xbeta - zu %*% sdu_tild
        for(k in 1:q_rand){
          y_tild <- tmp_y + zu[, k] * sdu_tild[k]
          # Sigma_sk <- 1/(crossprod(zu[, k], Gamma.inv %*% zu[, k])/se2 + 1/(se2 * tau2[2]*omega2_sdu[k]))
          Sigma_sk <- 1/(zuGIzu_FM[k]/se2 + 1/(se2 * tau2[2]*omega2_sdu[k]))
          # sdu_tild[k] <- truncnorm::rtruncnorm(1, as.numeric(Sigma_sk *crossprod(zu[, k], Gamma.inv %*% y_tild))/se2,
          #                                                sqrt(as.numeric(Sigma_sk)), a=0)
          sdu_tild[k] <- rtnorm(1, as.numeric(Sigma_sk *crossprod(zu[, k], Gamma.inv %*% y_tild))/se2,
                                sqrt(as.numeric(Sigma_sk)), A=0)
          if(is.na(sdu_tild[k]) | is.infinite(sdu_tild[k])){
            # print(it)
            sdu_tild[k] <- 1e-4
          }
          sdu[k] <- sdu_tild[k]
        }
        # print(c(round(min(sdu), 6), round(se2, 2)))
        
      }
    }
    zus <- zu %*% sdu_tild
    par <- c(beta, sdu_tild)
    S <- diag(sdu_tild)
    
    # --------------------- tau2 update -----
    # fixed part
    if(prior_fixed=="ss"){
      tau2 <- pmax(1e-6, extraDistr::rinvgamma(1, sum(g_b)/2 + 1/2, 1/2 + 0.5*sum(beta^2)/se2))
    }else{
      tau2[1] <- max(1e-6, extraDistr::rinvgamma(1, 0.5*(q_fix + 1),  1/xi2[1] + 0.5*sum(beta^2/omega2_beta)/se2))
      xi2[1] <- extraDistr::rinvgamma(1, 1, 1 + 1/tau2[1])
    }
    
    # random part
    if(prior_rand == "SS"){
      tau2[2] <- pmax(1e-6, extraDistr::rinvgamma(1, sum(g_s)/2 + 1/2, 1/2 + 0.5*sum(sdu^2)/se2) )
    }
    
    if(prior_rand == "L"){
      tau2[2] <- extraDistr::rdgamma(1, shape = 1/2 + q_rand, scale = 1/2 + sum(omega2_sdu)/2)
    }
    
    if(prior_rand == "HS"){
      tau2[2] <- pmax(1e-6, extraDistr::rinvgamma(1, 0.5*(q_rand + 1), 1/xi2[2] + 0.5*sum(sdu^2 /omega2_sdu )/se2) )
      xi2[2] <- extraDistr::rinvgamma(1, 1, 1 + 1/tau2[2])
    }
    
    if(prior_rand == "NIG"){
      tau2[2] <- pmax(1e-6, extraDistr::rinvgamma(1, 0.5*(q_rand) + 1e-2, 1e-2 + 0.5*sum(sdu^2 /omega2_sdu )/se2) )
    }
    
    if(prior_rand == "cauchy"){}
    
    if(prior_rand == "DG"){}
    
    # --------------------- omega2 for beta and s -----
    # fixed part
    if(prior_fixed=="HS"){
      omega2_beta <- pmax(1e-3, extraDistr::rinvgamma(q_fix, 1, 1/nu2_beta + beta^2/(2*se2*tau2[1])))
      nu2_beta <- extraDistr::rinvgamma(q_fix, 1, 1 + 1/omega2_beta)
    }
    
    # random part
    if(prior_rand == "HS"){
      omega2_sdu <- pmax(1e-3, extraDistr::rinvgamma(q_rand, 1, 1/nu2_sdu + sdu^2/(2*se2*tau2[2])))
      nu2_sdu <- extraDistr::rinvgamma(q_rand, 1, 1 + 1/omega2_sdu)
    }
    
    if(prior_rand == "cauchy"){
      omega2_sdu <- pmax(1e-3, extraDistr::rinvgamma(q_rand, 1, 1/nu2_sdu + sdu^2/(2*se2*tau2[2])))
      # nu2_sdu <- extraDistr::rinvgamma(q_rand, 1, 1 + 1/omega2_sdu)
    }
    
    if(prior_rand == "L"){
      omega2_sdu  <- statmod::rinvgauss(q_rand, mean = tau2[2] * sqrt(se2)/abs(sdu), shape = tau2[2])
    }
    
    if(prior_rand == "NIG"){}
    
    if(prior_rand == "DG"){
      omega2_sdu <- GIGrvg::rgig(q_rand, a_s-1/2, a_s*kappa_2, sdu^2)
    }
    
    # --------------------- a_s and kappa2 for double gamma prior-----
    if(prior_rand == "DG"){
      # # MH step for a_s ______________
      before <- a_s
      # from shrinkTVP package Bitto & Frühwirth-Schnatter 2018
      a_s <- MH_step(a_s, epsilon_a_s, q_rand, kappa_2, sdu^2, b_w, nu=5, d1, d2)
      if(before != a_s) counter_a_s <- counter_a_s+1
      
      #_________________________
      
      kappa_2 <- rgamma(1, shape = d1+q_rand*a_s, rate = d2 + a_s/2*sum(omega2_sdu))
    }
    
    
    # ----- u  -----
    # Lambda <- kronecker(diag(n_niv), S)
    # iR <- kronecker(diag(n_niv), chol2inv(t(B)))
    # iS_part <- solve(crossprod(Z %*% Lambda, Gamma.inv %*% Z %*% Lambda)/se2 + iR)
    # mu_part <- crossprod(iS_part, crossprod(Z %*% Lambda, Gamma.inv %*% (y - zz_alpha -X%*%beta)/se2))
    # u <- mvnfast::rmvn(1, mu_part, iS_part)
    # U <- matrix(u, q_rand,  n_niv)
    
    
    # u <- svd_Sigma_u_cpp(Gamma_inv = Gamma_0_inv, B = B, sdu = sdu_tild, X = x, sigma2 = se2, y_tilde = y - zz_alpha -X%*%beta)
    u_tmp <- svd_Sigma_u_cpp_2(Gamma_inv =Gamma_0_inv, B = B, sdu = sdu_tild, X = x, sigma2 = se2, y_tilde = y - zz_alpha -X%*%beta)
    u <- u_tmp$u_t
    # u <- toto$u
    U <- matrix(u, q_rand,  n_niv)
    # u_tt <- u_tmp$u_tt
    
    
    zu <- crossprod(tZ, crossprod(kron_cpp(U, diag_q_rand_2), L))
    xzu <- cbind(X, zu)
    
    s_K <- 2
    r_K <- 0.03
    
    if(correlation_rand){
      i=2
      for(i in 1:length(id_bloc_R)){
        if(id_bloc_R[[i]] != 1){
          # ----- K ---------
          # K_i <- K[i]
          # l_p_K_i <- dgamma(K_i, s_K, r_K, log = TRUE)
          # 
          # l_p_R0_K_i <- 2*K_i*sum(log(sin(theta[[i]])))  + lcst_Norm(K_i, id_bloc_R[i]) #l_p_R0_K_fct(q_rand, theta, K_i)
          # 
          # if(adaptive){
          #   H <- 200
          #   if(it>H+10){
          #     K <- K_hist[[i]][(it-H-1) : (it-1), 1]
          #     K_tilde <- K - mean(K) 
          #     epsilon_K[k] <-  max(2.4/(H+1) * c(crossprod(K_tilde)), 1)
          #     K_hist_epsilon[[k]][it, 1] <- epsilon_K[k]
          #   }
          # }
          # 
          # b_inf <- 0
          # K_star <- truncnorm::rtruncnorm(1, a=b_inf, b=100, mean = K_i, sd=sqrt(epsilon_K[i]))
          # l_p_star_old <- log(truncnorm::dtruncnorm(K_star, a = b_inf, b=100, mean = K_i, sd=sqrt(epsilon_K[i])))
          # l_p_old_star <- log(truncnorm::dtruncnorm(K_i, a = b_inf, b=100, mean = K_star, sd=sqrt(epsilon_K[i])))
          # 
          # l_p_K_star  <- dgamma(K_star, s_K, r_K, log = TRUE)
          # l_p_R0_K_star <-  2*K_star*sum(log(sin(theta[[i]])))  + lcst_Norm(K_star, id_bloc_R[i]) #l_p_R0_K_fct(q_rand, theta, K_star)
          # 
          # log_ratio <- min(0, l_p_R0_K_star + l_p_K_star - l_p_R0_K_i - l_p_K_i + l_p_old_star - l_p_star_old)
          # 
          # #browser()
          # if(runif(1) < exp(log_ratio)){
          #   # print(round(K_star, 2))
          #   K[i] <- K_star
          #   if(it > burnin) counter_1_K <- counter_1_K+1
          # }
          # if(it > burnin) counter_tot_K <- counter_tot_K+1
          # 
          
          
          # print(K)
          # ----- R0 -----
          #correlation R0 = B B'     B= f(theta)
          # if(it==2){
          l_p_theta_i     <- 0 # sum(index * log(sin(theta )))
          # l_p_Y_theta_i    <- sum(mvnfast::dmvn(t(U), rep(0, q_rand), t(B), log=TRUE, isChol=TRUE))
          l_p_Y_theta_i <- l_p_Y_theta_star_cpp_3(n_niv, id_bloc_R[i], U[id_bloc_vec == i, ], theta[[i]], B_list[[i]], K[i])$log_density
          
          # browser()
          # l_p_Y_theta_i <- l_p_Y_theta_star_cpp_4(n_niv, id_bloc_R[i], S[id_bloc_vec == i, id_bloc_vec == i] %*% U[id_bloc_vec == i, ], theta[[i]], B_list[[i]], K[i])$log_density
          
          k=1
          for(k in sample(1:length(theta[[i]]))){
            # print(c(i, j))
            theta_star <- theta[[i]]
            
            if(adaptive){
              H <- 200
              if(it>H+10) {
                K <- theta_hist[[i]][(it-H-1) : (it-1), k]
                K_tilde <- K - mean(K) 
                epsilon_angle[i] <-  max(2.4/(H+1) * c(crossprod(K_tilde)), 1)
                theta_hist_epsilon[[i]][it, k] <- epsilon_angle[i]
              }
            }
            # theta_star[k] <- runif(1, max(0, theta_star[k]-6*epsilon_angle), min(pi, theta_star[k]+6*epsilon_angle))
            # l_p_star_old <- l_p_old_star <- 0
            
            theta_star[k] <- truncnorm::rtruncnorm(1, a=0, b=pi, mean = theta[[i]][k], sd=sqrt(epsilon_angle[i]))
            l_p_star_old <- log(truncnorm::dtruncnorm(theta_star[k], a = 0, b=pi, mean = theta[[i]][k], sd=sqrt(epsilon_angle[i])))
            l_p_old_star <- log(truncnorm::dtruncnorm(theta[[i]][k], a = 0, b=pi, mean = theta_star[k], sd=sqrt(epsilon_angle[i])))
            
            l_p_theta_star  <- 0 # sum(index * log(sin(theta_star)))
            
            theta_mat_star <- theta_mat[[i]]; theta_mat_star[df_id_theta[[i]]$id_row[k], df_id_theta[[i]]$id_col[k]] <-  theta_star[k] 
            
            B_star  <- corrPourhamandi_new(B_list[[i]], theta_mat_star, df_id_theta[[i]]$id_row[k], df_id_theta[[i]]$id_col[k])
            # l_p_Y_theta_star <- sum(mvnfast::dmvn(t(U), rep(0, q_rand), t(B_star), log=TRUE, isChol=TRUE))
            # l_p_Y_theta_star <- l_p_Y_theta_star_cpp_2(n_niv, q_rand, U, theta_star[lower.tri(theta_star)], B_star, K)
            tmp <- l_p_Y_theta_star_cpp_3(n_niv, id_bloc_R[i], U[id_bloc_vec == i, ], theta_star[lower.tri(theta_star)], B_star, K[i])
            l_p_Y_theta_star <- tmp$log_density
            iR0_star <- tmp$iR
            # l_p_Y_theta_star <- matrixNormal::dmatnorm(matrix(y, m,  n_niv), matrix(zz_alpha - X %*% beta - kronecker(diag_niv, x%*%S%*%B_star)%*%u_tt, m,  n_niv), diag_ind, Gamma_0*se2, log=TRUE)
            # l_p_Y_theta_star <- log_d_Y_cpp(matrix(y, m,  n_niv), matrix(zz_alpha - X %*% beta - kronecker(diag_niv, x%*%S%*%B_star)%*%u_tt, m,  n_niv), diag_ind, Gamma_0_inv/se2)
            # l_p_Y_theta_star <- log_d_Y_cpp_2(matrix(y, m,  n_niv), matrix(zz_alpha - X %*% beta - kronecker(diag_niv, x%*%S%*%B_star)%*%u_tt, m,  n_niv), diag_ind, Gamma_0_inv/se2)
            
            # if(it == 20) browser()
            
            
            #________________________________
            
            log_ratio <- min(0, l_p_Y_theta_star + l_p_theta_star - l_p_Y_theta_i - l_p_theta_i + l_p_old_star - l_p_star_old)
            
            if(runif(1) < exp(log_ratio)){
              theta[[i]][k] <-  theta_star[k]
              theta_mat[[i]][df_id_theta[[i]]$id_row[k], df_id_theta[[i]]$id_col[k]] <-  theta_star[k] 
              l_p_Y_theta_i <- l_p_Y_theta_star
              l_p_theta_i <- l_p_theta_star
              B_list[[i]] <- B_star
              iR0_list[[i]] <- iR0_star # my_chol2inv(B)
              if(it > burnin) counter_1 <- counter_1+1
            }
            if(it > burnin) counter_tot <- counter_tot+1
            
          }
        }
      }
      B <- as.matrix(Matrix::bdiag(B_list))
      iR0 <- as.matrix(Matrix::bdiag(iR0_list))
    }
    
    
    # print(round(K, 2))
    
    # ----- rho  -----
    # rho <- 0
    # diag_0_Gamma =  arma::ones(n_niv, 1) * (1+rho*rho)
    # diag_0_Gamma(0) = diag_0_Gamma(n_niv-1) = 1
    # diag_1_Gamma = -rho* arma::ones(n_niv-1, 1)
    # Gamma_0_inv = (arma::diagmat(diag_1_Gamma, -1) + arma::diagmat(diag_0_Gamma, 0) + arma::diagmat(diag_1_Gamma, 1))/(1-rho*rho);
    
    # rho <- 0
    # Gamma_0 <- diag(n_niv)
    # Gamma_0_inv <- diag(n_niv)
    # Gamma.inv <- diag(n_niv*m)
    
    if(AR_residual){
      # Metropolis-Hastings
      tp <- (y - zz_alpha - xzu %*% par)
      tp_mat <- matrix(tp, m, n_niv)
      tmp_MH_rho <- crossprod(tp, crossprod(Gamma.inv, tp))
      det_gamma <- determinant(Gamma_0, logarithm = TRUE)$modulus
      
      rho_star <- runif(1, max(-1, rho-epsilon_rho),min(1, rho+epsilon_rho))
      log_Qnew = dunif(rho_star,max(-1, rho-epsilon_rho), min(1,rho+epsilon_rho), log = T) # proposal distribution Q
      log_Qold = dunif(rho,max(-1, rho_star-epsilon_rho), min(1,rho_star+epsilon_rho), log = T)
      
      Gamma_0_star <- matrix(rho_star, n_niv, n_niv); for(i in 1:n_niv) for(j in 1:n_niv) Gamma_0_star[i, j] <- rho_star^abs(i-j)
      Gamma_0_star_inv <- update_Gamma_inv(diag_niv, rho_star)
      # Gamma_star.inv <- kron_cpp(Gamma_0_star_inv, diag_ind)
      # tmp_star <- crossprod(tp, Gamma_star.inv) %*% tp
      tmp_star <- y_Gamma_inv_y_cpp(tp_mat, Gamma_0_star_inv)
      
      det_gamma_0_star <- determinant(Gamma_0_star, logarithm = TRUE)$modulus
      
      log_ratio <- -ni[1]/2 *  det_gamma_0_star -
        1/(2*se2) *  tmp_star + log_Qold +
        ni[1]/2 * det_gamma +
        1/(2*se2)*  tmp_MH_rho - log_Qnew
      
      
      if(log(runif(1)) < (as.numeric(log_ratio))) {
        tmp_MH_rho <- tmp_star
        det_gamma <- det_gamma_0_star
        rho <- rho_star
        # Gamma <- Gamma_star
        Gamma_0 <- Gamma_0_star
        Gamma_0_inv <- Gamma_0_star_inv
        Gamma.inv <- kron_cpp_s_s(Gamma_0_star_inv, diag_ind)
        # Gamma.inv <- kronecker(Gamma_0_star_inv, diag_ind)
        counter_rho <- counter_rho+1
      }
    }
    
    y_tmp <- y - zz_alpha - X %*% beta - zu %*% sdu_tild
    
    if(!AR_residual){
      # browser()
      for(k in 1:(n_niv-1)){
        sd_t[k] <- sqrt(extraDistr::rinvgamma(1, 1e-3 + m/2, as.numeric(1e-3 + 0.5*crossprod(y_tmp[niv == k])/se2)))
        if(is.null(sd_t[k]) | is.na(sd_t[k]) | sd_t[k] < sqrt(1e-3)) sd_t[k] <- sqrt(1e-3)
      }
      Gamma_0 <- diag(sd_t^2)
      Gamma_0_inv <- Diagonal(x=1/sd_t^2)
      Gamma.inv <- kron_cpp_s_s(Gamma_0_inv, diag_ind)
    }
    
    
    
    # ----- se2  -----
    # se2 <- 1
    # tp_mat <- matrix(y - zz_alpha - X %*% beta - zu %*% sdu_tild, m, n_niv)
    # ytilde <- 0.5 * y_Gamma_inv_y_cpp(tp_mat, Gamma_0_inv) 
    
    ytilde <- 0.5 * crossprod(y_tmp, Gamma.inv %*% y_tmp)
    
    # if(individual_effect){
    #   if(prior_rand == "L"){
    #     se2 <- extraDistr::rinvgamma(1, 1e-3 + 0.5*(no + sum(g_b) + sum(g_s)+1),
    #                                  as.numeric(1e-3 + 0.5*crossprod(alpha)/sd_alpha^2 + ytilde +
    #                                               0.5*sum(beta^2/omega2_beta) +
    #                                               0.5*sum( sdu^2/omega2_sdu ) ))
    #     
    #   }else{
    #     se2 <- extraDistr::rinvgamma(1, 1e-3 + 0.5*(no + sum(g_b) + sum(g_s)+1),
    #                                  as.numeric(1e-3 + 0.5*crossprod(alpha)/sd_alpha^2 + ytilde +
    #                                               0.5*sum(beta^2/omega2_beta)/tau2[1] +
    #                                               0.5*sum( sdu^2/omega2_sdu )/tau2[2] ))
    #   }
    # }else{
    if(prior_rand == "L"){
      se2 <- extraDistr::rinvgamma(1, 1e-3 + 0.5*(no + sum(g_b) + sum(g_s)),
                                   as.numeric(1e-3 + ytilde +
                                                0.5*sum(beta^2/omega2_beta) +
                                                0.5*sum( sdu^2/omega2_sdu ) ))
      
    }else{
      se2 <- extraDistr::rinvgamma(1, 1e-3 + 0.5*(no + sum(g_b) + sum(g_s)),
                                   as.numeric(1e-3 + ytilde +
                                                0.5*sum(beta^2/omega2_beta)/tau2[1] +
                                                0.5*sum( sdu^2/omega2_sdu )/tau2[2] ))
    }
    # }
    
    # print(se2)
    if(it > burnin & it %% thinin == 0){
      # chain$mu[ii] <- mu
      if(individual_effect) chain$alpha[ii, ] <- alpha
      if(individual_effect) chain$sd_alpha[ii] <- sd_alpha
      chain$beta[ii, ] <- beta
      if(prior_fixed=="SS") chain$g_b[ii, ] <- g_b
      if(prior_fixed=="SS") chain$pi_beta[ii] <- pi_beta
      if(method == "folded") chain$s[ii, ] <- s
      chain$sdu[ii, ] <- sdu
      chain$sdu_tild[ii, ] <- sdu_tild
      if(prior_rand=="SS")  chain$g_s[ii, ] <- g_s
      if(prior_rand=="SS")  chain$pi_sdu[ii] <- pi_sdu
      chain$se2[ii] <- se2
      chain$rho[ii] <- rho
      if(prior_rand=="DG")  chain$a_s[ii] <- a_s
      if(prior_rand=="DG")  chain$kappa2[ii] <- kappa_2
      chain$omega2_beta[ii, ] <- omega2_beta
      chain$omega2_sdu[ii, ] <- omega2_sdu
      if(correlation_rand) for(i in 1:length(id_bloc_R)) chain$theta[[i]][ii, ] <- theta[[i]]
      if(correlation_rand) chain$R0[ii, ] <- tcrossprod(B)[lower.tri(B)]
      if(correlation_rand) chain$K[ii, ] <- K
      chain$u[ii, ] <- c(U)
      chain$tau2[ii, ] <- tau2
      # chain$sd_t[ii, ] <- se2_t
      ii <- ii+1
    }
    if(adaptive){
      for(i in 1:length(id_bloc_R)) theta_hist[[i]][it, ] <- theta[[i]]; 
      for(i in 1:length(id_bloc_R)) K_hist[[i]] [it, ] <- K[i]}
  } 
  if(adaptive) {chain$theta_hist_epsilon <- theta_hist_epsilon; chain$K_hist_epsilon <- K_hist_epsilon}
  
  if(correlation_rand){
    chain$acc_rate_theta <- counter_1/counter_tot
    print("Acceptance ratio of theta:")
    print(counter_1/counter_tot)
  }
  if(prior_rand=="DG"){
    print("Acceptance ratio of a_s:")
    print(counter_a_s/niter)
    chain$acc_rate_a_s <- counter_a_s/niter
  }
  print("Acceptance ratio of rho:")
  print(counter_rho/niter)
  chain$acc_rate_rho <- counter_rho/niter
  return(chain)
}













#  Old version

# ris <- function(y,
#                 x,
#                 ind,
#                 prior_fixed = "HS",
#                 prior_rand="HS",
#                 correlation = TRUE,
#                 individual_effect = FALSE,
#                 niter=10000,
#                 burnin=round(niter / 2),
#                 thinin=1,
#                 epsilon=0.1,
#                 epsilon_a_s=1,
#                 K=1)
# {
# 
#   if(!prior_fixed  %in% c("SS", "HS")) stop("prior_fixed choice must be 'SS'or 'HS'." )
#   if(!prior_rand  %in% c("SS", "HS", "DG", "cauchy", "NIG", "L")) stop("prior_rand choice must be 'SS', 'HS', 'DG', 'cauchy', 'L' or 'NIG'." )
#   if(burnin >= niter ) stop("'burnin' must be lower than 'niter'")
# 
#   print("Variance components selection in RIS models")
#   # if(prior_rand=="HS") print("Horseshoe prior for variance component in RIS models")
#   # if(prior_rand=="cauchy") print("Cauchy prior for variance component in RIS models")
#   # if(prior_rand=="NIG") print("Gaussian prior for variance component in RIS models")
# 
#   # Settings ________________________________________________________________
#   levels(ind) <- 1:length(unique(levels(ind)))
#   n_niv         <- length(levels(ind))
#   no          <- length(y)
#   ni          <- table(ind)
#   q_fix       <- ncol(x)
#   q_rand      <- q_fix + 1
# 
#   index <- unlist(lapply(1:(q_fix), function(l) 2*K+rep(q_rand - l, q_rand - l)))
# 
#   #sort by increasing factor
#   y       <- y[order(ind)]
#   x       <- x[order(ind), ]
#   ind <- sort(ind)
#   J       <- t(as(ind,  Class = "sparseMatrix"))
#   L       <- matrix(0, (q_rand)^2, q_rand)
#   pos     <- cbind(which(diag(q_rand)>0), 1:(q_rand))
#   L[pos]  <- 1
# 
#   if(individual_effect){
#     table_ind <- table(ind)
#     zz <- matrix(0, length(y), max(table_ind))
#     for(k in 1:n_niv) zz[ind==k, 1:table_ind[k]] <- diag(table_ind[k])
#   }
# 
#   # Outputs _________________________________________________________________
#   chain             <- list()
#   if(individual_effect) chain$alpha       <- matrix(0, floor(niter - burnin)/thinin, ncol(zz))
#   if(individual_effect) chain$sd_alpha    <- rep(NA, floor(niter - burnin)/thinin)
#   # chain$mu          <- rep(NA, floor(niter - burnin)/thinin)
#   chain$beta        <- matrix(0, floor(niter - burnin)/thinin, q_fix)
#   chain$g_b         <- matrix(0, floor(niter - burnin)/thinin, q_fix)
#   chain$pi_beta     <- rep(NA, floor(niter - burnin)/thinin)
#   chain$sdu         <- matrix(0, floor(niter - burnin)/thinin, q_rand)
#   chain$sdu_tild    <- matrix(0, floor(niter - burnin)/thinin, q_rand)
#   chain$g_s         <- matrix(0, floor(niter - burnin)/thinin, q_rand)
#   chain$pi_sdu     <- rep(NA, floor(niter - burnin)/thinin)
#   chain$tau2 <- matrix(NA, floor(niter - burnin)/thinin, 2) # by convention tau[, 1]=tau_beta et tau[, 2]=tau_sdu
#   chain$omega2_beta <- matrix(1, floor(niter - burnin)/thinin, q_fix)
#   chain$omega2_sdu  <- matrix(1, floor(niter - burnin)/thinin, q_rand)
#   chain$theta       <- matrix(0, floor(niter - burnin)/thinin, 0.5*q_rand*(q_rand-1))
#   chain$se2         <- rep(NA, floor(niter - burnin)/thinin)
#   chain$u           <- matrix(NA, floor(niter - burnin)/thinin, n_niv*q_rand)
#   if(prior_rand=="DG") chain$a_s        <- rep(NA, floor(niter - burnin)/thinin)
#   if(prior_rand=="DG") chain$kappa2     <- rep(NA, floor(niter - burnin)/thinin)
# 
# 
#   # Initialisation __________________________________________________________
#   mu <- 0 # mean(y)
#   if(individual_effect) alpha <- rep(0, ncol(zz))
#   if(individual_effect) sd_alpha <- 10000
#   U                 <- matrix(rnorm(q_rand*n_niv), q_rand, n_niv)
#   u                 <- c(U)
#   # u <- init$u
#   # U               <- matrix(u, q_rand, n_niv)
#   beta <- chain$beta[1, ] <- rnorm(q_fix)
#   sdu_tild <- rnorm(q_rand)
#   sdu <- abs(sdu_tild)
#   se2 <- chain$se2[1] <- 1 # rgamma(1, 2, 2)
#   pi_beta <- pi_sdu <- 0.5
# 
#   if(correlation) theta <- sample_thetaPourhamandi(q_rand) else theta <- matrix(0, q_rand, q_rand); theta[lower.tri(theta)] <- pi/2
#   chain$theta[1, ] <- theta[lower.tri(theta)]
#   if(correlation){
#     B <- corrPourhamandi(theta)
#     R0 <- tcrossprod(B)
#   }else{
#     B <- R0 <- diag(q_rand)
#   }
# 
# 
#   tau2 <- chain$tau2[1, ] <- c(1, 1)
#   if(prior_fixed=="SS"){
#     g_b <- chain$g_b[1, ] <- sample(c(0, 1), q_fix, replace = TRUE)
#     beta[g_b==0] <- 0
#     chain$beta[1, ] <- beta
#   }else{ # Horseshoe
#     g_b <- chain$g_b[1, ] <- rep(1, q_fix)
#     tau2 <- chain$tau2[1, ] <- c(1e-5, 1)
#   }
#   xi2 <- rep(1, 2)
#   omega2_beta <- chain$omega2_beta[1, ] <- rep(1, q_fix)
#   nu2_beta <- rep(1, q_fix)
#   omega2_sdu <- chain$omega2_sdu[1, ] <- rep(1, q_rand)
#   nu2_sdu <- rep(1, q_rand)
#   if(prior_rand == "cauchy"){
#     nu2_sdu <- rep(2, q_rand)
#   }
#   if(prior_rand=="SS"){
#     g_s <- chain$g_s[1, ] <- sample(c(0, 1), q_rand, replace = TRUE)
#     sdu[g_s == 0] <- 0
#   }else{
#     g_s <- chain$g_s[1, ] <- rep(1, q_rand)
#   }
#   if(prior_rand == "DG"){
#     d1 <- 0.001
#     d2 <- 0.001
#     kappa_2 <- 20 #rgamma(1, shape = d1, rate = d1)
# 
#     b_w <- 10
#     a_s <- 0.1 # rexp(1, b_w)
#     l_p_omega_sdu_i <- sum(dgamma(omega2_sdu, a_s, a_s*kappa_2/2, log = TRUE))
#     l_p_a_s_i <- dexp(a_s, b_w, log=TRUE)
#   }
# 
#   #___________________________________________________________
# 
#   diag_niv <- Diagonal(n_niv)
#   diag_q_rand <- Diagonal(q_rand)
#   diag_q_rand_2 <- diag(q_rand)
#   tZ <- KhatriRao(t(J),  t(cbind(1, x)))
#   Z <- t(tZ)
#   zu <- crossprod(tZ, crossprod(kron_cpp(U, diag_q_rand_2), L))
#   zus <- zu %*% sdu_tild
#   xzu <- cbind(x, zu)
#   iR <- chol2inv(t(B))
#   id_lower.tri <- lower.tri(theta)
#   counter_tot <- counter_1 <- counter_a_s <- 0
#   z_list <- list(); for(k in 1:n_niv) z_list[[k]] <- cbind(1, x[ind==k, ])
#   if(individual_effect){
#     zz_alpha <- zz %*% alpha
#     zztzz <- crossprod(zz)
#     svd_zztzz <- svd(zztzz)
#   }else{
#     zz_alpha <- rep(0, length(y))
#   }
#   par <- c(beta, sdu_tild)
# 
#   ii <- 1
# 
#   # MCMC ___________________________________________________________
# 
#   print("0 %")
#   for(it in 2:niter){
#     # print(it)
#     if(it %% (niter/10) == 0) print(paste((it %/% (niter/10))*10, "%"))
#     # if(it == 100) browser()
# 
#     # --------------------- update mu
#     # mu <- rnorm(1, mean(y-zz_alpha - xzu %*% par), sd = sqrt(se2/no))
# 
#     # --------------------- update alpha (individual effect)
#     if(individual_effect){
#       iS_alpha <- diag(1/diag(zztzz/se2 + diag(ncol(zz))/sd_alpha^2))
#       mu_alpha <- crossprod(iS_alpha, crossprod(zz, y - xzu %*% par)/se2)
#       alpha <- c(mvnfast::rmvn(1, mu_alpha, iS_alpha))
#       zz_alpha <- zz %*% alpha
#       sd_alpha <-  sqrt(1/rgamma(1, 1/2+ncol(zz)/2, rate=1/2 + crossprod(alpha)/2))
#     }
# 
#     # browser()
#     # --------------------- update beta
#     if(prior_fixed=="SS"){
#       for(k in 1:q_fix){
#         y_tild <- y - zz_alpha - x[, -k] %*% beta[-k] - zus
#         Sigma_bk <- 1/(crossprod(x[, k])/se2 + 1/(se2 * tau2[1]))
#         log_R <- as.numeric(log(pi_beta)- log(1-pi_beta) - 0.5 *log(crossprod(x[, k]) + tau2[1]) -0.5 * log(tau2[1]) + (0.5 * t(y_tild) %*% x[, k] %*% Sigma_bk %*% t(x[, k]) %*% y_tild / (se2)^2 ))
#         g_b[k] <- 1*(log(runif(1)) < log_R - sumlogs(log_R, 0))
# 
#         if(g_b[k]==0) {
#           beta[k] <- 0
#         }else{
#           beta[k] <- rnorm(1, as.numeric(Sigma_bk *crossprod(x[, k], y_tild))/se2, sqrt(Sigma_bk))
#         }
#       }
#       pi_beta <- rbeta(1+sum(g_b), 1+q_fix-sum(g_b));
#     }else{ # horseshoe
#       if(prior_rand == "L"){ tmp_V <- c(omega2_beta) }else{ tmp_V <- c(tau2[1]*omega2_beta) }
#       iV <- Diagonal(x = 1/tmp_V)
#       iS_part <- solve(crossprod(x) + iV)
#       mu_part <- crossprod(iS_part, crossprod(x, y - zz_alpha - zus))
#       beta <- as.vector(mvnfast::rmvn(1, mu_part, se2*iS_part))
#     }
#     xbeta <- x %*% beta
# 
# 
# 
#     # --------------------- update  sdu
#     if(prior_rand=="SS"){
#       for(k in 1:q_rand){
#         y_tild <- y - zz_alpha - xbeta - zu[, -k] %*% sdu_tild[-k]
#         Sigma_sk <- 1/(crossprod(zu[, k])/se2 + 1/(se2 * tau2[2]))
#         log_R <- as.numeric(log(pi_sdu)- log(1-pi_sdu) - 0.5 *log(crossprod(zu[, k]) + tau2[2]) -0.5 * log(tau2[2]) + (0.5 * t(y_tild) %*% zu[, k] %*% Sigma_sk %*% t(zu[, k]) %*% y_tild / (se2)^2 ))
#         g_s[k] <- 1*(log(runif(1)) < log_R- sumlogs(log_R, 0))
#         if(g_s[k]==0) {
#           sdu[k] <- sdu_tild[k] <- 0
#         }else{
#           sdu_tild[k] <- rnorm(1, as.numeric(Sigma_sk *crossprod(zu[, k], y_tild))/se2, sqrt(Sigma_sk))
#           sdu[k] <- abs(sdu_tild[k])
#         }
#       }
#       pi_sdu <- rbeta(1+sum(g_s), 1+q_rand-sum(g_s));
# 
#     }else{ # horseshoe
#       if(prior_rand == "L"){ tmp_V <- omega2_sdu }else{ tmp_V <- tau2[2]*omega2_sdu }
#       iV <- Diagonal(x = 1/tmp_V)
#       iS_part <- solve(crossprod(zu) + iV)
#       mu_part <- crossprod(iS_part, crossprod(zu, y - zz_alpha - xbeta))
#       sdu_tild <- as.vector(mvnfast::rmvn(1, mu_part, se2*iS_part))
#       sdu <- abs(sdu_tild)
# # sdu_tild <- abs(sdu_tild)  # <<<< !!!!!!!!!!!
#     }
#     zus <- zu %*% sdu_tild
#     # zus <- zu %*% sdu
#     par <- c(beta, sdu_tild)
# 
# 
#     # --------------------- update u
#     # Lambda <- kron_cpp_s_s(diag_niv, Diagonal(x = sdu_tild))
#     # iR <- kron_cpp(diag(n_niv), chol2inv(t(B)))
#     # iS_part <- solve(crossprod(Z %*% Lambda)/se2 + iR)
#     # mu_part <- crossprod(iS_part, crossprod(Z %*% Lambda, (y - zz_alpha -x%*%beta)/se2))
#     # u <- mvnfast::rmvn(1, mu_part, iS_part)
#     # U <- matrix(u, q_rand,  n_niv)
# 
#     k=1
#     # iR <- chol2inv(t(B))
#     for(k in 1:n_niv){
#       tmp_0 <- z_list[[k]] * rep(sdu_tild, each = nrow(z_list[[k]] )) # element wise product
# 
#       tmp_1 <- crossprod(tmp_0)/se2
#       iS_part <- solve(tmp_1 + iR)
#       # mu_part <- crossprod(iS_part, crossprod(tmp_0, (y[ind==k]  - zz_alpha[ind==k] -x[ind==k,]%*%beta)/se2))
#       #  mu_part <- iS_part %*% t(tmp_0) %*% (y[ind==k] - zz_alpha[ind==k] -x[ind==k,]%*%beta)/se2
#       if(individual_effect){ # to speed up
#         mu_part <- my_mu_part(iS_part, tmp_0, (y[ind==k] - zz_alpha[ind==k] -x[ind==k,]%*%beta)/se2)
#       }else{
#         mu_part <- my_mu_part(iS_part, tmp_0, (y[ind==k] -x[ind==k,]%*%beta)/se2)
#       }
# 
#       u_tmp <- mvnfast::rmvn(1, mu_part, iS_part)
#       U[, k] <- c(u_tmp)
#     }
#     zu <- crossprod(tZ, crossprod(kron_cpp(U, diag_q_rand_2), L))
#     xzu <- cbind(x, zu)
# 
# 
# 
#     # --------------------- update tau2
#     # fixed part
#     if(prior_fixed=="ss"){
#       tau2 <- pmax(1e-6, extraDistr::rinvgamma(1, sum(g_b)/2 + 1/2, 1/2 + 0.5*sum(beta^2)/se2))
#     }else{
#       tau2[1] <- max(1e-6, extraDistr::rinvgamma(1, 0.5*(q_fix + 1),  1/xi2[1] + 0.5*sum(beta^2/omega2_beta)/se2))
#       xi2[1] <- extraDistr::rinvgamma(1, 1, 1 + 1/tau2[1])
#     }
# 
#     # random part
#     if(prior_rand == "SS"){
#       tau2[2] <- pmax(1e-6, extraDistr::rinvgamma(1, sum(g_s)/2 + 1/2, 1/2 + 0.5*sum(sdu^2)/se2) )
#     }
# 
#     if(prior_rand == "L"){
#       tau2[2] <- extraDistr::rdgamma(1, shape = 1/2 + q_rand, scale = 1/2 + sum(omega2_sdu)/2)
#     }
# 
#     if(prior_rand == "HS"){
#       tau2[2] <- pmax(1e-6, extraDistr::rinvgamma(1, 0.5*(q_rand + 1), 1/xi2[2] + 0.5*sum(sdu^2 /omega2_sdu )/se2) )
#       xi2[2] <- extraDistr::rinvgamma(1, 1, 1 + 1/tau2[2])
#     }
# 
#     if(prior_rand == "NIG"){
#       tau2[2] <- pmax(1e-6, extraDistr::rinvgamma(1, 0.5*(q_rand) + 1e-2, 1e-2 + 0.5*sum(sdu^2 /omega2_sdu )/se2) )
#     }
# 
#     if(prior_rand == "cauchy"){}
# 
#     if(prior_rand == "DG"){}
# 
#     # --------------------- omega2 for beta and s
#     # fixed part
#     if(prior_fixed=="HS"){
#       omega2_beta <- pmax(1e-3, extraDistr::rinvgamma(q_fix, 1, 1/nu2_beta + beta^2/(2*se2*tau2[1])))
#       nu2_beta <- extraDistr::rinvgamma(q_fix, 1, 1 + 1/omega2_beta)
#     }
# 
#     # random part
#     if(prior_rand == "HS"){
#       omega2_sdu <- pmax(1e-3, extraDistr::rinvgamma(q_rand, 1, 1/nu2_sdu + sdu^2/(2*se2*tau2[2])))
#       nu2_sdu <- extraDistr::rinvgamma(q_rand, 1, 1 + 1/omega2_sdu)
#     }
# 
#     if(prior_rand == "cauchy"){
#       omega2_sdu <- pmax(1e-3, extraDistr::rinvgamma(q_rand, 1, 1/nu2_sdu + sdu^2/(2*se2*tau2[2])))
#       # nu2_sdu <- extraDistr::rinvgamma(q_rand, 1, 1 + 1/omega2_sdu)
#     }
# 
#     if(prior_rand == "L"){
#       omega2_sdu  <- statmod::rinvgauss(q_rand, mean = tau2[2] * sqrt(se2)/abs(sdu), shape = tau2[2])
#     }
# 
#     if(prior_rand == "NIG"){}
# 
#     if(prior_rand == "DG"){
#       omega2_sdu <- GIGrvg::rgig(q_rand, a_s-1/2, a_s*kappa_2, sdu^2)
#     }
# 
#     # --------------------- a_s and kappa2 for double gamma prior
# 
#     if(prior_rand == "DG"){
#       # # MH step for a_s ______________
#       before <- a_s
#       # from shrinkTVP package Bitto & Frühwirth-Schnatter 2018
#       a_s <- MH_step(a_s, epsilon_a_s, q_rand, kappa_2, sdu^2, b_w, nu=5, d1, d2)
#       if(before != a_s) counter_a_s <- counter_a_s+1
# 
#       #_________________________
# 
#       kappa_2 <- rgamma(1, shape = d1+q_rand*a_s, rate = d2 + a_s/2*sum(omega2_sdu))
#     }
# 
# 
#     # --------------------- correlation R0 = B B'     B= f(theta)
#     if(correlation){
#       if(it==2){
#         l_p_theta_i     <- sum(index * log(sin(theta[id_lower.tri] )))
#         l_p_Y_theta_i    <- sum(mvnfast::dmvn(t(U), rep(0, q_rand), t(B), log=TRUE, isChol=TRUE))
#         # l_p_Y_theta_i <- l_p_Y_theta_star_cpp_2(n_niv, q_rand, U, theta[lower.tri(theta)], B, 0)
#       }
#       for(i in sample(2:q_rand, q_rand-1)){
#         for(j in sample(1:(i-1), i-1)){
#           # print(c(i, j))
#           theta_star <- theta
#           # theta_star[i, j] <- runif(1, max(0, theta_star[i, j]-epsilon), min(pi, theta_star[i, j]+epsilon))
#           # l_p_star_old <- l_p_old_star <- 0
#           theta_star[i, j] <- truncnorm::rtruncnorm(1, a=0, b=pi, mean = theta[i, j], sd=epsilon)
#           l_p_star_old <- log(truncnorm::dtruncnorm(theta_star[i, j], a = 0, b=pi, mean = theta[i, j], sd=epsilon))
#           l_p_old_star <- log(truncnorm::dtruncnorm(theta[i, j], a = 0, b=pi, mean = theta_star[i, j], sd=epsilon))
# 
#           l_p_theta_star  <- sum(index * log(sin(theta_star[id_lower.tri])))
#           B_star  <- corrPourhamandi_new(B, theta_star, i, j)
#           l_p_Y_theta_star <- sum(mvnfast::dmvn(t(U), rep(0, q_rand), t(B_star), log=TRUE, isChol=TRUE))
# 
#           # iR_star <- chol2inv(t(B_star))
#           # iR_star <- my_chol2inv(B_star)
#           # l_p_Y_theta_star <- l_p_Y_theta_star_cpp_2(n_niv, q_rand, U, theta_star[lower.tri(theta_star)], B_star, 0)
# 
#           #________________________________
# 
#           log_ratio <- min(0, l_p_Y_theta_star + l_p_theta_star - l_p_Y_theta_i - l_p_theta_i + l_p_old_star - l_p_star_old)
#           # log_ratio <- min(0, l_p_Y_theta_star - l_p_Y_theta_i) #  + l_p_theta_star - l_p_theta_i + l_p_old_star - l_p_star_old)
# 
#           #browser()
#           if(runif(1) < exp(log_ratio)){
#             theta[i, j] <-  theta_star[i, j]
#             l_p_Y_theta_i <- l_p_Y_theta_star
#             l_p_theta_i <- l_p_theta_star
#             B <- B_star
#             iR <- my_chol2inv(B)
#             counter_1 <- counter_1+1
#           }
#           counter_tot <- counter_tot+1
#         }
#       }
#     }
# 
#     # ---------------------update se2
#     # se2 <- 1
#     ytilde <- 0.5 * crossprod(y - zz_alpha - xzu %*% par)
#     if(prior_rand == "L"){
#       se2 <- extraDistr::rinvgamma(1, 1e-3 + 0.5*(no + sum(g_b) + sum(g_s)),
#                                    as.numeric(1e-3 + ytilde +
#                                                 0.5*sum(beta^2/omega2_beta) +
#                                                 0.5*sum( sdu^2/omega2_sdu ) ))
# 
#     }else{
#       se2 <- extraDistr::rinvgamma(1, 1e-3 + 0.5*(no + sum(g_b) + sum(g_s)),
#                                    as.numeric(1e-3 + ytilde +
#                                                 0.5*sum(beta^2/omega2_beta)/tau2[1] +
#                                                 0.5*sum( sdu^2/omega2_sdu )/tau2[2] ))
#     }
#     # }
# 
#     if(it > burnin & it %% thinin == 0){
#       # chain$mu[ii] <- mu
#       if(individual_effect) chain$alpha[ii, ] <- alpha
#       if(individual_effect) chain$sd_alpha[ii] <- sd_alpha
#       chain$beta[ii, ] <- beta
#       if(prior_fixed=="SS") chain$g_b[ii, ] <- g_b
#       if(prior_fixed=="SS") chain$pi_beta[ii] <- pi_beta
#       chain$sdu[ii, ] <- sdu
#       chain$sdu_tild[ii, ] <- sdu_tild
#       if(prior_rand=="SS") chain$g_s[ii, ] <- g_s
#       if(prior_rand=="SS") chain$pi_sdu[ii] <- pi_sdu
#       chain$se2[ii] <- se2
#       if(prior_rand=="DG") chain$a_s[ii] <- a_s
#       if(prior_rand=="DG") chain$kappa2[ii] <- kappa_2
#       # chain$k[ii] <- k_i
#       chain$omega2_beta[ii, ] <- omega2_beta
#       chain$omega2_sdu[ii, ] <- omega2_sdu
#       chain$theta[ii, ] <- theta[id_lower.tri]
#       chain$u[ii, ] <- c(U)
#       chain$tau2[ii, ] <- tau2
#       ii <- ii+1
#     }
#   }
# 
# 
#   if(correlation){
#     chain$acc_rate_theta <- counter_1/counter_tot
#     print("Acceptance ratio of theta:")
#     print(counter_1/counter_tot)
#   }
#   if(prior_rand=="DG"){
#     print("Acceptance ratio of a_s:")
#     print(counter_a_s/niter)
#     chain$acc_rate_a_s <- counter_a_s/niter
#   }
#   return(chain)
# }




#__________________________________________________________________________________________________________
#__________________________________________________________________________________________________________






#__________________________________________________________________________________________________________
#__________________________________________________________________________________________________________


#' MCMC algorithm
#'
#' @param y n-vector of the centered response variable
#' @param x n x q covariable matrix
#' @param ind n-vector indicating the level of each observation
#' @param prior_fixed character indicating which prior must be used to achieve selection of fixed effects. "SS" denote the Spike and Slab prior, "HS" denote the Horseshoe prior. The default value is "SS".
#' @param prior_rand character indicating which prior must be used to achieve selection of variance componantes. "SS" denote the Spike and Slab prior, "HS" denote the Horseshoe prior, "cauchy" denote the cauchy prior, "DG" denote the double-gamma prior, "L" denote the Laplace prior and "NIG" denote the Normal Inverse-Gamma prior. The default value is "SS".
#' @param correlation logical value indicating whether to learn the correlation matrix of the random effects. The default value is TRUE.
#' @param individual_effect logical value indicating whether to add a random individual effect. The default value is FALSE.
#' @param niter positive integer, indicating the number of MCMC iterations to perform, including the burn-in. Has to be larger than  nburn + 1. The default value is 10000.
#' @param burnin non-negative integer, indicating the number of iterations discarded as burn-in. Has to be smaller than to niter - 1. The default value is round(niter / 2).
#' @param thinin positive integer, indicating the degree of thinning to be performed. Every nthin draw is kept and returned. The default value is 1, implying that every draw is kept.
#' @param epsilon positive, real number. Determines the standard deviation of the proposal distribution for the Metropolis Hastings step for theta. The default value is 0.3.
#' @param epsilon_a_s positive, real number. Only available for prior_rand = "DG". Determines the standard deviation of the proposal distribution for the Metropolis Hastings step for a_xi. The default value is 1.
#' @param rho_epsilon positive, real number. Determines the standard deviation of the proposal distribution for the Metropolis Hastings step for rho. The default value is 0.05.
#'
#' @details For details concerning the algorithm please refer to the paper by ...
#' @return The value returned is a list object containing samples for parameters:
#'
#' \itemize{
#' \item  alpha, only if individual_effect==TRUE
#' \item  sd_alpha, only if individual_effect==TRUE
#' \item  beta
#' \item  g_b, only if prior_fixed=="SS"
#' \item  pi_beta, only if prior_fixed=="SS"
#' \item  sdu
#' \item  sdu_tild
#' \item  g_s, only if prior_fixed=="SS"
#' \item  pi_sdu, only if prior_fixed=="SS"
#' \item  se2
#' \item  a_s, only if prior_rand=="DG"
#' \item  kappa2, only if prior_rand=="DG"
#' \item  omega2_beta
#' \item  omega2_sdu
#' \item  theta
#' \item  u
#' \item  tau2
#' \item rho
#' }
#'  The list contains also:
#' \itemize{
#' \item  acc_rate_theta: acceptence rate of the MH step for sample theta, only if correlation == TRUE.
#' \item  acc_rate_a_s:  acceptence rate of the MH step for sample a_s, only if prior_rand=="DG".
#' }
#'
#' @note     To cite this package please reference:
#'
#' @export
#'
#' @examples
#'
# ris_AR <- function(y, X, ind, prior_fixed = "HS", prior_rand="HS", correlation = TRUE,
#                    individual_effect = FALSE, niter, burnin, thinin, epsilon=0.3,
#                    epsilon_a_s=1, rho_epsilon=0.05, K=0)
# {
#   if(!prior_fixed  %in% c("SS", "HS")) stop("prior_fixed choice must be 'SS'or 'HS'." )
#   if(!prior_rand  %in% c("SS", "HS", "DG", "cauchy", "NIG", "L")) stop("prior_rand choice must be 'SS', 'HS', 'DG', 'cauchy', 'L' or 'NIG'." )
#   if(burnin >= niter ) stop("'burnin' must be lower than 'niter'")
# 
#   print("Variance components selection in RIS models")
#   # if(prior_rand=="HS") print("Horseshoe prior for variance component in RIS models")
#   # if(prior_rand=="cauchy") print("Cauchy prior for variance component in RIS models")
#   # if(prior_rand=="NIG") print("Gaussian prior for variance component in RIS models")
# 
#   # Settings _________________________________________________________________
#   levels(ind) <- 1:length(unique(levels(ind)))
#   n_niv       <- length(levels(ind))
#   no          <- length(y)
#   ni          <- table(ind)
#   q_fix       <- ncol(X)
#   q_rand      <- q_fix + 1
# 
#   index <- unlist(lapply(1:(q_fix), function(l) 2*K+rep(q_rand - l, q_rand - l)))
# 
#   #sort by increasing factor
#   y       <- y[order(ind)]
#   X       <- X[order(ind), ]
#   ind <- sort(ind)
#   J       <- t(as(ind,  Class = "sparseMatrix"))
#   L       <- matrix(0, (q_rand)^2, q_rand)
#   pos     <- cbind(which(diag(q_rand)>0), 1:(q_rand))
#   L[pos]  <- 1
# 
#   if(individual_effect){
#     table_ind <- table(ind)
#     zz <- matrix(0, length(y), max(table_ind))
#     for(k in 1:n_niv) zz[ind==k, 1:table_ind[k]] <- diag(table_ind[k])
#   }
# 
#   # Outputs _________________________________________________________________
#   chain             <- list()
#   if(individual_effect) chain$alpha       <- matrix(0, floor(niter - burnin)/thinin, ncol(zz))
#   if(individual_effect) chain$sd_alpha    <- rep(NA, floor(niter - burnin)/thinin)
#   # chain$mu          <- rep(NA, floor(niter - burnin)/thinin)
#   chain$beta        <- matrix(0, floor(niter - burnin)/thinin, q_fix)
#   if(prior_fixed=="SS") chain$g_b         <- matrix(0, floor(niter - burnin)/thinin, q_fix)
#   if(prior_fixed=="SS") chain$pi_beta     <- rep(NA, floor(niter - burnin)/thinin)
#   chain$sdu         <- matrix(0, floor(niter - burnin)/thinin, q_rand)
#   chain$sdu_tild    <- matrix(0, floor(niter - burnin)/thinin, q_rand)
#   chain$g_s         <- matrix(0, floor(niter - burnin)/thinin, q_rand)
#   if(prior_rand=="SS") chain$pi_sdu         <- rep(NA, floor(niter - burnin)/thinin)
#   chain$tau2 <- matrix(NA, floor(niter - burnin)/thinin, 2) # by convention tau[, 1]=tau_beta et tau[, 2]=tau_sdu
#   chain$omega2_beta <- matrix(1, floor(niter - burnin)/thinin, q_fix)
#   chain$omega2_sdu  <- matrix(1, floor(niter - burnin)/thinin, q_rand)
#   chain$theta       <- matrix(0, floor(niter - burnin)/thinin, 0.5*q_rand*(q_rand-1))
#   chain$se2         <- rep(NA, floor(niter - burnin)/thinin)
#   chain$rho <- rep(NA, floor(niter - burnin)/thinin)
#   chain$u           <- matrix(NA, floor(niter - burnin)/thinin, n_niv*q_rand)
#   if(prior_rand=="DG") chain$a_s        <- rep(NA, floor(niter - burnin)/thinin)
#   if(prior_rand=="DG") chain$kappa2     <- rep(NA, floor(niter - burnin)/thinin)
#   # chain$k         <- rep(NA, floor(niter - burnin)/thinin)
# 
#   # Initialisation __________________________________________________________
#   pi_beta <- pi_sdu <- 0.5
#   mu <- 0 # mean(y)
#   if(individual_effect) alpha <- rep(0, ncol(zz))
#   if(individual_effect) sd_alpha <- 10000
#   U                 <- matrix(rnorm(q_rand*n_niv), q_rand, n_niv)
#   u                 <- c(U)
#   # u <- init$u
#   # U               <- matrix(u, q_rand, n_niv)
#   beta <- chain$beta[1, ] <- rnorm(q_fix)
#   sdu_tild <- rnorm(q_rand)
#   sdu <- abs(sdu_tild)
#   se2 <- chain$se2[1] <- 1 # rgamma(1, 2, 2)
# 
#   if(correlation) theta <- sample_thetaPourhamandi(q_rand) else theta <- matrix(0, q_rand, q_rand); theta[lower.tri(theta)] <- pi/2
#   chain$theta[1, ] <- theta[lower.tri(theta)]
#   if(correlation){
#     B <- corrPourhamandi(theta)
#     R0 <- tcrossprod(B)
#   }else{
#     B <- R0 <- diag(q_rand)
#   }
# 
# 
#   tau2 <- chain$tau2[1, ] <- c(1, 1)
#   if(prior_fixed=="SS"){
#     g_b <- chain$g_b[1, ] <- sample(c(0, 1), q_fix, replace = TRUE)
#     beta[g_b==0] <- 0
#     chain$beta[1, ] <- beta
#   }else{ # Horseshoe
#     g_b <- rep(1, q_fix)
#     tau2 <- chain$tau2[1, ] <- c(1e-5, 1)
#   }
#   xi2 <- rep(1, 2)
#   omega2_beta <- chain$omega2_beta[1, ] <- rep(1, q_fix)
#   nu2_beta <- rep(1, q_fix)
#   omega2_sdu <- chain$omega2_sdu[1, ] <- rep(1, q_rand)
#   nu2_sdu <- rep(1, q_rand)
#   if(prior_rand == "cauchy"){
#     nu2_sdu <- rep(2, q_rand)
#   }
#   if(prior_rand=="SS"){
#     g_s <- chain$g_s[1, ] <- sample(c(0, 1), q_rand, replace = TRUE)
#     sdu[g_s == 0] <- 0
#   }else{
#     g_s <- chain$g_s[1, ] <- rep(1, q_rand)
#   }
#   if(prior_rand == "DG"){
#     d1 <- 0.001
#     d2 <- 0.001
#     kappa_2 <- 20 #rgamma(1, shape = d1, rate = d1)
# 
#     b_w <- 10
#     a_s <- 0.1 # rexp(1, b_w)
#     l_p_omega_sdu_i <- sum(dgamma(omega2_sdu, a_s, a_s*kappa_2/2, log = TRUE))
#     l_p_a_s_i <- dexp(a_s, b_w, log=TRUE)
#   }
# 
#   #___________________________________________________________
# 
#   diag_niv <- Diagonal(n_niv)
#   diag_q_rand <- Diagonal(q_rand)
#   diag_q_rand_2 <- diag(q_rand)
#   tZ <- KhatriRao(t(J),  t(cbind(1, X)))
#   Z <- t(tZ)
#   zu <- crossprod(tZ, crossprod(kron_cpp(U, diag_q_rand_2), L))
#   zus <- zu %*% sdu_tild
#   xzu <- cbind(X, zu)
#   iR <- chol2inv(t(B))
#   id_lower.tri <- lower.tri(theta)
#   counter_tot <- counter_1 <- counter_a_s <- counter_rho <- 0
#   z_list <- list(); for(k in 1:n_niv) z_list[[k]] <- cbind(1, X[ind==k, ])
#   if(individual_effect){
#     zz_alpha <- zz %*% alpha
#     zztzz <- crossprod(zz)
#     svd_zztzz <- svd(zztzz)
#   }else{
#     zz_alpha <- rep(0, length(y))
#   }
#   par <- c(beta, sdu_tild)
# 
#   ii <- 1
# 
# 
#   rho <- runif(1)
#   # diag_0_Gamma =  arma::ones(n_niv, 1) * (1+rho*rho)
#   # diag_0_Gamma(0) = diag_0_Gamma(n_niv-1) = 1
#   # diag_1_Gamma = -rho* arma::ones(n_niv-1, 1)
#   # Gamma_inv = (arma::diagmat(diag_1_Gamma, -1) + arma::diagmat(diag_0_Gamma, 0) + arma::diagmat(diag_1_Gamma, 1))/(1-rho*rho);
# 
#   Gamma_0 <- matrix(rho, n_niv, n_niv); for(i in 1:n_niv) for(j in 1:n_niv) Gamma_0[i, j] <- rho^abs(i-j) # matrice intervenant dans la variance résiduelle
#   diag_ind <- diag(ni[1])
#   # Gamma.inv <- solve(kron_cpp(Gamma, diag_ind))
#   Gamma_0_inv <- solve(Gamma_0)
#   Gamma.inv <- kron_cpp(solve(Gamma_0), diag_ind)
#   # Gamma_k <- kron_cpp((Gamma_0[-k, -k]), diag_ind)
# 
#   x <- cbind(1, X[ind == 1, ])
# 
#   # MCMC ___________________________________________________________
# 
#   print("0 %")
#   for(it in 2:niter){
#     # print(it)
#     # if(it %% (niter/100) == 0) print(paste((it %/% (niter/100))*1, "%"))
#     if(it %% (niter/10) == 0) print(paste((it %/% (niter/10))*10, "%"))
#     # if(it == 100) browser()
# 
#     # --------------------- update mu
#     # mu <- rnorm(1, mean(y-zz_alpha - xzu %*% par), sd = sqrt(se2/no))
# 
#     # --------------------- update alpha (individual effect)
#     if(individual_effect){
#       iS_alpha <- diag(1/diag(crossprod(zz, Gamma.inv%*% zz)/se2 + diag(ncol(zz))/sd_alpha^2))
#       mu_alpha <- crossprod(iS_alpha, crossprod(zz, Gamma.inv %*% (y - X %*% beta - zu %*% sdu_tild))/se2)
#       alpha <- c(mvnfast::rmvn(1, mu_alpha, iS_alpha))
#       zz_alpha <- zz %*% alpha
#       sd_alpha <-  sqrt(1/rgamma(1, 1/2+ncol(zz)/2, rate=1/2 + crossprod(alpha)/2))
#     }
# 
#     # browser()
#     # --------------------- update beta
#     if(prior_fixed=="SS"){
#       for(k in 1:q_fix){
#         y_tild <- y - zz_alpha - X[, -k] %*% beta[-k] - zus
#         Sigma_bk <- 1/(crossprod(X[, k], Gamma.inv %*% X[, k])/se2 + 1/(se2 * tau2[1]))
#         log_R <- as.numeric(log(pi_beta)- log(1-pi_beta) - 0.5 *log(crossprod(X[, k], Gamma.inv %*% X[, k]) + tau2[1]) -0.5 * log(tau2[1]) + (0.5 * t(y_tild) %*% Gamma.inv %*% X[, k] %*% Sigma_bk %*% t(X[, k]) %*% Gamma.inv %*% y_tild / (se2)^2 ))
#         g_b[k] <- 1*(log(runif(1)) < log_R - sumlogs(log_R, 0))
# 
#         if(g_b[k]==0) {
#           beta[k] <- 0
#         }else{
#           beta[k] <- rnorm(1, as.numeric(Sigma_bk * crossprod(X[, k], Gamma.inv %*% y_tild))/se2, sqrt(as.numeric(Sigma_bk)))
#         }
#       }
#       pi_beta <- rbeta(1, 1+sum(g_b), 1+q_fix-sum(g_b));
#     }else{ # horseshoe
#       if(prior_rand == "L"){ tmp_V <- c(omega2_beta) }else{ tmp_V <- c(tau2[1]*omega2_beta) }
#       iV <- Diagonal(x = 1/tmp_V)
#       iS_part <- solve(crossprod(X, Gamma.inv %*% X) + iV)
#       mu_part <- crossprod(iS_part, crossprod(X, Gamma.inv %*% (y - zz_alpha - zus)))
#       beta <- as.vector(mvnfast::rmvn(1, mu_part, se2*iS_part))
#     }
#     xbeta <- X %*% beta
# 
#     # --------------------- update  sdu
# 
#     if(prior_rand=="SS"){
#       for(k in 1:q_rand){
#         y_tild <- y - zz_alpha - xbeta - zu[, -k] %*% sdu_tild[-k]
#         Sigma_sk <- 1/(crossprod(zu[, k], Gamma.inv %*% zu[, k])/se2 + 1/(se2 * tau2[2]))
#         log_R <- as.numeric(log(pi_sdu)- log(1-pi_sdu) - 0.5 *log(crossprod(zu[, k], Gamma.inv %*% zu[, k]) + tau2[2]) -0.5 * log(tau2[2]) + (0.5 * t(y_tild) %*% Gamma.inv %*% zu[, k] %*% Sigma_sk %*% t(zu[, k]) %*% Gamma.inv %*% y_tild / (se2)^2 ))
#         g_s[k] <- 1*(log(runif(1)) < log_R- sumlogs(log_R, 0))
#         if(g_s[k]==0) {
#           sdu[k] <- sdu_tild[k] <- 0
#         }else{
#           sdu_tild[k] <- rnorm(1, as.numeric(Sigma_sk *crossprod(zu[, k], Gamma.inv %*% y_tild))/se2, sqrt(as.numeric(Sigma_sk)))
#           sdu[k] <- abs(sdu_tild[k])
#         }
#       }
#       pi_sdu <- rbeta(1, 1+sum(g_s), 1+q_rand-sum(g_s));
#     }else{ # horseshoe
#       if(prior_rand == "L"){ tmp_V <- omega2_sdu }else{ tmp_V <- tau2[2]*omega2_sdu }
#       iV <- Diagonal(x = 1/tmp_V)
#       iS_part <- solve(crossprod(zu, Gamma.inv %*% zu) + iV)
#       mu_part <- crossprod(iS_part, crossprod(zu, Gamma.inv %*%(y - zz_alpha - xbeta)))
#       sdu_tild <- as.vector(mvnfast::rmvn(1, mu_part, se2*iS_part))
#       sdu <- abs(sdu_tild)
#     }
#     zus <- zu %*% sdu_tild
#     par <- c(beta, sdu_tild)
#     # S <- Diagonal(x=sdu_tild)
#     S <- diag(sdu_tild)
# 
#     # --------------------- update u
#     Lambda <- kronecker(diag(n_niv), S)
#     iR <- kronecker(diag(n_niv), chol2inv(t(B)))
#     iS_part <- solve(crossprod(Z %*% Lambda, Gamma.inv %*% Z %*% Lambda)/se2 + iR)
#     mu_part <- crossprod(iS_part, crossprod(Z %*% Lambda, Gamma.inv %*% (y - zz_alpha -X%*%beta)/se2))
#     u <- mvnfast::rmvn(1, mu_part, iS_part)
#     U <- matrix(u, q_rand,  n_niv)
# 
# 
# 
#     # iS_part <- solve( kronecker(Gamma_0_inv, crossprod(x %*% S %*% B))/se2 + diag(q_rand*n_niv))
#     # mu_part <- crossprod(iS_part,  kronecker(Gamma_0_inv, S %*% t(x)) %*% (y - zz_alpha -X%*%beta)/se2 )
#     # u <- mvnfast::rmvn(1, mu_part, iS_part)
#     # U <- matrix(u, q_rand,  n_niv)
# 
# 
#     # k=1
#     # # iR <- chol2inv(t(B))
#     # for(k in 1:niv){
#     #   tmp_0 <- z_list[[k]] * rep(sdu_tild, each = nrow(z_list[[k]] )) # element wise product
#     #
#     #   tmp_1 <- crossprod(tmp_0)/se2
#     #   iS_part <- solve(tmp_1 + iR)
#     #   # mu_part <- crossprod(iS_part, crossprod(tmp_0, (y[ind==k]  - zz_alpha[ind==k] -X[ind==k,]%*%beta)/se2))
#     #   #  mu_part <- iS_part %*% t(tmp_0) %*% (y[ind==k] - zz_alpha[ind==k] -X[ind==k,]%*%beta)/se2
#     #   if(individual_effect){ # to speed up
#     #     mu_part <- my_mu_part(iS_part, tmp_0, (y[ind==k] - zz_alpha[ind==k] -X[ind==k,]%*%beta)/se2)
#     #   }else{
#     #     mu_part <- my_mu_part(iS_part, tmp_0, (y[ind==k] -X[ind==k,]%*%beta)/se2)
#     #   }
#     #
#     #   u_tmp <- mvnfast::rmvn(1, mu_part, iS_part)
#     #   U[, k] <- c(u_tmp)
#     # }
#     zu <- crossprod(tZ, crossprod(kron_cpp(U, diag_q_rand_2), L))
#     xzu <- cbind(X, zu)
# 
# 
# 
#     # --------------------- update tau2
#     # fixed part
#     if(prior_fixed=="ss"){
#       tau2 <- pmax(1e-6, extraDistr::rinvgamma(1, sum(g_b)/2 + 1/2, 1/2 + 0.5*sum(beta^2)/se2))
#     }else{
#       tau2[1] <- max(1e-6, extraDistr::rinvgamma(1, 0.5*(q_fix + 1),  1/xi2[1] + 0.5*sum(beta^2/omega2_beta)/se2))
#       xi2[1] <- extraDistr::rinvgamma(1, 1, 1 + 1/tau2[1])
#     }
# 
#     # random part
#     if(prior_rand == "SS"){
#       tau2[2] <- pmax(1e-6, extraDistr::rinvgamma(1, sum(g_s)/2 + 1/2, 1/2 + 0.5*sum(sdu^2)/se2) )
#     }
# 
#     if(prior_rand == "L"){
#       tau2[2] <- extraDistr::rdgamma(1, shape = 1/2 + q_rand, scale = 1/2 + sum(omega2_sdu)/2)
#     }
# 
#     if(prior_rand == "HS"){
#       tau2[2] <- pmax(1e-6, extraDistr::rinvgamma(1, 0.5*(q_rand + 1), 1/xi2[2] + 0.5*sum(sdu^2 /omega2_sdu )/se2) )
#       xi2[2] <- extraDistr::rinvgamma(1, 1, 1 + 1/tau2[2])
#     }
# 
#     if(prior_rand == "NIG"){
#       tau2[2] <- pmax(1e-6, extraDistr::rinvgamma(1, 0.5*(q_rand) + 1e-2, 1e-2 + 0.5*sum(sdu^2 /omega2_sdu )/se2) )
#     }
# 
#     if(prior_rand == "cauchy"){}
# 
#     if(prior_rand == "DG"){}
# 
#     # --------------------- omega2 for beta and s
#     # fixed part
#     if(prior_fixed=="HS"){
#       omega2_beta <- pmax(1e-3, extraDistr::rinvgamma(q_fix, 1, 1/nu2_beta + beta^2/(2*se2*tau2[1])))
#       nu2_beta <- extraDistr::rinvgamma(q_fix, 1, 1 + 1/omega2_beta)
#     }
# 
#     # random part
#     if(prior_rand == "HS"){
#       omega2_sdu <- pmax(1e-3, extraDistr::rinvgamma(q_rand, 1, 1/nu2_sdu + sdu^2/(2*se2*tau2[2])))
#       nu2_sdu <- extraDistr::rinvgamma(q_rand, 1, 1 + 1/omega2_sdu)
#     }
# 
#     if(prior_rand == "cauchy"){
#       omega2_sdu <- pmax(1e-3, extraDistr::rinvgamma(q_rand, 1, 1/nu2_sdu + sdu^2/(2*se2*tau2[2])))
#       # nu2_sdu <- extraDistr::rinvgamma(q_rand, 1, 1 + 1/omega2_sdu)
#     }
# 
#     if(prior_rand == "L"){
#       omega2_sdu  <- statmod::rinvgauss(q_rand, mean = tau2[2] * sqrt(se2)/abs(sdu), shape = tau2[2])
#     }
# 
#     if(prior_rand == "NIG"){}
# 
#     if(prior_rand == "DG"){
#       omega2_sdu <- GIGrvg::rgig(q_rand, a_s-1/2, a_s*kappa_2, sdu^2)
#     }
# 
#     # --------------------- a_s and kappa2 for double gamma prior
# 
#     if(prior_rand == "DG"){
#       # # MH step for a_s ______________
#       before <- a_s
#       # from shrinkTVP package Bitto & Frühwirth-Schnatter 2018
#       a_s <- MH_step(a_s, epsilon_a_s, q_rand, kappa_2, sdu^2, b_w, nu=5, d1, d2)
#       if(before != a_s) counter_a_s <- counter_a_s+1
# 
#       #_________________________
# 
#       kappa_2 <- rgamma(1, shape = d1+q_rand*a_s, rate = d2 + a_s/2*sum(omega2_sdu))
#     }
# 
# 
#     # --------------------- correlation R0 = B B'     B= f(theta)
#     # if(correlation){
#     #   if(it == 2){
#     #     # l_p_theta_i     <- sum(index * log(sin(theta[id_lower.tri] )))
#     #     # l_p_Y_theta_i    <- sum(mvnfast::dmvn(t(U), rep(0, q_rand), t(B), log=TRUE, isChol=TRUE))
#     #     l_p_Y_theta_i <- l_p_Y_theta_star_cpp_2(n_niv, q_rand, U, theta[lower.tri(theta)], B, K)
#     #   }
#     #   for(i in sample(2:q_rand, q_rand-1)){
#     #     for(j in sample(1:(i-1), i-1)){
#     #       # print(c(i, j))
#     #       theta_star <- theta
#     #       theta_star[i, j] <- runif(1, max(0, theta_star[i, j]-epsilon), min(pi, theta_star[i, j]+epsilon))
#     #       l_p_star_old <- l_p_old_star <- 0
#     #       # theta_star[i, j] <- truncnorm::rtruncnorm(1, a=0, b=pi, mean = theta[i, j], sd=epsilon)
#     #       # l_p_star_old <- log(truncnorm::dtruncnorm(theta_star[i, j], a = 0, b=pi, mean = theta[i, j], sd=epsilon))
#     #       # l_p_old_star <- log(truncnorm::dtruncnorm(theta[i, j], a = 0, b=pi, mean = theta_star[i, j], sd=epsilon))
#     #
#     #       # l_p_theta_star  <- sum(index * log(sin(theta_star[id_lower.tri])))
#     #       B_star  <- corrPourhamandi_new(B, theta_star, i, j)
#     #       # l_p_Y_theta_star <- sum(mvnfast::dmvn(t(U), rep(0, q_rand), t(B_star), log=TRUE, isChol=TRUE))
#     #
#     #       # iR_star <- chol2inv(t(B_star))
#     #       # iR_star <- my_chol2inv(B_star)
#     #       l_p_Y_theta_star <- l_p_Y_theta_star_cpp_2(n_niv, q_rand, U, theta_star[lower.tri(theta_star)], B_star, K)
#     #
#     #       #________________________________
#     #
#     #       # log_ratio <- min(0, l_p_Y_theta_star + l_p_theta_star - l_p_Y_theta_i - l_p_theta_i + l_p_old_star - l_p_star_old)
#     #       log_ratio <- min(0, l_p_Y_theta_star  - l_p_Y_theta_i  ) # + l_p_theta_star - l_p_theta_i + l_p_old_star - l_p_star_old)
#     #
#     #       #browser()
#     #       if(runif(1) < exp(log_ratio)){
#     #         theta[i, j] <-  theta_star[i, j]
#     #         l_p_Y_theta_i <- l_p_Y_theta_star
#     #         # l_p_theta_i <- l_p_theta_star
#     #         B <- B_star
#     #         iR <- my_chol2inv(B)
#     #         counter_1 <- counter_1+1
#     #       }
#     #       counter_tot <- counter_tot+1
#     #     }
#     #   }
#     # }
# 
# 
#     # --------------------- correlation R0 = B B'     B= f(theta)
#     if(correlation){
#       if(it==2){
#         l_p_theta_i     <- sum(index * log(sin(theta[id_lower.tri] )))
#         l_p_Y_theta_i    <- sum(mvnfast::dmvn(t(U), rep(0, q_rand), t(B), log=TRUE, isChol=TRUE))
#         # l_p_Y_theta_i <- l_p_Y_theta_star_cpp_2(n_niv, q_rand, U, theta[lower.tri(theta)], B, 0)
#       }
#       for(i in sample(2:q_rand, q_rand-1)){
#         for(j in sample(1:(i-1), i-1)){
#           # print(c(i, j))
#           theta_star <- theta
#           # theta_star[i, j] <- runif(1, max(0, theta_star[i, j]-epsilon), min(pi, theta_star[i, j]+epsilon))
#           # l_p_star_old <- l_p_old_star <- 0
#           theta_star[i, j] <- truncnorm::rtruncnorm(1, a=0, b=pi, mean = theta[i, j], sd=epsilon)
#           l_p_star_old <- log(truncnorm::dtruncnorm(theta_star[i, j], a = 0, b=pi, mean = theta[i, j], sd=epsilon))
#           l_p_old_star <- log(truncnorm::dtruncnorm(theta[i, j], a = 0, b=pi, mean = theta_star[i, j], sd=epsilon))
# 
#           l_p_theta_star  <- sum(index * log(sin(theta_star[id_lower.tri])))
#           B_star  <- corrPourhamandi_new(B, theta_star, i, j)
#           l_p_Y_theta_star <- sum(mvnfast::dmvn(t(U), rep(0, q_rand), t(B_star), log=TRUE, isChol=TRUE))
# 
#           # iR_star <- chol2inv(t(B_star))
#           # iR_star <- my_chol2inv(B_star)
#           # l_p_Y_theta_star <- l_p_Y_theta_star_cpp_2(n_niv, q_rand, U, theta_star[lower.tri(theta_star)], B_star, 0)
# 
#           #________________________________
# 
#           log_ratio <- min(0, l_p_Y_theta_star + l_p_theta_star - l_p_Y_theta_i - l_p_theta_i + l_p_old_star - l_p_star_old)
#           # log_ratio <- min(0, l_p_Y_theta_star - l_p_Y_theta_i) #  + l_p_theta_star - l_p_theta_i + l_p_old_star - l_p_star_old)
# 
#           #browser()
#           if(runif(1) < exp(log_ratio)){
#             theta[i, j] <-  theta_star[i, j]
#             l_p_Y_theta_i <- l_p_Y_theta_star
#             l_p_theta_i <- l_p_theta_star
#             B <- B_star
#             iR <- my_chol2inv(B)
#             counter_1 <- counter_1+1
#           }
#           counter_tot <- counter_tot+1
#         }
#       }
#     }
# 
# 
#     # ----------------------- update rho
#     # rho <- 1
#     # diag_0_Gamma =  arma::ones(T, 1) * (1+rho*rho)
#     # diag_0_Gamma(0) = diag_0_Gamma(T-1) = 1
#     # diag_1_Gamma = -rho* arma::ones(T-1, 1)
#     # Gamma_inv = (arma::diagmat(diag_1_Gamma, -1) + arma::diagmat(diag_0_Gamma, 0) + arma::diagmat(diag_1_Gamma, 1))/(1-rho*rho);
# 
# 
#     # Metropolis-Hastings
#     tp <- (y - zz_alpha - xzu %*% par)
#     tmp_MH_rho <- crossprod(tp, crossprod(Gamma.inv, tp))
#     det_gamma <- determinant(Gamma_0, logarithm = TRUE)$modulus
# 
#     rho_star <- runif(1, max(-1, rho-rho_epsilon),min(1, rho+rho_epsilon))
#     log_Qnew = dunif(rho_star,max(-1, rho-rho_epsilon), min(1,rho+rho_epsilon), log = T) # proposal distribution Q
#     log_Qold = dunif(rho,max(-1, rho_star-rho_epsilon), min(1,rho_star+rho_epsilon), log = T)
# 
#     Gamma_star <- matrix(rho_star, n_niv, n_niv); for(i in 1:n_niv) for(j in 1:n_niv) Gamma_star[i, j] <- rho_star^abs(i-j)
#     Gamma_star.inv <- kron_cpp(solve(Gamma_star), diag_ind)
#     tmp_star <- crossprod(tp, Gamma_star.inv) %*% tp
# 
#     det_gamma_star <- determinant(Gamma_star, logarithm = TRUE)$modulus
# 
#     log_ratio <- -ni[1]/2 *  det_gamma_star -
#       1/(2*se2) *  tmp_star + log_Qold +
#       ni[1]/2 * det_gamma +
#       1/(2*se2)*  tmp_MH_rho - log_Qnew
# 
# 
#     if(log(runif(1)) < (as.numeric(log_ratio))) {
#       tmp_MH_rho <- tmp_star
#       det_gamma <- det_gamma_star
#       rho <- rho_star
#       # Gamma <- Gamma_star
#       Gamma_0 <- Gamma_star
#       Gamma_0_inv <- solve(Gamma_0)
#       Gamma.inv <- Gamma_star.inv
#       counter_rho <- counter_rho+1
#     }
#     # print(rho)
# 
# 
# 
# 
#     # ---------------------update se2
#     # se2 <- 1
#     ytilde <- 0.5 * crossprod(y - zz_alpha - X %*% beta - zu %*% sdu_tild, Gamma.inv %*% (y - zz_alpha - X %*% beta - zu %*% sdu_tild))
#     # if(individual_effect){
#     #   if(prior_rand == "L"){
#     #     se2 <- extraDistr::rinvgamma(1, 1e-3 + 0.5*(no + sum(g_b) + sum(g_s)+1),
#     #                                  as.numeric(1e-3 + 0.5*crossprod(alpha)/sd_alpha^2 + ytilde +
#     #                                               0.5*sum(beta^2/omega2_beta) +
#     #                                               0.5*sum( sdu^2/omega2_sdu ) ))
#     #
#     #   }else{
#     #     se2 <- extraDistr::rinvgamma(1, 1e-3 + 0.5*(no + sum(g_b) + sum(g_s)+1),
#     #                                  as.numeric(1e-3 + 0.5*crossprod(alpha)/sd_alpha^2 + ytilde +
#     #                                               0.5*sum(beta^2/omega2_beta)/tau2[1] +
#     #                                               0.5*sum( sdu^2/omega2_sdu )/tau2[2] ))
#     #   }
#     # }else{
#     if(prior_rand == "L"){
#       se2 <- extraDistr::rinvgamma(1, 1e-3 + 0.5*(no + sum(g_b) + sum(g_s)),
#                                    as.numeric(1e-3 + ytilde +
#                                                 0.5*sum(beta^2/omega2_beta) +
#                                                 0.5*sum( sdu^2/omega2_sdu ) ))
# 
#     }else{
#       se2 <- extraDistr::rinvgamma(1, 1e-3 + 0.5*(no + sum(g_b) + sum(g_s)),
#                                    as.numeric(1e-3 + ytilde +
#                                                 0.5*sum(beta^2/omega2_beta)/tau2[1] +
#                                                 0.5*sum( sdu^2/omega2_sdu )/tau2[2] ))
#     }
#     # }
# 
#     if(it > burnin & it %% thinin == 0){
#       # chain$mu[ii] <- mu
#       if(individual_effect) chain$alpha[ii, ] <- alpha
#       if(individual_effect) chain$sd_alpha[ii] <- sd_alpha
#       chain$beta[ii, ] <- beta
#       if(prior_fixed=="SS") chain$g_b[ii, ] <- g_b
#       if(prior_fixed=="SS") chain$pi_beta[ii] <- pi_beta
#       chain$sdu[ii, ] <- sdu
#       chain$sdu_tild[ii, ] <- sdu_tild
#       if(prior_rand=="SS")  chain$g_s[ii, ] <- g_s
#       if(prior_rand=="SS")  chain$pi_sdu[ii] <- pi_sdu
#       chain$se2[ii] <- se2
#       chain$rho[ii] <- rho
#       if(prior_rand=="DG")  chain$a_s[ii] <- a_s
#       if(prior_rand=="DG")  chain$kappa2[ii] <- kappa_2
#       chain$omega2_beta[ii, ] <- omega2_beta
#       chain$omega2_sdu[ii, ] <- omega2_sdu
#       chain$theta[ii, ] <- theta[id_lower.tri]
#       chain$u[ii, ] <- c(U)
#       chain$tau2[ii, ] <- tau2
#       ii <- ii+1
#     }
#   }
#   if(correlation){
#     chain$acc_rate_theta <- counter_1/counter_tot
#     print("Acceptance ratio of theta:")
#     print(counter_1/counter_tot)
#   }
#   if(prior_rand=="DG"){
#     print("Acceptance ratio of a_s:")
#     print(counter_a_s/niter)
#     chain$acc_rate_a_s <- counter_a_s/niter
#   }
#   print("Acceptance ratio of rho:")
#   print(counter_rho/niter)
#   chain$acc_rate_rho <- counter_rho/niter
#   return(chain)
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
