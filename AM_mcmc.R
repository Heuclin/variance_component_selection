



# Rcpp::sourceCpp('fct.cpp')

mode <- function(x) density(x)$x[which.max(density(x)$y)]

sdumlogs <- function(m){
  M <- max(m)
  return(M+log(sum(exp(m-M))))
}

proba_jointe <- function(gamma, idx, n = 3){
  proba.g.jointe <- sort(table(apply(1*gamma, 1, function(l) paste(l, collapse=""))) / nrow(gamma), decreasing = T ) [1:n]
  for(i in 1:length(proba.g.jointe)){
    st <- names(proba.g.jointe)[i]
    names(proba.g.jointe)[i] <- paste(rownames(idx)[str_locate_all(st, "1")[[1]][, 1]], collapse = ", ")
  }
  return(proba.g.jointe)
}

animal_model <- function(y, svdA, prior = "HS", niter = 5000, burnin = 2000, thinin = 5, method = "folded", ...){
  
  if(prior != "NG" & prior != "HS" & prior != "C" & prior !="SS" & prior !="NIG") stop("Prior choice must be one of 'HS', 'NIG', 'C', 'NG' or 'SS'" )
  
  if(prior %in% c("HS", "C", "NIG") ){
    time <- system.time(chain <- am(y = y, svdA = svdA, prior = prior, niter = niter, burnin = burnin, thinin = thinin, method = method, ...))
  }
  
  if(prior == "NG"){
    time <- system.time(chain <- am.ng(y = y, svdA = svdA, niter = niter, burnin = burnin, thinin = thinin, ...))
  }
  
  if(prior == "SS"){
    time <- system.time(chain <- am.ss(Y=y, MAT = svdA, niter = niter, burnin = burnin, thinin = thinin, method = method, ...))
  }
  
  chain$time <- time
  return(chain)
}

# Normal-gamma ------------------------------------------------------------------

#' @title Half-Cauchy prior for variance component in Animal models
#' @author  M. Denis, B. Heuclin and F. Mortier
am.ng <- function(y, svdA, s0=1e-3, r0=1e-3, stau = 0.5, rtau = 0.5, niter, burnin = 0, thinin = 1, epsilon_a_s=1){
  print("Normal Gamma prior for variance component in Animal models")
  
  require(extraDistr)
  #require(GIGrvg)
  no <- length(y)
  q <- length(svdA)
  mu.out <- se2.out <- zeta.out <- rep(0, niter)
  tau2.out <- rep(1, niter)
  su.out <- xi.out <-  omega.out <- nu.out <- matrix(0,niter,q)
  mu.out[1] <- mean(y)
  u <- matrix(rnorm(q*no),no)
  mu.y <- mean(y)
  se2.out[1] <- var(y-mu.out[1]-rowMeans(u))
  # tau2.out[1] <- 1
  xi.out[1, ] <- (rnorm(q))
  su.out[1, ] <- abs(xi.out[1, ])
  
  omega.out[1, ] <- 1
  logLik <- matrix(0, niter, no)
  zeta.out <- 2
  nu.out <- rep(2, q)
  
  d1 <- 0.001
  d2 <- 0.001
  kappa_2 <- 20 #rgamma(1, shape = d1, rate = d1)
  b_w <- 10
  a_s <- 0.1 # rexp(1, b_w)
  counter_a_s <- 0
  
  print("0 %")
  for(it in 2:niter){
    if(it %% (niter/10) == 0) print(paste((it %/% (niter/10))*10, "%"))
    
    # update mu
    ytilde <- y-u%*%su.out[it-1,]
    mu.out[it] <- rnorm(1,mean(ytilde),sqrt(se2.out[it-1]/no))
    
    su.tmp <- su.out[it-1,]
    xi.tmp <- xi.out[it-1,]
    
    for(j in 1:q){
      ytilde <- y-mu.out[it]-u[,-j,drop=FALSE]%*%su.tmp[-j]
      
      # update su
      sSDj <- 1/(1/(omega.out[it-1,j]*tau2.out[it-1])+crossprod(u[,j]))
      mSUj <- crossprod(sign(xi.tmp[j]) * u[,j], ytilde)*sSDj
      xi.tmp[j] <- rnorm(1,mean = mSUj, sd = sqrt(se2.out[it-1]*sSDj))
      su.tmp[j] <- abs(xi.tmp[j])
      
      # update u 
      sDj <- 1/(su.tmp[j]^2/se2.out[it-1]+1/svdA[[j]]$d)
      mUj <- su.tmp[j]/se2.out[it-1]*sDj*crossprod(svdA[[j]]$u, ytilde)
      u[,j] <- svdA[[j]]$u%*%rnorm(no, mUj, sqrt(sDj))
    }
    
    su.out[it,] <- su.tmp
    xi.out[it,] <- xi.tmp
    ytilde <- y-mu.out[it]- u%*%su.out[it,]
    
    #_________________________
    # print(max(su.out[it,]))
    
    omega.out[it, ] <- GIGrvg::rgig(q, a_s-1/2, a_s*kappa_2, su.out[it,]^2)
    
    # # MH step for a_s ______________
    before <- a_s
    # from shrinkTVP package Bitto & Frühwirth-Schnatter 2018
    a_s <- MH_step(a_s, epsilon_a_s, q, kappa_2, su.out[it,]^2, b_w, nu=5, d1, d2)
    if(before != a_s) counter_a_s <- counter_a_s+1
    
    kappa_2 <- rgamma(1, shape = d1+q*a_s, rate = d2 + a_s/2*sum(omega.out[it,]))
    
    # print(round(c(a_s, kappa_2), 2))
    # print(tau2.out[it])
    #_________________________
    
    #se2
    se2.out[it] <- extraDistr::rinvgamma(1, no/2+q/2+s0, 0.5*sum(ytilde^2) + (0.5/tau2.out[it-1])*sum(su.out[it,]^2/omega.out[it-1,])+r0)
    
    # # update omegaj
    # omega.out[it,] <- extraDistr::rinvgamma(q, 1, su.out[it,]^2/(2*se2.out[it]*tau2.out[it-1])+ 1/nu.out)
    # nu.out <- extraDistr::rinvgamma(q, 1, 1 + 1/omega.out[it, ])
    
    # # update tau2
    # tau2.out[it] <- extraDistr::rinvgamma(1, (1+q)/2, 0.5/se2.out[it]*sum(su.out[it,]^2/omega.out[it,]) + 1/zeta.out)
    # zeta.out <- extraDistr::rinvgamma(1, 1, 1 + 1/tau2.out[it])
    
    logLik[it,] <- dnorm(y,mu.out[it]+u%*%su.out[it,],sqrt(se2.out[it]),log = TRUE)
  }
  iter <- seq(burnin+1, niter, thinin)
  chain <- list(sdu=su.out[iter, ], 
                se2=se2.out[iter], 
                mu = mu.out[iter], 
                tau2 = tau2.out[iter],
                omega = omega.out[iter, ], 
                s=xi.out[iter, ], 
                logLik=logLik[iter, ])
  
  print("Acceptance ratio of a_s:")
  print(counter_a_s/niter)
  chain$acc_rate_a_s <- counter_a_s/niter
  
  return(chain)
}


# Horseshoe ------------------------------------------------------------------

#' @title Half-Cauchy prior for variance component in Animal models
#' @author  M. Denis, B. Heuclin and F. Mortier
am <- function(y, svdA, prior = "HS", s0=1e-3, r0=1e-3, stau = 0.5, rtau = 0.5, niter, burnin = 0, thinin = 1, method = "folded"){
  if(prior == "HS") print("Horseshoe prior for variance component in Animal models")
  if(prior == "C") print("Cauchy prior for variance component in Animal models")
  if(prior == "NG") print("Normal gamma prior for variance component in Animal models")
  if(prior == "NIG") print("Normal inverse-gamma prior for variance component in Animal models")
  
  require(extraDistr)
  #require(GIGrvg)
  no <- length(y)
  q <- length(svdA)
  
  chain <- list()
  chain$mu <- chain$tau2 <- chain$se2 <- rep(NA, floor(niter - burnin)/thinin)
  chain$sdu <- chain$s <- chain$omega <- matrix(NA, floor(niter - burnin)/thinin, q)
  chain$logLik <- matrix(NA, floor(niter - burnin)/thinin, no)
  
  mu <-  mean(y)
  # if(is.null(sim$U)) u <- matrix(rnorm(q*no),no) else u <- sim$U
  u <- matrix(rnorm(q*no),no)
  se2 <-  c(var(y-mu-rowMeans(u)))
  if(method == "folded") xi <- rnorm(q)
  if(method == "truncated") xi <- abs(rnorm(q))
  sdu <- abs(xi)
  tau2 <- 100
  omega <- rinvgamma(q, 1, 1)
  nu <- rep(2, q)
  zeta <- 2
  
  if(prior == "NG"){
    d1 <- 0.001
    d2 <- 0.001
    kappa_2 <- 20 #rgamma(1, shape = d1, rate = d1)
    b_w <- 10
    a_s <- 0.1 # rexp(1, b_w)
    counter_a_s <- 0
    epsilon_a_s <- 1
    tau2 <- 1
  }
  
  # if(prior == "SS"){
  #   g <- sample(c(TRUE, FALSE), q, replace = TRUE)
  #   xi <- rep(0, q)
  #   xi[g] <- rnorm(sum(g))
  #   sdu <- abs(xi)
  #   
  #   chain$g <- matrix(NA, floor(niter - burnin)/thinin, q)
  #   chain$pi <- rep(NA, floor(niter - burnin)/thinin)
  # }
  
  print("0 %")
  ii <- 1
  for(it in 1:niter){
    if(it %% (niter/10) == 0) print(paste((it %/% (niter/10))*10, "%"))
    
    # update mu
    ytilde <- y - u %*% sdu
    mu <- rnorm(1, mean(ytilde), sqrt(se2/no))
    # mu <- sim$mu
    
    for(j in sample(1:q)){
      ytilde <- y - mu - u[, -j, drop = FALSE] %*% sdu[-j]
      
      # update sdu
      sSDj <- 1/( 1/(omega[j]*tau2) + crossprod(u[,j]) )
      mSUj <- crossprod(sign(xi[j]) * u[,j], ytilde)*sSDj
      if(method == "folded"){
        xi[j] <- rnorm(1, mean = mSUj, sd = sqrt(se2*sSDj))
      }
      if(method == "truncated"){
        xi[j] <-  extraDistr::rtnorm(1, mean = mSUj , sd = sqrt(se2*sSDj), a = 0, b = Inf)
      }
      sdu[j] <- abs(xi[j])
      
      # update u 
      sDj <- 1/(sdu[j]^2/se2+1/svdA[[j]]$d)
      mUj <- sdu[j]/se2*sDj*crossprod(svdA[[j]]$u, ytilde)
      u[,j] <- svdA[[j]]$u%*%rnorm(no, mUj, sqrt(sDj))
      
      # v <- diag( 1 / (sdu[j]^2/se2 + 1/ svdA[[j]]$d) )
      # u[, j] <- tcrossprod(svdA[[j]]$u , rmvn(1, sdu[j]/se2* crossprod(v, crossprod(svdA[[j]]$u, ytilde)), sqrt(v), isChol = TRUE))
    }
    
    
    if(prior=="HS"){
      # update omegaj
      omega <- extraDistr::rinvgamma(q, 1, sdu^2/(2*se2*tau2) + 1/nu)
      nu <- extraDistr::rinvgamma(q, 1, 1 + 1/omega)
      
      # update tau2
      tau2 <- extraDistr::rinvgamma(1, (1+q)/2, 1/2/se2*sum(sdu^2/omega) + 1/zeta)
      zeta <- extraDistr::rinvgamma(1, 1, 1 + 1/tau2)
      # tau2 <- 1e-5
    }
    
    if(prior == "C"){
      omega <- extraDistr::rinvgamma(q, 1, sdu^2/(2*se2*tau2) + 0.5)
      tau2 <- 1
    }
    
    
    
    
    if(prior=="NG"){
      omega <- GIGrvg::rgig(q, lambda = a_s-1/2, chi = a_s*kappa_2, psi = sdu^2)
      # # MH step for a_s ______________
      before <- a_s
      # from shrinkTVP package Bitto & Frühwirth-Schnatter 2018
      a_s <- MH_step(a_s, epsilon_a_s, q, kappa_2, sdu^2, b_w, nu=5, d1, d2)
      if(before != a_s) counter_a_s <- counter_a_s+1
      kappa_2 <- rgamma(1, shape = d1+q*a_s, rate = d2 + a_s/2*sum(omega))
      tau2 <- 1
      print(round(c(a_s, kappa_2), 2))
      
    }
    
    
    # #_________________________
    # omega.out[it, ] <- GIGrvg::rgig(q, a_s-1/2, a_s*kappa_2, su.out[it,]^2)
    # # # MH step for a_s ______________
    # before <- a_s
    # # from shrinkTVP package Bitto & Frühwirth-Schnatter 2018
    # a_s <- MH_step(a_s, epsilon_a_s, q, kappa_2, su.out[it,]^2, b_w, nu=5, d1, d2)
    # if(before != a_s) counter_a_s <- counter_a_s+1
    # 
    # kappa_2 <- rgamma(1, shape = d1+q*a_s, rate = d2 + a_s/2*sum(omega.out[it,]))
    # #_________________________
    
    
    
    if(prior == "NIG"){
      omega <- rep(1, q)
      tau2 <- pmax(1e-6, extraDistr::rinvgamma(1, 0.5*q + 1e-2, 1e-2 + 0.5*sum(sdu^2 /omega)/se2) )
    }
    
    # if(prior == "SS"){
    #   omega <- rep(1, q)
    #   tau2 <- pmax(1e-6, extraDistr::rinvgamma(1, sum(g)/2 + 1/2, 1/2 + 0.5*sum(sdu^2)/se2) )
    # }
    
    #se2
    ytilde <- y - mu - u %*% sdu
    se2 <- extraDistr::rinvgamma(1, no/2 + q/2 + s0, 0.5*sum(ytilde^2) + (0.5/tau2)*sum(sdu^2/omega) + r0)
    # se2 <- sim$se2
    
    
    if(it > burnin & it %% thinin == 0){
      chain$mu[ii] <- mu
      chain$se2[ii] <- se2
      chain$sdu[ii, ] <- sdu
      chain$s[ii, ] <- xi
      chain$tau2[ii] <- tau2
      chain$omega[ii, ] <- omega
      chain$logLik[ii, ] <- dnorm(y, mu + u%*%sdu , sqrt(se2), log = TRUE)
      # if(prior=="SS"){
      #   chain$g[ii, ] <- g
      #   chain$pi[ii] <- pi_sdu
      # }
      ii <- ii+1
    }
  }
  
  return(chain)
}




# Spike and Slab ----------------------------------------------------------


my.sample <- function(a)  ifelse(length(a) == 1, a, sample (a, 1))

am.ss <- function(Y, MAT, a=1, b=1, niter = 5000, burnin = 0, thinin = 1, method = "folded"){   # Dirac sdur le spike, slab independant
  print("Spike and slab prior for variance component in Animal models")
  
  require(mvnfast)
  # require(msm)
  require(statmod)
  require(stringr)
  # require(MCMCglmm)
  
  no <- length(Y); q <- length(MAT)
  # Valeurs initiales 
  pi_sdu <- 0.5
  g <- sample(x = c(TRUE, FALSE), size = q, replace = T, prob = c(pi_sdu, 1-pi_sdu)) # init$g				# gamma
  mu <- mean(Y)	
  # if(is.null(sim$U)) U <- matrix(rnorm(q*no),no) else U <- sim$U
  U <- matrix(NA, no, q)
  for(i in 1:q){
    U[, i] <- tcrossprod(MAT[[i]]$u, rmvn(1, rep(0, no), diag(sqrt(MAT[[i]]$d)), isChol = T)   )
  }
  xi <- rep(0, q)
  xi[g] <- rnorm(sum(g))
  if(method == "folded") xi[g] <- rnorm(sum(g))
  if(method == "truncated") xi[g] <- abs(rnorm(sum(g)))
  sdu <- abs(xi)
  se2 <- extraDistr::rinvgamma(1,  a,  b)		# sigma^2
  tau2 <- extraDistr::rinvgamma(1,  a,  b) # 100
  omega <- rep(1, q)
  
  # declaration des sorties
  chain <- list()
  chain$mu <- chain$tau2 <- chain$se2 <- rep(NA, floor(niter - burnin)/thinin)
  chain$sdu <- chain$s <- chain$omega <- matrix(NA, floor(niter - burnin)/thinin, q)
  chain$logLik <- matrix(NA, floor(niter - burnin)/thinin, no)
  chain$g <- matrix(NA, floor(niter - burnin)/thinin, q)
  chain$pi <- rep(NA, floor(niter - burnin)/thinin)
  
  
  print("0 %")
  ii <- 0
  for (it in 1:niter){
    if(it %% (niter/10) == 0) print(paste((it %/% (niter/10))*10, "%"))
    
    # 1. Mise a jour de mu:
    ytilde <- Y - U %*% sdu
    mu <- rnorm(1, mean(ytilde), sqrt(se2/no))
    Ymu <- Y - mu
    
    # 2. Mise a jour de gamma et sdu avec integration univariee sdur  sdu_i
    for(i in sample(1:q)){
      
      # if(!g[i]) U[, i] <- tcrossprod( MAT[[i]]$u, (rmvn(1, rep(0, no), diag(MAT[[i]]$d))) )
      uu <- crossprod(U[, i])
      var <- se2/(uu + se2/tau2)
      Ytild <- Y - mu - as.matrix(U[, -i]) %*% sdu[-i]
      lp <- pnorm(0, mean = var * crossprod(U[, i], Ytild)/ se2 , sd = sqrt(var), lower.tail = F, log.p = T)
      
      p1 <-  log(pi_sdu) + lp + log(2) -1/2*log(1 + tau2*uu/se2) + 1/2*crossprod(Ytild, U[, i])^2 * var/se2^2
      p0 <-  log(1-pi_sdu)
      ratio <- p1 - sdumlogs(c(p1, p0))
      
      if(runif(1) < exp(ratio)){
        g[i] <- TRUE
        
        # var <- solve(1/tau2 + crossprod(U[, i])/ se2)
        var <- solve(1/tau2/se2 + crossprod(U[, i])/ se2)
        m <- crossprod(sign(xi[i])* U[, i], Ymu - as.matrix(U[, -i]) %*% sdu[-i]) * var / se2
        if(method == "folded"){
          xi[i] <- rnorm(1, mean = m, sd = sqrt(var))
        }
        if(method == "truncated"){
          xi[i] <-  extraDistr::rtnorm(1, mean = m , sd = sqrt(var), a = 0, b = Inf)
        }
        sdu[i] <- abs(xi[i])
        
        v <- diag(1/(sdu[i]^2/se2 + 1/ MAT[[i]]$d))
        # if(is.null(sim$U)){
        U[, i] <- tcrossprod(MAT[[i]]$u , rmvn(1, sdu[i]/se2* crossprod(v, crossprod(MAT[[i]]$u, (Ymu - as.matrix(U[, -i]) %*% sdu[-i] ))), sqrt(v), isChol = T))
        # }
      }else{
        g[i] <- FALSE
        sdu[i] <- xi[i] <- 0
        # if(is.null(sim$U)){
        U[, i] <- tcrossprod( MAT[[i]]$u, (rmvn(1, rep(0, no), diag(MAT[[i]]$d))) )
        # }
      }
    }
    pi_sdu <- rbeta(1, 1+sum(g), 1+q-sum(g));
    
    
    # 3. Mise a jour de U:
    if(sum(g) != 0){
      # for(j in my.sample(which(g))){  # 1:q
      #   v <- diag(1/(sdu[j]^2/se2 + 1/ MAT[[j]]$d))
      #   U[, j] <- tcrossprod(MAT[[j]]$u , rmvn(1, sdu[j]/se2* crossprod(v, crossprod(MAT[[j]]$u, (Ymu - as.matrix(U[, -j]) %*% sdu[-j] ))), sqrt(v), isChol = T))
      # }
      Usdu <- as.matrix(U[, g]) %*% sdu[g]
    }else{
      Usdu <- 0
    }
    
    
    # 4. Mise a jour de tau2 :
    tau2 <- pmax(1e-6, extraDistr::rinvgamma(1, sum(g)/2 + 1/2, 1/2 + 0.5*sum(sdu^2)/se2) )
    # tau2 <- 10
    
    # 5. Mise a jour de sigma^2:
    # se2 <- extraDistr::rinvgamma(1,  a + no/2,  b + crossprod(Ymu - Usdu) /2)
    se2 <- extraDistr::rinvgamma(1,  a + no/2 + sum(g)/2,  b + crossprod(Ymu - Usdu) /2 + sum(sdu^2)/2/tau2)
    
    
    if(it > burnin & it %% thinin == 0){
      chain$mu[ii] <- mu
      chain$se2[ii] <- se2
      chain$sdu[ii, ] <- sdu
      chain$s[ii, ] <- xi
      chain$tau2[ii] <- tau2
      chain$omega[ii, ] <- omega
      chain$logLik[ii, ] <- dnorm(Y, mu + U%*%sdu , sqrt(se2), log = TRUE)
      chain$g[ii, ] <- g
      chain$pi[ii] <- pi_sdu
      ii <- ii+1
    }
    
  }
  return(chain)
}




