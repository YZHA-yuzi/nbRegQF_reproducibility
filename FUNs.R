########################################################################
# ------ Self-defined R functions used in the simulation study --------#
########################################################################

ff <- function(x, epsilon){
  x[x <= epsilon] <- epsilon
  return(x)
}


## Basis functions given in Equations (3) and (4) in Reich (2012) ##
Bl <- function(tau.vec, l, L, basis.fun, shape = 5){
  k.vec = seq(0, 1, length.out = L+1)
  if(basis.fun == "Gau"){
    B.vec <- c()
    if(k.vec[l] < 0.5){
      for(tau in tau.vec){
        if(tau < k.vec[l]){
          B = qnorm(k.vec[l]) - qnorm(k.vec[l+1]) 
        }else if(k.vec[l] <= tau & tau < k.vec[l+1]){
          B = qnorm(tau) - qnorm(k.vec[l+1]) 
        }else{
          B = 0
        }
        B.vec <- c(B.vec, B)
      }
    }else{
      for(tau in tau.vec){
        if(tau < k.vec[l]){
          B = 0
        }else if(k.vec[l] <= tau & tau < k.vec[l+1]){
          B = qnorm(tau) - qnorm(k.vec[l]) 
        }else{
          B = qnorm(k.vec[l+1]) - qnorm(k.vec[l])
        }
        B.vec <- c(B.vec, B)
      }
    }
    return(B.vec)
  }
  if(basis.fun == "Gamma"){
    B.vec <- c()
    if(k.vec[l] < 0.5){
      for(tau in tau.vec){
        if(tau < k.vec[l]){
          B = qgamma(k.vec[l], shape = shape) - qgamma(k.vec[l+1], shape = shape) 
        }else if(k.vec[l] <= tau & tau < k.vec[l+1]){
          B = qgamma(tau, shape = shape) - qgamma(k.vec[l+1], shape = shape) 
        }else{
          B = 0
        }
        B.vec <- c(B.vec, B)
      }
    }else{
      for(tau in tau.vec){
        if(tau < k.vec[l]){
          B = 0
        }else if(k.vec[l] <= tau & tau < k.vec[l+1]){
          B = qgamma(tau, shape = shape) - qgamma(k.vec[l], shape = shape) 
        }else{
          B = qgamma(k.vec[l+1], shape = shape) - qgamma(k.vec[l], shape = shape)
        }
        B.vec <- c(B.vec, B)
      }
    }
    return(B.vec)
  }
  
}


## define quantile function based on those basis functions ##
quan.fun <- function(tau, L, theta.vec, alpha, basis.fun, shape = 5){
  if(basis.fun == "Gau"){
    B.mat <- matrix(NA, nrow = L, ncol = length(tau))
    for(i in 1:L){
      B.mat[i, ] <- sapply(tau, Bl, l = i, L = L, 
                           basis.fun = basis.fun,
                           shape = shape,
                           simplify = T)
    }
    theta.vec <- matrix(theta.vec, nrow = 1)
    Q.sim <- (alpha + theta.vec%*%B.mat)
    return(Q.sim)
  }
  if(basis.fun == "Gamma"){
    B.mat <- matrix(NA, nrow = L, ncol = length(tau))
    for(i in 1:L){
      B.mat[i, ] <- sapply(tau, Bl, l = i, L = L, 
                           basis.fun = basis.fun,
                           shape = shape,
                           simplify = T)
    }
    theta.vec <- matrix(theta.vec, nrow = 1)
    Q.sim <- (alpha + theta.vec%*%B.mat)
    return(Q.sim)
  }
}


## A function to evaluate density function ##
den.fun <- function(x, L, theta.vec, alpha, basis.fun, shape = 5){
  k.vec = seq(0, 1, length.out = L+1)
  if(basis.fun == "Gau"){
    q.vec <- quan.fun(tau = k.vec, L = L, 
                      theta.vec = theta.vec,
                      alpha = alpha,
                      basis.fun = basis.fun,
                      shape = shape)
    interval.mat <- cbind(q.vec[1:L], q.vec[2:(L+1)])
    a.vec <- get_basis_a(k.vec, q.vec, theta.vec)
    den.vec <- get_den(x, interval.mat, avec = a.vec,
                       thetavec = theta.vec)
    return(den.vec)
  }
  if(basis.fun == "Gamma"){
    q.vec <- quan.fun(tau = k.vec, L = L, 
                      theta.vec = theta.vec,
                      alpha = alpha, 
                      basis.fun = basis.fun,
                      shape = shape)
    interval.mat <- cbind(q.vec[1:L], q.vec[2:(L+1)])
    a.vec <- get_basis_a_gamma(k.vec, q.vec, theta.vec, shape)
    den.vec <- get_den_gamma(x, interval.mat, avec = a.vec,
                             thetavec = theta.vec, shape = shape)
    return(den.vec)
  }
}


### orthnormal Bernstein basis polynomials of degree n ###
### Eqn. (7)
bernstein_orth <- function(tau, j, n){
  comp1 = sqrt(2*(n-j)+1)
  comp2 = (1-tau)^(n-j)
  k.vec = 0:j
  comp3 <- function(x){
    sum((-1)^(k.vec)*choose(2*n+1-k.vec, j-k.vec)*choose(j, k.vec)*x^(j-k.vec))
  } 
  comp4  = sapply(tau, comp3, simplify = T)
  return(comp1*comp2*comp4)
}

get_integration_bernstein <- function(L, n, basis.fun, shape = 5){
  ## INPUTS:
  ## L: number of basis functions in representating quantile functions
  ## n: degrees of Bernstein basis functions
  ## basis.fun: name of the basis function used for estimating quan
  ## shape: shape parameter for gamma basis fun used for estimating quan
  ## OUTPUE: (L+1)xK (B(tau)K(.,tau)^T) matrix 
  K = n+1
  integration.mat <- matrix(NA, ncol = K, nrow = L + 1)
  for(ll in 1:(L+1)){
    if(ll == 1){
      for(k in 0:(K-1)){
        intergral <- function(x, df){
          bernstein_orth(tau = x, j = k, n = df)
        }
        integration.mat[ll, k+1] <- integrate(intergral, lower = 0, upper = 1,
                                              df = n)$value
      }
    }else{
      for(k in 0:(K-1)){
        intergral <- function(x, l, L, df){
          Bl(tau.vec = x, l = l, L = L, 
             basis.fun = basis.fun,
             shape = shape)*bernstein_orth(tau = x, j = k, n = df)
        }
        integration.mat[ll, k+1] <- integrate(intergral, lower = 0, upper = 1,
                                              l = ll-1, L = L,
                                              df = n)$value
      }
    }
  }
  return(integration.mat)
}


fit.health.Bernstein <- function(y,
                                 X.design,
                                 niter, ntimes, 
                                 burn_in = 20000,
                                 returnpost = FALSE){
  ## INPUTs:
  ## y: observed counts 
  ## X.design: design matrix
  ## niter (# MCMC iterations), ntimes (# time points)
  
  # ------------ priors ----------- #
  # priors #
  ## beta 
  nbeta.input = ncol(X.design)
  c.beta = matrix(rep(0, nbeta.input), ncol = 1)
  C.beta.pre = diag(rep(1/100, nbeta.input))
  
  # tuning parameters #
  xi.tun = 0.2
  
  ## initialize data frame to store posterior samples ##
  omega.mat <- matrix(NA, ncol = length(y), nrow = niter + 1)
  beta.mat <- matrix(NA, ncol = nbeta.input, nrow = niter + 1)
  parm.mat <- as.data.frame(matrix(NA, ncol = 1, nrow = niter + 1))
  colnames(parm.mat) <- c("xi")
  accept.mat <- matrix(NA, ncol = 1, nrow = niter + 1)
  colnames(accept.mat) <- c("xi")
  
  parm.mat[1, ] <- c(5)
  
  ## dimension of BK.mat = int B(tau)K(tau)^T dtau  = (L+1)xK
  dat.sim = data.frame(y = y, X.design)
  colnames(dat.sim) <- c("y", "Int", paste0("x", 1:(nbeta.input-1)))
  beta.mat[1, ] <- rep(0, nbeta.input)
  
  pg.parm1 <- dat.sim$y + parm.mat[1, "xi"]
  pg.parm2 <- as.matrix(dat.sim[,c("Int", paste0("x", 1:(nbeta.input-1)))])%*%
    matrix(beta.mat[1, ], ncol = 1)
  ntotal = ntimes
  omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
  
  for(i in 1:niter){
    # 1. update omega 
    omega.vec <- omega.mat[i, ]
    Omega <- Diagonal(x = omega.vec)
    
    ## compute latent normal variables ##
    z.vec <- (dat.sim$y - parm.mat[i, "xi"])/(2*omega.vec)
    
    # 2. update  beta (conjugate)
    ## posterior mean and variance of beta ##
    M <- solve(crossprod(X.design*sqrt(omega.vec)) + C.beta.pre)
    m <- M%*%(C.beta.pre%*%c.beta + 
                t(sqrt(omega.vec)*X.design)%*%
                (sqrt(omega.vec)*matrix(z.vec, ncol = 1)))
    beta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = m, Sigma = M) )
    
    # 3. update xi (over dispersion) 
    eta.vec.i1 = X.design%*%matrix(beta.mat[i+1, ], ncol = 1)
    q.nb = 1/(1+exp(eta.vec.i1)) # 1 - Pr(success)
    xi.star = rtruncnorm(1, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun)
    
    ll.star = sum(dnbinom(dat.sim$y, size = xi.star, prob = q.nb, log = TRUE))
    ll.curr = sum(dnbinom(dat.sim$y, size = parm.mat[i, "xi"], prob = q.nb, log = TRUE))
    q.star = log( dtruncnorm(parm.mat[i, "xi"], a = 0, mean = xi.star, sd = xi.tun) )
    q.curr = log( dtruncnorm(xi.star, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun) )
    
    ratio = min(1, exp(ll.star + q.star - ll.curr - q.curr))
    if(ratio > runif(1)){
      parm.mat[i+1, "xi"] <- xi.star
      accept.mat[i+1, 1] = 1
    }else{
      parm.mat[i+1, "xi"] <- parm.mat[i, "xi"]
      accept.mat[i+1, 1] = 0
    }
    
    omega.mat[i+1, ] <- rpg(ntotal, dat.sim$y + parm.mat[i+1, "xi"], eta.vec.i1)
    
    if(i%%1000 == 0){cat(i)}
    
    ## tuning parameters ##
    if(i <= 5000 & i%%100 == 0){
      
      accept.rate <- mean(accept.mat[2:i, 1])
      if(accept.rate > 0.6){xi.tun = 1.1*xi.tun}
      if(accept.rate < 0.2){xi.tun = 0.8*xi.tun}
      
    }
    
  }# END MCMC iterations 
  
  ## compute WAIC ##
  ## WAIC = -2*(lppd - pWAIC); 
  ## lppd (log pointwise predictive density); 
  ## pWAIC (effective number of parameters)
  compq <- function(x, X.design){
    q = 1/(1+exp(X.design%*%matrix(x[1:(nbeta.input)], ncol = 1)))
  }
  q.vec = apply(beta.mat[-(1:burn_in), ], 1, compq, X.design = X.design)
  
  compll <- function(x, xi.post){
    ll.i = dnbinom(x[1], size = xi.post, prob = x[-1], log = TRUE)
    return(ll.i)
  }
  log.pd = t( apply(cbind(dat.sim$y, q.vec), 1, compll, 
                    xi.post = parm.mat[-c(1:burn_in),"xi"]) )
  
  lppd = sum(log(apply(exp(log.pd), 1, mean)))
  pWAIC = sum(apply(log.pd, 1, var))
  
  if(returnpost){
    re <- list(beta = beta.mat[-c(1:burn_in), ], 
               parm = parm.mat[-c(1:burn_in), ], 
               accept.rate = mean(accept.mat[-c(1:burn_in), ]), 
               WAIC = -2*(lppd - pWAIC), lppd = lppd, pWAIC = pWAIC)
  }
  if(! returnpost){
    re <- list(beta = beta.mat, parm = parm.mat, 
               accept.rate = mean(accept.mat[-c(1:burn_in), ]), 
               WAIC = -2*(lppd - pWAIC), lppd = lppd, pWAIC = pWAIC)
  }
  
  return(re)
  
}


### A function to find initial values of thetas using quadratic programming ###
find.initial <- function(y.obs, Bl.mat, L, tau.vec){
  Dmat <- t(Bl.mat) %*% Bl.mat
  dvec <- t(Bl.mat) %*% matrix(y.obs, ncol = 1)
  Amat <- diag(rep(1, L+1))
  bvec <- matrix(0, ncol = L + 1)
  theta.initial <- solve.QP(Dmat = Dmat, dvec = dvec,
                            Amat = Amat, bvec = bvec)$solution
  return(theta.initial)
}


##### introduce a "prior" for alpha_t,0 and theta_t,l ######
fit.health.Bernstein.var <- function(y,
                                     BK.mat,
                                     theta.mat.ini,
                                     theta.pri.mu,
                                     theta.pri.pre,
                                     niter, ntimes, 
                                     nbeta.input,
                                     burn_in = 2500){
  ## INPUTs:
  ## y: observed counts 
  ## BK.mat: int B(tau)phi(tau) dtau (L+1)x(K)
  ## theta.mat.ini: initialized values of 
  ## (alpha0, theta_1, ... ,theta_L) (Tx(L+1))
  ## theta.pri.mu: mean of the vectorized matrix 
  ## (alpha0, theta_1, ... ,theta_L) stacked by row
  ## theta.pri.pre: precision matrix of the vectorized matrix 
  ## (alpha0, theta_1, ... ,theta_L)
  ## nbeta.input: # of betas
  ## niter (# MCMC iterations), ntimes (# time points)
  
  # ------------ priors ----------- #
  # priors #
  ## beta 
  # nbeta.input = ncol(X.design)
  c.beta = matrix(rep(0, nbeta.input), ncol = 1)
  C.beta.pre = diag(rep(1/100, nbeta.input))
  
  # tuning parameters #
  xi.tun = 0.2
  
  ## initialize data frame to store posterior samples ##
  omega.mat <- matrix(NA, ncol = length(y), nrow = niter + 1)
  beta.mat <- matrix(NA, ncol = nbeta.input, nrow = niter + 1)
  parm.mat <- as.data.frame(matrix(NA, ncol = 1, nrow = niter + 1))
  colnames(parm.mat) <- c("xi")
  accept.mat <- matrix(NA, ncol = 1, nrow = niter + 1)
  colnames(accept.mat) <- c("xi")
  theta.mat.list <- replicate(1 + niter, 
                              matrix(NA, ncol = L + 1, nrow = ntimes), 
                              simplify = F)
  
  parm.mat[1, ] <- c(5)
  theta.mat.list[[1]] <- theta.mat.ini
  
  ## dimension of BK.mat = int B(tau)K(tau)^T dtau  = (L+1)xK
  X.design.curr = cbind(1, theta.mat.list[[1]]%*%BK.mat)
  dat.sim = data.frame(y = y, X.design.curr)
  colnames(dat.sim) <- c("y", "Int", paste0("x", 1:(nbeta.input-1)))
  beta.mat[1, ] <- rep(0, nbeta.input)
  
  pg.parm1 <- dat.sim$y + parm.mat[1, "xi"]
  pg.parm2 <- as.matrix(dat.sim[,c("Int", 
                                   paste0("x", 1:(nbeta.input-1)))])%*%
    matrix(beta.mat[1, ], ncol = 1)
  ntotal = ntimes
  omega.mat[1, ] <- rpg(ntotal, pg.parm1, pg.parm2)
  
  theta.pri.prod = theta.pri.pre%*%theta.pri.mu
  
  L = nrow(BK.mat) - 1
  theta.pri.pre.list <- NULL
  theta.pri.mu.list <- NULL
  for(t in 1:ntimes){
    indext = ((L+1)*(t-1) + 1):((L+1)*t)
    theta.pri.pre.list[[t]] <- as.matrix(theta.pri.pre[indext, indext])
    theta.pri.mu.list[[t]] <- theta.pri.mu[indext]
  }
  theta.pri.prod.list <- lapply(1:ntimes, 
                          function(x) 
                          theta.pri.pre.list[[x]]%*%theta.pri.mu.list[[x]])
  
  
  for(i in 1:niter){
    # 1. update omega 
    omega.vec <- omega.mat[i, ]
    Omega <- Diagonal(x = omega.vec)
    
    ## compute latent normal variables ##
    z.vec <- (dat.sim$y - parm.mat[i, "xi"])/(2*omega.vec)
    
    # 2. update  beta (conjugate)
    ## posterior mean and variance of beta ##
    M <- solve(crossprod(X.design.curr*sqrt(omega.vec)) + C.beta.pre)
    m <- M%*%(C.beta.pre%*%c.beta + 
                t(sqrt(omega.vec)*X.design.curr)%*%
                (sqrt(omega.vec)*matrix(z.vec, ncol = 1)))
    beta.mat[i+1, ] <- as.numeric( mvrnorm(1, mu = m, Sigma = M) )
    
    # 3. update vectorized theta.mat 
    z.theta.vec <- z.vec - beta.mat[i+1, 1]
    BK.beta.i <- matrix(t(BK.mat%*%matrix(beta.mat[i+1, -1], ncol = 1)), 
                        nrow = 1)
    ## posterior mean and variance of vectorized theta.mat ##
    M.list <- get_Mlist(BK_beta = BK.beta.i,
                        omega_vec = omega.vec,
                        theta_pre = theta.pri.pre.list)
    m.list <- get_mulist(theta_prod = theta.pri.prod.list,
                         M_list = M.list,
                         omega_vec = omega.vec,
                         z = z.theta.vec,
                         BK_beta = BK.beta.i)
    theta.mat.list[[i+1]] <- sample_basiscoef(mu_list = m.list,
                                              sigma_list = M.list,
                                              ntheta = L+1)
    
    # 3. update xi (over dispersion) 
    X.design.curr = cbind(1, theta.mat.list[[i+1]]%*%BK.mat)
    eta.vec.i1 = X.design.curr%*%matrix(beta.mat[i+1, ], ncol = 1)
    q.nb = 1/(1+exp(eta.vec.i1)) # 1 - Pr(success)
    xi.star = rtruncnorm(1, a = 0, mean = parm.mat[i, "xi"], sd = xi.tun)
    
    ll.star = sum(dnbinom(dat.sim$y, size = xi.star, 
                          prob = q.nb, log = TRUE))
    ll.curr = sum(dnbinom(dat.sim$y, size = parm.mat[i, "xi"], 
                          prob = q.nb, log = TRUE))
    q.star = log( dtruncnorm(parm.mat[i, "xi"], a = 0,
                             mean = xi.star, sd = xi.tun) )
    q.curr = log( dtruncnorm(xi.star, a = 0, 
                             mean = parm.mat[i, "xi"], sd = xi.tun) )
    
    ratio = min(1, exp(ll.star + q.star - ll.curr - q.curr))
    if(ratio > runif(1)){
      parm.mat[i+1, "xi"] <- xi.star
      accept.mat[i+1, 1] = 1
    }else{
      parm.mat[i+1, "xi"] <- parm.mat[i, "xi"]
      accept.mat[i+1, 1] = 0
    }
    
    omega.mat[i+1, ] <- rpg(ntotal, 
                            dat.sim$y + parm.mat[i+1, "xi"], 
                            eta.vec.i1)
    
    if(i%%1000 == 0){cat(i)}
    
    ## tuning parameters ##
    if(i <= 5000 & i%%100 == 0){
      
      accept.rate <- mean(accept.mat[2:i, 1])
      if(accept.rate > 0.6){xi.tun = 1.1*xi.tun}
      if(accept.rate < 0.2){xi.tun = 0.8*xi.tun}
      
    }
    
  }# END MCMC iterations 
  

  ## compute WAIC ##
  ## WAIC = -2*(lppd - pWAIC); 
  ## lppd (log pointwise predictive density); 
  ## pWAIC (effective number of parameters)
  compq <- function(beta.vec, theta.mat, BK.mat){
    X.design = cbind(1, theta.mat%*%BK.mat)
    eta.vec = X.design%*%matrix(beta.vec, ncol = 1)
    q = 1/(1+exp(eta.vec))
    return(as.numeric(q))
  }
  beta.post <- beta.mat[-(1:burn_in), ]
  theta.mat.post <- theta.mat.list[-c(1:burn_in)]
  q.vec <- matrix(NA, ncol = nrow(beta.post), nrow = ntimes)
  for(isamp in 1:nrow(beta.post)){
    q.vec[,isamp] <- compq(beta.vec = beta.post[isamp, ],
                           theta.mat = theta.mat.post[[isamp]],
                           BK.mat = BK.mat)
  }
  
  compll <- function(x, xi.post){
    ll.i = dnbinom(x[1], size = xi.post, prob = x[-1], log = TRUE)
    return(ll.i)
  }
  log.pd = t( apply(cbind(dat.sim$y, q.vec), 1, compll, 
                    xi.post = parm.mat[-c(1:burn_in),"xi"]) )
  
  lppd = sum(log(apply(exp(log.pd), 1, mean)))
  pWAIC = sum(apply(log.pd, 1, var))
  
  re <- list(beta = beta.mat, parm = parm.mat,
             theta.mat.list = theta.mat.list[-c(1:burn_in)],
             accept.rate = mean(accept.mat[-c(1:burn_in), ]),
             WAIC = -2*(lppd - pWAIC), lppd = lppd, pWAIC = pWAIC)
  
  return(re)
  
}


## A function to compute coefs of orthnormal Bernstein basis functions used for
## approximating given function 
## Eqn. (22)
get_beta_bernstein <- function(j, n, f){
  ## INPUT:
  ## n: degree of orthnormal Bernstein basis functions
  ## f: function to be approximated 
  integral <- function(x){
    bernstein_orth(x, j, n)*f(x)
  }
  re = integrate(integral, lower = 0, upper = 1)$value
  return(re)
}


## compute int_beta(tau)dtau based on orthnormal Bernstein basis polynomials ##
get_int_beta_bernstein <- function(n, coef){
  intergral <- function(x, n, coef){
    X <- NULL
    for(j in 0:n){
      X <- cbind(X, bernstein_orth(tau = x, j = j, n = n))
    }
    return( X%*%matrix(coef, ncol = 1) )
  }
  re <- integrate(intergral, lower = 0, upper = 1,
                  n = n, coef = coef)$value
  return(re)
}




########### FUNCTIONS TO SUMMARIZE RESULTS ###########
## A function to get the matrix containing basis functions
get_bernstein_mat <- function(tau, n){
  X <- NULL
  for(j in 0:n){
    X <- cbind(X, bernstein_orth(tau = tau, j = j, n = n))
  }
  return(X)
}

### A function to summarize results 6:96 ###
get_res_QFknown <- function(re, tau.vec, n, burn_in){
  beta.post <- re$beta[-c(1:burn_in), ]
  parm.post <- re$parm[-c(1:burn_in), ]
  
  beta.hat <- colMeans(beta.post)
  xi.hat <- mean(parm.post)
  
  ## compute estimated beta(tau) = sum_i(K(tau_i, tau)*alpha_i)
  Kmat <- get_bernstein_mat(tau = tau.vec, n = n)
  beta.tau.hat <- Kmat%*%matrix(beta.hat[-1], ncol = 1)
  num.beta = ncol(beta.post)
  beta.tau.post <- apply(beta.post, 1,
                         function(x) Kmat%*%matrix(x[2:num.beta], ncol=1))
  beta.tau.CI <- apply(beta.tau.post, 1, quantile, c(0.025, 0.975))
  
  beta.tau.est.mat = data.frame(quan.level = tau.vec,
                                est = beta.tau.hat,
                                sd = apply(beta.tau.post, 1, sd),
                                lower = beta.tau.CI[1, ],
                                upper = beta.tau.CI[2, ])
  
  beta.int.hat = get_int_beta_bernstein(n = n, coef = beta.hat[-1])
  get_int_bernstein <- function(n){
    re.vec <- c()
    for(j in 0:n){
      re.vec <- c(re.vec, integrate(bernstein_orth, lower = 0, upper = 1,
                                    j = j, n = n)$value)
    }
    return(re.vec)
  }
  inter = get_int_bernstein(n = n)
  
  beta.int.post = apply(beta.post[,-1], 1, function(x) sum(inter*x))
  beta.int.est.mat = data.frame(est = mean(beta.int.post),
                                sd  = sd(beta.int.post),
                                lower = quantile(beta.int.post, 0.025),
                                upper = quantile(beta.int.post, 0.975))
  
  re <- list(n.fit = n,
             beta.int.est.mat = beta.int.est.mat,
             beta.tau.est.mat = beta.tau.est.mat)
}

get_res_mean <- function(re, burn_in){
  beta.post <- re$beta[-c(1:burn_in), ]
  parm.post <- re$parm[-c(1:burn_in), ]
  parm.post <- cbind(beta.post, parm.post)
  tab = data.frame(est = colMeans(parm.post),
                   sd = apply(parm.post, 2, sd),
                   lower = apply(parm.post, 2, quantile, 0.025),
                   upper = apply(parm.post, 2, quantile, 0.975))
  return(tab)
}
