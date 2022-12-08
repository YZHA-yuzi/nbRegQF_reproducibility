########################################################################
# --------------- R codes used in the simulation study --------------- #
# NOTE:
# Estimating quantile functions
# please make sure the working directory is the folder where 
# R script named "FUNs.R" and simulated data are located at 
########################################################################

# Load packages 
library(MASS); library(igraph); 
library(mvtnorm)
library(Matrix)
library(Rcpp)
library(dplyr)
library(tidyr)

# source self-defined R functions 
sourceCpp("FUNs_quan_Rcpp.cpp")
source("FUNs.R")

# read in simulated data
load("data_sim.rda")
x.sim.mat <- dat.all$exp$x.sim.mat

################## BEGIN ESTIMATION ################
set.seed(12345)

L = 4 # number of basis functions
shape = 5 # shape parameter specified in Gamma basis functions
ntimes = nrow(x.sim.mat) # number of time points

## get adjacency matrix ##
adj.list = as.matrix(cbind(c(1:(ntimes-1)), c(2:ntimes)))
W.t = get.adjacency(graph.edgelist(adj.list, directed=FALSE))
Dw.t <- Diagonal(x = apply(W.t, 1, sum))

## priors ##
## alpha0, theta_l_bar
c.prior = 0; C.prior = 100
## tau12, tau22
a = 0.1; b = 0.1
## rho1, rho2
# compute quantities beforehand to facilitate the updating of rho1, and rho2
lambda.t = eigen(solve(Dw.t)%*%W.t, only.values = TRUE)$values
rho1.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
ll.rho1.1 = sapply(rho1.prior.val, function(x) 0.5*sum(log(1-x*lambda.t)),
                   simplify = TRUE)

rho2.prior.val = qbeta(seq(1e-4, 1-1e-4, length.out = 1000), 1, 1)
ll.rho2.1 = sapply(rho2.prior.val, function(x) 0.5*L*sum(log(1-x*lambda.t)),
                   simplify = TRUE)

## tuning parameters
alph.tun = rep(0.1, ntimes)
theta.star.tun = matrix(rep(0.1, ntimes*L),
                        nrow = L, ncol = ntimes)

## the number of MCMC iterations 
niter = 10000
alph <- matrix(NA, ncol = ntimes, nrow = niter+1);
theta.star.list <- replicate(L,
                             matrix(NA, ncol = ntimes, nrow = niter+1),
                             simplify = F)
theta.bar.mat <- matrix(NA, ncol = L, nrow = niter + 1)
parm.mat <- matrix(NA, ncol = 5, nrow = niter + 1)
colnames(parm.mat) <- c("alpha", "tau12", "tau22",
                        "rho1", "rho2")

for(ll in 1:L){
  theta.star.list[[ll]][1, ] <- 1
}
alph[1, ] <- 0
theta.bar.mat[1, ] <- 1
parm.mat[1,  ] <- rep(0.5, ncol(parm.mat))
accept.mat <- matrix(NA, ncol = ntimes + ntimes, nrow = niter+1)
epsilon = 0.01


for(i in 1:niter){
  
  # 1. update alpha_t,0 (ntimes by 1) M-H
  theta.start.mat.i <- matrix(NA, nrow = L, ncol = ntimes)
  for(ll in 1:L){
    theta.start.mat.i[ll, ] <- theta.star.list[[ll]][i, ]
  }
  theta.mat.i <- apply(theta.start.mat.i, 2, ff, epsilon = epsilon)
  Omega1.i <- as.matrix(1/parm.mat[i, "tau12"]*(Dw.t - parm.mat[i,"rho1"]*W.t))

  ## current value of alpha = (alpha_10, ..., alpha_T0) (T x 1)
  alph.i <- alph[i, ]
  alph.val.i.vec <- rep(parm.mat[i, "alpha"], ntimes)
  
  ll.curr.vec <- comp_ll_alpha(alph_c = alph.i,
                               x_sim_mat = x.sim.mat,
                               theta_mat_i = theta.mat.i,
                               f = den.fun, L = L, shape = shape,
                               ntimes = ntimes, basis_fun = "Gamma")
  ll.curr.vec[!is.finite(ll.curr.vec)] <- log(1e-200)

  ### pre-compute the conditional mean and variance for alpha_t|alpha_-t ###
  alph.cond.mu <- alph.cond.sd <- rep(NA, ntimes)
  alph.bar.i = parm.mat[i, "alpha"]
  for(a.ii.cond in 1:ntimes){
    index2 = (1:ntimes)[-a.ii.cond]
    ### using precision matrix ###
    omega12.ii = Omega1.i[a.ii.cond, index2]
    alph.cond.mu[a.ii.cond] = alph.bar.i - (1/Omega1.i[a.ii.cond, a.ii.cond])*
      sum(omega12.ii*(alph.i[-a.ii.cond] - alph.bar.i))
    alph.cond.sd[a.ii.cond] = sqrt(1/Omega1.i[a.ii.cond, a.ii.cond])
  }

  for(a.ii in 1:ntimes){
    ## propose a candidate alpha_t0*
    alph.star.ii <- rnorm(1, alph.i[a.ii], sd = alph.tun[a.ii])
    alph.star <- alph.i
    alph.star[a.ii] <- alph.star.ii
    ## evaluate likelihood
    ## likelihood for the proposed alpha
    ll.star.t <- sum(log(den.fun(x = na.omit(x.sim.mat[a.ii, ]), L = L,
                                 theta.vec = theta.mat.i[ ,a.ii],
                                 alpha = alph.star[a.ii],
                                 shape = shape, basis.fun = "Gamma")))
    if(!is.finite(ll.star.t)){ll.star.t <- log(1e-200)}
    
    ll.star <- ll.star.t + dnorm(alph.star.ii,
                                 alph.cond.mu[a.ii],
                                 alph.cond.sd[a.ii], log = T)
    
    ll.curr <- ll.curr.vec[a.ii] + dnorm(alph.i[a.ii],
                                         alph.cond.mu[a.ii],
                                         alph.cond.sd[a.ii], log = T)
    
    ratio = min(1, exp(ll.star - ll.curr))
    if(ratio >= runif(1)){
      alph.i[a.ii] <- alph.star.ii
      ll.curr.vec[a.ii] <- ll.star.t
      accept.mat[i+1, a.ii] <- 1
    }else{
      accept.mat[i+1, a.ii] <- 0
    }
  } # end loop over alpha_t0, t = 1, ..., T
  alph[i+1, ] <- alph.i
  
  # 2. update theta_tl* t = 1, ..., ntimes, l = 1, ..., L
  ## update theta_t as a Lx1 vector
  Omega2.i.sparse <- 1/parm.mat[i, "tau22"]*(Dw.t - parm.mat[i,"rho2"]*W.t)
  Omega2.i <- as.matrix(Omega2.i.sparse)
  theta.start.mat.i <- matrix(NA, nrow = L, ncol = ntimes)
  for(ll in 1:L){
    theta.start.mat.i[ll, ] <- theta.star.list[[ll]][i, ]
  }
  
  ### pre-compute the conditional variance for theta_tl|theta_-tl ###
  theta.cond.sd <- sqrt(1/diag(Omega2.i))
  theta.bar.i = theta.bar.mat[i, ]
  for(the.ii in 1:ntimes){
    ## propose a new theta_lt*|current values
    theta.star.star.ii <- rnorm(L, theta.start.mat.i[,the.ii],
                                theta.star.tun[,the.ii])
    theta.star.ii <- ff(theta.star.star.ii, epsilon = epsilon)
    ll.star.t <- sum(log(den.fun(x = na.omit(x.sim.mat[the.ii, ]), L = L,
                                 theta.vec = theta.star.ii,
                                 alpha = alph.i[the.ii],
                                 shape = shape,
                                 basis.fun = "Gamma")))
    if(!is.finite(ll.star.t)){ll.star.t <- log(1e-200)}
    
    ### compute conditional mean of theta_tl|theta_-tl using precision matrix ###
    ### compute it within loop,
    ### since the conditional mean depends on current values of theta_-tl
    index2 = (1:ntimes)[-the.ii]
    omega12.ii = Omega2.i[the.ii, index2]
    H.inv = 1/Omega2.i[the.ii, the.ii]
    theta.cond.mu.t <- rep(NA, L)
    for(lll in 1:L){
      theta.cond.mu.t[lll] = theta.bar.i[lll] -
        H.inv*sum(omega12.ii*(theta.start.mat.i[lll, -the.ii] - theta.bar.i[lll]))
    }
    
    p.star = sum(dnorm(theta.star.star.ii, theta.cond.mu.t,
                       rep(theta.cond.sd[the.ii], L), log = T))
    p.curr = sum(dnorm(theta.start.mat.i[,the.ii], theta.cond.mu.t,
                       rep(theta.cond.sd[the.ii], L), log = T))
    ll.star <- ll.star.t + p.star
    ll.curr <- ll.curr.vec[the.ii] + p.curr
    ratio = min(1, exp(ll.star - ll.curr))
    if(ratio >= runif(1)){
      theta.start.mat.i[,the.ii] <- theta.star.star.ii
      accept.mat[i+1, ntimes + the.ii] <- 1
      ll.curr.vec[the.ii] <- ll.star.t
    }else{
      accept.mat[i+1, ntimes + the.ii] <- 0
    }
  } # end loop over time for theta_lt*
  for(ll in 1:L){
    theta.star.list[[ll]][i+1, ] <- theta.start.mat.i[ll, ]
  }

  # 3. update theta_bar_l, l = 1,..., L conjugate
  X.design.theta = matrix(rep(1, ntimes), ncol = 1)
  pre.comp = t(X.design.theta) %*% Omega2.i
  for(ll in 1:L){
    M = 1/(pre.comp %*% X.design.theta + 1/C.prior)
    m = M*(pre.comp%*%matrix(theta.start.mat.i[ll, ], ncol = 1))
    theta.bar.mat[i+1, ll] <- rnorm(1, m, sqrt(M))
  }
  
  # 4. update alpha_bar conjugate
  X.design.alph = matrix(rep(1, ntimes), ncol = 1)
  pre.comp = t(X.design.alph) %*% Omega1.i
  M = 1/(pre.comp %*% X.design.alph + 1/C.prior)
  m = M*(pre.comp%*%matrix(alph[i+1, ], ncol = 1))
  parm.mat[i+1, "alpha"] <- rnorm(1, m, sqrt(M))
  
  # 5. update tau12, conjugate, inverse gamma
  # variance parameter of the intercept
  res.alph.i <- alph[i+1, ] - parm.mat[i+1, "alpha"]
  parm.mat[i+1, "tau12"] <- 1/rgamma(1, ntimes/2 + a,
                                     b + as.numeric(t(matrix(res.alph.i, ncol = 1))%*%
                                                      (Dw.t - parm.mat[i, "rho1"]*W.t)%*%
                                                      matrix(res.alph.i, ncol = 1))/2)
  
  # 6. update tau22, conjugate, inverse gamma
  # variance parameter of coefs of basis functions
  innter.mat = (Dw.t - parm.mat[i, "rho2"]*W.t)
  sum.ig = 0
  for(ll in 1:L){
    res.theta.ll = theta.star.list[[ll]][i+1, ] - theta.bar.mat[i+1, ll]
    re = as.numeric(t(matrix(res.theta.ll, ncol = 1))%*%
                      innter.mat%*%
                      matrix(res.theta.ll, ncol = 1))
    sum.ig = sum.ig + re
  }
  parm.mat[i+1, "tau22"] <- 1/rgamma(1, L*ntimes/2 + a,
                                     b + sum.ig/2)
  
  # 7. update rho1 (M-H)
  inter = matrix(alph[i+1, ] - parm.mat[i+1, "alpha"], ncol = 1)
  inter1 = as.numeric( t(inter)%*%W.t%*%inter )
  ll.rho1 = ll.rho1.1 + rho1.prior.val/(2*parm.mat[i+1, "tau12"])*inter1
  parm.mat[i+1, "rho1"] <- sample(x = rho1.prior.val, size = 1,
                                  prob = exp(ll.rho1 - max(ll.rho1)))
  
  # 8. update rho2 (M-H)
  inter2 = 0
  for(ll in 1:L){
    inter = matrix(theta.star.list[[ll]][i+1, ] -
                     theta.bar.mat[i+1, ll], ncol = 1)
    inter1 = as.numeric( t(inter)%*%W.t%*%inter )
    inter2 = inter2 + inter1
  }
  ll.rho2 = ll.rho2.1 + rho2.prior.val/(2*parm.mat[i+1, "tau22"])*inter2
  parm.mat[i+1, "rho2"] <- sample(x = rho2.prior.val, size = 1,
                                  prob = exp(ll.rho2 - max(ll.rho2)))
  
  ## tuning parameters ##
  if(i <= 1000 & i%%100 == 0){
    accept.rate.all <- colMeans(accept.mat[2:i, ])
    for(acc.ii in 1:ntimes){
      accept.rate = accept.rate.all[acc.ii]
      if(accept.rate > 0.45){alph.tun[acc.ii] = 1.2*alph.tun[acc.ii]}
      if(accept.rate < 0.25){alph.tun[acc.ii] = 0.8*alph.tun[acc.ii]}
      acc.the.ii = acc.ii+ntimes
      accept.rate.the = accept.rate.all[acc.the.ii]
      if(accept.rate.the > 0.45){
        theta.star.tun[,acc.ii] = 1.2*theta.star.tun[,acc.ii]
      }
      if(accept.rate.the < 0.25){
        theta.star.tun[,acc.ii] = 0.8*theta.star.tun[,acc.ii]
      }
    }
  } # end tunning
  
} # end MCMC procedure

re  <- list(alph = alph,
            theta.bar = theta.bar.mat,
            parm = parm.mat,
            theta.star = theta.star.list,
            accept = accept.mat,
            alph.tun = alph.tun,
            theta.star.tun = theta.star.tun)
save(re, file = paste0("./inter_res/Res_estquan.rda"))
### !!!! NOTES !!!! 
### This estimation takes very long to run (~4 hours) 
### on a computer with 2.9 GHz 6-Core Intel Core i9 and 32 Gb memory
### SAVE the results as intermediate results ####
### ------- summarize results ---------- ###

## read in estimation results using relative path ###
load(paste0("./inter_res/Res_estquan.rda"))

### get var-cov of MVN priors for alpha, theta from its posterior samples ###
burn_in = 5000
ntimes = 1000; 
epsilon = 0.01

alpha.t.post <- re$alph[-c(1:burn_in), ]
theta.star.list.post <- lapply(re$theta.star, function(x) x[-c(1:burn_in), ])
alpha.theta.vec.post <- matrix(NA, ncol = ntimes*(L+1), nrow = nrow(alpha.t.post))
for(isamp in 1:nrow(alpha.t.post)){
  alpha.t.i = alpha.t.post[isamp, ]
  theta.star.mat.i = rbind(theta.star.list.post[[1]][isamp, ],
                           theta.star.list.post[[2]][isamp, ],
                           theta.star.list.post[[3]][isamp, ],
                           theta.star.list.post[[4]][isamp, ])
  theta.mat.i = apply(theta.star.mat.i, 2, ff, epsilon = epsilon)
  inter = rbind(alpha.t.i, theta.star.mat.i)
  alpha.theta.vec.post[isamp, ] <- c(inter)
  if(isamp %% 1000 == 0){cat(isamp, ",")}
}
mu.post = colMeans(alpha.theta.vec.post)

Sigma.list <- list()
for(itime in 1:ntimes){
  alpha.theta.post.itime <- cbind(alpha.t.post[,itime],
                                  theta.star.list.post[[1]][,itime],
                                  theta.star.list.post[[2]][,itime],
                                  theta.star.list.post[[3]][,itime],
                                  theta.star.list.post[[4]][,itime])
  Sigma.list[[itime]] <- cov(alpha.theta.post.itime)
}
Sigma.post <- bdiag(Sigma.list)
matrix_norm_post <- list(M = mu.post, Sigma.bdiag = Sigma.post)
save(matrix_norm_post, file = paste0("./inter_res/mu_cov_MVNprior.rda"))


