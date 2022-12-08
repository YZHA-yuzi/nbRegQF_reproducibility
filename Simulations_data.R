########################################################################
# --------------- R codes used in the simulation study --------------- #
# NOTE:
# generating data used in simulation studies
# please make sure the working directory is the folder where 
# files "FUNs.R" and "FUNs_quan_Rcpp.cpp" are
########################################################################

# Load packages 
library(igraph); library(mvtnorm); library(Matrix)
library(spam); library(Rcpp)

# source self-defined R functions 
sourceCpp("FUNs_quan_Rcpp.cpp")
source("FUNs.R")

# 1. DATA GENERATION 
## 1.1 generate exposure data 
## following quantile functions and health model given in Section 4 

### set random seed 
set.seed(13579)

# the number of basis functions used to define exposure quantile functions
L = 4
# the shape parameter of used piecewise Gamma functions 
shape = 5
# number of time points
ntimes = 1000 

# specify parameters 
alph.bar.true = 7.2 # theta_0 given in Section 4
tau12 = 1 # sigma_0^2 given in Section 4
rho1 = 0.9 
theta.bar.true = c(0.9,0.9,0.9,0.9) # theta_l, l=1,...,4 given in Section 4
tau22 = 0.02 # sigma_0^2 given in Section 4
rho2 = 0.9

# get an adjacency matrix  
adj.list = as.matrix(cbind(c(1:(ntimes-1)), c(2:ntimes)))
W.t = get.adjacency(graph.edgelist(adj.list, directed=FALSE))
Dw.t <- Diagonal(x = apply(W.t, 1, sum))
Sigma1 = tau12*solve(Dw.t - rho1*W.t)
Sigma2 = tau22*solve(Dw.t - rho2*W.t)

# generate theta_0,i, i=1,...,1000 given in Section 4
alph.vec.true = rmvnorm(1, mean = rep(alph.bar.true, ntimes),
                        sigma = as.matrix(Sigma1))
# first generate theta_l,i* 
theta.mat.true = matrix(NA, ncol = ntimes, nrow = L)
for(i in 1:L){
  inter.vec = rmvnorm(1, mean = rep(theta.bar.true[i], ntimes),
                      sigma = as.matrix(Sigma2))
  theta.mat.true[i, ] <- inter.vec
}
## theta_l,i cannot be small than 0 
theta.mat.true[theta.mat.true < 0] = 0.01

### randomly select 9 time points, 
### to check how time-specific quantile functions deviate from the overall one  
prob.vec.true <- quan.fun(tau = seq(0, 1, 0.01),
                          L = L, theta.vec = theta.bar.true,
                          alpha = alph.bar.true,
                          basis.fun = "Gamma",
                          shape = shape)
index = sample(1:ntimes, 9)
mean.vec <- c()
prob.vec.sim <- NULL
for(iii in index){
  prob.vec.sim <- rbind(prob.vec.sim,
                        quan.fun(tau = seq(0, 1, 0.01),
                                 L = L, theta.vec = theta.mat.true[,iii],
                                 alpha = alph.vec.true[iii],
                                 basis.fun = "Gamma", shape = shape))
  mean.vec <- c(mean.vec, integrate(quan.fun,
                                    L = L, theta.vec = theta.mat.true[,iii],
                                    alpha = alph.vec.true[iii],
                                    basis.fun = "Gamma", shape = shape,
                                    lower = 0, upper = 1)$value)
}


# number of observations at time point t
nsize = 100  
# specify ranges of individual-level exposures
x.vec = seq(0, 100, 0.005)
# generate individual-level exposures at each time point 
x.sim.mat <- matrix(NA, nrow = ntimes, ncol = nsize)
for(i in 1:ntimes){
  prob.vec = den.fun(x = x.vec, L = L,
                     theta.vec = theta.mat.true[, i],
                     alpha = alph.vec.true[i],
                     basis.fun = "Gamma",
                     shape = shape)
  x.sim = sample(x.vec, nsize, prob = prob.vec)
  x.sim.mat[i, ] <- x.sim
}
dat.exp <- list(alph.vec.true = alph.vec.true,
                theta.mat.true = theta.mat.true,
                x.sim.mat = x.sim.mat)
dat.all <- NULL
dat.all[["exp"]] <- dat.exp




## 1.2 generate health data 
# the number of simulations
nsims = 100
# intercept in the health model
beta0 = -3.5
# the over-dispersion parameter
xi = 20 

### Scenario 1 beta(tau) is a constant = 0.5 ###
intergral <- function(x, L, theta.vec, alpha, shape){
  re.beta = 0.5
  re.quan = as.numeric(quan.fun(x, L = L,
                                theta.vec = theta.vec,
                                alpha = alpha,
                                basis.fun = "Gamma",
                                shape = shape))
  re = re.beta*re.quan
  return(re)
}
# compute int beta(tau)Qt(tau)dtau
integration.vec = rep(NA, ntimes)
for(i in 1:ntimes){
  integration.vec[i] <- integrate(intergral, lower = 0, upper = 1,
                                  L = L, theta.vec = dat.exp$theta.mat.true[,i],
                                  alpha = dat.exp$alph.vec.true[i],
                                  shape = shape)$value
  if(i%%100 == 0){cat(i, ",")}
}
eta.vec.true <- beta0 + integration.vec
q.nb = 1/(1+exp(eta.vec.true))
set.seed(1234)
y.sim = matrix(NA, ncol = nsims, nrow = ntimes)
for(i.sim  in 1:nsims){
  y.sim[,i.sim] <- apply(cbind(rep(xi, ntimes), q.nb), 1,
                         function(x) rnbinom(1, size = x[1], prob = x[2]))
}
dat.all[["health1"]] <- y.sim  # each column is a simulation


### Scenario 2 beta(tau) = x ####
# compute int beta(tau)Qt(tau)dtau
intergral <- function(x, L, theta.vec, alpha, shape){
  re.beta = x
  re.quan = as.numeric(quan.fun(x, L = L,
                                theta.vec = theta.vec,
                                alpha = alpha,
                                basis.fun = "Gamma",
                                shape = shape))
  re = re.beta*re.quan
  return(re)
}
integration.vec = rep(NA, ntimes)
for(i in 1:ntimes){
  integration.vec[i] <- integrate(intergral, lower = 0, upper = 1,
                                  L = L, theta.vec = dat.exp$theta.mat.true[,i],
                                  alpha = dat.exp$alph.vec.true[i],
                                  shape = shape)$value
  if(i%%100 == 0){cat(i, ",")}
}
eta.vec.true <- beta0 + integration.vec
q.nb = 1/(1 + exp(eta.vec.true))
set.seed(1234)
y.sim = matrix(NA, ncol = nsims, nrow = ntimes)
for(i.sim  in 1:nsims){
  y.sim[,i.sim] <- apply(cbind(rep(xi, ntimes), q.nb), 1,
                         function(x) rnbinom(1, size = x[1], prob = x[2]))
}
dat.all[["health2"]] <- y.sim  # each column is a simulation


### Scenario 3 beta(tau) quadratic effects 1.5x^2 ###
# compute int beta(tau)Qt(tau)dtau
intergral <- function(x, L, theta.vec, alpha, shape){
  re.beta = 1.5*x^2 
  re.quan = as.numeric(quan.fun(x, L = L,
                                theta.vec = theta.vec,
                                alpha = alpha,
                                basis.fun = "Gamma",
                                shape = shape))
  re = re.beta*re.quan
  return(re)
}
integration.vec = rep(NA, ntimes)
for(i in 1:ntimes){
  integration.vec[i] <- integrate(intergral, lower = 0, upper = 1,
                                  L = L, theta.vec = dat.exp$theta.mat.true[,i],
                                  alpha = dat.exp$alph.vec.true[i],
                                  shape = shape)$value
  if(i%%100 == 0){cat(i, ",")}
}
eta.vec.true <- beta0 + integration.vec
q.nb = 1/(1 + exp(eta.vec.true))
set.seed(1234)
y.sim = matrix(NA, ncol = nsims, nrow = ntimes)
for(i.sim  in 1:nsims){
  y.sim[,i.sim] <- apply(cbind(rep(xi, ntimes), q.nb), 1,
                         function(x) rnbinom(1, size = x[1], prob = x[2]))
}
dat.all[["health3"]] <- y.sim  # each column is a simulation


### Scenario 4 beta(tau) linear and then constant ###
### beta(tau) = 4/3*x if x < 0.5, 2/3 if x >= 0.5
betatau.true <- function(x){
  re = rep(NA, length(x))
  index = x < 0.5
  re[index] = (4/3)*x[index]
  re[!index] = (2/3)
  return(re)
}
# compute int beta(tau)Qt(tau)dtau
intergral <- function(x, L, theta.vec, alpha, shape){
  re.beta = betatau.true(x)
  re.quan = as.numeric(quan.fun(x, L = L,
                                theta.vec = theta.vec,
                                alpha = alpha,
                                basis.fun = "Gamma",
                                shape = shape))
  re = re.beta*re.quan
  return(re)
}
integration.vec = rep(NA, ntimes)
for(i in 1:ntimes){
  integration.vec[i] <- integrate(intergral, lower = 0, upper = 1,
                                  L = L, theta.vec = dat.exp$theta.mat.true[,i],
                                  alpha = dat.exp$alph.vec.true[i],
                                  shape = shape)$value
  if(i%%100 == 0){cat(i, ",")}
}
eta.vec.true <- beta0 + integration.vec
q.nb = 1/(1 + exp(eta.vec.true))
set.seed(1234)
y.sim = matrix(NA, ncol = nsims, nrow = ntimes)
for(i.sim  in 1:nsims){
  y.sim[,i.sim] <- apply(cbind(rep(xi, ntimes), q.nb), 1,
                         function(x) rnbinom(1, size = x[1], prob = x[2]))
}
dat.all[["health4"]] <- y.sim  # each column is a simulation


### Scenario 5: decreasing ####
### beta(tau) exp(-tau^2/0.328)
betatau.true <- function(x, d = 0.328){
  exp(-x^2/d)
}
# compute int beta(tau)Qt(tau)dtau
intergral <- function(x, L, theta.vec, alpha, shape){
  re.beta = betatau.true(x)
  re.quan = as.numeric(quan.fun(x, L = L,
                                theta.vec = theta.vec,
                                alpha = alpha,
                                basis.fun = "Gamma",
                                shape = shape))
  re = re.beta*re.quan
  return(re)
}
integration.vec = rep(NA, ntimes)
for(i in 1:ntimes){
  integration.vec[i] <- integrate(intergral, lower = 0, upper = 1,
                                  L = L, theta.vec = dat.exp$theta.mat.true[,i],
                                  alpha = dat.exp$alph.vec.true[i],
                                  shape = shape)$value
  if(i%%100 == 0){cat(i, ",")}
}
eta.vec.true <- beta0 + integration.vec
q.nb = 1/(1 + exp(eta.vec.true))
set.seed(1234)
y.sim = matrix(NA, ncol = nsims, nrow = ntimes)
for(i.sim  in 1:nsims){
  y.sim[,i.sim] <- apply(cbind(rep(xi, ntimes), q.nb), 1,
                         function(x) rnbinom(1, size = x[1], prob = x[2]))
}
dat.all[["health5"]] <- y.sim  # each column is a simulation


### Scenario 7: decreasing ####
### beta(tau)=-tau+1
betatau.true <- function(x){
  -1*x+1
}
# compute int beta(tau)Qt(tau)dtau
intergral <- function(x, L, theta.vec, alpha, shape){
  re.beta = betatau.true(x)
  re.quan = as.numeric(quan.fun(x, L = L,
                                theta.vec = theta.vec,
                                alpha = alpha,
                                basis.fun = "Gamma",
                                shape = shape))
  re = re.beta*re.quan
  return(re)
}
integration.vec = rep(NA, ntimes)
for(i in 1:ntimes){
  integration.vec[i] <- integrate(intergral, lower = 0, upper = 1,
                                  L = L, theta.vec = dat.exp$theta.mat.true[,i],
                                  alpha = dat.exp$alph.vec.true[i],
                                  shape = shape)$value
  if(i%%100 == 0){cat(i, ",")}
}
eta.vec.true <- beta0 + integration.vec
q.nb = 1/(1 + exp(eta.vec.true))
set.seed(1234)
y.sim = matrix(NA, ncol = nsims, nrow = ntimes)
for(i.sim  in 1:nsims){
  y.sim[,i.sim] <- apply(cbind(rep(xi, ntimes), q.nb), 1,
                         function(x) rnbinom(1, size = x[1], prob = x[2]))
}
dat.all[["health6"]] <- y.sim  # each column is a simulation
save(dat.all, file = paste0("data_sim.rda"))
