########################################################################
# --------------- R codes used in the simulation study --------------- #
# NOTE:
# Summarize results when estiamted quantile functions are used
# please make sure the working directory is the folder where 
# R script named "FUNs.R" and simulated data are located at 
########################################################################

args <-  commandArgs(trailingOnly = TRUE)
sce.index = eval( parse(text=args[1]) )

# sce.index = 1, ..., 6

cat("THE current scenario is", sce.index)

# Load packages 
library(MASS)
library(igraph); library(mvtnorm); library(Matrix)
library(spam); library(Rcpp)
library(BayesLogit) # rpg, draw samples from PG 
library(truncnorm)

# source self-defined R functions 
sourceCpp("FUNs_quan_Rcpp.cpp")
source("FUNs.R")

# read in simulated data
load("data_sim.rda")

betatau.true.list <- list()
betatau.true.list[["health1"]] <- function(x){ rep(0.5, length(x)) }
betatau.true.list[["health2"]] <- function(x){ x }
betatau.true.list[["health3"]] <- function(x){ 1.5*x^2 }
betatau.true.list[["health4"]] <- function(x){
  re = rep(NA, length(x))
  index = x < 0.5
  re[index] = (4/3)*x[index]
  re[!index] = (2/3)
  return(re)
}
betatau.true.list[["health5"]] <- function(x){
  re = exp(-x^2/0.328)
  return(re)
}
betatau.true.list[["health6"]] <- function(x){
  re = -x+1
  return(re)
}

beta.true.list <- NULL
count = 1
for(n in 2){
  beta.true.mat <- matrix(NA, ncol = n + 1, nrow = 6)
  for(i in 1:6){
    beta.true.mat[i, ] <- 
      sapply(0:n, get_beta_bernstein, n = n, 
             f = betatau.true.list[[i]], simplify = T)
  }
  beta.true.list[[count]] <- beta.true.mat
  count = count + 1
}

beta.int.true.mat <- matrix(NA, nrow = 1, ncol = 6)
for(i in 1:6){
  beta.int.true.mat[1, i] <- 
    get_int_beta_bernstein(n = 2, coef = beta.true.list[[1]][i, ])
}

tab.sim.list <- list()
for(i in 1){
  tab.sim <- t(data.frame(beta0 = rep(-3.5, 6),
                          beta.true.list[[i]],
                          xi = rep(20, 6)))
  rownames(tab.sim) <- c(paste0("beta", 0:(i+1+1)), "xi")
  colnames(tab.sim) <- paste0("S", 1:6)
  tab.sim.list[[i]] <- tab.sim
}

# compute int beta(tau)Qt(tau)dtau
L = 4; shape = 5
re1 <- dat.all$exp
ntimes = nrow(dat.all$exp$x.sim.mat)

load("./inter_res/eta_true.rda") # eta.true.mat


### load in simulation results and store them into three lists, 
### each element contains results from ONE simulation; 
### these three lists have length of 100 ###
nsims = 100
n.fit = 2

re.all <- list()
for(i.sim in 1:nsims){
  load(paste0('./inter_res/Res_QF_errors_S', 
              sce.index, "_", i.sim, ".rda"))
  re.all[[i.sim]] <- re.fit.quan.errors
  if(i.sim%%10 == 0){cat(i.sim)}
}
cat("simulation results have been loaded")

## A function to get the matrix containing kernel functions
get_bernstein_mat <- function(tau, n){
  X <- NULL
  for(j in 0:n){
    X <- cbind(X, bernstein_orth(tau = tau, j = j, n = n))
  }
  return(X)
}

### 6:96 ###
get_res <- function(re, tau.vec, n, burn_in){
  
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

##################################################################
#### initialize a list to save all results ####
tab.list <- list()
### summarize results across simulations ###
beta.tau.est.list <- list()
beta.int.est.list <- list()
for(ii in 1:nsims){
  res.ii = get_res(re = re.all[[ii]], tau.vec = seq(0, 1, 0.01),
                   n = n.fit, burn_in = 2500)
  ## integrated beta(tau) ##
  beta.int.est.list[[ii]] <- res.ii$beta.int.est.mat
  beta.tau.est.list[[ii]] <- res.ii$beta.tau.est.mat
  if(ii %% 10){cat(ii, ",")}
}
save(beta.tau.est.list,
     file = paste0("./inter_res/Sum_betatau_QFerrors_S",
                   sce.index, ".rda"))

##### 1. get int_beta(tau) in the proposed and regular models #####
#### int_beta(tau)dtau: bias, RMSE, coverage #####
beta.int.est.mat <- do.call(rbind.data.frame, beta.int.est.list)
sum_res_intbeta <- function(beta.int.est.mat, beta.int.true){
  bias = mean(beta.int.est.mat$est - beta.int.true)
  bias.rel = mean((beta.int.est.mat$est - beta.int.true)/beta.int.true)
  RMSE = sqrt(mean((beta.int.est.mat$est - beta.int.true)^2))
  MSE = mean((beta.int.est.mat$est - beta.int.true)^2)
  coverage = mean(beta.int.est.mat$lower <= beta.int.true &
                    beta.int.true <= beta.int.est.mat$upper)
  return(data.frame(est.avg = mean(beta.int.est.mat$est),
                    sd = sd(beta.int.est.mat$est),
                    se.avg = mean(beta.int.est.mat$sd),
                    bias = bias, bias.rel = bias.rel,
                    RMSE = RMSE, MSE = MSE,
                    coverage = coverage))
}
tab.betaint = sum_res_intbeta(beta.int.est.mat = beta.int.est.mat, 
                              beta.int.true = 0.5)
rownames(tab.betaint) <- c("quan")
tab.list[["tab.betaint"]] <- tab.betaint


#### 2. compute int_beta(tau)Qt(tau) ######
get_covariatevals_eta_mu_X <- function(re, n.fit, BK.mat,
                                       burn_in, model = "quan"){
  
  if(model == "quan"){
    theta.mat.list <- re$theta.mat.list
    get_X <- function(theta.alpha.mat, BK.mat){
      return(theta.alpha.mat%*%BK.mat)
    }
    X.design.list <- lapply(theta.mat.list, get_X, BK.mat = BK.mat)
    beta.post <- re$beta[-c(1:burn_in), 2:(n.fit+1+1)]
    parm.post <- re$parm[-c(1:burn_in), 1]
    pred.post <- matrix(NA, nrow = nrow(X.design.list[[1]]), 
                        ncol = nrow(beta.post))
    for(i in 1:nrow(beta.post)){
      pred.post[,i] <- X.design.list[[i]] %*% matrix(beta.post[i, ], ncol = 1)
    }
    tab = data.frame(est = rowMeans(pred.post),
                     sd = apply(pred.post, 1, sd),
                     lower = apply(pred.post, 1, quantile, 0.025),
                     upper = apply(pred.post, 1, quantile, 0.975))
    eta.post <- t(apply(pred.post, 1, 
                        function(x) x + re$beta[-c(1:burn_in), 1]))
    tab.eta = data.frame(est = rowMeans(eta.post),
                         sd = apply(eta.post, 1, sd),
                         lower = apply(eta.post, 1, quantile, 0.025),
                         upper = apply(eta.post, 1, quantile, 0.975))
    mu.post = t(apply(eta.post, 1, 
                      function(x) x + log(parm.post)))
    tab.mu = data.frame(est = rowMeans(mu.post),
                        sd = apply(mu.post, 1, sd),
                        lower = apply(mu.post, 1, quantile, 0.025),
                        upper = apply(mu.post, 1, quantile, 0.975))
    return(list(tab.vals = tab, 
                tab.eta = tab.eta,
                tab.mu = tab.mu,
                pred.post = pred.post))
  }else if(model == "mean"){
    beta.post <- re$beta[-c(1:burn_in), 2]
    parm.post <- re$parm[-c(1:burn_in), 1]
    pred.post <- do.call(cbind, lapply(1:length(beta.post), function(x)
      beta.post[x]*re$X.bar.post[x, ]))
    tab = data.frame(est = rowMeans(pred.post),
                     sd = apply(pred.post, 1, sd),
                     lower = apply(pred.post, 1, quantile, 0.025),
                     upper = apply(pred.post, 1, quantile, 0.975))
    eta.post <- t(apply(pred.post, 1,
                        function(x) x + re$beta[-c(1:burn_in), 1]))
    tab.eta = data.frame(est = rowMeans(eta.post),
                         sd = apply(eta.post, 1, sd),
                         lower = apply(eta.post, 1, quantile, 0.025),
                         upper = apply(eta.post, 1, quantile, 0.975))
    
    mu.post = t(apply(eta.post, 1,
                      function(x) x + log(parm.post)))
    tab.mu = data.frame(est = rowMeans(mu.post),
                        sd = apply(mu.post, 1, sd),
                        lower = apply(mu.post, 1, quantile, 0.025),
                        upper = apply(mu.post, 1, quantile, 0.975))
    
    return(list(tab.vals = tab,
                tab.eta = tab.eta,
                tab.mu = tab.mu,
                pred.post = pred.post))
  }
  
}

BK.mat = get_integration_bernstein(L = L, n = n.fit,
                                   basis.fun = "Gamma",
                                   shape = shape)
pred.vals <- pred.eta <- pred.mu <- list()
for(i in 1:nsims){
  inter <- get_covariatevals_eta_mu_X(re = re.all[[i]],
                                      n.fit = 2, 
                                      BK.mat = BK.mat,
                                      burn_in = 2500,
                                      model = "quan")
  pred.vals[[i]] <- inter$tab.vals
  pred.eta[[i]] <- inter$tab.eta
  pred.mu[[i]] <- inter$tab.mu
  cat(i, ",")
}

sum_covariatevals <- function(tab, covariate.true){
  bias = (tab$est - covariate.true)
  bias.rel = (tab$est - covariate.true)/covariate.true
  difsquare = (tab$est - covariate.true)^2
  coverage.ind = (tab$lower <= covariate.true) &
    (covariate.true <= tab$upper)
  return(list(bias = bias,
              bias.rel = bias.rel,
              difsquare = difsquare, 
              coverage.ind = coverage.ind))
}

results.vals <- results.eta <- results.mu <- list()
for(i in 1:nsims){
  results.vals[[i]] <- sum_covariatevals(tab = pred.vals[[i]],
                                         covariate.true = 
                                           eta.true.mat[,sce.index] + 3.5)
  results.eta[[i]] <- sum_covariatevals(tab = pred.eta[[i]],
                                        covariate.true = 
                                          eta.true.mat[,sce.index])
  results.mu[[i]] <- sum_covariatevals(tab = pred.mu[[i]],
                                       covariate.true = 
                                         eta.true.mat[,sce.index] + log(20))
}

bias.mat <- do.call(rbind, lapply(results.vals, function(x) x$bias))
tab.list[["bias.mat.vals"]] <- bias.mat

difsquare.cov <- do.call(c, lapply(results.vals, 
                                   function(x) x$difsquare))
bias.cov = mean(do.call(c, lapply(results.vals, 
                                  function(x) x$bias)))
bias.rel.cov = mean(do.call(c, lapply(results.vals, 
                                      function(x) x$bias.rel)))
RMSE.cov = sqrt( mean( difsquare.cov ) )
MSE.cov = mean( difsquare.cov )
coverage.cov = mean( rowMeans(do.call(cbind, 
                                      lapply(results.vals, 
                                             function(x) x$coverage.ind))) )
tab.covariatevals = data.frame(avg.true =
                                 rep(mean(eta.true.mat[,sce.index] + 3.5), 1),
                               bias = c(bias.cov),
                               bias.rel = c(bias.rel.cov),
                               RMSE = c(RMSE.cov),
                               MSE = c(MSE.cov),
                               coverage = c(coverage.cov))
rownames(tab.covariatevals) <- c("quan")
tab.list[["tab.covariatevals"]] <- tab.covariatevals
cat("complete the computation of predicted covariates")


######## predicted eta ########
bias.mat <- do.call(rbind, lapply(results.eta, function(x) x$bias))
tab.list[["bias.mat.eta"]] <- bias.mat
difsquare.eta <- do.call(c, lapply(results.eta, 
                                   function(x) x$difsquare))
bias.eta = mean(do.call(c, lapply(results.eta, 
                                  function(x) x$bias)))
bias.rel.eta = mean(do.call(c, lapply(results.eta, 
                                      function(x) x$bias.rel)))
RMSE.eta = sqrt( mean( difsquare.eta ) )
MSE.eta =  mean( difsquare.eta ) 
coverage.eta = mean( rowMeans(do.call(cbind, 
                                      lapply(results.eta, 
                                             function(x) x$coverage.ind))) )

tab.eta = data.frame(avg.true =
                       rep(mean(eta.true.mat[,sce.index]), 1),
                     bias = c(bias.eta),
                     bias.rel = c(bias.rel.eta),
                     RMSE = c(RMSE.eta),
                     MSE = c(MSE.eta),
                     coverage = c(coverage.eta))
rownames(tab.eta) <- c("quan")
tab.list[["tab.eta"]] <- tab.eta

######## predicted mu ########
bias.mat <- do.call(rbind, lapply(results.mu, function(x) x$bias))
tab.list[["bias.mat.mu"]] <- bias.mat
difsquare.mu <- do.call(c, lapply(results.mu, 
                                  function(x) x$difsquare))
bias.mu = mean(do.call(c, lapply(results.mu, 
                                 function(x) x$bias)))
bias.rel.mu = mean(do.call(c, lapply(results.mu, 
                                     function(x) x$bias.rel)))
RMSE.mu = sqrt( mean( difsquare.mu ) )
MSE.mu = mean( difsquare.mu )
coverage.mu = mean( rowMeans(do.call(cbind, 
                                     lapply(results.mu, 
                                            function(x) x$coverage.ind))) )

tab.mu = data.frame(avg.true =
                      rep(mean(eta.true.mat[,sce.index]+log(20)), 1),
                    bias = c(bias.mu),
                    bias.rel = c(bias.rel.mu),
                    RMSE = c(RMSE.mu),
                    MSE = c(MSE.mu),
                    coverage = c(coverage.mu))
rownames(tab.mu) <- c("quan")
tab.list[["tab.mu"]] <- tab.mu


######### 3. compute attributable ED counts #####
get_postED_X <- function(re, burn_in, BK.mat, n.fit, model = "quan"){
  if(model == "quan"){
    get_X <- function(theta.alpha.mat, BK.mat, n.fit){
      X <- cbind(1, theta.alpha.mat%*%BK.mat)
      X.null <- X
      X.null[,2:(n.fit+1+1)] <- 0
      return(list(X = X, X.null = X.null))
    }
    inter.list <- lapply(re$theta.mat.list, get_X, 
                         BK.mat = BK.mat, n.fit = n.fit)
    X.design.list <- lapply(inter.list, function(x) x$X)
    X.null.list <- lapply(inter.list, function(x) x$X.null)
    beta.post = re$beta[-c(1:burn_in), ]
    ED.post <- ED.post.null <- matrix(NA, nrow = nrow(re$theta.mat.list[[1]]),
                                      ncol = nrow(beta.post))
    for(i in 1:nrow(beta.post)){
      beta.vec.i = matrix(beta.post[i, ], ncol = 1)
      ED.post[,i] <-  X.design.list[[i]] %*% beta.vec.i
      ED.post.null[,i] <- X.null.list[[i]] %*% beta.vec.i
    }
    ED.attri <- exp(ED.post) - exp(ED.post.null)
    xi.post = re$parm[-c(1:burn_in),  1]
    ED.attri.all.post <- rowSums(apply(ED.attri, 1, function(x) x*xi.post))
    ED.attri.mat <- data.frame(est = mean(ED.attri.all.post),
                               sd = sd(ED.attri.all.post),
                               lci = quantile(ED.attri.all.post, 0.025),
                               uci = quantile(ED.attri.all.post, 0.975))
    return(ED.attri.mat)
  }else if(model == "mean"){
    get_X <- function(X.bar){
      X <- cbind(1, X.bar)
      X.null <- X
      X.null[,2] <- 0
      return(list(X = X, X.null = X.null))
    }
    inter.list <- lapply(1:nrow(re$X.bar.post), function(x)
      get_X(X.bar = re$X.bar.post[x, ]))
    X.design.list <- lapply(inter.list, function(x) x$X)
    X.null.list <- lapply(inter.list, function(x) x$X.null)
    beta.post = re$beta[-c(1:burn_in), ]
    ED.post <- ED.post.null <- matrix(NA, nrow = ncol(re$X.bar.post),
                                      ncol = nrow(beta.post))
    for(i in 1:nrow(beta.post)){
      beta.vec.i = matrix(beta.post[i, ], ncol = 1)
      ED.post[,i] <-  X.design.list[[i]] %*% beta.vec.i
      ED.post.null[,i] <- X.null.list[[i]] %*% beta.vec.i
    }
    ED.attri <- exp(ED.post) - exp(ED.post.null)
    xi.post = re$parm[-c(1:burn_in),  1]
    ED.attri.all.post <- rowSums(apply(ED.attri, 1, function(x) x*xi.post))
    ED.attri.mat <- data.frame(est = mean(ED.attri.all.post),
                               sd = sd(ED.attri.all.post),
                               lci = quantile(ED.attri.all.post, 0.025),
                               uci = quantile(ED.attri.all.post, 0.975))
    return(ED.attri.mat)
  }
  
}
ED.mat <- data.frame(est = rep(NA, nsims), sd = rep(NA, nsims),
             lower = rep(NA, nsims), upper = rep(NA, nsims))
for(ii in 1:nsims){
  ED.mat[ii, ] <- get_postED_X(re = re.all[[ii]],
                               burn_in = 2500,
                               BK.mat = BK.mat, n.fit = n.fit,
                               model = "quan")
  if(ii %% 10){cat(ii, ",")}
}
tab.ED.list <- list(
  quan = ED.mat
)
ED.true <- sum(20*(exp(eta.true.mat[,sce.index]) - exp(-3.5)))
tab.ED.list[["true"]] <- ED.true
tab.list[["tab.ED.list"]] <- tab.ED.list
cat("complete the computation of attributable ED visits")


###### 4. get beta(tau) ######
tau.vec = seq(0, 1, 0.01)
tau.vec.1 = seq(0.05, 0.95, 0.01)  # 6:96
tau.vec.2 = seq(0.01, 0.99, 0.01)  # 2:100

betatau.true.val = betatau.true.list[[sce.index]](tau.vec)
betatau.true.val.1 = betatau.true.list[[sce.index]](tau.vec.1)
betatau.true.val.2 = betatau.true.list[[sce.index]](tau.vec.2)

beta.tau.avg <- Reduce("+", beta.tau.est.list)/length(beta.tau.est.list)

bias.rel.betatau = mean(do.call(c, 
                                lapply(beta.tau.est.list, 
                                       function(x) (x$est - betatau.true.val)/betatau.true.val)))
RMSE.betatau = sqrt(mean(do.call(c, 
                                 lapply(beta.tau.est.list, 
                                        function(x) (x$est - betatau.true.val)^2))))
MSE.betatau = RMSE.betatau^2



bias.rel.betatau.1 = mean(do.call(c, 
                                  lapply(beta.tau.est.list, 
                                         function(x) (x$est[6:96] - betatau.true.val.1)/betatau.true.val.1)))
RMSE.betatau.1 = sqrt(mean(do.call(c, 
                                   lapply(beta.tau.est.list, 
                                          function(x) (x$est[6:96] - betatau.true.val.1)^2))))
MSE.betatau.1 = RMSE.betatau.1^2


bias.rel.betatau.2 = mean(do.call(c, 
                                  lapply(beta.tau.est.list, 
                                         function(x) (x$est[2:100] -betatau.true.val.2)/betatau.true.val.2)))
RMSE.betatau.2 = sqrt(mean(do.call(c, 
                                   lapply(beta.tau.est.list, 
                                          function(x) (x$est[2:100] -betatau.true.val.2)^2))))
MSE.betatau.2 = RMSE.betatau.2^2

coverage.betatau <- mean(rowMeans(do.call(cbind, lapply(beta.tau.est.list,
                                                        function(x) x$lower <= betatau.true.val &
                                                          betatau.true.val <= x$upper))))

coverage.betatau.1 <- mean(rowMeans(do.call(cbind, lapply(beta.tau.est.list,
                                                          function(x) x$lower[6:96] <= betatau.true.val.1 &
                                                            betatau.true.val.1 <= x$upper[6:96]))))

coverage.betatau.2 <- mean(rowMeans(do.call(cbind, lapply(beta.tau.est.list,
                                                          function(x) x$lower[2:100] <= betatau.true.val.2 &
                                                            betatau.true.val.2 <= x$upper[2:100]))))

tab.betatau = data.frame(bias.rel = c(bias.rel.betatau, bias.rel.betatau.1,
                                      bias.rel.betatau.2),
                         RMSE = c(RMSE.betatau, RMSE.betatau.1, RMSE.betatau.2),
                         MSE = c(MSE.betatau, MSE.betatau.1, MSE.betatau.2),
                         coverage = c(coverage.betatau, coverage.betatau.1,
                                      coverage.betatau.2))
rownames(tab.betatau) <- c("quan0_1", "quan5_95", "quan1_99")
tab.list[["tab.betatau"]] <- tab.betatau
cat("COMPLETE computations of beta(tau)")


##### 6. get estimated betas #######
get_res_bets <- function(re, burn_in){
  beta.post <- re$beta[-c(1:burn_in), ]
  parm.post <- re$parm[-c(1:burn_in), ]
  parm.post <- cbind(beta.post, parm.post, beta.post[,1] + log(parm.post))
  tab = data.frame(est = colMeans(parm.post),
                   sd = apply(parm.post, 2, sd),
                   lower = apply(parm.post, 2, quantile, 0.025),
                   upper = apply(parm.post, 2, quantile, 0.975))
  return(tab)
}
burn_in = 2500
re.inter <- lapply(re.all, function(x) get_res_bets(re = x, burn_in = burn_in))
est.mat.inter <- do.call(rbind, lapply(re.inter, function(x) x$est))
sd.mat.inter <- do.call(rbind, lapply(re.inter, function(x) x$sd))
lower.mat.inter <- do.call(rbind, lapply(re.inter, function(x) x$lower))
upper.mat.inter <- do.call(rbind, lapply(re.inter, function(x) x$upper))
sum_res <- function(est.mat, sd.mat, 
                    lower.mat, upper.mat, true.vec){
  coverage.vec <- rep(NA, ncol(est.mat))
  for(i in 1:ncol(est.mat)){
    coverage.vec[i] <- mean(lower.mat[,i] <= true.vec[i] &
                              true.vec[i] <= upper.mat[,i])
  }
  re = data.frame(true = true.vec,
                  est.avg = colMeans(est.mat),
                  sd = apply(est.mat, 2, sd),
                  se.avg = colMeans(sd.mat),
                  coverage = coverage.vec)
}
tab.basiscoef <- sum_res(est.mat = est.mat.inter,
                         sd.mat = sd.mat.inter,
                         lower.mat = lower.mat.inter,
                         upper.mat = upper.mat.inter,
                         true.vec = c(tab.sim.list[[n.fit-1]][,sce.index],
                                      tab.sim.list[[n.fit-1]][1,sce.index] + 
                                        log(tab.sim.list[[n.fit-1]]["xi",sce.index])))
tab.list[["tab.basiscoef"]] <- tab.basiscoef

save(tab.list, file = paste0("./inter_res/Sum_res_errors_S", 
                             sce.index, ".rda"))




