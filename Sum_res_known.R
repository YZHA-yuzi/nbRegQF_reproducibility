########################################################################
# --------------- R codes used in the simulation study --------------- #
# NOTE:
# Summarize results
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

# start = proc.time()[3]

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

### get eta for all scenarios ####
intergral <- function(x, f.betatau, L, theta.vec, alpha, shape){
  re.beta = f.betatau(x)
  re.quan = as.numeric(quan.fun(x, L = L,
                                theta.vec = theta.vec,
                                alpha = alpha,
                                basis.fun = "Gamma",
                                shape = shape))
  re = re.beta*re.quan
  return(re)
}

# compute int beta(tau)Qt(tau)dtau
L = 4; shape = 5
re1 <- dat.all$exp
ntimes = nrow(dat.all$exp$x.sim.mat)

#### ONLY NEED TO RUN THE FOLLOWING CHUNK OF CODES ONCE ####
# eta.true.mat <- matrix(NA, nrow = ntimes, ncol = 6)
# start = proc.time()[3]
# for(ii in 1:6){
#   integration.vec.ii = rep(NA, ntimes)
#   for(i in 1:ntimes){
#     integration.vec.ii[i] <- integrate(intergral, lower = 0, upper = 1,
#                                        f.betatau = betatau.true.list[[ii]],
#                                        L = L, 
#                                        theta.vec = re1$theta.mat.true[,i],
#                                        alpha = re1$alph.vec.true[i],
#                                        shape = shape)$value
#   }
#   eta.true.mat[,ii] <- tab.sim.list[[1]][1, ii] + integration.vec.ii
#   cat(ii, ",")
# }
# cat("complete the computation of true eta")
# proc.time()[3] - start
### !!!! NOTES:
### It takes about 2min 
### on a computer with 2.9 GHz 6-Core Intel Core i9 and 32 Gb memory
### So, save this matrix as an intermediate results
# save(eta.true.mat, file = "./inter_res/eta_true.rda")

load("./inter_res/eta_true.rda") # eta.true.mat


### load in simulation results and store them into TWO lists, 
### each element contains results from ONE simulation; 
### these TWO lists have length of 100 ###
nsims = 100

# re.all <- re.all.mean <- list()
# for(i.sim in 1:nsims){
#   
#   load(paste0("./inter_res/Res_QF_known_S", ss, "_", i.sim, ".rda"))
#   re.all[[i.sim]] <- re.fit.quan
#   
#   load(paste0("./inter_res/Res_Mean_known_S", ss, "_", i.sim, ".rda"))
#   re.all.mean[[i.sim]] <- re.fit.mean
# 
#   save(re.all,
#        file = paste0("./inter_res/Res_QF_known_S", sce.index, "_all.rda"))
#   save(re.all.mean,
#        file = paste0("./inter_res/Res_Mean_known_S", sce.index, "_all.rda"))
# 
# }

load(paste0("./inter_res/Res_QF_known_S", sce.index, "_all.rda"))
load(paste0("./inter_res/Res_Mean_known_S", sce.index, "_all.rda"))


##################################################################
n.fit = 2
#### initialize a list to save all results ####
tab.list <- list()
### summarize results across simulations ###
beta.tau.est.list <- list()
beta.int.est.list <- list()
for(ii in 1:nsims){
  res.ii = get_res_QFknown(re = re.all[[ii]], 
                           tau.vec = seq(0, 1, 0.01),
                           n = n.fit, burn_in = 2500)
  ## integrated beta(tau) ##
  beta.int.est.list[[ii]] <- res.ii$beta.int.est.mat
  beta.tau.est.list[[ii]] <- res.ii$beta.tau.est.mat
}
save(beta.tau.est.list,
     file = paste0("./inter_res/Sum_betatau_QFknown_S",
                   sce.index, ".rda"))

##### 1. get int_beta(tau) in the proposed and regular models #####
#### int_beta(tau)dtau: bias, RMSE, coverage #####
beta.int.est.mat <- do.call(rbind.data.frame, beta.int.est.list)
beta1.mean.mat <- do.call(rbind.data.frame,
                          lapply(lapply(re.all.mean,
                                        get_res_mean, burn_in = 2500),
                                 function(x) x[2, ]))
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
                    MSE = MSE, RMSE = RMSE,
                    coverage = coverage))
}
tab.betaint = rbind.data.frame(
  sum_res_intbeta(beta.int.est.mat = beta.int.est.mat,
                  beta.int.true = 0.5),
  sum_res_intbeta(beta1.mean.mat, beta.int.true = 0.5))
rownames(tab.betaint) <- c("quan", "mean")
tab.list[["tab.betaint"]] <- tab.betaint


##### WAIC comparison ######
WAIC.mat.comp <- data.frame(quan = sapply(re.all, function(x) x$WAIC),
                            mean = sapply(re.all.mean, function(x) x$WAIC))
colnames(WAIC.mat.comp) <- c("quan", "mean")
tab.list[["WAIC.mat.comp"]] <- WAIC.mat.comp


#### 2. compute int_beta(tau)Qt(tau) ######
theta.alpha.mat <- t(rbind(re1$alph.vec.true, re1$theta.mat.true))
BK.mat = get_integration_bernstein(L = L, n = n.fit,
                                   basis.fun = "Gamma",
                                   shape = shape)
X.design = cbind(1, theta.alpha.mat%*%BK.mat)

mean.gamma.vec = rep(NA, ntimes)
for(i in 1:ntimes){
  mean.gamma.vec[i] <- integrate(quan.fun, L = 4,
                                 theta.vec = re1$theta.mat.true[,i],
                                 alpha = re1$alph.vec.true[i],
                                 basis.fun = "Gamma",
                                 shape = shape,
                                 lower = 0, upper = 1)$value
}
X.design.mean <- as.matrix(cbind(1, mean.gamma.vec))


get_covariatevals_eta_mu <- function(re, n.fit, X, model, burn_in){
  if(model == "mean"){
    beta.post <- re$beta[-c(1:burn_in), 2]
    parm.post <- re$parm[-c(1:burn_in), 1]
    pred.post <- matrix(X, ncol = 1) %*% matrix(beta.post, nrow = 1)
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
                tab.mu = tab.mu))

  }else if(model == "quan"){
    beta.post <- re$beta[-c(1:burn_in), 2:(n.fit+1+1)]
    parm.post <- re$parm[-c(1:burn_in), 1]
    pred.post <- X %*% t(beta.post)
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
                tab.mu = tab.mu))
  }
}

pred.vals <- pred.vals.mean <- list()
pred.eta <- pred.eta.mean <- list()
pred.mu <- pred.mu.mean <- list()
for(i in 1:nsims){
  burn_in = 2500
  inter <- get_covariatevals_eta_mu(re = re.all[[i]],
                                    n.fit = 2, X = X.design[,-1],
                                    model = "quan", burn_in = burn_in)
  pred.vals[[i]] <- inter$tab.vals
  pred.eta[[i]] <- inter$tab.eta
  pred.mu[[i]] <- inter$tab.mu

  inter.mean <- get_covariatevals_eta_mu(re = re.all.mean[[i]],
                                         n.fit = 2, X = mean.gamma.vec,
                                         model = "mean", burn_in = burn_in)
  pred.vals.mean[[i]] <- inter.mean$tab.vals
  pred.eta.mean[[i]] <- inter.mean$tab.eta
  pred.mu.mean[[i]] <- inter.mean$tab.mu
  if(i%%10 == 0){cat(i)}
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

results.vals <- results.vals.mean <- list()
results.eta <- results.eta.mean <- list()
results.mu <- results.mu.mean <- list()
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

  results.vals.mean[[i]] <- sum_covariatevals(tab = pred.vals.mean[[i]],
                                              covariate.true =
                                                eta.true.mat[,sce.index] + 3.5)
  results.eta.mean[[i]] <- sum_covariatevals(tab = pred.eta.mean[[i]],
                                             covariate.true =
                                               eta.true.mat[,sce.index])
  results.mu.mean[[i]] <- sum_covariatevals(tab = pred.mu.mean[[i]],
                                            covariate.true =
                                              eta.true.mat[,sce.index] + log(20))
}


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

difsquare.cov.mean <- do.call(c, lapply(results.vals.mean,
                                        function(x) x$difsquare))
bias.cov.mean =  mean(do.call(c, lapply(results.vals.mean,
                                        function(x) x$bias)))
bias.rel.cov.mean =  mean(do.call(c, lapply(results.vals.mean,
                                        function(x) x$bias.rel)))
RMSE.cov.mean = sqrt( mean(difsquare.cov.mean) )
MSE.cov.mean = mean(difsquare.cov.mean)
coverage.cov.mean = mean( rowMeans(do.call(cbind,
                                           lapply(results.vals.mean,
                                                  function(x) x$coverage.ind))) )
tab.covariatevals = data.frame(avg.true =
                                 rep(mean(eta.true.mat[,sce.index] + 3.5), 2),
                               bias = c(bias.cov,
                                        bias.cov.mean),
                               bias.rel = c(bias.rel.cov, bias.rel.cov.mean),
                               RMSE = c(RMSE.cov,
                                        RMSE.cov.mean),
                               MSE = c(MSE.cov, MSE.cov.mean),
                               coverage = c(coverage.cov,
                                            coverage.cov.mean))
rownames(tab.covariatevals) <- c("quan", "mean")
tab.list[["tab.covariatevals"]] <- tab.covariatevals
cat("complete the computation of predicted covariates")


######## predicted eta ########
difsquare.eta <- do.call(c, lapply(results.eta,
                                   function(x) x$difsquare))
bias.eta = mean(do.call(c, lapply(results.eta,
                                  function(x) x$bias)))
bias.rel.eta = mean(do.call(c, lapply(results.eta,
                                      function(x) x$bias.rel)))
RMSE.eta = sqrt( mean( difsquare.eta ) )
MSE.eta = mean( difsquare.eta )
coverage.eta = mean( rowMeans(do.call(cbind,
                                      lapply(results.eta,
                                             function(x) x$coverage.ind))) )



difsquare.eta.mean <- do.call(c, lapply(results.eta.mean,
                                        function(x) x$difsquare))
bias.eta.mean =  mean(do.call(c, lapply(results.eta.mean,
                                        function(x) x$bias)))
bias.rel.eta.mean = mean(do.call(c, lapply(results.eta.mean,
                                           function(x) x$bias.rel)))
RMSE.eta.mean = sqrt( mean(difsquare.eta.mean) )
MSE.eta.mean = mean(difsquare.eta.mean)
coverage.eta.mean = mean( rowMeans(do.call(cbind,
                                           lapply(results.eta.mean,
                                                  function(x) x$coverage.ind))) )
tab.eta = data.frame(avg.true =
                       rep(mean(eta.true.mat[,sce.index]), 2),
                     bias = c(bias.eta,
                              bias.eta.mean),
                     bias.rel = c(bias.rel.eta, bias.rel.eta.mean),
                     RMSE = c(RMSE.eta,
                              RMSE.eta.mean),
                     MSE = c(MSE.eta, MSE.eta.mean),
                     coverage = c(coverage.eta,
                                  coverage.eta.mean))
rownames(tab.eta) <- c("quan", "mean")
tab.list[["tab.eta"]] <- tab.eta

######## predicted mu ########
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

difsquare.mu.mean <- do.call(c, lapply(results.mu.mean,
                                       function(x) x$difsquare))
bias.mu.mean =  mean(do.call(c, lapply(results.mu.mean,
                                       function(x) x$bias)))
bias.rel.mu.mean = mean(do.call(c, lapply(results.mu.mean,
                                          function(x) x$bias.rel)))
RMSE.mu.mean = sqrt( mean(difsquare.mu.mean) )
MSE.mu.mean = mean(difsquare.mu.mean)
coverage.mu.mean = mean( rowMeans(do.call(cbind,
                                          lapply(results.mu.mean,
                                                 function(x) x$coverage.ind))) )
tab.mu = data.frame(avg.true =
                      rep(mean(eta.true.mat[,sce.index]+log(20)), 2),
                    bias = c(bias.mu,
                             bias.mu.mean),
                    bias.rel = c(bias.rel.mu, bias.rel.mu.mean),
                    RMSE = c(RMSE.mu,
                             RMSE.mu.mean),
                    MSE = c(MSE.mu, MSE.mu.mean),
                    coverage = c(coverage.mu,
                                 coverage.mu.mean))
rownames(tab.mu) <- c("quan", "mean")
tab.list[["tab.mu"]] <- tab.mu


######### 3. compute attributable ED counts #####
get_postED <- function(re, burn_in, X, X.null){
  beta.post = re$beta[-c(1:burn_in), ]
  ED.post = X %*% t(beta.post)
  ED.post.null = X.null %*% t(beta.post)
  ED.attri <- exp(ED.post) - exp(ED.post.null)
  # xi.post = re$parm[-c(1:burn_in),  "xi"]
  xi.post = re$parm[-c(1:burn_in),  1]
  ED.attri.all.post <- rowSums(apply(ED.attri, 1, function(x) x*xi.post))
  ED.attri.mat <- data.frame(est = mean(ED.attri.all.post),
                             sd = sd(ED.attri.all.post),
                             lci = quantile(ED.attri.all.post, 0.025),
                             uci = quantile(ED.attri.all.post, 0.975))
  return(ED.attri.mat)
}

ED.mat <- ED.mean.mat <-
  data.frame(est = rep(NA, nsims), sd = rep(NA, nsims),
             lower = rep(NA, nsims), upper = rep(NA, nsims))
for(ii in 1:nsims){

  X.null <- X.design
  X.null[,2:(n.fit+1+1)] <- 0

  # X.design.mean = as.matrix(cbind(1, mean.gamma.vec))
  X.null.mean <- X.design.mean
  X.null.mean[,2] <- 0

  burn_in = 2500

  ED.mat[ii, ] <- get_postED(re = re.all[[ii]],
                             burn_in = burn_in,
                             X = X.design,
                             X.null = X.null)
  ED.mean.mat[ii, ] <- get_postED(re = re.all.mean[[ii]],
                                  burn_in = burn_in,
                                  X = X.design.mean,
                                  X.null = X.null.mean)
  # if(ii %% 10){cat(ii, ",")}
}
tab.ED.list <- list(
  quan = ED.mat,
  mean = ED.mean.mat
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

bias.rel.betatau = mean(do.call(c, lapply(beta.tau.est.list,
                          function(x) (x$est - betatau.true.val)/betatau.true.val)))
RMSE.betatau = sqrt(mean(do.call(c,
                                 lapply(beta.tau.est.list,
                                        function(x) (x$est - betatau.true.val)^2))))
MSE.betatau = RMSE.betatau^2


bias.rel.betatau.1 = mean(do.call(c, lapply(beta.tau.est.list,
                                          function(x) (x$est[6:96] - betatau.true.val.1)/
                                            betatau.true.val.1)))
RMSE.betatau.1 = sqrt(mean(do.call(c,
                                   lapply(beta.tau.est.list,
                                          function(x) (x$est[6:96] - betatau.true.val.1)^2))))
MSE.betatau.1 = RMSE.betatau.1^2


bias.rel.betatau.2 = mean(do.call(c, lapply(beta.tau.est.list,
                                            function(x) (x$est[2:100] - betatau.true.val.2)/
                                              betatau.true.val.2)))
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

re.inter <- lapply(re.all.mean, function(x) get_res_bets(re = x, burn_in = burn_in))
est.mat.inter <- do.call(rbind, lapply(re.inter, function(x) x$est))
sd.mat.inter <- do.call(rbind, lapply(re.inter, function(x) x$sd))
lower.mat.inter <- do.call(rbind, lapply(re.inter, function(x) x$lower))
upper.mat.inter <- do.call(rbind, lapply(re.inter, function(x) x$upper))
sum_res_1 <- function(est.mat, sd.mat,
                      lower.mat, upper.mat){
  re = data.frame(est.avg = colMeans(est.mat),
                  sd = apply(est.mat, 2, sd),
                  se.avg = colMeans(sd.mat))
}
tab.coef.mean <- sum_res_1(est.mat = est.mat.inter,
                           sd.mat = sd.mat.inter,
                           lower.mat = lower.mat.inter,
                           upper.mat = upper.mat.inter)
tab.list[["tab.coef.mean"]] <- tab.coef.mean

save(tab.list, file = paste0("./inter_res/Sum_res_known_S",
                             sce.index, ".rda"))

# proc.time()[3] - start



