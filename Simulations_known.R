########################################################################
# --------------- R codes used in the simulation study --------------- #
# NOTE:
# fitting health models assuming quantile functions are known
# please make sure the working directory is the folder where 
# R script named "FUNs.R" and simulated data are located at 
########################################################################

args <-  commandArgs(trailingOnly = TRUE)
simindex = eval( parse(text=args[1]) )
sceindex = eval( parse(text=args[2]) )

## simindex = 1, ..., 100 (index of the simulation)
## sceindex = 1, ..., 6 (index of the simulation scenario)

# Load packages 
library(MASS)
library(igraph); library(mvtnorm); library(Matrix)
library(spam); library(Rcpp)
library(BayesLogit) # rpg, draw samples from PG 
library(truncnorm)

# source self-defined R functions 
sourceCpp("FUNs_quan_Rcpp.cpp")
source("FUNs.R")

dir.create("inter_res")

# read in simulated data
load("data_sim.rda")

ntimes = 1000; nsize = 100
nsims = 100; L = 4; shape = 5

set.seed(12345 + simindex)


re1 <- dat.all$exp
theta.mat = as.matrix(t(rbind(re1$alph.vec.true, 
                              re1$theta.mat.true)))

# prepare design matrix used in fitting "mean" models
### compute the averaged exposure at each time point ###
mean.gamma.vec = rep(NA, ntimes)
for(i in 1:ntimes){
  mean.gamma.vec[i] <- integrate(quan.fun, L = 4,
                                 theta.vec = re1$theta.mat.true[,i],
                                 alpha = re1$alph.vec.true[i],
                                 basis.fun = "Gamma",
                                 shape = shape,
                                 lower = 0, upper = 1)$value
}
X.design.mean = as.matrix(cbind(1, mean.gamma.vec))

sce.vec = 1:6
sce.i = sce.vec[sceindex]
cat("the current simulation is", simindex)
cat("the current sceanrio is", sce.i)

# get health data under the current simulation scenario
y.sim = dat.all[[paste0("health",sce.i)]]


## fit the proposed scalar-on-quantile models
BK.mat = get_integration_bernstein(L = L, n = 2,
                                   basis.fun = "Gamma",
                                   shape = shape)
X.design = cbind(1, theta.mat%*%BK.mat)
## ONE simulation takes about 20s 
## on a computer with 2.9 GHz 6-Core Intel Core i9 and 32 Gb memory
re.fit.quan <- fit.health.Bernstein(y = y.sim[,simindex],
                                    X.design = X.design,
                                    niter = 5000,
                                    ntimes = ntimes,
                                    burn_in = 2500)
save(re.fit.quan,
     file = paste0("./inter_res/Res_QF_known_S", sce.i, "_", simindex, ".rda"))

## fit the "mean" model 
## one simulation takes about 20s 
## on a computer with 2.9 GHz 6-Core Intel Core i9 and 32 Gb memory
re.fit.mean <- fit.health.Bernstein(y = y.sim[,simindex],
                                    X.design = X.design.mean,
                                    niter = 5000,
                                    ntimes = ntimes, 
                                    burn_in = 2500)
save(re.fit.mean,
     file = paste0("./inter_res/Res_Mean_known_S", sce.i, "_", simindex, ".rda"))
