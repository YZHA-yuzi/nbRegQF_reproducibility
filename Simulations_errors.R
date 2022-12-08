########################################################################
# --------------- R codes used in the simulation study --------------- #
# NOTE:
# fitting health models using estimated quantile functions
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
library(Rcpp)
library(BayesLogit) # rpg, draw samples from PG
library(truncnorm)
library(spam)

# sessionInfo()

# source self-defined R functions
sourceCpp("FUNs_quan_Rcpp.cpp")
source("FUNs.R")

# read in simulated data
load("data_sim.rda")

# read in prepared mu and covariance matrix of MVN priors
## using relative path
load("./inter_res/mu_cov_MVNprior.rda")

ntimes = 1000; nsize = 100
nsims = 100; L = 4; shape = 5

set.seed(12345 + simindex)

sce.vec = 1:6
sce.i = sce.vec[sceindex]
cat("the current simulation is", simindex)
cat("the current sceanrio is", sce.i)

theta.pri.mu = matrix(matrix_norm_post$M, ncol = 1)
theta.pri.pre = solve(matrix_norm_post$Sigma.bdiag)
theta.mat.ini = matrix(matrix_norm_post$M,
                       byrow = T, nrow = 1000, ncol = 5)

y.sim = dat.all[[paste0("health", sce.i)]]
BK.mat = get_integration_bernstein(L = L, n = 2,
                                   basis.fun = "Gamma",
                                   shape = shape)

## ONE simulation takes about 35s to run
## on a computer with 2.9 GHz 6-Core Intel Core i9 and 32 Gb memory
# start = proc.time()[3]
re.fit.quan.errors <- fit.health.Bernstein.var(y = y.sim[,simindex],
                                               BK.mat = BK.mat,
                                               theta.mat.ini = theta.mat.ini,
                                               theta.pri.mu = theta.pri.mu,
                                               theta.pri.pre = theta.pri.pre,
                                               niter = 5000,
                                               ntimes = ntimes,
                                               nbeta.input = 1+1+2,
                                               burn_in = 2500)
# proc.time()[3] - start
save(re.fit.quan.errors,
     file = paste0("./inter_res/Res_QF_errors_S", sce.i, "_", simindex, ".rda"))


