########################################################################
# R codes to generate the table containing simulation results (Table 1) #
# NOTE:
# Summarize results when estiamted quantile functions are used
# please make sure the working directory is the folder where 
# R script named "FUNs.R" and simulated data are located at 
########################################################################

library(ggplot2)
library(ggpubr)
library(ggsci)
library(sf)
library(tmap)

### true beta(tau) ###
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


tab.betaint.list <- tab.betatau.list <- 
  tab.covariatevals.list <-tab.WAIC.comp <- list()
tab.ED.list.known <- tab.ED.list.errors <- list()

betatau.bias.mat.errors <- betatau.bias.mat.known <- list()

sce.vec = 1:6
for(i in sce.vec){
  load(paste0("./inter_res/Sum_res_known_S", i, ".rda"))
  tab.list.known <- tab.list
  
  load(paste0("./inter_res/Sum_res_errors_S", i, ".rda"))
  tab.list.errors <- tab.list
  
  beta.tau.true.vec <- betatau.true.list[[i]](seq(0, 1, 0.01))
  
  load(paste0("./inter_res/Sum_betatau_QFerrors_S", i, ".rda"))
  beta.tau.mat.list.errors <- beta.tau.est.list

  load(paste0("./inter_res/Sum_betatau_QFknown_S", i, ".rda"))
  beta.tau.mat.list.known <- beta.tau.est.list
  
  ## int_beta(tau)dtau 
  tab.betaint.list[[i]] <- rbind.data.frame(tab.list.known$tab.betaint,
                                            tab.list.errors$tab.betaint)
  
  ## beta(tau)
  beta.tau.bias.known.i <- do.call(cbind, lapply(beta.tau.mat.list.known, 
                                   function(x) x$est - beta.tau.true.vec))
  beta.tau.bias.errors.i <- do.call(cbind, lapply(beta.tau.mat.list.errors, 
                                    function(x) x$est - beta.tau.true.vec))
  
  betatau.bias.mat.errors[[i]] <- beta.tau.bias.errors.i
  betatau.bias.mat.known[[i]] <- beta.tau.bias.known.i
  
  tab.betatau.list[[i]] <- data.frame(bias = c(mean(c(beta.tau.bias.known.i)),
                                               mean(c(beta.tau.bias.errors.i))),
                                      rbind.data.frame(tab.list.known$tab.betatau[1, ],
                                                       tab.list.errors$tab.betatau[1, ]))
  
  tab.covariatevals.list[[i]] <- 
    rbind.data.frame(tab.list.known$tab.covariatevals,
                     tab.list.errors$tab.covariatevals)
  
  tab.WAIC.comp[[i]] <- tab.list.known$WAIC.mat.comp
  tab.ED.list.known[[i]] <- tab.list.known$tab.ED.list
  tab.ED.list.errors[[i]] <- tab.list.errors$tab.ED.list

}

sum_ED <- function(x, model){
  tab.1 <- data.frame(true = x$true,
                      avg = mean(x[[model]]$est),
                      bias = mean(x[[model]]$est - x$true),
                      bias.rel = mean((x[[model]]$est - x$true)/x$true),
                      MSE = mean((x[[model]]$est - x$true)^2),
                      cover = mean(x[[model]]$lower <= x$true &
                                     x$true <= x[[model]]$upper))
  return(tab.1)
}
tab.sum.ED.known <- do.call(rbind.data.frame,
                           lapply(tab.ED.list.known, sum_ED, model = "quan"))
tab.sum.ED.mean.known <- do.call(rbind.data.frame,
                                lapply(tab.ED.list.known, sum_ED, model = "mean"))
tab.sum.ED.all.known <- rbind.data.frame(tab.sum.ED.known,
                                        tab.sum.ED.mean.known)[c(sapply(1:6, 
                                                                       function(x) c(x, x+6))), ]

tab.sum.ED.errors <- do.call(rbind.data.frame,
                          lapply(tab.ED.list.errors, sum_ED, model = "quan"))

sort.index = c(sapply(1:6, function(x) c(x, x+6, x+12)))
tab.sum.ED.print <- rbind.data.frame(tab.sum.ED.all.known[seq(2, 12, 2),],
                                     tab.sum.ED.all.known[seq(1, 12, 2),],
                                     tab.sum.ED.errors)[sort.index,]

num.sce = 6
sort.index.1 = c(seq(2, 3*num.sce, 3), 
                 seq(1, 3*num.sce, 3), 
                 seq(3, 3*num.sce, 3))
tab.sum.betaint.print <- do.call(rbind.data.frame, 
                                 tab.betaint.list)[sort.index.1, ][sort.index, ]

tab.betatau.list.1 <- lapply(tab.betatau.list, function(x) rbind(NA, x))
tab.sum.betatau.print <- do.call(rbind.data.frame, tab.betatau.list.1)

tab.covariatevals.print <- do.call(rbind.data.frame, 
                                   tab.covariatevals.list)[sort.index.1, ][sort.index, ]

get_relMSE <- function(tab){
  sort.index = c(sapply(1:6, function(x) c(x, x+6, x+12)))
  ref = tab$MSE[seq(1, 18, 3)]
  quan.true = tab$MSE[seq(2, 18, 3)]
  quan.est = tab$MSE[seq(3, 18, 3)]
  re = c(ref/ref, quan.true/ref, quan.est/ref)
  re = re[sort.index]
  return(re)
}
rel.MSE.betaint <- get_relMSE(tab = tab.sum.betaint.print)
rel.MSE.covari <- get_relMSE(tab = tab.covariatevals.print)
rel.MSE.ED <- get_relMSE(tab = tab.sum.ED.print)

col1 = formatC(tab.sum.betaint.print$bias.rel, digits = 3, format = "f")
col2 = formatC(rel.MSE.betaint, digits = 2, format = "f")
col3 = tab.sum.betaint.print$coverage*100

col4 = formatC(tab.sum.betatau.print$bias, digits = 3, format = "f")
col5 = formatC(tab.sum.betatau.print$MSE, digits = 3, format = "f")
col6 = formatC(tab.sum.betatau.print$coverage*100, digits = 2, format = "f")


col7 = formatC(tab.covariatevals.print$bias.rel, digits = 3, format = "f")
col8 = formatC(rel.MSE.covari, digits = 3, format = "f")
col9 = formatC(tab.covariatevals.print$coverage*100, digits = 2, format = "f")

col10 = formatC(tab.sum.ED.print$bias.rel, digits = 3, format = "f")
col11 = formatC(rel.MSE.ED, digits = 3, format = "f")
col12 = tab.sum.ED.print$cover*100

tab.print.all <- cbind(col1, col2, col3, col4, col5, col6,
                       col7, col8, col9, col10, col11, col12)
colnames(tab.print.all) <- c("betaint.relbias", "betaint.relMSE", "betaint.cover",
                             "betatau.bias", "betatau.MSE", "betatau.cover",
                             "coval.relbias", "coval.relMSE", "coval.cover",
                             "ED.relbias", "ED.relMSE", "ED.cover")
# library(kableExtra)
# library(dplyr)
# kbl(tab.print.all, format = "latex")

dir.create("TabsFigs")
write.csv(tab.print.all, 
          file = paste0("./TabsFigs/Tab1.csv"))


##### GENERATE FIG1 ######
cols.vec <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
              "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
png(paste0("./TabsFigs/Fig1.png"),
    height = 8, width = 8, units = "in", res = 600)
tau.vec = seq(0, 1, 0.01)
plot(tau.vec, betatau.true.list[[1]](tau.vec), type="l", lty = 1, 
     xlab = expression(paste("quantile levels ", tau)),
     ylab = expression(paste("true ", beta, "(", tau, ")")),
     lwd = 2, ylim = c(0, 1.8), col = cols.vec[2])
for(i in c(2:6)){
  if(i %in% 2:4){
    lines(tau.vec, betatau.true.list[[i]](tau.vec), type="l", lty = i, 
          lwd = 2, col = cols.vec[1])
  }else if(i %in% 5:6){
    lines(tau.vec, betatau.true.list[[i]](tau.vec), type="l",
          lwd = 2, col = cols.vec[3], lty = i-3)
  }
}
legend("topleft", 
       legend = c(expression(paste("S1: ",beta, "(", tau, ") = 0.5")), 
                  expression(paste("S2: ",beta, "(", tau, ") = ", tau)),
                  expression(paste("S3: ",beta, "(", tau, ") = 1.5", tau^2)),
                  expression(paste("S4: ",beta, "(", tau, ") = ",
                                   frac(4,3), "I(",tau,"<0.5)+",
                                   frac(2,3), "I(", tau, ">=", 0.5,")")),
                  expression(paste("S5: ",beta,"(", tau, ") = ", 
                                   exp, "(", frac(-tau^2,0.328), ")")),
                  expression(paste("S6: ",beta,"(", tau, ") = ", 
                                   -tau, "+1"))),
       lty = c(1, 2:4, 2:3), lwd = 2, bty = "n",
       col = c(cols.vec[2], rep(cols.vec[1], 3), rep(cols.vec[3], 2)),
       cex = 0.8)
dev.off()



###### GENERATE FIG2 ######
### !!!! NOTES: the following chunk of data takes about 5 minutes to run 
### on a computer with 2.9 GHz 6-Core Intel Core i9 and 32 Gb memory, 
### therefore, we have provided data frames that used
### for generating Fig 2 in the folder "inter_res"

# # estimated overall quantile functions for 9 consecutive days
# ### read in simulated data ####
# load("data_sim.rda")
# x.sim.mat <- dat.all$exp$x.sim.mat
# 
# ## read in estimation results ###
# load("./inter_res/Res_estquan.rda")
# 
# source("FUNs.R")
# 
# ### true parameters ###
# alph.bar.true = 7.2
# tau12 = 1; rho1 = 0.9
# L = 4; shape = 5
# theta.bar.true = c(0.9,0.9,0.9,0.9)
# tau22 = 0.02; rho2 = 0.9
# ntimes = 1000 # number of time points
# parm.true = c(alph.bar.true, tau12, tau22, rho1, rho2)
# 
# theta.bar.mat <- re$theta.bar
# alph <- re$alph
# parm.mat <- re$parm
# theta.star.list <- re$theta.star
# epsilon = 0.01
# 
# get_est <- function(re, burn_in = 5000){
#   theta.bar.mat <- re$theta.bar
#   alph <- re$alph
#   parm.mat <- re$parm
#   tab.est <- data.frame(est = c(colMeans(theta.bar.mat[-c(1:burn_in),]),
#                                 colMeans(parm.mat[-c(1:burn_in), ])),
#                         sd = apply(cbind(theta.bar.mat, parm.mat)[-c(1:burn_in), ],
#                                    2, sd),
#                         t(apply(cbind(theta.bar.mat, parm.mat)[-c(1:burn_in), ],
#                                 2, quantile, c(0.025, 0.975))))
#   colnames(tab.est) <- c("est","sd", "lower","upper")
#   rownames(tab.est) <- c(paste0("theta",1:L), "alpha", "tau12", "tau22",
#                          "rho1", "rho2")
#   return(tab.est)
# }
# 
# ff_new <- function(x, epsilon){
#   x[x <= epsilon] <- epsilon
#   return(x)
# }
# 
# get_pointest_t <- function(re, burn_in = 5000, L){
#   
#   alph <- re$alph[-c(1:burn_in), ]
#   theta.star.list <- re$theta.star
#   theta.star.post <- lapply(theta.star.list, function(x) x[-c(1:burn_in), ])
#   theta.t.post <- lapply(theta.star.post, ff_new, epsilon = 0.01)
#   theta.t.hat <- sapply(theta.t.post, colMeans)
#   est.mat <- data.frame(alphat = colMeans(alph), theta.t.hat)
#   colnames(est.mat)[-1] = paste0("theta_", 1:L)
#   
#   sd.mat = data.frame(alphat = apply(alph, 2, sd),
#                       sapply(theta.t.post, function(x) apply(x, 2, sd)))
#   colnames(sd.mat)[-1] = paste0("theta_", 1:L)
#   
#   lower.mat = data.frame(alphat = apply(alph, 2, quantile, 0.025),
#                          sapply(theta.t.post,
#                                 function(x) apply(x, 2, quantile, 0.025)))
#   colnames(lower.mat)[-1] = paste0("theta_", 1:L)
#   
#   upper.mat = data.frame(alphat = apply(alph, 2, quantile, 0.975),
#                          sapply(theta.t.post,
#                                 function(x) apply(x, 2, quantile, 0.975)))
#   colnames(upper.mat)[-1] = paste0("theta_", 1:L)
#   
#   return(list(est = est.mat, sd = sd.mat,
#               lower = lower.mat, upper = upper.mat))
# }
# 
# 
# #### get point estimate of alpha_t0, theta_tl #####
# tab_basiscoef <- get_pointest_t(re = re, burn_in = 5000, L = L)
# theta.mat.hat <- tab_basiscoef$est
# 
# burn_in = 5000
# tau.vec = seq(0, 1, length.out = 100)
# alph <- re$alph[-c(1:burn_in), ]
# theta.star.list <- re$theta.star
# theta.star.post <- lapply(theta.star.list, function(x) x[-c(1:burn_in), ])
# theta.t.post <- lapply(theta.star.post, ff_new, epsilon = 0.01)
# quan.mat.list <- list()
# sel.time.vec = 1:9
# count = 1
# for(itime in sel.time.vec){
#   alph.theta.mat.itime <- cbind(alph[,itime],
#                                 do.call(cbind, 
#                                         lapply(theta.star.post, function(x)
#                                           x[,itime])))
#   quan.itime.post <- matrix(NA, nrow = 5001, ncol = length(tau.vec))
#   for(j in 1:5001){
#     x = alph.theta.mat.itime[j, ]
#     quan.itime.post[j, ] <- quan.fun(tau = tau.vec, L = L, theta.vec = x[2:5],
#                                      alpha = x[1], basis.fun = "Gamma", shape = 5)
#     if(j%%1000 == 0){ cat(j) }
#   }
#   quan.true.itime <- quan.fun(tau = tau.vec, L = L, 
#                               theta.vec = dat.all$exp$theta.mat.true[,itime],
#                               alpha = dat.all$exp$alph.vec.true[,itime], 
#                               basis.fun = "Gamma", shape = 5)
#   
#   quan.itime.mat <- data.frame(true = as.numeric(quan.true.itime),
#                                est = colMeans(quan.itime.post),
#                                lower = apply(quan.itime.post, 2, quantile, 0.025),
#                                upper = apply(quan.itime.post, 2, quantile, 0.975))
#   quan.mat.list[[count]] <- quan.itime.mat
#   count = count + 1
#   cat(count)
# }
# df.pl.quan = data.frame(do.call(rbind.data.frame, quan.mat.list),
#                         quan = rep(tau.vec, 9),
#                         time = rep(paste0("time ", 1:9), each = 100))
# df.pl.quan <- subset(df.pl.quan, est != Inf & quan != 0)
# 
# df.pl.obs <- do.call(rbind.data.frame, 
#                      lapply(1:9, function(x) data.frame(y = x.sim.mat[x, ], 
#                                                         time = paste0("time ", x))))
# dat.fig2 <- list(df.pl.quan = df.pl.quan,
#                  df.pl.obs = df.pl.obs)
# save(dat.fig2, file = "./inter_res/dat_fig2.rda")

load("./inter_res/dat_fig2.rda")

df.pl.obs <- dat.fig2$df.pl.obs
df.pl.quan <- dat.fig2$df.pl.quan

fig2A <- ggplot(data = df.pl.obs, aes(group = time, x = time, y = y)) + 
  geom_boxplot(notch = FALSE, outlier.shape = 1) + 
  stat_summary(fun = mean, geom="point", shape=17, size = 2, col = "black") +
  theme_bw() + labs(x = "time points", y = "observed exposures") + 
  theme(text = element_text(size=10))
col.vec = c("true" = "red", "estimated (95% CI)" = "black")
fig2B <- ggplot(df.pl.quan) + 
  geom_line(aes(x = quan, y = true, color = "true")) + 
  geom_line(aes(x = quan, y = est, color = "estimated (95% CI)")) + 
  geom_ribbon(aes(x = quan, y = est, 
                  ymin = lower, ymax = upper), alpha=0.25, color = "lightgrey") + 
  facet_wrap(~time) + 
  labs(x = "quantile levles", 
       y = paste0("true/estimated exposure quantile functions"),
       color = "") + 
  scale_color_manual(values = col.vec) + 
  theme_bw() + theme(legend.position = "bottom", 
                     text = element_text(size=10),
                     axis.text.x = element_text(size=7.5))
fig2 <- ggarrange(fig2A, fig2B, ncol = 2, labels = c("(a)", "(b)"),
                  font.label=list(color="black",size=11,face="plain"))

ggsave(filename = paste0("./TabsFigs/Fig2.png"),
       width = 12, height = 6, dpi = 600, 
       units = "in", plot = fig2,
       bg = "white")



####### GENERATE FIG 4 #######
load("./inter_res/dat_fig4.rda")
WAIC.mat <- dat.fig4$WAIC
df.betatau.forplot <- dat.fig4$betatau

pl.betatau.sce = data.frame(expand.grid(c("as_whz_any", "resp_any", "cvd_any"),
                                        c("pm25", "nox", "co", "EC")),
                            sce = 1:(3*4))
pl.betatau.list <- list()
count = 1
title.vec = c(expression(paste(PM[2.5], ", ", "as_whz")),
              expression(paste(PM[2.5], ", ", "resp")),
              expression(paste(PM[2.5], ", ", "cvd")),
              expression(paste(NO[x], ", ", "as_whz")),
              expression(paste(NO[x], ", ", "resp")),
              expression(paste(NO[x], ", ", "cvd")),
              expression(paste(CO, ", ", "as_whz")),
              expression(paste(CO, ", ", "resp")),
              expression(paste(CO, ", ", "cvd")),
              expression(paste(EC, ", ", "as_whz")),
              expression(paste(EC, ", ", "resp")),
              expression(paste(EC, ", ", "cvd")))
for(airpol in c("pm25", "nox", "co", "EC")){
  
  for(outcome in c("as_whz_any", "resp_any", "cvd_any")){
    
    waic.i = subset(WAIC.mat, exp == airpol & out == outcome)
    sub.title.i = paste0("WAIC: ", 
                         "mean (", 
                         format(round(min(waic.i[5])), big.mark = ",") ,
                         ")", 
                         " vs. quantile (",
                         format(round(min(waic.i[3:4])), big.mark = ","), 
                         ")")
    title.i = title.vec[count]
    
    df.betatau.i = subset(df.betatau.forplot,
                          exposure %in% airpol & disease %in% outcome & 
                            quan.level %in% seq(0, 1, 0.01))
    pl.betatau.i <- ggplot(df.betatau.i, 
                           aes(color = model, fill = model)) +
      geom_line(aes(x = quan.level, y = est)) + 
      geom_ribbon(aes(x = quan.level, y = est, 
                      ymin = lower, ymax = upper), 
                  alpha=0.25, color = NA) + 
      labs(x = expression(paste("quantile levels ", tau)), 
           y = expression(paste(beta, "(", tau, ")")),
           color = "Exposure", fill = "Exposure",
           subtitle = sub.title.i,
           title = title.i) + 
      scale_color_nejm(labels = c("mean", "quantile function")) + 
      scale_fill_nejm(labels = c("mean", "quantile function")) + 
      theme_bw() + theme(legend.position = "bottom", 
                         text = element_text(size=9),
                         plot.title = element_text(size=10),
                         plot.subtitle = element_text(size=8))
    pl.betatau.list[[count]] <- pl.betatau.i
    count = count + 1
  }
}
pl.betatau <- ggarrange(pl.betatau.list[[7]], 
                        pl.betatau.list[[8]], 
                        pl.betatau.list[[9]], 
                        pl.betatau.list[[10]], 
                        pl.betatau.list[[11]], 
                        pl.betatau.list[[12]],
                        pl.betatau.list[[4]], 
                        pl.betatau.list[[5]], 
                        pl.betatau.list[[6]],
                        pl.betatau.list[[1]], 
                        pl.betatau.list[[2]], 
                        pl.betatau.list[[3]], 
                        nrow = 4, ncol = 3, common.legend = T, 
                        legend = "bottom",
                        labels = paste0("(",letters[1:12],")"),
                        font.label=list(color="black",
                                        size=11,face="plain"))
ggsave(filename = paste0("./TabsFigs/Fig4.png"),
       width = 10, height = 8, units = "in", dpi = 600, 
       plot = pl.betatau,
       bg = "white")




##### GENERATE FIG 6 #######
load("./inter_res/dat_fig6.rda")
dat.map <- st_read("./inter_res/cb_2015_us_zcta510_500k/cb_2015_us_zcta510_500k.shp")
dat.difAD.byzip <- dat_fig6$AD
zip.metro <- dat_fig6$zip
index = dat.map$ZCTA5CE10 %in% zip.metro
dat.map.sub <- dat.map[index, ]

df.pl.difAD.sub1 <- left_join(dat.map.sub, 
                              subset(dat.difAD.byzip, exposure == "co" & disease == "cvd_any"),
                              by = c("ZCTA5CE10" = "zip"))
pl.co.cvd <- tm_shape(df.pl.difAD.sub1) + 
  tm_polygons("dif", title = "Relative \n differences", 
              style = "cont", midpoint = NA, 
              palette = "RdYlBu", breaks = seq(-1, 1.4, 0.4)) + 
  tm_text("ZCTA5CE10", size = "AREA") + 
  tm_layout(inner.margins = 0.01, legend.title.size = 1,
            legend.position = c("right", "bottom"), frame = FALSE,
            title = "(a) CO, cvd", title.size = 1)

df.pl.difAD.sub2 <- left_join(dat.map.sub, 
                              subset(dat.difAD.byzip, exposure == "EC" & disease == "resp_any"),
                              by = c("ZCTA5CE10" = "zip"))
pl.EC.resp <- tm_shape(df.pl.difAD.sub2) + 
  tm_polygons("dif", title = "Relative \n differences", 
              style = "cont", midpoint = NA, 
              palette = "RdYlBu", breaks = seq(-1, 1.4, 0.4)) + 
  tm_text("ZCTA5CE10", size = "AREA") + 
  tm_layout(inner.margins = 0.01, legend.title.size = 1,
            legend.position = c("right", "bottom"), frame = FALSE,
            title = "(b) EC, resp", title.size = 1)

df.pl.difAD.sub3 <- left_join(dat.map.sub, 
                              subset(dat.difAD.byzip, 
                                     exposure == "pm25" & disease == "as_whz_any"),
                              by = c("ZCTA5CE10" = "zip"))
pl.pm25.aswhz <- tm_shape(df.pl.difAD.sub3) + 
  tm_polygons("dif", title = "Relative \n differences", 
              style = "cont", midpoint = NA,
              palette = "RdYlBu", breaks = seq(-1, 1.4, 0.4)) + 
  tm_text("ZCTA5CE10", size = "AREA") + 
  tm_layout(inner.margins = 0.01, legend.title.size = 1,
            legend.position = c("right", "bottom"), frame = FALSE,
            title = "(c) PM2.5, as_whz", title.size = 1)

pl.map.AD <- tmap_arrange(pl.co.cvd, pl.EC.resp, pl.pm25.aswhz, 
                          ncol = 3, outer.margins = rep(0.01, 4))
tmap_save(tm = pl.map.AD,
          filename = paste0("./TabsFigs/Fig6.png"),
          width = 15, height = 6, units = "in", dpi = 480)






