
######################################
###    Simulation experiments      ###
######################################

rm(list=ls())

# Required packages and auxiliary functions are loaded at the start
source('AuxFunctions.R')
library(MASS)
library(dplyr)
library(sae)
library(nlme)
library(betareg)


###########################
### Simulation scenario ###
###########################

# Simulation parameters (B, R, D, H) are clearly defined
B=200 # number of iterations
R=200 # bootstrap replicates
D=40  # number of areas
H=10  # MC iterations for the FHMQ (Marchetti et al., 2025)

# Population sizes (Nd) and sample sizes (nd) are calculated based on the target sampling fraction
Nd=rep(100,D)
N=sum(Nd)
n=N*0.05 # p
regioncode=rep(1:D,Nd)

# Continuous covariates (xdj1, xdj2, xdj3) are simulated from lognormal and normal distributions
xdj1 <- rlnorm(N,  0.75, 0.5)
xdj2 <- rlnorm(N,  1,   0.25)
xdj3 <- rnorm(N,   2,  0.75)

# Sampling weights and probabilities (fdj, pdj, wdj) are computed for unequal probability sampling
fdj <- runif(N, 1, 5)       
pdj <- n * fdj / sum(fdj)   
wdj <- 1/pdj

# Sample and non-sample subsets
s_list <- lapply(split(1:N, regioncode), function(idx) {
  sample(idx, 1, prob = pdj[idx]) # guarantee at least 1 from each region
})

# remaining units sampled as before
s <- sort(c(unlist(s_list), sample(setdiff(1:N, unlist(s_list)), 
            n - length(s_list), prob = pdj[setdiff(1:N, unlist(s_list))])))

x.s<-xdj1[s]
regioncode.s<-regioncode[s]
wdj.s<-wdj[s]

x.r<-xdj1[-s]
regioncode.r<-regioncode[-s]
wdj.r<-wdj[-s]

nd=as.numeric(table(regioncode.s))

# Quantile grids (tau and tau.large) are generated for M-quantile estimation
tau <- sort(unique(c(seq(0.0015,0.9985,0.05), 
                     1- seq(0.0015,0.9985,0.05)), 0.5))
tau.large <- sort(unique(c(seq(0.0015,0.9985,0.02), 
                           1- seq(0.0015,0.9985,0.02)), 0.5))
GALI01 <- lapply(tau.large, simu.GALI01)

# Outliers 
mean.ui<-9
mean.e<-20
ki<-as.integer(.1 * D) 

# Both single-covariate and multi-covariate scenarios are implemented (more.covariates flag)
  more.covariates = T
# more.covariates = T if Xd=(Xd1, Xd2, Xd3) from census files
# more.covariates = F if Xd=Xd1 from census files

if(more.covariates==T){ 
  Xd <- sapply(list(xdj1, xdj2, xdj3), function(x) {
        aggregate(x, by = list(regioncode), mean)[, 2]})
} else { Xd <- aggregate(xdj1, by=list(regioncode), mean)[, 2]
}
  
Xd.j <- merge(data.frame('regioncode.s' =regioncode.s),
              data.frame('regioncode.s'=1:D, Xd), 
                   by='regioncode.s')[, -1] 

res.s <- pred.Hajek <-  pred.FH <- pred.FH.R <- pred.FH2 <- pred.NER <- 
         pred.MQ <- pred.AMQ <- pred.AMQ2 <- Yd <- pred.AMQ.ell1 <- 
         sd.Hajek <- sd.Hajek2 <- sd <- MSE.AMQ <- MSE.AMQ2 <- list()

# Estimate the MSE or not
do.MSE <- TRUE

#############################
### Start the simulations ###
#############################

for (b in 1:B){
  set.seed(b)
  print(b)
  
  # Outliers in area effects (ud) and residuals (edj) are introduced optionally 
  ud <- rnorm(D,0,sqrt(3))
  # % contamination in ui
  out.ui <- 0
  uui<-ud
  u1 <- rnorm(ki, mean.ui, sqrt(20))
  uui[(D-ki+1):D]<-u1
  out.ui <- rep(out.ui, D)
  ud <- ifelse(out.ui > 0, uui, ud)
  
  # % contamination in e
  out.e <- 0.00
  n1<-rbinom(N,1,out.e)
  edj <- (1-n1)*rnorm(N,0,sqrt(6))+n1*rnorm(N, mean.e, sqrt(150))
  
  if(more.covariates==T){
    ydj <- 100 + 5*xdj1 + 3*xdj2 + 2*xdj3 + rep(ud,Nd) + edj
  } else {
    ydj <- 100 + 5*xdj1 + rep(ud,Nd) + edj
  }
  
  Yd[[b]] <- aggregate(ydj, by=list(regioncode), mean)[, 2]
  
  # Sample observations
  y.s <- ydj[s]
  data.s <- data.frame(regioncode.s, y.s, x.s, wdj.s)
  sum.s <- aggregate(y.s, by=list(regioncode.s), sum)$x
  
  
  #####################################
  ###   1. Hajek direct estimates   ###
  #####################################
  # Direct Hajek estimators and their variances are computed per area
  
  Y.d.Hajek <- as.data.frame(data.s %>% mutate(w = wdj.s) %>% 
                            group_by(regioncode.s) %>%
                    summarise(Y.dir = sum(y.s * w) / sum(w),
                    Y.var = sum(w * (w - 1) * (y.s - Y.dir)^2) / (sum(w)^2),
                    Nd.hat =(sum(w)^2),
      .groups = "drop")) %>% mutate(Y.var = ifelse(Y.var == 0, 
                                    min(Y.var[Y.var != 0]), Y.var))
  pred.Hajek[[b]] <- Y.d.Hajek$Y.dir
  sd.Hajek[[b]] <- sqrt(Y.d.Hajek$Y.var)
  
  covHajek <- cov.Hajek(regioncode.s, y.s, wdj.s)
  
  # Multiple small area predictors are computed: FH, robust FH, GVF-adjusted FH, 
  # NER, unit-level MQ, and area-level MQ (AMQ l=1,2, and FHMQ Marchetti et al. 2025)
  
  ################################
  ###   2. Fay-Herriot model   ###
  ################################
  
  FHmod <- eblupFH(Y.d.Hajek$Y.dir ~ Xd, vardir = Y.d.Hajek$Y.var)
  pred.FH[[b]] <- FHmod$eblup   # EBLUP
  
  gamma.d <- FHmod$fit$refvar/(FHmod$fit$refvar+Y.d.Hajek$Y.var)
  # u.d <- gamma.d * (Y.d.Hajek$Y.dir - cbind(1, Xd) %*% FHmod$fit$estcoef$beta)
  # EBLUP <- cbind(1, Xd) %*% FHmod$fit$estcoef$beta + u.d
  
  # REBLUP
  u.d.rob <- gamma.d*hub.psi(Y.d.Hajek$Y.dir - cbind(1, Xd) %*% FHmod$fit$estcoef$beta, k=3)
  pred.FH.R[[b]] <- cbind(1, Xd) %*% FHmod$fit$estcoef$beta + u.d.rob
  
  # GVF method
  GVF.mod <- lm(log(Y.d.Hajek$Y.var) ~  Y.d.Hajek$Y.dir + nd + Y.d.Hajek$Y.dir*nd)
  GVF.sel <- step(GVF.mod, direction = "backward", trace = 0)
  
  v <- deviance(GVF.sel)/df.residual(GVF.sel)
  sd.Hajek2[[b]] <- sqrt(exp(v/2)*exp(predict(GVF.sel)))
  
  FHmod2 <- eblupFH(Y.d.Hajek$Y.dir ~ Xd, vardir = sd.Hajek2[[b]]^2)
  pred.FH2[[b]] <- FHmod2$eblup
  
  
  ############################################
  ###   3. Nested error regression model   ###
  ############################################
  
  modelo.NER <- lme(fixed = y.s~x.s, random=~1|regioncode.s, data = data.s)
  
  pred.NER[[b]] <- 1/Nd *(sum.s + sapply(1:D, function(d) 
    sum(as.matrix(cbind(1,x.r)[regioncode.r==d, ]) %*% 
          t(as.matrix(coef(modelo.NER)[d, ])))))

  
  ##################################
  ###   4. Unit-level MQ model   ###
  ##################################
  
  mod <- QRLM(x=cbind(1,x.s), y=y.s, q=tau, maxit=30, k = 1.345)
  qo <- matrix(c(gridfitinter(y.s,mod$fitted.values,
                              mod$q.values)),nrow=n,ncol=1)
  qmat <- data.frame(qo, regioncode.s)
  mqo <- aggregate(qmat[,1],by=list(qmat[,2]),mean)[,2]
  
  mod.SAE <- QRLM(x=cbind(1,x.s), y=y.s, q=mqo, maxit=30, k = 1.345)
  
  pred.MQ[[b]] <- 1/Nd *(sum.s + sapply(1:D, function(d) 
    sum(as.matrix(cbind(1,x.r)[regioncode.r==d, ]) %*% mod.SAE$coef[,d])))
  
  
  #########################################
  ###   5. Area-level MQ model \ell=1   ###
  #########################################
  
  mod.auxi <- QRLM(x=cbind(1, Xd), y=Y.d.Hajek$Y.dir, q=tau, maxit=30, k = 1.345)
  mqol1 <- matrix(c(gridfitinter(Y.d.Hajek$Y.dir, mod.auxi$fitted.values,
                                    mod.auxi$q.values)), nrow=D, ncol=1)
  mqol1 <- 0.5 * (tanh(2 * (mqol1 - 0.5)) + 1)

  mod.mqol1 <- QRLM(x=cbind(1, Xd), y=Y.d.Hajek$Y.dir, q=mqol1, maxit=30, k = 1.345)
  pred.AMQ.ell1[[b]] <- sapply(1:D, function(d) cbind(1,Xd)[d, ] %*% mod.mqol1$coef[, d])
  
  
  #########################################
  ###   6. Area-level MQ model \ell=2   ###
  #########################################
  
  mod.area <- QRLM(x=cbind(1, Xd.j), y=y.s, q=tau, maxit=30, k = 1.345)
  mqo.area <- matrix(c(gridfitinter(y.s, mod.area$fitted.values,
                                    mod.area$q.values)), nrow=n, ncol=1)
  qmat.area <- data.frame(mqo.area, regioncode.s)
  mqo.area  <- aggregate(qmat.area[,1], by=list(qmat.area [,2]), mean)[,2]
  
  mod.area.50 <-  QRLM(x=cbind(1, Xd), y=Y.d.Hajek$Y.dir, q=0.50, maxit=30, k = 1.345)
  mod.area.SAE <- QRLM(x=cbind(1, Xd), y=Y.d.Hajek$Y.dir, q=mqo.area, maxit=30, k = 1.345)
  
  res.s <- lapply(1:D, function(d){ Y.d.Hajek$Y.dir-cbind(1,Xd)%*%mod.area.SAE$coef[, d] } )
  sd[[b]] <- sapply(res.s, fun.MAD)
  
  pred.AMQ[[b]] <- sapply(1:D, function(d) cbind(1,Xd)[d, ] %*% mod.area.SAE$coef[, d])
  
  u.MQ <- sapply(1:D, function(d) cbind(1,Xd)[d, ]%*%(mod.area.SAE$coef[, d] - mod.area.50$coef))
  u.MQ.C <- u.MQ - mean(u.MQ)
  
  # Analytical and bootstrap MSE estimators for AMQ are computed within the loop

  # MSE estimation for AMQ models
  if(do.MSE == TRUE){
    ## Analytical MSE estimation
    MSE.AMQ[[b]] <- MQpred.mse(Y.d.Hajek$Y.dir, cbind(1,Xd), Y.d.Hajek$Y.var, 
                               covHajek, mod.area.SAE, D)
    
    ## Bootstrap MSE estimation
    MSE.AMQ2[[b]] <- MSE.boot(mqo.area, mod.area.50, u.MQ.C, tau.large, 
                              GALI01, MSE.AMQ[[b]], Y.d.Hajek$Y.dir, cbind(1,Xd),
                              sd[[b]], D, R, mod.area.SAE)
  }
  
  ################################################################
  ###   7. Area-level MQ model from Marchetti et al., (2025)   ###
  ################################################################
  
  pred.AMQ2[[b]] <- compute_AMQ2(Xd, Y.d.Hajek, tau, D, H)
  
}

#########################################################
###   PERFORMANCE OF THE POPULATION MEAN PREDICTORS   ###
#########################################################

# Performance metrics (RRMSE, RBIAS, ARBIAS) are calculated across all simulation replicates

Yd.mean <- Reduce('+', Yd)/B

rrmse.Hajek    <- sqrt(Reduce('+', lapply(Map('-', pred.Hajek, Yd), function(x){ x^2 }))/B)/Yd.mean
rrmse.FH       <- sqrt(rowMeans(sapply(1:B, function(b) (pred.FH[[b]]   - Yd[[b]])^2), na.rm = TRUE))/Yd.mean
rrmse.FH.R     <- sqrt(rowMeans(sapply(1:B, function(b) (pred.FH.R[[b]] - Yd[[b]])^2), na.rm = TRUE))/Yd.mean
rrmse.FH2      <- sqrt(rowMeans(sapply(1:B, function(b) (pred.FH2[[b]]  - Yd[[b]])^2), na.rm = TRUE))/Yd.mean
rrmse.NER      <- sqrt(Reduce('+', lapply(Map('-', pred.NER,   Yd), function(x){ x^2 }))/B)/Yd.mean
rrmse.MQ       <- sqrt(Reduce('+', lapply(Map('-', pred.MQ,    Yd), function(x){ x^2 }))/B)/Yd.mean
rrmse.AMQ.ell1 <- sqrt(Reduce('+', lapply(Map('-', pred.AMQ.ell1,   Yd), function(x){ x^2 }))/B)/Yd.mean
rrmse.AMQ      <- sqrt(Reduce('+', lapply(Map('-', pred.AMQ,   Yd), function(x){ x^2 }))/B)/Yd.mean
rrmse.AMQ2     <- sqrt(Reduce('+', lapply(Map('-', pred.AMQ2,  Yd), function(x){ x^2 }))/B)/Yd.mean

arbias.Hajek    <- (Reduce('+', lapply(Map('-', pred.Hajek, Yd), function(x){ abs(x) }))/B)/Yd.mean
arbias.FH       <- rowMeans(sapply(1:B, function(b) abs(pred.FH[[b]]   - Yd[[b]])), na.rm = TRUE)/Yd.mean
arbias.FH.R     <- rowMeans(sapply(1:B, function(b) abs(pred.FH.R[[b]] - Yd[[b]])), na.rm = TRUE)/Yd.mean
arbias.FH2      <- rowMeans(sapply(1:B, function(b) abs(pred.FH2[[b]]  - Yd[[b]])), na.rm = TRUE)/Yd.mean
arbias.NER      <- (Reduce('+', lapply(Map('-', pred.NER,   Yd), function(x){ abs(x) }))/B)/Yd.mean
arbias.MQ       <- (Reduce('+', lapply(Map('-', pred.MQ,    Yd), function(x){ abs(x) }))/B)/Yd.mean
arbias.AMQ      <- (Reduce('+', lapply(Map('-', pred.AMQ,   Yd), function(x){ abs(x) }))/B)/Yd.mean
arbias.AMQ.ell1 <- (Reduce('+', lapply(Map('-', pred.AMQ,   Yd), function(x){ abs(x) }))/B)/Yd.mean
arbias.AMQ2     <- (Reduce('+', lapply(Map('-', pred.AMQ2,  Yd), function(x){ abs(x) }))/B)/Yd.mean

rbias.Hajek    <- (Reduce('+', lapply(Map('-', pred.Hajek, Yd), function(x){ x }))/B)/Yd.mean
rbias.FH       <- rowMeans(sapply(1:B, function(b) (pred.FH[[b]]   - Yd[[b]])), na.rm = TRUE)/Yd.mean
rbias.FH.R     <- rowMeans(sapply(1:B, function(b) (pred.FH.R[[b]] - Yd[[b]])), na.rm = TRUE)/Yd.mean
rbias.FH2      <- rowMeans(sapply(1:B, function(b) (pred.FH2[[b]]  - Yd[[b]])), na.rm = TRUE)/Yd.mean
rbias.NER      <- (Reduce('+', lapply(Map('-', pred.NER,   Yd), function(x){ x }))/B)/Yd.mean
rbias.MQ       <- (Reduce('+', lapply(Map('-', pred.MQ,    Yd), function(x){ x }))/B)/Yd.mean
rbias.AMQ      <- (Reduce('+', lapply(Map('-', pred.AMQ,   Yd), function(x){ x }))/B)/Yd.mean
rbias.AMQ.ell1 <- (Reduce('+', lapply(Map('-', pred.AMQ.ell1,   Yd), function(x){ x }))/B)/Yd.mean
rbias.AMQ2     <- (Reduce('+', lapply(Map('-', pred.AMQ2,  Yd), function(x){ x }))/B)/Yd.mean

####################
###   Boxplots   ### 550 x 500
####################

# Boxplots visualize relative RMSE and bias for all models

box.rrmse  <- data.frame(rrmse.Hajek, rrmse.FH, rrmse.FH.R, rrmse.FH2, rrmse.AMQ.ell1,
                        rrmse.AMQ, rrmse.AMQ2, rrmse.NER, rrmse.MQ)

box.rbias <- data.frame(rbias.Hajek, rbias.FH, rbias.FH.R, rbias.FH2, rbias.AMQ.ell1,
                         rbias.AMQ, rbias.AMQ2, rbias.NER,  rbias.MQ)

par(mfrow = c(1, 1),
    mar = c(4, 4, 2, 2),
    oma = c(0, 0, 0, 0),
    mgp = c(1.6, 0.5, 0))

boxplot(box.rrmse, names = c("dir", "fh1", "fhr", "fh2", "amq1", "amq2", 
                             "fhmq", "ner", "mq"), ylab='',
        las=2, main = "RRMSE (%)", xlab = "", pch=19, 
        cex.axis = 1.35, cex.lab = 1.35, ylim=c(0.00, 0.09))

boxplot(box.rbias, names = c("dir", "fh1", "fhr", "fh2", "amq1", "amq2", 
                             "fhmq", "ner", "mq"),
        las=2, main = "RBIAS (%)", xlab = "", pch=19,  ylab='',
        cex.axis = 1.35, cex.lab = 1.35, ylim=c(-0.08, 0.04))
abline(h=0, col='red', lty=2, lwd=2)


###########################
###   Tabular results   ###
###########################

# Reference values: AMQ
rrmse_AMQ <- round(mean(rrmse.AMQ), 4)
arbias_AMQ <- round(mean(arbias.AMQ), 4)

k.sec <- D
# k.sec <- as.integer(.1 * D)

rrmse_AMQ <- round(mean(rrmse.AMQ[(D-k.sec+1):D]), 4) # or [1:(D-k.sec)]
results <- data.frame(Model = c("Hajek", "FH", "FH.R", "FH2", "AMQell1", "AMQ", "AMQ2", "NER",  "MQ"),
           ARBIAS = round(c(mean(arbias.Hajek[(D-k.sec+1):D]), mean(arbias.FH[(D-k.sec+1):D]), 
                            mean(arbias.FH.R[(D-k.sec+1):D]),  mean(arbias.FH2[(D-k.sec+1):D]),   
                            mean(arbias.AMQ.ell1[(D-k.sec+1):D]),
                            mean(arbias.AMQ[(D-k.sec+1):D]),   mean(arbias.AMQ2[(D-k.sec+1):D]),  
                            mean(arbias.NER[(D-k.sec+1):D]),   mean(arbias.MQ[(D-k.sec+1):D])), 4),   
           RRMSE = round(c(mean(rrmse.Hajek[(D-k.sec+1):D]),   mean(rrmse.FH[(D-k.sec+1):D]),
                            mean(rrmse.FH.R[(D-k.sec+1):D]),   mean(rrmse.FH2[(D-k.sec+1):D]), 
                            mean(rrmse.AMQ.ell1[(D-k.sec+1):D]),
                            mean(rrmse.AMQ[(D-k.sec+1):D]),    mean(rrmse.AMQ2[(D-k.sec+1):D]),  
                            mean(rrmse.NER[(D-k.sec+1):D]),    mean(rrmse.MQ[(D-k.sec+1):D])), 4))
results$RRMSE_INCR <- round((results$RRMSE - rrmse_AMQ) / rrmse_AMQ * 100, 2)
results


######################################################
###   PERFORMANCE OF THE PROPOSED MSE ESTIMATORS   ###
######################################################

rmse.AMQ   <- sqrt(Reduce('+', lapply(Map('-', pred.AMQ,     Yd), function(x){ x^2 }))/B)
rmse.Hajek <- sqrt(Reduce('+', lapply(Map('-', pred.Hajek,   Yd), function(x){ x^2 }))/B)

RMSE.mse <- sqrt(apply(sapply(MSE.AMQ, function(x){ (sqrt(x[[1]]) - rmse.AMQ)^2 }), 1, mean))
BIAS.mse <- apply(sapply(MSE.AMQ, function(x){ (sqrt(x[[1]]) - rmse.AMQ) }), 1, mean)

RMSE.mse2 <- sqrt(apply(sapply(MSE.AMQ2, function(x){ (sqrt(x[[1]]) - rmse.AMQ)^2 }), 1, mean))
BIAS.mse2 <- apply(sapply(MSE.AMQ2, function(x){ (sqrt(x[[1]]) - rmse.AMQ) }), 1, mean)

RMSE.mse3 <- sqrt(apply(sapply(sd.Hajek, function(x){ (x - rmse.Hajek)^2 }), 1, mean))
BIAS.mse3 <- apply(sapply(sd.Hajek, function(x){ (x - rmse.Hajek) }), 1, mean)

RMSE <- cbind(RMSE.mse, RMSE.mse2, RMSE.mse3)
BIAS <- cbind(BIAS.mse, BIAS.mse2, BIAS.mse3)

# Results are saved to CSV for further analysis and reproducibility
# write.csv(data.frame(RMSE, BIAS), row.names = FALSE, file = paste0("MSEresultsB/MSEresults_R", R, ".csv"))

####################
###   Boxplots   ###
####################

par(mfrow = c(1, 1),
    mar = c(4, 4, 1, 1),
    oma = c(0, 0, 0, 0),
    mgp = c(1.6, 0.5, 0))

boxplot(100*RMSE/rmse.AMQ, names = c("AMQ analytic", "AMQ boot", "Hajek"), ylim=c(10,200),
        ylab = "RRMSE (%)", xlab = "", pch=19, cex.axis = 1.35, cex.lab = 1.35)

boxplot(100*BIAS/rmse.AMQ, names = c("AMQ analytic", "AMQ boot", "Hajek"), ylim=c(-150,200),
        ylab = "RBIAS (%)", xlab = "", pch=19, cex.axis = 1.35, cex.lab = 1.35)
abline(h=0, col='red', lty=2, lwd=2)

