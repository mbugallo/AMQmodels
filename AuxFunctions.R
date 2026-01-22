
###########################
###   HAJEK COVARIANCE  ###
###########################

cov.Hajek <- function(regioncode, y, w) {
  areas <- sort(unique(regioncode))
  D <- length(areas)
  
  # total weights and Hajek averages by area
  Wtot <- tapply(w, regioncode, sum)[areas]
  Ydir <- tapply(w*y, regioncode, sum)[areas] / Wtot
  
  # residuals weighted by area
  residuales <- lapply(areas, function(a) {
    idx <- regioncode == a
    w[idx]*(y[idx] - Ydir[which(areas==a)])
  })
  
  # internal function for each pair
  cov_fun <- function(d, g) {
    if (d == g) {
      idx <- regioncode == d
      sum(w[idx]*(w[idx]-1)*(y[idx]-Ydir[which(areas==d)])^2)/(Wtot[which(areas==d)]^2)
    } else {
      -sum(outer(residuales[[which(areas==d)]],
                 residuales[[which(areas==g)]]))/
        (Wtot[which(areas==d)]*Wtot[which(areas==g)])
    }
  }
  
  outer(areas, areas, Vectorize(cov_fun))
}


##########################################
###   MQ REGRESSION FITTING ALGORITHM  ###
##########################################

QRLM <-function (x, y, case.weights = rep(1, nrow(x)), k=1.345, 
                 var.weights = rep(1, nrow(x)), ..., w = rep(1, nrow(x)), 
                 init = "ls", psi = psi.huber, 
                 scale.est = c("MAD", "Huber", "proposal 2"), 
                 k2 = 1.345, method = c("M", "MM"), maxit = 20, 
                 acc = 1e-04, test.vec = "resid", q = 0.5)
{
  irls.delta <- function(old, new) sqrt(sum((old - new)^2)/max(1e-20, sum(old^2)))
  irls.rrxwr <- function(x, w, r) {
    w <- sqrt(w)
    max(abs((matrix(r*w,1,length(r)) %*% x)/sqrt(matrix(w,1,length(r)) %*% (x^2))))/sqrt(sum(w*r^2))
  }
  method <- match.arg(method)
  nmx <- deparse(substitute(x))
  if (is.null(dim(x))) {
    x <- as.matrix(x)
    colnames(x) <- nmx
  }
  else x <- as.matrix(x)
  if (is.null(colnames(x)))
    colnames(x) <- paste("X", seq(ncol(x)), sep = "")
  if (qr(x)$rank < ncol(x))
    stop("x is singular: singular fits are not implemented in rlm")
  if (!(any(test.vec == c("resid", "coef", "w", "NULL")) || is.null(test.vec)))
    stop("invalid testvec")
  if (length(var.weights) != nrow(x))
    stop("Length of var.weights must equal number of observations")
  if (any(var.weights < 0))
    stop("Negative var.weights value")
  if (length(case.weights) != nrow(x))
    stop("Length of case.weights must equal number of observations")
  w <- (w * case.weights)/var.weights
  if (method == "M") {
    scale.est <- match.arg(scale.est)
    if (!is.function(psi))
      psi <- get(psi, mode = "function")
    arguments <- list(...)
    if (length(arguments)) {
      pm <- pmatch(names(arguments), names(formals(psi)), nomatch = 0)
      if (any(pm == 0))
        warning(paste("some of ... do not match"))
      pm <- names(arguments)[pm > 0]
      formals(psi)[pm] <- unlist(arguments[pm])
    }
    if (is.character(init)) {
      if (init == "ls")
        temp <- lm.wfit(x, y, w, method = "qr")
      else if (init == "lts")
        temp <- lqs.default(x, y, intercept = FALSE, nsamp = 200)
      else stop("init method is unknown")
      coef <- temp$coef
      resid <- temp$resid
    }
    else {
      if (is.list(init))
        coef <- init$coef
      else coef <- init
      resid <- y - x %*% coef
    }
  }
  else if (method == "MM") {
    scale.est <- "MM"
    temp <- lqs.default(x, y, intercept = FALSE, method = "S", k0 = 1.548)
    coef <- temp$coef
    resid <- temp$resid
    psi <- psi.bisquare
    if (length(arguments <- list(...)))
      if (match("c", names(arguments), nomatch = FALSE)) {
        c0 <- arguments$c
        if (c0 > 1.548) {
          psi$c <- c0
        }
        else warning("c must be at least 1.548 and has been ignored")
      }
    scale <- temp$scale
  }
  else stop("method is unknown")
  done <- FALSE
  conv <- NULL
  n1 <- nrow(x) - ncol(x)
  if (scale.est != "MM")
    scale <- mad(resid/sqrt(var.weights), 0)
  
  qest <- matrix(0, nrow = ncol(x), ncol = length(q))
  qwt <- matrix(0, nrow = nrow(x), ncol = length(q))
  qfit <- matrix(0, nrow = nrow(x), ncol = length(q))
  qres <- matrix(0, nrow = nrow(x), ncol = length(q))
  qvar.matrix <- array(rep(0,ncol(x)*ncol(x)),dim=c(ncol(x),ncol(x),length(q)))
  qscale <- NULL
  for(i in 1:length(q)) {
    for (iiter in 1:maxit) {
      if (!is.null(test.vec))
        testpv <- get(test.vec)
      if (scale.est != "MM") {
        if (scale.est == "MAD")
          scale <- median(abs(resid/sqrt(var.weights)))/0.6745
        else {gamma<- 4*k2^2*(1-pnorm(k2))*((1-q[i])^2+q[i]^2) - 4*k2*dnorm(k2)*((1-q[i])^2+q[i]^2) + 
          4*(1-q[i])^2*(pnorm(0)-(1-pnorm(k2))) + 4*q[i]^2*(pnorm(k2)-pnorm(0))
        scale <- sqrt(sum(pmin(resid^2/var.weights,(k2*scale)^2))/(n1*gamma))
        }
        if (scale == 0) {
          done <- TRUE
          break
        }
      }
      w <- psi(resid/(scale * sqrt(var.weights)),k=k) * case.weights
      ww <- 2 * (1 - q[i]) * w
      ww[resid > 0] <- 2 * q[i] * w[resid > 0]
      w <- ww
      temp <- lm.wfit(x, y, w, method = "qr")
      coef <- temp$coef
      resid <- temp$residuals
      if (!is.null(test.vec))
        convi <- irls.delta(testpv, get(test.vec))
      else convi <- irls.rrxwr(x, wmod, resid)
      conv <- c(conv, convi)
      done <- (convi <= acc)
      if (done)
        break
    }
    if (!done)
      warning(paste("rlm failed to converge in", maxit, "steps at q = ", q[i]))
    qest[, i] <- coef
    qscale[i]<-scale
    qwt[, i] <- w
    qfit[, i] <- temp$fitted.values
    qres[,i] <- resid
    
    tmp.res.mq<-qres[,i]/qscale[i]
    Epsi2<-(sum((qwt[,i]*tmp.res.mq)^2)/(nrow(x)-ncol(x)))
    Epsi<-(1/qscale[i])*(sum(2*(q[i]*(0<=tmp.res.mq & tmp.res.mq<= k)+
                                  (1-q[i])*(-k <=tmp.res.mq & tmp.res.mq<0)))/nrow(x))
    qvar.matrix[,,i]<- (((Epsi2)/Epsi^2)*solve(t(x)%*%x))
    
  }
  list(fitted.values = qfit, residuals = qres, q.values = q, q.weights = qwt, coef= qest,
       qscale=qscale,var.beta=qvar.matrix)
}


###########################################
###   COMPUTING OF THE QUANTILE-ORDERS  ###
###########################################

"zerovalinter"<-function(y, x)
{
  if(min(y) > 0) {
    xmin <- x[y == min(y)]
    if(length(xmin) > 0)
      xmin <- xmin[length(xmin)]
    xzero <- xmin
  }
  else {
    if(max(y) < 0) {
      xmin <- x[y == max(y)]
      if(length(xmin) > 0)
        xmin <- xmin[1]
      xzero <- xmin
    }
    else {
      y1 <- min(y[y > 0])
      if(length(y1) > 0)
        y1 <- y1[length(y1)]
      y2 <- max(y[y < 0])
      if(length(y2) > 0)
        y2 <- y2[1]
      x1 <- x[y == y1]
      if(length(x1) > 0)
        x1 <- x1[length(x1)]
      x2 <- x[y == y2]
      if(length(x2) > 0)
        x2 <- x2[1]
      xzero <- (x2 * y1 - x1 * y2)/(y1 - y2)
      xmin <- x1
      if(abs(y2) < y1)
        xmin <- x2
    }
  }
  resu <-  xzero
  resu
}


###########################################
###    LINEAR INTERPOLATION FUNCTION    ###
###########################################

"gridfitinter"<-function(y,expectile,Q)
{
  nq<-length(Q)
  diff <- y %*% t(as.matrix(rep(1, nq))) - expectile        
  vectordest <- apply(diff, 1, zerovalinter,Q)    
}


#######################################
###    MEDIAN ABSOLUTE DEVIATION    ###
#######################################

fun.MAD <- function(x){ median(abs( x-median(x)) )/0.6745 }


###################################################
###  MODIFIED MQ REGRESSION FITTING ALGORITHM   ###
###################################################

# Huber's psi function (symmetric)
hub.psi <- function(x, k){ ifelse(abs(x) <= k, x, sign(x) * k) }

# Asymmetric influence function for quantile q
hub.psi.q <- function(x, k, q){  2 * hub.psi(x, k) * (ifelse(x > 0, q, 1 - q)) }

QRLMmod <- function(Xdd, Y.dir, var.dir, q) {
  
  beta <- coef(lm(Y.dir ~ -1 + Xdd))
  m <- nrow(Xdd)
  
  tol = 1e-6
  max_iter = 200
  
  for (iter in 1:max_iter) {
    # Initialize estimating equation (score function)
    eq <- rep(0, length(beta))
  
    resid <- Y.dir - Xdd %*% beta            
    psi_vals <- hub.psi.q(resid / sqrt(var.dir), k = 1.345, q)  
    eq <- t(Xdd) %*% psi_vals     
    
    # Newton-Raphson update
    beta_new <- beta + solve(t(Xdd) %*% Xdd) %*% eq
    
    # Check convergence
    if ((max(abs(beta_new - beta)) < tol)|(iter==max_iter)) {
      return(list(fitted.values=Xdd %*%beta_new, beta=beta_new, q.values=q))
    }
  }
  
}


######################################################
###  AMQ REGRESSION FROM MARCHETTI ET AL. (2025)   ###
######################################################

compute_AMQ2 <- function(Xd, Y.d.Hajek, tau, D, H) {
  
  # Step 1: Fit the MQ models for all tau values
  QRLMmod.result <- lapply(tau, function(q) {
          QRLMmod(cbind(1, Xd), Y.d.Hajek$Y.dir, Y.d.Hajek$Y.var, q = q)})
  
  # Step 2: Interpolate M-quantiles to get tau for area-level SAE
  mqo.QRLMmod <- matrix(gridfitinter(Y.d.Hajek$Y.dir,
                        sapply(QRLMmod.result, function(x) x[[1]]),
                        tau), nrow = D, ncol = 1)
  
  # Step 3: Fit MQ SAE models at interpolated tau
  QRLM.SAE <- lapply(mqo.QRLMmod, function(q) {
      QRLMmod(cbind(1, Xd), Y.d.Hajek$Y.dir, Y.d.Hajek$Y.var, q = q) })
  
  QRLM.SAE.coef <- sapply(QRLM.SAE, function(x) x[[2]])
  
  # Step 4: Compute area-level predictions
  Y.d.pred <- sapply(1:D, function(d) cbind(1, Xd)[d, ] %*% QRLM.SAE.coef[, d])
  
  # Step 5: Monte Carlo loop to estimate expected M-quantiles
  mqo.QRLM.h <- vector("list", H)
  for(h in 1:H) {
    Y.d.pred.h <- Y.d.pred + sapply(sqrt(Y.d.Hajek$Y.var), rnorm, n = 1, mean = 0)
    
    QRLMmod.h <- lapply(tau, function(q) {
      QRLMmod(cbind(1, Xd), Y.d.pred.h, Y.d.Hajek$Y.var, q = q)
    })
    
    mqo.QRLM.h[[h]] <- matrix(gridfitinter(Y.d.Hajek$Y.dir,
                              sapply(QRLMmod.h, function(x) x[[1]]),
      tau), nrow = D, ncol = 1)
  }
  
  # Step 6: Average over H Monte Carlo replicates
  mqo.H <- Reduce('+', mqo.QRLM.h) / H
  
  # Step 7: Fit MQ SAE model at averaged tau
  QRLM.SAE.h <- lapply(mqo.H, function(q) {
    QRLMmod(cbind(1, Xd), Y.d.Hajek$Y.dir, Y.d.Hajek$Y.var, q = q)
  })
  
  QRLM.SAE.h <- sapply(QRLM.SAE.h, function(x) x[[2]])
  
  # Step 8: Compute final predictions for each area 
  pred <- sapply(1:D, function(d) cbind(1, Xd)[d, ] %*% QRLM.SAE.h[, d])
  
  return(pred)
}


############################################################################
###  A Working Likelihood for M-quantiles: The Generalised Asymmetric  #####
###                Least Informative (GALI) distribution               #####
############################################################################

# Preliminary functions
rho <- function(x, k){ 2*ifelse(abs(x) <= k, x^2/2, abs(x)*k-k^2/2) }
rho.q <- function(u, q){ ifelse(u < 0, (q-1)*u, q*u) * rho(u, k=1.345) }

### SIMULATE GALI(0,1,theta.d) RANDOM VARIABLES

simu.GALI01 <- function(theta.d){
  
  grid <- seq(-30,30,by=.01)
  integrando <- function(u, theta.d){ exp(-rho.q(u, theta.d)) }
  B <- integrate(integrando, theta.d=theta.d, lower = -Inf, upper = Inf)[[1]] 
  integrando <- function(u, theta.d, B){ exp(-rho.q(u, theta.d))/B }
  
  # GALI01.distribucion
  GALI01.distribucion <- function(x, theta.d, B){ integrate(integrando, theta.d=theta.d, 
                                  B=B, lower = -Inf, upper = x)[[1]] }
  
  GALI01.dist.emp <- sapply(grid,  GALI01.distribucion, theta.d=theta.d, B=B)
  
  index <- which(GALI01.dist.emp == max(GALI01.dist.emp))[1]
  
  GALI01.dist.emp[index:length(GALI01.dist.emp)] <- 1
  
  return(list(GALI01.dist.emp, grid))
}  

simu.GALI01.values <- function(GALI01.dist.emp, grid, L){
  
  simu.GALI01.emp <- sapply(runif(L), 
          function(x){ grid[ which.min(abs(GALI01.dist.emp-x))] } )
  
  return(simu.GALI01.emp)
}


###################################################
###  MSE ESTIMATION FROM AMQ MODELS: ANALYTICAL ###
###################################################

MQpred.mse <- function(Ydir, Xd, Y.var, covHajek, mod.area.SAE, D){
  
  n <- dim(Xd)[1]
  p <- dim(Xd)[2]
  
  AMQ<-mse1<-var1<-var2<-bias<-H<-list()
  
  for (d in 1:D){
    weig.d <- diag(mod.area.SAE$q.weights[, d])
    
    XtW   <- t(Xd) %*% weig.d
    invXtWX <- solve(XtW %*% Xd)
    ed    <- invXtWX %*% XtW
    
    # Hat matrix for domain d
    H[[d]] <- Xd %*% ed
    
    # Row d of the projection matrix
    ad <- Xd[d, ] %*% ed
    
    ones.d <- rep(0, D); ones.d[d] <- 1
      
    AMQ[[d]] <- ad  %*% Ydir
    
    var1[[d]]   <- (ad - ones.d) %*% covHajek %*% t(ad - ones.d)
    # var1[[d]]   <- (ad - ones.d) %*% diag(Y.var) %*% t(ad - ones.d)
    aux <- sum(sapply(1:D, function(g){ ad[g]*Xd[g, ]%*%mod.area.SAE$coef[, g] }))
    
    bias[[d]] <- aux - AMQ[[d]]
    
    mse1[[d]] <- var1[[d]] + Y.var[d] + 2*((ad - ones.d) %*% covHajek[d,]) + bias[[d]]^2
    
  }
    return(list(mse = unlist(mse1), H=H))
}
   
mean_upper_trim <- function(x, upper.trim, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  if (upper.trim <= 0 || upper.trim >= 1) stop("upper.trim must be between 0 and 1")
  
  q <- quantile(x, probs = 1 - upper.trim, na.rm = na.rm)
  mean(x[x <= q], na.rm = na.rm)
}


###################################################
###  MSE ESTIMATION FROM AMQ MODELS: BOOTSTRAP  ###
###################################################

MSE.boot <- function(mqo.area, mod.area.50, u.MQ.C, tau.large, GALI01, MSE.AMQ,
                     Y.dir, Xd, sd, D, R, mod.area.SAE){
  
  mqo.pos <- sapply(mqo.area, function(x) which.min(abs(tau.large - x)))
  
  GALI01.pos <- lapply(mqo.pos, function(x){ GALI01[[x]]} )
  
  pred.50 <- Xd%*%mod.area.50$coef
  
  pred.AMQ <- Y.boot <- list()
  for(r in 1:R){
    Y.boot[[r]] <- pred.50 + hub.psi(sd*sample(u.MQ.C, replace = TRUE, D), k=2)
    
    GALI.random <- list()
    for(d in 1:D){
      GALI.random[[d]] <- simu.GALI01.values(GALI01.dist.emp=GALI01.pos[[d]][[1]], 
                                   grid=GALI01.pos[[d]][[2]], D)
      
      GALI.random[[d]] <- (diag(1, D, D)-MSE.AMQ$H[[d]])%*%(Xd%*%mod.area.SAE$coef[,d] 
                                                            + sd[d] * GALI.random[[d]])
      GALI.random[[d]] <- GALI.random[[d]][d]
    }
    
    Y.dir.boot <- Y.boot[[r]] + unlist(GALI.random)
    
    mod.area.SAE.r <- QRLM(x=Xd, y=Y.dir.boot, q=mqo.area, maxit=30, k = 1.345)
    
    pred.AMQ[[r]] <- sapply(1:D, function(d) Xd[d, ] %*% mod.area.SAE.r$coef[, d])
  }
  
  mse  <- apply(sapply(1:R, function(r) (pred.AMQ[[r]] - Y.boot[[r]])^2), 
                1, mean_upper_trim, 0.05)
  
 return(list(mse))
  
}

