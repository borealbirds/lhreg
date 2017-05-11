## ----eval=FALSE,results='hide'-------------------------------------------
#  devtools::install_github("borealbirds/lhreg")
#  devtools::build_vignettes("~/repos/lhreg")

## ------------------------------------------------------------------------
library(lhreg)
data(lhreg_data)
str(lhreg_data)

## ----eval=FALSE----------------------------------------------------------
#  library(ape)
#  mph <- read.nexus("11960.tre") # 1000 trees with Ericson backbone
#  cph <- consensus(mph)
#  table(sapply(mph, function(z) length(z$tip.label)))
#  CORR <- TRUE
#  vv <- list()
#  vv[[1]] <- vcv(mph[[1]], corr=CORR)
#  for (i in 2:length(mph)) {
#      v <- vcv(mph[[i]], corr=CORR)
#      v <- v[rownames(vv[[1]]), colnames(vv[[1]])]
#      vv[[i]] <- v
#  }
#  vvv <- v
#  for (i in 1:length(v)) {
#      vvv[i] <- mean(sapply(vv, function(z) z[i]))
#  }
#  spp <- intersect(rownames(lhreg_data), rownames(vvv))
#  vvv <- vvv[spp,spp]
#  cor_matrix <- as.matrix(nearPD(vvv, corr=TRUE)$mat)

## ------------------------------------------------------------------------
data(cor_matrix)
str(cor_matrix)
heatmap(cor_matrix)

## ------------------------------------------------------------------------
y <- lhreg_data$logphi
m4 <- lm(y ~ logmass + Mig2 + Nesthm + xMaxFreqkHz + Hab4, lhreg_data)
m3 <- lm(y ~ logmass + Mig2 + Nesthm + xMaxFreqkHz + Hab3, lhreg_data)
m2 <- lm(y ~ logmass + Mig2 + Nesthm + xMaxFreqkHz + Hab2, lhreg_data)
m4s <- step(m4, trace=0)
m3s <- step(m3, trace=0)
m2s <- step(m2, trace=0)
AIC(m4,m3,m2,m4s,m3s,m2s)
summary(m2s)
mp <- update(m2s, .~.-Nesthm)
mp0 <- lm(y ~ 1)
summary(mp0) # null model for SR
summary(mp) # full model for SR
## evaluating interactions
mpx <- lm(y ~ logmass + Mig2 + Nesthm + xMaxFreqkHz + Hab2 +
    logmass:xMaxFreqkHz + Hab2:xMaxFreqkHz, lhreg_data)
mpx2 <- step(mpx)
summary(mpx2)

## ------------------------------------------------------------------------
y <- lhreg_data$logtau
m4 <- lm(y ~ logmass + Mig2 + Nesthm + xMaxFreqkHz + Hab4, lhreg_data)
m3 <- lm(y ~ logmass + Mig2 + Nesthm + xMaxFreqkHz + Hab3, lhreg_data)
m2 <- lm(y ~ logmass + Mig2 + Nesthm + xMaxFreqkHz + Hab2, lhreg_data)
m4s <- step(m4, trace=0)
m3s <- step(m3, trace=0)
m2s <- step(m2, trace=0)
AIC(m4,m3,m2,m4s,m3s,m2s)
summary(m2s)
mt <- update(m2s, .~.-Nesthm)
mt0 <- lm(y ~ 1)
summary(mp0) # null model for DD
summary(mp) # full model for DD
## evaluating interactions
mtx <- lm(y ~ logmass + Mig2 + Nesthm + xMaxFreqkHz + Hab2 +
    logmass:xMaxFreqkHz + Hab2:xMaxFreqkHz, lhreg_data)
mtx2 <- step(mtx)
summary(mtx2)

## ----eval=FALSE----------------------------------------------------------
#  met <- "DE" # can use "SANN", "Nelder-Mead" is quickest
#  x <- lhreg_data
#  vc <- cor_matrix
#  
#  ## model matrix definitions
#  X0 <- matrix(1, nrow(x), 1) # null (for both)
#  colnames(X0) <- "Intercept"
#  Xt <- model.matrix(mt) # DD
#  Xp <- model.matrix(mp) # SR
#  
#  ## --- model selection
#  
#  
#  ## tau model selection
#  
#  tM0 <- lhreg(Y=x$logtau, X=X0, SE=x$logtauSE, V=vc, lambda=0,
#      hessian=TRUE, method=met)
#  tMp0 <- lhreg(Y=x$logtau, X=X0, SE=x$logtauSE, V=vc, lambda=1,
#      hessian=TRUE, method=met)
#  tMl0 <- lhreg(Y=x$logtau, X=X0, SE=x$logtauSE, V=vc, lambda=NA,
#      hessian=TRUE, method=met)
#  tMx <- lhreg(Y=x$logtau, X=Xt, SE=x$logtauSE, V=vc, lambda=0,
#      hessian=TRUE, method=met)
#  tMpx <- lhreg(Y=x$logtau, X=Xt, SE=x$logtauSE, V=vc, lambda=1,
#      hessian=TRUE, method=met)
#  tMlx <- lhreg(Y=x$logtau, X=Xt, SE=x$logtauSE, V=vc, lambda=NA,
#      hessian=TRUE, method=met)
#  
#  ## phi model selection
#  
#  pM0 <- lhreg(Y=x$logphi, X=X0, SE=x$logphiSE, V=vc, lambda=0,
#      hessian=TRUE, method=met)
#  pMp0 <- lhreg(Y=x$logphi, X=X0, SE=x$logphiSE, V=vc, lambda=1,
#      hessian=TRUE, method=met)
#  pMl0 <- lhreg(Y=x$logphi, X=X0, SE=x$logphiSE, V=vc, lambda=NA,
#      hessian=TRUE, method=met)
#  pMx <- lhreg(Y=x$logphi, X=Xp, SE=x$logphiSE, V=vc, lambda=0,
#      hessian=TRUE, method=met)
#  pMpx <- lhreg(Y=x$logphi, X=Xp, SE=x$logphiSE, V=vc, lambda=1,
#      hessian=TRUE, method=met)
#  pMlx <- lhreg(Y=x$logphi, X=Xp, SE=x$logphiSE, V=vc, lambda=NA,
#      hessian=TRUE, method=met)
#  
#  aict <- AIC(tM0, tMp0, tMl0, tMx, tMpx, tMlx)
#  aict$dAIC <- aict$AIC-min(aict$AIC)
#  
#  aicp <- AIC(pM0, pMp0, pMl0, pMx, pMpx, pMlx)
#  aicp$dAIC <- aicp$AIC-min(aicp$AIC)
#  
#  cbind(aict, t(sapply(list(tM0, tMp0, tMl0, tMx, tMpx, tMlx),
#      function(z) z$summary[c("sigma","lambda"), 1])))
#  
#  cbind(aicp, t(sapply(list(pM0, pMp0, pMl0, pMx, pMpx, pMlx),
#      function(z) z$summary[c("sigma","lambda"), 1])))

## ----eval=FALSE----------------------------------------------------------
#  library(parallel)
#  
#  lam <- seq(0, 2, by=0.01)
#  
#  cl <- makeCluster(4)
#  
#  tmp <- clusterEvalQ(cl, library(lhreg))
#  
#  object <- tMl0
#  clusterExport(cl, "object")
#  pl_tMl0 <- pbsapply(lam, function(z, ...) profile_lambda1(object, z, ...), cl=cl,
#      method=met)
#  
#  object <- tMlx
#  clusterExport(cl, "object")
#  pl_tMlx <- pbsapply(lam, function(z, ...) profile_lambda1(object, z, ...), cl=cl,
#      method=met)
#  
#  object <- pMl0
#  clusterExport(cl, "object")
#  pl_pMl0 <- pbsapply(lam, function(z, ...) profile_lambda1(object, z, ...), cl=cl,
#      method=met)
#  
#  object <- pMlx
#  clusterExport(cl, "object")
#  pl_pMlx <- pbsapply(lam, function(z, ...) profile_lambda1(object, z, ...), cl=cl,
#      method=met)
#  
#  stopCluster(cl)

## ----eval=FALSE----------------------------------------------------------
#  ## makes sense to use lm for EDR LOO
#  tM0$coef
#  c(coef(mt0), log(summary(mt0)$sigma))
#  tMx$coef
#  c(coef(mt), log(summary(mt)$sigma))
#  
#  pM0$coef
#  c(coef(mp0), log(summary(mp0)$sigma))
#  pMx$coef
#  c(coef(mp), log(summary(mp)$sigma))
#  
#  n <- nrow(x)
#  cl <- makeCluster(4)
#  
#  pr_tM0  <- t(pbsapply(1:n, loo1, mm=tM0))
#  pr_tMp0 <- t(pbsapply(1:n, loo2, mm=tMp0))
#  pr_tMl0 <- t(pbsapply(1:n, loo2, mm=tMl0))
#  pr_tMx  <- t(pbsapply(1:n, loo1, mm=tMx))
#  pr_tMpx <- t(pbsapply(1:n, loo2, mm=tMpx))
#  pr_tMlx <- t(pbsapply(1:n, loo2, mm=tMlx))
#  
#  pr_pM0  <- t(pbsapply(1:n, loo1, mm=pM0))
#  pr_pMp0 <- t(pbsapply(1:n, loo2, mm=pMp0))
#  pr_pMl0 <- t(pbsapply(1:n, loo2, mm=pMl0))
#  pr_pMx  <- t(pbsapply(1:n, loo1, mm=pMx))
#  pr_pMpx <- t(pbsapply(1:n, loo2, mm=pMpx))
#  pr_pMlx <- t(pbsapply(1:n, loo2, mm=pMlx))
#  
#  stopCluster(cl)
#  
#  SSEt <- c(
#      tM0 = sum((pr_tM0 - tM0$Y)^2),
#      tMp0 = sum((pr_tMp0 - tMp0$Y)^2),
#      tMl0 = sum((pr_tMl0 - tMl0$Y)^2),
#      tMx = sum((pr_tMx - tMx$Y)^2),
#      tMpx = sum((pr_tMpx - tMpx$Y)^2),
#      tMlx = sum((pr_tMlx - tMlx$Y)^2))
#  SSEp <- c(
#      pM0 = sum((pr_pM0 - pM0$Y)^2),
#      pMp0 = sum((pr_pMp0 - pMp0$Y)^2),
#      pMl0 = sum((pr_pMl0 - pMl0$Y)^2),
#      pMx = sum((pr_pMx - pMx$Y)^2),
#      pMpx = sum((pr_pMpx - pMpx$Y)^2),
#      pMlx = sum((pr_pMlx - pMlx$Y)^2))
#  MSEt <- SSEt / n
#  MSEp <- SSEp / n

## ------------------------------------------------------------------------
#load(system.file("extdata", "lhreg-results-DE.Rdata", package = "lhreg"))

