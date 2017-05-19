## ----update-stuff,eval=FALSE,results='hide',echo=FALSE-------------------
#  devtools::install_github("borealbirds/lhreg")
#  devtools::build_vignettes("~/repos/lhreg")

## ----install,eval=FALSE--------------------------------------------------
#  devtools::install_github("borealbirds/lhreg")

## ----vignette,eval=FALSE-------------------------------------------------
#  vignette(topic = "lhreg", package = "lhreg")

## ----trait-data,message=FALSE--------------------------------------------
library(lhreg)
data(lhreg_data)
str(lhreg_data)
with(lhreg_data, plot(exp(logphi), exp(logtau),
    cex=logmass*0.5, col=Mig2, pch=c(21, 22)[Hab2]))
legend("topright", bty="n", pch=c(21, 21, 22, 22), col=c(1,2,1,2),
    legend=c("Migratory/Closed", "Resident/Closed",
    "Migratory/Open", "Resident/Open"))

## ----phylo-corr,eval=FALSE-----------------------------------------------
#  library(ape)
#  mph <- read.nexus("11960.tre") # 1000 trees with Ericson backbone
#  lhreg_tree <- consensus(mph)
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

## ----heatmap-------------------------------------------------------------
library(ape)
data(cor_matrix)
str(cor_matrix)
heatmap(cor_matrix)

## ----phylo-tree,fig.height=12,fig.width=7--------------------------------
plot(compute.brlen(vcv2phylo(cor_matrix)), cex=0.5)

## ----screening-sr--------------------------------------------------------
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

amp <- anova(mp)
round(structure(100 * amp[,"Sum Sq"] / sum(amp[,"Sum Sq"]),
    names=rownames(amp)), 1)

mp0 <- lm(y ~ 1)
summary(mp0) # null model for SR
summary(mp) # full model for SR
## evaluating interactions
mpx <- lm(y ~ logmass + Mig2 + Nesthm + xMaxFreqkHz + Hab2 +
    logmass:xMaxFreqkHz + Hab2:xMaxFreqkHz, lhreg_data)
mpx2 <- step(mpx)
summary(mpx2)

## ----screening-dd--------------------------------------------------------
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

amt <- anova(mt)
round(structure(100 * amt[,"Sum Sq"] / sum(amt[,"Sum Sq"]),
    names=rownames(amt)), 1)

mt0 <- lm(y ~ 1)
summary(mp0) # null model for DD
summary(mp) # full model for DD
## evaluating interactions
mtx <- lm(y ~ logmass + Mig2 + Nesthm + xMaxFreqkHz + Hab2 +
    logmass:xMaxFreqkHz + Hab2:xMaxFreqkHz, lhreg_data)
mtx2 <- step(mtx)
summary(mtx2)

## ----models,eval=FALSE---------------------------------------------------
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

## ----profile-lik,eval=FALSE----------------------------------------------
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

## ----loo,eval=FALSE------------------------------------------------------
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

## ----load-results--------------------------------------------------------
load(system.file("extdata", "lhreg-results-DE.rda", package = "lhreg"))

Crit <- -0.5*qchisq(0.95, 1)
ltmp <- seq(0, 1, by=0.0001)
## red-yl-blue
Col <- c("#2C7BB6", "#6BAACF", "#ABD9E9", "#D4ECD3", "#FFFFBF", "#FED690",
    "#FDAE61", "#EA633E", "#D7191C")
Col1 <- Col[1]
Col2 <- rgb(171/255, 217/255, 233/255, 0.5) # Col[3]
Col3 <- Col[9]
Col4 <- rgb(253/255, 174/255, 97/255, 0.5) # Col[7]
prt <- exp(pr_tMx)
prp <- exp(pr_pMlx)

## ----table-1-------------------------------------------------------------
ff <- function(z) {
    zz <- z$summary

    pcut <- function(p) {
        factor(c("***", "**", "*", "+", "ns")[as.integer(cut(p,
            c(1, 0.1, 0.05, 0.01, 0.001, 0), include.lowest=TRUE, right=FALSE))],
            levels=c("***", "**", "*", "+", "ns"))
    }
    structure(sapply(1:nrow(zz), function(i)
        paste0(round(zz[i,1], 3), " (SE +/- ", round(zz[i,2], 3), pcut(zz[i,4]), ")")),
        names=rownames(zz))
}
zzz <- lapply(list(pM0, pMp0, pMl0, pMx, pMpx, pMlx), ff)
m <- matrix("", length(zzz[[6]]), 6)
rownames(m) <- names(zzz[[6]])
colnames(m) <- c("M0", "Mp0", "Ml0", "Mx", "Mpx", "Mlx")
for (i in 1:6) {
    j <- match(names(zzz[[i]]), rownames(m))
    m[j,i] <- zzz[[i]]
}
ii <- grep("(SE +/- NANA)", m,fixed=TRUE)
m[ii] <- gsub("(SE +/- NANA)", "(fixed)", m[ii], fixed=TRUE)
m[m==""] <- "n/a"
m1 <- m

zzz <- lapply(list(tM0, tMp0, tMl0, tMx, tMpx, tMlx), ff)
m <- matrix("", length(zzz[[6]]), 6)
rownames(m) <- names(zzz[[6]])
colnames(m) <- c("M0", "Mp0", "Ml0", "Mx", "Mpx", "Mlx")
for (i in 1:6) {
    j <- match(names(zzz[[i]]), rownames(m))
    m[j,i] <- zzz[[i]]
}
ii <- grep("(SE +/- NANA)", m,fixed=TRUE)
m[ii] <- gsub("(SE +/- NANA)", "(fixed)", m[ii], fixed=TRUE)
m[m==""] <- "n/a"
m2 <- m

m1 <- rbind(m1, df=aicp$df, dAIC=round(aicp$dAIC, 3),
    XV_MSE=round(MSEp, 3))
m2 <- rbind(m2, df=aict$df, dAIC=round(aict$dAIC, 3),
    XV_MSE=round(MSEt, 3))

print.default(m1, quote=FALSE) # SR results
print.default(m2, quote=FALSE) # DD results

## ----fig-1,width=14,height=6---------------------------------------------
#pdf("Fig1.pdf", width=14, height=6)
op <- par(mfrow=c(1,2))

Res1 <- pl_pMl0
Res2 <- pl_pMlx
yv <- exp(Res1-max(Res1))
plot(lam, yv, type="n", lwd=1, ylim=exp(c(2*Crit, 0)),
    xlab=expression(lambda), ylab="Profile Likelihood Ratio")
tmp1 <- splinefun(lam, yv)(ltmp)
i1 <- tmp1 > exp(Crit)
yv <- exp(Res2-max(Res2))
tmp2 <- splinefun(lam, yv)(ltmp)
i2 <- tmp2 > exp(Crit)
polygon(c(ltmp[i1], rev(ltmp[i1])),
    c(tmp1[i1], rep(-1, sum(i1))),
    col=Col2, border=NA)
polygon(c(ltmp[i2], rev(ltmp[i2])),
    c(tmp2[i2], rep(-1, sum(i2))),
    col=Col4, border=NA)
abline(h=exp(Crit), col="grey")
lines(ltmp, tmp1, col=Col1, lwd=1)
lines(ltmp[i1], tmp1[i1], col=Col1, lwd=3)
lines(rep(ltmp[which.max(tmp1)], 2), c(1, -1), col=Col1)
lines(ltmp, tmp2, col=Col3, lwd=1)
lines(ltmp[i2], tmp2[i2], col=Col3, lwd=3)
lines(rep(ltmp[which.max(tmp2)], 2), c(1, -1), col=Col3)
box()
legend(0.1, 1, bty="n", border=c(Col1, Col3), fill=c(Col2, Col4), pt.cex=2, pt.lwd=2,
    legend=c(expression(M[lambda*0]-SR), expression(M[lambda*beta]-SR)))
round(c(Max_Ml0=ltmp[which.max(tmp1)], CI=range(ltmp[i1])), 3)
round(c(Max_Mlx=ltmp[which.max(tmp2)], CI=range(ltmp[i2])), 3)

Res1 <- pl_tMl0
Res2 <- pl_tMlx
yv <- exp(Res1-max(Res1))
plot(lam, yv, type="n", lwd=1, ylim=exp(c(2*Crit, 0)),
    xlab=expression(lambda), ylab="Profile Likelihood Ratio")
tmp1 <- splinefun(lam, yv)(ltmp)
i1 <- tmp1 > exp(Crit)
yv <- exp(Res2-max(Res2))
tmp2 <- splinefun(lam, yv)(ltmp)
i2 <- tmp2 > exp(Crit)
polygon(c(ltmp[i1], rev(ltmp[i1])),
    c(tmp1[i1], rep(-1, sum(i1))),
    col=Col2, border=NA)
polygon(c(ltmp[i2], rev(ltmp[i2])),
    c(tmp2[i2], rep(-1, sum(i2))),
    col=Col4, border=NA)
abline(h=exp(Crit), col="grey")
lines(ltmp, tmp1, col=Col1, lwd=1)
lines(ltmp[i1], tmp1[i1], col=Col1, lwd=3)
lines(rep(ltmp[which.max(tmp1)], 2), c(1, -1), col=Col1)
lines(ltmp, tmp2, col=Col3, lwd=1)
lines(ltmp[i2], tmp2[i2], col=Col3, lwd=3)
lines(rep(ltmp[which.max(tmp2)], 2), c(1, -1), col=Col3)
box()
legend(0.1, 1, bty="n", border=c(Col1, Col3), fill=c(Col2, Col4), pt.cex=2, pt.lwd=2,
    legend=c(expression(M[lambda*0]-DD), expression(M[lambda*beta]-DD)))
round(c(Max_Ml0=ltmp[which.max(tmp1)], CI=range(ltmp[i1])), 3)
round(c(Max_Mlx=ltmp[which.max(tmp2)], CI=range(ltmp[i2])), 3)

par(op)

#dev.off()

## ----fig-2,width=14,height=7---------------------------------------------
#pdf("Fig2.pdf", width=14, height=7)

op <- par(mfrow=c(1,2))

Max <- 0.7
plot(prp, xaxs = "i", yaxs = "i", type="n",
    ylim=c(0, Max), xlim=c(0, Max),
    xlab="Time-removal SR Estimate", ylab="LOO SR Estimate")
abline(0,1,lty=1, col=1)
abline(0,2,lty=2, col="grey")
abline(0,1/2,lty=2, col="grey")
abline(0,1.5,lty=2, col="grey")
abline(0,1/1.5,lty=2, col="grey")
points(prp, xaxs = "i", yaxs = "i",
    col=c(Col1,Col3)[as.integer(x$Mig2)],
    cex=0.2+2*x$MaxFreqkHz/10)
box()
legend("topleft", pch=21, col=c(Col1,Col3,1,1),
    pt.cex=c(1.5,1.5,2,1), bty="n",
    legend=c("Migrant", "Resident", "Song Picth: High", "Song Pitch: Low"))
text(0.9*Max, 0.05*Max, expression(M[lambda*beta]-SR))

Max <- 2.1
plot(prt, xaxs = "i", yaxs = "i", type="n",
    ylim=c(0, Max), xlim=c(0, Max),
    xlab="Distance-sampling DD Estimate", ylab="LOO DD Estimate")
abline(0,1,lty=1, col=1)
abline(0,2,lty=2, col="grey")
abline(0,1/2,lty=2, col="grey")
abline(0,1.5,lty=2, col="grey")
abline(0,1/1.5,lty=2, col="grey")
points(prt, xaxs = "i", yaxs = "i",
    cex=0.2+2*x$logmass/5,
    col=c(Col1,Col3)[as.integer(x$Hab2)])
box()
legend("topleft", pch=21, col=c(Col1,Col3,1,1), pt.cex=c(1.5,1.5,2,1), bty="n",
    legend=c("Habitat: Closed", "Habitat: Open", "Body Mass: Large", "Body Mass: Small"))
text(0.9*Max, 0.05*Max, expression(M[0*beta]-DD))

par(op)

#dev.off()

