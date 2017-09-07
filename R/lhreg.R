lhreg <-
function(Y, X, SE, V, init=NULL, lambda=NA, method="Nelder-Mead",
hessian=FALSE, DElimit=10, eval=FALSE)
{
    .solvenear <-
    function(x)
    {
        xinv <- try(solve(x), silent = TRUE)
        if (inherits(xinv, "try-error"))
            xinv <- as.matrix(solve(Matrix::nearPD(x)$mat))
        xinv
    }

    colnames(X) <- gsub("[[:punct:]]", "", colnames(X))
    if (is.null(init)) {
        m <- lm(Y ~ .-1, data.frame(X))
        init <- c(coef(m), log_sigma=log(summary(m)$sigma)) # lambda is fixed
        if (is.na(lambda))
            init <- c(init, log_lambda=0) # lambda is optimized over
    }
    n <- length(Y)
    nn <- seq_len(n)
    D <- diag(1, n, n)
    diag(D) <- SE^2
    nllfun <- function(par) {
        cf <- par[1:ncol(X)]
        sigma_sq <- exp(par[ncol(X)+1])^2
        if (is.na(lambda)) {
            #lambda <- plogis(par[ncol(X)+2]) # lambda in 0-1
            lambda <- exp(par[ncol(X)+2]) # real valued lambda
        }
        VV <- V
        VV[lower.tri(VV)] <- lambda * V[lower.tri(VV)]
        VV[upper.tri(VV)] <- lambda * V[upper.tri(VV)]
        mu <- drop(X %*% cf)
        Sigma <- sigma_sq * VV + D
        out <- -sum(mvtnorm::dmvnorm(Y, mu, Sigma, log=TRUE))
        if (is.na(out) || is.infinite(out))
            10^15 else out
    }
    if (eval)
        return(nllfun(init))
    if (method == "DE") {
        up <- rep(DElimit, length(init))
        lo <- -up
        opt <- DEoptim(fn=nllfun, lower=lo, upper=up,
            control=list(trace=FALSE, itermax=length(init)*200))
        cf <- opt$optim$bestmem
        ll <- -opt$optim$bestval
        if (hessian) {
            hess <- optimHess(opt$optim$bestmem, nllfun)
            S <- .solvenear(hess)
        } else {
            S <- NULL
        }
    } else {
        opt <- optim(init, nllfun, hessian=hessian, method=method)
        cf <- opt$par
        ll <- -opt$value
        S <- if (hessian)
            .solvenear(opt$hessian) else NULL
    }
    if (hessian) {
        mvn <- rmvnorm(10^4, cf, S)
        mvn[,ncol(X)+1] <- exp(mvn[,ncol(X)+1])
        if (is.na(lambda)) {
            mvn[,ncol(X)+2] <- exp(mvn[,ncol(X)+2])#plogis(mvn[,ncol(X)+2])
        }
        trcf <- colMeans(mvn)
        trse <- apply(mvn, 2, sd)
        if (!is.na(lambda)) {
            trcf <- c(trcf, lambda)
            trse <- c(trse, NA)
        }
    } else {
        trcf <- cf
        trcf[ncol(X)+1] <- exp(trcf[ncol(X)+1])
        if (is.na(lambda)) {
            trcf[ncol(X)+2] <- exp(trcf[ncol(X)+2])#plogis(trcf[ncol(X)+2])
        } else {
            trcf[ncol(X)+2] <- lambda
        }
        trse <- rep(NA, length(trcf))
    }
    tstat <- trcf/trse
    pval <- 2 * pnorm(-abs(tstat))
    coefs <- cbind(trcf, trse, tstat, pval)
    colnames(coefs) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(coefs) <- c(paste0("beta_", colnames(X)), "sigma", "lambda")
    out <- list(coef=cf, loglik=ll, vcov=S, nobs=length(Y), df=length(cf),
        summary=coefs,
        Y=Y, X=X, SE=SE, V=V, init=init, method=method, lambda=lambda)
    class(out) <- "lhreg"
    out
}

logLik.lhreg <-
function (object, ...)
{
    structure(object$loglik, df = object$df, class = "logLik")
}

summary.lhreg <-
function(object, ...)
{
    printCoefmat(object$summary, signif.legend = TRUE, ...)
}

profile_lambda1 <-
function(object, value, ...)
{
    if (!is.na(object$lambda))
        stop("lambda must be fixed!")
    logLik(lhreg(Y=object$Y, X=object$X, SE=object$SE, V=object$V,
        lambda=value, hessian=FALSE, init=object$init, ...))
}

loo1 <-
function(i, object, return_coefs=TRUE)
{
    yy <- object$Y[-i]
    xx <- object$X[-i,,drop=FALSE]
    mod <- lm(yy ~ xx-1)
    pr <- drop(object$X[i,,drop=FALSE] %*% coef(mod))
    out <- c(obs=unname(object$Y[i]), pred=unname(pr))
    if (return_coefs)
        out <- c(out, est=c(coef(mod), sigma=summary(mt0)$sigma))
    out
}

loo2 <-
function(i, object, return_coefs=TRUE, method=NULL)
{
    if (is.null(method))
        method <- object$method
    remod <- lhreg(Y=object$Y[-i], X=object$X[-i,,drop=FALSE], SE=object$SE[-i],
        V=object$V[-i,-i,drop=FALSE], lambda=object$lambda, hessian=FALSE,
        method=method, init=object$coef)

    VV <- remod$summary["sigma",1]^2 * object$V
    VV[lower.tri(VV)] <- remod$summary["lambda",1] * VV[lower.tri(VV)]
    VV[upper.tri(VV)] <- remod$summary["lambda",1] * VV[upper.tri(VV)]

    n <- length(object$Y)
    y <- object$Y
    mu <- drop(object$X %*% object$coef[1:ncol(object$X)])

    y2 <- y[-i]
    mu1 <- mu[i]
    mu2 <- mu[-i]

    #Sig11 <- VV[i,i,drop=FALSE]
    Sig12 <- VV[i,-i,drop=FALSE]
    #Sig21 <- VV[-i,i,drop=FALSE]
    Sig22 <- VV[-i,-i,drop=FALSE]

    ## interested in mu12 only (due to observation error)
    mu12 <- drop(mu1 + Sig12 %*% solve(Sig22) %*% (y2-mu2))
    out <- c(obs=unname(y[i]), pred=unname(mu12))
    if (return_coefs)
        out <- c(out, est=remod$summary[,1])
    out
}

simulate.lhreg <-
function(object, nsim = 1, seed = NULL, lambda = NA, obs_error = FALSE, ...)
{
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    if (nsim < 1)
        stop("'nsim' must be at least 1")
    lam <- if (is.na(lambda))
        object$summary["lambda",1] else lambda
    VV <- object$summary["sigma",1]^2 * object$V
    VV[lower.tri(VV)] <- lam * VV[lower.tri(VV)]
    VV[upper.tri(VV)] <- lam * VV[upper.tri(VV)]
    if (obs_error)
        diag(VV) <- diag(VV) + object$SE
    mu <- drop(object$X %*% object$coef[1:ncol(object$X)])
    t(mvtnorm::rmvnorm(nsim, mu, VV, ...))
}
