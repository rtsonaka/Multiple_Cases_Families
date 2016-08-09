gauher <- function (n) {
    m <- trunc((n + 1)/2)
    x <- w <- rep(-1, n)
    for (i in seq_len(m)) {
        z <- if (i == 1) {
            sqrt(2 * n + 1) - 1.85575 * (2 * n + 1)^(-0.16667)
        }
        else if (i == 2) {
            z - 1.14 * n^0.426/z
        }
        else if (i == 3) {
            1.86 * z - 0.86 * x[1]
        }
        else if (i == 4) {
            1.91 * z - 0.91 * x[2]
        }
        else {
            2 * z - x[i - 2]
        }
        for (its in seq_len(10)) {
            p1 <- 0.751125544464943
            p2 <- 0
            for (j in seq_len(n)) {
                p3 <- p2
                p2 <- p1
                p1 <- z * sqrt(2/j) * p2 - sqrt((j - 1)/j) *
                  p3
            }
            pp <- sqrt(2 * n) * p2
            z1 <- z
            z <- z1 - p1/pp
            if (abs(z - z1) <= 3e-14)
                break
        }
        x[i] <- z
        x[n + 1 - i] <- -z
        w[i] <- 2/(pp * pp)
        w[n + 1 - i] <- w[i]
    }
    list(x = x, w = w)
}

dmvnorm <- function (x, mu, Sigma, log = FALSE) {
    if (!is.matrix(x))
        x <- rbind(x)
    p <- length(mu)
    if (p == 1) {
        dnorm(x, mu, sqrt(Sigma), log = log)
    }
    else {
        t1 <- length(mu) == length(Sigma)
        t2 <- all(abs(Sigma[lower.tri(Sigma)]) < sqrt(.Machine$double.eps))
        if (t1 || t2) {
            if (!t1)
                Sigma <- diag(Sigma)
            nx <- nrow(x)
            ff <- rowSums(dnorm(x, rep(mu, each = nx), sd = rep(sqrt(Sigma),
                each = nx), log = TRUE))
            if (log)
                ff
            else exp(ff)
        }
        else {
            ed <- eigen(Sigma, symmetric = TRUE)
            ev <- ed$values
            evec <- ed$vectors
            if (!all(ev >= -1e-06 * abs(ev[1])))
                stop("'Sigma' is not positive definite")
            ss <- x - rep(mu, each = nrow(x))
            inv.Sigma <- evec %*% (t(evec)/ev)
            quad <- 0.5 * rowSums((ss %*% inv.Sigma) * ss)
            fact <- -0.5 * (p * log(2 * pi) + sum(log(ev)))
            if (log)
                as.vector(fact - quad)
            else as.vector(exp(fact - quad))
        }
    }
}

cd <- function (x, f, ..., eps = 0.001) {
    n <- length(x)
    res <- numeric(n)
    ex <- pmax(abs(x), 1)
    for (i in 1:n) {
        x1 <- x2 <- x
        x1[i] <- x[i] + eps * ex[i]
        x2[i] <- x[i] - eps * ex[i]
        diff.f <- c(f(x1, ...) - f(x2, ...))
        diff.x <- x1[i] - x2[i]
        res[i] <- diff.f/diff.x
    }
    res
}



GHpoints <- function (y, id, X, thetas, kinmat, n.mems, GHk = 5) {
    betasM <- thetas[-length(thetas)]
    sigma.b <- exp(thetas[length(thetas)])
    Sigma <- sigma.b^2 * kinmat # change dim Sigma
    if(marginalized){
        margProbs <- plogis(c(X %*% betasM))
        Delta.ij <- numeric(nrow(X))
        ind <- rep(seq_len(n.mems), length.out = nrow(X))
        for (l in 1:nrow(X)) {
            Delta.ij[l] <- Delta.ijFun(0, margProbs[l], sigma.b = sqrt(Sigma[ind[l], ind[l]])) # change dim Sigma
        }
    }else{
        Delta.ij <- c(X %*% betasM)
    }
    n <- length(unique(id))
    GH <- gauher(GHk)
    b <- as.matrix(expand.grid(rep(list(GH$x), n.mems)))
    k <- nrow(b)
    wGH <- as.matrix(expand.grid(rep(list(GH$w), n.mems)))
    wGH <- 2^(n.mems/2) * apply(wGH, 1, prod) * exp(rowSums(b * b))
    b <- sqrt(2) * b
    ###########
    fn <- function (b, y.i, delta.ij) {
        log.p.b <- dmvnorm(b, rep(0, n.mems), Sigma, log = TRUE)
        log.p.yb <- sum(dbinom(y.i, 1, plogis(delta.ij + b), log = TRUE))
        - log.p.yb - log.p.b
    }
    gr <- function (b, y.i, delta.ij) {
        cd(b, fn, y.i = y.i, delta.ij = delta.ij)
    }
    scaled.b <- vector("list", n)
    det.Bs <- numeric(n)
    for (i in seq_len(n)) {
        id.i <- id == i
        opt <- optim(rep(0, n.mems), fn, gr, y.i = y[id.i], delta.ij = Delta.ij[id.i],
            method = "BFGS", hessian = TRUE)
        invChol <- solve(chol(opt$hessian))
        scaled.b[[i]] <- t(opt$par + tcrossprod(invChol, b))
        det.Bs[i] <- det(invChol)
    }
    list(scaled.b = scaled.b, det.Bs = det.Bs, wGH = wGH, unscaled.b = b)
}



Delta.ijFun <- function (startDelta, margProb, sigma.b,
        lowerInt = -25, upperInt = 25) {
    fn.AGH <- function(startDelta) {
        f <- function (b) {
            exp(plogis(startDelta + b, log.p = TRUE) +
                dnorm(b, sd = sigma.b, log = TRUE))
        }
        margProb - integrate(f, lower = lowerInt, upper = upperInt, rel.tol = 1e-05)$value
    }
    try.n <- 1
    unr <- try(uniroot(fn.AGH, c(startDelta - 10, startDelta + 10))$root, TRUE)
    while (inherits(unr, "try-error") && try.n < 10) {
        try.n <- try.n + 1
        unr <- try(uniroot(fn.AGH, c(startDelta - 10*try.n, startDelta + 10*try.n))$root, TRUE)
    }
    unr
}
