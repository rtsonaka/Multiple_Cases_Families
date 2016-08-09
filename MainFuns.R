logLikFam <- function (thetas, y, yPatterns, X, kinmat, n.mems, id, GHs, marginalized = FALSE,
                       ascert = FALSE, lambda, truePrev, sig.x, MAF, pseudo = FALSE,
                       penalty.glmm = FALSE, penalty.glmmN = FALSE) {
  # parameters
  betasM <- thetas[-length(thetas)]
  sigma.b <- exp(thetas[length(thetas)])
  Sigma <- sigma.b^2 * kinmat
  
  GHs <- GHpoints(y, id, X, c(betasM, log(sigma.b)), kinmat, n.mems = n.mems, GHk = 15)
  
  #####
  # calculate Dij
  if (marginalized) {
    margProbs <- plogis(c(X %*% betasM))
    Delta.ij <- numeric(nrow(X))
    ind <- rep(seq_len(n.mems), length.out = nrow(X))
    for (k in 1:nrow(X)) {
      Delta.ij[k] <- Delta.ijFun(0, margProbs[k], sigma.b = sqrt(Sigma[ind[k], ind[k]])) # change dim Sigma
    }
  } else {
    Delta.ij <- c(X %*% betasM)
  }
  ####
  # calculate log Lik
  log.p.b <- if(pseudo){
    t(sapply(GHs$scaled.b, dmvnorm, mu = rep(0, n.mems), Sigma = Sigma, log = TRUE))
  }else{
    t(dmvnorm(GHs$unscaled.b, mu = rep(0, n.mems), Sigma = Sigma, log = TRUE))
  }
  probs <- if(pseudo){
    plogis(Delta.ij + do.call(rbind, lapply(GHs$scaled.b, t)))
  }else{
    plogis(Delta.ij + t(GHs$unscaled.b)[rep(1:n.mems, n), ])
  }
  ff <- function (yy) {
    log.p.y.b <- rowsum(dbinom(yy, 1, probs, log = TRUE), id)
    if(pseudo){
      c((GHs$det.Bs * exp(log.p.y.b + log.p.b)) %*% GHs$wGH)
    }else{
      c((exp(log.p.y.b + log.p.b[rep(1, n), ])) %*% GHs$wGH)
    }
  }
  p.y <- ff(y)
  p.yPatterns <- rowSums(apply(yPatterns, 2, ff), na.rm = TRUE)
  p.yPatterns <- pmax(pmin(0.9999, p.yPatterns), 0.0001)
  ####
  # Penalty
  if(penalty.glmm){
    
    intgrAge <- function (SNP, sig.x) {
      f <- function (age, SNP, sig.x) {
        exp(plogis(c(cbind(1, SNP, age) %*% betasM/sqrt(1 + 0.346 * sigma.b^2)), log.p = TRUE)
            + dnorm(age, mean = 0, sd = sig.x, log = TRUE))
      }
      integrate(f, lower = -20, upper = 20, SNP = SNP, sig.x = sig.x)$value
    }
    modelPrev <- sum(c((1-MAF)^2, 2*MAF*(1-MAF), MAF^2) * sapply(0:2, intgrAge, sig.x = sig.x))
    
  }else{
    
    intgrAge <- function (SNP, sig.x) {
      f <- function (age, SNP, sig.x) {
        exp(plogis(c(cbind(1, SNP, age) %*% betasM), log.p = TRUE)
            + dnorm(age, mean = 0, sd = sig.x, log = TRUE))
      }
      integrate(f, lower = -20, upper = 20, SNP = SNP, sig.x = sig.x)$value
    }
    modelPrev <- sum(c((1-MAF)^2, 2*MAF*(1-MAF), MAF^2) * sapply(0:2, intgrAge, sig.x = sig.x))
  }
  
  if(penalty.glmmN){
    
    GHk <- 21
    GH <- gauher(GHk)
    b <- as.matrix(expand.grid(rep(list(GH$x), 1)))
    
    k <- nrow(b)
    wGH <- as.matrix(expand.grid(rep(list(GH$w), 1)))
    wGH <- 2^(n.mems/2) * apply(wGH, 1, prod) * exp(rowSums(b * b))
    b <- sqrt(2) * b
    integr.grid <- expand.grid(b, 0:2)
    
    modelPrev <- sum(rowsum(cbind(apply(integr.grid, 1, function(x){
      
      f <- function (age, SNP, betasM, sig.x, b) {
        
        exp(plogis(c(cbind(1, SNP, age) %*% betasM) + b, log.p = TRUE) +
              dnorm(age, mean = 0, sd = sig.x, log = TRUE))
        
      }
      integrate(f, lower = -25, upper = 25, SNP = x[2], betasM = betasM, sig.x = sig.x, b = x[1])$value
      
    }) * (dnorm(b, sd = sigma.b) * wGH)[rep(1:GHk, 3)]), rep(1:3, each = GHk)) * c((1-MAF)^2, 2 * MAF * (1-MAF), MAF^2))
    
    
    #print(modelPrev)
  }
  
  penalty <- lambda * (qlogis(truePrev) - qlogis(modelPrev))^2
  #penalty <- lambda * (log(trueHerit) - log(sigma.b))^2
  
  ####
  #cat("\n", c(betasM, sigma.b))
  if(!ascert) {
    -sum(log(p.y), na.rm = TRUE) + penalty
  } else {
    -sum(log(p.y) - log(1 - p.yPatterns), na.rm = TRUE) + penalty
  }
}





scoreFam <- function(thetas, y, yPatterns, X, kinmat, n.mems, id, GHs, marginalized, 
                     ascert, lambda, truePrev, sig.x, MAF, pseudo, penalty.glmm, penalty.glmmN) {
  cd(thetas, logLikFam, y = y, yPatterns = yPatterns, X = X, kinmat = kinmat, n.mems = n.mems,
     id = id, GHs = GHs, marginalized = marginalized, ascert = ascert, lambda = lambda, 
     truePrev = truePrev, sig.x = sig.x, MAF = MAF, pseudo = pseudo,
     penalty.glmm = penalty.glmm, penalty.glmmN = penalty.glmmN)
}



#logLikFam(initThet, y, as.matrix(rep(0, length(y))), X, kinmat, n.mems, id, GHs, marginalized = TRUE,
#    lambda = 2, truePrev = 0.2, sig.x = 1, MAF = 0.2)




Hessian.logLikFam <- function(thetas, y, yPatterns, X, kinmat, n.mems, id, GHs, 
                              marginalized, ascert, lambda, truePrev, sig.x, 
                              MAF, pseudo, penalty.glmm, penalty.glmmN){
  JM:::cd.vec(thetas, f = scoreFam, y, yPatterns, X, kinmat, n.mems, 
              id, GHs, marginalized, ascert, lambda, 
              truePrev, sig.x, MAF, pseudo, penalty.glmm, penalty.glmmN)
}



if(FALSE){
penalty.glmm <- FALSE
penalty.glmmN <- FALSE


if(penalty.glmm){
  
  intgrAge <- function (SNP, sig.x, betasM) {
    f <- function (age, SNP, sig.x) {
      exp(plogis(c(cbind(1, SNP, age) %*% betasM), log.p = TRUE)
          + dnorm(age, mean = 0, sd = sig.x, log = TRUE))
    }
    integrate(f, lower = -20, upper = 20, SNP = SNP, sig.x = sig.x)$value
  }
  truePrev <- sum(c((1-MAF)^2, 2 * MAF * (1-MAF), MAF^2) * sapply(0:2, intgrAge, sig.x = sig.x, betasM = betasM.0/sqrt(1+0.346*sigma.b^2)))
  print(truePrev)
  
  
}else{
  intgrAge <- function (SNP, sig.x, betasM) {
    f <- function (age, SNP, sig.x) {
      exp(plogis(c(cbind(1, SNP, age) %*% betasM), log.p = TRUE)
          + dnorm(age, mean = 0, sd = sig.x, log = TRUE))
    }
    integrate(f, lower = -20, upper = 20, SNP = SNP, sig.x = sig.x)$value
  }
  truePrev <- sum(c((1-MAF)^2, 2 * MAF * (1-MAF), MAF^2) * sapply(0:2, intgrAge, sig.x = sig.x, betasM = betasM.0))
  print(truePrev)
}



if(penalty.glmmN){
  GHk <- 21
  GH <- gauher(GHk)
  b <- as.matrix(expand.grid(rep(list(GH$x), 1)))
  
  k <- nrow(b)
  wGH <- as.matrix(expand.grid(rep(list(GH$w), 1)))
  wGH <- 2^(n.mems/2) * apply(wGH, 1, prod) * exp(rowSums(b * b))
  b <- sqrt(2) * b
  integr.grid <- expand.grid(b, 0:2)
  
  truePrev. <- function(betasM, sigma.b){
    
    sum(rowsum(cbind(apply(integr.grid, 1, function(x){
      
      f <- function (age, SNP, betasM, sig.x, b) {
        
        exp(plogis(c(cbind(1, SNP, age) %*% betasM) + b, log.p = TRUE) +
              dnorm(age, mean = 0, sd = sig.x, log = TRUE))
        
      }
      integrate(f, lower = -25, upper = 25, SNP = x[2], betasM = betasM, sig.x = sig.x, b = x[1])$value
      
    }) * (dnorm(b, sd = sigma.b) * wGH)[rep(1:GHk, 3)]), rep(1:3, each = GHk)) * c((1-MAF)^2, 2 * MAF * (1-MAF), MAF^2))
    
  }
  truePrev <- truePrev.(betasM, sigma.b)
}
}