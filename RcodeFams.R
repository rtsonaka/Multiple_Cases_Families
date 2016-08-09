#########################################################################################################
# Author: Roula Tsonaka (s.tsonaka@lumc.nl)
#         Leiden University Medical Center (LUMC)
#########################################################################################################
# Title: RcodeFams.R
# Aim: R code to analyse case-control family data with various kinship matrices. 
#      Thus it can handle twins, sibships, etc.
# Notes: 
#      This code has been used to run the simulation study descibed in Tsonaka et.al. (2015).
#      The code runs only for families with equal size. 
#      Results can be derived for: (1) the marginalized mixed-effects models (marginalized = TRUE),  
#                                  (2) mixed-effects logistic regression (marginalized = FALSE),
#                                  (3) without ascertainment, i.e. randomly sampled families (ascert = FALSE),
#                                  (4) with ascertainment (i.e. at least 1 or at least 2 cases, ascert = TRUE).
#      
# Reference: Tsonaka, R., van der Woude, D. and Houwing-Duistermaat, J. (2015). 
#            Marginal genetic effects estimation in family and twin studies using random-effects models. 
#            Biometrics, 71, 1130 - 1138.
#
# Date: 04AUGUST2016
###########################################################################################################


# Load main and supporting functions
source("./SupportFuns.R")
source("./MainFuns.R")

# Load data
data. <- read.table("C:/RoulaTsonaka/Papers/MargzedFamilies2/GIT/MargFams/Dataset.txt")


asc.size <- 1

# Prepare function initialization
n <- length(unique(data.$id))
n.mems <- table(data.$id)[1]

# kinship matrix - here for siblings
if(n.mems > 2){
  kin.mat <- matrix(1/2, ncol = n.mems, nrow = n.mems)
  diag(kin.mat) <- rep(1, n.mems)
  #kin.mat[c(1:2), c(1:2)] <- diag(2) 
}
if(n.mems == 2){
  kin.mat <- matrix(1/2, ncol = n.mems, nrow = n.mems)
  diag(kin.mat) <- rep(1, n.mems)
  #kin.mat[c(1:2), c(1:2)] <- diag(2)
}

# it has to be adjusted for diff ascertainment - here for 1 or 2
pos.pats <- function(x, asc.size = 2){
  pat.len <- length(x)
  x1 <- rep(0, pat.len)
  xx <- expand.grid(lapply(rep(1, pat.len), seq, from = 0))
  x2 <- xx[which(apply(xx, 1, function(x) sum(x) == 1)), ]
  x3 <- xx[which(apply(xx, 1, function(x) sum(x) == 2)), ]
  if(asc.size == 1) as.matrix(x1)
  if(asc.size == 2) rbind(x1, x2)
  #if(asc.size == 3) rbind(x1, x2, x3)
}
y <- data.$Resp
id <- data.$id
pats. <- split(y, id)
pats. <- lapply(pats., function(x) t(pos.pats(x, asc.size = 2)))

yPatterns <- if(asc.size == 1) {
  as.matrix(rep(0, length(y)))
} else {
  pats. <- do.call("rbind", pats.)
  dimnames(pats.) <- NULL
  pats.
}




# Example
#--------

# Simulated data can be found in Dataset.txt
# n = 48 multiple-cases sibships of size 3 have been simulated with at least 1 case.
# true coefficients: c(-2.5, 0, -0.5); intercept, SNP, age
# true genetic variance: 3
# true prevalence 0.09

# Design matrix
X <- model.matrix(~ SNP + age, data = data.) 

# Run unpenalized ascertainment-corrected marginalized model
initThet <- c(2.5, 0.0, 0.0, log(sqrt(4))) # initial values
optFam.0 <- try(optim(initThet, logLikFam, scoreFam, 
                      y = y, yPatterns = yPatterns, X = X,
                      kinmat = kin.mat, n.mems = n.mems, id = id, 
                      GHs = 11, lambda = 0, truePrev = 0.1, sig.x = 1, MAF = 0.3,
                      marginalized = TRUE, ascert = TRUE, pseudo = FALSE,
                      penalty.glmm = FALSE, penalty.glmmN = FALSE,
                      method = "BFGS", hessian = FALSE,
                      control = list(trace = 10, maxit = 500, 
                                     parscale = c(0.1, 0.1, 0.1, 0.01))))

#Output
if(!inherits(optFam.0, "try-error")){
  if(optFam.0$convergence == 0 & optFam.0$value > 0){
    parsP0 <- c(optFam.0$par[-length(optFam.0$par)], exp(optFam.0$par[length(optFam.0$par)])) 
    cat("\nParameter Estimates (betas, sigmab):", parsP0)
  }
  if(optFam.0$convergence == 0 & optFam.0$value > 0){
    library(JM)
    hes0 <- Hessian.logLikFam(optFam.0$ par, y, yPatterns = yPatterns,
                              X, kinmat = kin.mat, n.mems, id, GHs = 11,
                              lambda = 0, truePrev = 0.1, sig.x = 1, MAF = 0.3,
                              marginalized = TRUE, ascert = TRUE, pseudo = FALSE,
                              penalty.glmm = FALSE, penalty.glmmN = FALSE)
    cat("\nStd Errors:", sqrt(diag(solve(hes0)))) 

    # Using the prevalence information would be helpful - this is what we do next

  }
}



# Run penalized ascertainment-corrected marginalized model
initThet <- c(2.5, 0.0, 0.0, log(sqrt(4))) # initial values
optFam <- try(optim(initThet, logLikFam, scoreFam, 
                        y = y, yPatterns = yPatterns, X = X,
                        kinmat = kin.mat, n.mems = n.mems, id = id, 
                        GHs = 11, lambda = 100, truePrev = 0.1, sig.x = 1, MAF = 0.3,
                        marginalized = TRUE, ascert = TRUE, pseudo = FALSE,
                        penalty.glmm = FALSE, penalty.glmmN = FALSE,
                        method = "BFGS", hessian = FALSE,
                        control = list(trace = 10, maxit = 500, 
                                       parscale = c(0.1, 0.1, 0.1, 0.01))))
  
#Output
if(!inherits(optFam, "try-error")){
    if(optFam$convergence == 0 & optFam$value > 0){
      parsP <- c(optFam$par[-length(optFam$par)], exp(optFam$par[length(optFam$par)])) 
      cat("\nParameter Estimates (betas, sigmab):", parsP)
      
    }
    if(optFam$convergence == 0 & optFam$value > 0){
      library(JM)
      hes <- Hessian.logLikFam(optFam$ par, y, yPatterns = yPatterns,
                                X, kinmat = kin.mat, n.mems, id, GHs = 11,
                                lambda = 100, truePrev = 0.1, sig.x = 1, MAF = 0.3,
                                marginalized = TRUE, ascert = TRUE, pseudo = FALSE,
                                penalty.glmm = FALSE, penalty.glmmN = FALSE)
      cat("\nStd Errors:", sqrt(diag(solve(hes))))
    }
}    





gc()
