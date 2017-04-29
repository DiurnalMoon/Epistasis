Lasso_A <- function(Y, X) {
  
  Regress <- cv.glmnet(X, Y)
  min.cvm <- min(Regress$cvm)
  upper.cvm <- min.cvm + Regress$cvsd[match(min.cvm, Regress$cvm)]
  candi <- intersect(which(Regress$cvm < upper.cvm),
                     which(Regress$nzero > 1))
  lambda <- Regress$lambda[candi]
  m <- length(lambda)
  Summer <- data.frame(No.SNP = integer(m),
                       Adj.R.squared = numeric(m),
                       Lambda = lambda,
                       CV.Err = Regress$cvm[candi])
  Trace <- vector("list", length = m)
  for (i in 1:length(lambda)) {
    lam <- lambda[i]
    coef <- as.matrix(coef(Regress, s=lam))[, 1]
    coef <- coef[which(coef != 0)]
    coef <- coef[-1]
    
    reg <- lm(Y ~ X[, names(coef)])
    loci <- unlist(strsplit(names(coef), split="-A"))
    Summer$No.SNP[i] <- length(coef)
    Summer$Adj.R.squared[i] <- summary(reg)$adj.r.squared
    betaA <- reg$coefficients[-1]
    Trace[[i]] <- data.frame(Loci=loci, BetaA=betaA)
    rownames(Trace[[i]]) <- loci
  }
  
  lambdamax <- min(Regress$lambda.1se, max(lambda))
  if (Regress$lambda.min %in% lambda == FALSE) {
    lambdamin <- lambda[round(length(lambda)/2)]
  } else {
    lambdamin <- Regress$lambda.min
  }
  return(list(Summary=Summer, Coefficients=Trace,
              Lambda1se=c(lambdamax, lambdamin, min(lambda))))
}

Lasso <- function(Y, X) {
  
  Regress <- cv.glmnet(X, Y)
  min.cvm <- min(Regress$cvm)
  upper.cvm <- min.cvm + Regress$cvsd[match(min.cvm, Regress$cvm)]
  candi <- intersect(which(Regress$cvm < upper.cvm),
                     which(Regress$nzero > 1))
  lambda <- Regress$lambda[candi]
  m <- length(lambda)
  Summer <- data.frame(No.SNP = integer(m),
                       Adj.R.squared = numeric(m),
                       Lambda = lambda,
                       CV.Err = Regress$cvm[candi])
  Trace <- vector("list", length = m)
  for (i in 1:length(lambda)) {
    lam <- lambda[i]
    coef <- as.matrix(coef(Regress, s=lam))[, 1]
    coef <- coef[which(coef != 0)]
    coef <- coef[-1]
    Kept_SNP <- sapply(strsplit(names(coef), split = "-"), function(x) x[1])
    Kept_SNP <- unique(Kept_SNP)
    Kept_Geno <- X[, paste0(rep(Kept_SNP, each=2),
                            rep(c("-A", "-D"), length(Kept_SNP))) ]
    reg <- lm(Y ~ Kept_Geno)
    Summer$No.SNP[i] <- length(coef)
    Summer$Adj.R.squared[i] <- summary(reg)$adj.r.squared
    betaA <- reg$coefficients[grep('-A', names(reg$coefficients))]
    betaD <- reg$coefficients[grep('-D', names(reg$coefficients))]
    Trace[[i]] <- data.frame(Loci=Kept_SNP, BetaA=betaA, BetaD=betaD)
    rownames(Trace[[i]]) <- Trace[[i]]$Loci
  }
  
  lambdamax <- min(Regress$lambda.1se, max(lambda))
  if (Regress$lambda.min %in% lambda == FALSE) {
    lambdamin <- lambda[round(length(lambda)/2)]
  } else {
    lambdamin <- Regress$lambda.min
  }
  return(list(Summary=Summer, Coefficients=Trace,
              Lambda1se=c(lambdamax, lambdamin, min(lambda))))
}

LassoII <- function(X, Y, penalty) {
  
  Regress <- cv.glmnet(x=X, y=Y, penalty.factor = penalty)
  print("Reg_LassoII")
  
  nm <- length(which(penalty == 0))
  force_in <- 1:nm
  min.cvm <- min(Regress$cvm)
  upper.cvm <- min.cvm + Regress$cvsd[match(min.cvm, Regress$cvm)]
  candi <- intersect(which(Regress$cvm < upper.cvm),
                     which(Regress$nzero > 1))
  lambda <- Regress$lambda[candi]
  m <- length(lambda)
  Summer <- data.frame(No.Pair = integer(m),
                       Adj.R.squared = numeric(m),
                       Lambda = lambda,
                       CV.Err = Regress$cvm[candi])
  Trace <- vector("list", length = m)
  for (i in 1:length(lambda)) {
    lam <- lambda[i]
    coef <- as.matrix(coef(Regress, s=lam))[, 1]
    coef <- coef[which(coef != 0)]
    coef <- coef[-1]
    Kept_All <- names(coef)
    Kept_Geno <- X[, Kept_All]
    reg <- lm(Y ~ Kept_Geno)
    
    beta <- reg$coefficients[-1]
    names(beta) <- Kept_All
    beta <- beta[which(is.na(beta) == FALSE)]
    Summer$No.Pair[i] <- length(beta) - nm
    Summer$Adj.R.squared[i] <- summary(reg)$adj.r.squared
    Trace[[i]] <- data.frame(Eff=names(beta), Beta=beta)
    rownames(Trace[[i]]) <- Trace[[i]]$Eff
    Trace[[i]]$Eff <- as.character(Trace[[i]]$Eff)
  }
  
  lambdamax <- min(Regress$lambda.1se, max(lambda))
  if (Regress$lambda.min %in% lambda == FALSE) {
    lambdamin <- lambda[round(length(lambda)/2)]
  } else {
    lambdamin <- Regress$lambda.min
  }
  Step2 <- list(Summary=Summer, Coefficients=Trace,
                Lambda1se=c(lambdamax, lambdamin, min(lambda)))
  return(Step2)
}



