source("/mnt/gt4sp_1/ysi2/SimuEpi/Script/Fun_network.R")

list.of.packages <- c("sommer", "glmnet")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

MarginalAA <- function(pheno, add, dom, dataDir, Scenario, thin, seed, boot,
                       penals=c(1, 0.8, 0.5), pcut=0.995) {

  set.seed(seed + 1984)

  # Step I

  stepfile <- paste0(dataDir, 'StepI_', Scenario ,'.rds')
  if (file.exists(stepfile)) {
    Step1 <- readRDS(stepfile)
  } else {
    Step1 <- Lasso_A(X=add, Y=pheno)
    saveRDS(Step1, stepfile)
  }

  Summer <- Step1$Summary[, 2:4]
  index <- match(Step1$Lambda1se, Summer$Lambda)
  Sset <- Step1$Coefficients[[index[2]]]
  if (length(which(is.na(Sset[, 2]) == TRUE)) != 0) {
    Sset <- Sset[-which(is.na(Sset[, 2]) == TRUE), ]
  }
  loci <- as.character(Sset$Loci)
  main <- add[, loci]

  GRM <- function(M1, M2=NULL){

    if (is.null(M2)) {
      K <- tcrossprod(M1)
    } else {
      K <- tcrossprod(M1, M2)
    }
    K <- K/mean(diag(K))
    return(K)

  }

  bootstrapAnova <- function(mA, m0, GRMi, B=1000){
    y <- m0$model$pheno
    sigmaN <- var(m0$residuals)
    sigmaA <- var(mA$cond.residuals)
    V <- mA$var.comp[1] * GRMi +
         mA$var.comp[2] * diag(nrow=nrow(GRMi))
    Vinv <- chol2inv(chol(V))
    n <- nrow(V)

    oneBootstrap <- function(m0, mA){
      d <- drop(simulate(m0))[, 1]

      l0 <- (-n/2) * log(sigmaN) - sum((d-m0$fitted.values)^2)/sigmaN/2
      l1 <- (-n/2) * log(sigmaA) - t(d-mA$fitted.y) %*%
        Vinv %*% (d-mA$fitted.y) /2
      return(l0-l1)
    }

    nulldist <- replicate(B, oneBootstrap(m0, mA))

    l0 <- (-n/2) * log(sigmaN) - sum((y-m0$fitted.values)^2)/sigmaN/2
    l1 <- (-n/2) * log(sigmaA) - t(y-mA$fitted.y) %*%
      Vinv %*% (y-mA$fitted.y) /2
    tao <- as.numeric(l0-l1)


    ret <- length(which(nulldist > tao)) / length(nulldist)
    return(ret)
  }

  set <- seq(1, ncol(add), by=thin)

  stepfile <- paste0(dataDir,'MargAA_', Scenario, '.rds')
  if (file.exists(stepfile)) {
  result <- readRDS(stepfile)
  } else {
  Stats <- sapply(set, function(i) {
    if (length(grep(paste0('^', colnames(add)[i], '$'), colnames(main))) == 0 ) {
      X <- cbind(main, add[,i], dom[, i])
    } else {
      X <- main
    }
    m0 <- lm(pheno~X)
    mat <- apply(add[, -i], 2, function(x) {
      x * add[, i]
    })
    GRMi <- GRM(mat)
    I <- diag(nrow = nrow(add))
    RE <- list(MA=list(Z=I, K=GRMi))
    mA <- mmer(Y=pheno, X=X, Z=RE,
               method="EMMA", silent = TRUE)
    p <- bootstrapAnova(mA, m0, GRMi, B=boot)
    return(p)
  })
  result <- rbind(set, unlist(Stats))
  saveRDS(result, paste0(dataDir,'MargAA_', Scenario, '.rds'))

  }

  # Steo II
  print("Sart step II")

  Interact <- function(X, G1, G2){
    return(G1[,X[1]] * G2[,X[2]])
  }

  set.seed(seed + 1984)

  Summer <- Step1$Summary[, 2:4]
  index <- match(Step1$Lambda1se, Summer$Lambda)
  for (a in 1:length(index)) {
    if (is.na(index[a])) {
      index[a] <- which.min(abs(Summer$Lambda - Step1$Lambda1se[a]))
  }}
  Sset <- Step1$Coefficients[[index[2]]]
  if (length(which(is.na(Sset[, 2]) == TRUE)) != 0) {
    Sset <- Sset[-which(is.na(Sset[, 2]) == TRUE), ]
  }

  loci <- as.character(Sset$Loci)
  force_in <- paste0(rep(loci, each=2), c('-A', '-D'))

  locs <- as.integer(result[1, ])

  epicandy <- locs[which(result[2, ] > pcut)]

  Lset <- Step1$Coefficients[[index[3]]]
    if (length(which(is.na(Lset[, 2]) == TRUE)) != 0) {
    Lset <- Lset[-which(is.na(Lset[, 2]) == TRUE), ]
  }

  epicandy <- union(match(as.character(Lset$Loci), colnames(add)), epicandy)
  saa_index <- combn(as.vector(epicandy), 2)
  Effcandy <- paste(saa_index[1, ], saa_index[2, ], sep = '_')

  sA <- add[, epicandy]
  sD <- dom
  sAA <- apply(saa_index, 2, Interact, G1=add, G2=add)
  colnames(sA) <- paste(colnames(sA), 'A', sep='-')
  colnames(sD) <- paste(colnames(sD), 'D', sep='-')
  colnames(sAA) <- paste0(Effcandy, '-AA')

  Table <- matrix(0, nrow = length(penals), ncol = 5)
  colnames(Table) <- c("Penalty.Ratio", "Sum.D.Hat", "Sum.AA.Hat",
                       "V.D.Hat", "V.AA.Hat")

  for ( i in 1:length(penals) ) {

    print(penals[i])

    FinalGeno <- cbind(sA, sD, sAA)

    pety <- penals[i]
    Table[i, "Penalty.Ratio"] <- pety
    penalty <- rep(1, ncol(FinalGeno))
    penalty[colnames(FinalGeno) %in% colnames(dom) == TRUE] <- pety
    penalty[colnames(FinalGeno) %in% force_in == TRUE] <- 0

    SubScenario <- paste0('R', pety, '-', Scenario)
    stepfile <- paste0(dataDir, 'StepII_', SubScenario, '.rds')

    if (file.exists(stepfile)) {
      Step2 <- readRDS(stepfile)
    } else {
      Step2 <- LassoII(X=FinalGeno, Y=pheno, penalty=penalty)
      saveRDS(Step2, stepfile)
    }

    lambdamin <- Step2$Lambda1se[3]
    pin <- which(Step2$Summary$Lambda == lambdamin)
    if (lambdamin < max(Step2$Summary$Lambda)) {
      pin <- which.min(Step2$Summary$Lambda)
    }
    FinalGeno <- FinalGeno[, Step2$Coefficients[[pin]]$Eff]
    coef <- Step2$Coefficients[[pin]]$Beta
    names(coef) <- colnames(FinalGeno)

    Table[i, 'Sum.D.Hat'] <- sum(coef[grep("-D$", names(coef))])
    Table[i, 'Sum.AA.Hat'] <- sum(coef[grep("-AA$", names(coef))])

    dind <- grep("-D$", names(coef))
    aind <- grep("-AA$", names(coef))
    gDhat <- as.matrix(FinalGeno[, dind]) %*% as.vector(coef[dind])
    gAAhat <- as.matrix(FinalGeno[, aind]) %*% as.vector(coef[aind])
    Table[i, "V.D.Hat"] <- var(gDhat, na.rm=TRUE)
    Table[i, "V.AA.Hat"] <- var(gAAhat, na.rm=TRUE)
  }

  saveRDS(Table,
      file=paste0(dataDir, 'SumEff_', Scenario, '.rds'))

  return(Table)
}














