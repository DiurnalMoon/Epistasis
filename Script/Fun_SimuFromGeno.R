# Dependent package
library(igraph)
library(MASS)

Simu_pheno <- function(add, dom, n.main, n.size, n.orlp, p.neg,
                    dataDir, scenario,
                    var_cop=NULL, seed=1984,
                    powerrate=1, exp=1, rho=0.7) {
  
  set.seed(seed+1984)
  EffType <- c("A", "D", "AA", "AD", "DD", "Residual")

  ### Simulate scare-free graphs.
  
  # Simulate three graphs with the same size and the same set of nodes.
  # represent one single gene network, 
  # but with three types of epistatic effect.

  negP <- as.list(p.neg)  # Percentage of negative
  names(negP) <- EffType[1:5]
  size <- rep(n.size, 3) # Size of network
  mainsize <- n.orlp # Number of QTLs with both main and add-add effect
  Graphs <- SimuGraph(n=length(size), rate=powerrate,
                      size=size, mainsize=mainsize)
  
  ### Generate Effect Size
  
  # Set correlation between effect size and connectivity
  if (length(rho) == 1) rho <- rep(rho, 3)
  rho_add <- rho[1]
  rho_dom <- rho[2]
  rho_epi <- rho[3]
  
  # Epistatic effect
  if (length(exp) == 1) exp <- rep(exp, 5)
  exprate <- exp[3:5]
  tol <- 0.02
  Epi_Index <- data.frame(Node1=c(),
                          Node2=c(),
                          Group=c(),
                          Coef=c(),
                          degr1=c(),
                          degr2=c(), 
                          Type=c())
  for (i in 1:length(Graphs)) {
    graph <- Graphs[[i]]
    degr1 <- graph$degr[graph$edge[, 1]]
    degr2 <- graph$degr[graph$edge[, 2]]
    weight <- degr1 + degr2
    weight <- weight / sqrt(var(weight) * exprate[i]^2)
    coef <- SimuEff(weight=weight, rate=exprate[i],
                    rho0=rho_epi, bound=tol, negP=negP[[i+2]])
    Epi_Index <- rbind(Epi_Index,
                       data.frame(Node1=graph$edge[, 1],
                                  Node2=graph$edge[, 2],
                                  Group=rep(i, length(coef)),
                                  Coef=coef,
                                  degr1=degr1,
                                  degr2=degr2, 
                                  Type=rep(EffType[i+2], length(coef)) ))
  }
  
  # Main effect
  exprate <- exp[1:2]
  tol <- 0.02
  Mai_Index <- data.frame(Node=c(),
                          Group=c(),
                          Coef=c(),
                          Type=c())
  i <- 1
  # Currently only consider the overlap between main effect QTL
  # and QTL with additive-by-additive epistatic effect
  # Sample n.orlp QTLs from the AA network
  # Randomly sample the rest (n.main - n.orlp) from whole genome
  graph <- Graphs[[i]]
  mnod <- graph$mainnode
  degr <- graph$degr[mnod]
  
  # Additive
  coef <- rexp(n=n.main, rate=exp[1])
  Mai_Index <- rbind(Mai_Index,
                     data.frame(Node=mnod,
                                Group=rep(i, n.orlp),
                                Coef=coef[1:n.orlp],
                                Type=rep('A', n.orlp)))
  Mai_Index <- rbind(Mai_Index,
                     data.frame(Node=rep(0, n.main-n.orlp),
                                Group=rep(0, n.main-n.orlp),
                                Coef=coef[(n.orlp+1):length(coef)],
                                Type=rep('A', n.main-n.orlp)))
  # Dominant
  coef <- rexp(n=n.main, rate=exp[2])
  Mai_Index <- rbind(Mai_Index,
                     data.frame(Node=mnod,
                                Group=rep(i, n.orlp),
                                Coef=coef[1:n.orlp],
                                Type=rep('D', n.orlp)))
  Mai_Index <- rbind(Mai_Index,
                     data.frame(Node=rep(0, n.main-n.orlp),
                                Group=rep(0, n.main-n.orlp),
                                Coef=coef[(n.orlp+1):length(coef)],
                                Type=rep('D', n.main-n.orlp)))
  
  
  ### Generate QTLs and genotypic design matrices
  ### for selected loci/pairs
  
  # Select QTLs from all snps
  # Randomly assign snps to nodes
  pools <- sample(1:ncol(add), size = size[1])
  map <- matrix(nrow = size[1], ncol = 2)
  map[, 1] <- 1:size[1]
  map[, 2] <- sample(pools, size =size[1])
  
  Epi_Index <- cbind(Epi_Index,
                     Loc1=rep(0, nrow(Epi_Index)),
                     Loc2=rep(0, nrow(Epi_Index)))
  Epi_Index$Loc1 <- map[Epi_Index$Node1, 2]
  Epi_Index$Loc2 <- map[Epi_Index$Node2, 2]
  
  Mai_Index <- cbind(Mai_Index,
                     Loc=rep(0, nrow(Mai_Index)))
  index <- which(Mai_Index$Node != 0)
  Mai_Index$Loc[index] <- map[Mai_Index$Node[index], 2]
  LonelyMain <- sample((1:ncol(add))[-pools], size = n.main-n.orlp)
  Mai_Index$Loc[-index] <- LonelyMain
  MainQTL <- unique(Mai_Index$Loc)
  
  # Calculate genotypic design matrices
  A <- add[, MainQTL]
  colnames(A) <- paste(MainQTL, "A", sep = "-")
  D <- dom[, MainQTL]
  colnames(D) <- paste(MainQTL, "D", sep = "-")
  
  Interact <- function(X, G1, G2){
    return(G1[,X[1]] * G2[,X[2]])
  }
  
  index <- as.matrix(Epi_Index[which(Epi_Index$Type == "AA"), c('Loc1', 'Loc2')])
  AA <- apply(index, 1, Interact, G1=add, G2=add)
  colnames(AA) <- sapply(1:nrow(index), function(i) {
    paste(paste(sort(index[i, ]), collapse = "_"), "AA", sep = "-")
  })
  index <- as.matrix(Epi_Index[which(Epi_Index$Type == "AD"), c('Loc1', 'Loc2')])
  AD <- apply(index, 1, Interact, G1=add, G2=dom)
  colnames(AD) <- sapply(1:nrow(index), function(i) {
    paste(paste(sort(index[i, ]), collapse = "_"), "AD", sep = "-")
  })
  index <- as.matrix(Epi_Index[which(Epi_Index$Type == "DD"), c('Loc1', 'Loc2')])
  DD <- apply(index, 1, Interact, G1=dom, G2=dom)
  colnames(DD) <- sapply(1:nrow(index), function(i) {
    paste(paste(sort(index[i, ]), collapse = "_"), "DD", sep = "-")
  })
  
  sub_geno <- list(A=A, D=D, AA=AA, AD=AD, DD=DD)
  
  # Calculate genetic value based on simulated coefficients
  gen_eff <- list() # each is a vector of length # individual
  gen_var <- c()    # variance of each genetic effects among population
  Coef <- list(A=Mai_Index$Coef[which(Mai_Index$Type == 'A')],
               D=Mai_Index$Coef[which(Mai_Index$Type == 'D')])
  for ( i in 1:3 ) {
    Coef[[i+2]] <- Epi_Index$Coef[which(Epi_Index$Type == names(sub_geno)[i+2])]
  }
  names(Coef) <- names(sub_geno)
  for ( i in 1:length(sub_geno) ) {
    gen_eff[[i]] <- sub_geno[[i]] %*% Coef[[i]]
    gen_var[i] <- var(apply(gen_eff[[i]], 1, sum))
  }
  
  # Scale variance components
  if (is.null(var_cop)) {
    var_cop <- c(0.25, 0.15, 0.2, 0.1, 0.1, 0.2)
  }
  if (length(var_cop) < 6) {
    var_cop <- c(var_cop, rep((1-sum(var_cop))/(6-length(var_cop)), 
                              6-length(var_cop)))
    cat("Simulated variance components: \n")
    print(var_cop)
  }
  gen_var <- c(gen_var, 1)
  var_scl <- gen_var
  coef <- Coef
  for ( i in 1:length(var_cop) ) {
    var_scl[i] <- sqrt(gen_var[1] / var_cop[1] * var_cop[i] /gen_var[i])
  }
  error <- rnorm(nrow(add), 0, var_scl[6])
  coef_new <- lapply(1:5, function(i) {
    coef[[i]] <- coef[[i]] * var_scl[i]
  })
  
  for ( i in 1:5 ) {
    gen_eff[[i]] <- sub_geno[[i]] %*% coef_new[[i]]
    gen_var[i] <- var(apply(gen_eff[[i]], 1, sum))
  }
  gen_eff[[6]] <- error
  gen_var[6] <- var(error)
  var_cop <- sapply(gen_var, function(x) x/sum(gen_var))
  names(var_cop) <- EffType
  
  # Generate final phenotypic value
  pheno <- gen_eff[[1]]
  for ( i in 2:length(gen_eff) ) pheno <- pheno + gen_eff[[i]]
  
  ### Save/return simulated data
  EffectIndex <- list(Mai_Index=Mai_Index,
                      Epi_Index=Epi_Index,
                      Graphs = Graphs,
                      Coefficients=coef_new,
                      Subgeno=sub_geno)
  saveRDS(pheno, paste0(dataDir, 'pheno_',
                        scenario, '_', 'S', seed, '.rds'))
  saveRDS(EffectIndex, paste0(dataDir, 'EffectIndex_',
                              scenario, '_', 'S', seed, '.rds'))
  return(pheno)
}


SimuGraph <- function(n, rate=1, size, mainsize=NULL) {
  
  if (!is.null(mainsize) & length(mainsize) < n) {
    mainsize <- c(mainsize, rep(0, n-length(mainsize)))
  }
  Graphs <- list()
  if (length(rate) == 1) rate <- rep(rate, n)
  if (length(size) == 1) size <- rep(size, n)
  for ( i in 1:n ) {
    nets <- barabasi.game(n = size[i], directed = FALSE)
    edge <- get.edgelist(nets)
    degr <- degree(nets)
    if (!is.null(mainsize)) {
      mnot <- sample(x=1:size[i], size=mainsize[i], prob = degr)
      Graphs[[i]] <- list(network=nets, edge=edge, degr=degr, mainnode=mnot)
    } else {
      Graphs[[i]] <- list(network=nets, edge=edge, degr=degr)
    }
  }
  return(Graphs)
}


SimuEff <- function(weight, rate, rho0, bound, negP=0.5) {
  
  rho <- 0
  n <- length(weight)
  V <- matrix(c(1, rho0, rho0, 1), nrow = 2)
  while (rho < rho0-bound | rho > rho0+bound) {
    z <- mvrnorm(n, mu = rep(0, 2), Sigma=V, empirical=T)
    u <- pnorm(z)
    orind <- order(u[, 2])
    x2 <- weight
    so <- sort(weight)
    for ( i in 1:n ) {
      x2[orind[i]] <- so[i]
    }
    x1 <- qexp(u[, 1], rate = 1)
    rho <- cor(x1, x2)
  }
  x3 <- x1
  for (i in 1:n) {
    index <- which(weight == x2[i])[1]
    x3[index] <- x1[i]
    weight[index] <- -1
  }
  negindex <- sample(1:length(x3), size=negP*length(x3))
  x3[negindex] <- (-1) * x3[negindex]
  return(x3)
}