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
