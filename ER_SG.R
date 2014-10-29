# Generates Erdos Renyi graphs
# Uses subgraph method to calculate DMAPs
#
library(igraph)

t0 <- proc.time()
# create the erdos-renyi graphs
# define the parameters
m <- 1000
n <- 100
p <- runif(m, 0, 1)

G <- list(erdos.renyi.game(n, p[[1]], "gnp"))
if (m > 1) {
  for (i in 2:m) {
    G <- c(G, list(erdos.renyi.game(n, p[[i]], "gnp")))
  }
}

# find the motifs of each graph
rho <- matrix(,9, m)
for (i in 1:m) {

  if (n > 1) {
    rho[1,i] <- sum(degree(G[[i]])*(degree(G[[i]]) == 1))/(2*choose(n,2))#(ecount(G[[i]])/choose(n,2))
  } else
    {
      rho[1,i] = 0
    }

  if (n > 2) {
    d <- graph.motifs(G[[i]], 3)
    d <- d[!is.na(d)]
    rho[2:3,i] <- (d/choose(n,3)) }
  else {
    rho[2:3,i] <- rep(0, 2)
}

if (n > 3) {
d <- graph.motifs(G[[i]], 4)
d <- d[!is.na(d)]
rho[4:9,i] <- (d/choose(n, 4))
} else {
rho[4:9,i] <- rep(0,6)
}
}

# Perform the DMAP calculation
dist = matrix(0, m, m)
epsilon = 0.5

for (i in 1:m) {
  for (j in 1:m) {
    dist[i,j] <- sum((rho[,i] - rho[,j])^2)
  }
}

W = exp(-dist/(epsilon^2))
A = matrix(0,m,m)

for (i in 1:m) {
  for (j in 1:m) {
    A[i,] = W[i,]/sum(W[i,])
  }
}

V <- eigen(A, TRUE)
t1 <- proc.time()
print(t1 - t0)



plot(log(Re(V$values[1:10])), xlab = "Index", ylab = "log(Eigenvalues)" )
title("Log(Eigenvalues) vs. Index; Motif Metric")

plot(V$vectors[,2], p, xlab = "v2", ylab = "p" )
title("Eigenvector 2; Motif Metric")

plot(p, V$vectors[,3], ylab = "v3", xlab = "p" )
title("Eigenvector 3; Motif Metric")

plot(V$vectors[,2], V$vectors[,3], xlab = "v2", ylab = "v3" )
plot(V$vectors[,2], V$vectors[,4], xlab = "v2", ylab = "v4" )
plot(V$vectors[,2], V$vectors[,5], xlab = "v2", ylab = "v5" )
