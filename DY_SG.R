#!/usr/bin/env Rscript
library(igraph)

# Simulates dynamic process
# Subgraph approach

# Recreate dynamic process
t0 <- proc.time()
# set parameters

n = 100
t = 1000
r = 0.1
p = 0.3

Gdummy = erdos.renyi.game(n, p, directed=FALSE)
G = list(Gdummy)
distribution = matrix(0, n, t)

for (i in 1:t) {
  a = degree.distribution(Gdummy)
  distribution[1:length(a),i] = a
  pair = sample(1:n, 2, replace=F)
  if (!are.connected(Gdummy, pair[[1]], pair[[2]])) {
    Gdummy = add.edges(Gdummy,pair)
  }

  if (runif(1,0,1) < r) {
    Gdummy = delete.edges(Gdummy,E(Gdummy)[sample(1:length(E(Gdummy)), 1)])
  }
G = c(G, list(Gdummy))
}

m = length(G)

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
epsilon = 10

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
save(list = ls(all = TRUE), file = "dynamic_DMAP.RData")
