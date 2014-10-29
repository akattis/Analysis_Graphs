# Uses Erdos Renyi models
# Using the spectral approach, calculates DMAPs

library(igraph)

t0 <- proc.time()
# create the erdos-renyi graphs
# define the parameters
m <- 1000
n <- 10
p <- runif(m, 0, 1)

G <- list(erdos.renyi.game(n, p[[1]], "gnp"))
if (m > 1) {
  for (i in 2:m) {
    G <- c(G, list(erdos.renyi.game(n, p[[i]], "gnp")))
  }
}

# find the spectral metric
lam_lo = 0.00001
lam_hi = 0.001
seqlength = 100
lambda = seq(lam_lo, lam_hi, length.out = seqlength)

pee = matrix(1/n, n,1)
q = pee

S <- matrix(,length(lambda), m)
for (i in 1:m) {
DC = eigen(get.adjacency(G[[i]]))
D = diag(DC$values)
P = DC$vectors

  for (j in 1:length(lambda)) {
    S[j,i] = t(q)%*%P%*%exp(lambda[[j]]*D)%*%t(P)%*%pee
  }
}

# Perform the DMAP calculation
dist = matrix(0, m, m)
epsilon = 5

for (i in 1:m) {
  for (j in 1:m) {
    dist[i,j] <- sum((S[,i] - S[,j])^2)
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
