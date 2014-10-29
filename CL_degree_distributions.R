library(igraph)

## create the Chung-Lu graphs

# define the parameters
m <- 1000
n <- 100
pA <- seq(0.7, 1, 0.1)
rA <- seq(0, 0.09, 0.03)

dstrb <- array(0, c(length(pA), length(rA), m, n))

for (aa in 1:length(pA)) {
for (bb in 1:length(rA)) {
p = pA[[aa]]
r = rA[[bb]]

# create the first graph
w = matrix(1:n,1,n)
w = n*p*((w/n)^r)
B = (t(w)%*%w)/sum(w)
B = B*(B < 1) + (B > 1)
Q = matrix(runif(n*n, 0, 1),n,n)
Q = matrix(1,n,n)*(B > Q)
Q = upper.tri(Q)*Q
Q = Q + t(Q)


G <- list(graph.adjacency(Q, "undirected"))
dstrb[aa,bb,1,1:length(degree.distribution(G[[1]]))] = degree.distribution(G[[1]])

if (m > 1) {
  for (i in 2:m) {
    w = matrix(1:n,1,n)
    w = n*p*((w/n)^r)
    B = (t(w)%*%w)/sum(w)
    B = B*(B < 1) + (B > 1)
    Q = matrix(runif(n*n, 0, 1),n,n)
    Q = matrix(1,n,n)*(B > Q)
    Q = upper.tri(Q)*Q
    Q = Q + t(Q)
    G <- c(G, list(graph.adjacency(Q, "undirected")))

dstrb[aa,bb,i,1:length(degree.distribution(G[[i]]))] = degree.distribution(G[[i]])
  }
}


}
}
