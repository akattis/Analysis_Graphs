library(igraph)
library(scatterplot3d)
library(ggplot2)
t0 <- proc.time()

## create the Chung-Lu graphs

# define the parameters
m <- 1000
n <- 100
p <- runif(m, 0.5, 1)
r <- runif(m, 0, 0.5)

# create the first graph
w = matrix(1:n,1,n)
w = n*p[[1]]*((w/n)^r[[1]])
B = (t(w)%*%w)/sum(w)
B = B*(B < 1) + (B > 1)
Q = matrix(runif(n*n, 0, 1),n,n)
Q = matrix(1,n,n)*(B > Q)
Q = upper.tri(Q)*Q
Q = Q + t(Q)


G <- list(graph.adjacency(Q, "undirected"))

if (m > 1) {
  for (i in 2:m) {
    w = matrix(1:n,1,n)
    w = n*p[[i]]*((w/n)^r[[i]])
    B = (t(w)%*%w)/sum(w)
    B = B*(B < 1) + (B > 1)
    Q = matrix(runif(n*n, 0, 1),n,n)
    Q = matrix(1,n,n)*(B > Q)
    Q = upper.tri(Q)*Q
    Q = Q + t(Q)
    G <- c(G, list(graph.adjacency(Q, "undirected")))
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
print(t1 - t0)

# Plot everything

# Figure 5
t1 <- qplot(1:10, log(V$values[1:10], data = data.frame(1:10, V$values[1:10]))
t2 <- qplot(p, r, data = data.frame(p,r,V$vectors[,2]), color = V$vectors[,2]) + scale_color_gradientn(colours=rainbow(10))
t3 <- qplot(p, r, data = data.frame(p,r,V$vectors[,3]), color = V$vectors[,3]) + scale_color_gradientn(colours=rainbow(10))

#Only problematic one
t4 <- qplot(p, r, data = data.frame(p,r,V$vectors[,4]), color = V$vectors[,4]) + scale_color_gradientn(colours=rainbow(10))


t5 <- qplot(p, r, data = data.frame(p,r,V$vectors[,5]), color = V$vectors[,5]) + scale_color_gradientn(colours=rainbow(10), limits = c(-0.15, 0.15))

#Figure 6
t6 <- qplot(V$vectors[,2], V$vectors[,3], data = data.frame(V$vectors[,2], V$vectors[,3], p), color = p) + scale_color_gradientn(colours=rainbow(10))
t7 <- qplot(V$vectors[,2], V$vectors[,3], data = data.frame(V$vectors[,2], V$vectors[,3], r), color = r) + scale_color_gradientn(colours=rainbow(10))

# Figure 7
t8 <- scatterplot3d(V$vectors[,2], V$vectors[,3], V$vectors[,4])
