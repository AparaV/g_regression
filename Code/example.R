##
## Simple example to work through
##

library(pcalg)
library(igraph)
library(RBGL)

source("bucket_decomposition.R")
source("mle.R")


set.seed(3)

num_nodes <- 10
poss_degrees <- c(2, 3, 4, 5)
avg_degree <- 3

prob_edge <- avg_degree / (num_nodes - 1)

# G <- erdos.renyi.game(num_nodes, prob_edge, type="gnp")
# dag <- r.gauss.pardag(num_nodes, prob_edge, lbe=1, ube=1)
# dag <- randomDAG(num_nodes, prob_edge, lB=1, uB=1)

# Function to sample edge weights
sample_edge_weights <- function(m) {
    abs_weight <- runif(m, min=0.1, max=2)
    sign_weight <- -1 + 2*rbinom(m, 1, 0.5)
    weight <- abs_weight * sign_weight
    return(weight)
}

dag <- randDAG(num_nodes, avg_degree, weighted=TRUE, wFUN=list(sample_edge_weights))

# Topologically sort the DAG so it is easier to generate data
sorted_dag <- tsort(dag)
dag_matrix <- as(dag, "matrix")
dag_matrix <- dag_matrix[sorted_dag, sorted_dag]

# Error distribution 1 in Guo + Perkovic (2022)
error_dist_1 <- function(n) {
    v <- runif(n, min=0.5, max=6)
    epsilon <- rnorm(n, mean=0, sd=v)
    return(epsilon)
}


# Method to sample data from DAG
sample_data <- function(n, dag, err) {
    # n: Number of data points
    # dag: Topologically sorted dag
    # err: Function that generates error
    p <- nrow(dag)
    X <- matrix(0, nrow=n, ncol=p)
    X[, p] <- err(n)
    for (i in (p-1):1) {
        X[, i] <- dag[i, ] %*% t(X) + err(n)
    }
    return(X)
}

X <- sample_data(100, dag_matrix, err=error_dist_1)

# Need to work with the CPDAG for bucket decomposition
cpdag <- dag2cpdag(dag)
cpdag_matrix <- as(cpdag, "matrix")

buckets <- ordered_bucket_decomposition(cpdag_matrix)

lambda_mle <- find_lambda_mle(X, cpdag_matrix, buckets)


