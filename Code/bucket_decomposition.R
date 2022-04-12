##
## Function to get ordered bucket decomposition
## Algorithm 2 of Guo and Perkovic 2022
##

library(igraph)
library(RBGL)

bucket_decomposition <- function(G) {
    # G is an MPDAG
    
    # Find the undirected subgraph
    G_undirected <- (G + t(G)) / 2
    G_undirected[G_undirected <= 0.5] <- 0
    G_U <- graph.adjacency(G_undirected, mode="undirected")
    G_U <- as_graphnel(G_U)
    
    # Get all connected components
    connected_components <- connectedComp(G_U)
    
    return(connected_components)
}

ordered_bucket_decomposition <- function(G) {
    # G is an MPDAG
    
    # Get bucket decomposition
    connected_components <- bucket_decomposition(G)
    
    # Construct the partial causal ordering decomposition
    buckets <- list()
    n <- length(connected_components)
    remaining_vertices <- V
    while (n > 0) {
        C <- connected_components[[1]]
        C_prime <- setdiff(remaining_vertices, C)
        
        # Check whether all edges between C and C_prime are into C
        into <- TRUE
        for (x in C_prime) {
            for (y in C) {
                if (G[y, x] != 0) {
                    into <- FALSE
                    break
                }
            }
            if (!into) {
                break
            }
        }
        
        # If they are, add C to beginning of bucket list
        if (into) {
            buckets <- c(list(C), buckets)
            connected_components <- connected_components[-1]
            n <- n - 1
            remaining_vertices <- C_prime
        }
        # Otherwise, try another connected component in the next iteration
        else {
            connected_components <- c(connected_components[-1], list(C))
        }
        
    }
    
    return(buckets)
}

V <- colnames(G)

G <- cpdag_matrix
G_undirected <- (G + t(G)) / 2
G_undirected[G_undirected <= 0.5] <- 0


G_U <- graph.adjacency(G_undirected, mode="undirected")
G_U <- as_graphnel(G_U)
connected_components <- connectedComp(G_U)


