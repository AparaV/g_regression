##
## Check identifiability condition for total effect
##

directed_path <- function(G, source, dest, forbidden=c()) {
    # cat("From", source, "to", dest, "\n")
    
    vertex_names <- colnames(G)
    source_id <- which(vertex_names == source)[1]
    dest_id <- which(vertex_names == dest)[1]
    
    if (is.na(source_id) || is.na(dest_id)) {
        stop("Source and destination not found in G. Provide vertex names.")
    }
    
    
    nbrs <- which(G[source_id, ] == 1, arr.ind=TRUE)
    # cat("Neighbors", vertex_names[nbrs], "\n")
    if (length(nbrs) == 0) {
        return(NULL)
    }
    
    for (i in 1:length(nbrs)) {
        X <- nbrs[i]
        
        # Check if X is in forbidden list
        if (X %in% forbidden) {
            next
        }
        
        # Check if edge is undirected
        if (G[X, source_id] == 1) {
            next
        }
        
        # Check if we reached destination
        # cat("X", X, "dest", dest, "\n")
        if (X == dest_id) {
            return(c(source, dest))
        }
        
        # Recursively find a causal path from X to destination
        # Update forbidden nodes with the source (avoid cycles)
        path <- directed_path(G, vertex_names[X], dest, forbidden)
        if (!is.null(path)) {
            path <- c(source, path)
            return(path)
        }
    }
    
    return(NULL)
}

ancestors <- function(X, G) {
    if (length(X) != 1) {
        stop("X needs to be a  single node")
    }
    
    an <- c(X)
    
    vertices <- colnames(G)
    n <- ncol(G)
    for (i in 1:n) {
        if (X == i) {
            next
        }
        # print(i)
        path <- directed_path(G, vertices[i], X)
        if (!is.null(path)) {
            an <- c(an, vertices[i])
        }
    }
    
    return(an)
}


possibly_causal_path <- function(G, source, dest, forbidden=c(), undirected=TRUE, seen_nodes=c()) {
    
    nbrs <- which(G[source, ] == 1, arr.ind=TRUE)
    if (length(nbrs) == 0) {
        return(NULL)
    }
    
    for (i in 1:length(nbrs)) {
        X <- nbrs[i]
        
        # Check if X is in forbidden list
        if (X %in% forbidden) {
            next
        }
        
        # Check if the first edge is undirected
        if (undirected && G[X, source] == 0) {
            next
        }
        
        # Check if we reached destination
        if (X == dest) {
            return(c(source, dest))
        }
        
        # Check if there is a directed edge to previously visited nodes
        possibly_causal <- TRUE
        for (V in seen_nodes) {
            if (G[X, V] == 1 && G[V, X] == 0) {
                possibly_causal <- FALSE
            }
        }
        
        # Recursively find a possibly causal path from X to destination
        # Update forbidden nodes with the source (avoid undirected cycles)
        # Update seen nodes with source to check for directed edges to source
        if (possibly_causal) {
            path <- possibly_causal_path(G, X, dest, forbidden=c(forbidden, source),
                                         undirected=FALSE, seen_nodes=c(seen_nodes, source))
            if (!is.null(path)) {
                path <- c(source, path)
                return(path)
            }
        }
    }
    
    return(NULL)
}


## The total causal effect A on Y is identified given an MPDAG G if and only if
## there is no proper, possibly causal path from A to Y in G that starts with an undirected edge.
is_identifiable <- function(G, A, Y) {
    
    if (length(Y) != 1) {
        stop("Outcome variable should be univariate")
    }
    
    n <- length(A)
    for (i in 1:n) {
        path <- possibly_causal_path(G, A[i], Y, forbidden=A[-i])
        if (!is.null(path)) {
            return(FALSE)
        }
    }
    
    return(TRUE)
    
}
