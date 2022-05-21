##
## Sample DAG and data
##

library(pcalg)
library(igraph)
library(RBGL)

source("cpdag.R")


# Function to sample edge weights
sample_edge_weights <- function(m) {
    abs_weight <- runif(m, min=0.1, max=2)
    sign_weight <- -1 + 2*rbinom(m, 1, 0.5)
    weight <- abs_weight * sign_weight
    return(weight)
}


# Error distribution 1 in Guo + Perkovic (2022)
error_dist_1 <- function(n) {
    v <- runif(n, min=0.5, max=6)
    epsilon <- rnorm(n, mean=0, sd=v)
    return(epsilon)
}

error_dist_2 <- function(n) {
    v <- runif(n, min=0.5, max=1.5)
    epsilon <- sqrt(v) * rt(n, 5)
    return(epsilon)
}

error_dist_3 <- function(n) {
    s <- runif(n, min=0.4, max=0.7)
    epsilon <- rlogis(n, location=0, scale=s)
    return(epsilon)
}

error_dist_4 <- function(n) {
    a <- runif(n, min=1.2, max=2.1)
    epsilon <- runif(n, min=-a, max=a)
    return(epsilon)
}


compute_error <- function(est, true) {
    return(norm(true-est, type="2"))
}

sample_dag <- function(n, k, homebrew=TRUE) {
    dag <- pcalg::randDAG(n, k, weighted=TRUE, wFUN=list(sample_edge_weights))
    dag_matrix <- as(dag, "matrix")
    
    # Need to work with the CPDAG for bucket decomposition
    if (homebrew) {
        cpdag <- NULL
        cpdag_matrix <- get_cpdag(dag_matrix, verbose=T)
    }
    else {
        cpdag <- pcalg::dag2cpdag(dag)
        cpdag_matrix <- as(cpdag, "matrix")
    }
    
    # cpdag_matrix <- cpdag_matrix[sorted_dag, sorted_dag]
    
    # Get gamma matrix by sorting rows and columns in DAG
    # cols <- colnames(dag_matrix)
    # sorted_names <- as.character(sort(as.integer(cols)))
    # gamma <- dag_matrix[sorted_names, sorted_names]
    
    result <- list(
        dag=dag,
        cpdag=cpdag,
        dag_amat=dag_matrix,
        cpdag_amat=cpdag_matrix
    )
    return(result)
}

sample_treatment_outcome <- function(dag_amat, cpdag_amat, A_len, cpdag_est=NULL) {
    
    # Find nodes in DAG that have non-empty descendants
    dag_amat_og <- dag_amat
    dag_amat[which(dag_amat != 0)] <- 1
    A_candidates <- which(rowSums(dag_amat) != 0)
    
    # Search through all subsets of size A_len to find a pair (A, Y)
    # such that the effect of A on Y is identified
    possible_treatments <- combn(A_candidates, A_len, simplify=FALSE)
    indices_treatments <- seq(1:length(possible_treatments))
    
    found_A <- FALSE
    while(!found_A) {
        
        if (length(indices_treatments) == 0) {
            break
        }
        
        # We sample a treatment set and make sure we never pick it again
        idx <- sample(length(indices_treatments), 1)
        A_idx <- indices_treatments[idx]
        A <- possible_treatments[[A_idx]]
        indices_treatments <- indices_treatments[-idx]
        
        # Find all possible descendants of treatment set
        Y_candidates <- c()
        for (a in A) {
            a_descendants <- descendants(dag_amat_og, a)
            Y_candidates <- c(Y_candidates, a_descendants)
        }
        Y_candidates <- unique(Y_candidates)
        if (length(Y_candidates) == 0) {
            next
        }
        
        # For non-empty descendant sets, find an outcome so that the effect
        # is identified
        identified <- FALSE
        for (i in 1:length(Y_candidates)) {
            Y <- Y_candidates[i]
            if (Y %in% A) {
                next
            }
            # cat("A", A, ", Y", Y, "\n")
            identified <- is_identifiable(cpdag_amat, A, Y)
            if(!is.null(cpdag_est)) {
                identified <- identified & is_identifiable(cpdag_est, A, Y)
            }
            if (identified) {
                found_A <- TRUE
                break
            }
        }
    }
    
    if (found_A) {
        return(list(A=A, Y=Y))
    }
    else {
        return(NULL)
    }
}




sort_dag <- function(dag) {
    # Topologically sort the DAG so it is easier to generate data
    sorted_nodes <- rev(RBGL::tsort(dag))
    dag_mat <- as(dag, "matrix")
    dag_mat <- dag_mat[sorted_nodes, sorted_nodes]
    return(dag_mat)
}

dag_to_SEM <- function(dag_amat) {
    return(t(dag_amat))
}


# Method to sample data from DAG
sample_data <- function(n, dag, err) {
    # n: Number of data points
    # dag: Topologically sorted dag
    # err: Function that generates error
    
    sorted_dag_amat <- sort_dag(dag)
    gamma <- dag_to_SEM(sorted_dag_amat)
    
    p <- nrow(gamma)
    X <- matrix(0, nrow=n, ncol=p)
    X[, p] <- err(n)
    for (i in (p-1):1) {
        X[, i] <- gamma[i, ] %*% t(X) + err(n)
    }
    
    # The columns in dag are topologically sorted
    # So unsort the columns in the data matrix
    X <- data.frame(X)
    cols <- colnames(gamma)
    colnames(X) <- cols
    sorted_names <- as.character(sort(as.integer(cols)))
    X <- X[, sorted_names]
    X <- as(X, "matrix")
    
    return(X)
}