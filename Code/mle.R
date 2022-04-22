##
## Transform from Gamma to Lambda space
##

require(pracma)
source("identifiablity.R")
source("bucket_decomposition.R")

sample_covariance <- function(X) {
    n <- nrow(X)
    sigma <- t(X) %*% X / n
    return(sigma)
}


find_lambda_mle <- function(X, cpdag, buckets) {
    p <- ncol(X)
    sigma_n <- sample_covariance(X)
    K <- length(buckets)
    lambda <- matrix(0, nrow=p, ncol=p)
    
    dimnames(sigma_n) <- dimnames(lambda) <- dimnames(cpdag)
    for (k in 2:K) {
        B_k <- buckets[[k]]
        B_k <- as.integer(B_k)
        
        # Set all undirected edges within the bucket to 0
        # This helps us quickly find parents
        cpdag_temp <- cpdag
        cpdag_temp[B_k, B_k] <- 0
        parents_B_k <- which(cpdag_temp[, B_k] != 0, arr.ind=TRUE)
        
        # Compute MLE of lambda block
        if (length(parents_B_k) != 0) {
            size_parents <- size(parents_B_k)
            if (size_parents[1] == 1) {
                parents_B_k <- as.integer(parents_B_k)
            }
            else {
                parents_B_k <- as.integer(parents_B_k[, 1])
            }
            # cat(B_k, "<-", parents_B_k, "\n")
            
            if (length(parents_B_k) == 1) {
                inv <- 1 / sigma_n[parents_B_k, parents_B_k]
            }
            else {
                inv <- pracma::pinv(sigma_n[parents_B_k, parents_B_k])
            }
            lambda[parents_B_k, B_k] <- inv %*% sigma_n[parents_B_k, B_k]
        }
    }
    
    return(lambda)
}


find_true_lambda <- function(gamma, cpdag, buckets) {
    p <- ncol(gamma)
    K <- length(buckets)
    lambda <- matrix(0, nrow=p, ncol=p)
    
    dimnames(lambda) <- dimnames(cpdag)
    for (k in 2:K) {
        B_k <- buckets[[k]]
        B_k <- as.integer(B_k)
        
        # Set all undirected edges within the bucket to 0
        # This helps us quickly find parents
        cpdag_temp <- cpdag
        cpdag_temp[B_k, B_k] <- 0
        parents_B_k <- which(cpdag_temp[, B_k] != 0, arr.ind=TRUE)
        
        # Compute MLE of lambda block
        if (length(parents_B_k) != 0) {
            size_parents <- size(parents_B_k)
            if (size_parents[1] == 1) {
                parents_B_k <- as.integer(parents_B_k)
            }
            else {
                parents_B_k <- as.integer(parents_B_k[, 1])
            }
            # cat(B_k, "<-", parents_B_k, "\n")
            
            # if (length(parents_B_k) == 1) {
            #     inv <- 1 / sigma_n[parents_B_k, parents_B_k]
            # }
            # else {
            #     inv <- pracma::pinv(sigma_n[parents_B_k, parents_B_k])
            # }
            
            
            I <- diag(length(B_k))
            gamma_B_k <- gamma[B_k, B_k]
            
            # print(B_k)
            
            # print(gamma_B_k)
            
            right_inv <- pracma::pinv(I - gamma_B_k)
            
            # print(right_inv)
            
            # print(parents_B_k)
            # print(gamma[parents_B_k, B_k])
            
            lambda[parents_B_k, B_k] <- gamma[parents_B_k, B_k] %*% right_inv
        }
    }
    
    return(lambda)
}



G_estimator <- function(G, lambda, A, Y) {
    V <- seq(1, ncol(G), by=1)
    V2 <- setdiff(V, A)
    D <- ancestors(Y, G[V2, V2])
    # print(D)
    # D_idx <- which(colnames(G) %in% D)
    D_idx <- as.integer(D)
    # print(D_idx)
    
    left <- lambda[A, D_idx]
    
    if (length(D) == 1) {
        tau_AY <- matrix(left, ncol=1)
        colnames(tau_AY) <- Y
        rownames(tau_AY) <- names(left)
        return(tau_AY)
    }
    
    I <- diag(1, nrow=length(D))
    right_inner <- pinv(I - lambda[D_idx, D_idx])
    dimnames(right_inner) <- dimnames(lambda[D_idx, D_idx])
    # print(right_inner)
    
    D_idx_right <- which(colnames(right_inner) %in% D)
    Y_idx_right <- which(colnames(right_inner) %in% Y)
    right <- matrix(right_inner[D_idx_right, Y_idx_right], nrow=length(D))
    
    # print(left)
    # print(right)
    rownames(right) <- D
    colnames(right) <- Y
    # print(right)
    # right[D, ] <- right[colnames(left), ]
    # print(right)
    
    tau <- left %*% right
    return(tau)
}

estimate_causal_effect <- function(data, A, Y, mpdag) {
    if (length(Y) != 1) {
        stop("The outcome variable Y should be one-dimensional")
    }
    
    identifiable <- is_identifiable(mpdag, A, Y)
    if (!identifiable) {
        print("Effect is not identifiable!")
        return(NULL)
    }
    
    buckets <- ordered_bucket_decomposition(mpdag)
    lambda_mle <- find_lambda_mle(data, mpdag, buckets)
    # print(lambda_mle)
    tau_AY <- G_estimator(mpdag, lambda_mle, A, Y)
    
    return(t(tau_AY))
}

find_true_causal_effect <- function(gamma, A, Y, mpdag) {
    if (length(Y) != 1) {
        stop("The outcome variable Y should be one-dimensional")
    }
    
    identifiable <- is_identifiable(mpdag, A, Y)
    if (!identifiable) {
        print("Effect is not identifiable!")
        return(NULL)
    }
    
    buckets <- ordered_bucket_decomposition(mpdag)
    lambda_mle <- find_true_lambda(gamma, mpdag, buckets)
    # print(lambda_mle)
    tau_AY <- G_estimator(mpdag, lambda_mle, A, Y)
    
    return(t(tau_AY))
}


