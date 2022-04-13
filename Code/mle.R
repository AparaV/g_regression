##
## Transform from Gamma to Lambda space
##

require(pracma)

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
    
    dimnames(sigma_n) <- dimnames(lambda) <- dimnames(cpdag_matrix)
    
    for (k in 2:K) {
        B_k <- buckets[[k]]
        
        # Set all undirected edges within the bucket to 0
        # This helps us quickly find parents
        cpdag_temp <- cpdag_matrix
        cpdag_temp[B_k, B_k] <- 0
        parents_B_k <- which(cpdag_temp[, B_k] != 0, arr.ind=TRUE)
        
        # Compute MLE of lambda block
        if (length(parents_B_k) != 0) {
            lambda[parents_B_k, B_k] <- pracma::pinv(
                sigma_n[parents_B_k, parents_B_k]) %*% sigma_n[parents_B_k, B_k]
        }
    }
    
    return(lambda)
}



