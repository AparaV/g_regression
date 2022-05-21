##
## DAG to CPDAG

require(pcalg)

get_cpdag <- function(dag_amat, verbose=T) {
    dag_amat[dag_amat != 0] <- 1
    p <- nrow(dag_amat)
    
    pdag <- dag_amat + abs(t(dag_amat))
    pdag[pdag != 0] <- 1
    
    # Orient v-structures
    edges <- which(dag_amat == 1, arr.ind=T)
    num_edges <- nrow(edges)
    for (i in seq(1, num_edges)) {
        A <- edges[i, 1]
        B <- edges[i, 2]
        
        edges_to_B <- which(dag_amat[, B] == 1, arr.ind=T)
        for (j in 1:length(edges_to_B)) {
            C <- edges_to_B[j]
            if (C == A) {
                next
            }
            
            if (dag_amat[A, C] == 0 && dag_amat[C, A] == 0) {
                if (verbose) {
                    cat("Orienting", A, "->", B, "<-", C, "\n")
                }
                pdag[B, A] <- 0
                pdag[B, C] <- 0
            }
        }
    }
    
    # Complete Meek's rules
    cpdag <- addBgKnowledge(t(pdag), verbose=verbose, checkInput=F)
    
    return(t(cpdag))
}