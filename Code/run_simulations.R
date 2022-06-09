##
## Simulations script
##

library(optparse)

source("sampler.R")
source("mle.R")
source("identifiablity.R")

option_list = list(
    make_option(c("--num_nodes"), type="integer", default=10, 
                help="Size of the graph", metavar="integer"),
    make_option(c("--treatment_size"), type="integer", default=2, 
                help="Number of variables to intervene", metavar="integer"),
    make_option(c("--n"), type="integer", default=1000, 
                help="Number of data points", metavar="integer"),
    make_option(c("--num_reps"), type="integer", default=100, 
                help="Number of repetitions", metavar="integer"),
    make_option(c("--output"), type="character", default=NULL,
                help="Path to store output", metavar="character")
)
arg.parser = OptionParser(option_list=option_list)
args = parse_args(arg.parser)


# Parameters
num_nodes <- args$num_nodes
A_len <- args$treatment_size
n <- args$n
num_reps <- args$num_reps
output <- args$output


# num_nodes <- 50
# A_len <- 3
# n <- 1000


results_mat <- matrix(NA, nrow=num_reps, ncol=5)

seed_counter <- 0
set.seed(seed_counter)
iters <- 0
while (iters < num_reps) {
    
    seed_counter <- seed_counter + 1
    # if (seed_counter == 443) {
    #     next
    # }
    # seed_counter <- 1
    set.seed(seed_counter)
    
    cat("Iteration", iters + 1, "( seed =", seed_counter, "):\n")
    
    out <- tryCatch({
            # Sample a DAG
            cat("\tSampling a DAG\n")
            poss_degrees <- c(2, 3, 4, 5)
            avg_degree <- sample(poss_degrees, 1)
            dag_object <- sample_dag(num_nodes, avg_degree)
            
            # Sample data from DAG
            cat("\tSampling data\n")
            X <- sample_data(n, dag_object$dag, err=error_dist_1)
            sigma_n <- sample_covariance(X)
            
            # Sample treatments and outcomes
            cat("\tSampling treatment and outcome variables\n")
            AY <- sample_treatment_outcome(dag_object$dag_amat, dag_object$cpdag_amat, A_len)
            if (is.null(AY)) {
                cat("\t\tNo identified pair found!")
                cat("\t\tRepeating iteration", iters + 1, "\n")
                cat("--------------------\n")
                # There was no identifiable treatment, outcome pair found
                # return(FALSE)
                next
            }
            A <- AY$A
            Y <- AY$Y
            
            # Estimate effects using various methods
            cat("\tEstimating effects\n")
            cause_eff <- estimate_causal_effect(X, A, Y, dag_object$cpdag_amat)
            true_eff <- find_true_causal_effect(dag_object$dag_amat, A, Y, dag_object$cpdag_amat)
            adjustment_eff <- pcalg::ida(A, Y, sigma_n, dag_object$cpdag_amat, method="optimal", type="cpdag")
            rrc_eff <- pcalg::jointIda(A, Y, sigma_n, dag_object$cpdag_amat, technique="RRC", type="cpdag")
            mcd_eff <- pcalg::jointIda(A, Y, sigma_n, dag_object$cpdag_amat, technique="MCD", type="cpdag")
            
            # Compute and store errors
            cat("\tComputing L2 errors\n")
            results_iter <- rep(NA, 5)
            results_iter[1] <- seed_counter
            results_iter[2] <- compute_error(cause_eff, true_eff)
            if (is_effect_same(adjustment_eff)) {
                results_iter[3] <- compute_error(adjustment_eff[, 1], true_eff)
            }
            if (is_effect_same(rrc_eff)) {
                results_iter[4] <- compute_error(rrc_eff[, 1], true_eff)
            }
            if (is_effect_same(mcd_eff)) {
                results_iter[5] <- compute_error(mcd_eff[, 1], true_eff)
            }
            
            cat("\tFinished iteration", iters + 1, "\n")
            cat("--------------------\n")
            iters <- iters + 1
            results_mat[iters,] <- results_iter
            
            # return(TRUE)
        },
        error=function(cond){
            cat("\tRan into the following error:\n\t")
            print(cond)
            cat("\tTrying next seed\n")
            return(FALSE)
        }
    )
}

cat("Finished all simulations.\n")
cat("Saving results in", output, "\n")

results_df <- data.frame(results_mat)
colnames(results_df) <- c("Seed", "G_err", "adj_err", "rrc_err", "mcd_err")
results <- list(
    errors=results_df,
    num_nodes=num_nodes,
    A_len=A_len,
    n=n,
    num_reps=num_reps
)

saveRDS(results, output)

