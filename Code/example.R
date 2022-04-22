##
## Simple example to work through
##

source("sampler.R")
source("mle.R")

set.seed(3)

num_nodes <- 10
poss_degrees <- c(2, 3, 4, 5)
avg_degree <- 3

dag_object <- sample_dag(num_nodes, avg_degree)
X <- sample_data(5000, dag_object$dag, err=error_dist_1)
sigma_n <- sample_covariance(X)

A <- c(4, 6)
Y <- 10

cause_eff <- estimate_causal_effect(X, A, Y, dag_object$cpdag_amat)
true_eff <- find_true_causal_effect(dag_object$dag_amat, A, Y, dag_object$cpdag_amat)

eff2_eff <- eff2::estimateEffect(X, A, Y, t(dag_object$cpdag_amat))
adjustment_eff <- pcalg::ida(A, Y, sigma_n, dag_object$cpdag_amat, method="optimal", type="cpdag")
rrc_eff <- pcalg::jointIda(A, Y, sigma_n, dag_object$cpdag_amat, technique="RRC", type="cpdag")
mcd_eff <- pcalg::jointIda(A, Y, sigma_n, dag_object$cpdag_amat, technique="MCD", type="cpdag")


print(cause_eff)
print(true_eff)
print(eff2_eff)
print(adjustment_eff)
print(rrc_eff)
print(mcd_eff)

print(compute_error(cause_eff, true_eff))
print(compute_error(adjustment_eff[, 1], true_eff))
print(compute_error(rrc_eff[, 1], true_eff))
print(compute_error(mcd_eff[, 1], true_eff))

