#
# Working with DREAM4 data
#

library(tidyverse)
library(pcalg)

source("mle.R")
source("sampler.R")

estimate_s_ij <- function(s, tau, A, Y) {
    s_ij <- s[Y] - (s[A[1]] * tau[1] + s[A[2]] * tau[2])
    return(s_ij)
}

observed_data_dir <- "../Data/DREAM4_InSilico_Size10/insilico_size10_5/"
gold_standard_dir <- "../Data/DREAM4_Challenge2_GoldStandards/Size 10/Size 10 bonus round/"

X1 <- as.matrix(read_tsv(paste0(observed_data_dir, "time_1.tsv")))
X2 <- as.matrix(read_tsv(paste0(observed_data_dir, "time_2.tsv")))
X3 <- as.matrix(read_tsv(paste0(observed_data_dir, "time_3.tsv")))
X4 <- as.matrix(read_tsv(paste0(observed_data_dir, "time_4.tsv")))
X5 <- as.matrix(read_tsv(paste0(observed_data_dir, "time_5.tsv")))

s <- read_tsv(paste0(observed_data_dir, "insilico_size10_5_wildtype.tsv"))
s <- data.frame(s)
double_knockout <- read_tsv(paste0(gold_standard_dir, "insilico_size10_5_nonoise_dualknockouts.tsv"))
double_knockout <- data.frame(double_knockout)

network <- matrix(0, nrow=10, ncol=10)
network[1, 2] <- 1
network[3, 4] <- 1
network[4, 2] <- 1
network[5, 6] <- 1
network[6, 8] <- 1
network[6, 9] <- 1
network[7, 2] <- 1
network[7, 5] <- 1
network[8, 7] <- 1
network[8, 9] <- 1
network[9, 1] <- 1
network[10, 3] <- 1
rownames(network) <- c(1:10)
colnames(network) <- c(1:10)

network_1 <- network
network_1[5, 6] <- 0

network_2 <- network
network_2[6, 8] <- 0

network_3 <- network
network_3[8, 7] <- 0

network_4 <- network
network_4[7, 5] <- 0

networks <- list(network_1, network_2, network_3, network_4)

A <- list(c(6, 8), c(7, 8), c(8, 10), c(8, 5), c(8, 9))
Y <- c(1:10)


reps <- 100
G_na_errors <- data.frame(matrix(nrow=length(networks), ncol=reps))
adj_na_errors <- data.frame(matrix(nrow=length(networks), ncol=reps))
G_errors <- data.frame(matrix(nrow=length(networks), ncol=reps))
rrc_errors <- data.frame(matrix(nrow=length(networks), ncol=reps))

for (rep_id in 1:reps) {

cov_1 <- cov(X1[sample(nrow(X1), replace=T), ])
cov_2 <- cov(X2[sample(nrow(X2), replace=T), ])
cov_3 <- cov(X3[sample(nrow(X3), replace=T), ])
cov_4 <- cov(X4[sample(nrow(X4), replace=T), ])
cov_5 <- cov(X5[sample(nrow(X5), replace=T), ])
cov <- (cov_1 + cov_2 + cov_3 + cov_4 + cov_5) / 5
rownames(cov) <- c(1:10)
colnames(cov) <- c(1:10)


errors <- data.frame(matrix(nrow=0, ncol=7))
colnames(errors) <- c("Network", "G", "eff2", "adj", "rrc", "baseline", "denominator")

for (k in 1:length(networks)) {
    for (i in 1:length(A)) {
        for (j in 1:length(Y)) {
            
            if (Y[j] %in% A[[i]]) {
                s_ij_est_G <- 0
                s_ij_est_adj <- 0
                s_ij_est_rrc <- 0
                # next
            }
            else {
                cause_eff <- estimate_causal_effect(NULL, A[[i]], Y[j], networks[[k]], cov)
                eff2_eff <- eff2:::.estimateEffect(cov, A[[i]], Y[j], networks[[k]])
                adjustment_eff <- pcalg::ida(A[[i]], Y[j], cov, networks[[k]], method="optimal", type="pdag")
                rrc_eff <- pcalg::jointIda(A[[i]], Y[j], cov, networks[[k]], technique="RRC", type="pdag")[, 1]
                
                s_ij_est_G <- estimate_s_ij(s, cause_eff, A[[i]], Y[j])#[[1]]
                s_ij_est_eff2 <- estimate_s_ij(s, eff2_eff, A[[i]], Y[j])#[[1]]
                s_ij_est_adj <- estimate_s_ij(s, adjustment_eff, A[[i]], Y[j])#[[1]]
                s_ij_est_rrc <- estimate_s_ij(s, rrc_eff, A[[i]], Y[j])#[[1]]
            }
            
            baseline_est <- s[j]
            true <- double_knockout[i, j]
            
            errors_k_G <- (s_ij_est_G - true)^2
            errors_k_eff2 <- (s_ij_est_eff2 - true)^2
            errors_k_adj <- (s_ij_est_adj - true)^2
            errors_k_rrc <- (s_ij_est_rrc - true)^2
            errors_k_baseline <- (baseline_est - true)^2
            normalizing <- true^2
            
            errors[nrow(errors) + 1, ] <- c(k, errors_k_G, errors_k_eff2, errors_k_adj, errors_k_rrc, errors_k_baseline, normalizing)
            
        }
    }
}

    for (k in 1:length(networks)) {
        results_k <- errors[which(errors$Network == k), ]
        
        results_k_na <- results_k[-which(is.na(results_k$adj)), ]
        G_na <- sum(results_k_na$G) / sum(results_k_na$denominator)
        adj_na <- sum(results_k_na$adj) / sum(results_k_na$denominator)
        
        G <- sum(results_k$G) / sum(results_k$denominator)
        rrc <- sum(results_k$rrc) / sum(results_k$denominator)
        
        G_na_errors[k, rep_id] <- G_na
        adj_na_errors[k, rep_id] <- adj_na
        G_errors[k, rep_id] <- G
        rrc_errors[k, rep_id] <- rrc
        
    }

}

N <- nrow(X1) + nrow(X2) + nrow(X3) + nrow(X4) + nrow(X5)
z <- qnorm(0.975)

G_na_se <- sqrt(apply(G_na_errors, 1, var)) * z / sqrt(N) * 100
adj_na_se <- sqrt(apply(adj_na_errors, 1, var)) * z / sqrt(N) * 100
G_se <- sqrt(apply(G_errors, 1, var)) * z / sqrt(N) * 100
rrc_se <- sqrt(apply(rrc_errors, 1, var)) * z / sqrt(N) * 100

print(G_na_se)
print(adj_na_se)
print(G_se)
print(rrc_se)
