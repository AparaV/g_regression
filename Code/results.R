##
## Analyze results
##

gm_mean = function(x, na.rm=TRUE){
    if (na.rm) {
        x <- x[which(!is.na(x))]
    }
    exp(sum(log(x[x > 0])) / length(x))
}


poss_V <- c(20, 50, 100)
poss_A <- c(1, 2, 3, 4)
poss_n <- c(100, 1000)
num_reps <- 1000

for (V in poss_V) {
    for (A in poss_A) {
        for (n in poss_n) {
            
            V <- 20
            A <- 4
            n <- 1000
            
            path_name <- paste0("../Results/V", V, "_A", A, "_n", n, ".rds")
            
            cat("V = ", V, "A = ", A, "n = ", n, "\n")
            
            results <- readRDS(path_name)
            errors <- results$errors
            
            errors$adj_ratio <- errors$adj_err^2 / errors$G_err^2
            errors$mcd_ratio <- errors$mcd_err^2 / errors$G_err^2
            errors$rrc_ratio <- errors$rrc_err^2 / errors$G_err^2
            
            adj_ratio_gm <- gm_mean(errors$adj_ratio, na.rm=T)
            mcd_ratio_gm <- gm_mean(errors$mcd_ratio, na.rm=T)
            rrc_ratio_gm <- gm_mean(errors$rrc_ratio, na.rm=T)
            
            adj_missing <- length(which(is.na(errors$adj_err)))
            mcd_missing <- length(which(is.na(errors$mcd_err)))
            rrc_missing <- length(which(is.na(errors$rrc_err)))
            
            cat("Errors:\n")
            cat("\tadj.0:", adj_ratio_gm, "\n")
            cat("\tida.M:", mcd_ratio_gm, "\n")
            cat("\tida.R:", rrc_ratio_gm, "\n")
            
            cat("\nMissing:\n")
            cat("\tadj.0:", adj_missing / num_reps * 100, "%\n")
            cat("\tida.M:", rrc_missing / num_reps * 100, "%\n")
            cat("\tida.R:", mcd_missing / num_reps * 100, "%\n")
            
            cat("========\n")
        }
    }
}


library(tidyverse)
library(ggplot2)


errors_ratio <- errors[c("adj_ratio", "rrc_ratio", "mcd_ratio")]
long_data <- errors_ratio %>% gather(key="Method", value="Ratio")
# long_data_na <- long_data[which(!is.na(long_data$Ratio)), ]


p <- ggplot(long_data, aes(Method, Ratio)) +
    geom_violin(trim=T, na.rm=T, scale="width", kernel="triangular")
p
