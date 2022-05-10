##
## Analyze results
##


library(tidyverse)
library(ggplot2)
library(latex2exp)
library(ggh4x)

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

full_df <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(full_df) <- c("V", "A", "n", "adj_ratio", "mcd_ratio", "rrc_ratio")


ggplot_theme <- theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          strip.background = element_rect(colour="black", fill="white"),
          plot.title = element_text(hjust = 0.5))

for (V in poss_V) {
    for (A in poss_A) {
        for (n in poss_n) {
            
            path_name <- paste0("../Results/V", V, "_A", A, "_n", n, ".rds")
            
            cat("V = ", V, "A = ", A, "n = ", n, "\n")
            
            results <- readRDS(path_name)
            errors <- results$errors
            
            errors$adj_ratio <- errors$adj_err^2 / errors$G_err^2
            errors$mcd_ratio <- errors$mcd_err^2 / errors$G_err^2
            errors$rrc_ratio <- errors$rrc_err^2 / errors$G_err^2
            
            
            errors_ratio <- errors[c("adj_ratio", "rrc_ratio", "mcd_ratio")]
            long_data <- errors_ratio %>% gather(key="Method", value="Ratio")
            
            # this_df <- errors[, c("adj_ratio", "mcd_ratio", "rrc_ratio")]
            this_df <- long_data
            this_df$V <- V
            this_df$A <- A
            this_df$n <- n
            
            full_df <- rbind(full_df, this_df)
            
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
            
            plot_title <- TeX(paste0("$|V|$ = ", V, ", $|A|$ = ", A, " (", n, " points)"))
            y_breaks <- c(1e-3, 1, 1e3, 1e6, 1e9)
            y_labels <- c(TeX("$10^{-3}$"), TeX("$1$"), TeX("$10^3$"), TeX("$10^6$"), TeX("$10^9$"))
            
            p <- ggplot(this_df, aes(Method, Ratio)) +
                geom_hline(yintercept=1, color="black", size=0.7, linetype="dashed") +
                geom_violin(trim=T, na.rm=T, scale="width", kernel="triangular") +
                stat_summary(fun="mean", geom="point") +
                # These make a "+" sign at the median
                stat_summary(fun="median", geom="crossbar", width=0.2, size=0.2) +
                stat_summary(fun="median", geom="crossbar", width=0.02, size=2) +
                ggplot_theme +
                xlab("Method") + ylab("Relative error") + ggtitle(plot_title) +
                scale_y_log10(breaks=y_breaks, labels=y_labels) +
                scale_x_discrete(labels=c("Optimal\nAdjustment", "IDA\n(MCD)", "IDA\n(RRC)"))
            
            fname <- paste0("../Figures/violin_", V, "_", A, "_", n, ".pdf")
            ggsave(fname, units="in", width=5, height=5, dpi=200)
        }
    }
}


# long_data_na <- long_data[which(!is.na(long_data$Ratio)), ]

plot_title <- TeX(paste0("$|V|$ = ", V, ", $|A|$ = ", A, " (", n, " points)"))
y_breaks <- c(1e-3, 1, 1e3, 1e6, 1e9)
y_labels <- c(TeX("$10^{-3}$"), TeX("$1$"), TeX("$10^3$"), TeX("$10^6$"), TeX("$10^9$"))

p <- ggplot(full_df, aes(Method, Ratio)) +
    geom_hline(yintercept=1, color="black", size=0.7, linetype="dashed") +
    geom_violin(trim=T, na.rm=T, scale="width", kernel="triangular") +
    stat_summary(fun="mean", geom="point") +
    # These make a "+" sign at the median
    stat_summary(fun="median", geom="crossbar", width=0.2, size=0.2) +
    stat_summary(fun="median", geom="crossbar", width=0.02, size=2) +
    ggplot_theme +
    xlab("Method") + ylab("Relative error") +
    scale_y_log10(breaks=y_breaks, labels=y_labels) +
    scale_x_discrete(labels=c("Optimal\nAdjustment", "IDA\n(MCD)", "IDA\n(RRC)")) +
    facet_nested(A ~ V + n, labeller=label_both)

ggsave("../Figures/violin_ensemble.pdf", units="in", width=12, height=10, dpi=200)

