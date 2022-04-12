library(pcalg)
library(graph)


set.seed(33)

size <- 6
prob <- 0.5

d <- randomDAG(size, prob, lB=1, uB=1)

dev.off()
plot(d)

cpdag <- dag2cpdag(d)

dev.off()
plot(cpdag)


d2 <- matrix(c(c(0, 0, 0, 0),
               c(1, 0, 0, 0),
               c(1, 1, 0, 0),
               c(0, 1, 1, 0)),
             nrow=4, ncol=4)
d2 <- as(d2, "graphNEL")

dev.off()
plot(d2)

cpdag2 <- dag2cpdag(d2)

dev.off()
plot(cpdag2)
