## Testing data
##

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DREAM4")

library("DREAM4")
library("SummarizedExperiment")


data("dream4_010_05")
names(assays(dream4_010_01))
expressionData <- assays(dream4_010_01)$simulated
names(metadata(dream4_010_01))
goldStandardMatrix <- metadata(dream4_010_01)$goldStandardAdjacencyMatrix

goldStandardMatrix <- metadata(dream4_010_05)$goldStandardAdjacencyMatrix
