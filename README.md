# Efficient least squares for total causal effect estimation under linearity and causal sufficiency

Implementation of F. R. Guo, E. Perković (2022). Efficient Least Squares for Estimating Total Effects under Linearity and Causal Sufficiency.

This is a part of my Ph.D. Preliminary Examination at the Department of Statistics at University of Washington.

## Setup

All simulations were performed in R v4.1.2 using `pcalg` v2.7.5. Standard graphical libraries are required along with `tidyverse` and `ggplot2`.

```R
install.packages("pcalg")
install.packages("igraph")
install.packages("RBGL")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("latex2exp")
install.packages("ggh4x")
install.packages("optparse")
```

## Implementation

- `Code\bucket_decomposition.R` performs partial causal ordering.
- `Code\cpdag.R` is my own implementation of coverting a DAG into a CPDAG since the clusters owned by the Department of Statistics at UW had issues running `pcalg:::dag2cpdag`.
- `Code\identifiability.R` contains implementation of finding ancestors, descendants, directed paths, etc. Most importantly, it checks for the causal effect identifiability condition.
- `Code\mle.R` computes the G-regression estimator.
- `Code\sampler.R` samples random DAGs, treatment sets and outcome variable, and data points.


## Running

### Synthetic datset

All simulations are seeded. So they should be reproducible.

```shell
$ cd Code
$ Rscript run_simulations.R --num_nodes 10 --treatment_size 2 --n 1000 --num_reps 5 --output ../Results/test.rds
```

The above command uses the true CPDAG to estimate causal effects. To replicate results in Appendix G where the CPDAG is estimated using Greedy Equivalence Search, replace `run_simulations.R` with `run_simulations_unknown_cpdag.R`.

Arguments:
- `--num_nodes` - Number of nodes in the DAG
- `--treatment_size` - Number of variables in the treatment set
- `--n` - Number of iid data points
- `--num_reps` - Number of replications
- `--output` - File to store output

The cluster jobs in `Code\Slurm\synthetic_sims\` can automate replicating results in the report. The cluster jobs in `Code\Slurm\synthetic_sims_ges\` can help replicate results in Appendix G. Plots are generated using `Code\results.R`. You might need to modify filenamse to accomodate your results.

### DREAM4 dataset

The dataset was obtained from https://www.synapse.org/#!Synapse:syn3049712/wiki/74630. This is locally available in the `Data\` directory. Table 3 in the report can be generated using the following commands:

```shell
$ cd Code
$ Rscript dream4_results.R
$ Rscript dream4_bootstrap.R
```

`dream4_results.R` generates the normalized squared errors. `dream4_bootstrap.R` computes standard errors using a bootstrap approach described in Appendix E.