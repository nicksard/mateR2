
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mateR2

<!-- badges: start -->
<!-- badges: end -->

mateR2 is an MCMC-driven simulation engine for generating breeding
matrices and single-cohort pedigrees under precise demographic
constraints. It treats the mating pool of a dioecious species as a
bipartite network, which preserves the joint covariance between male and
female reproductive output. It is aimed at the skewed, “sweepstakes”
reproductive systems common in wild aquatic populations (e.g. salmon and
sturgeon), and at evaluating how pedigree reconstruction frameworks
behave under them.

The pipeline runs in three stages:

1.  **Metropolis-Hastings MCMC.** `generate_map_table()` searches a
    compressed state space of discrete mating-block configurations for a
    Mate-Pair Summary Table hitting your targets for parental abundance
    ($N_P$), operational sex ratio ($SR$) and mean number of mates
    ($MM$).
2.  **Degree-preserving edge swapping.** `mp_table_to_matrix()` expands
    that table into an individual-level binary matrix, and
    `randomize_mating_structure()` rewires it to dissolve the
    block-diagonal structure *without altering any individual’s mate
    count*.
3.  **Fecundity and sampling.** `brd.mat.fitness()` converts the binary
    matrix **M** into a weighted breeding matrix **K**, and
    `mat.sub.sample()` draws a juvenile sample without replacement to
    emulate field sampling effort.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("nicksard/mateR2", build_vignettes = TRUE)
```

Optional extras: `coda` for `run_mcmc_chains()`, `igraph` for
`count_network_components()`.

## A worked example

``` r
library(mateR2)

# Stage 1: 100 parents, 2:1 sex ratio, 2 mates each
fit <- generate_map_table(
  Np_target = 100, sr_target = 2.0, mean_mates_target = 2.0,
  max_males_per_female = 10, max_females_per_male = 10,
  np_weight = 50, sr_weight = 50, mm_weight = 50, decay_constant = -0.05,
  n_iter = 50000, burn_in = 5000, thin = 10, seed = 42
)
#> MCMC finished. Overall Acceptance Rate: 2.97%
unlist(fit$map_stats)
#>         Np         SR  MeanMates MaxLogProb 
#>  99.000000   2.000000   2.000000  -2.852517
```

``` r
# Stage 2: expand to individuals, then rewire at 25% mixing intensity
M <- mp_table_to_matrix(fit$map_table)
M <- randomize_mating_structure(M, I = 0.25)
M[1:4, 1:6]
#>      F028 F015 F020 F025 F011 F001
#> M066    1    0    0    0    0    0
#> M015    0    0    0    0    0    0
#> M040    0    0    0    0    0    0
#> M064    0    0    0    0    0    0
```

Rewiring dissolves the block structure while leaving every marginal
untouched:

``` r
mat.stats(M)[, c("num.males", "num.females", "overall.mean.mates")]
#>   num.males num.females overall.mean.mates
#> 1        66          33                  2
count_network_components(mp_table_to_matrix(fit$map_table))  # before
#> [1] 26
count_network_components(M)                                  # after
#> [1] 16
```

``` r
# Stage 3: allocate fecundity, then sample 500 juveniles
K <- brd.mat.fitness(M, min.fert = 1000, max.fert = 5000, type = "uniform")
ped <- convert2ped(mat.sub.sample(K, num_offspring = 500))
head(ped, 3)
#>    off  mom  dad
#> 1 off1 F028 M066
#> 2 off2 F028 M066
#> 3 off3 F028 M066
```

``` r
# Sibship structure and parental detection in the sample
sampled <- ped2mat(ped)
sib.stats(sampled)[, c("mean_FS", "mean_MHS", "mean_PHS")]
#>   mean_FS mean_MHS mean_PHS
#> 1   7.416   10.308    2.648
parent.class.stats(sampled)[, c("detected", "singletons", "doubletons")]
#>          detected singletons doubletons
#> maternal       33          0          0
#> paternal       65          4          5
#> overall        98          4          5
```

## Reproducibility

The sampler draws from R’s random number generator, so `set.seed()`
governs a run, and passing `seed` to `generate_map_table()` does the
same. Identical seeds produce identical chains across platforms.

`run_mcmc_chains()` runs several dispersed chains and returns `coda`
traces with Gelman-Rubin and effective sample size diagnostics.

## Getting started

The quick-start vignette walks through the same pipeline in more detail:

``` r
vignette("quickstart-breeding-sim", package = "mateR2")
```

## Citation

``` r
citation("mateR2")
```
