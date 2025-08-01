# Stacked SVD or SVD Stacked?

## R package functionality

To install the R package, run `devtools::install_github("phillipnicol/stackedSVD")`. The following demonstrates the functionality using a simulated example with 3 matrices and a rank 1 signal:

```
set.seed(134609)
n <- 500
d <- 1000

library(stackedSVD)

M <- 3

theta <- c(1, 1.25, 1.5)

X.list <- list() 
v <- rnorm(n=d)
v <- v / sqrt(sum(v^2))
for(m in 1:M) {
  u <- rnorm(n=n)
  u <- u / sqrt(sum(u^2))
  
  X.list[[m]] <- theta[m]*outer(u,v) + matrix(rnorm(n*d,sd=1/sqrt(d)),
                                         nrow=n,ncol=d)
}

##Run methods 
v.stacksvd <- stackedSVD(X.list, rank=1)
cor.stacksvd <- cor(v.stacksvd, v)^2

v.svdstack <- SVDstacked(X.list, rank=1)
cor.svdstack <- cor(v.svdstack, v)^2

v.stacksvdw <- stackedSVDWeighted(X.list, rank=1)
cor.stacksvdw <- cor(v.stacksvdw, v)^2

v.svdstackw <- SVDstackedWeighted(X.list, rank=1)
cor.svdstackw <- cor(v.svdstackw, v)^2

all.cor <- c(cor.stacksvd, cor.svdstack, cor.stacksvdw, cor.svdstackw)
names(all.cor) <- c("stackedSVD", "SVDstacked", "stackedSVD (weighted)", "SVDstacked (weighted)")

print(all.cor)
```

## Replicating analysis

To generate the results of Figure 3, navigate to `analysis/sim/code` and run `different_c_theta_final.R` (a sample `sbatch` script is also provided). Then `plotting.R` generates the plots. 

To replicate the results in Figures 4 and 5, navigate to `analysis/python_plots` and run the Python file `simulations.py`. Note the optional flags `--test` and `--help` to test the code using a smaller number of parameter choices. 

To replicate the single-cell RNA-seq results, navigate to `analysis/zhengmix_ambient/code` and run `single-cell_full.R`. 

## Reference

Baharav, T.Z., Nicol, P.B., Irizarry, R.A., and Ma, R. (2025). Stacked SVD or SVD Stacked? A Random Matrix Theory perspective on data integration. [ArXiv](https://arxiv.org/abs/2507.22170)
