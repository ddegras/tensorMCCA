# tensorMCCA
Tensor-based Multiple Canonical Correlation Analysis

`R` package for the multiple canonical correlation analysis (MCCA) of data tensors (multidimensional arrays) measured on the same individuals/objects. The goal is to find canonical weight tensors such that the sum of covariances/correlations between canonical scores is maximized across all pairs of datasets. Various scaling contraints and orthogonality constraints can be imposed on canonical weights during the optimization. The package implements several initialization methods (CCA-based, HOSVD, random) and optimization methods (block coordinate ascent, gradient-based with retraction by scaling, gradient-based with retraction by rotation).  Permutation tests can also be run to select the number of canonical components to retain. The package also features bootstrap methods, including basic-, percentile-, and normal bootstrap confidence intervals.   

**To install the `R` package:**
```
library(devtools)
install_github("https://github.com/ddegras/tensorMCCA", subdir = "tensorMCCA")
```

**Main functions:**

-  MCCA optimization: `mcca.cov`, `mcca.cor`
-  Initialization: `mcca.init.cca`, `mcca.init.svd`, `mcca.init.random` 
-  Permutation tests: `permutation.test`
-  Bootstrap confidence intervals: `mcca.boot`, `mcca.boot.ci`
