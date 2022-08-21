# tensorMCCA
Tensor-based Multiple Canonical Correlation Analysis

`R` package for the multiple canonical correlation analysis (MCCA) of data tensors measured on the same individuals/objects. The goal is to find canonical tensors that maximize the sum of covariances/correlations between canonical scores across all pairs of datasets. Various scaling contraints and orthogonality constraints can be imposed on canonical tensors during the optimization. The package implements several initialization methods (CCA-based, HOSVD, random) and optimization methods (block coordinate ascent, gradient-based with retraction by scaling, gradient-based with retraction by rotation).  Permutation tests can also be run to select the number of canonical components to retain. 

**To install the `R` package:**
```
library(devtools)
install_github("https://github.com/ddegras/tensorMCCA", subdir = "tensorMCCA")
```

**Main functions:**

-  `mcca.cov`, `mcca.cor`
-  `mcca.init.cca`, `mcca.init.svd`, `mcca.init.random` 
-  `permutation.test`
