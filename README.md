# tensorMCCA
Tensor-based Multiple Canonical Correlation Analysis

`R` package for the multiple canonical correlation analysis (MCCA) of data tensors measured on the same individuals/objects. The goal is to find canonical tensors that maximize the sum of covariances/correlations between canonical scores across all pairs of datasets. Various orthogonality constraints can be imposed to find higher-order canonical components. The functions of the package implement initialization algorithms (random, HOSVD), optimization algorithms (block coordinate ascent, gradient), as well as permutation tests to select the number of canonical components to retain. 

**To install the `R` package:**
```
library(devtools)
install_github("https://github.com/ddegras/tensorMCCA", subdir = "tensorMCCA")
```
