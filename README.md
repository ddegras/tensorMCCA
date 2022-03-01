# tensorMCCA
Tensor-based Multiple Canonical Correlation Analysis

A package for multiple canonical correlation analysis (MCCA) with data tensors measured on the same individuals/objects. The goal is to find canonical tensors that maximize sums of covariances or sums of correlations between pairs of datasets under orthogonality constraints. The functions in the package implement initialization algorithms (random, HOSVD), optimization algorithms (block coordinate ascent, gradient), as well as permutation tests to select the number of canonical components to retain. 
