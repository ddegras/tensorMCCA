##########################################
# FUNCTIONS TO APPROXIMATE RANK-1 TENSORS 
# UNDER ORTHOGONALITY CONSTRAINTS 
##########################################


# min sum_i,k || v_ik - v0_ik ||^2

# sum_i < v_i, vll_i > = 0, ll = 1,...,l
# (1/m) sum_i || v_i ||^2 = 1



# sum_i w_ik(i)  < v_ik(i) , vll_ik(i) > = 0 for ll = 1, ..., l
# w_ik(i) = prod_(k!=k(i)) < v_ik , vll_ik > 

# approx.rk1.tsnr.rk1.ortho <- function(v, ortho, maxit = 1000L, tol = 1e-6)
# {
# stopifnot(is.list(v) && is.list(ortho))
# ortho <- as.matrix(ortho)
# stopifnot(length(v) == nrow(ortho))
# m <- length(v)
# northo <- ncol(ortho)
# blocks <- expand.grid(lapply(v, seq_along))
# nblocks <- nrow(blocks)
# v0 <- v





	
# }