tnsr.rk1 <- function(x, maxit = 100, tol = 1e-6)
{
stopifnot(is.numeric(x))
stopifnot(maxit >= 0 && tol >= 0)
if (is.vector(x)) return(list(x))
if (all(x == 0)) return(lapply(dim(x), numeric))
d <- length(dim(x)) 
if (d == 3L) return(tnsr3d.rk1(x, maxit, tol))
svdx <- hosvd(x, 1) 
nrmv <- as.numeric(svdx$core)
if (d == 2L || maxit == 0) {
	v <- lapply(svdx$vectors, "*", y = nrmv^(1/d))
	return(v)
}
v <- svdx$vectors
v[[d]] <- v[[d]] * nrmv
kronv <- Reduce(kronecker, v)
for (it in 1:maxit) {
	kronv.old <- kronv
	for (k in 1:d) {
		v[[k]] <- tvec.prod(x, v[-k], (1:d)[-k])
		if (k < d) v[[k]] <- v[[k]] / sqrt(sum(v[[k]]^2))
	}	
	kronv <- Reduce(kronecker, v)
	e <- sqrt(sum((kronv - kronv.old)^2))
	if (e <= tol * max(nrmv, nrmv.old, 1)) break
}

s <- c(rep.int(nrmv^(1/d), d - 1), nrmv^(1/d - 1))
for (k in 1:d) v[[k]] <- v[[k]] * s[k] 
v	
}

