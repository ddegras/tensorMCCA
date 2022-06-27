tnsr.rk1 <- function(x, scale = FALSE, maxit = 100, tol = 1e-6)
{
stopifnot(is.numeric(x))
stopifnot(maxit >= 0 && tol >= 0)
eps <- 1e-14
if (is.vector(x)) {
	if (scale) {
		nrm <- sqrt(sum(x^2))
		if (nrm > eps) {
			return(list(x / nrm))
		} else {
			return(list(rep.int(1 / sqrt(length(x)), length(x))))
		}
	} else {
		return(list(x))
	}
}
if (all(x == 0)) {
	if (scale) {
		return(lapply(dim(x), 
			function(xx) rep.int(1/sqrt(length(x)), length(x))))		
	} else {
		return(lapply(dim(x), numeric))
	}
}
d <- length(dim(x)) 
if (d == 3L) { return(tnsr3d.rk1(x, scale, maxit, tol)) }
svdx <- hosvd(x, 1) 
nrmv <- if (scale) 1 else as.numeric(svdx$core)
v <- svdx$vectors
if (d == 2L || maxit == 0) {
	if (scale) {
		return(v)
	} else {
		return(lapply(v, "*", y = nrmv^(1/d)))		
	}
}
if (!scale) v[[d]] <- v[[d]] * nrmv
kronv <- Reduce(kronecker, v)
for (it in 1:maxit) {
	nrmv.old <- nrmv
	kronv.old <- kronv
	for (k in 1:d) {
		v[[k]] <- tvec.prod(x, v[-k], (1:d)[-k])
		if (k < d) v[[k]] <- v[[k]] / sqrt(sum(v[[k]]^2))
	}
	if (scale) {
		v[[d]] <- v[[d]] / sqrt(sum(v[[d]]^2))
	} else {
		nrmv <- sqrt(sum(v[[d]]^2))
	}
	kronv <- Reduce(kronecker, v)
	e <- sqrt(sum((kronv - kronv.old)^2))
	if (e <= tol * max(nrmv, nrmv.old, 1)) break
}

if (scale) {
	return(v)
}
s <- c(rep.int(nrmv^(1/d), d - 1), nrmv^(1/d - 1))
for (k in 1:d) v[[k]] <- v[[k]] * s[k] 
v	
}

