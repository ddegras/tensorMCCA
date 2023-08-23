
hosvd <- function(x, r = NULL)
{
stopifnot(is.numeric(x))
if (length(dim(x)) <= 1) {
	p <- length(x) 
	r <- 1L
} else {
	p <- dim(x)
}
d <- length(p) 
if (is.null(r)) {
	r <- p 
} else {
	stopifnot(all(r > 0))
	r <- rep_len(as.integer(r), d)
	r <- pmin(r, p)
}

## Trivial cases 
if (d == 1L) {
	nrm <- sqrt(sum(x^2))
	if (nrm == 0) {
		return(list(factors = list(rep.int(1/sqrt(p), p)), 
			core = 0))
	} else {
		return(list(factors = list(x/nrm), core = nrm))
	}
}
if (d == 2L) {
	RSpectra.flag <- all(p >= 3)
	svdx <- if (RSpectra.flag) {
		tryCatch(svds(x, max(r), r[1], r[2]), 
			error = function(e) NULL)	
	} else NULL
	if (is.null(svdx)) 
		svdx <- svd(x, r[1], r[2])
	return(list(factors = list(svdx$u, svdx$v), 
		core = diag(svdx$d, r[1], r[2])))
}

## Outputs
core <- x
u <- vector("list", d)

for (k in 1:d) {
	## Mode-k flattening
	if (k > 1) {
		perm <- 1:d
		perm[c(1,k)] <- perm[c(k,1)]
		core <- aperm(core, perm)
	}
	dimc <- dim(core)
	dim(core) <- c(p[k], prod(dimc[-1]))
	
	## SVD
	RSpectra.flag <- all(dim(core) > max(2, r[k]))
	svdk <- if (RSpectra.flag) {
		tryCatch(svds(core, r[k]), error = function(e) NULL)
	} else NULL
	if (is.null(svdk))
		svdk <- svd(core, nu = r[k], nv = r[k]) }
	u[[k]] <- svdk$u
	
	## Reshape and permute core
	core <- svdk$d[1:(r[k])] * t(svdk$v)
	dim(core) <- c(r[k], dimc[-1])
	if (k > 1) core <- aperm(core, perm)	
}	

list(factors = u, core = core)
	
}


