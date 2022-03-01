permutation.test <- function(x, obj, nperm = 100, parallel = TRUE, ncores = NULL)
{
parallel.flag <- parallel && require(parallel)

f <- function(idx, x, obj) {
m <- length(x)
dimx <- lapply(x, dim)
n <- tail(dimx[[1]], 1)
r <- length(obj$objective)
for (i in 1:m) {
	d <- length(dimx[[i]])
	dim(x[[i]]) <- c(prod(head(dimx[[i]],d-1)), n)
	x[[i]] <- x[[i]][,sample(n)]
	dim(x[[i]]) <- dimx[[i]]
}
mCCA.cor(x, r = r, c = obj$c, maxit = obj$maxit, tol = obj$tol, 
	init.type = obj$init.type, init.value = obj$init.value, 
	ortho = obj$ortho, verbose = FALSE)$objective
}

if (is.null(ncores)) 
	ncores <- getOption("mc.cores", 2L)
r <- length(obj$objective)

perm.objective <- if (parallel.flag) {
	mclapply(1:nperm, f, x = x, obj = obj, mc.cores = ncores)
	} else { lapply(1:nperm, f, x = x, obj = obj) }
perm.objective <- unlist(perm.objective)
dim(perm.objective) <- c(r, nperm)
pvals <- rowMeans(perm.objective >= obj$objective)
list(pvals = pvals, perm.objective = t(perm.objective))
}