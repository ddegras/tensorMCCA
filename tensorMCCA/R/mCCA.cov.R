mCCA.cov <- function(x, r, c = 1, init.type = c("svd", "ones", "random"), 
	init.value = NULL, cnstr = c("block", "global"), ortho = c("block.score", 
	"global.score", "canon.t.1", "canon.t.all"), balance = TRUE, maxit = 1000, 
	tol = 1e-6, sweep = c("cyclical", "random"), verbose = FALSE)
{

## Check arguments
test <- check.arguments(x, init.value)

## Data dimensions
m <- length(x) # number of datasets
dimx <- lapply(x, dim) # tensor dimensions for each dataset
dim.img <- lapply(dimx, function(x) x[-length(x)]) 
# tensor dimensions for each dataset, last mode (instances) omitted
ndim.img <- sapply(dim.img, length) 
# numbers of image dimensions for each dataset
n <- tail(dimx[[1]], 1) # numbers of instances per dataset

## Objective weights
stopifnot(length(c) == 1 || (is.matrix(c) && all(dim(c) == length(x))))
stopifnot(all(c >= 0) && any(c > 0))
if (length(c) == 1) {
	c <- matrix(1/m^2, m, m)
} else {
	c <- (c + t(c)) / (2 * sum(c))
}
# c <- c / n
csum <- colSums(c)

## Match arguments
init.type <- match.arg(init.type)
ortho <- match.arg(ortho)
cnstr <- match.arg(cnstr)
sweep <- match.arg(sweep)

## Optimization method
mCCA.single.cov <- if (cnstr == "block") {
	mCCA.single.block.cov } else {
	mCCA.single.global.cov }

## Data centering
for(i in 1:m) {
    mu <- rowMeans(x[[i]], dims = ndim.img[i])
    x[[i]] <- x[[i]] - as.vector(mu)
}
		
## Adjust number of canonical components as needed
r0 <- as.integer(r)
r <- switch(ortho, 
	block.score = min(r0, n-1), global.score = min(r0, n-1),
	canon.t.1 = min(r0, unlist(lapply(dim.img, "[[", 1))),
	canon.t.all = min(r0, unlist(dim.img)))
if (verbose && r != r0)
	warning(paste("Argument 'r' set to", r,
		"to satisfy orthogonality constraints"))
				
## Create output objects 
v <- vector("list", m) # canonical vectors
for (i in 1:m) 
	v[[i]] <- lapply(dim.img[[i]], function(dd) matrix(0, dd, r))
block.score <- array(dim = c(n, m, r)) 
global.score <- matrix(nrow = n, ncol = r) 
input <- list(obj = "cov", r = r, c = c, init.type = init.type, 
	init.value = init.value, ortho = ortho, balance = balance, 
	maxit = maxit, tol = tol, sweep = sweep) 


## MAIN LOOP
objective <- iters <- numeric(r)
for (k in 1:r) {	
	## Initialize canonical vectors
	if (!is.null(init.value)) {
		v0 <- init.value
		test <- NCOL(v0[[1]][[1]]) == r
		if (test && k > 1) {
			for (i in 1:m)
				for (l in 1:ndim.img[i])
					v0[[i]][[l]] <- v0[[i]][[l]][,k] 
		} 
	} else {
		v0 <- switch(init.type, 
			svd = init.v.svd(x, type = "norm", cnstr = cnstr),
			ones = init.v.ones(x, type = "norm", cnstr = cnstr),
			random = init.v.random(x, type = "norm", cnstr = cnstr))
	}
	
	## Run MCCA and store results
	if (verbose) cat("\n\nMCCA: Component",k,"\n")
	out <- mCCA.single.cov(x, v0, c, balance, sweep, maxit, tol, verbose)
	objective[k] <- out$objective
	block.score[,,k] <- out$y # canonical scores
	iters[k] <- out$iters

	## Deflate canonical vectors (necessary?)
	vk <- out$v 
	if (k > 1 && ortho %in% c("canon.t.1", "canon.t.all")) {
		vk <- deflate.v(out$v, v, ortho)	
		vk <- scale.v(vk, x, type = "norm", cnstr = cnstr)
		block.score[,,k] <- image.scores(x, vk) 
		objective[k] <- sum(c * crossprod(block.score[,,k])) / n
	}
	
	global.score[,k] <- block.score[,,k] %*% csum
	
	## Add new canonical vectors to output	
	for (i in 1:m) 
		for (dd in 1:ndim.img[i]) 
			v[[i]][[dd]][,k] <- vk[[i]][[dd]] 
	
	## Deflate data matrix
	if (k < r) 	
		x <- deflate.x(x, vk, block.score[,,k], ortho, check.args = FALSE)	
	
} 

## Re-order results according to objective values if needed
o <- order(objective, decreasing = TRUE)
if (!identical(o,1:r)) {
	for (i in 1:m)
		for (k in seq_along(v[[i]]))
			v[[i]][[k]] <- v[[i]][[k]][,o] 
	block.score <- block.score[,,o]
	global.score <- global.score[,o]
	objective <- objective[o]
	iters <- iters[o]
}


list(v = v, block.score = block.score, global.score = global.score,
	objective = objective, iters = iters, input = input)
}


