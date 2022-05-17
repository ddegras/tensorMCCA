
mCCA.cor <- function(x, r, c = 1, ortho = c("block.score", "global.score", "canon.tnsr"), 
	init.type = c("svd", "cca", "random"), init.value = NULL, maxit = 1000, tol = 1e-6,
	sweep = c("cyclical", "random"), verbose = FALSE)
{

## Check arguments
test <- check.arguments(x, init.value)

## Data dimensions
m <- length(x) # number of datasets
dimx <- lapply(x, dim) # tensor dimensions for each dataset
p <- lapply(dimx, function(x) x[-length(x)]) 
# tensor dimensions for each dataset, last mode (instances) omitted
d <- sapply(p, length) 
# numbers of image dimensions for each dataset
n <- tail(dimx[[1]], 1) # numbers of instances per dataset

## Objective weights
stopifnot(length(c) == 1 || (is.matrix(c) && all(dim(c) == length(x))))
stopifnot(all(c >= 0) && any(c > 0))
c <- if (length(c) == 1) { 	matrix(1/m^2, m, m) 
	} else { (c + t(c)) / (2 * sum(c)) }

## Match arguments
init.type <- match.arg(init.type)
ortho <- match.arg(ortho)
sweep <- match.arg(sweep)

## Data centering
for(i in 1:m) {
    mu <- rowMeans(x[[i]], dims = d[i])
    x[[i]] <- x[[i]] - as.vector(mu)
}
		
## Adjust number of canonical components as needed
r0 <- as.integer(r)
pp <- sapply(p, prod)
r <- switch(ortho, 
	block.score = min(n-1, max(pp), r0), 
	global.score = min(n-1, sum(pp), r0),
	canon.tnsr = min(max(pp), r0))
if (verbose && r != r0)
	warning(paste("Argument 'r' set to", r,
		"to satisfy orthogonality constraints"))

## Prepare optimization if orthogonality 
## constraints are on canonical tensors
if (ortho == "canon.tnsr" && r > 1) {
	## Set tensor modes on which constraints apply
	ortho.mode <- array(dim = c(r, r, m))	
	for (i in 1:m) {
		mat <- matrix(0L, r, r)
		mod <- unlist(lapply((r-1):1, 
			function(len) rep_len(1:d[i], len)))
		if (sweep == "random") mod <- sample(mod)
		mat[lower.tri(mat)] <- mod
		ortho.mode[,,i] <- mat + t(mat)	
	}
}

## Create output objects 
v <- vector("list", m * r) # canonical vectors
dim(v) <- c(m, r)
block.score <- array(dim = c(n, m, r)) 
global.score <- matrix(nrow = n, ncol = r) 
input <- list(obj = "cor", r = r, c = c, init.type = init.type, 
	init.value = init.value, ortho = ortho, maxit = maxit, 
	tol = tol, sweep = sweep) 

## Trivial case: constant datasets  
test <- sapply(x, function(a) all(abs(a) <= 1e-15))
test <- outer(test, test, "&")
if (all(test | c == 0)) {
	r <- 1
	input$r <- r
	v <- vector("list", m)
	for (i in 1:m) 
		v[[i]] <- lapply(p[[i]], function(len) matrix(0, len, r))
	return(list(v = v, block.score = array(0, c(n, m, r)), 
		global.score = matrix(0, n, r), objective = 0, 
		iters = 0, input = input))
}


## MAIN LOOP
objective <- iters <- numeric(r)
vprev <- NULL
for (l in 1:r) {	
	## Initialize canonical vectors
	v0 <- if (!is.null(init.value) && NCOL(init.value) == 1) {
		init.value 
	} else if (!is.null(init.value) && NCOL(init.value) >= l) {
		init.value[, l]	
	} else {
		switch(init.type, 
			svd = init.mcca.svd(x, objective = "cor", center = FALSE),
			cca = init.mcca.cca(x, objective = "cor", center = FALSE),
			random = init.mcca.random(x, objective = "cor", center = FALSE))
	}
	
	## Run MCCA and store results
	if (verbose) cat("\n\nMCCA: Component",l,"\n")
	out <- mCCA.single.cor(x, v0, c, sweep, maxit, tol, verbose)
	objective[l] <- out$objective
	block.score[,,l] <- out$y # canonical scores
	iters[l] <- out$iters	
	v[,l] <- out$v 
	
	## Calculate global scores
	global.score[,l] <- rowMeans(block.score[,,l])	
	
	## Prepare orthogonality constraints for next stage
	if (ortho  == "canon.tnsr" && l < r) {
		for (i in 1:m) {
			k <- ortho.mode[l, l+1, i]
			vprev[[i]][[k]] <- v[[i,l]][[k]]
		}
	}

	## Deflate data matrix
	if (l < r) 	
		x <- deflate.x(x, vprev, block.score[,,l], ortho, check.args = FALSE)	

	## Monitor objective value
	if (objective[l] <= 1e-14) break 
	
} 

## Re-order results according to objective values if needed
objective <- objective[objective > 0]
o <- order(objective, decreasing = TRUE)
if (!identical(o,1:r)) {
	v <- v[,o]
	block.score <- block.score[,,o]
	global.score <- global.score[,o]
	objective <- objective[o]
	iters <- iters[o]
}

list(v = v, block.score = block.score, global.score = global.score,
	objective = objective, iters = iters, input = input)
}


