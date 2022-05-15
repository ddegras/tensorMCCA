mCCA.cor <- function(x, r, c = 1, init.type = c("svd", "cca", "random"), 
	init.value = NULL, ortho = c("block.score", "global.score", "canon.tnsr"), 
	maxit = 1000, tol = 1e-6, sweep = c("cyclical", "random"), verbose = FALSE)
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
if (length(c) == 1) {
	c <- matrix(1/m^2, m, m)
} else {
	c <- (c + t(c)) / (2 * sum(c))
}

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
r <- switch(ortho, 
	block.score = min(n-1, sapply(p, prod), r0), 
	global.score = min(n-1, sum(sapply(p, prod)), r0),
	canon.tnsr = min(sapply(p, prod), r0))
if (verbose && r != r0)
	warning(paste("Argument 'r' set to", r,
		"to satisfy orthogonality constraints"))

## Set up sequence of constraints if orthogonality 
## constraints are on canonical tensors
if (ortho == "canon.tnsr") {
	ortho.mode <- array(0L, c(r, r, m))
	
	# , i=1:m, l=1:r, 
	
	
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



## Main loop
objective <- iters <- numeric(r)
for (l in 1:r) {	
	## Initialize canonical vectors
	if (!is.null(init.value) && NCOL(init.value) >= l) {
		v0 <- init.value[,l] 
	} else {
		v0 <- switch(init.type, 
			svd = init.mcca.svd(x, objective = "cor"),
			cca = init.mcca.cca(x, objective = "cor"),
			random = init.mcca.random(x, objective = "cor"))
	}
	
	## Run MCCA and store results
	if (verbose) cat("\n\nMCCA: Component",l,"\n")
	out <- mCCA.single.cor(x, v0, c, sweep, maxit, tol, verbose) 
	objective[l] <- out$objective
	block.score[,,l] <- out$y # canonical scores
	iters[l] <- out$iters
	
	## Deflate canonical vectors 
	vl <- out$v 
	if (l > 1 && ortho  == "canon.tnsr") {
		vl <- deflate.v(out$v, v, ortho)	
		vl <- scale.v(vl, x, type = "var")
		block.score[,,l] <- image.scores(x, vl) 
		objective[l] <- sum(c * crossprod(block.score[,,l])) / n
	}
	
	## Calculate global scores
	nrmv <- numeric(m)
	for (i in 1:m) 
		nrmv[i] <- prod(sapply(vl[[i]], function(x) sum(x^2)))	
	global.score[,l] <- rowSums(block.score[,,l]) / sum(nrmv)
	
	## Rescale block scores
	block.score[,,l] <- block.score[,,l] %*% diag(1/nrmv)
	
	## Add new canonical vectors to output	
	for (i in 1:m) 
		for (k in 1:d[i]) 
			v[i,l][[k]] <- vl[[i]][[k]] 
	
	## Deflate data matrix
	if (l < r) 	
		x <- deflate.x(x, vl, block.score[,,l], ortho, check.args = FALSE)	
	
} 

## Re-order results according to objective values if needed
o <- order(objective, decreasing = TRUE)
if (!identical(o,1:r)) {
	v <- v[,o]
	# for (i in 1:m)
		# for (l in seq_along(v[[i]]))
			# v[[i]][[l]] <- v[[i]][[l]][,o] 
	block.score <- block.score[,,o]
	global.score <- global.score[,o]
	objective <- objective[o]
	iters <- iters[o]
}

list(v = v, block.score = block.score, global.score = global.score,
	objective = objective, iters = iters, input = input)
}


