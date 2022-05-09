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
	block.score = min(r0, n-1), global.score = min(r0, n-1),
	canon.t.1 = min(r0, unlist(lapply(p, "[[", 1))),
	canon.t.all = min(r0, unlist(p)))
if (verbose && r != r0)
	warning(paste("Argument 'r' set to", r,
		"to satisfy orthogonality constraints"))
				
## Create output objects 
v <- vector("list", m * r) # canonical vectors
if (r > 1) dim(v) <- c(m, r)
# for (i in 1:m) 
	# v[[i]] <- lapply(p[[i]], function(dd) matrix(0, dd, r))
block.score <- array(dim = c(n, m, r)) 
global.score <- matrix(nrow = n, ncol = r) 
input <- list(obj = "cor", r = r, c = c, init.type = init.type, 
	init.value = init.value, ortho = ortho, maxit = maxit, 
	tol = tol, sweep = sweep) 




## Main loop
objective <- iters <- numeric(r)
for (l in 1:r) {	
	## Initialize canonical vectors
	if (!is.null(init.value)) {
		v0 <- init.value
		test <- NCOL(v0[[1]][[1]]) == r
		if (test && l > 1) {
			for (i in 1:m)
				for (k in 1:d[i])
					v0[[i]][[k]] <- v0[[i]][[k]][, l, drop = FALSE] 
		} 
	} else {
		v0 <- switch(init.type, 
			svd = init.mcca.svd(x, objective = "cor"),
			cca = init.mcca.cca(x, objective = "cor"),
			random = init.mcca.random(x, objective = "cor"))
	}
	
	## Run MCCA and store results
	if (verbose) cat("\n\nMCCA: Component",l,"\n")
	# Insert here code to identify orthogonality constraints 
	# w.r.t. previous canonical vectors ONLY needed if ortho == "canon.tnsr"
	if (ortho == "canon.tnsr") {
		# ... 
	} else vprev <- NULL
	out <- mCCA.single.cor(x, v0, c, sweep, maxit, tol, verbose) # ortho = vprev
	objective[l] <- out$objective
	block.score[,,l] <- out$y # canonical scores
	iters[l] <- out$iters
	
	## Deflate canonical vectors (necessary?)
	# Not anymore: this will be done (if needed) inside the optimization
	# vl <- out$v 
	# if (l > 1 && ortho %in% c("canon.t.1", "canon.t.all")) {
		# vl <- deflate.v(out$v, v, ortho)	
		# vl <- scale.v(vl, x, type = "var")
		# block.score[,,l] <- image.scores(x, vl) 
		# objective[l] <- sum(c * crossprod(block.score[,,l])) / n
	# }
	
	## Calculate global scores
	nrmv <- numeric(m)
	for (i in 1:m) 
		nrmv[i] <- prod(sapply(vl[[i]], function(x) sum(x^2)))	
	global.score[,l] <- rowSums(blocl.score[,,l]) / sum(nrmv)
	
	## Rescale block scores
	block.score[,,l] <- block.score[,,l] %*% diag(1/nrmv)
	
	## Add new canonical vectors to output	
	for (i in 1:m) 
		for (k in 1:d[i]) 
			v[i,l][[k]] <- vl[[i]][[k]] 
	
	## Deflate data matrix
	if (l < r && ortho %in% c("block.score", "global.score")) 	
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


