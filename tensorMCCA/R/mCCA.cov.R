mcca.cov <- function(x, r = 1, w = 1, norm = c("block", "global"), 
	ortho = c("score", "canon.tnsr"), init = c("cca", "svd", "random"), 
	maxit = 1000, tol = 1e-6, sweep = c("cyclical", "random"), 
	control = list(), verbose = FALSE)
{

## Check arguments
test <- if (is.list(init)) { 
	check.arguments(x, init, w)
} else check.arguments(x, NULL, w)
eps <- 1e-14
r0 <- r 
r <- as.integer(r)

## Data dimensions
m <- length(x) 
dimx <- lapply(x, dim) 
d <- sapply(dimx, length) - 1L 
p <- mapply(head, dimx, d, SIMPLIFY = FALSE) 
n <- tail(dimx[[1]], 1) 

## Objective weights
w <- if (length(w) == 1) {
	matrix(1 / (m^2), m, m)
} else {
	(w + t(w)) / (2 * sum(w)) 
}

## Match other arguments
ortho <- match.arg(ortho)
norm <- match.arg(norm)
sweep <- match.arg(sweep)
if (is.character(init)) init <- match.arg(init)

## Set initialization method
mcca.init <- if (is.character(init)) {
	switch(init, cca = mcca.init.cca, 
		svd = mcca.init.svd, random = mcca.init.random)
} else NULL

## Center data
if (ortho == "canon.tnsr" && norm == "block") {
	xbar <- mapply(rowMeans, x = x, dims = d, SIMPLIFY = FALSE)	
} else {	
	for(i in 1:m) {
	    xbar <- rowMeans(x[[i]], dims = d[i])
	    x[[i]] <- x[[i]] - as.vector(xbar)
	}
}

## Adjust number of canonical components 
r <- min(r, n - 1L)
r0 <- r
pp <- sapply(p, prod)
r <- if (ortho == "score" && norm == "block") { 
	min(n - 1, max(pp), r0)
} else if (ortho == "score" && norm == "global") {
	min(n - 1, sum(pp), r0)
} else if (ortho == "canon.tnsr" && norm == "block") {
	min(max(pp), r0)
} else {
	min(sum(pp), r0)
}
# if (verbose && r != r0)
	# warning(paste("\nArgument 'r' set to", r,
		# "to satisfy orthogonality constraints"))

## Prepare optimization if orthogonality 
## constraints are on canonical tensors
if (ortho == "canon.tnsr" && r > 1) {
	## Set tensor modes on which constraints apply
	ortho.mode <- set.ortho.mode(x, r, 
		cnstr = control$ortho$cnstr, 
		method = control$ortho$method)
}

## Initialization parameters
test <- (is.list(control) && !is.null(control$init))
if (identical(init, "cca")) {
	init.args <- list(k = NULL, w = w, 
		objective = "cov", norm = "block", center = FALSE, 
		search = ifelse(m <= 5, "exhaustive", "approximate"))
	if (test) {
		names. <- intersect(names(control$init), 
			c("k", "norm", "search"))
		init.args[names.] <- control$init[names.]
	}
} else if (identical(init, "svd")) {
	init.args <- list(objective = "cov", 
		norm = "block", center = FALSE)
	if (test) {
		names. <- intersect(names(control$init), c("norm"))
		init.args[names.] <- control$init[names.]
	}
} else if (identical(init, "random")) {
	init.args <- list(r = 1L, objective = "cov")
}

## Create output objects 
v <- vector("list", m * r) # canonical vectors
dim(v) <- c(m, r)
block.score <- array(dim = c(n, m, r)) # canonical scores
global.score <- matrix(nrow = n, ncol = r) 
objective <- iters <- numeric(r)
input <- list(objective = "covariance", r = r0, w = w, 
	init = init, ortho = ortho, maxit = maxit, tol = tol, 
	sweep = sweep) 

## Trivial case: constant datasets  
test <- sapply(x, function(a) all(abs(a) <= eps))
test <- outer(test, test, "&")
if (all(test | w == 0)) {
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
vzero <- vector("list", m)
for (i in 1:m) vzero[[i]] <- lapply(d, numeric) 
vortho <- NULL
# browser()
if (ortho == "canon.tnsr") x0 <- x

for (l in 1:r) {	
	
	
	## Deflate data matrix
	if (l > 1) { 	
		if (ortho == "score" && norm == "block") { 
			x <- deflate.x(x, score = block.score[,,l], 
				check.args = FALSE)
		} else if (ortho == "score" && norm == "global") { 
			x <- deflate.x(x, score = global.score[,l], 
				check.args = FALSE)
		} else if (ortho == "canon.tnsr" && norm == "block") {
			cnstr <- set.ortho.mat(v = v[, 1:(l-1)], 
				modes = t(ortho.mode[1:(l-1), l, ]))
			for(i in 1:m) 
				x[[i]] <- tnsr.mat.prod(
					x = x0[[i]] - as.vector(xbar[[i]]), 
					mat = cnstr$mat[[i]], modes = cnstr$modes[[i]])
		} 
	}
	
	## Initialize canonical vectors
	if (is.character(init)) {
		init.args$x <- x
		v0 <- do.call(mcca.init, init.args)
	} else if (is.list(init)) {
		v0 <- if (is.vector(init)) {
			init } else { init[, min(l, ncol(init))] }
		if (ortho == "canon.tnsr" && norm == "block")
			v0
	} else {
		v0 <- 
	} 
	
	## Run MCCA and store results
	if (verbose) cat("\n\nMCCA: Component", l, "\n")
	out <- if (norm == "block") {
		mcca.single.block.cov(x = x, v = v0, w = w, sweep = sweep, 
			maxit = maxit, tol = tol, verbose = verbose)
	} else {
		mcca.single.global.cov(x = x, v = v0, w = w, ortho = vortho, 
			sweep = sweep, maxit = maxit, tol = tol, verbose = verbose)
	}
	objective[l] <- out$objective
	block.score[,,l] <- out$y 
	global.score[,l] <- rowMeans(block.score[,,l]) 
	iters[l] <- out$iters
	v[,l] <- out$v

	## Monitor objective value
	if (objective[l] <= eps) break 

	## Prepare orthogonality constraints for next stage
	if (ortho == "canon.tnsr" && l < r) {
		if (norm == "block") {
			vprev <- lapply(d, function(len) vector("list", len))
			for (i in 1:m) 
				vprev[[i]][ortho.mode[[i]]] <- v[[i,l]][ortho.mode[[i]]]
		} else { # norm == "global"
			vprev <- v[, 1:l, drop = FALSE]
			for (i in 1:m) {
				idx <- setdiff(1:d[i], ortho.mode[[i]])
				for (ll in 1:l)
					vprev[[i,ll]][idx] <- lapply(p[[i]][idx], numeric)
			}
		}
	} 
		
} 

## Re-order results according to objective values if needed
o <- order(objective, decreasing = TRUE)
o <- o[objective[o] > eps]
if (!identical(o, 1:r)) {
	v <- v[,o]
	block.score <- block.score[,,o]
	global.score <- global.score[,o]
	objective <- objective[o]
	iters <- iters[o]
}

r <- length(o)
if (r == 1) {
	dim(block.score) <- c(n, m)
	dim(global.score) <- NULL
} 

list(v = v, block.score = block.score, global.score = global.score,
	objective = objective, iters = iters, input = input)
}


