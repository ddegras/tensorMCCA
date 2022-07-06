mcca.cov <- function(x, r, w = 1, cnstr = c("block", "global"), 
	ortho = c("score", "canon.tnsr"), init = c("cca", "svd", "random"), 
	maxit = 1000, tol = 1e-6, sweep = c("cyclical", "random"), 
	control = list(), verbose = FALSE)
{

## Check arguments
test <- if (is.list(init)) { 
	check.arguments(x, init, w)
} else check.arguments(x, NULL, w)
eps <- 1e-14
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
cnstr <- match.arg(cnstr)
sweep <- match.arg(sweep)
if (is.character(init)) init <- match.arg(init)

## Optimization method
mcca.single.cov <- if (cnstr == "block") {
	mcca.single.block.cov } else {
	mcca.single.global.cov }
mcca.init <- if (is.character(init)) {
	switch(init, cca = mcca.init.cca, 
		svd = mcca.init.svd, random = mcca.init.random)
} else NULL

## Data centering
for(i in 1:m) {
    xbar <- rowMeans(x[[i]], dims = d[i])
    x[[i]] <- x[[i]] - as.vector(xbar)
}

## Prepare optimization if orthogonality 
## constraints are on canonical tensors
if (ortho == "canon.tnsr" && r > 1) {
	## Set tensor modes on which constraints apply
	if (length(control) > 0 && (!is.null(control$ortho))) {
		ortho.mode <- if (is.integer(control$ortho)) {
			as.list(rep_len(control$ortho, m))
		} else if (is.list(control$ortho)) {
			rep_len(control$ortho, m)
		} else if (identical(control$ortho, "all")) {
			mapply(seq.int, from = rep(1L, m), to = d, 
				SIMPLIFY = FALSE)
		} else {
			as.list(rep_len(1L, m))
		}
	} else {
		ortho.mode <- as.list(rep_len(1L, m))
	}
}
	
## Adjust number of canonical components as needed
## (Not strictly needed when calculating canonical components iteratively)
r0 <- r
if (cnstr == "block" && ortho == "score") { 
	pp <- sapply(p, prod)
	r <- min(n - 1, pp, r0)
} else if (cnstr == "global" && ortho == "score") {
	pp <- sapply(p, prod)
	r <- min(n - 1, sum(pp), r0)
} else if (cnstr == "block" && ortho == "canon.tnsr") {
	rr <- unlist(mapply("[", p, ortho.mode))
	r <- min(n - 1, rr, r0)
} else {
	rr <- sapply(1:m, function(i) min(p[ortho.mode[[i]]]))
	r <- min(n - 1, sum(rr), r0)
}
if (verbose && r != r0)
	warning(paste("\nArgument 'r' set to", r,
		"to satisfy orthogonality constraints"))

## Initialization parameters
test <- (is.list(control) && !is.null(control$init))
if (is.character(init) && init == "cca") {
	init.args <- list(r = 1, k = NULL, w = w, 
		objective = "cov", cnstr = "block", center = FALSE, 
		search = ifelse(m <= 5, "exhaustive", "approximate"))
	if (test) {
		names. <- intersect(names(control), 
			c("k", "cnstr", "search"))
		init.args[names.] <- control[names.]
	}
}
if (is.character(init) && init == "svd") {
	init.args <- list(r = 1, objective = "cov", 
		cnstr = "block", center = FALSE)
	if (test) {
		names. <- intersect(names(control), c("cnstr"))
		init.args[names.] <- control[names.]
	}
}
if (is.character(init) && init == "random") 
	init.args <- list(r = 1, scale = "norm")

## Create output objects 
v <- vector("list", m * r) # canonical vectors
dim(v) <- c(m, r)
block.score <- array(dim = c(n, m, r)) # canonical scores
global.score <- matrix(nrow = n, ncol = r) 
objective <- iters <- numeric(r)
input <- list(objective = "covariance", r = r, w = w, 
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
vprev <- NULL

for (l in 1:r) {	
	
	## Initialize canonical vectors
	if (is.character(init)) {
		init.args$x <- x
		v0 <- do.call(mcca.init, init.args)
	} else if (is.list(init) && is.vector(init)) {
		v0 <- init
	} else {
		v0 <- init[, min(l, NCOL(init))]
	} 
	
	## Run MCCA and store results
	if (verbose) cat("\n\nMCCA: Component", l, "\n")
	out <- mcca.single.cov(x, v0, w, sweep, maxit, tol, verbose)
	objective[l] <- out$objective
	block.score[,,l] <- out$y 
	global.score[,l] <- rowMeans(block.score[,,l]) 
	iters[l] <- out$iters
	v[,l] <- out$v

	## Prepare orthogonality constraints for next stage
	if (l < r && ortho == "canon.tnsr") {
		if (cnstr == "block") {
			vprev <- lapply(d, function(len) vector("list", len))
			for (i in 1:m) 
				vprev[[i]][ortho.mode[[i]]] <- v[[i,l]][ortho.mode[[i]]]
		} else {
			vprev <- v[, 1:l]
			for (i in 1:m) {
				idx <- setdiff(1:d[i], ortho.mode[[i]])
				for (ll in 1:l)
					vprev[[i,ll]][idx] <- lapply(p[[i]][idx], numeric)
			}
		}
	} 
	x	
	## Deflate data matrix
	if (l < r) { 	
		x <- if (cnstr == "block" && ortho == "score") { 
			deflate.x(x, score = block.score[,,l])
		} else if (cnstr == "global" && ortho == "score") { 
			deflate.x(x, score = global.score[,l])
		 else { 
			deflate.x(x, v = vprev, scope = cnstr)	
	}
	
	## Monitor objective value
	if (objective[l] <= eps) break 
		
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

if (r == 1) {
	dim(block.score) <- c(n, m)
	dim(global.score) <- NULL
} 

list(v = v, block.score = block.score, global.score = global.score,
	objective = objective, iters = iters, input = input)
}


