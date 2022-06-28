
mcca.cor <- function(x, r, w = 1, ortho = c("block.score",
	"global.score", "canon.tnsr"), init = c("cca", "svd", "random"), 
	maxit = 1000, tol = 1e-6, sweep = c("cyclical", "random"), 
	control = list(), verbose = FALSE)
{

## Check arguments
test <- if (is.list(init)) {
	check.arguments(x, init, w)
} else check.arguments(x, NULL, w) 
eps <- 1e-14

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
	warning(paste("\nArgument 'r' set to", r,
		"to satisfy orthogonality constraints"))

## Initialization parameters
test <- (is.list(control) && !is.null(control$init))
if (is.character(init) && init == "cca") {
	init.args <- list(r = 1, k = NULL, w = w, 
		objective = "cor", cnstr = "block", center = FALSE, 
		search = ifelse(m <= 5, "exhaustive", "approximate"))
	if (test) {
		names. <- intersect(names(control), 
			c("k", "cnstr", "search"))
		init.args[names.] <- control[names.]
	}
}
if (is.character(init) && init == "svd") {
	init.args <- list(r = 1, objective = "cor", 
		cnstr = "block", center = FALSE)
	if (test) {
		names. <- intersect(names(control), c("cnstr"))
		init.args[names.] <- control[names.]
	}
}
if (is.character(init) && init == "random") 
	init.args <- list(r = 1, scale = "var")

## Prepare optimization if orthogonality 
## constraints pertain to canonical tensors
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
block.score <- array(0, c(n, m, r)) 
global.score <- matrix(0, n, r) 
input <- list(objective = "correlation", r = r, w = w, 
	init = init, ortho = ortho, maxit = maxit, tol = tol, 
	sweep = sweep) 
objective <- iters <- numeric(r)

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
vprev <- vector("list", m)
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
	if (verbose) cat("\n\nMCCA: Component",l,"\n")
	out <- mcca.single.cor(x, v0, w, sweep, maxit, tol, verbose)
	objective[l] <- out$objective
	block.score[,,l] <- out$y # canonical scores
	global.score[,l] <- rowMeans(block.score[,,l])	
	iters[l] <- out$iters	
	v[,l] <- out$v 
		
	## Prepare orthogonality constraints for next stage
	if (ortho  == "canon.tnsr" && l < r) {
		vprev <- lapply(d, function(len) vector("list", len))
		for (i in 1:m) {
			k <- ortho.mode[l, l+1, i]
			vprev[[i]][[k]] <- v[[i,l]][[k]]
		}
	}

	## Deflate data matrix
	if (l < r) { 	
		x <- switch(ortho, 
			block.score = deflate.x(x, score = block.score[,,l], 
				check.args = FALSE),
			global.score = deflate.x(x, score = global.score[,l], 
				check.args = FALSE),  
			canon.tnsr = deflate.x(x, v = vprev, check.args = FALSE))	
	}

	## Monitor objective value
	if (objective[l] <= eps) break 
	
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


