
mcca.cor <- function(x, r, w = 1, ortho = c("score", "weight"), 
	optim = c("bca", "grad.scale"), init = c("cca", "svd", "random"), 
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
sweep <- match.arg(sweep)
optim <- match.arg(optim)
if (is.character(init)) init <- match.arg(init)

## Set initialization method
if (is.character(init)) 
	mcca.init <- switch(init, cca = mcca.init.cca, 
		svd = mcca.init.svd, random = mcca.init.random)

## Center data
for(i in 1:m) {
    mu <- rowMeans(x[[i]], dims = d[i])
    x[[i]] <- x[[i]] - as.vector(mu)
}
if (ortho == "weight") x0 <- x

## Adjust number of canonical components as needed
r0 <- as.integer(r)
pp <- sapply(p, prod)
r <- switch(ortho, score = min(n - 1, max(pp), r0),
	weight = min(max(pp), r0)) 

## Set orthogonality constraints
if (r > 1) {
	if (ortho == "weight") {
		## Set tensor modes on which constraints apply
		ortho.mode <- set.ortho.mode(x, r, 
			cnstr = control$ortho$cnstr, 
			method = control$ortho$method)
	} else {
		ortho.cnstr <- vector("list", m * r)
		dim(ortho.cnstr) <- c(m, r)
	}
}

## Set up initialization parameters
test <- (is.list(control) && !is.null(control$init))
if (identical(init, "cca")) {
	init.args <- list(k = NULL, w = w, 
		objective = "cor", scale = "block", center = FALSE, 
		search = ifelse(m <= 5, "exhaustive", "approximate"))
	if (test) {
		names. <- intersect(names(control), 
			c("k", "search"))
		init.args[names.] <- control[names.]
	}
} else if (identical(init, "svd")) {
	init.args <- list(objective = "cor", center = FALSE)
} else if (identical(init, "random")) { 
	init.args <- list(objective = "cor")
}

## Create output objects 
v <- vector("list", m * r) # canonical vectors
dim(v) <- c(m, r)
block.score <- array(0, c(n, m, r)) 
global.score <- matrix(0, n, r) 
objective <- iters <- numeric(r)
input <- list(objective = "cor", r = r0, w = w, 
	ortho = ortho, init = init, maxit = maxit, 
	tol = tol, sweep = sweep, control = control)


## MAIN LOOP

for (l in 1:r) {	
	
	## Prepare orthogonality constraints after first round
	if (l > 1) {
		if (ortho == "weight") {
			## Deflate data
			ortho.cnstr <- set.ortho.mat(v = v[, 1:(l-1)], 
				modes = ortho.mode[, 1:(l-1), l])
			x <- deflate.x(x = x0, v = v[,1:(l-1)],  
				ortho.mode = ortho.mode[,1:(l-1),l], 
				check.args = FALSE)
		} else { # ortho == "score"
			## Calculate tensors that form orthogonality constraints		
			for (i in 1:m)
				ortho.cnstr[[i,l-1]] <- tnsr.vec.prod(x[[i]], 
					block.score[,i,l-1], d[i]+1)
		}
	}
		
	## Initialize canonical vectors
	if (is.character(init)) {
		init.args$x <- x
		v0 <- do.call(mcca.init, init.args)
	} else if (is.list(init)) {
		v0 <- if (is.vector(init)) {
			init } else { init[, min(l, ncol(init))] }
		if (l > 1) {
			v0 <- if (ortho == "weight") {
				tnsr.rk1.mat.prod(v0, mat = ortho.cnstr$mat, 
					modes = ortho.cnstr$modes, transpose.mat = FALSE)
			} else {
				tnsr.rk1.ortho(v0, ortho.cnstr[,1:(l-1)], 
					maxit = 100L, tol = 1e-6)
			}
		}
		v0 <- scale.v(v0, type = "var", x = x, check.args = FALSE)
	}
	
	## Run MCCA and store results
	if (verbose) cat("\n\nMCCA: Component",l,"\n")
	out <- mcca.single.cor(x = x, v = v0, w = w, 
		ortho = if (ortho == "weight") ortho.cnstr[,1:(l-1)] else NULL,
		sweep = sweep, maxit = maxit, tol = tol, verbose = verbose)
	v[,l] <- out$v
	objective[l] <- out$objective
	block.score[,,l] <- out$y 
	global.score[,l] <- rowMeans(block.score[,,l])	
	iters[l] <- out$iters	
	
	## Transform back canonical tensors to original search space 
	## if needed	
	if (l > 1 && ortho == "weight") {
		v[,l] <- tnsr.rk1.mat.prod(v = v[,l], 
			mat = ortho.cnstr$mat, modes = ortho.cnstr$modes, 
			transpose.mat = TRUE)
	}

	## Monitor objective value
	if (objective[l] <= eps) break 
	
} 

## Re-order results according to objective values if needed
objective <- objective[objective > 0]
o <- order(objective, decreasing = TRUE)
if (!identical(o, 1:r)) {
	v <- v[,o]
	block.score <- block.score[,,o]
	global.score <- global.score[,o]
	objective <- objective[o]
	iters <- iters[o]
}

## Clean up scores
block.score[abs(block.score) < eps] <- 0
global.score[abs(global.score) < eps] <- 0

list(v = v, block.score = block.score, global.score = global.score,
	objective = objective, iters = iters, input = input)
}


