mcca.boot <- function(x, object, method = c("raw", "residual"), nboot = 100, parallel = TRUE, 
	target = c("v", "score.cov", "noise.cov", "score.cor", "noise.cor"))
{
## Determine bootstrap targets
lookup <- c("v", "score.cov", "noise.cov", "score.cor", "noise.cor")
target <- lookup[pmatch(target, lookup)]
stopifnot(length(target) > 0)

## Match other arguments
method <- match.arg(method)

## Determine if bootstrap should be computed in parallel 
parallel.flag <- FALSE
if (parallel) {
	parallel.flag <- require(foreach) && getDoParRegistered()
	if (!parallel.flag) 
		warning("Executing the function sequentially. For parallel execution,", 
			"please load the package 'foreach' and register a backend, e.g.,",
			"'doParallel', before running 'mcca.boot'.")
}

## Data dimensions
m <- length(x)
dimx <- lapply(x, dim)
p <- lapply(dimx, function(dims) dims[-length(dims)])
n <- tail(dimx[[1]], 1) 
d <- sapply(p, length)

## Type of maximization in fitted object
r <- object$call.args$r
v <- object$v
if (is.null(dim(v))) dim(v) <- c(m, r)
call.args <- object$call.args
objective.type <- call.args$objective.type

## Prepare resampling in case of residual-based bootstrap
if (method == "residual") {
	resid <- calculate.residuals(x, object, matrix.out = TRUE)
	fitted <- lapply(1:m, function(i) as.numeric(x[[i]]) - resid[[i]])
}

## Perform bootstrap in parallel or sequentially 
if (parallel.flag) {
	boot.out <- if (method == "raw") {
		foreach(b = 1:nboot) %dopar% mcca.boot.single(x, method, target, call.args)
	} else foreach(b = 1:nboot) %dopar% mcca.boot.single(list(fitted, resid), method, 
		target, call.args)
} else {
	boot.out <- vector("list", nboot)
	if (method == "raw") {
		for (b in 1:nboot) boot.out[[b]] <- mcca.boot.single(x, method, target, call.args)
	} else {
		for (b in 1:nboot) boot.out[[b]] <- mcca.boot.single(list(fitted, resid), method, 
			target, call.args)
	}
}

## Estimate inference targets from original data
mcca.out <- list()
score.target <- NULL
if ("score.cov" %in% target) score.target <- "cov"
if ("score.cor" %in% target) score.target <- c(score.target, "cor")
if (!is.null(score.target)) 
	mcca.out <- c(mcca.out, calculate.score.cov(object, score.target, 
		matrix.out = FALSE))
noise.target <- NULL
if ("noise.cov" %in% target) noise.target <- "cov"
if ("noise.cov" %in% target) noise.target <- c(noise.target, "cor")
if (!is.null(noise.target)) 
	mcca.out <- c(mcca.out, calculate.noise.cov(x, object, type = score.target, 
		matrix.out = FALSE))

## Realign the bootstrap estimates with original estimates
if ("v" %in% target) {
	mcca.out$v <- object$v
	for (b in 1:nboot) {
		idx <- which(colMeans(cosine.dist(boot.out[[b]]$v, mcca.out$v)) > 1)
		for (l in idx) { 
			for (i in 1:m) {
				boot.out[[b]]$v[[i,l]][[1]] <- - boot.out[[b]]$v[[i,l]][[1]]
			}
		}	
	}
}

list(bootstrap = boot.out, original = mcca.out)
}


###################
# Helper functions
###################


calculate.score.cov <- function(object, type = c("cov", "cor"), matrix.out = TRUE)
{
type <- intersect(type, c("cov", "cor")) 
if (length(type) == 0) 
	stop("Please specify 'type' as 'cov', 'cor', or both these values.")
out <- list()
block.score <- object$block.score
dims <- dim(block.score)
n <- dims[1]
m <- dims[2]
r <- dims[3]
global.score <- object$global.score
if ("cov" %in% type) {	
	out$score.cov.block <- array(apply(block.score, 3, cov), dim = c(m, m, r))
	out$score.cov.global <- cov(global.score)
	if (!matrix.out) {
		idx <- which(lower.tri(matrix(0, m, m), TRUE))
		dim(out$score.cov.block) <- c(m^2, r)
		out$score.cov.block <- out$score.cov.block[idx,,drop = FALSE]
		idx <- which(lower.tri(matrix(0, r, r), TRUE))
		out$score.cov.global <- out$score.cov.global[idx]	
	}
}
if ("cor" %in% type) {
	 out$score.cor.block <- array(apply(block.score, 3, cor), dim = c(m, m, r))
	out$score.cor.global <- cor(global.score)
	if (!matrix.out) {
		idx <- which(lower.tri(matrix(0, m, m)))
		dim(out$score.cor.block) <- c(m^2, r)
		out$score.cor.block <- out$score.cor.block[idx,,drop = FALSE]
		idx <- which(lower.tri(matrix(0, r, r)))
		out$score.cor.global <- out$score.cor.global[idx]	
	}
}
out
}


calculate.residuals <- function(x, object, matrix.out = FALSE)
{
m <- length(x)
dimx <- lapply(x, dim)
p <- lapply(dimx, function(dims) dims[-length(dims)])
pp <- sapply(p, prod)
n <- tail(dimx[[1]], 1)
r <- object$call.args$r	
objective.type <- object$call.args$objective.type
scale.cnstr <- object$call.args$scale
ortho.type <- object$call.args$ortho.type
test <- (objective.type == "cov" && scale.cnstr == "block" && ortho.type == "weight") 
resid <- vector("list", m)
for (i in 1:m) { 
	resid[[i]] <- matrix(x[[i]], pp[i], n)
	vi <- matrix(unlist(tnsr.rk1.expand(object$v[i,])), pp[i], r) 
	score <- if (test) {
		object$block.score[,i,] 
	} else {
		crossprod(resid[[i]], qr.Q(qr(vi)))
	}
	resid[[i]] <- resid[[i]]  - tcrossprod(vi, score)
	if (!matrix.out) dim(resid[[i]]) <- dimx[[i]]	
}
resid	
}


calculate.noise.cov <- function(x, object, type = c("cov", "cor"), matrix.out = TRUE)
{
type <- intersect(type, c("cov", "cor")) 
if (length(type) == 0) 
	stop("Please specify 'type' as 'cov', 'cor', or both these values.")
out <- list()
if (!matrix.out) {
	m <- length(x)
	dimx <- lapply(x, dim)
	pp <- sapply(dimx, function(dims) prod(dims[-length(dims)]))
}
resid <- calculate.residuals(x, object, matrix.out = TRUE)
if ("cov" %in% type) {
	out$noise.cov <- lapply(resid, cov)	
	if (!matrix.out) {
		for (i in 1:m) {
			idx <- which(lower.tri(matrix(0, pp[i], pp[i]), diag = TRUE))
			out$noise.cov[[i]] <- out$noise.cov[[i]][idx]	
		}
	}		
}
if ("cor" %in% type) {
	out$noise.cor <- lapply(resid, cor)	
	if (!matrix.out) {
		for (i in 1:m) {
			idx <- which(lower.tri(matrix(0, pp[i], pp[i]), diag = FALSE))
			out$noise.cor[[i]] <- out$noise.cor[[i]][idx]	
		}
	}		
}
out
}


mcca.boot.single <- function(x, method, target, call.args) {
	
## Data dimensions
if (method == "raw") {
	m <- length(x)
	dimx <- lapply(x, dim)
} else { # method == "residual"
	m <- length(x[[1]]) 
	dimx <- lapply(x[[1]], dim)	
}
n <- dimx[[1]][length(dimx[[1]])] 
p <- lapply(dimx, function(dims) dims[-length(dims)])
d <- sapply(p, length)
r <- call.args$r
out <- list()

## Resample data
idx <- sample(n, n, replace = TRUE)
if (method == "raw") {
	xstar <- x
	for (i in 1:m) {
		dim(xstar[[i]]) <- c(prod(p[[i]]), n)
		xstar[[i]] <- xstar[[i]][,idx]
		dim(xstar[[i]]) <- dimx[[i]]
	}
} else { # method == "residual"
	xstar <- x$fitted
	for (i in 1:m) {
		xstar[[i]] <- xstar[[i]] + x$resid[[i]][,idx]
		dim(xstar[[i]]) <- dimx[[i]]
	}
}

# Run mCCA on bootstrapped data
objective.type <- call.args$objective.type
init.args <- call.args$init.args
result <- if (objective.type == "cor") {
	mcca.cor(xstar, r = call.args$r, w = init.args$w, ortho = call.args$ortho.type,
		optim = call.args$optim, init = call.args$init.method, sweep = call.args$sweep,
		maxit = call.args$maxit, tol = call.args$tol, control = call.args$control)
} else {
	mcca.cov(xstar, r = call.args$r, w = init.args$w, scale = init.args$scale, 
		ortho = call.args$ortho.type, optim = call.args$optim, init = call.args$init.method, 
		sweep = call.args$sweep, maxit = call.args$maxit, tol = call.args$tol, 
		control = call.args$control)
}

## Add canonical weights to output
if ("v" %in% target) out$v <- result$v

## Calculate covariance and/or correlation of scores
score.target <- NULL
if ("score.cov" %in% target) score.target <- "cov"
if ("score.cor" %in% target) score.target <- c(score.target, "cor")
if (!is.null(score.target)) 
	out <- c(out, calculate.score.cov(result, score.target, 
		matrix.out = FALSE))

## Calculate covariance and/or correlation of residuals
noise.target <- NULL
if ("noise.cov" %in% target) noise.target <- "cov"
if ("noise.cov" %in% target) noise.target <- c(noise.target, "cor")
if (!is.null(noise.target)) 
	out <- c(out, calculate.noise.cov(x, result, type = score.target, 
		matrix.out = FALSE))

out
} 

