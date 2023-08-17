mcca.boot <- function(x, object, 
	target = c("v", "score.cov", "noise.cov", "score.cor", "noise.cor"),
	resample = c("data", "residuals"), nboot = 100, parallel = TRUE, 
	control = list())
{
## Determine bootstrap targets
lookup <- c("v", "score.cov", "noise.cov", "score.cor", "noise.cor")
target <- lookup[pmatch(target, lookup)]
stopifnot(length(target) > 0)

## Match other arguments
resample <- match.arg(resample)
.errorhandling <- match.arg(control[[".errorhandling"]],  
	c("stop", "remove", "pass"))
unbiased.score.cov <- as.logical(control[["unbiased.score.cov"]])
if (length(unbiased.score.cov) == 0) unbiased.score.cov <- TRUE
unbiased.noise.cov <- as.logical(control[["unbiased.noise.cov"]])
if (length(unbiased.noise.cov) == 0) unbiased.noise.cov <- TRUE

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
if (resample == "data") {
	xx <- x
} else { # resample == "residuals"
	resid <- calculate.residuals(x, object, matrix.out = TRUE)
	fitted <- lapply(1:m, function(i) as.numeric(x[[i]]) - resid[[i]])
	xx <- list(fitted = fitted, resid = resid, dimx = dimx)
}

## Perform bootstrap in parallel or sequentially 
if (parallel.flag) {
	boot.out <- foreach(b = 1:nboot, .errorhandling = .errorhandling) %dopar% {
			mcca.boot.single(xx, resample, target, call.args, unbiased.score.cov) 
	}
} else {
	boot.out <- vector("list", nboot)
	if (.errorhandling == "stop") {
		for (b in 1:nboot) 
			boot.out[[b]] <- mcca.boot.single(xx, resample, target, call.args, 
				unbiased.score.cov)
	} else {
		for (b in 1:nboot) 
			boot.out[[b]] <- tryCatch(
				mcca.boot.single(xx, resample, target, call.args, unbiased.score.cov), 
				error = function(e) e)
	}
	if (.errorhandling == "remove") {
		error.flag <- sapply(boot.out, inherits, what = "error")
		if (any(error.flag)) boot.out <- boot.out[!error.flag]
	}
}

## Estimate inference targets from original data
mcca.out <- list()
if ("v" %in% target) mcca.out$v <- object$v
score.target <- NULL
if ("score.cov" %in% target) score.target <- "cov"
if ("score.cor" %in% target) score.target <- c(score.target, "cor")
if (!is.null(score.target)) 
	mcca.out <- c(mcca.out, calculate.score.cov(object, score.target, 
		unbiased.score.cov, matrix.out = FALSE))
noise.target <- NULL
if ("noise.cov" %in% target) noise.target <- "cov"
if ("noise.cor" %in% target) noise.target <- c(noise.target, "cor")
if (!is.null(noise.target)) 
	mcca.out <- c(mcca.out, calculate.noise.cov(x, object, type = noise.target, 
		unbiased.noise.cov, matrix.out = FALSE))

## Realign the bootstrap estimates with original estimates
if ("v" %in% target) {
	for (b in 1:length(boot.out)) {
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


calculate.score.cov <- function(object, type = c("cov", "cor"), 
	unbiased = TRUE, matrix.out = TRUE)
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
	block.score.cov <- 
		array(apply(block.score, 3, cov), dim = c(m, m, r))
	out$block.score.cov <- block.score.cov
	if (unbiased) {
		adjust <- rowMeans(apply(block.score, 2, cov))
		global.score.cov <- (m * cov(global.score) - adjust) / (m-1)
	} else {
		global.score.cov <- cov(global.score)
	}
	out$global.score.cov <- global.score.cov
	if (!matrix.out) {
		idx <- which(lower.tri(matrix(0, m, m), TRUE))
		dim(out$block.score.cov) <- c(m^2, r)
		out$block.score.cov <- out$block.score.cov[idx,,drop = FALSE]
		idx <- which(lower.tri(matrix(0, r, r), TRUE))
		out$global.score.cov <- out$global.score.cov[idx]	
	}
}
if ("cor" %in% type) {
	if ("cov" %in% type) {
		out$block.score.cor <- 
			array(apply(block.score.cov, 3, cov2cor), dim = c(m, m, r))
		out$global.score.cor <- cov2cor(global.score.cov)
	} else {
		out$block.score.cor <- 
			array(apply(block.score, 3, cor), dim = c(m, m, r))
		if (unbiased) {
			adjust <- rowMeans(apply(block.score, 2, cov))
			out$global.score.cor <- cov2cor(m * cov(global.score) - adjust)
		} else {
			out$global.score.cor <- cor(global.score)
		} 
	}
	if (!matrix.out) {
		idx <- which(lower.tri(matrix(0, m, m)))
		dim(out$block.score.cor) <- c(m^2, r)
		out$block.score.cor <- out$block.score.cor[idx,,drop = FALSE]
		idx <- which(lower.tri(matrix(0, r, r)))
		out$global.score.cor <- out$global.score.cor[idx]	
	}
}
out
}



# # calculate.score.cov.no.bias <- function(score)
# {
# stopifnot(length(dim(score)) == 3)
# n <- dim(score)[1]
# m <- dim(score)[2]
# r <- dim(score)[3]

# global.score.cov <- matrix(0, r, r)
# mask <- lower.tri(matrix(0, m, m))
# for (l1 in 1:r) {
	# for (l2 in 1:l1) {
		# cov12 <- cov(score[,,l1], score[,,l2])
		# global.score.cov[l1,l2] <- mean(cov12[mask]) 
		# global.score.cov[l2,l1] <- global.score.cov[l1,l2]
	# }
# }
# global.score.cov
# }




calculate.residuals <- function(x, object, matrix.out = FALSE)
{
m <- length(x)
dimx <- lapply(x, dim)
n <- tail(dimx[[1]], 1)
r <- object$call.args$r	
objective.type <- object$call.args$objective.type
scale.cnstr <- object$call.args$scale
ortho.type <- object$call.args$ortho.type
test <- ((objective.type == "cov") && (scale.cnstr == "block") && 
	(ortho.type == "weight"))
resid <- vector("list", m)
for (i in 1:m) { 
	resid[[i]] <- matrix(x[[i]], ncol = n)
	vi <- matrix(unlist(tnsr.rk1.expand(object$v[i,])), ncol = r) 
	if (test) {
		score <- object$block.score[,i,] 
	} else {
		vi <- qr.Q(qr(vi))
		score <- crossprod(resid[[i]], vi)
	}
	resid[[i]] <- resid[[i]]  - tcrossprod(vi, score)
	if (!matrix.out) dim(resid[[i]]) <- dimx[[i]]	
}
resid	
}





calculate.noise.cov <- function(x, object, type = c("cov", "cor"), 
	unbiased = TRUE, matrix.out = TRUE)
{
type <- intersect(type, c("cov", "cor")) 
if (length(type) == 0) 
	stop("Please specify 'type' as 'cov', 'cor', or both these values.")
out <- list()
m <- length(x)
n <- tail(dim(x[[1]]), 1)
r <- NCOL(object$v)
if (unbiased) {
	noise.cov <- vector("list", m)
	for (i in 1:m) {
		vi <- matrix(unlist(tnsr.rk1.expand(object$v[i,])), ncol = r)
		vi <- sweep(vi, 2, sqrt(global.score.var), "*")
		noise.cov[[i]] <- cov(t(matrix(x[[i]], ncol = n))) - tcrossprod(vi)					}
	if ("cov" %in% type) 
		out$noise.cov <- noise.cov
	if ("cor" %in% type) 
		out$noise.cor <- lapply(noise.cov, cov2cor)
	rm(vi, noise.cov)
} else {
	resid <- calculate.residuals(x, object, matrix.out = TRUE)
	if ("cov" %in% type) 
		out$noise.cov <- lapply(resid, function(x) cov(t(x)))
	if ("cor" %in% type) 
		out$noise.cor <- lapply(resid, function(x) cor(t(x)))
	rm(resid)
}
if (!matrix.out) {
	for (i in 1:m) {
		if ("cov" %in% type) {
			idx <- which(lower.tri(matrix(0, pp[i], pp[i]), diag = TRUE))
			out$noise.cov[[i]] <- out$noise.cov[[i]][idx]	
		}
		if ("cor" %in% type) {
			idx <- which(lower.tri(matrix(0, pp[i], pp[i]), diag = FALSE))
			out$noise.cor[[i]] <- out$noise.cor[[i]][idx]	
		}
	}
}
out
}


mcca.boot.single <- function(x, resample, target, call.args, 
	unbiased.score.cov, unbiased.noise.cov) {
	
## Data dimensions
if (resample == "data") {
	m <- length(x)
	dimx <- lapply(x, dim)
} else { # resample == "residuals"
	dimx <- x$dimx
	m <- length(dimx) 
}
n <- tail(dimx[[1]], 1) 
out <- list()

## Resample data
idx <- sample(n, n, replace = TRUE)
if (resample == "data") {
	xstar <- x
	p <- lapply(dimx, function(dims) dims[-length(dims)])
	for (i in 1:m) {
		dim(xstar[[i]]) <- c(prod(p[[i]]), n)
		xstar[[i]] <- xstar[[i]][,idx]
		dim(xstar[[i]]) <- dimx[[i]]
	}
} else { # resample == "residuals"
	xstar <- vector("list", m)
	for (i in 1:m) {
		xstar[[i]] <- x$fitted[[i]] + x$resid[[i]][,idx]
		dim(xstar[[i]]) <- dimx[[i]]
	}
}

# Run mCCA on bootstrapped data
result <- if (call.args$objective.type == "cor") {
	mcca.cor(xstar, r = call.args$r, w = call.args$w, ortho = call.args$ortho.type,
		optim = call.args$optim, init = call.args$init.method, sweep = call.args$sweep,
		maxit = call.args$maxit, tol = call.args$tol, control = call.args$control)
} else {
	mcca.cov(xstar, r = call.args$r, w = call.args$w, scale = call.args$scale, 
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
	out <- c(out, calculate.score.cov(result, score.target, unbiased.score.cov,
		matrix.out = FALSE))

## Calculate covariance and/or correlation of residuals
noise.target <- NULL
if ("noise.cov" %in% target) noise.target <- "cov"
if ("noise.cor" %in% target) noise.target <- c(noise.target, "cor")
if (!is.null(noise.target)) 
	out <- c(out, calculate.noise.cov(xstar, result, noise.target, 
		unbiased.noise.cov, matrix.out = FALSE))

out
} 

