permutation.test <- function(x, obj, nperm = 100, parallel = TRUE)
{
## Prepare MCCA
test <- check.arguments(x)
eps <- 1e-14
m <- length(x)
n <- tail(dim(x[[1]]), 1)
dimx <- lapply(x, dim)
d <- sapply(dimx, length) - 1
fields <- c("init.val", "objective.type", "ortho.type", 
	"ortho.cnstr", "r", "scale")
for (name in fields) 
	assign(name, obj$call.args[[name]])
ortho <- ortho.type
if (ortho == "weight" && scale == "block") {
	ortho.mode <- ortho.cnstr
	if (!is.null(init.val)) {
		if (length(init.val) < m * r) 
			obj$call.args$init.val <- rep_len(init.val, m * r)
		dim(obj$call.args$init.val) <- c(m, r)
	}
}
v <- obj$v

## Center data
for(i in 1:m) {
    xbar <- rowMeans(x[[i]], dims = d[i])
    if (all(abs(xbar) < eps))
	    x[[i]] <- x[[i]] - as.vector(xbar)
}

## Build deflated datasets
# Note: deflation is not strictly needed if orthogonality constraints
# are on weights: it's already explicitly enforced during optimization.
# However, it helps finding better starting points. 
if (r > 1) {
	xdfl <- vector("list", m * (r - 1))
	dim(xdfl) <- c(m, r - 1)
	for (l in 2:r) {
		if (ortho == "weight" && scale == "block") {
			ortho.cnstr <- set.ortho.mat(v = v[,1:(l-1)], 
				modes = ortho.mode[, 1:(l-1), l])
			xdfl[,l-1] <- mapply(tnsr.mat.prod, x = x, 
				mat = ortho.cnstr$mat, modes = ortho.cnstr$modes, 
				SIMPLIFY = FALSE) 
			if (!is.null(init.val)) 
				obj$call.args$init.val[,l] <- tnsr.rk1.mat.prod(
					v = obj$call.args$init.val[,l], transpose.mat = FALSE,
					mat = ortho.cnstr$mat, modes = ortho.cnstr$modes)
		} else if (ortho == "weight" && scale == "global") {
			xdfl[,l-1] <- deflate.x(x, v = v[, 1:(l-1)], 
				check.args = FALSE)
		} else if (ortho == "score" && scale == "block") {
			xdfl[,l-1] <- deflate.x(x, 
				score = obj$block.score[,,1:(l-1)],
				scope = "block", check.args = FALSE)	
		} else if (ortho == "score" && scale == "global") {
			xdfl[,l-1] <- deflate.x(x, 
				score = obj$global.score[,1:(l-1)],
				scope = "global", check.args = FALSE)			
		}
	}	
}

## Scale initial canonical weights if provided
objective.cov <- (objective.type == "cov")
if (!is.null(init.val)) {
	for (l in 1:ncol(init.val)) 
		obj$call.args$init.val[,l] <- scale.v(init.val[,l], 
			type = if (objective.cov) "norm" else "var",
			scale = scale, x = if (objective.cov) { NULL 
			} else if (l == 1) { x } else { xdfl[,l-1] },
			check.args = FALSE)
}


## MAIN LOOP
if (parallel && require(foreach) && getDoParRegistered()) {
	rhoperm <- foreach(j = 1:nperm, .combine = rbind) %dopar% {
		set.seed(j)
		rhopermj <- numeric(r) # canonical correlations
		for (l in 1:r) {
			xx <- if (l == 1) x else xdfl[,l-1]
			for (i in 1:m) {
				p <- head(dim(xx[[i]]), d[i])
				dim(xx[[i]]) <- c(prod(p), n)
				xx[[i]] <- xx[[i]][, sample(n)]
				dim(xx[[i]]) <- c(p, n)
			}
			rhopermj[l] <- mcca.perm(xx, obj, l)
		}
		rhopermj
	}
	rownames(rhoperm) <- NULL
} else {
	rhoperm <- matrix(, nperm, r)
	for(j in 1:nperm) {
		set.seed(j)
		for (l in 1:r) {
			xx <- if (l == 1) x else xdfl[,l-1]
			for (i in 1:m) {
				p <- head(dim(xx[[i]]), d[i])
				dim(xx[[i]]) <- c(prod(p), n)
				xx[[i]] <- xx[[i]][, sample(n)]
				dim(xx[[i]]) <- c(p, n)
			}	
			rhoperm[j, l] <- mcca.perm(xx, obj, l)
		}
	}
}

## p-values
pval.uncorrected <- rowMeans(t(abs(rhoperm)) > abs(obj$objective))
pval.corrected <- Reduce(max, pval.uncorrected, accumulate = TRUE)

list(pval.corrected = pval.corrected, 
	pval.uncorrected = pval.uncorrected,
	objective = obj$objective, 
	objective.perm = rhoperm)
}






mcca.perm <- function(x, obj, l)
{
## Define variables with shorter names for easier handling
for (name in names(obj$call.args))
	assign(name, obj$call.args[[name]])
v <- obj$v
normvar <- switch(objective.type, cov = "norm", cor = "var")
if (ortho.type == "weight" && scale == "block") 
	ortho.cnstr <- set.ortho.mat(v = v[,1:(l-1)], 
		modes = ortho.cnstr[, 1:(l-1), l])

## Initialize canonical weights
if (!is.null(init.val)) {
	v0 <- if (is.matrix(init.val)) {
		init.val[, min(l, ncol(init.val))]
	} else init.val
	# if (ortho.type == "weight" && scale == "block" && l > 1) 
		# v0 <- tnsr.rk1.mat.prod(v0, mat = ortho.cnstr$mat, 
			# modes = ortho.cnstr$modes, transpose.mat = FALSE)
		# v0 <- scale.v(v0, type = "var", x = x, check.args = FALSE)
} else {
	mcca.init <- switch(init.method, 
		cca = mcca.init.cca, 
		svd = mcca.init.svd, 
		random = mcca.init.random)
	init.args$x <- x
	v0 <- do.call(mcca.init, init.args)
}
if (ortho.type == "score" && l > 1) {
	v0 <- tnsr.rk1.ortho(v0, ortho.cnstr[, 1:(l-1), drop = FALSE], 
		maxit = 100L, tol = 1e-6)
	v0 <- scale.v(v0, type = normvar, x = x, check.args = FALSE)
}

## Run MCCA
out <- if (objective.type == "cov" && optim == "bca" && 
	scale == "block") {
	mcca.cov.bca.block(x = x, v = v0, w = w, 
		ortho = if (ortho.type == "score" && l > 1) {
		ortho.cnstr[,1:(l-1)] } else NULL, 
		sweep = sweep, maxit = maxit, tol = tol, 
		verbose = FALSE)
} else if (objective.type == "cov" && optim == "bca" && 
	scale == "global") { 
	mcca.cov.bca.global(x = x, v = v0, w = w, 
		ortho = if (l == 1) { NULL 
		} else if (ortho.type == "weight") { v[, 1:(l-1)] 
		} else ortho.cnstr[,1:(l-1)], sweep = sweep, 
		maxit = maxit, tol = tol, verbose = FALSE)
} else if (objective.type == "cor" && optim == "bca") {
	mcca.cor.bca(x = x, v = v0, w = w, 
		ortho = if (ortho.type == "score" && l > 1) {
			ortho.cnstr[, 1:(l-1), drop = FALSE] } else NULL,
		sweep = sweep, maxit = maxit, tol = tol, verbose = FALSE)	
} else if (optim == "grad.scale") {
	mcca.gradient.scale(x = x, v = v0, w = w, scale = scale, 
		type = normvar, maxit = maxit, tol = tol, verbose = FALSE)
} else if (optim == "grad.rotate") {  
	mcca.gradient.rotate(x = x, v = v0, w = w, 
		maxit = maxit, tol = tol, verbose = FALSE)
}	

out$objective

}


