permutation.test <- function(x, object, nperm = 100, parallel = TRUE)
{
## Prepare MCCA
test <- check.arguments(x)
dimx <- lapply(x, dimfun)
d <- sapply(dimx, length) - 1L
m <- length(x)
n <- tail(dimx[[1]], 1)
p <- mapply(head, x = dimx, d, SIMPLIFY = FALSE)
p[d == 0] <- 1

fields <- c("init.val", "objective.type", "ortho.type", 
	"ortho.cnstr", "r", "scope")
for (name in fields) 
	assign(name, object$call.args[[name]])
ortho <- ortho.type
if (ortho == "weight" && scope == "block") 
	ortho.mode <- ortho.cnstr
v0 <- init.val
if (!is.null(v0) && is.vector(v0)) 
	dim(v0) <- c(m,1)
v <- object$v

## Center data
for(i in 1:m) {
    xbar <- if (d[i] == 0) { mean(x[[i]]) 
    	} else rowMeans(x[[i]], dims = d[i])
    if (any(abs(xbar) > 1e-16))
	    x[[i]] <- x[[i]] - as.vector(xbar)
}

## Build deflated datasets 
# Note: deflation is not strictly needed if orthogonality constraints
# are on weights: it's already explicitly enforced during optimization.
# However, it helps finding better starting points. 
if (r > 1) {
	xdfl <- vector("list", m * (r-1))
	dim(xdfl) <- c(m, r-1)
	for (l in 2:r) {
		if (ortho == "weight" && scope == "block") {
			ortho.cnstr <- set.ortho.mat(v = v[,1:(l-1)], 
				modes = ortho.mode[, 1:(l-1), l])
			xdfl[,l-1] <- mapply(tnsr.mat.prod, x = x, 
				mat = ortho.cnstr$mat, modes = ortho.cnstr$modes, 
				SIMPLIFY = FALSE) 
		} else if (ortho == "weight" && scope == "global") {
			xdfl[,l-1] <- deflate.x(x, scope = "global",
				ortho = v[,1:(l-1), drop=FALSE])
		} else if (ortho == "score" && scope == "block") {
			xdfl[,l-1] <- deflate.x(x, 
				ortho = ortho.cnstr[,1:(l-1), drop=FALSE])	
		} else if (ortho == "score" && scope == "global") {
			xdfl[,l-1] <- deflate.x(x, scope = "global",
				ortho = ortho.cnstr[,1:(l-1), drop=FALSE])			
		}
	}
}

## Enforce scaling and orthogonality constraints on 
## initial canonical weights if provided
if (!is.null(v0)) {
	scale.type <- switch(objective.type, cov = "norm", cor = "var")
	v0[,1] <- scale.v(v0[,1], scale.type, scope, x, FALSE)
	idx <- (1:ncol(v0))[-1]
	for (l in idx) {
		if (ortho == "weight" && scope == "block") {
			cnstr <- set.ortho.mat(v = v[,1:(l-1)], 
				modes = ortho.mode[, 1:(l-1), l])
			v0[,l] <- tnsr.rk1.mat.prod(v = v0[,l],
				mat = ortho.cnstr$mat, modes = ortho.cnstr$modes)
		} else {	
			cnstr <- switch(ortho, 
				weight = v[,1:(l-1),drop=FALSE],
				score = ortho.cnstr[,1:(l-1),drop=FALSE])
			v0[,l] <- tnsr.rk1.ortho(v0[,l], cnstr, 100L, 1e-6)
		}
		v0[,l] <- scale.v(v0[,l], scale.type, scope, xdfl[,l-1], FALSE)
	}
	object$call.args$init.val <- v0
}


## MAIN LOOP
if (parallel && require(foreach) && getDoParRegistered()) {
	rhoperm <- foreach(j = 1:nperm, .combine = rbind) %dopar% {
		set.seed(j)
		rhopermj <- numeric(r) # canonical correlations
		for (l in 1:r) {
			xx <- if (l == 1) x else xdfl[,l-1]
			for (i in 1:m) {
				p <- head(dimfun(xx[[i]]), d[i])
				dim(xx[[i]]) <- c(prod(p), n)
				xx[[i]] <- xx[[i]][, sample(n)]
				dim(xx[[i]]) <- c(p, n)
			}
			rhopermj[l] <- mcca.perm(xx, object, l)
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
				p <- head(dimfun(xx[[i]]), d[i])
				dim(xx[[i]]) <- c(prod(p), n)
				xx[[i]] <- xx[[i]][, sample(n)]
				dim(xx[[i]]) <- c(p, n)
			}	
			rhoperm[j, l] <- mcca.perm(xx, object, l)
		}
	}
}

## p-values
pval.uncorrected <- rowMeans(t(abs(rhoperm)) > abs(object$objective))
pval.corrected <- Reduce(max, pval.uncorrected, accumulate = TRUE)

list(pval.corrected = pval.corrected, 
	pval.uncorrected = pval.uncorrected,
	objective = object$objective, 
	objective.perm = rhoperm)
}






mcca.perm <- function(x, object, l)
{
## Define variables with shorter names for easier handling
for (name in names(object$call.args))
	assign(name, object$call.args[[name]])
v <- object$v
scale.type <- switch(objective.type, cov = "norm", cor = "var")
ortho.cnstr <- if (l == 1) {
	NULL
} else if (ortho.type == "weight" && scope == "block") {
	NULL
} else if (ortho.type == "weight" && scope == "global") {
	v[, 1:(l-1), drop = FALSE]
} else {
	ortho.cnstr[, 1:(l-1), drop = FALSE]
}

## Initialize canonical weights
if (!is.null(init.val) && ncol(init.val) >= l) {
	v0 <- init.val[,l] 
	# Already satisfies scaling & orthogonality constraints
} else {
	mcca.init <- switch(init.method, 
		cca = mcca.init.cca, 
		svd = mcca.init.svd, 
		random = mcca.init.random)
	init.args$x <- x
	v0 <- do.call(mcca.init, init.args)
	v0 <- tnsr.rk1.ortho(v0, ortho.cnstr, maxit = 100L, tol = 1e-6)
	## Scale initial solution
	v0 <- scale.v(v0, scale.type, scope, x, FALSE)
}


## Run TMCCA
out <- if (objective.type == "cov" && optim == "bca" && 
	scope == "block") {
	mcca.cov.bca.block(x = x, v = v0, w = w, ortho = ortho.cnstr, 
		sweep = sweep, maxit = maxit, tol = tol, verbose = FALSE)
} else if (objective.type == "cov" && optim == "bca" && 
	scope == "global") { 
	mcca.cov.bca.global(x = x, v = v0, w = w, ortho = ortho.cnstr, 
		sweep = sweep, maxit = maxit, tol = tol, verbose = FALSE)
} else if (objective.type == "cor" && optim == "bca") {
	mcca.cor.bca(x = x, v = v0, w = w, ortho = ortho.cnstr,
		sweep = sweep, maxit = maxit, tol = tol, verbose = FALSE)	
} else if (optim == "grad.scale") {
	mcca.gradient.scale(x = x, v = v0, w = w, scope = scope, 
		type = scale.type, maxit = maxit, tol = tol, verbose = FALSE)
} else if (optim == "grad.rotate") {  
	mcca.gradient.rotate(x = x, v = v0, w = w, 
		maxit = maxit, tol = tol, verbose = FALSE)
}	

out$objective

}


