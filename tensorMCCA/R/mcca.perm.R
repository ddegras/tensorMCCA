mcca.perm <- function(x, obj, l)
{
## Define variables with shorter names for easier handling
for (name in names(obj$call.args))
	assign(obj$call.args[[name]], name)
v <- obj$v
normvar <- switch(objective, cov = "norm", cor = "var")

## Deflate data as needed
if (l > 1 && ortho.type == "weight" && scale == "block") {
	ortho.mode <- ortho.cnstr
	ortho.cnstr <- set.ortho.mat(v = obj$v[,1:(l-1)], 
		modes = ortho.mode[, 1:(l-1), l])
	x <- deflate.x(x, v = v[, 1:(l-1)], 
		ortho.mode = ortho.mode[, 1:(l-1), l], 
		check.args = FALSE)	
}

## Initialize canonical weights
if (!is.null(init.val)) {
	v0 <- if (is.matrix(init.val)) {
		init.val[, min(l, ncol(init.val))]
	} else init.val
	if (l > 1 && ortho.type == "weight" && scale == "block") 
		v0 <- tnsr.rk1.mat.prod(v0, mat = ortho.cnstr$mat, 
			modes = ortho.cnstr$modes, transpose.mat = FALSE)
	v0 <- scale.v(v0, type = normvar, scale = scale, 
		check.args = FALSE)
} else {
	mcca.init <- switch(init.method, 
		cca = mcca.init.cca, 
		svd = mcca.init.svd, 
		random = mcca.init.random)
	init.args$x <- x
	v0 <- do.call(mcca.init, init.args)
}
if (ortho == "score" && l > 1) {
	v0 <- tnsr.rk1.ortho(v0, ortho.cnstr[, 1:(l-1), drop = FALSE], 
		maxit = 100L, tol = 1e-6)
	v0 <- scale.v(v0, type = normvar, x = x, check.args = FALSE)
}

## Run MCCA
out <- if (objective == "cov" && optim == "bca" && scale == "block") {
	mcca.cov.bca.block(x = x, v = v0, w = w, 
		ortho = if (ortho.type == "score" && l > 1) {
		ortho.cnstr[,1:(l-1)] } else NULL, 
		sweep = sweep, maxit = maxit, tol = tol, 
		verbose = FALSE)
} else if (objective == "cov" && optim == "bca" && scale == "global") { 
	mcca.cov.bca.global(x = x, v = v0, w = w, 
		ortho = if (l == 1) { NULL 
		} else if (ortho == "weight") { v[, 1:(l-1)] 
		} else ortho.cnstr[,1:(l-1)], sweep = sweep, 
		maxit = maxit, tol = tol, verbose = FALSE)
} else if (objective == "cor" && optim == "bca") {
	mcca.cor.bca(x = x, v = v0, w = w, 
		ortho = if (ortho == "score" && l > 1) {
			ortho.cnstr[, 1:(l-1), drop = FALSE] } else NULL,
		sweep = sweep, maxit = maxit, tol = tol, verbose = verbose)	
} else if (optim == "grad.scale") {
	mcca.gradient.scale(x = x, v = v0, w = w, scale = scale, 
		type = normvar, maxit = maxit, tol = tol, verbose = verbose)
} else if (optim == "grad.rotate") {  
	mcca.gradient.rotate(x = x, v = v0, w = w, 
		maxit = maxit, tol = tol, verbose = verbose)
}	

out$objective

}


