###########################################
# Function to calculate gradients of 
# objective and constraint functions, 
# Lagrange multipliers, and KKT residuals 
##########################################

kkt <- function(x, fit)
{
## Check compatibility of arguments
test <- check.arguments(x, fit$v, fit$call.args$w)

## Extract relevant quantities
v <- fit$v
score <- fit$block.score
w <- fit$call.args$w
scale <- fit$call.args$scale
ortho <- fit$call.args$ortho.type
objective <- fit$call.args$objective.type

## Calculate data dimensions
m <- length(x)
dimfun <- function(x) {
	if (length(dim(x)) == 0) c(1L, length(x)) else dim(x) }
dimx <- lapply(x, dimfun)
d <- sapply(dimx, length) - 1L 
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
n <- tail(dimx[[1]], 1L)
r <- NCOL(v)

## Auxiliary function for gradient calculation
gradfun <- function(x, y = NULL) {
if (is.null(y)) y <- x
d <- length(x)
if (d == 1) return(y[[1]])
cp <- numeric(d)
for (k in 1:d) cp[k] <- sum(x[[k]] * y[[k]])
out <- vector("list", d)
for (k in 1:d) out[[k]] <- prod(cp[-k]) * y[[k]]
unlist(out)	
}


## Calculate gradients
grad.obj <- grad.cnstr <- vector("list", m * r)
dim(grad.obj) <- dim(grad.cnstr) <- c(m, r)
sump <- sapply(p, sum)
score <- score / n
global.score <- if (ortho == "score" && scale = "global") {
	fit$global.score / n } else NULL
for (i in 1:m) {
	for (l in 1:r) {
		## Gradient of objective function
		tvprod <- vector("list", d[i])
		scorei <- score[,,l] %*% w[,i] 
		for (k in 1:d[i]) {
			tvprod[[k]] <- tnsr.vec.prod(x = x[[i]], 
				v = v[[i,l]][-k], modes = (1:d[i])[-k])
			if (p[[i]][k] == 1) dim(tvprod[[k]]) <- c(1,n)
		}
		tvprod <- do.call(rbind, tvprod)
		dim(tvprod) <- c(sump[i], n)
		# tvprod <- tvprod - rowMeans(tvprod)
		grad.obj[[i,l]] <- tvprod %*% scorei
		
		## Gradient of scaling constraint
		mat <- matrix(nrow = sump[i], ncol = l)
		mat[,l] <- if (objective == "cor") {
			tvprod %*% score[,i,l]	
		} else if (objective == "cov" && scale == "block") {
			unlist(v[[i,l]])
		} else {
			gradfun(v[[i,l]])
		}
		if (l == 1) {
			grad.cnstr[[i,l]] <- mat
			next
		}	
			
		## Gradient of orthogonality constraints
		if (ortho == "weight") {
			for (ll in 1:(l-1))
				mat[,ll] <- gradfun(v[[i,l]], v[[i,ll]])			
		} else if (ortho == "score" && scale == "block") {
			mat[,1:(l-1)] <- tvprod %*% score[,i,1:(l-1)]
		} else { # ortho == "score" && scale == "global"
			mat[,1:(l-1)] <- tvprod %*% global.score[,i,1:(l-1)] 
		}		
		grad.cnstr[[i,l]] <- mat			
	}
}

if (scale == "global") {
	for (l in 1:r) {
		grad.obj[[1,l]] <- do.call(rbind, grad.obj[,l])
		grad.cnstr[[1,l]] <- do.call(rbind, grad.cnstr[,l])
	}
	grad.obj <- grad.obj[1,]
	grad.cnstr <- grad.cnstr[1,]
}

## Regress objective gradient against constraints gradients
if (scale == "block") {
	coefs <- vector("list", m * r)
	dim(coefs) <- c(m, r)
	residuals <- coefs
	for (i in 1:m) {
		for (l in 1:r) {
			reg <- suppressWarnings(lsfit(grad.cnstr[[i,l]], 
				grad.obj[[i,l]], NULL, FALSE))
			coefs[[i,l]] <- reg$coefficients
			residuals[[i,l]] <- reg$residuals
		}
	}
	for (l in 1:r) 
		residuals[[1,l]] <- unlist(residuals[,l])	
	residuals <- residuals[1,]
} else { # scale == "global"
	coefs <- residuals <- vector("list", r)
	for (l in 1:r) {
		reg <- suppressWarnings(lsfit(grad.cnstr[[l]], 
			grad.obj[[l]], NULL, FALSE))
		coefs[[l]] <- reg$coefficients
		residuals[[l]] <- reg$residuals
	}
}

list(grad.objective = grad.obj, grad.constraints = grad.cnstr, 
	multipliers = coefs, residuals = residuals)	
}


