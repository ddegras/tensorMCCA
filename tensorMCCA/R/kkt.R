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
sump <- sapply(p, sum)
cumsump <- cumsum(c(0, sump))
pp <- cumsump[m+1]
grad.obj <- matrix(0, pp, r)
grad.cnstr <- vector("list", r)
score <- score / n
if (ortho == "score" && scale == "global") 
	global.score <- fit$global.score / n
for (l in 1:r) {
	grad.cnstr[[l]] <- switch(scale, 
		block = matrix(0, pp, l*m), 
		global = matrix(0, pp, l))
	for (i in 1:m) {
		idxrow <- (cumsump[i]+1):cumsump[i+1]
		idxcol <- switch(scale, 
			block = seq(i, by = m, len = l),
			global = 1:l)
	
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
		grad.obj[idxrow,l] <- tvprod %*% scorei
		
		## Gradient of scaling constraint
		mat <- matrix(nrow = sump[i], ncol = l)
		mat[,l] <- if (objective == "cor") {
			tvprod %*% score[,i,l]	
		} else if (objective == "cov" && scale == "block") {
			unlist(v[[i,l]])
		} else {
			gradfun(v[[i,l]])
		}
			
		## Gradient of orthogonality constraints
		if (l > 1) {
			if (ortho == "weight") {
				for (ll in 1:(l-1))
					mat[,ll] <- gradfun(v[[i,l]], v[[i,ll]])			
			} else if (ortho == "score" && scale == "block") {
				mat[,1:(l-1)] <- tvprod %*% score[,i,1:(l-1)]
			} else { # ortho == "score" && scale == "global"
				mat[,1:(l-1)] <- tvprod %*% global.score[,i,1:(l-1)] 
			}		
		}
		
		grad.cnstr[[l]][idxrow,idxcol] <- mat
	}
}

## Regress objective gradient against constraints gradients
coefs <- vector("list", r)
residuals <- matrix(0, pp, r)
for (l in 1:r) {
	reg <- suppressWarnings(
		lsfit(grad.cnstr[[l]], grad.obj[,l], NULL, FALSE))
	coefs[[l]] <- reg$coefficients
	if (scale == "block") dim(coefs[[l]]) <- c(m,l)
	residuals[,l] <- reg$residuals
}


list(grad.objective = grad.obj, grad.constraints = grad.cnstr, 
	multipliers = coefs, residuals = residuals)	
}


