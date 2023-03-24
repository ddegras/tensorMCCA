######################################
# Function to calculate KKT residuals
######################################

kkt <- function(x, object)
{
## Check compatibility of arguments
test <- check.arguments(x, object$v, object$call.args$w)

## Extract fitted objects
v <- object$v[,1]
score <- object$block.score[,,1]
w <- object$call.args$w

## Calculate data dimensions
m <- length(x)
dimx <- lapply(x, dim)
d <- sapply(dimx, length) - 1L
n <- tail(dimx[[1]], 1L)

## Precalculate Lagrange multipliers
objective.type <- object$call.args$objective.type
scale <- object$call.args$scale
if (objective.type == "cov" && scale == "global") {
	lambda <- vector("list", m)
	for (i in 1:m) {
		if (d[i] == 1) {
			lambda[[i]] <- rep(objective / m, d[i])
		} else {
			sqnrmv <- sapply(v[[i]], function(x) sum(x^2))
			for (k in 1:d[i]) 
				lambda[[i]][k] <- objective * prod(sqnrmv[-k]) / m 
		}
	}			
} else {
	lambda <- colSums(crossprod(score) * w)  / n
	lambda <- mapply(rep, lambda, d, SIMPLIFY = FALSE)
}
	
## Calculate KKT residuals 
residuals <- vector("list", m)
for (i in 1:m) {
	for (k in 1:d[i]) {
		## Calculate gradient of objective function
		tvprod <- tnsr.vec.prod(x[[i]], v[[i]][-k], (1:d[i])[-k])
		grad <- tvprod %*% score %*% (w[,i] / n) 		

		residuals[[i]][[k]] <- if (objective.type == "cov") {
			grad - lambda[[i]][k] * v[[i]][[k]]
		} else {
			grad - lambda[[i]][k] * (tvprod %*% score[,i] / n)	
		}		
	}	
}

residuals
}