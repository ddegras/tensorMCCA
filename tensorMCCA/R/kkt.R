######################################
# Function to calculate KKT residuals
######################################

kkt <- function(x, object)
{
## Check compatibility of arguments
test <- check.arguments(x, object$v, object$call.args$w)

## Extract relevant quantities
v <- object$v
score <- object$block.score
w <- object$call.args$w
scale <- object$call.args$scale
ortho <- object$call.args$ortho.type
obj <- object$call.args$objective.type

## Calculate data dimensions
m <- length(x)
dimx <- lapply(x, dim)
d <- sapply(dimx, length) - 1L 
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
n <- tail(dimx[[1]], 1L)
r <- NCOL(v)

## Auxiliary function for gradient calculation
if (ortho == "weight" && r > 1) {
	gradfun <- function(x, y) {
	d <- length(x)
	if (d == 1) return(y[[1]])
	cp <- numeric(d)
	for (k in 1:d) cp[k] <- sum(x[[k]] * y[[k]])
	out <- vector("list", d)
	for (k in 1:d) out[[k]] <- prod(cp[-k]) * y[[k]]
	unlist(out)	
	}
}

## Calculate gradients
grad.obj <- grad.cnstr <- vector("list", m * r)
dim(grad.obj) <- dim(grad.cnstr) <- c(m, r)
sump <- sapply(p, sum)
for (i in 1:m) {
	for (l in 1:r) {
		## Gradient of objective function
		tvprod <- vector("list", d[i])
		for (k in 1:d[i]) 
			tvprod[[k]] <- tnsr.vec.prod(x = x[[i]], 
				v = v[[i,l]][-k], modes = (1:d[i])[-k])
		tvprod <- do.call(rbind, tvprod)
		# dim(tvprod) <- c(sump[i], n)
		tvprod <- tvprod - rowMeans(tvprod)
		grad.obj[[i,l]] <- tvprod %*% score[,,l] %*% (w[,i] / n)
		## Gradient of scaling constraint
		mat <- matrix(nrow = sump[i], ncol = l)
		mat[,l] <- if (obj == "cor") {
			tvprod %*% score[,i,l]	
		} else {
			unlist(v[[i,l]])
		}
		if (l == 1) {
			grad.cnstr[[i,l]] <- mat
			next
		}		
		## Gradient of orthogonality constraints
		if (ortho == "score") {
			mat[,1:(l-1)] <- tvprod %*% score[,i,1:(l-1)]
		} else {
			for (ll in 1:(l-1))
				mat[,ll] <- gradfun(v[[i,l]], v[[i,ll]])
		}		
		grad.cnstr[[i,l]] <- mat			
	}
}

## Calculate regression residuals
residuals <- vector("list", m * r)
dim(residuals) <- c(m, r)
for (i in 1:m)
for (l in 1:r)
residuals[[i,l]] <- suppressWarnings(lsfit(grad.cnstr[[i,l]], 
	grad.obj[[i,l]], NULL, FALSE)$residuals)
	
list(grad.objective = grad.obj, grad.constraints = grad.cnstr, 
	residuals = residuals)	
}


