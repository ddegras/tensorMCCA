


########################
# Random initialization 
# of canonical vectors 
########################


mcca.init.random <- function(x, objective = c("covariance", "correlation"))
{
test <- check.arguments(x)
eps <- 1e-14
## Data dimensions
m <- length(x) # number of datasets 
dimx <- lapply(x, dim) # full data dimensions
p <- lapply(dimx, function(idx) idx[-length(idx)]) 
# image dimensions (ignore last dimension of datasets = replications)
d <- sapply(p, length)

## Scaling constraints
objective <- match.arg(objective) # objective function to maximize

## Initialize canonical vectors randomly
v <- vector("list", m)
for (i in 1:m)
	v[[i]] <- lapply(p[[i]], runif, a = -1, b = 1)

## Scale canonical vectors
if (objective == "covariance") {
	v <- scale.v(v, cnstr = "block")
} else {
	y <- canon.scores(x, v)
	ybar <- colMeans(y)
	nrm <- sqrt(colMeans(y^2) - ybar^2)
	nrm[is.nan(nrm)] <- 0
	nz <- which(nrm > eps)
	for (i in 1:m) {
		v[[i]] <- if (nz[i]) {
			lapply(v[[i]], "/", y = nrm[i]^(1/d[i]))
		} else { lapply(p[[i]], numeric) }
	}
}

v
}


