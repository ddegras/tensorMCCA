


########################
# Random initialization 
# of canonical vectors 
########################


init.mcca.random <- function(x, objective = c("covariance", "correlation"))
{
## Check argument x if required
test <- check.arguments(x)

## Data dimensions
m <- length(x) # number of datasets 
dimx <- lapply(x, dim) # full data dimensions
p <- lapply(dimx, function(idx) idx[-length(idx)]) 
# image dimensions (ignore last dimension of datasets = replications)
d <- sapply(p, length)

## Scaling constraints
objective <- match.arg(objective) # objective function to maximize

## Initialize canonical vectors randomly
v <- vector("list",m)
for (i in 1:m)
for (k in 1:d[i])
	v[[i]][[k]] <- runif(p[[i]][k], -1, 1)

## Scale canonical vectors
v <- scale.v(v, cnstr = "block")
if (objective == "correlation") {
	y <- image.scores(x, v)
	y <- y - rowMeans(y)
	nrm <- sqrt(colMeans(y^2))
	nrm[nrm <= 1e-14] <- 1
	for (i in 1:m)
	for (k in 1:d[i])
		v[[i]][[k]] <- v[[i]][[k]] / nrm[k]^(1/d[i])	
}

v
}


