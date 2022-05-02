


########################
# Random initialization 
# of canonical vectors 
########################


init.mcca.random <- function(x, objective = c("cov", "cor"), 
	cnstr = c("block", "global"), balance = TRUE)
{
## Check argument x if required
test <- check.arguments(x)

## Data dimensions
m <- length(x) # number of datasets 
dimx <- lapply(x, dim) # full data dimensions
dim.img <- lapply(dimx, function(idx) idx[-length(idx)]) 
# image dimensions (ignore last dimension of datasets = replications)
ndim.img <- sapply(dim.img, length)

## Scaling constraints
objective <- match.arg(objective) # objective function to maximize
type <- switch(objective, cov = "norm", cor = "var") # norm or variance constraints
cnstr <- match.arg(cnstr) # block or global constraints
if (objective == "cor") cnstr <- "block"

## Initialize canonical vectors randomly
v <- vector("list",m)
for (i in 1:m)
	v[[i]] <- lapply(dim.img[[i]], function(len) runif(len))

## Scale canonical vectors
v <- scale.v(v, x, type, cnstr, balance)
return(v)
}


