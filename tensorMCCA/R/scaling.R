


#########################################
# Function for scaling canonical vectors
#########################################

scale.v <- function(v, x = NULL, type = c("norm", "var"), 
	cnstr = c("block", "global"), balance = TRUE)
{
type <- match.arg(type)   # norm or variance constraints
cnstr <- match.arg(cnstr) # block or global constraints
m <- length(v)
d <- sapply(v, length)
eps <- 1e-15 # numerical tolerance for zero

## If variance constraints, calculate variances and rescale
if (type == "var") {
	y <- image.scores(x, v)
	sdy <- sqrt(colMeans(y^2))
	for (i in 1:m) {
		if (sdy[i] < eps) next
		for (k in 1:d[i]) 
			v[[i]][[k]] <- v[[i]][[k]] / sdy[i]^(1/d[i])
	}
	if (!balance) return(v)	
}

## Calculate Frobenius norms of canonical tensors
## (= products of Euclidean norms of canonical vectors)
nrmv <- vector("list", m) 
nrmt <- numeric(m)
nrmfun <- function(vv) sqrt(sum(vv^2))
for (i in 1:m) {
	nrmv[[i]] <- sapply(v[[i]], nrmfun)
	nrmt[i] <- prod(nrmv[[i]])
}

## If variance constraints and balancing requirement,
## make norms of canonical vectors equal in each block
## without changing their product  
if (type == "var") {
	for (i in 1:m) {
		s <- if (nrmt[i] < eps) { numeric(d[i])
		} else { nrmt[i]^(1/d[i]) / nrmv[[i]] }
		for (k in 1:d[i]) 
			v[[i]][[k]] <- s[k] * v[[i]][[k]] 
	}
}

## If tensor norm constraints at the block level,
## rescale each canonical vector by its norm
if (type == "norm" && cnstr == "block") {
	for (i in 1:m) {
		s <- if (nrmt[i] < eps) { numeric(d[i]) 
			} else if (balance) { 1 / nrmv[[i]]
			} else { rep(1 / nrmt[i]^(1/d[i]), d[i]) } 
		for (k in 1:d[i]) 
			v[[i]][[k]] <- s[k] * v[[i]][[k]]	
	}
} 

## If tensor norm constraints at the global level,
## rescale each canonical vector in a given mode/dimension
## by the global norm of all canonical vectors in this mode
if (type == "norm" && cnstr == "global") {
	global.s <- mean(nrmt^2)
	global.s <- if (global.s < eps) 0 else (1/sqrt(global.s))
	for (i in 1:m) {
		s <- if (balance) {
			(global.s * nrmt[i])^(1/d[i]) / nrmv[[i]]
		} else { rep(global.s^(1/d[i]), d[i]) }
		s[is.nan(s)] <- 0
	for (k in 1:d[i]) 
		v[[i]][[k]] <- s[k] * v[[i]][[k]]	
	}
} 

return(v)
}




