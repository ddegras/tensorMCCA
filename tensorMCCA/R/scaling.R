


#########################################
# Function for scaling canonical vectors
#########################################

scale.v <- function(v, cnstr = c("block", "global"))
{
cnstr <- match.arg(cnstr) # block or global constraints
m <- length(v)
d <- sapply(v, length)
eps <- 1e-15 # numerical tolerance for zero

## Calculate Frobenius norms of canonical tensors
## (= products of Euclidean norms of canonical vectors)
nrmv <- vector("list", m) 
nrmt <- numeric(m)
nrmfun <- function(vv) sqrt(sum(vv^2))
for (i in 1:m) {
	nrmv[[i]] <- sapply(v[[i]], nrmfun)
	nrmt[i] <- prod(nrmv[[i]])
}

## If tensor norm constraints at the block level,
## rescale each canonical vector by its norm
if (cnstr == "block") {
	for (i in 1:m) {
		for (k in 1:d[i]) {
			v[[i]][[k]] <- if (nrmt[i] < eps) { 
				numeric(length(v[[i]][[k]])) 
			} else { v[[i]][[k]] / nrmv[[i]][k] }
		}
	}
} 

## If tensor norm constraints at the global level,
## rescale each canonical tensor by the same factor 
## in a way that its canonical vectors are balanced 
if (cnstr == "global") {
	s <- sqrt(mean(nrmt^2)) # grand scaling factor
	for (i in 1:m) {
		for (k in 1:d[i]) {
			sik <- (nrmt[i] / s)^(1/d[i]) / nrmv[[i]][k]
			v[[i]][[k]] <- if (nrmt[i] < eps) { 
				numeric(length(v[[i]][[k]])) 
			} else { v[[i]][[k]] * sik } 
		}
	}
} 

v
}




