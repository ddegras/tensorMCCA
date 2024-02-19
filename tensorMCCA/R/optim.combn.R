optim.combn.exact <- function(score, w)
{
n <- dim(score)[1]
m <- dim(score)[2]

## Simpler case m = 2
if (m == 2) {
	a <- matrix(0, m, m)
	for (i in 1:m)
	for (j in 1:m)
		a[i,j] <- w[1,1] * sum(score[,1,i] * score[,1,i]) +
			w[1,2] * sum(score[,1,i] * score[,2,j]) +
			w[2,1] * sum(score[,2,j] * score[,1,i]) +
			w[2,2] * sum(score[,2,j] * score[,2,j])
	idx <- arrayInd(which.max(a), dim(a))
	return(list(idx = c(idx), sign = c(1,1), objective = a[idx]/n))
} 


## Cross-products between scores
cp <- array(0, rep(m,4))
for (i in 1:m) { 
	for (j in i:m) {
		cp[,,i,j] <- if (i == j) {
			diag(w[i,i] * colSums(score[,i,]^2))
		} else {
			2 * w[i,j] * crossprod(score[,i,], score[,j,])
		}
	}
}

## Set of sign flips
flip <- lapply(0:(m-1), 
	function(size) as.list(data.frame(combn(2:m,size))))
flip <- unlist(flip, FALSE, FALSE)

objective.best <- -Inf
iperm <- integer(m)

for (f in 1:length(flip)) {
	anyflip <- (length(flip[[f]]) > 0)
	if (anyflip) {
		idx <- flip[[f]] 
		cp[,,idx,-idx] <- - cp[,,idx,-idx]
		cp[,,-idx,idx] <- - cp[,,-idx,idx]
	}
	a <- array(0, rep(m,m))
	for (i in 1:m) { 
		for (j in i:m) {
			if (i == j) {
				cpij <- diag(cp[,,i,j])
				perm <- c(i,(1:m)[-i])
			} else {
				cpij <- as.vector(cp[,,i,j])
				perm <- c(i,j,(1:m)[-c(i,j)])
			}
			a <- aperm(a, perm)
			a <- a + cpij
			iperm[perm] <- 1:m
			a <- aperm(a, iperm)	
		}
	}
	objective <- max(a)
	if (objective > objective.best) {
		objective.best <- objective
		idx.best <- arrayInd(which.max(a), dim(a))
		flip.best <- flip[[f]] 
	}
	if (anyflip) {
		idx <- flip[[f]] 
		cp[,,idx,-idx] <- - cp[,,idx,-idx]
		cp[,,-idx,idx] <- - cp[,,-idx,idx]
	}	
}

sign.best <- rep(1,m)
sign.best[flip.best] <- -1
return(list(idx = c(idx.best), sign = sign.best, objective = objective.best))

}




optim.combn.greedy <- function(score, w)
{
n <- dim(score)[1]
m <- dim(score)[2]

## Simpler case m = 2
if (m == 2) {
	a <- matrix(0, m, m)
	for (i in 1:m)
	for (j in 1:m)
		a[i,j] <- w[1,1] * sum(score[,1,i] * score[,1,i]) +
			w[1,2] * sum(score[,1,i] * score[,2,j]) +
			w[2,1] * sum(score[,2,j] * score[,1,i]) +
			w[2,2] * sum(score[,2,j] * score[,2,j])
	idx <- arrayInd(which.max(a), dim(a))
	return(list(idx = c(idx), sign = c(1,1), objective = a[idx]/n))
} 

## Cross-products between scores
for (i in 1:m)
for (j in 1:m)
if (w[i,j] == 0) score[,i,j] <- NA
cp <- array(NA, rep(m,4))
for (i in 1:m) { 
	for (j in i:m) {
		if (w[i,j] == 0) next
		if (i == j) {
			mat <- matrix(NA, m, m)
			diag(mat) <- w[i,i] * colSums(score[,i,]^2)
			cp[,,i,i] <- mat
		} else {
			cp[,,i,j] <- w[i,j] * crossprod(score[,i,], score[,j,])
		}
	}
}

## Greedy search
idx.best <- integer(m)
cpcopy <- cp
while (!all(idx.best > 0)) {
	idx <- arrayInd(which.max(cp), dim(cp))
	dset <- idx[3:4]
	item <- idx[1:2]
	cp[,,dset[1],dset[2]] <- NA
	cp[-item[1],,dset[1],] <- cp[,-item[1],,dset[1]] <- NA
	cp[-item[2],,dset[2],] <- cp[,-item[2],,dset[2]] <- NA
	idx <- which(idx.best[dset] == 0)	
	idx.best[dset[idx]] <- item[idx]
}
cp <- cpcopy
cp[is.na(cp)] <- 0

objective <- matrix(,m,m)
for (i in 1:m) 
for (j in i:m) 
	objective[i,j] <- objective[j,i] <- 
		cp[idx.best[i], idx.best[j], i, j]
objective.best <- objective.sum <- sum(objective)
sign.best <- rep(1,m)

## Check if sign flips can improve objective
flip <- lapply(1:(m-1), 
	function(size) as.list(data.frame(combn(2:m,size))))
flip <- unlist(flip, FALSE, FALSE)
for (f in 1:length(flip)) {
	idx <- flip[[f]]
	objective.tmp <- objective.sum - 4 * sum(objective[idx,-idx])	
	if (objective.tmp > objective.best) { 
		objective.best <- objective.tmp
		sign.best <- rep(1,m)
		sign.best[idx] <- -1
	}
}

return(list(idx = c(idx.best), sign = sign.best, objective = objective.best))

}
