set.ortho.mode <- function(x, r, cnstr = NULL,  
	method = c("cyclical", "random", "alldim", "maxdim"))
{
# test <- check.arguments(x)
stopifnot(r >= 1)
r <- as.integer(r)
method <- match.arg(method)

m <- length(x)
dimx <- lapply(x, dim)
p <- lapply(dimx, function(idx) idx[-length(idx)])	
d <- sapply(p, length)
n <- tail(dimx[[1]], 1)
pp <- sapply(p, prod)

ortho.mode <- vector("list", m * r^2)
dim(ortho.mode) <- c(m, r, r)
if (r == 1) return(ortho.mode)

if (!is.null(cnstr)) {
	if (!(length(cnstr) %in% c(1,m)))
		stop(paste("'cnstr' must be a single integer", 
			"or a vector of integers with same length as 'x'")) 
	cnstr <- as.integer(cnstr)
	cnstr <- rep_len(cnstr, m)
	if (!(all(cnstr > 0) && all(cnstr <= d)))
		stop(paste("The values in 'cnstr' must be positive integers",
		"less or equal to the number of modes/dimensions in the",
		"corresponding component of 'x'"))
	for (i in 1:m) {
		ortho.mode[i,,] <- list(cnstr[i])
		diag(ortho.mode[i,,]) <- vector("list", r)
	}
} else if (method %in% c("cyclical", "random")) {
	mask.up <- which(upper.tri(matrix(, r, r)), arr.ind = TRUE)
	mask.lo <- mask.up[,2:1]
	rc2 <- choose(r, 2)
	mat <- vector("list", r^2)
	dim(mat) <- c(r, r)
	for (i in 1:m) {
		vals <- switch(method, cyclical = rep_len(1:d[i], rc2),
			random = sample.int(d[i], rc2, replace = TRUE))
		mat[mask.up] <- mat[mask.lo] <- as.list(vals)
		ortho.mode[i,,] <- mat
	}
} else {
	for (i in 1:m) {
		vals <- switch(method, alldim = 1:d[i], 
			maxdim = which.max(p[[i]]))
		ortho.mode[i,,] <- list(vals)
		diag(ortho.mode[i,,]) <- vector("list", r)		
	}
}

ortho.mode
}


set.ortho.mat <- function(v, modes)
{
if (!is.matrix(v)) v <- as.matrix(v)
if (!is.matrix(modes)) 
	modes <- as.matrix(modes)
stopifnot(all(dim(v) == dim(modes)))
m <- nrow(v)
ortho.mat <- out.modes <- vector("list", m)
vzero <- logical(m)
for (i in 1:m) {
	idxi <- modes[i, ]
	di <- length(v[[i]]) 
	mat <- vector("list", di)
	for (k in 1:di) {
		idxik <- which(sapply(idxi, is.element, el = k))
		if (length(idxik) == 0) next
		vik <- sapply(v[i,idxik], "[[", k)
		qrvik <- qr(vik)
		pik <-  nrow(vik)
		rik <- qrvik$rank
		if (rik > 0 && rik < pik) {
			mat[[k]] <- t(qr.Q(qrvik, TRUE)[, -(1:rik)])
		} else if (rik == pik) {
			mat[[k]] <- matrix(0, 1, nrow(vik))
			vzero[i] <- TRUE
		}
	}
	modei <- which(!sapply(mat, is.null))
	ortho.mat[[i]] <- mat[modei]
	out.modes[[i]] <- modei
}

list(mat = ortho.mat, modes = out.modes, vzero = vzero)
}
