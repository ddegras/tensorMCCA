


#########################################
# Function for scaling canonical weights
#########################################

scale.v <- function(v, type = c("norm", "var"), 
	scale = c("block", "global"), x = NULL, check.args = TRUE)
{
type <- match.arg(type) 
scale <- match.arg(scale)
stopifnot(type == "norm" || !is.null(x))
stopifnot(is.list(v))
islist.v1 <- is.list(v[[1]])
if (!islist.v1) v <- list(v)
m <- NROW(v)
r <- NCOL(v)
isvec.v <- is.null(dim(v))
if (isvec.v) dim(v) <- c(m, r)
if (!is.null(x) && !is.list(x)) x <- list(x) 
if (check.args && !is.null(x)) test <- check.arguments(x, v)
d <- sapply(v[, 1], length)
eps <- 1e-15 # numerical tolerance for zero
p <- vector("list", m)
for (i in 1:m) p[[i]] <- sapply(v[[i,1]], length)

## Calculate norms of canonical tensors
## (= products of Euclidean norms of canonical vectors)
nrmv <- vector("list", m * r) 
dim(nrmv) <- c(m, r)
nrmt <- matrix(nrow = m, ncol = r)
nrmfun <- function(vv) sqrt(sum(vv^2))
for (i in 1:m) {
	for (l in 1:r) {
		nrmv[[i,l]] <- sapply(v[[i,l]], nrmfun)
		nrmt[i,l] <- prod(nrmv[[i,l]])
	}
}

if (type == "norm" && scale == "block") {
	## rescale each canonical vector by its norm
	zero <- (nrmt < eps)
	scalefun <- function(x,y) as.numeric(x) / y
	for (i in 1:m) {
		for (l in 1:r) {
			v[[i,l]] <- if (zero[i,l]) { 
				lapply(p[[i]], numeric) 
			} else { 
				mapply(scalefun, v[[i,l]], nrmv[[i,l]], SIMPLIFY = FALSE) 
			}
		}
	}
} else if (type == "norm" && scale == "global") {
	## rescale each canonical tensor by the same factor 
	## in a way that its canonical vectors are balanced 
	sl <- sqrt(colMeans(nrmt^2)) # grand scaling factor
	zero <- (nrmt < eps)
	scalefun <- function(x,y) as.numeric(x) * y
	for (i in 1:m) {
		for (l in 1:r) {
			if (zero[i,l]) { 
				v[[i,l]] <- lapply(p[[i]], numeric) 			
			} else { 
				sil <- (nrmt[i,l] / sl[l])^(1/d[i]) / nrmv[[i,l]]
				v[[i,l]] <- mapply(scalefun, v[[i,l]], sil, SIMPLIFY = FALSE) 
			}
		}
	}
} else { # type == "var"
	score <- canon.scores(x, v)
	n <- dim(score)[1]
	dim(score) <- c(n, m, r)
	sd.score <- colMeans(score^2) - colMeans(score)^2
	sd.score <- sqrt(pmax(sd.score, 0))
	zero <- (nrmt < eps | sd.score < eps)
	scalefun <- function(x,y) as.numeric(x) * y
	for (i in 1:m) {
		for (l in 1:r) {
			if (zero[i,l]) { 
				v[[i,l]] <- lapply(p[[i]], numeric) 
			} else { 
				sil <- (nrmt[i,l] / sd.score[i,l])^(1/d[i]) / nrmv[[i,l]]
				v[[i,l]] <- mapply(scalefun, v[[i,l]], sil, SIMPLIFY = FALSE) 
			}
		}
	}	
}

if (isvec.v) dim(v) <- NULL
if (!islist.v1) v <- v[[1]]
v
}




