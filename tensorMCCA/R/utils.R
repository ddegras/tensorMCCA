tnsr.rk1.nrm  <- function(v, norm = c("block", "global"))
{
stopifnot(is.list(v))
stopifnot(is.vector(v) || is.matrix(v))
norm <- match.arg(norm)
lenv <- length(v) 
out <- numeric(lenv)
sqnrmfun <- function(x) sum(x^2)
for (i in 1:lenv) 
	out[i] <- prod(sapply(v[[i]], sqnrmfun))
dim(out) <- dim(v)
if (norm == "block") return(sqrt(out))
if (is.vector(out)) return(sqrt(mean(out))) 
sqrt(rowMeans(out))
}



tnsr.rk1.cp <- function(x, y = NULL)
{
cpfun <- function(x, y) sum(x * y)
	
if (is.null(y)) {
	if (!is.matrix(x)) x <- as.matrix(x)
	m <- nrow(x)
	r <- ncol(x)
	out <- array(dim = c(r, r, m))
	for (i in 1:m) {
		for (k in 1:r) {
			for (l in 1:k) {
				out[k,l,i] <- prod(mapply(cpfun, x[[i,k]], x[[i,l]]))
				out[l,k,i] <- out[k,l,i]
			}
		}
	}
} else {
	stopifnot(length(x) == length(y))
	m <- length(x)
	out <- numeric(m)
	for (i in 1:m)
		out[i] <- prod(mapply(cpfun, x[[i]], x[[i]]))
}

out
}
