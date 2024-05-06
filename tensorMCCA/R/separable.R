###########################################################
# This internal function tests whether a symmetric matrix 
# with nonnegative entries and at least one postive entry 
# is separable
###########################################################

# If 'freediag' is set to TRUE, the diagonal elements of the matrix 'x'
# are variable, i.e., they can be adjusted to try to make 'x' separable


separable <- function(x, freediag = FALSE)
{
stopifnot(is.numeric(x))
stopifnot(all(!is.na(x)) && all(!is.nan(x)))
stopifnot(all(x >= 0) && any(x > 0))
stopifnot(length(x) == 1 || (is.matrix(x) && nrow(x) == ncol(x)))
stopifnot(x == t(x))

## Trivial case: scalar argument
if (length(x) == 1) {
	out <- list(separable = TRUE, a = sqrt(x))
	return(out)
}
	
m <- ncol(x)

## Case: all matrix elements are fixed
if (!freediag) {
	eigx <- eigen(x, TRUE)
	eigvals <- eigx$values 
	out <- if (all(eigvals[-1] <= max(1e-15, 1e-12 * eigvals[1]))) {
		list(separable = TRUE, a = sqrt(eigvals[1]) * abs(eigx$vectors[,1]))
	} else {
		list(separable = FALSE, a = NULL)
	}
		return(out)
}	

## Case: nondiagonal matrix elements are fixed, diagonal elements are free


if (m == 2) {
	out <- list(separable = TRUE, a = rep(sqrt(x[1,2]),2))
	return(out)
}
	
diag(x) <- NA
elemzero <- which(x == 0, TRUE)
colzero <- apply((x == 0) | is.na(x), 2, all)

if (length(elemzero) > 0) {
	## If element (i,j) is zero, then i-th or j-th column must be zero
	## for symmetric matrix to be separable
	test <- (colzero[elemzero[,1]] | colzero[elemzero[,2]])
	if (!test || all(colzero)) {
		out <- list(separable = FALSE, a = NULL)
		return(out)	
	}
	## Reduce matrix to row/cols with nonzero (fixed) elements
	x <- x[!colzero,!colzero]
}

mm <- NCOL(x) # number of rows/cols with nonzero (fixed) elements

## Case: single row/column with nonzero (fixed) elements
if (mm == 1) {
	a <- rep(0, m)
	a[which(!colzero)] <- sqrt(x)
	out <- list(separable = TRUE, a = a)
	return(out)
}

## General case: log-transform positive entries and check if they
## satisfy associated (over-)determined system of linear equations 
idx <- t(combn(mm,2))
logx <- log(x[idx])
cmm2 <- nrow(idx)
mat <- matrix(0, cmm2, mm)
mat[cbind(rep(1:cmm2, each = 2), as.integer(t(idx)))] <- 1
logaa <- tryCatch(solve(mat, logx), error = function(e) NULL)
if (is.null(logaa)) {
	out <- list(separable = FALSE, a = NULL)
	return(out)
}

a <- numeric(m)
a[!colzero] <- exp(logaa)
out <- list(separable = TRUE, a = a)
out

}