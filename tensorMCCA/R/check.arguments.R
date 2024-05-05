check.arguments <- function(x, v = NULL, w = NULL)
{
## Check argument types
stopifnot(is.list(x))
if (!(is.null(v) || is.list(v)))
	stop("Argument 'v' must be left unspecified (NULL) or be a list.")
stopifnot(is.null(w) || is.numeric(w))
	
## Data dimensions
if (!(is.null(v) || length(x) == NROW(v)))
	stop(paste("If specified, 'v' must be a vector of lists",
	"of same length as 'x' or a matrix of lists with as many rows",
	"as the length of 'x'."))
dimx <- lapply(x, dimfun) # full data dimensions
m <- length(x)

## Check that all datasets have same number of instances/individuals
n <- sapply(dimx, tail, 1L) # numbers of instances per dataset
if (!all(n == n[1L]))
	stop(paste("All components of 'x' must have", 
		"the same dimension in the last mode."))

## Check that there are at least 2 datasets and 2 instances 
if (length(x) < 2L)
	stop("Make sure that 'x' has length at least 2 (number of datasets).")
if (n[1] < 2L) 
	stop(paste("Make sure that the arrays in 'x'",
		"have their last dimension at least 2 (number of instances)."))

## Check that the data contain no missing values (NA) or NaN
test <- sapply(x, function(xx) any(is.na(xx) | is.nan(xx)))
if (any(test)) 
	stop("The data arrays in 'x' must not contain NA or NaN values.")

## Check optional argument 'v'
if (length(v) > 0) {
	r <- NCOL(v)
	if (!is.matrix(v)) dim(v) <- c(m, r)
	p <- lapply(dimx, function(idx) idx[-length(idx)]) 
	d <- sapply(p, length)
	for (i in 1:m) {
		for (l in 1:r) {
			if (length(v[[i, l]]) != d[i])
				stop("Component v[[", i, ",", l, "]] ",
				"must a list of length the number of dimensions ",
				"of x[[", i, "]] minus 1.")
			if (any(sapply(v[[i,l]], length) != p[[i]]))
				stop("Component v[[", i, ",", l, "]][[", k, "]] ",
				"must be a numerical vector of length ",
				"dim(x[[", i, "]])[", k, "].")
		}
	}		
	if (any(is.na(unlist(v))) || any(is.nan(unlist(v))))
		stop(paste("The canonical vectors in 'v' must not contain",
			"NA or NaN values."))
}

## Check optional argument 'w'
if (!is.null(w)) {
	test1 <- (length(w) == 1)
	test2 <- (is.matrix(w) && all(dim(w) == length(x)))
	if (!(test1 || test2)) 
		stop(paste("If specified, 'w' must either be a single number", 
		"or a square matrix with numbers of rows and columns equal", 
		"to the length of 'x' (the number of datasets)."))
	stopifnot(all(w >= 0) && any(w > 0))	
}

NULL
}
