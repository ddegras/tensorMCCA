check.arguments <- function(x, v = NULL, w = NULL)
{
## Check argument types
stopifnot(is.list(x))
if (!(is.null(v) || is.list(v)))
	stop("Argument 'v' must be left unspecified (NULL) or be a list.")
stopifnot(is.null(w) || is.numeric(w))
	
## Data dimensions
if (!(is.null(v) || length(x) == length(v)))
	stop("If specified, 'v' must have the same length as 'x'")
dimx <- lapply(x, dim) # full data dimensions

## Check that all datasets have same number of instances/individuals
n <- sapply(dimx, tail, 1) # numbers of instances per dataset
if (!all(n == n[1]))
	stop(paste("All components of 'x' must have", 
		"the same dimension in the last mode."))

## Check that there are at least 2 datasets and 2 instances 
if (length(x) < 2)
	stop("Make sure that 'x' has length at least 2 (number of datasets).")
if (n[1] < 2) 
	stop(paste("Make sure that the arrays in 'x'",
		"have their last dimension at least 2 (number of instances)."))

## Check that the data arrays have at most 4 dimensions
## (i.e. at most 3 image dimensions + 1 dimension for instances)
dim.img <- lapply(dimx, function(idx) idx[-length(idx)]) 
ndim.img <- sapply(dim.img, length)
if (any(ndim.img > 3)) 
	stop(paste("Make sure that the arrays in 'x' have",
		"at most 4 dimensions: 1, 2, or 3 for images and 1",
	  	"for instances/individuals"))

## Check that the data contain no missing values (NA) or NaN
test <- sapply(x, function(xx) any(is.na(xx) | is.nan(xx)))
if (any(test)) 
	stop("The data arrays in 'x' must not contain NA or NaN values.")

## Check optional argument 'v'
if (length(v) > 0) {
	dimv <- lapply(v, function(l) sapply(l, NROW))
	if (!identical(dim.img, dimv))
		stop(paste("Make sure that each component of 'v' is a list",
			"of vectors or column matrices whose length match the",
			"corresponding dimensions of 'x'."))
			
	r <- unlist(lapply(v, function(l) lapply(l, NCOL)))
	if (!all(r == r[1]))
		stop(paste("The matrices contained in 'v' should all",
			"have the same number of columns."))
			
	if (any(is.na(unlist(v))) || any(is.nan(unlist(v))))
		stop(paste("The canonical vectors in 'v' must not contain",
			"NA or NaN values."))
}

## Check optional argument 'w'
if (!is.null(w)) {
	test1 <- (length(w) == 1)
	test2 <- (is.matrix(w) && all(dim(w) == length(x)))
	if (!(test1 || test2)) 
		stop(paste("If specified, 'w' must either be:\n* A single number, or", 
		"\n* A square matrix with dimensions equal", 
		"to the length of 'x' (the number of datasets)"))
	stopifnot(all(w >= 0) && any(w > 0))	
}

NULL
}
