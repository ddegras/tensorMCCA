\name{check.arguments}
\alias{check.arguments}

\title{
Check Data and Initial Canonical Tensors for Multiple CCA
}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}

\usage{check.arguments(x, v = NULL)}

\arguments{
  \item{x}{
Data: list of arrays. 
}
\item{v}{
Canonical vectors (optional): list of lists of vectors/matrices. 
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}

\value{
0 in case of successful test, error message otherwise
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, v = NULL) 
{
    if (!is.list(x)) 
        stop("Argument 'x' must be a list.")
    if (!(is.null(v) || is.list(v))) 
        stop("Argument 'v' must be left unspecified (NULL) or be a list.")
    if (!(is.null(v) || length(x) == length(v))) 
        stop("If specified, 'v' must have the same length as 'x'")
    dimx <- lapply(x, dim)
    n <- sapply(dimx, tail, 1)
    if (!all(n == n[1])) 
        stop(paste("All components of 'x' must have", "the same dimension in the last mode."))
    if (length(x) < 2) 
        stop("Make sure that 'x' has length at least 2 (number of datasets).")
    if (n[1] < 2) 
        stop(paste("Make sure that the arrays in 'x'", "have their last dimension at least 2 (number of instances)."))
    dim.img <- lapply(dimx, function(idx) idx[-length(idx)])
    ndim.img <- sapply(dim.img, length)
    if (any(ndim.img > 3)) 
        stop(paste("Make sure that the arrays in 'x' have", "at most 4 dimensions: 1, 2, or 3 for images and 1", 
            "for instances/individuals"))
    test <- sapply(x, function(xx) any(is.na(xx) | is.nan(xx)))
    if (any(test)) 
        stop("The data arrays in 'x' must not contain NA or NaN values.")
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
    return(0)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
