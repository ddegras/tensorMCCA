





###########################################
# Internal function to calculate objective
###########################################

# This function calculates the -weighted- sum of covariances.
# For speed, it does not check the dimensions of arguments. 
# To calculate the sum of correlations with this function, 
# the canonical vectors should be scaled beforehand to have 
# unit variance, e.g., with function 'scale.v'.

objective.internal <- function(x, v, w)
{
n <- tail(dim(x[[1]]), 1) 
score <- canon.scores(x, v)
r <- NCOL(v)
if (r == 1) {
	return(sum(w * crossprod(score)) / n)	
}
out <- numeric(r)
for (l in 1:r) 
	out[l] <- sum(w * crossprod(score[, , l])) / n
sum(out)
}




############################################
# Function to calculate sum of correlations
############################################


objective.cor <- function(x, v, w = 1)
{
test <- check.arguments(x, v, w)
m <- length(x)
n <- tail(dim(x[[1]]), 1)
r <- NCOL(v)
w <- if (length(w) == 1) { 1 / (m^2)
} else { w <- (w + t(w)) / (2 * sum(w)) }
score <- canon.scores(x, v)
if (r == 1) {
	return(sum(w * cor(score)))
}
out <- numeric(r)
for (l in 1:r) 
	out[l] <- sum(w * cor(score[, , l]))
sum(out)
}




############################################
# Function to calculate sum of covariances
############################################


objective.cov <- function(x, v, w = 1)
{
test <- check.arguments(x, v, w)
m <- length(x)
n <- tail(dim(x[[1]]), 1)
r <- NCOL(v)
w <- if (length(w) == 1) { 1 / (m^2) 
} else { w <- (w + t(w)) / (2 * sum(w)) }
score <- canon.scores(x, v)
if (r == 1) {
	return(sum(w * cov(score)) * (n - 1) / n)
} 
out <- numeric(r)
for (l in 1:r) 
	out[l] <- sum(w * cov(score[,,l])) * (n - 1) / n
sum(out)
}

