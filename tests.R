library(devtools)
path <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/2021/mCCA/code/tensorMCCA/tensorMCCA"
install.packages(path, type = "source", repos = NULL)
library(tensorMCCA)
load_all(path)


n <- 10
m <- 3
dimx <- list(c(10,n), c(10,15,n), c(10,5,5,n))
w <- matrix(1, m, m)
w2 <- w
diag(w2) <- 0

x <- vector("list", m)
for (i in 1:m) {
	x[[i]] <- array(runif(prod(dimx[[i]])), dimx[[i]])
	v[[i]] <- lapply(dimx[[i]][-length(dimx[[i]])], runif)
}
# list(matrix(runif(5*n), 5, n), array(runif(6*n), c(2, 3, n)), 
	# array(runif(8*n), c(2, 2, 2, n))) 
# v <- list(list(runif(5)), list(runif(2), runif(3)), 
	# list(runif(2), runif(2), runif(2)))
# control <- list(ortho = 1:3)

# options(error = NULL)
# debug(mcca.cov)
# undebug(mcca.cov)	
test <- mcca.cov(x, r = 10, init = "cca", verbose = TRUE)
test <- mcca.cov(x, r = 10, init = "svd", verbose = TRUE)
test <- mcca.cov(x, r = 10, init = "random", verbose = TRUE)
test <- mcca.cov(x, r = 10, init = v, verbose = TRUE)

test <- mcca.cov(x, r = 10, ortho = "canon", init = "cca", verbose = TRUE)

undebug(mcca.single.global.cov)
test <- mcca.cov(x, r = 10, w = w2, norm = "global",
	ortho = "canon.tnsr", control = control, verbose = TRUE)
test <- mcca.cov(x, r = 10, init = "svd", verbose = TRUE)
test <- mcca.cov(x, r = 10, init = "random", verbose = TRUE)
test <- mcca.cov(x, r = 10, init = v, verbose = TRUE)

# debug(tnsr.vec.prod)
# debug(objective.gradient)

objective.cov(x, v)
objective.gradient(x, v, w)

vv <- cbind(v, v)
dimnames(vv) <- NULL
test <- scale.v(vv, scale = "norm", cnstr == "block")
test <- scale.v(vv, scale = "norm", cnstr = "global")
test <- scale.v(v, x, scale = "var")
test <- mcca.init.random(x, r = 2, ortho.mode = 1:2)
i <- 2; k <- 2; crossprod(sapply(test[i,], "[[", k))

<x, v1> ortho <x, v2>