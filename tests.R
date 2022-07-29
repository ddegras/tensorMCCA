library(devtools)
path <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/2021/mCCA/code/tensorMCCA/tensorMCCA"
install.packages(path, type = "source", repos = NULL)
library(tensorMCCA)
load_all(path)


n <- 10
w <- matrix(1/9, 3, 3)
w2 <- matrix(1/6, 3, 3)
diag(w2) <- 0
x <- list(matrix(runif(5*n), 5, n), array(runif(6*n), c(2, 3, n)), 
	array(runif(8*n), c(2, 2, 2, n))) 
v <- list(list(runif(5)), list(runif(2), runif(3)), 
	list(runif(2), runif(2), runif(2)))
control <- list(ortho = 1:3)

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