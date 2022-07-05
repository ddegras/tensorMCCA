library(devtools)
path <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/2021/mCCA/code/tensorMCCA/tensorMCCA"
install.packages(path, type = "source", repos = NULL)
library(tensorMCCA)
load_all(path)



w <- matrix(1/9, 3, 3)
x <- list(matrix(runif(50), 5, 10), array(runif(60), c(2, 3, 10)), 
	array(runif(80), c(2, 2, 2, 10))) 
v <- list(list(runif(5)), list(runif(2), runif(3)), 
	list(runif(2), runif(2), runif(2)))
	
options(error = recover)
debug(mcca.cov)
undebug(mcca.cov)	
test <- mcca.cov(x, r = 10, verbose = TRUE)

debug(tnsr.vec.prod)
debug(objective.gradient)

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