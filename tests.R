library(devtools)
path <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/2021/mCCA/code/tensorMCCA/tensorMCCA"
install.packages(path, type = "source", repos = NULL)
library(tensorMCCA)
load_all(path)


## Generate artificial data
n <- 10
m <- 3
dimx <- list(c(10,n), c(10,15,n), c(10,5,5,n))
w <- matrix(1, m, m)
w2 <- w
diag(w2) <- 0
v <- vector("list", m)
x <- vector("list", m)
for (i in 1:m) {
	x[[i]] <- array(runif(prod(dimx[[i]])), dimx[[i]])
	xbar <- rowMeans(x[[i]], dims = length(dimx[[i]]) - 1)
	x[[i]] <- x[[i]] - as.vector(xbar)
	v[[i]] <- lapply(dimx[[i]][-length(dimx[[i]])], runif)
}

## Make sure that the main functions run without errors
args0 <- list(x = x, r = 5, w = w, verbose = FALSE)
combs <- expand.grid(
	init = c("cca", "svd", "random", "custom"), 
	scale = c("block", "global"), 
	ortho = c("score", "weight"), 
	sweep = c("cyclical", "random"),
	optim = c("bca", "grad.scale", "grad.rotate"),
	stringsAsFactors = FALSE)

# Test mcca.cov 
idx <- which(combs$optim == "bca")
for (i in idx) {
args <- c(args0, as.list(combs[i,]))
if (combs[i,"init"] == "custom") args$init <- v
test <- do.call(mcca.cov, args)
}

idx <- which(combs$optim == "grad.scale" & 
	combs$init == "cca")
for (i in idx) {
args <- c(args0, as.list(combs[i,]))
test <- do.call(mcca.cov, args)
}

idx <- which(combs$optim == "grad.rotate" & 
	combs$init == "cca" & combs$scale == "block")
for (i in idx) {
args <- c(args0, as.list(combs[i,]))
test <- do.call(mcca.cov, args)
}


# Test mcca.cor
idxr <- which(combs$optim == "bca" & 
	combs$scale == "block")
idxc <- which(colnames(combs) != "scale")
for (i in idxr) {
args <- c(args0, as.list(combs[i,idxc]))
if (combs[i,"init"] == "custom") args$init <- v
test <- do.call(mcca.cor, args)
}
idx <- which(combs$optim == "grad.scale" & 
	combs$init == "cca")
for (i in idxr) {
args <- c(args0, as.list(combs[i,idxc]))
test <- do.call(mcca.cov, args)
}

idxr <- which(combs$optim == "grad.rotate" & 
	combs$init == "cca" & combs$scale == "block")
for (i in idxr) {
args <- c(args0, as.list(combs[i,idxc]))
test <- do.call(mcca.cov, args)
}


