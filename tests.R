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
	v[[i]] <- lapply(dimx[[i]][-length(dimx[[i]])], runif)
}


## Make sure that the main functions run without errors
combs <- expand.grid(
	init = c("cca", "svd", "random", "custom"), 
	norm = c("block", "global"), 
	ortho = c("score", "canon.tnsr"), 
	sweep = c("cyclical", "random"),
	stringsAsFactors = FALSE)
ncombs <- nrow(combs)

# Test mcca.cov
for (i in 1:ncombs) {
init. <- if (combs[i,"init"] == "custom") {
	v } else combs[i,"init"] 
test <- mcca.cov(x, r = 5, w = w, init = init., 
	norm = combs[i, "norm"], 
	ortho = combs[i, "ortho"],
	sweep = combs[i, "sweep"], verbose = FALSE)
}

# Test mcca.cor
idx <- which(combs$norm == "block" | combs$ortho == "score")
for (i in idx) {
init. <- if (combs[i,"init"] == "custom") {
	v } else combs[i,"init"] 
test <- mcca.cor(x, r = 5, w = w, init = init., 
	norm = combs[i, "norm"], 
	ortho = combs[i, "ortho"],
	sweep = combs[i, "sweep"], verbose = TRUE)
}


