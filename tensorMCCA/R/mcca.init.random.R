mcca.init.random <- function(x, objective = c("cov", "cor"), 
	scale = c("block", "global"))
{
test <- check.arguments(x)
objective <- match.arg(objective) 
scale <- match.arg(scale)
m <- length(x) 
dimx <- lapply(x, dim) 
d <- sapply(dimx, length) - 1L
p <- mapply(head, dimx, d, SIMPLIFY = FALSE) 
v <- vector("list", m)
for (i in 1:m) 
	v[[i]] <- lapply(p[[i]], runif, min = -1, max = 1)
scale.v(v, x = x, check.args = FALSE, scale = scale,
	type = switch(objective, cov = "norm", cor = "var"))
}
