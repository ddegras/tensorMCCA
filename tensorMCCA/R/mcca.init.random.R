


########################
# Random initialization 
# of canonical vectors 
########################


mcca.init.random <- function(x, r = 1L, objective = c("cov", "cor"), 
	ortho.mode = NULL)
{
test <- check.arguments(x)
objective <- match.arg(objective) 

m <- length(x) 
dimx <- lapply(x, dim) 
d <- sapply(dimx, length) - 1L
p <- mapply(head, dimx, d, SIMPLIFY = FALSE) 
r <- as.integer(r)
if (r == 1) ortho.mode <- NULL

## Set up orthogonality constraints
len.ortho <- length(ortho.mode)
if (len.ortho > 0) {
	if (is.numeric(ortho.mode)) 
		ortho.mode <- replicate(m, list(ortho.mode))
	if (is.list(ortho.mode) && len.ortho < m) 
		ortho.mode <- ortho.mode[rep_len(1:len.ortho, m)]
	if (identical(ortho.mode, "all")) 
		ortho.mode <- mapply(seq.int, from = rep(1L, m),
			to = d, SIMPLIFY = FALSE)
	for (i in 1:m)
		ortho.mode[[i]] <- intersect(ortho.mode[[i]], 1:d[i])
}

## Initialize canonical vectors randomly
v <- vector("list", m * r)
dim(v) <- c(m, r)
for (i in 1:m) {
	for (l in 1:r) {
		v[[i,l]] <- lapply(p[[i]], runif, min = -1, max = 1)
	}
}

## Orthogonalize tensors
if (len.ortho > 0) {
	for (i in 1:m) {
		for (k in ortho.mode[[i]]) {
			vk <- sapply(v[i,], "[[", k)
			vk <- qr.Q(qr(vk))
			for (l in 1:r) v[[i,l]][[k]] <- vk[, l]
		}
	}
}

## Scale and balance canonical vectors 
v <- scale.v(v, x = x, check.args = FALSE,
	scale = switch(objective, cov = "norm", cor = "var"))
if (r == 1) dim(v) <- NULL

v
}


