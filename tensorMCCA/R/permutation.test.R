permutation.test <- function(x, obj, nperm = 100, parallel = TRUE)
{
## Prepare MCCA
test <- check.arguments(x)
m <- length(x)
n <- tail(dim(x[[1]]), 1)
dimx <- lapply(x, dim)
d <- sapply(dimx, length) - 1
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
r <- obj$call.args$r
rhohat <- obj$objective
fields <- c("ortho.type", "ortho.cnstr", "scope")
for (name in fields) 
	assign(obj$call.args[[name]], name)
v <- obj$v

## Center data
for(i in 1:m) {
	d <- length(dim(x[[i]])) - 1
    xbar <- rowMeans(x[[i]], dims = d)
    x[[i]] <- x[[i]] - as.vector(xbar)
}

## Build deflated datasets
if (r > 1) {
	xdfl <- vector("list", m * (r - 1))
	dim(xdfl) <- c(m, r - 1)
	for (l in 2:r) {
		if (ortho.type == "weight" && scope == "block") {
			ortho.cnstr <- set.ortho.mat(v = v[,1:(l-1)], 
				modes = ortho.mode[, 1:(l-1), l])
			xfltd[,l-1] <- mapply(tnsr.mat.prod, x = x, 
				mat = ortho.cnstr$mat, modes = ortho.cnstr$modes, 
				SIMPLIFY = FALSE) 
			# if (is.list(obj$init)) { ... }
		} else if (ortho.type == "weight" && scope == "global") {
			xdfl[,l-1] <- deflate.x(x, v = v[, 1:(l-1)], 
				check.args = FALSE)
			# if (is.list(obj$init)) { ... }
		} else if (ortho.type == "score" && scope == "block") {
			xdfl[,l-1] <- deflate.x(x, 
				score = obj$block.score[,,1:(l-1)],
				scope = "block", check.args = FALSE)	
			# if (is.list(obj$init)) { ... }
		} else if (ortho.type == "score" && scope == "global") {
			xdfl[,l-1] <- deflate.x(x, 
				score = obj$global.score[,1:(l-1)],
				scope = "global", check.args = FALSE)			
			# if (is.list(obj$init)) { ... }
		}
	}	
}

## MAIN LOOP
if (parallel && require(foreach) && getDoParRegistered()) {
	rhostar <- foreach(j = 1:nperm, .combine = rbind) %dopar% {
		set.seed(j)
		perm <- sample(n)
		rho <- numeric(r) # canonical correlations
		for (l in 1:r) {
			xx <- if (l == 1) x else xdfl[,l-1]
			for (i in 1:m) {
				dim(xx[[i]]) <- c(prod(p[[i]]), n)
				xx[[i]] <- xx[[i]][, perm]
				dim(x[[i]]) <- dimx[[i]]
			}
			
			# ...
			# rho[l] <- mcca.perm()
		}
		rho
	}
} else {
	rhostar <- matrix(, nperm, r)
	for(j in 1:nperm) {
		set.seed(j)
		perm <- sample(n)
		for (l in 1:r) {
			xx <- if (l == 1) x else xdfl[,l-1]
			for (i in 1:m) {
				dim(xx[[i]]) <- c(prod(p[[i]]), n)
				xx[[i]] <- xx[[i]][, perm]
				dim(x[[i]]) <- dimx[[i]]
			}
			
			# ...
			# rhostar[j,l] <- mcca.perm()
		}
}

uncorrected.pval <- rowMeans(t(abs(rhostar)) > abs(rho))
corrected.pval <- uncorrected.pval



list(pvals = pvals, perm.objective = t(perm.objective))
}