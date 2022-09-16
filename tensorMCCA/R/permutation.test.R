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
fields <- c("init.val", "objective.type", "ortho.type", 
	"ortho.cnstr", "scale", "scope")
for (name in fields) 
	assign(obj$call.args[[name]], name)
if ((!is.null(init.val)) && (!is.matrix(init.val))) 
	dim(init.val) <- c(m, 1)
v <- obj$v

## Center data
for(i in 1:m) {
    xbar <- rowMeans(x[[i]], dims = d[i])
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
			xdfl[,l-1] <- mapply(tnsr.mat.prod, x = x, 
				mat = ortho.cnstr$mat, modes = ortho.cnstr$modes, 
				SIMPLIFY = FALSE) 
			if ((!is.null(init.val)) && l <= ncol(init.val)) 
				init.val[,l] <- tnsr.rk1.mat.prod(init.val[,l], 
					mat = ortho.cnstr$mat, modes = ortho.cnstr$modes, 
					transpose.mat = FALSE)
		} else if (ortho.type == "weight" && scope == "global") {
			xdfl[,l-1] <- deflate.x(x, v = v[, 1:(l-1)], 
				check.args = FALSE)
		} else if (ortho.type == "score" && scope == "block") {
			xdfl[,l-1] <- deflate.x(x, 
				score = obj$block.score[,,1:(l-1)],
				scope = "block", check.args = FALSE)	
		} else if (ortho.type == "score" && scope == "global") {
			xdfl[,l-1] <- deflate.x(x, 
				score = obj$global.score[,1:(l-1)],
				scope = "global", check.args = FALSE)			
		}
	}	
}

## Scale initial canonical weights if provided
objective.cov <- (objective.type == "cov")
if (!is.null(init.val)) {
	for (l in 1:ncol(init.val)) 
		obj$call.args$init.val[,l] <- scale.v(init.val[,l], 
			type = if (objective.cov) "norm" else "var",
			scale = scale, x = if (objective.cov) NULL else x,
			check.args = FALSE)
}


## MAIN LOOP
if (parallel && require(foreach) && getDoParRegistered()) {
	rhoperm <- foreach(j = 1:nperm, .combine = rbind) %dopar% {
		set.seed(j)
		perm <- sample(n)
		rhopermj <- numeric(r) # canonical correlations
		for (l in 1:r) {
			xx <- if (l == 1) x else xdfl[,l-1]
			for (i in 1:m) {
				dim(xx[[i]]) <- c(prod(p[[i]]), n)
				xx[[i]] <- xx[[i]][, perm]
				dim(x[[i]]) <- dimx[[i]]
			}
			rhopermj[l] <- mcca.perm(xx, obj, l)
		}
		rhopermj
	}
} else {
	rhoperm <- matrix(, nperm, r)
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
			rhoperm[j, l] <- mcca.perm(xx, obj, l)
		}
	}
}

## p-values
pval.uncorrected <- rowMeans(t(abs(rhoperm)) > abs(obj$objective))
pval.corrected <- Reduce(max, uncorrected.pval, accumulate = TRUE)

list(pval.corrected = pval.corrected, 
	pval.uncorrected = pval.uncorrected,
	objective = obj$objective, 
	objective.perm = t(rhoperm))
}
