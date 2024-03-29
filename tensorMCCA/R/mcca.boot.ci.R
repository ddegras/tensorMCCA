#################################
# Main function for building 
# bootstrap confidence intervals
#################################

mcca.boot.ci <- function(object, level = 0.95, 
	type = c("basic", "normal", "normal-debiased", "percentile"))
{
lookup <- c("basic", "normal", "normal-debiased", "percentile")
type <- lookup[pmatch(type, lookup)]
if (all(is.na(type))) stop("Please specify 'type' as 'basic', 'normal',", 
	"'normal-debiased', 'percentile', or any combination thereof",
	"(partial match allowed).")
target <- names(object$original)
out <- list()
nboot <- length(object$bootstrap)


if ("v" %in% target) {
	v <- object$original$v
	m <- NROW(v)
	r <- NCOL(v)
	out$v <- vector("list", m * r)
	dim(out$v) <- c(m, r)
	d <- sapply(v[1:m], length) 
	for (i in 1:m) {
		pi <- sapply(v[[i]], length)
		for (l in 1:r) {
			if (d[i] == 1) {
				bootmat <- lapply(object$bootstrap, function(obj) obj$v[[i,l]][[1]])
				stat <- v[[i,l]][[1]]
			} else {
				bootmat <- lapply(object$bootstrap, 
					function(obj) outer.prod.nodim(obj$v[[i,l]]))
				stat <- outer.prod.nodim(v[[i,l]]) 
				# bootmat <- lapply(object$bootstrap, 
					# function(obj) Reduce(kronecker, rev(obj$v[[i,l]])))
				# stat <- Reduce(kronecker, rev(v[[i,l]])) 
			}
			bootmat <- matrix(unlist(bootmat), ncol = nboot)
			ci <- build.ci(bootmat, stat, level, type)
			for (name in type) {
				dim(ci[[name]]) <- c(pi, 2)
				dimn <- vector("list", d[i] + 1L)
				names(dimn) <- c(paste0("dim",1:d[i]), "bound")
				dimn$bound <- c("lwr", "upr")
				dimnames(ci[[name]]) <- dimn
			}
			out$v[[i,l]] <- ci
		}
	}
}

if ("block.score.cov" %in% target) {
	r <- ncol(object$original$block.score.cov)
	out$block.score.cov <- out$global.score.cov <- vector("list", r)
	for (l in 1:r) {
		bootmat <- sapply(object$bootstrap, function(obj) obj$block.score.cov[,l])
		if(!is.matrix(bootmat)) 
			bootmat <- matrix(bootmat, ncol = nboot)
		stat <- object$original$block.score.cov[,l]
		ci <- build.ci(bootmat, stat, level, type)
		out$block.score.cov[[l]] <- lapply(ci, vec2cov)	 
		names(out$block.score.cov[[l]]) <- names(ci)
	}	
	bootmat <- sapply(object$bootstrap, "[[", "global.score.cov")
	stat <- object$original$global.score.cov
	ci <- build.ci(bootmat, stat, level, type)
	out$global.score.cov <- lapply(ci, vec2cov)
	names(out$global.score.cov) <- names(ci)
}

if ("score.cor.block" %in% target) {
	r <- ncol(object$original$score.cor.block) 
	out$score.cor.block <- out$score.cor.global <- vector("list", r)
	for (l in 1:r) {
		bootmat <- sapply(object$bootstrap, 
			function(obj) obj$score.cor.block[,l])
		stat <- object$original$score.cor.block[,l]
		ci <- build.ci(bootmat, stat, level, type)	
		out$score.cor.block[[l]] <- lapply(ci, vec2cor)
		names(out$score.cor.block[[l]]) <- names(ci)
	}	
	# if (r > 1) {
		# out$score.cor.global <- replicate(length(type), 
			# c(lwr = 1, upr = 1), FALSE)		
		# names(out$score.cor.global) <- type
	# } else {
		bootmat <- sapply(object$bootstrap, "[[", "score.cor.global")
		stat <- object$original$score.cor.global
		ci <- build.ci(bootmat, stat, level, type)	
		out$score.cor.global <- lapply(ci, vec2cor)
	# }
}

if ("noise.cov" %in% target) {
	m <- length(object$original$noise.cov)
	out$noise.cov <- vector("list", m)
	for (i in 1:m) {
		bootmat <- sapply(object$bootstrap, 
			function(obj) obj$noise.cov[[i]])
		if(!is.matrix(bootmat)) 
			bootmat <- matrix(bootmat, ncol = nboot)
		stat <- object$original$noise.cov[[i]]
		ci <- build.ci(bootmat, stat, level, type)
		out$noise.cov[[i]] <- lapply(ci, vec2cov)	 
		names(out$noise.cov[[i]]) <- names(ci)
	}	
}

if ("noise.cor" %in% target) {
	m <- length(object$original$noise.cor)
	out$noise.cor <- vector("list", m)
	for (i in 1:m) {
		bootmat <- sapply(object$bootstrap, 
			function(obj) obj$noise.cor[[i]])
		if(!is.matrix(bootmat)) 
			bootmat <- matrix(bootmat, ncol = nboot)
		stat <- object$original$noise.cor[[i]]
		ci <- build.ci(bootmat, stat, level, type)
		out$noise.cor[[i]] <- lapply(ci, vec2cor)	 
		names(out$noise.cor[[i]]) <- names(ci)
	}	
}

out

}




###################################	
# Helper function to build CIs 
# and calculate SEs
###################################	

build.ci <- function(bootmat, stat, level, type)
{
out <- list()
if (!is.matrix(bootmat)) 
	bootmat <- matrix(bootmat, nrow = length(stat))
bootmean <- rowMeans(bootmat) 
out$bias <- bootmean - stat
out$se <- sqrt(rowMeans(bootmat^2) - bootmean^2) # standard error
test <- is.nan(out$se)
if (any(test)) out$se[test] <- 0
alpha <- min(level, 1 - level) / 2
if (any(c("basic", "percentile") %in% type)) {
	percent.ci <- t(apply(bootmat, 1, quantile, 
		probs = c(alpha, 1 - alpha), names = FALSE))
	if ("percentile" %in% type) 
		out$percentile <- percent.ci
	if ("basic" %in% type) 
		out$basic <- matrix(2 * bootmean - percent.ci[,2:1], ncol = 2)
}
if 	("normal" %in% type) {
	z <- qnorm(1 - alpha)
	halfwidth <- z * out$se
	out$normal <- cbind(stat - halfwidth, stat + halfwidth)  
}
if 	("normal-debiased" %in% type) {
	stat <- stat - out$bias 
	z <- qnorm(1 - alpha)
	halfwidth <- z * out$se
	out[["normal-debiased"]] <- cbind(stat - halfwidth, stat + halfwidth)  
}

out
}


##################################
# Helper function to convert vector 
# to covariance matrix
##################################

# Input: vectorized half-covariance matrix 
# or matrix of vectorized half-covariance matrices
# (one per column)

vec2cov <- function(vec) {
	mat <- as.matrix(vec)
	n <- ncol(mat)
	stopifnot(n %in% c(1,2))
	m <- as.integer((sqrt(8 * nrow(mat) + 1) - 1) / 2)
	if (m == 1) {
		out <- vec
	} else {
		idxlo <- which(lower.tri(matrix(0, m, m)))
		idxrow <- rep(1:(m-1), (m-1):1)
		idxcol <- rep(2:m, 1:(m-1)) 
		idxup <- (idxcol - 1L) * m + idxrow
		idxdiag <- cumsum(c(0L, m:2)) + 1L
		idxdiag2 <- seq(0, by = m, len = m) + (1:m)
		out <- matrix(0, m^2, n)
		out[idxlo,] <- mat[-idxdiag,]
		out[idxup,] <- out[idxlo,]
		out[idxdiag2,] <-  mat[idxdiag,]
	}
	if (n == 1) {
		dim(out) <- c(m, m)			
	} else {
		dim(out) <- c(m, m, n)	
		dimnames(out) <- list(data1 = 1:m, data2 = 1:m, bound = c("lwr", "upr"))			
	}
	out
}


####################################
# Helper function to convert vector 
# to correlation matrix
####################################

# Input: vectorized half-correlation matrix 
# or matrix of vectorized half-correlation matrices
# (one per column)

vec2cor <- function(vec) {
	mat <- as.matrix(vec)
	n <- max(1, ncol(mat))
	m <- as.integer((sqrt(8 * nrow(mat) + 1) + 1) / 2)
	if (m == 1) {
		out <- 1
	} else {
		idxlo <- which(lower.tri(matrix(0, m, m)))
		idxrow <- rep(1:(m-1), (m-1):1)
		idxcol <- rep(2:m, 1:(m-1)) 
		idxup <- (idxcol - 1L) * m + idxrow
		idxdiag <- seq(0, by = m, len = m) + (1:m)
		out <- matrix(1, m^2, n)
		out[idxlo,] <- mat
		out[idxup,] <- mat
		out[idxdiag,] <- 1
	}
	if (n == 1) {
		dim(out) <- c(m, m)
	} else {
		dim(out) <- c(m, m, n)
		dimnames(out) <- list(data1 = 1:m, data2 = 1:m, bound = c("lwr", "upr"))	
	}
	out
}
