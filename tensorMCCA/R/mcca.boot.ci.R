#################################
# Main function for building 
# bootstrap confidence intervals
#################################

mcca.boot.ci <- function(object, level = 0.95, 
	type = c("basic", "normal", "normal-debiased", "percentile"))
{
if (!is.numeric(level) || is.complex(level) || length(level) != 1L || !is.finite(level) ||
    level <= 0 || level >= 1)
    stop("'level' must be a finite number strictly between 0 and 1.")
lookup <- c("basic", "normal", "normal-debiased", "percentile")
if (!is.character(type) || length(type) == 0L)
    stop("'type' must specify at least one supported interval type.")
type <- lookup[pmatch(type, lookup)]
if (anyNA(type)) stop("Please specify 'type' as 'basic', 'normal',",
	"'normal-debiased', 'percentile', or any combination thereof",
	"(partial match allowed).")
type <- unique(type)
target <- names(object$original)
out <- list()
nboot <- length(object$bootstrap)
if (nboot == 0L) stop("At least one bootstrap replicate is required.")
for (name in target) {
    values <- unlist(object$original[[name]])
    if (length(values) && (!is.numeric(values) || is.complex(values) || any(!is.finite(values))))
        stop("Original '", name, "' statistics must be finite real numbers.")
}

# Failed or incomplete refits must not silently disappear during unlisting.
same.shape <- function(a, b) {
    if (!identical(dim(a), dim(b)) || length(a) != length(b) ||
        is.list(a) != is.list(b)) return(FALSE)
    if (is.list(a)) return(all(mapply(same.shape, a, b)))
    TRUE
}
for (b in seq_len(nboot)) {
    replica <- object$bootstrap[[b]]
    if (inherits(replica, "error"))
        stop("Bootstrap replicate ", b, " failed: ", conditionMessage(replica))
    for (name in target) {
        if (!name %in% names(replica) ||
            !same.shape(object$original[[name]], replica[[name]]))
            stop("Bootstrap replicate ", b, " has a missing or incompatible '", name, "' target.")
        values <- unlist(replica[[name]])
        if (length(values) && (!is.numeric(values) || is.complex(values) || any(!is.finite(values))))
            stop("Bootstrap replicate ", b, " has non-finite or non-real '", name, "' statistics.")
    }
}

# Correlation diagonals are fixed: their bias and SE are zero, their bounds one.
restore.cor.ci <- function(ci) {
    result <- lapply(ci, vec2cor)
    for (name in intersect(c("bias", "se"), names(result))) diag(result[[name]]) <- 0
    result
}


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

if ("block.score.cor" %in% target) {
	r <- ncol(object$original$block.score.cor)
	out$block.score.cor <- out$global.score.cor <- vector("list", r)
	for (l in 1:r) {
		bootmat <- sapply(object$bootstrap, 
			function(obj) obj$block.score.cor[,l])
		stat <- object$original$block.score.cor[,l]
		ci <- build.ci(bootmat, stat, level, type)	
		out$block.score.cor[[l]] <- restore.cor.ci(ci)
	}	
	bootmat <- sapply(object$bootstrap, "[[", "global.score.cor")
	stat <- object$original$global.score.cor
	ci <- build.ci(bootmat, stat, level, type)
	out$global.score.cor <- restore.cor.ci(ci)
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
		out$noise.cor[[i]] <- restore.cor.ci(ci)
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
stat <- as.numeric(stat)
if (length(stat) == 0L) {
    out$bias <- out$se <- numeric(0)
    for (name in type) out[[name]] <- matrix(numeric(0), 0L, 2L)
    return(out)
}
if (!is.matrix(bootmat)) 
	bootmat <- matrix(bootmat, nrow = length(stat))
if (nrow(bootmat) != length(stat) || ncol(bootmat) == 0L ||
    any(!is.finite(bootmat)) || any(!is.finite(stat)))
    stop("Bootstrap statistics must be finite, with matching dimensions and at least one replicate.")
bootmean <- rowMeans(bootmat) 
out$bias <- bootmean - stat
# Retain denominator B while avoiding cancellation between large second moments.
out$se <- sqrt(rowMeans((bootmat - bootmean)^2))
alpha <- (1 - level) / 2
if (any(c("basic", "percentile") %in% type)) {
	percent.ci <- t(apply(bootmat, 1, quantile, 
		probs = c(alpha, 1 - alpha), names = FALSE))
	if ("percentile" %in% type) 
		out$percentile <- percent.ci
	if ("basic" %in% type) 
		out$basic <- 2 * stat - percent.ci[,2:1, drop = FALSE]
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
		out <- mat
	} else {
		grid <- matrix(seq_len(m^2), m, m)
		idxlo <- which(lower.tri(grid))
		idxup <- t(grid)[idxlo]
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
		out <- rep(1, n)
	} else {
		grid <- matrix(seq_len(m^2), m, m)
		idxlo <- which(lower.tri(grid))
		idxup <- t(grid)[idxlo]
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
