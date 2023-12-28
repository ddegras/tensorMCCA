#######################################
# Main function to calculate Hessian 
# matrix of the Lagrange function
# and its projection on the orthogonal
# of the constraint gradients 
#######################################


hessian <- function(x, fit, grad = NULL)
{
## Constraint gradients and Lagrange multipliers 
if (is.null(grad)) grad <- kkt(x, fit)
lambda <- grad$multipliers
grad <- grad$grad.constraints
pp <- nrow(grad[[1]])

## Calculate Hessian matrices of Lagrange function
H <- hessian.internal(x, fit, lambda)		
		
## Project Hessian matrices on orthogonal of constraint gradients
r <- fit$call.args$r
projH <- vector("list", r)
for (l in 1:r) {
	svdl <- svd(grad[[l]], nu = pp, nv = 0)
	keep <- 1:pp
	rmv <- which(svdl$d > max(1e-14, svdl$d[1] * 1e-8))
	if (length(rmv) > 0) keep <- keep[-rmv]
	U <- if (length(keep) > 0) {
		svdl$u[,keep] } else numeric(pp) 	
	projH[[l]] <- crossprod(U, H[,,l] %*% U)
}	

list(H = H, projH = projH)
}






#################################
# Internal function to calculate 
# Hessian of Lagrange function 
#################################


hessian.internal <- function(x, fit, lambda)
{
m <- length(x)
dimx <- lapply(x, dim)
p <- lapply(dimx, function(x) if (is.null(x)) 1 else x[-length(x)])
d <- sapply(p, length) 
n <- if (is.null(dimx[[1]])) length(x[[1]]) else tail(dimx[[1]], 1)
r <- fit$call.args$r

## Canonical weights + scores and objective weights
objective <- fit$call.args$objective.type
ortho <- fit$call.args$ortho.type
scale <- fit$call.args$scale
v <- fit$v
score <- fit$block.score / n
ortho.score.block <- (ortho == "score" && scale == "block")
ortho.score.global <- (ortho == "score" && scale == "global")
if (ortho.score.global) global.score <- fit$global.score / n
w <- fit$call.args$w

## Block indices for output matrices
cumsump <- cumsum(c(0, unlist(p)))
idxH <- vector("list", m)
count <- 0
for (i in 1:m) {
for (k in 1:d[i]) {
	count <- count + 1
	idxH[[i]][[k]] <- (cumsump[count]+1):cumsump[count+1]
}}
pp <- cumsump[count+1]

## Auxiliary function
cpfun <- function(x,y) sum(x*y)

## Build Hessian matrices
H <- array(0, c(pp, pp, r))
for (l in 1:r) {
	tvprod <- vector("list", m)
	
	## First pass
	for (i in 1:m) {
		modes <- 1:d[i]
		for (k in 1:d[i])
			tvprod[[i]][[k]] <- tnsr.vec.prod(x[[i]], v[[i,l]][-k], modes[-k])
	}
	
	## Second pass
	if (scale == "global") lam <- lambda[[l]]
	for (i in 1:m) {	
		vil <- v[[i,l]]
		if (objective == "cov") {
			cp <- mapply(cpfun, x = unlist(v[i,1:l], FALSE), y = vil)
			dim(cp) <- c(d[i], l) }
		if (scale == "block") lam <- lambda[[l]][i,]
		# Right-multiplier vector for matrix terms x(i,t) x_(-(k,kk)) v(i,l)
		scorei <- score[,,l] %*% w[,i]
		if (objective == "cor")	
			scorei <- scorei - lam[l] * score[,i,l]
		if (ortho.score.block && l > 1)
			scorei <- scorei - as.matrix(score[,i, 1:(l-1)]) %*% lam[-l]
		if (ortho.score.global && l > 1)
			scorei <- scorei - as.matrix(global.score[,1:(l-1)]) %*% lam[-l]
		scorei <- list(scorei)			
		for (j in 1:i) {
			# Multiplicative constant for matrix terms
			# x(i,t) x_(-k) v(i,l) (x(j,t) x_(-kk) v(j,l))^T
			test <- (objective == "cor" && i == j)
			cst <- (w[i,j] - ifelse(test, lam[l], 0)) / n
			for (k in 1:d[i]) {
			for (kk in 1:d[j]) {
				idxrow <- idxH[[i]][[k]]
				idxcol <- idxH[[j]][[kk]]
				Acc <- 0
				if (cst != 0) { # obj + scale (cor)
					tcp <- tcrossprod(tvprod[[i]][[k]], tvprod[[j]][[kk]])
					Acc <- cst * tcp
				}
				if (i > j) {
					H[idxrow,idxcol,l] <- Acc
					next
				}
				if (k != kk) { # obj + scale (cor) + ortho (score)
					vs <- c(vil[-c(k,kk)], scorei)
					modes <- (1:(d[i]+1))[-c(k,kk)]
					tvprod2 <- tnsr.vec.prod(x[[i]], vs, modes)
					if (k > kk) tvprod2 <- t(tvprod2)
					Acc <- Acc + tvprod2	
				}
				if (k == kk && objective == "cov") { # scale (cov) 
						val <- prod(lam[l], cp[-k,l])
						Acc <- Acc - diag(val, p[[i]][k])
				} 
				if (k != kk && ortho == "weight") { # ortho (weight)
					for (ll in 1:l) {
						val <- prod(lam[ll], cp[-c(k,kk),ll])
						Acc <- Acc - tcrossprod(
							val * v[[i,ll]][[k]], v[[i,ll]][[kk]])
					}
				}
				H[idxrow,idxcol,l] <- Acc
			}} # end k, kk loops			
			if (i > j) {
				idxrow <- unlist(idxH[[i]])
				idxcol <- unlist(idxH[[j]])
				H[idxcol,idxrow,l] <- t(H[idxrow,idxcol,l])	
			}	
		} # end j loop
	} # end i loop
} # end l loop

H
}



