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

## Auxiliary quantities
r <- fit$call.args$r
scale <- fit$call.args$scale
if (scale == "block") {
	m <- length(x)
	dimx <- lapply(x, dim)
	p <- lapply(dimx, function(x) if (is.null(x)) 1 else x[-length(x)])
	sump <- sapply(p, sum) 
	cumsump <- cumsum(c(0, sump))
	idxH <- mapply(seq, from = cumsump[1:m] + 1, to = cumsump[-1], 
		SIMPLIFY = FALSE)
}

## Calculate Hessian matrices of Lagrange function
H <- hessian.internal(x, fit, lambda)		
		
## Project Hessian matrices on orthogonal of constraint gradients
projH <- vector("list", r)
# if (scale == "block") {
	# for (l in 1:r) {
		# UHU <- vector("list", m^2)
		# dim(UHU) <- c(m,m)
		# U <- vector("list",m)
		# for (i in 1:m) {
			# svdil <- svd(grad[[i,l]], nu = sump[i], nv = 0)
			# keep <- 1:sump[i]
			# rmv <- which(svdil$d > max(1e-14, svdil$d[1] * 1e-8))
			# if (length(rmv) > 0) keep <- keep[-rmv]
			# U[[i]] <- if (length(keep) > 0) {
				# svdil$u[,keep] } else numeric(sump[i])
		# }
		# for (i in 1:m) {
		# for (j in 1:i) {
			# UHU[[i,j]] <- crossprod(U[[i]], H[idxH[[i]], idxH[[j]], l] %*% U[[j]])
			# if (i > j) UHU[[j,i]] <- t(UHU[[i,j]])		
		# }}
		# for (j in 1:m)
			# UHU[[1,j]] <- do.call(rbind, UHU[,j])
		# nr <- nrow(UHU[[1]])
		# projH[[l]] <- matrix(unlist(UHU[1,]), nr, nr)	
	# }
# } else {
pp <- cumsump[m+1]
	# nr <- nrow(grad[[1]])
for (l in 1:r) {
	# svdl <- svd(grad[[l]], nu = nr, nv = 0)
	# keep <- 1:nr
	svdl <- svd(grad[[l]], nu = pp, nv = 0)
	keep <- 1:pp
	rmv <- which(svdl$d > max(1e-14, svdl$d[1] * 1e-8))
	if (length(rmv) > 0) keep <- keep[-rmv]
	U <- if (length(keep) > 0) {
		svdl$u[,keep] } else numeric(pp) # numeric(nr)	
	projH[[l]] <- crossprod(U, H[,,l] %*% U)
}	
# }	

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
				idxik <- idxH[[i]][[k]]
				idxjkk <- idxH[[j]][[kk]]
				Acc <- 0
				if (cst != 0) { # obj + scale (cor)
					tcp <- tcrossprod(tvprod[[i]][[k]], tvprod[[j]][[kk]])
					Acc <- cst * tcp
				}
				if (i > j) {
					H[idxik,idxjkk,l] <- Acc
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
				H[idxik,idxjkk,l] <- Acc
			}} # end k, kk loops			
			if (i > j) {
				idxi <- unlist(idxH[[i]])
				idxj <- unlist(idxH[[j]])
				H[idxj,idxi,l] <- t(H[idxi,idxj,l])	
			}	
		} # end j loop
	} # end i loop
} # end l loop

H
}



