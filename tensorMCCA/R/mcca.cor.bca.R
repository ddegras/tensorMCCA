###################################################
# Internal function for MCCA with tensor data
# (maximize sum of correlations one block at a time)
###################################################

# This function is not meant to be called by the user
# It calculates a single set of canonical vectors
# The optimization is conducted one data block at a time
	
mcca.cor.bca <- function(x, v, w, ortho, sweep, maxit, 
	tol, verbose)
{
	
## Data dimensions
dimx <- lapply(x, dimfun)
d <- sapply(dimx, length) - 1L
m <- length(x)
n <- tail(dimx[[1]], 1) 
p <- mapply(head, dimx, d, SIMPLIFY = FALSE)
p[d == 0] <- 1L

objective <- numeric(maxit + 1L)
objective[1] <- objective.internal(x, v, w)
if (verbose) 
	cat("\nIteration", 0, "Objective", objective[1])
	
## The next two lines assume that the TMCCA starting point 
## is feasible, which is enforced in mcca.cor
vbest <- v
objective.best <- objective[1]

score <- matrix(0, n, m) # canonical scores <X_it, v_i>
if (sweep == "cyclical") idxi <- 1:m
if (!is.null(ortho) && !is.matrix(ortho)) 
	dim(ortho) <- c(m, 1)

for (it in 1:maxit) {
	if (sweep == "random") idxi <- sample(m)
	for (ii in 1:m) { 	
		i <- idxi[ii]	
		## Calculate the inner products <X_jt, v_j> 
		## After the first algorithm iteration (it = 1), in each 
		## iteration of the i loop, only the inner products associated 
		## with the previous value of i need being updated
		idxj <- if (it == 1 && ii == 1) { idxi[-1] 
			} else if (ii == 1) { lastidx } else idxi[ii-1]
		for (j in idxj) 
			score[, j] <- if (d[j] == 0) {
				x[[j]] * unlist(v[[j]])
			} else {
				tnsr.vec.prod(x[[j]], v[[j]], 1:d[j]) 
			}
		
		## Set up linear program
		a <- if (m == 1) { 
			array(0, p[[i]]) 
		} else if (d[i] == 0) {
			sum(x[[i]] * (score[, -i, drop = FALSE] %*% w[-i, i])) / n
		} else {
			tnsr.vec.prod(x = x[[i]], modes = d[i] + 1L,
				v = score[, -i, drop = FALSE] %*% (w[-i, i] / n))
		}
		
		## Update canonical vectors
		v[[i]] <- optim.block.cor(v = v[[i]], obj = a, 
			scale = x[[i]], maxit = maxit, tol = tol,
			ortho = ortho[i,])		
	}	
	lastidx <- idxi[m]

	## Calculate objective value
	objective[it + 1L] <- objective.internal(x, v, w)	
	if (objective[it + 1L] > objective.best) {
		vbest <- v
		objective.best <- objective[it + 1L]
	}

	if (verbose) 
		cat("\nIteration",it,"Objective",objective[it+1])
	
	## Check convergence 
	if (it > 1 && abs(objective[it + 1L] - objective[it]) <= 
	    	tol * max(1, objective[it])) break
}

list(v = vbest, score = canon.scores(x, vbest), 
	objective = objective.best, iters = it, 
	trace = objective[1:(it+1)])
}



