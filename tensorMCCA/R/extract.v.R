## Extract single vectors from list of rank-1 tensors

# Inputs 
# v:	list of length m (datasets)
# 		Each component v[[i]] is a list of length d[i] (i=1:m)
#		Each component v[[i]][[k]] is a matrix of dimensions p[[i]][k] x r (k=1:d[i])
# orders: component order(s) to extract (integer vector, subset of 1:r)
# modes:  mode(s) to extract (list of length m)
#		Each component modes[[i]] has same length as orders 
#		Each 		

extract.v.internal <- function(v, orders, modes)
{
	m <- length(v)
	d <- sapply(v, length)
	for (i in 1:m)	{
		for (k in 1:d[i]) {
			idx <- which(modes[[i]] == k)
			if (length(idx) > 0) 
				v[[i]][[k]] <- vc[[i]][[k]][,orders[idx]]
		}
	}
	v
}



# extract.v <- function(v, idx, sort = c("order", "block"))
# {
	# test <- check.arguments(v = v)
	# r <- NCOL(v[[1]][[1]])
	# idx <- as.integer(idx)
	# if (!all(idx %in% (1:r)))
		# stop(paste("Values in 'idx' must be integers between 1",
			# "and the common number of columns of each matrix in 'v'"))
	# sort <- match(sort)
	# m <- length(v)
	# d <- sapply(v, length)
	# if (sort == "block") {
		# for (i in 1:m) 
			# for (k in 1:d[i])
				# v[[i]][[k]] <- v[[i]][[k]][, idx, drop = FALSE]
	# } else {
		# vout <- replicate(r, vector("list", r))
	# }
		
	# }
	# vout <- vector("list", ifelse(, m, length(idx)))
	
# }