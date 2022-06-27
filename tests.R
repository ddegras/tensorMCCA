library(devtools)
path <- "~/Library/Mobile Documents/com~apple~CloudDocs/Research/2021/mCCA/code/tensorMCCA/tensorMCCA"
install.packages(path, type = "source", repos = NULL)
load_all(path)




x <- list(matrix(runif(50), 5, 10), array(runif(60), c(2, 3, 10)), 
	array(runif(80), c(2, 2, 2, 10))) 
	
options(error = recover)
debug(mcca.cov)
undebug(mcca.cov)	
test <- mcca.cov(x, r = 10, verbose = TRUE)

