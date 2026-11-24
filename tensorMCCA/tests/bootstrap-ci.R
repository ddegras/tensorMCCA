library(tensorMCCA)

# The expected values below use interval definitions, explicit triangle loops,
# and ordinary outer products, independently of the package helpers.
run.bootstrap.ci.tests <- function(ns, stop.on.failure = TRUE) {
    ci <- get("mcca.boot.ci", envir = ns)
    build <- get("build.ci", envir = ns)
    restore.cov <- get("vec2cov", envir = ns)
    restore.cor <- get("vec2cor", envir = ns)
    types <- c("basic", "normal", "normal-debiased", "percentile")
    results <- list()
    check <- function(name, expr) {
        failure <- tryCatch({force(expr); NULL}, error = conditionMessage)
        results[[name]] <<- if (is.null(failure)) "passed" else failure
    }
    equal <- function(x, y, tolerance = 1e-12) {
        stopifnot(length(x) == length(y),
                  isTRUE(all.equal(as.numeric(x), as.numeric(y),
                                   tolerance = tolerance)))
    }
    expected.ci <- function(stat, draws, level = .95) {
        stat <- as.numeric(stat)
        draws <- matrix(draws, nrow = length(stat))
        average <- apply(draws, 1L, mean)
        error <- apply(draws, 1L, sd) * sqrt((ncol(draws)-1)/ncol(draws))
        quantiles <- t(apply(draws, 1L, quantile,
                            probs = c((1-level)/2, (1+level)/2), names = FALSE))
        halfwidth <- qnorm((1+level)/2) * error
        list(bias = average-stat, se = error,
             basic = cbind(2*stat-quantiles[,2], 2*stat-quantiles[,1]),
             normal = cbind(stat-halfwidth, stat+halfwidth),
             "normal-debiased" = cbind(2*stat-average-halfwidth,
                                        2*stat-average+halfwidth),
             percentile = quantiles)
    }
    triangle <- function(values, size, diagonal = TRUE, fixed.diagonal = 1) {
        answer <- diag(fixed.diagonal, size)
        position <- 0L
        for (column in seq_len(size)) {
            for (row in seq_len(size)) {
                if (row > column || (diagonal && row == column)) {
                    position <- position+1L
                    answer[row,column] <- answer[column,row] <- values[position]
                }
            }
        }
        answer
    }

    check("four interval definitions and separate type requests", {
        expected <- expected.ci(10, matrix(1:4, 1L))
        together <- build(matrix(1:4, 1L), 10, .95, types)
        for (name in names(expected)) equal(together[[name]], expected[[name]])
        for (name in types) {
            separate <- build(matrix(1:4, 1L), 10, .95, name)
            equal(separate[[name]], expected[[name]])
        }
        equal(together$basic, c(16.075, 18.925))
        equal(together$bias, -7.5)
        equal(together$se, sqrt(1.25))
    })
    check("confidence levels below one half", {
        widths <- numeric()
        for (level in c(.1, .2, .5, .9, .95)) {
            actual <- build(matrix(1:4, 1L), 10, level, types)
            expected <- expected.ci(10, matrix(1:4, 1L), level)
            for (name in types) equal(actual[[name]], expected[[name]])
            widths <- c(widths, diff(actual$percentile[1,]))
        }
        stopifnot(all(diff(widths) > 0))
    })
    check("large offsets and constant bootstrap statistics", {
        offsets <- c(-2,-1,0,1,2,-2,-1,0,1,2,-1,1)
        actual <- build(matrix(1e8+offsets, 1L), 1e8, .95, types)
        equal(actual$se, sqrt(22/12))
        equal(actual$normal[1,]-1e8, qnorm(c(.025,.975))*sqrt(22/12), 1e-8)
        constant <- build(matrix(rep(7, 8), 2L), c(7,7), .95, types)
        equal(constant$se, c(0,0))
        for (name in types) equal(constant[[name]], rep(7,4))
    })
    check("one dimensional array estimates", {
        draws <- outer(1:6, c(.7,.9,1.1,1.3))
        actual <- build(draws, array(1:6, 6L), .95, types)
        expected <- expected.ci(1:6, draws)
        for (name in names(expected)) equal(actual[[name]], expected[[name]])
    })
    check("partial type names and a single deterministic replicate", {
        original <- list(noise.cov = list(10))
        object <- list(original = original, bootstrap = list(list(noise.cov = list(7))))
        actual <- ci(object, type = c("bas", "normal", "normal-deb", "perc"))$noise.cov[[1]]
        stopifnot(all(types %in% names(actual)))
        equal(actual$bias, -3); equal(actual$se, 0)
        equal(actual$basic, c(13,13)); equal(actual$percentile, c(7,7))
        equal(actual$normal, c(10,10)); equal(actual[["normal-debiased"]], c(13,13))
    })

    for (size in 1:8) for (columns in 1:2) {
        check(paste("covariance matrix reconstruction", size, columns), {
            draws <- matrix(seq_len(size*(size+1)/2*columns), ncol = columns)
            actual <- restore.cov(draws)
            stopifnot(identical(dim(actual), if (columns == 1L) c(size,size) else c(size,size,columns)))
            for (column in seq_len(columns)) {
                expected <- triangle(draws[,column], size)
                part <- if (columns == 1L) actual else actual[,,column]
                equal(part, expected)
                equal(part, t(part))
            }
        })
        check(paste("correlation matrix reconstruction", size, columns), {
            draws <- matrix(seq_len(size*(size-1)/2*columns)/100, ncol = columns)
            actual <- restore.cor(draws)
            stopifnot(identical(dim(actual), if (columns == 1L) c(size,size) else c(size,size,columns)))
            for (column in seq_len(columns)) {
                expected <- triangle(draws[,column], size, FALSE)
                part <- if (columns == 1L) actual else actual[,,column]
                equal(part, expected)
                equal(part, t(part))
            }
        })
    }

    for (components in 1:2) {
        check(paste("vector and tensor weight intervals", components), {
            weights <- matrix(vector("list", 3L*components), 3L, components)
            for (component in seq_len(components)) {
                weights[[1,component]] <- list(c(1,2,3)*component)
                weights[[2,component]] <- list(c(.5,1)*component, c(.2,.4,.6))
                weights[[3,component]] <- list(c(.1,.2), c(.3,.4)*component, c(.5,.6))
            }
            replicas <- lapply(c(-.3,-.1,.1,.3), function(delta) {
                current <- weights
                for (component in seq_len(components)) for (block in 1:3)
                    current[[block,component]][[1]] <- current[[block,component]][[1]]+delta
                list(v = current)
            })
            object <- list(original = list(v = weights), bootstrap = replicas)
            actual <- ci(object)
            expand <- function(factors) as.numeric(Reduce(function(a,b) outer(a,b), factors))
            for (component in seq_len(components)) for (block in 1:3) {
                stat <- expand(weights[[block,component]])
                draws <- vapply(replicas, function(replica) expand(replica$v[[block,component]]), stat)
                expected <- expected.ci(stat, draws)
                shape <- c(vapply(weights[[block,component]], length, 1L), 2L)
                for (name in types) {
                    stopifnot(identical(dim(actual$v[[block,component]][[name]]), shape))
                    equal(actual$v[[block,component]][[name]], expected[[name]])
                }
                equal(actual$v[[block,component]]$bias, expected$bias)
                equal(actual$v[[block,component]]$se, expected$se)
            }
        })
        check(paste("score and noise matrix intervals", components), {
            original <- list(block.score.cov = matrix(seq_len(10L*components),10L,components),
                             global.score.cov = seq_len(components*(components+1L)/2),
                             block.score.cor = matrix(seq_len(6L*components)/100,6L,components),
                             global.score.cor = seq_len(components*(components-1L)/2)/100,
                             noise.cov = list(3, 1:10), noise.cor = list(numeric(0), (1:6)/100))
            replicas <- lapply(c(-.03,-.01,.01,.03), function(delta) {
                lapply(original, function(target) {
                    if (is.list(target)) lapply(target, function(x) x+delta) else target+delta
                })
            })
            actual <- ci(list(original = original, bootstrap = replicas))
            stopifnot(setequal(names(actual), names(original)))
            check.matrix <- function(target, vectors, dimension, correlation) {
                if (length(target)) {
                    expected <- expected.ci(target, do.call(cbind, vectors))
                } else {
                    expected <- list(bias = numeric(0), se = numeric(0))
                    for (name in types) expected[[name]] <- matrix(numeric(0),0L,2L)
                }
                answer <- list()
                for (name in names(expected)) {
                    values <- as.matrix(expected[[name]])
                    diagonal <- if (correlation && name %in% c("bias","se")) 0 else 1
                    matrices <- lapply(seq_len(ncol(values)), function(column)
                        triangle(values[,column], dimension, !correlation, diagonal))
                    answer[[name]] <- if (length(matrices) == 1L) matrices[[1]] else
                        array(unlist(matrices), c(dimension,dimension,2L))
                }
                answer
            }
            for (kind in c("cov","cor")) {
                block.name <- paste0("block.score.",kind)
                global.name <- paste0("global.score.",kind)
                noise.name <- paste0("noise.",kind)
                for (component in seq_len(components)) {
                    expected <- check.matrix(original[[block.name]][,component],
                        lapply(replicas,function(x) x[[block.name]][,component]),4L,kind == "cor")
                    for (name in names(expected)) equal(actual[[block.name]][[component]][[name]],expected[[name]])
                }
                expected <- check.matrix(original[[global.name]],
                    lapply(replicas,function(x) x[[global.name]]),components,kind == "cor")
                for (name in names(expected)) equal(actual[[global.name]][[name]],expected[[name]])
                for (block in 1:2) {
                    expected <- check.matrix(original[[noise.name]][[block]],
                        lapply(replicas,function(x) x[[noise.name]][[block]]),c(1L,4L)[block],kind == "cor")
                    for (name in names(expected)) equal(actual[[noise.name]][[block]][[name]],expected[[name]])
                }
            }
        })
    }
    check("failed and malformed bootstrap replicates are identified", {
        original <- list(noise.cov = list(c(1,0,1)))
        valid <- list(original = original, bootstrap = list(original, original))
        expect.error <- function(object, pattern, ...) {
            error <- tryCatch({ci(object,...); NULL}, error = conditionMessage)
            stopifnot(!is.null(error), grepl(pattern,error))
        }
        bad <- valid; bad$bootstrap[[2]] <- simpleError("refit failed")
        expect.error(bad,"[Rr]eplicate 2")
        bad <- valid; bad$bootstrap[[2]] <- list()
        expect.error(bad,"[Rr]eplicate 2")
        bad <- valid; bad$bootstrap[[2]]$noise.cov[[1]] <- c(1,2)
        expect.error(bad,"[Rr]eplicate 2")
        bad <- valid; bad$bootstrap[[2]]$noise.cov[[1]][1] <- Inf
        expect.error(bad,"[Rr]eplicate 2")
        bad <- valid; bad$bootstrap[[2]]$noise.cov[[1]][1] <- 1+1i
        expect.error(bad,"[Rr]eplicate 2")
        bad <- valid; bad$original$noise.cov[[1]][1] <- 1+1i
        expect.error(bad,"Original")
        bad <- valid; bad$original$noise.cov[[1]][1] <- Inf
        expect.error(bad,"Original")
        bad <- valid; bad$bootstrap <- list()
        expect.error(bad,"replicate")
        for (level in list(NA_real_,Inf,0,1,c(.9,.95),1+1i)) expect.error(valid,"level",level = level)
        expect.error(valid,"type",type = c("normal","unknown"))
    })
    failures <- results[results != "passed"]
    cat(length(results), "bootstrap CI regression groups:",
        length(results)-length(failures), "passed,", length(failures), "failed\n")
    if (length(failures) && stop.on.failure)
        stop(paste(names(failures), unlist(failures), sep = ": ", collapse = "\n"))
    invisible(results)
}

run.bootstrap.ci.tests(asNamespace("tensorMCCA"))
