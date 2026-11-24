library(tensorMCCA)

## These tests use the installed package, base R, and deterministic data only.
ns <- asNamespace("tensorMCCA")
checks <- character()
failures <- character()
check <- function(name, code) {
    problem <- tryCatch({ force(code); NULL }, error = function(e) conditionMessage(e))
    checks <<- c(checks, name)
    if (is.null(problem)) {
        cat("PASS", name, "\n")
    } else {
        failures <<- c(failures, paste(name, problem, sep = ": "))
        cat("FAIL", name, ":", problem, "\n")
    }
}
equal <- function(a, b, tolerance = 1e-11)
    stopifnot(isTRUE(all.equal(a, b, tolerance = tolerance, check.attributes = TRUE)))

## Local copies permit deterministic capture/failure injection without changing
## the installed namespace or the package's fitting algorithms.
local.bootstrap <- function() {
    e <- new.env(parent = ns)
    for (name in c("mcca.boot", "mcca.boot.single", "calculate.score.cov",
                   "calculate.noise.cov", "calculate.residuals")) {
        f <- get(name, envir = ns)
        environment(f) <- e
        assign(name, f, envir = e)
    }
    e
}
tensor.vector <- function(row) Reduce(kronecker, rev(row))
fixture <- function(r = 1L, m = 3L, n = 40L) {
    set.seed(2101)
    signal <- matrix(rnorm(3*n), 3, n)
    x <- lapply(seq_len(m), function(i) {
        a <- signal + matrix(rnorm(3*n, sd = 0.3), 3, n)
        a - rowMeans(a)
    })
    v <- matrix(vector("list", m*r), m, r)
    for (i in seq_len(m)) for (j in seq_len(r))
        v[[i,j]] <- list(diag(3)[,j])
    score <- array(0, c(n, m, r))
    for (i in seq_len(m)) for (j in seq_len(r))
        score[,i,j] <- as.vector(crossprod(v[[i,j]][[1]], x[[i]]))
    global <- matrix(apply(score, c(1, 3), mean), n, r)
    list(x = x, fit = list(v = v, block.score = score, global.score = global,
        call.args = list(objective.type = "cov", r = r, w = 1-diag(m),
            scope = "block", ortho.type = "weight", optim = "bca",
            init.method = "svd", init.val = NULL, sweep = "cyclical",
            maxit = 40, tol = 1e-6, control = list())))
}

check("single-component score covariance and correlation", {
    z <- fixture()
    f <- get("calculate.score.cov", ns)
    both <- f(z$fit, c("cov", "cor"), unbiased = TRUE)
    expected <- matrix(0, 1, 1)
    for (i in 1:3) for (j in 1:3) if (i != j)
        expected <- expected + cov(matrix(z$fit$block.score[,i,], 40, 1),
                                   matrix(z$fit$block.score[,j,], 40, 1))/6
    equal(both$global.score.cov, expected)
    equal(both$global.score.cor, matrix(1, 1, 1))
    only.cor <- f(z$fit, "cor", unbiased = TRUE)
    equal(only.cor$global.score.cor, matrix(1, 1, 1))
    stopifnot(identical(dim(f(z$fit, "cov", TRUE, FALSE)$block.score.cov), c(6L, 1L)),
              length(f(z$fit, "cor", TRUE, FALSE)$global.score.cor) == 0L)
})

check("multiple-component score covariance preserves cross terms", {
    z <- fixture(2L)
    f <- get("calculate.score.cov", ns)
    expected <- matrix(0, 2, 2)
    for (i in 1:3) for (j in 1:3) if (i != j)
        expected <- expected + cov(matrix(z$fit$block.score[,i,], 40, 2),
                                   matrix(z$fit$block.score[,j,], 40, 2))/6
    both <- f(z$fit, c("cov", "cor"), TRUE)
    equal(both$global.score.cov, expected)
    equal(f(z$fit, "cov", FALSE)$global.score.cov, cov(z$fit$global.score))
    equal(f(z$fit, "cor", FALSE)$global.score.cor, cor(z$fit$global.score))
})

check("noise covariance lower-triangle serialization", {
    z <- fixture(2L)
    ## Keep the covariance subtraction strictly inside the positive-definite
    ## case; indefinite adjusted noise estimates are a separate model issue.
    z$fit$global.score <- z$fit$global.score / 3
    f <- get("calculate.noise.cov", ns)
    full <- f(z$x, z$fit, c("cov", "cor"), TRUE, TRUE)
    saved <- f(z$x, z$fit, c("cov", "cor"), TRUE, FALSE)
    for (i in 1:3) {
        expected <- cov(t(z$x[[i]]))
        for (j in 1:2) expected <- expected - var(z$fit$global.score[,j]) *
            tcrossprod(tensor.vector(z$fit$v[[i,j]]))
        equal(full$noise.cov[[i]], expected)
        equal(saved$noise.cov[[i]], expected[lower.tri(expected, TRUE)])
        equal(saved$noise.cor[[i]], cov2cor(expected)[lower.tri(expected)])
    }
})

check("residual bootstrap uses scope and retains one-component matrices", {
    f <- get("calculate.residuals", ns)
    for (r in 1:2) {
        z <- fixture(r)
        got <- f(z$x, z$fit, TRUE)
        for (i in 1:3) {
            expected <- z$x[[i]]
            expected[seq_len(r),] <- 0
            equal(got[[i]], expected)
            equal(f(z$x, z$fit, FALSE)[[i]], expected)
        }
    }
})

check("residual projection removes the loading span only", {
    z <- fixture(2L)
    z$fit$call.args$scope <- "global"
    for (i in 1:3) z$fit$v[[i,2]] <- z$fit$v[[i,1]]
    got <- get("calculate.residuals", ns)(z$x, z$fit, TRUE)
    for (i in 1:3) {
        expected <- z$x[[i]]
        expected[1,] <- 0
        equal(got[[i]], expected)
    }
})

check("one shared resampling index and explicit initial values", {
    for (kind in c("cov", "cor")) for (r in 1:2) {
        z <- fixture(r)
        z$fit$call.args$objective.type <- kind
        z$fit$call.args$init.method <- NULL
        z$fit$call.args$init.val <- z$fit$v
        e <- local.bootstrap()
        captured <- NULL
        capture <- function(x, init, ...) { captured <<- list(x = x, init = init); z$fit }
        e$mcca.cov <- e$mcca.cor <- capture
        set.seed(411)
        indices <- sample(40, 40, replace = TRUE)
        set.seed(411)
        e$mcca.boot.single(z$x, "data", "v", z$fit$call.args, TRUE, TRUE)
        equal(captured$init, z$fit$v)
        for (i in 1:3) equal(captured$x[[i]], z$x[[i]][,indices, drop = FALSE])
    }
})

check("named initializer and residual index schedule stay unchanged", {
    z <- fixture(2L)
    for (method in c("svd", "random", "cca")) {
        z$fit$call.args$init.method <- method
        residual <- lapply(z$x, function(a) { a[1:2,] <- 0; a })
        fitted <- Map("-", z$x, residual)
        e <- local.bootstrap()
        captured <- NULL
        e$mcca.cov <- function(x, init, ...) { captured <<- list(x = x, init = init); z$fit }
        set.seed(915)
        indices <- sample(40, 40, replace = TRUE)
        set.seed(915)
        e$mcca.boot.single(list(fitted = fitted, resid = residual,
            dimx = lapply(z$x, dim)), "residuals", "v", z$fit$call.args, TRUE, TRUE)
        stopifnot(identical(captured$init, method))
        for (i in 1:3) equal(captured$x[[i]], fitted[[i]] + residual[[i]][,indices, drop = FALSE])
    }
})

check("common tensor sign uses amplitudes across datasets", {
    z <- fixture(2L)
    for (i in 1:3) z$fit$v[[i,1]][[1]] <- c(c(10, 1, 1)[i], 0, 0)
    e <- local.bootstrap()
    e$mcca.boot.single <- function(...) {
        vv <- z$fit$v
        for (i in 2:3) vv[[i,1]][[1]] <- -vv[[i,1]][[1]]
        list(v = vv, global.score.cov = c(4, 0.5, 9), global.score.cor = 0.2)
    }
    expected <- e$mcca.boot.single()
    ## The unweighted cosine signs average to -1/3, but the true inner
    ## product in stacked weight space is 100 - 1 - 1 = 98 > 0.
    got <- e$mcca.boot(z$x, z$fit, target = "v", nboot = 1, parallel = FALSE)$bootstrap[[1]]
    equal(got, expected)
})

check("common sign also transforms cross-component score statistics", {
    z <- fixture(2L)
    e <- local.bootstrap()
    e$mcca.boot.single <- function(...) {
        vv <- z$fit$v
        for (i in 1:3) vv[[i,1]][[1]] <- -vv[[i,1]][[1]]
        list(v = vv, global.score.cov = c(4, 0.5, 9), global.score.cor = 0.2)
    }
    for (target in list(c("v", "score.cov", "score.cor"), c("score.cov", "score.cor"))) {
        got <- e$mcca.boot(z$x, z$fit, target = target, nboot = 1, parallel = FALSE,
                          control = list(unbiased.score.cov = FALSE))$bootstrap[[1]]
        equal(got$global.score.cov, c(4, -0.5, 9))
        equal(got$global.score.cor, -0.2)
        if ("v" %in% target) equal(got$v, z$fit$v) else stopifnot(is.null(got$v))
    }
})

check("pass and remove preserve their error-handling contracts", {
    z <- fixture()
    for (policy in c("pass", "remove")) {
        e <- local.bootstrap()
        count <- 0L
        e$mcca.boot.single <- function(...) {
            count <<- count + 1L
            if (count %% 2L == 1L) stop("deterministic bootstrap failure")
            list(v = z$fit$v)
        }
        got <- e$mcca.boot(z$x, z$fit, "v", nboot = 3, parallel = FALSE,
                          control = list(.errorhandling = policy))
        stopifnot(count == 3L)
        if (policy == "pass") {
            stopifnot(length(got$bootstrap) == 3L,
                      inherits(got$bootstrap[[1]], "error"),
                      inherits(got$bootstrap[[3]], "error"))
            equal(got$bootstrap[[2]]$v, z$fit$v)
        } else {
            stopifnot(length(got$bootstrap) == 1L)
            equal(got$bootstrap[[1]]$v, z$fit$v)
        }
    }
    e <- local.bootstrap()
    e$mcca.boot.single <- function(...) stop("all fail")
    got <- e$mcca.boot(z$x, z$fit, "v", nboot = 2, parallel = FALSE,
                      control = list(.errorhandling = "remove"))
    stopifnot(length(got$bootstrap) == 0L)
})

check("invalid counts and targets are rejected before resampling", {
    z <- fixture()
    e <- local.bootstrap()
    draws <- 0L
    e$mcca.boot.single <- function(...) { draws <<- draws+1L; list(v = z$fit$v) }
    for (bad in list(0, -1, 1.5, NA_real_, Inf, numeric(), c(1, 2), 1i)) {
        err <- tryCatch(e$mcca.boot(z$x, z$fit, "v", nboot = bad, parallel = FALSE),
                        error = identity)
        stopifnot(inherits(err, "error"))
    }
    err <- tryCatch(e$mcca.boot(z$x, z$fit, "not-a-target", nboot = 1, parallel = FALSE),
                    error = identity)
    stopifnot(inherits(err, "error"), draws == 0L)
})

check("public covariance bootstrap supports both resampling modes", {
    z <- fixture()
    for (r in 1:2) for (initial in c("named", "explicit")) {
        fit <- mcca.cov(z$x, r = r, ortho = "weight", init = "svd", maxit = 40)
        if (initial == "explicit")
            fit <- mcca.cov(z$x, r = r, ortho = "weight", init = fit$v, maxit = 40)
        for (resample in c("data", "residuals")) {
            set.seed(722)
            boot <- mcca.boot(z$x, fit, resample = resample, nboot = 2, parallel = FALSE,
                control = list(unbiased.score.cov = FALSE, unbiased.noise.cov = FALSE))
            stopifnot(length(boot$bootstrap) == 2L,
                      identical(dim(boot$bootstrap[[1]]$v), c(3L, as.integer(r))),
                      identical(names(boot$original), names(boot$bootstrap[[1]])))
        }
    }
})

cat(length(checks)-length(failures), "of", length(checks), "bootstrap groups passed\n")
if (length(failures)) stop(paste(failures, collapse = "\n"), call. = FALSE)
