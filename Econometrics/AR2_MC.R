# AR2 estimators - Montecarlo Simulation
# Written by Jorge Paredes on Nov 2026. 

# [0] Preamble ----
library(ggplot2)
nobs = 100 ; burn_in = 200
set.seed(190799)
R = 1000

fig_dir <- "Econometrics/figures"; res_dir <- "Econometrics/results"
dir.create(fig_dir, showWarnings = FALSE); dir.create(res_dir, showWarnings = FALSE)

# [1] Data ----
## True parameters
b = matrix(c(
    9, # b_0
    0.9, #b_1 
    -0.5 #b_2
    ), ncol = 1)
sig2 = 1.2
coefs <- c("b_0" = b[1], "b_1" = b[2], "b_2" = b[3], "sigma^2" = sig2)

## Stationarity Check [All have to be True]
b[2] + b[3] < 1; b[3] - b[2] < 1; abs(b[3]) <1
## Complex roots check
b[2]^2 + 4*b[3] < 0

## Simulate the data <- 
sim_ar2 <- function(.beta, .T, .burn_in = 200){
    e <- matrix(rnorm(.T + .burn_in, 0, sqrt(.beta[4])), ncol = 1)
    y <- NULL; y[1] <- .beta[1]; y[2] <- .beta[1]
    for(t in 3: (.T + .burn_in)){
        y[t] = .beta[1] + .beta[2] * y[t-1] + .beta[3] * y[t-2] + e[t] 
    }
    y <- matrix(tail(y, .T), ncol = 1)
    e <- matrix(tail(e, .T), ncol = 1)
    return(y)
}
y <- sim_ar2(.beta = coefs, .T = nobs)

# [2] Estimation functions ----
## [2.1] OLS
est_ols <- function(.Y){
    
    T <- length(.Y); .X <- matrix(c(rep(1, T), c(NA, .Y[-T]), c(NA, NA, .Y[1:(T-2)])), ncol = 3)
    .X <- tail(.X, -2); .Y <- tail(.Y, -2)
    
    beta_hat <- solve(t(.X) %*% .X) %*% t(.X) %*% .Y
    residuals_hat <- .Y - .X %*% beta_hat
    sig2_hat <- as.numeric((t(residuals_hat) %*% residuals_hat)/(nrow(.X)-ncol(.X)))
    VCV <- sig2_hat * solve(t(.X) %*% .X)
    se_hat <- sqrt(diag(VCV))
    se_sig2 <- sig2_hat * sqrt(2/(nrow(.X) - ncol(.X)))   # SE of sigma^2 under normal errors
    se_hat <- c(se_hat, se_sig2)

    gamma_0_hat <- ((1-beta_hat[3])*sig2_hat)/((1+beta_hat[3])*( (1-beta_hat[3])^2 -beta_hat[2]^2))
    R2_hat <- 1- sig2_hat/gamma_0_hat

    return(list("betas"= as.numeric(beta_hat), "sigma_2" = sig2_hat, "standard_errors" = se_hat, "R2" = R2_hat))
}

## [2.2] Exact MLE 
est_mle <- function(.Y){
    # compute ols to have starting values
    ols <- est_ols(.Y)

    # compute X matrix
    T <- length(.Y)
    .X <- matrix(c(rep(1, T), c(NA, .Y[-T]), c(NA, NA, .Y[1:(T-2)])), ncol = 3)
    .X <- tail(.X, -2)

    # compute negative log-likelihood
    negll <- function(theta, .Y, .X){
        # Unpack the components of the vector
        b <- theta[1:3]; sig2 <- theta[4]
        # Get the total number of periods (T)
        T <- length(.Y)
        # Make sure is a stationary process and has positive variance
        if(b[2] + b[3] >= 1| b[3] - b[2] >= 1| abs(b[3]) >= 1| sig2 <= 0) return(1e10)
        # compute the residuals from linear regression
        e = .Y[3:T] - .X %*% b
        # Compute mean, variance and autocovariance
        mu <- b[1]/(1-b[2]-b[3])
        gamma_0 <- ((1-b[3])*sig2)/((1+b[3])*( (1-b[3])^2 -b[2]^2))
        gamma_1 <- (b[2]*sig2)/((1+b[3])*( (1-b[3])^2 -b[2]^2))
        # Starting values of joint density
        Y_0 = .Y[1:2]; Sigma_0 = matrix(c(gamma_0, gamma_1, gamma_1, gamma_0), 2)
        mu_0 = matrix(c(mu, mu),ncol=1)
        # log-likelihood
        ll = -log(2*pi) - 0.5 * log(det(Sigma_0)) - 0.5*t(Y_0 - mu_0) %*% solve(Sigma_0)%*% (Y_0 - mu_0) - (T-2)/2* log(2* pi * sig2) - (t(e)%*%e)/(2*sig2) 
        # get it back
        return(-as.numeric(ll))
    }
    opt <- optim(c(ols$betas, ols$sigma_2), negll, .Y = .Y, .X = .X, method = "BFGS", hessian = TRUE)
    
    # Convergence check
    if (opt$convergence != 0) return(list("betas" = rep(NA, 3), "sigma_2" = NA, "standard_errors" = rep(NA, 4), "R2" = NA))

    # unpack and tidy up for retrieval
    beta_hat <- c("b_0" = opt$par[1], "b_1" = opt$par[2], "b_2" = opt$par[3])
    se_hat <- sqrt(diag(solve(opt$hessian)))
    sig2_hat <- opt$par[4]
    gamma_0_hat <- ((1-beta_hat[3])*sig2_hat)/((1+beta_hat[3])*( (1-beta_hat[3])^2 -beta_hat[2]^2))
    R2_hat <- 1- sig2_hat/gamma_0_hat
    return(list("betas" = as.numeric(beta_hat), "sigma_2" = sig2_hat, "standard_errors" = se_hat, "R2" = R2_hat))
}

## [2.3] Gibbs Sampling 
est_gibbs <- function(.Y,.draws = 10000, .burn_in = 1000, .lambda, .m, .V=NULL, .a, .d){
    ols <- est_ols(.Y)

    T <- length(.Y); .X <- matrix(c(rep(1, T), c(NA, .Y[-T]), c(NA, NA, .Y[1:(T-2)])), ncol = 3)
    .X <- tail(.X, -2); .Y <- tail(.Y, -2)
    ## Set up
    n_iters <- .draws + .burn_in
    draws <- matrix(NA, n_iters, 4, dimnames = list(NULL, c("b_0", "b_1", "b_2", "sigma^2")))
    if(is.null(.V)) .V <- diag(c(100, .lambda^2, (.lambda/2)^2))
    Vi <- solve(.V); XX <- t(.X) %*% .X; XY <- t(.X) %*% .Y
    sig2_g <- ols$sigma_2

    # Gibbs loop
    for(i in 1:n_iters){
        V1 <- solve(Vi + XX/sig2_g)
        m1 <- V1 %*% (Vi %*% .m + XY/sig2_g)
        repeat{
            beta_hat <- m1 + t(chol(V1)) %*% rnorm(3)
            if (beta_hat[2] + beta_hat[3] < 1 & beta_hat[3] - beta_hat[2] < 1 & abs(beta_hat[3]) < 1) break
        }
        residuals_g <- .Y - .X %*% beta_hat
        sig2_g <- 1/rgamma(1, shape = .a + nrow(.X)/2, rate = .d + sum(residuals_g^2)/2)  # sig2 ~ IG
        draws[i, ] <- c(beta_hat, sig2_g)                                                 # store
    }
    draws <- draws[-(1:.burn_in), ]
    gamma_0_g <- ((1-draws[,3])*draws[,4])/((1+draws[,3])*((1-draws[,3])^2 - draws[,2]^2))
    return(list("betas" = as.numeric(colMeans(draws)[1:3]), "sigma_2" = mean(draws[,4]), "standard_errors" = apply(draws, 2, sd), "R2" = mean(1 - draws[,4]/gamma_0_g)))
}

# [3] Monte Carlo Simulation ----
## Comparison statistics
MC_comp <- function(.est, .se, .coefs_truth = coefs){
    truth <- matrix(rep(.coefs_truth, nrow(.est)), ncol = 4, byrow = TRUE, dimnames = list(NULL, c("b_0", "b_1", "b_2", "sigma2")))
    errors <- .est - truth
    stats <- rbind(
        bias = colMeans(errors, na.rm = TRUE), 
        sd = apply(.est, 2, sd, na.rm = TRUE), 
        rmse = sqrt(colMeans(errors^2, na.rm =TRUE)),
        avg_se = colMeans(.se, na.rm = TRUE), 
        coverage = colMeans(abs(errors) <= 1.96 * .se, na.rm = TRUE)
        )
    colnames(stats) <- colnames(truth)
    return(stats)
}

## Matrix preallocation
ols_est <- ols_se <- matrix(data = NA, nrow = R, ncol = 4, dimnames = list(NULL, c("b_0", "b_1", "b_2", "sigma2")))
mle_est <- mle_se <- matrix(data = NA, nrow = R, ncol = 4, dimnames = list(NULL, c("b_0", "b_1", "b_2", "sigma2")))
gibbs_est <- gibbs_se <- matrix(data = NA, nrow = R, ncol = 4, dimnames = list(NULL, c("b_0", "b_1", "b_2", "sigma2")))
times <- matrix(data = NA, nrow = R, ncol = 4, dimnames = list(NULL, c("replic", "ols", "mle", "gibbs")))
sizes <- c(25, 50, 75, 100, 250)
stats <- list(); times_all <- list(); est_all <- list()

for(n in sizes){
    for(r in 1:R){
        # generate the data
        y_r <- sim_ar2(coefs, n)
        # estimation
        ## OLS
        t_ols <- system.time(ols <- est_ols(.Y = y_r))["elapsed"]
        ols_est[r,] <- c(ols$betas, ols$sigma_2); ols_se[r,] <- ols$standard_errors
        ## MLE
        t_mle <- system.time(mle <- est_mle(.Y = y_r))["elapsed"]
        mle_est[r,] <- c(mle$betas, mle$sigma_2); mle_se[r,] <- mle$standard_errors
        ## Gibbs
        t_gibbs <- system.time(gibbs <- est_gibbs(.Y = y_r, .lambda = 0.3, .m = c(0,1,0), .a = 0.01, .d = 0.01))["elapsed"]
        gibbs_est[r,] <- c(gibbs$betas, gibbs$sigma_2); gibbs_se[r,] <- gibbs$standard_errors
        ## Log
        times[r,] <- c(r, t_ols, t_mle, t_gibbs)
        cat(sprintf("n = %d | Replic: %d | OLS done in %.2f s | MLE done in %.2f s | Gibbs-Sample done in %.2f s | Total %.2f s |\n", n, r, t_ols, t_mle, t_gibbs, t_ols + t_mle + t_gibbs)); flush.console()
    }
    # Store results
    stats[[paste0("n_", n)]] <- list(
        OLS = MC_comp(.est = ols_est, .se = ols_se), 
        MLE = MC_comp(.est = mle_est, .se = mle_se),
        GIBBS = MC_comp(.est = gibbs_est, .se = gibbs_se)
    )
    times_all[[paste0("n_",n)]] <- times
    est_all[[paste0("n_", n)]] <- list(OLS = ols_est, MLE = mle_est, GIBBS = gibbs_est)
}
## Save Monte Carlo output
saveRDS(list(est = est_all, stats = stats, times = times_all, truth = coefs, R = R, sizes = sizes), file.path(res_dir, "AR2_MC_results.rds"))

# [4] Monte Carlo results ----

library(ggplot2)
cols <- c(OLS = "#A52A2A", MLE = "#104E8B", GIBBS = "#556B2F")                    # one color per method, names as in est_all / stats
pars <- colnames(est_all[[1]]$OLS); truth <- setNames(coefs, pars)
plab <- c(b_0 = "b[0]", b_1 = "b[1]", b_2 = "b[2]", sigma2 = "sigma^2")            # plotmath labels

## Sampling distributions
xr <- sapply(pars, function(p) quantile(unlist(lapply(est_all, function(e) sapply(e, function(m) m[, p]))), c(0.005, 0.995), na.rm = TRUE))  # common x range per parameter
grid <- expand.grid(n = names(est_all), parameter = pars, method = names(cols), stringsAsFactors = FALSE)
dens <- do.call(rbind, Map(function(n, p, m) {
    x <- na.omit(est_all[[n]][[m]][, p]); d <- density(x, from = xr[1, p], to = xr[2, p])
    data.frame(n = paste("n =", sub("n_", "", n)), parameter = plab[[p]], method = m, x = d$x, y = d$y, mean = mean(x), truth = truth[[p]])
}, grid$n, grid$parameter, grid$method))
dens$y <- dens$y / ave(dens$y, dens$n, dens$parameter, FUN = max)                 # each panel peaks at 1, so rows can share the y axis
dens$n <- factor(dens$n, unique(dens$n)); dens$parameter <- factor(dens$parameter, plab)

ggplot(dens, aes(x, y, colour = method)) +
    geom_line(linewidth = 0.8) +
    geom_vline(data = unique(dens[c("n", "parameter", "truth")]), aes(xintercept = truth), linetype = "dashed", linewidth = 1) +   # true value
    geom_segment(data = unique(dens[c("n", "parameter", "method", "mean")]), aes(x = mean, xend = mean, y = 0, yend = 0.12), linewidth = 1.2) +  # mean ticks
    facet_grid(n ~ parameter, scales = "free_x", switch = "y", labeller = labeller(parameter = label_parsed)) + 
    scale_colour_manual(values = cols, breaks = names(cols), labels = c("OLS", "MLE", "Gibbs")) +
    labs(x = NULL, y = NULL, colour = NULL, caption = "Dashed line: true value. Ticks at base: mean of each estimator.") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", axis.text.y = element_blank(), panel.grid = element_blank(), panel.border = element_rect(fill = NA, colour = "grey40"), strip.text = element_text(size = 11))
ggsave(file.path(fig_dir, "fig_densities.png"), width = 10, height = 10, dpi = 160, bg = "white")

## Performance
perf <- do.call(rbind, lapply(names(stats), function(n) do.call(rbind, lapply(names(stats[[n]]), function(m) {
    s <- stats[[n]][[m]]; tr <- abs(unname(coefs))
    data.frame(n = sub("n_", "", n), method = m, parameter = colnames(s),
               metric = rep(c("RMSPE (%)", "Bias (%)", "Coverage (%)"), each = ncol(s)),
               value = c(100 * s["rmse", ] / tr, 100 * s["bias", ] / tr, 100 * s["coverage", ]))
}))))
perf <- within(perf, {
    n <- factor(n, sub("n_", "", names(stats)))
    method <- factor(method, names(cols), c("OLS", "MLE", "Gibbs"))
    metric <- factor(metric, c("RMSPE (%)", "Bias (%)", "Coverage (%)"))
    parameter <- factor(parameter, pars, c("b[0]", "b[1]", "b[2]", "sigma^2"))
})

ggplot(perf, aes(n, value, fill = method)) +
    geom_col(position = position_dodge(0.8), width = 0.75) +
    geom_hline(data = data.frame(metric = factor(c("Bias (%)", "Coverage (%)"), levels(perf$metric)), y = c(0, 95)),
               aes(yintercept = y), linetype = "dashed", linewidth = 0.6, colour = "grey20") +       # zero bias, nominal 95%
    facet_grid(metric ~ parameter, scales = "free_y", switch = "y", labeller = labeller(parameter = label_parsed)) +
    scale_fill_manual(values = unname(cols)) +
    labs(x = "Sample size (n)", y = NULL, fill = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", strip.placement = "outside", panel.grid.major.x = element_blank(), panel.grid.minor = element_blank())
ggsave(file.path(fig_dir, "fig_performance.png"), width = 11, height = 7.5, dpi = 160, bg = "white")