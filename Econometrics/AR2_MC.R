# AR2 estimators - Montecarlo Simulation
# Written by Jorge Paredes on Nov 2026. 

# [0] Preamble ----
library(dplyr); library(ggplot2)
nobs = 100 ; burn_in = 200
set.seed(190799)
R = 1000

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
## Matrix preallocation
ols_est <- ols_se <- matrix(data = NA, nrow = R, ncol = 4, dimnames = list(NULL, c("b_0", "b_1", "b_2", "sigma2")))
mle_est <- mle_se <- matrix(data = NA, nrow = R, ncol = 4, dimnames = list(NULL, c("b_0", "b_1", "b_2", "sigma2")))
gibbs_est <- gibbs_se <- matrix(data = NA, nrow = R, ncol = 4, dimnames = list(NULL, c("b_0", "b_1", "b_2", "sigma2")))
times <- matrix(data = NA, nrow = R, ncol = 4, dimnames = list(NULL, c("replic", "ols", "mle", "gibbs")))

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

## MonteCarlo simulation
sizes <- c(25, 50, 75, 100, 250)
stats <- list(); times_all <- list()

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
}

