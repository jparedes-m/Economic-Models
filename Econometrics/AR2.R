# AR2 estimators - why and how they differ
# Written by Jorge Paredes on Nov 2026. 

# [0] Preamble ----
library(dplyr); set.seed(1907)
nobs = 100 ; burn_in = 200

# [1] Data ----

## So an AR(2) has the form: y_t = b_0 + b_1 y_{t-1} + b_2 y_{t-2} + e_t, e_t ~ N(0, sigma^2)

## So we need that y_t follows a stationary process (so roots of the characteristic polinomial are out of unitary circle)
## The characteristic polinomial for this process is: 1 - b_1 z -b_2 z^2. Make it equal to zero to find the roots. 

## Let
b = matrix(c(
    9, # b_0
    0.9, #b_1 
    -0.5 #b_2
    ), ncol = 1)

sig2 = 1.2

## Stationarity Check [All have to be True]
b[2] + b[3] < 1; b[3] - b[2] < 1; abs(b[3]) <1

## Complex roots check
b[2]^2 + 4*b[3] < 0

## Generate the errors and the time series
e = matrix(rnorm(nobs+burn_in, 0, sqrt(sig2)), ncol = 1)

y = NULL; y[1] = b[1];  y[2] = b[1]
for(t in 3:(nobs+burn_in)){
    y[t] = b[1] + b[2] *y[t-1] + b[3]*y[t-2] + e[t]
}
y <- matrix(tail(y, nobs), ncol = 1)
e <- matrix(tail(e, nobs), ncol = 1)
rm(burn_in)

coefs <- c("b_0" = b[1], "b_1" = b[2], "b_2" = b[3], "sigma^2" = sig2)

# [2] OLS ----
Y = tail(y, -2)
x = matrix(c(rep(1, nobs), c(NA, y[-nobs]), c(NA, NA, y[1:(nobs-2)])),ncol = 3); X = tail(x,-2)

b_ols = solve(t(X) %*% X) %*% t(X) %*% Y

e_ols = Y - X %*% b_ols
sig2_ols = as.numeric((t(e_ols) %*% e_ols)/(nrow(X)-ncol(X)))
var_covar = sig2_ols * solve(t(X) %*% X)
se_ols = sqrt(diag(var_covar))
gamma_0_ols <- ((1-b_ols[3])*sig2_ols)/((1+b_ols[3])*( (1-b_ols[3])^2 -b_ols[2]^2))
R2_ols <- 1- sig2_ols/gamma_0_ols

coefs_ols <- c("b_0" = b_ols[1], "b_1" = b_ols[2], "b_2" = b_ols[3], "sigma^2" = sig2_ols)

# [3] Maximum Likelihood ----
## This is the exact maximum likelihood

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
opt = optim(c(b_ols, sig2_ols), negll, .Y = y, .X = X, method = "BFGS", hessian = TRUE)
coefs_mle <- c("b_0" = opt$par[1], "b_1" = opt$par[2], "b_2" = opt$par[3], "sigma^2" = opt$par[4])
se_mle <- sqrt(diag(solve(opt$hessian)))
gamma_0_hat <- ((1-opt$par[3])*opt$par[4])/((1+opt$par[3])*( (1-opt$par[3])^2 -opt$par[2]^2))
R2_mle <- 1- opt$par[4]/gamma_0_hat

# [4] Bayesian estimation ----

## Set up
G <- 11000; burn <- 1000
draws <- matrix(NA, G, 4, dimnames = list(NULL, c("b_0", "b_1", "b_2", "sigma^2")))

## Define the prior
lam <- 0.3
m <- c(0,1,0)
V <- diag(c(100, lam^2, (lam/2)^2))
a <- 0.01; d <- 0.01
Vi <- solve(V); XX <- t(X) %*% X; XY <- t(X) %*% Y 
sig2_g <- sig2_ols

# Gibbs loop
for (g in 1:G) {
    V1 <- solve(Vi + XX/sig2_g)                             # posterior var of b
    m1 <- V1 %*% (Vi %*% m + XY/sig2_g)                     # posterior mean of b
    repeat {                                                # redraw until stationary
        b_g <- m1 + t(chol(V1)) %*% rnorm(3)
        if (b_g[2] + b_g[3] < 1 & b_g[3] - b_g[2] < 1 & abs(b_g[3]) < 1) break
    }
    e_g <- Y - X %*% b_g                                    # residuals at b_g
    sig2_g <- 1/rgamma(1, shape = a + nrow(X)/2, rate = d + sum(e_g^2)/2)  # sig2 ~ IG
    draws[g, ] <- c(b_g, sig2_g)                            # store
}
draws <- draws[-(1:burn), ]                                 # drop burn-in

## Posterior summaries
coefs_bayes <- colMeans(draws)                              # posterior mean
se_bayes <- apply(draws, 2, sd)                             # posterior sd
gamma_0_g <- ((1-draws[,3])*draws[,4])/((1+draws[,3])*((1-draws[,3])^2 - draws[,2]^2))
R2_bayes <- mean(1 - draws[,4]/gamma_0_g)                   # implied R2, averaged over draws

# Comparison ----
plot(1:nobs, y, type = 'l', col = 'black')
lines(3:nobs, X%*%b_ols, col = 'red')
lines(3:nobs, X%*%coefs_mle[1:3], col = 'blue')
lines(3:nobs, X%*%coefs_bayes[1:3], col ='green')

fmt <- function(b, se) c(rbind(sprintf("%.3f", b), ifelse(is.na(se), "", sprintf("(%.3f)", se))))  # coef row + (se) row
tab <- cbind(Truth = fmt(coefs, NA), OLS = fmt(coefs_ols, c(se_ols, NA)), MLE = fmt(coefs_mle, se_mle), Bayes = fmt(coefs_bayes, se_bayes))
rownames(tab) <- c(rbind(names(coefs), ""))
print(tab, quote = FALSE)
                                                                 