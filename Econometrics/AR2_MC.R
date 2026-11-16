# AR2 estimators - Montecarlo Simulation
# Written by Jorge Paredes on Nov 2026. 

# [0] Preamble ----
library(dplyr); set.seed(1907)
nobs = 100 ; burn_in = 200

# [1] Data ----
## True parameters
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