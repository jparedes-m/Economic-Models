# AR2 estimators - why and how they differ
# Written by Jorge Paredes on Nov 2026. 

# [0] Preamble ----
library(dplyr)
set.seed(1907)
nobs = 100 # number of observations

# [1] Data ----

## So an AR(2) has the form: y_t = b_0 + b_1 y_{t-1} + b_2 y_{t-2} + e_t
## e_t ~ N(0, sigma^2)

## So we need that y_t follows a stationary process (so roots of the characteristic polinomial are out of unitary circle)
## The characteristic polinomial for this process is: 1 - b_1 -_b2. Make it equal to zero to find the roots. 

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

## Generate y_t
e_t = matrix(rnorm(nobs, 0, sqrt(sig2)), ncol = 1)

y = NULL; y[1] = b[1]; y[2] = b[1]
y[3:nobs] = b[1] + b[2]*lag(y) + b[3]*lag(y, 2)
