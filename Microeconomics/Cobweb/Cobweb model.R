# Cobweb Model

# Written by Jorge Paredes on May, 2022

# This script aims to recreate the Cobweb model presented in Ch. 17 on Fundamental Methods of Mathematical Economics from Alpha Chiang


# Preamble ----
library(ggplot2)
library(dplyr)
setwd("C:/Users/jparedesm/iCloudDrive/Desktop/Papers/Small Models/Cobweb")

# Linear CobWeb ----

lin_cobweb <- function(alpha, beta, gamma, delta, p0, iter = 100, tol = 6e-07, return_df=FALSE){
  # Set the vectors for supply, demand and price dynamics
  Qd <- numeric(iter)
  Qs <- numeric(iter)
  p <- numeric(iter)
  p[1] <- p0
  
  # Equilibrium
  p_eq <- (alpha+gamma)/(beta+delta)
  q_eq <- alpha-beta*p_eq
  
  # Create the functions for the demand and supply
  for(t in 1:iter){
  # Demand 
  Qd[t] = alpha - beta*t
  # Supply
  Qs[t] = -gamma + delta*t
  }  
  
  # Price dynamics
  for(t in 1:(iter-1)){
    # Price dynamics 
    p[t+1] = (alpha+gamma)/beta - (delta/beta)*p[t]
  }
  
  # Dataframe computation 
  prices <- data.frame(prices= p) %>% 
    mutate(Demand = alpha - beta*prices) %>% 
    mutate(Supply = -gamma + delta*lag(prices))
  
  
  #plotting
  plot(Qd, 0:(iter-1), type="l", col="red", xlim = c(0,2*q_eq), ylim=c(0,2*p_eq),
       ylab = "Price", xlab = "Quantity", main = "Market representation")
  lines(Qs,0:(iter-1), type="l", col="blue", lty=2)
  legend("topleft", 
         legend=c("Demand", "Supply"), 
         col=c("red", "blue"),
         border = "black",
         lty=1:2)
  
  

  convergence <- prices %>% 
    mutate(p_eq = p_eq,
           dif= abs(p_eq-prices),
           niter=seq(1:iter),
           when = ifelse(dif<=tol,"yes", "no")) %>% 
    filter(when=="yes") %>% 
    select(niter)
  convergence <- convergence[1,1]
  
  
  # Print the equilibrium price and quantity
  cat(paste(" Results of the Cobweb Model:","\n",
            "------------------------------","\n",
            "Intertemporal Equilibrium Price:", p_eq,"\n",
            "Intertemporal Equilibrium Quantity:", q_eq, "\n"))
      
  if(is.na(convergence)==T){
    cat(" The model has a steeper Supply. The price rather oscilates or explotes")
  }else{
    cat(paste("This model converges at price and quantity at the", convergence, "iteration","\n"))
  }
  
  
  
  readline(prompt="Press [enter] to continue to the next graph:")
  
  plot(0:(iter-1), p, type="l", main="Price Convergence", xlab = "Period", ylab="Price")

if(return_df==TRUE){
 return(prices) 
}
}

# Converge
lin_cobweb(alpha = 18, beta = 3, gamma = 3, delta = 2, p0 =2, iter=100)

# Oscila
lin_cobweb(alpha = 18, beta = 3, gamma = 3, delta = 3, p0 =2, iter=100)

# Diverge
lin_cobweb(alpha = 18, beta = 3, gamma = 3, delta = 4, p0 =2, iter=100)

