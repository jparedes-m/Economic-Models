# Cobweb Model

Written by Jorge Paredes on May, 2022. 

## The model:

Suppose you have a market, were the supply function is dependent of the lagged price, whereas the demand is not.
Formally:

```math
Q_{d,t} = \alpha - \beta \cdot P_{t}\\
Q_{s,t} = -\gamma +\delta P_{t-1}
```

The equilibrium is defined when the supply quantity equals the demand quantity, formally:

```math
Q_{d,t} \equiv Q_{s,t}
```

Therefore the equilibrium becomes price dependant:
```math
\delta P_{t-1} + \beta P_{t} = \alpha +\gamma
```
Normalizing that as a difference equation:

```math
P_{t+1} + \frac{\delta}{\beta}P_{t} = \frac{\alpha + \gamma}{\beta}
```

Where the price dynamics are:
```math
P_{t} = \left(P_{0} - \frac{\alpha + \gamma}{\beta + \delta}\right)\cdot \left(-\frac{\delta}{\beta}\right)^{t}+\frac{\alpha + \gamma}{\beta + \delta}
```
where $`P_{0}`$ is the initial price.

However, the intertemporal equilibrium price is defined as:
```math
P^{*}=\frac{\alpha + \gamma}{\beta + \delta}
```
and the intertemporal equilibrium quantity is:

```math
Q^{*} = \alpha - \beta \frac{\alpha + \gamma}{\beta + \delta}
```
Beware that if $\delta > \beta$, the price trajectory will be explosive.
Furthermore, if$\delta = \beta$ the price dynamics will be oscillating ad infinitum, the only case where the price dynamics converge is when $\delta < \beta$.

As Alpha Chiang exemplifies:

![image](https://user-images.githubusercontent.com/103344273/170894992-2d449b55-2908-4b41-9faf-ff2c74d4823f.png)


## The code:

I made one function `lin_cobweb()`, where you put $`\alpha, \beta, \gamma, \delta`$, and the initial price $`P_{0}`$, furthermore, you can set the number of iterations (`iter`) to make, and the tolerance (`tol`), to get the iteration number. 

So if you input in R this:

```
lin_cobweb(alpha = 18, beta = 3, gamma = 3, delta = 2, p0 =2, iter=100)
```
this would be the results:

```
 Results of the Cobweb Model: 
 ------------------------------ 
 Intertemporal Equilibrium Price: 4.2 
 Intertemporal Equilibrium Quantity: 5.4 
This model converges at price and quantity at the 39 iteration 
```

![image](https://user-images.githubusercontent.com/103344273/170894885-8610e7e0-2edb-4de6-8258-db61a0559a18.png)
![image](https://user-images.githubusercontent.com/103344273/170894890-55d77fc1-5b46-4972-8f2e-c14690b0719b.png)




