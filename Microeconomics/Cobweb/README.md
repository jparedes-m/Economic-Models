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

However, the intertemporal equilibrium price is defined as:
```math
\bar{P}=\frac{\alpha + \gamma}{\beta + \delta}
```
and the intertemporal equilibrium quantity is:

```math
\bar{Q} = \alpha - \beta \frac{\alpha + \gamma}{\beta + \delta}
```

