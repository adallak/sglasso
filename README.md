# Simplified Graphical Lasso
## Introduction

The `sglasso` is an R package that implements the Graphical Lasso algorithmwith the modifications proposed in Dallakyan and 
Pourahmadi (2024). The paper is available at [https://arxiv.org/abs/2403.12357].

The main function, `sglasso()`, takes a sample covariance matrix of the observations as input and returns an estimate of the precision matrix. 

## Installation

To install the latest version from Github, use

```s
library(devtools)
devtools::install_github("adallak/sglasso")
```

## Usage
### Data generation
```s
library(sglasso)
library(huge)
set.seed(12)
p <- 50
n <- 200
prob = 0.6
X <- huge.generator(n =n, d = p, prob = prob)$data
S <- cov(X)
```

### Estimating precision matrix with a fixed tuning parameter

The `sglasso()` function takes the following parameters:

* `S` - sample covariance matrix
* `init_theta` - initial value of precision matrix
* `lambda` - controls the tuning parameter
* `pen.diag` - whether to penalize the diagonal elements. True by default.
* `niter`   - maximum number of iterations.
* `inner_niter` - number of iterations for constrained box quadratic optimization
* `tolerance` - tolerance value for convergence.
* `thresh` - threshold value.

```s
theta = sglasso(S, lambda = 0.1)$theta
```

### With cross-validation
To select the threshold via cross-validation use the function `CVsglasso()`.
```s
cv = CVsglasso(X, lambda = 10^seq(-3, 3, 0.1))
```
