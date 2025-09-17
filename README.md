PSor Package Vignette
================

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

The goal of **PSor** is to estimate principal causal effects under
principal stratification using a margin-free, variational-independent
odds ratio sensitivity parameter, allowing analysis when monotonicity
may not hold. The framework unifies the monotonicity assumption with the
counterfactual intermediate independence assumption. The framework also
assumes the mean principal ignorability. The package accompanies the
paper \`\`Semiparametric Principal Stratification Analysis Beyond
Monotonicity’’ and provides point estimates, standard errors, and
confidence intervals for both the conditionally doubly robust (CDR) and
debiased machine learning (DML) estimators.

## Installation

You can install the development version of PSor from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("deckardt98/PSor")
```

## Example

This example demonstrates how to use `PSor.fit` to estimate principal
causal effects with simulated data from our manuscript. To summarize,
the data will include a binary treatment `Z`, a binary intermediate
outcome `D`, a continuous final outcome `Y`, and baseline covariates,
**X**. Let `Y(z)` and `D(z)` respectively denote the potential
final outcome and potential intermediate outcome under treatment value
`Z=z`. Under principal stratification, the estimand of interest is the
principal causal effect:
$$\mu_{d_0d_1}=E[Y(1)-Y(0)|D(0)=d_0,D(1)=d_1].$$ For example, under a
noncompliance setup where $D$ denotes the actual treatment received, the
principal strata variable $G=(D(0),D(1))$ can be interpreted as follows:
$G=11$ represents always-takers, $G=01$ represents compliers, $G=00$
represents never-takers, and $G=10$ represents defiers. The primary
estimand of interest is the complier average causal effect (CACE) within
the stratum $G=01$. Next, we use a simulated dataset with four
covariates to illustrate an example application of the package.

### 1. Load Package and Generate Data

First, we load the necessary packages and define a function to generate
a simulated dataset. This simulated data will include a binary treatment
`Z`, a binary intermediate outcome `D`, a continuous final outcome `Y`,
and four covariates, `X1`-`X4`.

``` r
library(truncnorm)

expit <- function(x){return(exp(x)/(1+exp(x)))}
simu_full_data <- function(n, seed=20250917, theta){
  # Input:
  # n: sample size
  # theta: odds ratio sensitivity parameter; theta = 1 assumes independence and theta = Inf assumes monotonicity
  set.seed(seed)
  Z <- D <- c()
  # Simulate covariates
  X1 <- rtruncnorm(n, a=-20, b=20, mean = 0, sd = 1)
  X2 <- rtruncnorm(n, a=-20, b=20, mean = 0, sd = 1)
  X3 <- rtruncnorm(n, a=-20, b=20, mean = 0, sd = 1)
  X4 <- rbinom(n, size = 1, prob = 0.5)
  if (theta==Inf){
    X1 <- abs(X1)
    X2 <- abs(X2)
    X3 <- abs(X3)
    probZ <- expit(0.1*(X1+X2+X3)+0.5*X4)
    Z <- rbinom(n, size = 1, prob = probZ)
    
    # Simulate G=(D(0),D(1))
    probD1 <- expit(1.2*X4)
    probD0 <- expit(-0.4-0.2*X1-0.2*X2-0.2*X3-0.2*X4)
    prob11 <- probD0
    prob01 <- probD1 - probD0
    prob00 <- 1 - probD1
    prob_matrix <- cbind(prob00, prob01, prob11)
    # Simulate G from the categorical distribution
    G <- apply(prob_matrix, 1, function(p) sample(0:2, size = 1, prob = p))
    # Simulate D(1)
    D1 <- as.numeric(I(G!=0))
    # Simulate D(0)
    D0 <- as.numeric(I(G==2))
    # Compute observed intermediate outcome
    D <- D1*Z+(1-Z)*D0
  } else if (theta==1) {
    probZ <- expit(0.1*(X1+X2+X3)+0.5*X4)
    Z <- rbinom(n, size = 1, prob = probZ)
    # Simulate (D(0),D(1))
    probD1 <- expit(0.3*X1+0.4*X2+0.3*X3+0.5*X4)
    probD0 <- expit(0.4*X1+0.3*X2+0.4*X3+0.5*X4)
    prob11 <- probD0*probD1
    prob10 <- probD0 - prob11
    prob01 <- probD1 - prob11
    prob00 <- 1-prob11-prob10-prob01
    prob_matrix <- cbind(prob10, prob00, prob01, prob11)
    jointD <- apply(prob_matrix, 1, function(p) sample(1:4, size = 1, prob = p))
    # Simulate D(0)
    D0 <- as.numeric(I(jointD==1|jointD==4))
    # Simulate D(1)
    D1 <- as.numeric(I(jointD==3|jointD==4))
    # Compute observed intermediate outcome
    D <- D1*Z+(1-Z)*D0
  } else {
    probZ <- expit(0.1*(X1+X2+X3)+0.5*X4)
    Z <- rbinom(n, size = 1, prob = probZ)
    # Simulate (D(0),D(1))
    probD1 <- expit(0.3*X1+0.4*X2+0.3*X3+0.5*X4)
    probD0 <- expit(0.4*X1+0.3*X2+0.4*X3+0.5*X4)
    deltaX <- (1+(theta-1)*(probD1+probD0))^2-4*theta*(theta-1)*probD1*probD0
    prob11 <- (1+(theta-1)*(probD1+probD0)-sqrt(deltaX))/2/(theta-1)
    prob10 <- probD0 - prob11
    prob01 <- probD1 - prob11
    prob00 <- 1-prob11-prob10-prob01
    prob_matrix <- cbind(prob10, prob00, prob01, prob11)
    jointD <- apply(prob_matrix, 1, function(p) sample(1:4, size = 1, prob = p))
    # Simulate D(0)
    D0 <- as.numeric(I(jointD==1|jointD==4))
    # Simulate D(1)
    D1 <- as.numeric(I(jointD==3|jointD==4))
    # Compute observed intermediate outcome
    D <- D1*Z+(1-Z)*D0
  }
  # Simulate potential final outcome
  meanY1 <- -1+D1+X1+3*X2+3*X3+3*X4
  meanY0 <- 3-D0-1.5*X1+2*X2+2*X3-2*X4
  Y1 <- rnorm(n = n, mean = meanY1, sd = 1)
  Y0 <- rnorm(n = n, mean = meanY0, sd = 1)
  Y <- Y1*Z+(1-Z)*Y0
  
  return(as.data.frame(cbind(X1,X2,X3,X4,Z,D,Y)))
}
```

### 2. Run `PSor.fit`

Now, we simulate a sample data set assuming counterfactual intermediate
independence with $\theta(\mathbf{X})=1$, and then call the main
function to compute principal causal effects under either correctly
specified odds ratio or incorrectly assumed monotonicity.

``` r
library(PSor)
library(SuperLearner)
#> Loading required package: nnls
#> Loading required package: gam
#> Loading required package: splines
#> Loading required package: foreach
#> Loaded gam 1.22-6
#> Super Learner
#> Version: 2.0-29
#> Package created on 2024-02-06
# Generate a data set under independence; or = 1
n = 500
theta = 1
df <-  simu_full_data(n, theta=theta)

# Fit correctly specified odds ratio, or = 1
PSor.fit(
  out.formula = Y~X1+X2+X3+X4,
  ps.formula = D~X1+X2+X3+X4,
  pro.formula = Z~X1+X2+X3+X4,
  df = df,
  out.name = "Y",
  int.name = "D",
  trt.name = "Z",
  cov.names = c("X1","X2","X3","X4"),
  or = 1,
  SLmethods = c("SL.glm", "SL.rpart", "SL.nnet"),
  n.fold = 5,
  scale = "RD",
  alpha = 0.05
)
#>               CDR.Est CDR.SE CDR.ci.low CDR.ci.up DML.Est DML.SE DML.ci.low
#> Always-Takers   2.193  0.330      1.545     2.840   2.148  0.360      1.443
#> Compliers      -0.198  0.444     -1.067     0.672  -0.241  0.460     -1.143
#> Never-Takers   -3.278  0.397     -4.056    -2.500  -3.235  0.405     -4.028
#> Defiers        -0.226  0.399     -1.009     0.556  -0.273  0.402     -1.061
#>               DML.ci.up
#> Always-Takers     2.853
#> Compliers         0.662
#> Never-Takers     -2.441
#> Defiers           0.515

# Fit by incorrectly assuming monotonicity
PSor.fit(
  out.formula = Y~X1+X2+X3+X4,
  ps.formula = D~X1+X2+X3+X4,
  pro.formula = Z~X1+X2+X3+X4,
  df = df,
  out.name = "Y",
  int.name = "D",
  trt.name = "Z",
  cov.names = c("X1","X2","X3","X4"),
  or = Inf,
  SLmethods = c("SL.glm", "SL.rpart", "SL.nnet"),
  n.fold = 5,
  scale = "RD",
  alpha = 0.05
)
#>                    CDR.Est  CDR.SE CDR.ci.lower CDR.ci.upper DML.Est DML.SE
#> Always-Takers (11)   1.438   0.310        0.831        2.045   1.433  0.317
#> Compliers (01)     -35.204 522.863    -1059.997      989.589  48.371 21.161
#> Never-Takers (00)   -2.305   0.298       -2.889       -1.722  -2.277  0.315
#>                    DML.ci.lower DML.ci.upper
#> Always-Takers (11)        0.812        2.055
#> Compliers (01)            6.895       89.846
#> Never-Takers (00)        -2.895       -1.659
```

Here, the function computes estimates under monotonicity by setting
`or = Inf`. The `CDR` estimator uses linear regression for the
continuous outcome and logistic regression for the intermediate outcome
and treatment propensity. For the DML estimator, we will use the
`SuperLearner` package for nuisance function estimation, including the
outcome regression, principal score, and propensity score. The argument
`SLmethods = c("SL.glm", "SL.rpart", "SL.nnet")` specifies the machine
learning algorithms used to estimate the nuisance functions.
