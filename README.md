
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesWinRatio

<!-- badges: start -->
<!-- badges: end -->

This package performs Bayesian monitoring using the win ratio approach
based on the manuscript “The Win Ratio Approach in Bayesian Monitoring
for Two-Arm Phase II Clinical Trial Designs with Multiple Time-to-Event
Endpoints”.

## Installation

You can install the development version of BayesWinRatio from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("xinran-h/BayesWinRatio")
```

## BayesWinRatio in a nutshell

The function OCC.Table is the core function of BayesWinRatio. This
function first generates simulation data, perform Bayesian monitoring
using this simulation data, and returns the operating characteristics of
the Bayesian monitoring. The user inputs are:

- N.sim: Number of simulations.
- N.max: Maximum number of patients to enroll.
- mu.trt: True mean for treatment arm (log).
- Sigma.trt: True variance for treatment arm (log).
- mu.ctrl: True mean for control arm (log).
- Sigma.ctrl: True variance for control arm (log).
- cens_upper: Upper limit for the censoring time.The censoring time is
  generated from Uniform(0, cens_upper).
- design: A numeric value indicating the type of design. 1 = proposed
  design, 2 = time to recurrence design, 3= time to death design, 4 =
  time to first event design.
- cohort: Interim cohort.
- recruit.int: Recruitment interval.
- m0: Prior mean for mu.trt/mu.ctrl.
- L0: Prior variance for mu.trt/mu.ctrl.
- v0: For the proposed design, this is the prior degrees of freedom for
  Sigma.trt/Sigma.ctrl. For the traditional designs, v0/2 is the prior
  shape for Sigma.trt/Sigma.ctrl.
- S0:For the proposed design, this is the prior scale matrix for
  Sigma.trt/Sigma.ctrl. For the traditional designs, v0\*S0/2 is the
  prior scale for Sigma.trt/Sigma.ctrl.
- time_max: The upper limit for the recurrence and death time sampled
  from truncated normal. This will set the upper limit to to time_max
  rather than Inf.
- eta: A pre-specified lower bound of acceptable performance based on
  historical information.
- lambda: Cutoff parameter.
- thin_MCMC: Thinning degree.
- Niter: Number of iterations for gibbs sampler.

The output is a list with the following components:

- mysim1: A numeric vector. mysim1 = c(PRN, PEN, EN), where PRN is the
  average percentage of trials that are not stopped, PEN is the average
  percentage of early termination, and EN is the average number of
  patients.  
- mu: A numeric vector. mu = c(mu.trt,mu.ctrl).
- Sigma.trt: True covariance matrix for treatment arm (log), same as the
  input.
- Sigma.ctrl: True covariance matrix for control arm (log), same as the
  input.
- lambda: Cutoff parameter.

### Tuning lambda

To access operating characteristics, we first tune lambda, the cutoff
parameter. This parameter is calibrated by controlling the percentage of
early termination (PET) at a desirable level (say, 0.1) under the null
scenario. The following code uses a bisection method to find the lambda.

``` r
library(BayesWinRatio)
RNGkind("L'Ecuyer-CMRG")
set.seed(126)

# PET = 0.1 under the null scenario
Cutoff2Prob <- 0.1 


# Define the bisection parameters
a <- 3  # The lowest value of mylambda
b <- 100  # The highest value of mylambda
tolerance <- 0.1   # Tolerance for stopping the bisection

# Initialize the bisection loop
while (b - a > tolerance) {
  c <- (a + b) / 2
  result <- OCC.Table(N.sim = 10000,
                      N.max = 100,
                      mu.trt = c(log(2), log(5)),
                      Sigma.trt = matrix(c(1, 0.5, 0.5, 1),ncol=2, byrow = T),
                      mu.ctrl= c(log(2), log(5)),
                      Sigma.ctrl = matrix(c(1, 0.5, 0.5, 1), ncol=2, byrow = T), 
                      cens_upper = 25,
                      design = 1, 
                      cohort = c(40,60,80,100),
                      recruit.int  = 0.25,
                      m0 = c(0,0),
                      L0 = diag(10^6, 2),
                      v0 = 4,
                      S0 = diag(10^(-6), 2),
                      time_max = 40,
                      eta = 1.5,
                      lambda = c,
                      thin_MCMC = 5,
                      Niter = 100000);    
  
  if (result$mysim1[2] == Cutoff2Prob) {
    break
  } 
  else if (result$mysim1[2] < Cutoff2Prob) {
    b <- c
  } else {
    a <- c
  }
}

# The maximum lambda value that satisfies the condition is approximately 'b'
max_lambda = b
```

### Accessing operating characterstics

We run the OCC_Table function to obtain the operating characteristics
using the lambda just calibrated.

``` r
RNGkind("L'Ecuyer-CMRG")
set.seed(126)

mu.death = c(log(2.5), log(3), log(3.5), log(4), log(4.5))
#varying mu_death for treatment arm
for (i in 1:length(mu.death)){
  
  result <- OCC.Table(N.sim = 10000,
                      N.max = 100,
                      mu.trt = c(log(2), mu.death[i]),
                      Sigma.trt = matrix(c(1, 0.5, 0.5, 1), ncol=2, byrow = T),
                      mu.ctrl= c(log(2), log(5)),
                      Sigma.ctrl = matrix(c(1, 0.5, 0.5, 1), ncol=2, byrow = T), 
                      cens_upper = 25,
                      design = 1,
                      cohort = c(40,60,80,100),
                      recruit.int  = 0.25,
                      m0 = c(0,0),
                      L0 = diag(10^6, 2),
                      v0 = 4,
                      S0 = diag(10^(-6), 2),
                      time_max = 40,
                      eta = 1.5,
                      lambda =  max_lambda,
                      thin_MCMC = 5,
                      Niter = 100000);       
  
  saveRDS(result, paste0("~/simulation results/varying_death", i,".rds"))
  
}
```
