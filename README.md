# JMR: An R Package for Recurrent Event Survival Analysis and Joint Modeling

While the main focus of this package is recurrent event data analysis and joint modeling, it also includes fitting regression models for right censored (non-recurrent) time-to-event data. Available options are Weibull, log-logistic, log-normal, exponentiated Weibull and generalized gamma accelerated failure time models, and Weibull and generalized
log-logistic proportional hazards models. This package also includes residual analysis, data simulation and dynamic predictions. Maximum likelihood method is used to fit recurrent event models, whereas Bayesian approach (**implemented in STAN**) is used for joint modeling. The rstan package is used to parse, compile, test, estimate and analyze Stan models. **See the pdf manual, JMR_0.0.0.9000.pdf, for details**.

# Inatallation: 

Step 1: Install the following R packages: rstan (>= 2.18.1), bayesplot (>= 1.7.1), coda (>= 0.19-3), numDeriv (>= 2016.8-1.1), Rmpfr (>= 0.8-1)

Step 2: From R,

library(devtools)

install_github("sa4khan/JMR")

# Functions for Recurrent Event Data Analysis

**etime.reg**: Fit a parametric time-to-event regression model. It includes both standard (non-recurrent) survival data analysis and recurrent event
data analysis.

**etime.resid**: Cox-Snell residual analysis of the etime.reg fit.

**LR.test**: Likelihood ratio test for Weibull AFT as a submodel of the exponentiated Weibull AFT and generalized gamma AFT models.

**sim.rec**: Simulate recurrent event data for fixed covariates.

# Main functions for Joint Modeling

**jmreg.aft**: Fit a joint model

**jm.summary**: Summary of a joint model fit

**jm.resid.plot**: Cox-Snell residual plot for the event process of the joint model

**jm.reffects**: Posterior means/medians of the random effects from a joint model fit

**jm.DIC**: Computes DIC for Bayesian fit of the joint model

**jm.WAIC**: Computes WAIC for Bayesian fit of the joint model

**jm.surv**: Dynamic predictions of survival probabilites

**jm.sim**: Simulate from joint models
