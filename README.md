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
data analysis. Available options are Weibull, log-normal, log-logistic, exponentiated Weibull and generalized gamma AFT models, and Weibull and generalized log-logistic PH models.

**etime.resid**: Cox-Snell residual analysis of the etime.reg fit.

**LR.test**: Likelihood ratio test for Weibull AFT as a submodel of the exponentiated Weibull AFT and generalized gamma AFT models.

**sim.rec**: Simulate recurrent event data for fixed covariates.

# Main functions for Joint Modeling

**jm.reg**: Bayesian fit of a joint model. The time-to-event process is described using a parametric survival model (available options are Weibull, log-normal, log-logistic, exponentiated Weibull and generalized gamma AFT models, and Weibull and generalized log-logistic PH models), and the longitudinal process is modeled using the linear mixed-effects model. It is assumed that the association between the two submodels is induced by the longitudinal value. The Markov Chain Monte Carlo (MCMC) for Bayesian inference is implemented via the "rstan" package, which provides the R interface to STAN.

**jm.summary**: Produces a summary of the Bayesian joint model fit. It uses a jm.reg fit as its argument.

**jm.icriteria**: Computes DIC and WAIC for a Bayesian joint model fit. It uses a jm.reg fit as its argument.

**jm.resid**: Returns posterior summaries of the Cox-Snell residuals for a jm.reg fit. It also produces a residual plot.

**jm.plots**: This function produces plots for posterior analysis. It includes pair plots, rhat values as either points or a histogram, ratios of effective sample
size to total sample size as either points or a histogram, trace plots, density plots, and plots of uncertainty intervals.

**jm.reffects**: Returns posterior summary of the random effects.

**jm.surv**: Dynamic predictions of survival probabilites.

**jm.lppredict**: Given a jm.reg fit and a data frame, this function produces marginal predictions (population-level) for the longitudinal outcome.

**jm.lspredict**: Given a jm.reg fit and a data frame, this function produces subject-specific predictions for the longitudinal outcome.

**jm.sim**: This function simulates longitudinal responses and event times from joint models. Available options for the time-to-event submodel are Weibull AFT, log-logistic AFT, log-normal AFT, exponentiated Weibull AFT, generalized gamma AFT, Weibull PH and generalized log-logistic PH.
