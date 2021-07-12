#' The 'JMR' package.
#'
#' @description While the main focus of this package involves recurrent event data analysis and joint modeling, 
#'    it also includes fitting regression models for right censored (non-recurrent) time-to-event data. 
#'    For the time-to-event process, available options are Weibull, log-logistic, log-normal, exponentiated Weibull and generalized gamma 
#'    accelerated failure time models, and Weibull and generalized log-logistic proportional hazards models. 
#'    This package also includes several utility functions, including functions for residual analysis, data simulation and dynamic predictions. 
#'    Maximum likelihood method is used to fit recurrent event models, whereas Bayesian approach (implemented in Stan) 
#'    is used for joint modeling. The 'rstan' package is used to parse, compile, test, estimate and analyze Stan models. 
#'
#' @docType package
#' @name JMR-package
#' @aliases JMR
#' @useDynLib JMR, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @author Shahedul Khan <khan@math.usask.ca>
#'
#' @references
#' Cook RJ and Lawless J, The statistical analysis of recurrent events, Springer, 2007.
#'
#' Cox DR and Oakes D, Analysis of survival data, Chapman and Hall/CRC, 1984.
#'
#'  Gabry J and Mahr T, bayesplot: Plotting for Bayesian Models, R package version 1.7.2, 
#'    \url{https://mc-stan.org/bayesplot}.
#'
#' Gelman A, Hwang J, and Vehtari A, Understanding predictive information criteria for Bayesian
#'    models, Statistics and Computing, 24: 997-1016, 2014.
#'
#' Henderson R, Diggle P, and Dobson A, Joint modelling of longitudinal measurements and event
#'   time data, Biostatistics, 1: 465-480, 2000.
#'
#' Kalbfleisch JD and Prentice RL, The statistical analysis of failure time data, Wiley, 2002.
#'
#' Rizopoulos D, JMbayes: Joint Modeling of Longitudinal and Time-to-Event Data under a Bayesian
#'    Approach, R package version 0.8-85, 2020, \url{https://cran.r-project.org/web/packages/JMbayes}.
#'
#' Rizopoulos D, The R package JMbayes for fitting joint models for longitudinal and time-to-event
#'    data Using MCMC, Journal of Statistical Software, 72(7): 1-45, 2016.
#'
#' Rizopoulos D, Joint Models for Longitudinal and Time-to-Event Data: With Applications in R,
#'    Chapman and Hall/CRC, 2012.
#'
#' Rizopoulos D and Ghosh P, A Bayesian semiparametric multivariate joint model for multiple longitudinal
#'    outcomes and a time-to-event, Statistics in Medicine, 30: 1366-1380, 2011.
#'
#' Stan Development Team, Stan Modeling Language Users Guide and Reference Manual, 2020, Version
#'  \url{https://mc-stan.org}.
#'
#' Stan Development Team, RStan: the R interface to Stan, 2020, Version \url{http://mc-stan.org/}.
#'
#' Therneau T and Grambsch P, Modeling survival data: extending the Cox Model, Springer-Verlag, 2000.
#'
#' Zhou H and Hanson T, A unified framework for fitting Bayesian semiparametric models to arbitrarily
#'    censored survival data, including spatially referenced data, Journal of the American Statistical
#'    Association, 113(522), 571-581, 2018.
NULL
