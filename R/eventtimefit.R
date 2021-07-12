##############################################
#' Fit a parametric time-to-event regression model
#' @description Fit a parametric time-to-event regression model for right censored data using the maximum likelihood method. 
#'     It includes both standard (non-recurrent) survival data analysis and recurrent event data analysis. Available options are
#'     Weibull, log-logistic, log-normal, exponentiated Weibull and generalized gamma accelerated failure time (AFT) models, 
#'     and Weibull and generalized log-logistic proportional hazards (PH) models.
#' @keywords Recurrent event data analysis, Time-to-event regression
#' @importFrom Rmpfr mpfr
#' @importFrom numDeriv hessian
#' @param Formula a formula of the form \code{response ~ model}. The \code{response} is 
#'        \code{c(st, status)} for standard survival analysis with right censored data, and \code{c(start, stop, status)} for recurrent event
#'        data analysis, where
#'        \itemize{ 
#'            \item \code{st}: survival time for right censored data,
#'            \item \code{status}: censoring indicator,
#'            \item \code{start} and \code{stop}: the start and end time of each time interval for recurrent event data. 
#'          }
#'        The linear predictor is speficified by the \code{model} argument of the Formula expression. 
#'        It consists of a series of terms (covariates) 
#'        separated by + operators. 
#'        It also allows \code{:} and \code{*} operators to define interaction terms. 
#'        \itemize{
#'           \item Example 1: \code{c(st, status) ~ z1 + z2 + z3}
#'           \item Example 2: \code{c(start, stop, status) ~ z1 + z2 * z3}
#'        }
#' @param init initial values for the parameters (optional). It is a matrix with each row has one set
#'       of initial values. If there are \eqn{p} 
#'       regression coefficients \eqn{\beta_1}, \eqn{\beta_2}, ...,  
#'       \eqn{\beta_p}, then each row has initial values with the following sequence:
#'       \eqn{\beta_1}, \eqn{\beta_2}, ..., 
#'       \eqn{\beta_p}, log(\eqn{\kappa}), log(\eqn{\gamma}), log(\eqn{\rho})
#'       for the exponentiated Weibull, generalized 
#'       gamma and generalized log-logistic models, and \eqn{\beta_1}, \eqn{\beta_2}, ..., 
#'       \eqn{\beta_p}, log(\eqn{\kappa}), log(\eqn{\rho})
#'       for the Weibull, log-logistic and log-normal models, where \eqn{\kappa} and 
#'       \eqn{\gamma} are the
#'       shape parameters and \eqn{\rho} is the rate parameter (see below). That is, multiple sets of
#'       initial values can be given. For example, for the exponentiated Weibull model, 
#'       \code{init = rbind(c(rep(0, p), rep(0, 3)), c(rep(0, p), rep(0.5, 3)))} gives two sets of 
#'       initial values.
#' @param Data a data frame in which to interpret the variables named in the Formula.
#' @param surv.model assumed distribution for survival times. Available options are given as follows.
#'        \itemize{ 
#'           \item AFT Models: \code{"weibull"}, \code{"lnormal"}, \code{"llogistic"}, \code{"eweibull"} and \code{"ggamma"}
#'                 for Weibull, log-normal, log-logistic, exponentiated Weibull and generalized gamma distributions, respectively. 
#'           \item PH Models: \code{"weibullph"} and \code{"gllph"} for Weibull and generalized log-logistic distributions, respectively.
#'         }
#'       Default is \code{"weibull"}. The probability density functions of these distributions are:  
#'        \itemize{ 
#'           \item Weibull: \eqn{f(t)=\kappa\rho(\rho}\eqn{t)^{\kappa-1}} exp\{\eqn{-(\rho}\eqn{t)^\kappa}\} 
#'           \item Log-logistic: \eqn{f(t) = \kappa\rho(\rho}\eqn{t)^{\kappa-1}/[1+(\rho}\eqn{t)^\kappa]^2}
#'           \item Log-normal: \eqn{f(t) = \kappa}
#'                              exp\eqn{\{-(\kappa}log\eqn{(\rho}\eqn{t))^2/2\}/[t(2\pi)^{1/2}]}
#'            \item Generalized gamma: \eqn{f(t) = \kappa\rho(\rho}\eqn{t)^{\kappa\gamma-1}} 
#'                              exp\eqn{\{-(\rho}\eqn{t)^\kappa\}/\Gamma(\gamma)}
#'            \item Exponentiated Weibull: \eqn{f(t) =\kappa\gamma\rho(\rho}\eqn{t)^{\kappa-1}
#'                       (1-} exp\eqn{\{-(\rho}\eqn{t)^\kappa\})^{\gamma-1}} exp\eqn{\{-(\rho}\eqn{t)^\kappa\}} 
#'            \item Generalized log-logistic: \eqn{f(t) = \kappa\rho(\rho}\eqn{t)^{\kappa-1}/[1+(\gamma}\eqn{t)^\kappa]^{(\rho/\gamma)^\kappa+1}}
#'         }
#' @param conf.int confidence level (default is 0.95).
#' @param iter.max maximum number of iterations for optimization (default is 250).
#' @details The AFT and the PH models are specified using the hazard functions 
#'       \eqn{h(t; z) = h_0[t} exp\eqn{(-z'\beta)]} exp\eqn{(-z'\beta)} and
#'       \eqn{h(t; z) = h_0(t)} exp\eqn{(z'\beta)}, respectively,
#'       where \eqn{h_0(.)} is the baseline hazard function of the assumed distribution, 
#'      \eqn{z} is the \eqn{p} x 1 vector of covariates, and \eqn{\beta} is the corresponding vector
#'       of regression coefficients (see Section 2.3 of Kalbfleisch and 
#'       Prentice, and Section 3.2 of Cook and Lawless). The baseline hazard functions of the 
#'        distributions are:
#'        \itemize{ 
#'           \item Weibull: \eqn{h(t) = \kappa\rho(\rho}\eqn{t)^{\kappa-1}}
#'           \item Log-logistic: \eqn{h(t) = \kappa\rho(\rho}\eqn{t)^{\kappa-1}/[1+(\rho}\eqn{t)^\kappa]}
#'           \item Log-normal: \eqn{h(t) = f(t)/S(t)}, whete \eqn{S(t)} is the survivor function of the log-normal distribution
#'           \item Generalized gamma: \eqn{h(t) = f(t)/S(t)}, whete \eqn{S(t)} is the survivor function of the generalized gamma distribution
#'           \item Exponentiated Weibull: \eqn{h(t) =\kappa\gamma\rho(\rho}\eqn{t)^{\kappa-1}
#'                       (1-} exp\eqn{\{-(\rho}\eqn{t)^\kappa\})^{\gamma-1}} exp\eqn{\{-(\rho}\eqn{t)^\kappa\}/
#'                       [1-(1-} exp\eqn{\{-(\rho}\eqn{t)^\kappa\})^\gamma]}
#'          \item Generalized log-logistic: \eqn{h(t) = \kappa\rho(\rho}\eqn{t)^{\kappa-1}/[1+(\gamma}\eqn{t)^\kappa]} 
#'          }
#'         For recurrent event analysis, each line of data for a given subject must include the start time and 
#'         stop time for each interval of follow-up.
#' @return \code{model:} the survival model.
#' @return \code{data summary:} number of observations, number of events, number of censored observations, and number of predictors.
#' @return \code{fit:} estimate, standard error, z, p value, and confidence interval.
#' @return \code{exp(coef)} and \code{exp(-est):} estimate and confidence interval.
#' @return \code{fit criteria:} log-likelihood, deviance, and AIC.
#' @return \code{optimizer:} nlminb or optim.
#' @return \code{cov:} covariance matrix.
#' @return \code{st:} survival times (for recurrent event, an n x 2 matrix for start and stop times).
#' @return \code{status:} censoring indicator.
#' @return \code{design.mat:} design matrix.
#' @references Cook RJ and Lawless J, The statistical analysis of recurrent events, Springer, 2007.
#' @references Kalbfleisch JD and Prentice RL, The statistical analysis of failure time data, Wiley, 2002.
#' @author Shahedul Khan <khan@math.usask.ca>
#' @seealso \code{\link{LR.test}}, \code{\link{etime.resid}}, 
#'          \code{\link[survival]{survreg}}, \code{\link[survival]{coxph}}.
#' @examples 
#' # Example 1: Recurrent Event Analysis
#'   library(frailtypack)
#'   data(readmission)
#'   fit.gg <- etime.reg(c(t.start, t.stop, event) ~ sex + chemo + charlson,
#'       surv.model = "ggamma",  Data = readmission)
#'   fit.gg
#'   fit.gg$cov
#'
#' # Example 2: Non-recurrent Data (Right Censored)
#'   library(survival)
#'   fit.ew <- etime.reg(c(time, status) ~ karno + diagtime + age + prior + celltype + trt,
#'       Data = veteran, surv.model = "eweibull")
#'   fit.ew
#'   fit.ew$design.mat
#' @export

etime.reg<-function(Formula,init=NULL,Data,surv.model="weibull",conf.int=0.95,iter.max=250){
if(!inherits(Formula,"formula")){
stop("Formula must be a formula of the form
'c(time, status) ~ covariates separated by +' for survival analysis, and 
'c(start, stop, status) ~ covariates separated by +' for recurrent event analysis.")
}
if(!(length(all.vars(Formula[1:2]))==2 || length(all.vars(Formula[1:2]))==3)){
stop("Formula must be a formula of the form
'c(time, status) ~ covariates separated by +' for survival analysis, and 
'c(start, stop, status) ~ covariates separated by +' for recurrent event analysis.")
}
if(!is.null(init)){
if(!is.matrix(init)){
stop("init must be a matrix with one set of initial values in each row")
}
}
allvar<-all.vars(Formula)
if(isFALSE(all(allvar %in% colnames(Data)))){
stop("variable names in 'Formula' did not match
     with the variable names in 'Data'")
}
if(isFALSE(surv.model=="weibull" || surv.model=="llogistic" || 
   surv.model=="eweibull" || surv.model=="ggamma" || 
   surv.model=="lnormal" || surv.model=="weibullph" || 
   surv.model=="gllph")){
stop("surv.model must be 'weibull' or 'llogistic' or 'eweibull' or
     'ggamma' or 'lnormal' or 'weibullph' or 'gllph'")
}
ind<-length(all.vars(Formula[1:2]))
if(surv.model=="weibull"){
if(ind>2){
fit<-Rwaft.fit(Formula,init=init,data=Data,conf.int=conf.int,iter.max=iter.max)
} else{
fit<-waft.fit(Formula,init=init,data=Data,conf.int=conf.int,iter.max=iter.max)
}
}
if(surv.model=="llogistic"){
if(ind>2){
fit<-Rllaft.fit(Formula,init=init,data=Data,conf.int=conf.int,iter.max=iter.max)
} else{
fit<-llaft.fit(Formula,init=init,data=Data,conf.int=conf.int,iter.max=iter.max)
}
}
if(surv.model=="eweibull"){
if(ind>2){
fit<-Rewaft.fit(Formula,init=init,data=Data,conf.int=conf.int,iter.max=iter.max)
} else{
fit<-ewaft.fit(Formula,init=init,data=Data,conf.int=conf.int,iter.max=iter.max)
}
}
if(surv.model=="ggamma"){
if(ind>2){
fit<-Rggamma.fit(Formula,init=init,data=Data,conf.int=conf.int,iter.max=iter.max)
} else{
fit<-ggamma.fit(Formula,init=init,data=Data,conf.int=conf.int,iter.max=iter.max)
}
}
if(surv.model=="lnormal"){
if(ind>2){
fit<-Rlnaft.fit(Formula,init=init,data=Data,conf.int=conf.int,iter.max=iter.max)
} else{
fit<-lnaft.fit(Formula,init=init,data=Data,conf.int=conf.int,iter.max=iter.max)
}
}
if(surv.model=="weibullph"){
if(ind>2){
fit<-Rwph.fit(Formula,init=init,data=Data,conf.int=conf.int,iter.max=iter.max)
} else{
fit<-wph.fit(Formula,init=init,data=Data,conf.int=conf.int,iter.max=iter.max)
}
}
if(surv.model=="gllph"){
if(ind>2){
fit<-Rgll.fit(Formula,init=init,data=Data,conf.int=conf.int,iter.max=iter.max)
} else{
fit<-gll.fit(Formula,init=init,data=Data,conf.int=conf.int,iter.max=iter.max)
}
}
return(fit)
}
