########################################################################
#' Fit joint model
#' @description Bayesian fit of a joint model. The time-to-event process is described using a 
#'       parametric survival model (both AFT and PH are available),
#'       and the longitudinal process is modeled using the linear mixed-effects model. It is assumed that
#'       the association between the two submodels is induced by the longitudinal value.
#'       This version does not allow association induced by the random intercepts and/or slopes.
#'       The Markov Chain Monte Carlo (MCMC) for Bayesian inference is implemented via the 
#'       'rstan' package (Stan Development Team, 2020), which provides the R interface to Stan. 
#' @keywords Joint modeling, MCMC
#' @param surv.fit an object of class \code{\link[survival]{coxph}} from the 'survival' package, representing the Cox PH fit. 
#'        It is required to specify \code{x = TRUE} in \code{coxph}.
#' @param lme.fit an object of class \code{\link[nlme]{lme}} from the 'nlme' package, 
#'        representing the linear mixed-effects model fit.
#' @param surv.model the survival model to be used to describe the event process. Available options are given as follow.
#'        \itemize{ 
#'           \item AFT Models: \code{"weibull"}, \code{"lnormal"}, \code{"llogistic"}, \code{"eweibull"} and \code{"ggamma"}
#'                 for Weibull, log-normal, log-logistic, exponentiated Weibull and generalized gamma distributions, respectively. 
#'           \item PH Models: \code{"weibullph"} and \code{"gllph"} for Weibull and generalized log-logistic distributions, respectively.
#'         }
#'       Default is \code{"weibullph"}. The parameterization of these distributions is described in \code{\link{etime.reg}}.
#' @param fixed.model a character string describing the form of the fixed-effects component in \code{lme}
#'       (this argument is only used when an AFT model is considered to describe the event process). For computational efficiency,
#'       we recommend to use \code{fixed.model = "simple"} if the form is 
#'       \eqn{x'(s)\alpha=\alpha_0+\alpha_1 s} or 
#'       \eqn{\alpha_0+\alpha_1 s + \alpha_2 x_2 + ... + \alpha_q x_q}, where \eqn{s} denotes the time points
#'       at which the longitudinal measurements are recorded, and \eqn{x_2, ... ,x_q} are 
#'       time-independent covariates; otherwise use \code{fixed.model = "nsimple"} (default).
#'      For example, use \code{fixed.model = "nsimple"} if \eqn{x'(s)\alpha} includes an interaction term involving \eqn{s}
#'      or includes \code{ns(s, df)}, where df is the degrees of freedom for a natural cubic spline.
#'      Note that \code{fixed.model = "nsimple"} also works if the form of the fixed-effects component
#'      is "simple" as described above. 
#' @param rand.model a character string describing the form of the random-effects component in \code{lme}
#'       (this argument is only used when an AFT model is considered to describe the event process). For computational efficiency,
#'       we recommend to use \code{rand.model = "simple"} if the form is 
#'       \eqn{w'(s)b=b_0+b_1 s}, where \eqn{b_0} and \eqn{b_1} are 
#'       random intercept and slope, respectively; otherwise use \code{rand.model = "nsimple"} (default).
#'      For example, use \code{rand.model = "nsimple"} if \eqn{w'(s)b} includes \code{ns(s, df)}.
#'      Note that \code{rand.model = "nsimple"} also works if the form of the random-effects component
#'      is "simple" as described above. 
#' @param timevar the name of the time variable in the linear mixed-effects model (a character string).
#' @param form a character string to describe the formulation of the joint model. Available options are
#'      \code{form = "henderson"} (default, can be abbreviated as "hen") for the formulation 
#'       proposed by Henderson et al. (2000), and \code{form = "rizopoulos"} (can be abbreviated as "riz") 
#'       for the formulation described by Rizopoulos (2012). See Details.
#' @param method estimation method to be used (a character string). Available options are 
#'      \code{method = "MCMC"} for Bayesian inference (implemented vis the 'rstan' package), and
#'      \code{method = "2stage"} for the two-stage estimation method. Note that the two-stage approach
#'      may produce biased results as demonstrated by many authors, and hence it is not recommended. 
#'      Since the MCMC approach is computatinally intensive
#'      and time consuming, the two-stage method is included for a quick but crude estimation of the parameters.
#' @param inits initial values of the parameters for MCMC; see Details.
#' @param warmup a positive integer specifying the number of burnin iterations per chain for MCMC (default is 1000).
#' @param iter a positive integer specifying the number of iterations for each chain including warmup (default is 3500). 
#' @param chains a positive integer specifying the number of Markov chains (default is 2). 
#' @param thin thinning interval for monitors (default is 1, the recommended value).
#' @param adapt_delta the target average proposal acceptance probability during Stan's adaptation period (default is 0.8).
#'      In general, you should not need to change \code{adapt_delta} unless you see a warning message about 
#'      divergent transitions. If there is such a warning message,
#'      increase \code{adapt_delta} from the default to a value closer to 1 (e.g., from 0.90 to 0.99).
#'      See the Stan manual for detail.
#' @param max_treedepth a positive integer specifying the maximum treedepth (see the Stan manual).
#'        If there is a warning about transitions exceeding the maximum treedepth,
#'        try increasing the \code{max_treedepth} parameter (e.g., \code{max_treedepth = 15}).
#' @param control a list of control values (see Details).
#' @details The MCMC algorithm for Bayesian inference is implemented in Stan. The \code{coxph} and \code{lme} fits 
#'     (the arguments \code{surv.fit} and \code{lme.fit}) are used to organize the data to be used in Stan. The event 
#'     process can be modeled by one of Weibull AFT, log-normal AFT, log-logistic AFT, exponentiated 
#'     Weibull AFT, generalized gamma AFT, Weibull PH and generalized log-logistic PH models, and the 
#'     longitudinal process is characterized by the linear mixed-effects model. The event time distribution
#'     is characterized by the parameters \eqn{\rho}, \eqn{\kappa} and \eqn{\gamma} as described in \code{\link{etime.reg}}.  
#'
#'     \strong{Longitudinal Process:} We model the longitudinal response \eqn{y_{ij}} at time \eqn{s_{ij}} by the relationship
#'     \eqn{y_{ij}=\mu_i(s_{ij})+U_i(s_{ij})+\epsilon_{ij}}, where \eqn{\mu_i(s_{ij})} is the mean response,
#'     \eqn{U_i(s_{ij})} incorporates subject-specific random effects, and \eqn{\epsilon_{ij}} ~ \eqn{N(0,\sigma^2)} 
#'     is a sequence of mutually independent measurement errors. We assume that the mean response at time \eqn{s}
#'     is characterized by a linear model \eqn{\mu_i(s)=x_i'(s)\alpha}, where \eqn{x_i(s)} is a vector of 
#'     covariates (possibly time-dependent) and \eqn{\alpha} is the corresponding vector of
#'     regression coefficients (fixed effects). For \eqn{U_i(s)}, we assume a linear random effects model
#'     \eqn{U_i(s)=w_i'(s)b_i}, where \eqn{w_i(s)} is a vector of covariates and \eqn{b_i} ~ \eqn{N(0,\Sigma_b)} 
#'     is the corresponding vector of random effects. 
#'
#'    \strong{AFT Event Process:} The event intensity process at time \eqn{t} can be expressed as 
#'     \eqn{\lambda_i(t)=\lambda_0[g_i(t)]} exp\eqn{[-z_i'\beta-V_i(t)]}, where \eqn{\lambda_0(.)} is the 
#'     baseline intensity function of the assumed distribution (see \code{\link{etime.reg}}),
#'     \eqn{g_i(t)=\int_0^t}exp\eqn{[-z_i'\beta-V_i(u)]du}, 
#'     \eqn{z_i} is a vector of baseline covariates, and \eqn{\beta} is the corresponding vector of regression
#'     coefficients (the specificantion of \eqn{V_i(t)} is described below). With this formulation,
#'     the survivor function can be written as \eqn{S_i(t)=S_0[g_i(t)]}, where \eqn{S_0(.)} is the 
#'     baseline survivor function of the assumed distribution (see Cox and Oakes (1984) for detail).
#'
#'    \strong{PH Event Process:} The event intensity process at time \eqn{t} can be expressed as 
#'     \eqn{\lambda_i(t)=\lambda_0(t)} exp\eqn{[z_i'\beta+V_i(t)]} (the specificantion of \eqn{V_i(t)} is described below). 
#'     With this formulation, the survivor function can be written as \eqn{S_i(t)=}exp\eqn{[-\Lambda_i(t)]}, where
#'     \eqn{\Lambda_i(t)=\int_0^t \lambda_0(u)}exp\eqn{[z_i'\beta+V_i(u)]du}.
#'
#'    \strong{Association Structure:} In our implementation, dependence between the longitudinal and the time-to-event 
#'     sub-models is captured through \eqn{V_i(t)}. In Rizopoulos (2012) (see also the 'JM' package), \eqn{V_i(t)} is defined as
#'    \eqn{V_i(t)=\phi[\mu_i(t)+U_i(t)]}, where \eqn{\phi} is the measure of association (regression coefficient for the
#'    time-dependent covariate of the event time process) induced by the fitted longitudinal values. On the other hand,
#'    Henderson et al. (2000) (see also the 'joineR' package) proposed \eqn{V_i(t)=\phi U_i(t)}. 
#'    Both these formulations have been implemented, with \code{form = "henderson"} and \code{form = "rizopoulos"}
#'    for the formulations described by Henderson et al. (2000) and Rizopoulos (2012), respectively.
#'
#'    \strong{Initial Values for MCMC:} Initial values of the parameters can be specified (optional) 
#'    using the \code{inits} argument. Use the "inits via list"
#'    format as described in the 'rstan' package: set inital values by providing a list equal in length to the number
#'    of chains. The elements of this list should themselves be named lists, where
#'    each of these named lists has the name of a parameter and is used to specify
#'    the initial values for that parameter for the corresponding chain. The form of each named list
#'    for a joint model with the survival process described by a two-parameter distribution (\eqn{\rho} and \eqn{\kappa}) is
#' 
#'    \preformatted{list(alpha = a vector, beta = a vector, phi = a scalar, squ_kappa = a scalar, 
#' inv_rho2 = a scalar, b_unscaled = a matrix, rand_cov = a matrix, inv_sigma2 = a scalar),}
#'
#'    where \code{alpha} = \eqn{\alpha}, \code{beta} = \eqn{\beta}, \code{phi} = \eqn{\phi},
#'    \code{squ_kappa} = \eqn{\kappa^2}, \code{inv_rho2} = \eqn{1/\rho^2},
#'    \code{b_unscaled} = unscaled random effects (a matrix with the number of columns equal to the
#'     number of random effects for each subject, and the number of rows equal to the number of subjects), 
#'     \code{rand_cov} = \eqn{\Sigma_b}, and \code{inv_sigma2} = \eqn{1/\sigma^2}.
#'    Note that \code{b_unscaled = solve(t(chol(rand_cov))) \%*\% t(b)}. 
#'    For a joint model with the survival process described by the exponentiated Weibull or the generalized gamma
#'    distribution (three parameters: \eqn{\rho}, \eqn{\kappa} and \eqn{\gamma}), 
#'    add \code{squ_gam}  (a scalar) in the above list, where \code{squ_gam} = \eqn{\gamma^2}.
#'    With the survival process described by the generalized log-logistic PH model, add  
#'    \code{inv_gam2} in the above list, where \code{inv_gam2} = \eqn{1/\gamma^2}. 
#'
#' \code{control}: a list of control values with components:
#'     \itemize{ 
#'     \item \code{alpha_mu}: \eqn{\mu_\alpha} (a vector) for the prior \eqn{\alpha} ~ \eqn{N(\mu_\alpha, \Sigma_\alpha)}, where
#'           \eqn{\alpha} is the vector of fixed-effects coefficients of the longitudinal model.      
#'     \item \code{alpha_sd}: standard deviations for \eqn{\alpha} (a vector), that is, square roots of 
#'           the diagonal elements of \eqn{\Sigma_\alpha}.
#'     \item \code{beta_mu}: \eqn{\mu_\beta} (a vector) for the prior \eqn{\beta} ~ \eqn{N(\mu_\beta, \Sigma_\beta)}, where
#'           \eqn{\beta} is the vector of coefficients for the baseline covariates 
#'           of the survival model (no intercept).      
#'     \item \code{beta_sd}: standard deviations for \eqn{\beta} (a vector), that is, square roots of 
#'          the diagonal elements of \eqn{\Sigma_\beta}.
#'     \item \code{phi_mu}: \eqn{\mu_\phi} (a scalar) for the prior \eqn{\phi} ~ \eqn{N(\mu_\phi, \sigma_\phi)}, where
#'           \eqn{\phi} is the association parameter (regression coefficient for the
#'           time-dependent covariate of the event time process).  
#'     \item \code{phi_sd}: the value of \eqn{\sigma_\phi} (a scalar).
#'     \item \code{A}: a positive definite matrix for the prior \eqn{\Sigma_b} ~ InvWishart(A, df), 
#'            where \eqn{\Sigma_b} is the covariance matrix of the random-effects vector \eqn{b_i}
#'            with dimension equal to the number of random effects for each \eqn{i}.
#'     \item \code{nu}: the value of the df for the inverse Wishart prior InvWishart(A, df).
#'     \item \code{a0}: the shape parameter for the prior \eqn{\sigma^{-2}} ~ gamma(shape, rate).
#'     \item \code{a1}: the rate parameter for the prior \eqn{\sigma^{-2}} ~ gamma(shape, rate).
#'     \item \code{b0}: the shape parameter for the prior \eqn{1/\rho^2} ~ gamma(shape, rate).
#'     \item \code{b1}: the rate parameter for the prior \eqn{1/\rho^2} ~ gamma(shape, rate).
#'     \item \code{c0}: the shape parameter for the prior \eqn{\kappa^2} ~ gamma(shape, rate).
#'     \item \code{c1}: the rate parameter for the prior \eqn{\kappa^2} ~ gamma(shape, rate).
#'     \item \code{d0}: the shape parameter for the prior \eqn{\gamma^2} ~ gamma(shape, rate) \cr
#'            (for the generalized log-logistic PH, \eqn{1/\gamma^2} ~ gamma(shape, rate)).
#'     \item \code{d1}: the rate parameter for the prior \eqn{\gamma^2} ~ gamma(shape, rate) \cr
#'            (for the generalized log-logistic PH, \eqn{1/\gamma^2} ~ gamma(shape, rate)).
#'     \item \code{qpoints}: number of quadrature points to approximate the integral involved
#'           in the survival model; available options are 3, 5, 7, 15, 21, 31 and 41.
#'    \item \code{seed:} the seed for random number generation in Stan.
#'          }
#' @return \code{Survival sub-model:} a description of the survival model used to describe the event process.
#' @return \code{Longitudinal sub-model:} a description the mixed-effects model used to describe the longitudinal process.
#' @return \code{Association structure:} the association structure linking the two processes.
#' @return \code{stan_fit:} posterior summaries of the parameters.
#' @return \code{time:} computational time.
#' @references Cox DR and Oakes D, Analysis of survival data, Chapman and Hall/CRC, 1984.
#' @references Henderson R, Diggle P, and Dobson A, Joint modelling of longitudinal measurements and event time data, 
#'           Biostatistics, 1: 465-480, 2000.
#' @references Rizopoulos D, Joint Models for Longitudinal and Time-to-Event Data: 
#'          With Applications in R, Chapman and Hall/CRC, 2012.
#' @references Stan Development Team, Stan Modeling Language Users Guide and Reference Manual, 2020, 
#'           Version \url{https://mc-stan.org}. 
#' @references Stan Development Team, RStan: the R interface to Stan, 2020, 
#'           Version \url{http://mc-stan.org/}. 
#' @author Shahedul Khan <khan@math.usask.ca>
#' @seealso \code{\link{jm.icriteria}}, \code{\link{jm.lppredict}} , \code{\link{jm.lspredict}}, 
#'      \code{\link{jm.plots}}, \code{\link{jm.resid}},
#'      \code{\link{jm.sim}}, \code{\link{jm.summary}}, \code{\link{jm.surv}}
#' @examples
#'   # Example 1: AIDS data from package 'JM'
#'   library(JM)
#'   data(aids.id)
#'   data(aids)
#'   surv.fit <- coxph(Surv(Time, death) ~ drug + gender + prevOI + AZT,
#'           data = aids.id, x = TRUE)
#'   lme.fit <- lme(CD4 ~ obstime + obstime:drug + gender + prevOI + AZT,
#'       random =  ~ obstime | patient, data = aids)
#'   # Weibull PH with form = "riz"
#'   jmfit.wph0 <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "weibullph", timevar = "obstime", form = "riz")
#'   jmfit.wph0
#'   # Generalized log-logistic PH with form = "riz"
#'   jmfit.gllph0 <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "gllph", timevar = "obstime", form = "riz")
#'   jmfit.gllph0
#'   # Weibull AFT with form = "riz"
#'   jmfit.w0 <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "weibull", timevar = "obstime", form = "riz")
#'   jmfit.w0
#'   # Exponentiated Weibull AFT with form = "hen"
#'   jmfit.ew0 <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "eweibull", timevar = "obstime", form = "hen")
#'   jmfit.ew0
#'   # Summary of jmfit.wph0
#'   jm.summary(jmfit.wph0)
#'   # WAIC and DIC
#'   jm.icriteria(jmfit.wph0)
#'   # MCMC Diagnostic plots
#'   jm.plots(jmfit.wph0, pairs = TRUE, rhat = TRUE, neff = TRUE,
#'     trace = TRUE, density = TRUE, uncertainty.intervals = TRUE)
#'   # Cox-Snell Residuals
#'   jm.resid(jmfit.wph0)
#'
#'   # Example 2: AIDS data from package 'JM'
#'   surv.fit <- coxph(Surv(Time, death) ~ drug + gender + prevOI + AZT,
#'           data = aids.id, x = TRUE)
#'   lme.fit <- lme(CD4 ~ obstime + drug + gender + prevOI + AZT,
#'       random =  ~ obstime | patient, data = aids)
#'   # For AFT models, use fixed.model = "simple" and rand.model = "simple"
#'   # (for computational efficiency)
#'   # Log-normal AFT with form = "hen"
#'   jmfit.ln1 <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "lnormal", fixed.model = "simple",
#'       rand.model = "simple",  timevar = "obstime", form = "hen")
#'   jm.summary(jmfit.ln1)
#'   # Generalized gamma AFT with form = "riz"
#'   jmfit.gg1 <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "ggamma",  fixed.model = "simple",
#'       rand.model = "simple",  timevar = "obstime",  form = "riz")
#'   jm.summary(jmfit.gg1)
#'   # For Weibull PH, use the default fixed.model = "nsimple" and rand.simple = "nsimple".
#'   jmfit.wph1 <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "weibullph", timevar = "obstime", form = "hen")
#'   jm.summary(jmfit.wph1)
#'
#'   # Example 3: pbc data
#'   lme.fit <- lme(log(bilirubin) ~ drug + ns(futime, 2),
#'     data = pbc.long, random = ~ ns(futime, 2) | id)
#'   surv.fit <- coxph(Surv(st, status2) ~ drug * age, data = pbc.surv, x = TRUE)
#'   # Two-stage estimation
#'   jmfit.gg2stage <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "ggamma", method = "2stage", timevar = "futime", form = "riz")
#'   jmfit.gg2stage
#'   # Bayesian estimation
#'   jmfit.ll2 <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "llogistic", timevar = "futime", form = "riz")
#'   jm.summary(jmfit.ll2)
#'   jmfit.wph2 <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "weibullph", timevar = "futime", form = "riz")
#'   jm.summary(jmfit.wph2)
#' @export

jm.reg<-function(surv.fit,lme.fit,surv.model="weibullph",
 fixed.model="nsimple",rand.model="nsimple",
 timevar,form="henderson",method="MCMC",inits=NULL,
 warmup=1000,iter=3500,chains=2,thin=1,
 adapt_delta=0.8,max_treedepth=10,control=list()){
# ---------- Checking -----------------------
if(!inherits(lme.fit, "lme")) stop("\nlme.fit must be from class lme.")
if(!inherits(surv.fit, "coxph")) stop("\nsurv.fit must be from class coxph.")
if(is.null(surv.fit$x)) stop("\nUse x = TRUE in coxph fit.")
if(!(surv.model=="weibull" || surv.model=="llogistic" || surv.model=="lnormal" || 
  surv.model=="eweibull" || surv.model=="ggamma" || surv.model=="weibullph" || surv.model=="gllph"))
stop("\nThe survival model must be one of 'weibull', 'llogistic', 'lnormal', 'eweibull'
  'ggamma', 'weibullph' and 'gllph'")
if(!(form=="riz" || form=="rizopoulos" || form=="hen" || form=="henderson"))
stop("\nThe form must be one of 'rizopoulos' or 'henderson'")
if(!exists("timevar"))
stop("\ntimevar is missing")
# ----------------------------------------
if(surv.model=="weibull" || surv.model=="llogistic" || surv.model=="lnormal" || surv.model=="weibullph"){
return(jmfit.2par(surv.fit=surv.fit,lme.fit=lme.fit,surv.model=surv.model,
 fixed.model=fixed.model,rand.model=rand.model,
 timevar=timevar,inits=inits,form=form,method=method,
 warmup=warmup,iter=iter,chains=chains,thin=thin,
 adapt_delta=adapt_delta,max_treedepth=max_treedepth,control=control))
} else{
return(jmfit.3par(surv.fit=surv.fit,lme.fit=lme.fit,surv.model=surv.model,
 fixed.model=fixed.model,rand.model=rand.model,
 timevar=timevar,inits=inits,form=form,method=method,
 warmup=warmup,iter=iter,chains=chains,thin=thin,
 adapt_delta=adapt_delta,max_treedepth=max_treedepth,control=control))
}
}

