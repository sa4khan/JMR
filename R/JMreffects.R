####################################################################
# Posterior summaries of the random effects
###################################################################
#' Posterior summaries of the random effects from a joint model fit
#' @description Returns posterior summary of the random effects. Note that the random effects are shared between the 
#'     longitudinal and the survival components, and the link between these two processes via the random effects 
#'     is commonly known as \emph{latent association}. The formulation of the joint model
#'     is described in \code{\link{jm.reg}}.
#' @keywords Joint modeling, MCMC, Posterior summary, Random effects
#' @param jmfit a \code{\link{jm.reg}} fit with \code{method = "MCMC"}.
#' @param rhat.plot if \code{TRUE}, produce a plot for rhats. 
#' @param neff.plot if \code{TRUE}, produce a plot for ratios of effective sample size to total sample size. 
#' @details The random effects are monitored in MCMC simulations in Stan. The summary method in rstan is used
#'        to produce the posterior summaries of the random effects.
#' @return \code{b:} posterior summaries of the random effects.
#' @return \code{plots:} plots if \code{rhat.plot = TRUE} and/or \code{neff.plot = TRUE}.
#' @author Shahedul Khan <khan@math.usask.ca>
#' @seealso \code{\link{jm.reg}}, \code{\link{jm.summary}}
#' @examples
#'   # Example: pbc data
#'   lme.fit <- lme(log(bilirubin) ~ drug + ns(futime, 2),
#'     data = pbc.long, random = ~ ns(futime, 2) | id)
#'   surv.fit <- coxph(Surv(st, status2) ~ drug * age, data = pbc.surv, x = TRUE)
#'   # Bayesian estimation
#'   jmfit.wph2 <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "weibullph", timevar = "futime", form = "riz")
#'   jm.reffects(jmfit.wph2, rhat.plot = TRUE, neff.plot = TRUE)
#' @export

jm.reffects<-function(jmfit,rhat.plot=FALSE,neff.plot=FALSE){
results<-jmfit$all.results
fit<-results$stan.fit
b<-rstan::summary(fit)[[1]][results$b.loc,]
if(rhat.plot){
r<-rhat(results$stan.fit,
  pars=names(results$stan.fit)[c(results$b.loc)])
windows()
print(mcmc_rhat(r))
}
if(neff.plot){
neff.ratio<-neff_ratio(fit, pars = names(fit)[results$b.loc])
windows()
print(mcmc_neff(neff.ratio) +
  ggtitle("Ratios of Effective Sample Size to Total Sample Size\n(lighter is better)\n")+ 
  theme(plot.title = element_text(hjust = 0.5)))
}
output<-list(b=b)
class(output)<-"JMR"
attr(output, "hidden")<-NULL
return(output)
}
###################

