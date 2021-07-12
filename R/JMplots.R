########################################################################
#' MCMC plots for posterior analysis
#' @description This function produces plots for posterior analysis. 
#'    The 'bayesplot' package is used for plotting. Available options are
#'    pair plots, rhat values as either points or a histogram, ratios of effective sample size to total sample size
#'     as either points or a histogram, trace plots, density plots, and plots of uncertainty intervals. See the
#'    'bayesplot' package for detail.
#' @keywords Joint modeling, MCMC, Pair plots, Rhat values, Effective sample size, Trace plot, Density plot, Uncertainty interval.
#' @param jmfit a \code{\link{jm.reg}} fit with \code{method = "MCMC"}. 
#' @param pairs if \code{TRUE}, produce pair plots. 
#' @param rhat if \code{TRUE}, produce a plot for rhats. 
#' @param neff if \code{TRUE}, produce a plot for ratios of effective sample size to total sample size. 
#' @param trace if \code{TRUE}, produce trace plots. 
#' @param density if \code{TRUE}, produce density plots. 
#' @param uncertainty.intervals if \code{TRUE}, produce uncertainty intervals computed from parameter draws. 
#' @details The 'bayesplot' package is used for plotting. See the
#'    'bayesplot' package for detail.
#' @return Plots for posterior analysis. 
#' @references Gabry J and Mahr T, bayesplot: Plotting for Bayesian Models, R
#'    package version 1.7.2, \url{https://mc-stan.org/bayesplot}
#' @author Shahedul Khan <khan@math.usask.ca>
#' @seealso \code{\link{jm.reg}}, \code{\link{jm.resid}}
#' @examples 
#' # Example: AIDS data from package JM
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
#'   jm.plots(jmfit.wph0, pairs = TRUE, rhat = TRUE, neff = TRUE,
#'     trace = TRUE, density = TRUE, uncertainty.intervals = TRUE)
#' @export

jm.plots<-function(jmfit,pairs=TRUE,rhat=FALSE,neff=FALSE,
  trace=FALSE,density=FALSE,uncertainty.intervals=FALSE){
if(class(jmfit)!="JMR"){
stop("jmfit must inherit from class JMR")
}
results<-jmfit$all.results
surv.data<-results$surv.data
surv.covariates<-colnames(surv.data)[-(1:3)]
long.data<-results$long.data
long.covariates<-colnames(long.data)[-(1:2)]
#-----
beta.loc<-results$beta.loc
phi.loc<-results$phi.loc
kappa.loc<-results$kappa.loc
if(results$surv.model=="eweibull" || results$surv.model=="ggamma" || results$surv.model=="gllph"){
gam.loc<-results$gam.loc
}
rho.loc<-results$rho.loc
alpha.loc<-results$alpha.loc
sigma.loc<-results$sigma.loc
cov.loc<-results$cov.loc
b.loc<-results$b.loc
#------
mat0<-matrix(cov.loc,ncol=sqrt(length(cov.loc)))
rsd.loc<-diag(mat0)
rcov.loc0<-upper.tri(mat0,diag=FALSE)
rcov.loc<-c(t(mat0))[c(t(rcov.loc0))]
#posterior[1:2,1,c(rsd.loc,rcov.loc)]
#------
vm<-sqrt(length(cov.loc))
mat <- matrix(0, vm, vm)
mat1<-which(mat ==0, arr.ind = T)
mat1[,1:2]<-mat1[,2:1]
mat2<-matrix(mat1[mat1[,"row"]<mat1[,"col"],],ncol=2)-1
gg<-NULL
for(k in 1:(vm*(vm-1)/2)){
gg<-c(gg,paste("cov(b",paste(mat2[k,1]),", b",paste(mat2[k,2]),")",sep = ''))
}
gg0<-(1:vm)-1
gg1<-NULL
for(k in 1:vm){
gg1<-c(gg1,paste("sd(b",paste(gg0[k]),")",sep = ''))
}
#-----
posterior <- as.array(results$stan.fit)
if(results$surv.model=="eweibull" || results$surv.model=="ggamma" || results$surv.model=="gllph"){
surv.par<-posterior[,,c(beta.loc,phi.loc,rho.loc,kappa.loc,gam.loc)]
surv.regpar<-posterior[,,c(beta.loc,phi.loc)]
surv.distpar<-posterior[,,c(rho.loc,kappa.loc,gam.loc)]
dimnames(surv.par)[[3]]<-c(surv.covariates,"association","rho","kappa","gamma")
dimnames(surv.regpar)[[3]]<-c(surv.covariates,"association")
dimnames(surv.distpar)[[3]]<-c("rho","kappa","gamma")
} else{
surv.par<-posterior[,,c(beta.loc,phi.loc,rho.loc,kappa.loc)]
surv.regpar<-posterior[,,c(beta.loc,phi.loc)]
surv.distpar<-posterior[,,c(rho.loc,kappa.loc)]
dimnames(surv.par)[[3]]<-c(surv.covariates,"association","rho","kappa")
dimnames(surv.regpar)[[3]]<-c(surv.covariates,"association")
dimnames(surv.distpar)[[3]]<-c("rho","kappa")
}
long.par<-posterior[,,c(alpha.loc,sigma.loc,rsd.loc,rcov.loc)]
long.regpar<-posterior[,,c(alpha.loc)]
long.covpar<-posterior[,,c(sigma.loc,rsd.loc,rcov.loc)]
dimnames(long.par)[[3]]<-c(long.covariates,"sd(resid)",gg1,gg)
dimnames(long.regpar)[[3]]<-c(long.covariates)
dimnames(long.covpar)[[3]]<-c("sd(resid)",gg1,gg)
#------
if(pairs){
np<-nuts_params(results$stan.fit)
windows()
print(mcmc_pairs(surv.regpar,np=np,
  grid_args = list(top="Regression Coefficients for the Survival Sub-Model\n")))
windows()
print(mcmc_pairs(surv.distpar,np=np,
  grid_args = list(top="Distributional Parameters for the Survival Sub-Model\n")))
windows()
print(mcmc_pairs(long.regpar,np=np,
  grid_args = list(top="Fixed-Effects Coefficients for the Longitudinal Sub-Model\n")))
windows()
print(mcmc_pairs(long.covpar,np=np,
  grid_args = list(top="Variance-Covariance Parameters for the Longitudinal Sub-Model\n")))
}
#-------------------
if(rhat){
if(results$surv.model=="eweibull" || results$surv.model=="ggamma" || results$surv.model=="gllph"){
r<-rhat(results$stan.fit,
  pars=names(results$stan.fit)[c(beta.loc,phi.loc,rho.loc,kappa.loc,gam.loc,
  alpha.loc,sigma.loc,rsd.loc,rcov.loc)])
} else{
r<-rhat(results$stan.fit,
  pars=names(results$stan.fit)[c(beta.loc,phi.loc,rho.loc,kappa.loc,
  alpha.loc,sigma.loc,rsd.loc,rcov.loc)])
}
windows()
print(mcmc_rhat(r))
}
#--------------
if(neff){
neff.ratio<-neff_ratio(results$stan.fit, pars = names(results$stan.fit)[-c(b.loc,max(b.loc)+1)])
print(mcmc_neff(neff.ratio)+ggtitle("Ratios of Effective Sample Size to Total Sample Size\n(lighter is better)\n")+ 
  theme(plot.title = element_text(hjust = 0.5)))
}
#-----
if(trace){
windows()
print(mcmc_trace(surv.par)+ggtitle("Survival Sub-Model\n")+ 
  theme(plot.title = element_text(hjust = 0.5)))
windows()
print(mcmc_trace(long.par)+ggtitle("Longitudinal Sub-Model\n")+ 
  theme(plot.title = element_text(hjust = 0.5)))
}
#-----
if(density){
windows()
print(mcmc_dens(surv.par)+ggtitle("Survival Sub-Model\n")+ 
  theme(plot.title = element_text(hjust = 0.5)))
windows()
print(mcmc_dens(long.par)+ggtitle("Longitudinal Sub-Model\n")+ 
  theme(plot.title = element_text(hjust = 0.5)))
}
#--------------------
if(uncertainty.intervals){
windows()
print(mcmc_intervals(surv.par)+ggtitle("Survival Sub-Model\n")+ 
  theme(plot.title = element_text(hjust = 0.5)))
windows()
print(mcmc_intervals(long.par)+ggtitle("Longitudinal Sub-Model\n")+ 
  theme(plot.title = element_text(hjust = 0.5)))
}
}
#-----------------