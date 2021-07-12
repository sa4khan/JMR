#--------------
# Compute survival probabilities for one MCMC simulation results
#--------------------
jm.surv0<-function(b0,st0,pred0,xalpha_fixed0,alphatime0,
  xalpha_nodes,nodesr,nodes0,weight,
  phi0,rho0,kappa0,gam0,model,method,lme.model){
if(method=="riz"){
if(model=="weibull"){
if(lme.model=="simple"){
psi = (exp(- pred0 - phi0 * (xalpha_fixed0 + b0[1])) *
  (1 - exp(- phi0 * (alphatime0 + b0[2]) * st0))) / (phi0 * (alphatime0 + b0[2]))
} else{
psi<-t(weight)%*%exp(-phi0*(xalpha_nodes + nodesr%*%b0))*st0*exp(-pred0)/2
}
logs<- -(rho0*psi)^kappa0
} else if(model=="llogistic"){
if(lme.model=="simple"){
psi = (exp(- pred0 - phi0 * (xalpha_fixed0 + b0[1])) *
  (1 - exp(- phi0 * (alphatime0 + b0[2]) * st0))) / (phi0 * (alphatime0 + b0[2]))
} else{
psi<-t(weight)%*%exp(-phi0*(xalpha_nodes + nodesr%*%b0))*st0*exp(-pred0)/2
}
logs<- -log(1+(rho0*psi)^kappa0)
} else if(model=="lnormal"){
if(lme.model=="simple"){
psi = (exp(- pred0 - phi0 * (xalpha_fixed0 + b0[1])) *
  (1 - exp(- phi0 * (alphatime0 + b0[2]) * st0))) / (phi0 * (alphatime0 + b0[2]))
} else{
psi<-t(weight)%*%exp(-phi0*(xalpha_nodes + nodesr%*%b0))*st0*exp(-pred0)/2
}
logs<- plnorm(psi, meanlog = -log(rho0), sdlog = 1/kappa0, lower.tail = FALSE, log.p = TRUE)
} else if(model=="eweibull"){
if(lme.model=="simple"){
psi = (exp(- pred0 - phi0 * (xalpha_fixed0 + b0[1])) *
  (1 - exp(- phi0 * (alphatime0 + b0[2]) * st0))) / (phi0 * (alphatime0 + b0[2]))
} else{
psi<-t(weight)%*%exp(-phi0*(xalpha_nodes + nodesr%*%b0))*st0*exp(-pred0)/2
}
logs<-log1p(-exp(gam0*log1mexp(-(rho0*psi)^kappa0)))
} else if(model=="ggamma"){
if(lme.model=="simple"){
psi = (exp(- pred0 - phi0 * (xalpha_fixed0 + b0[1])) *
  (1 - exp(- phi0 * (alphatime0 + b0[2]) * st0))) / (phi0 * (alphatime0 + b0[2]))
} else{
psi<-t(weight)%*%exp(-phi0*(xalpha_nodes + nodesr%*%b0))*st0*exp(-pred0)/2
}
logs<-pgamma((rho0*psi)^kappa0,shape=gam0,scale=1,lower.tail=FALSE,log.p=TRUE)
} else if(model=="weibullph"){
logs<- -c(weight%*%matrix(exp((kappa0-1)*log(nodes0)+phi0*(xalpha_nodes + nodesr%*%b0)),
   nrow=length(weight)))*st0*kappa0*rho0^kappa0*exp(pred0)/2
} else if(model=="gllph"){
logs<- -c(weight%*%matrix(exp((kappa0-1)*log(nodes0)-log(1+(gam0*nodes0)^kappa0)+
   phi0*(xalpha_nodes + nodesr%*%b0)),nrow=length(weight)))*st0*kappa0*rho0^kappa0*exp(pred0)/2
}
} else{
if(model=="weibull"){
if(lme.model=="simple"){
psi = (exp(- pred0 - phi0 * b0[1]) * (1 - exp(- phi0* b0[2] * st0))) / (phi0 * b0[2])
} else{
psi<-t(weight)%*%exp(-phi0*(nodesr%*%b0))*st0*exp(-pred0)/2
}
logs<--(rho0*psi)^kappa0
} else if(model=="llogistic"){
if(lme.model=="simple"){
psi = (exp(- pred0 - phi0 * b0[1]) * (1 - exp(- phi0* b0[2] * st0))) / (phi0 * b0[2])
} else{
psi<-t(weight)%*%exp(-phi0*(nodesr%*%b0))*st0*exp(-pred0)/2
}
logs<- -log(1+(rho0*psi)^kappa0)
} else if(model=="lnormal"){
if(lme.model=="simple"){
psi = (exp(- pred0 - phi0 * b0[1]) * (1 - exp(- phi0* b0[2] * st0))) / (phi0 * b0[2])
} else{
psi<-t(weight)%*%exp(-phi0*(nodesr%*%b0))*st0*exp(-pred0)/2
}
logs<- plnorm(psi, meanlog = -log(rho0), sdlog = 1/kappa0, lower.tail = FALSE, log.p = TRUE)
} else if(model=="eweibull"){
if(lme.model=="simple"){
psi = (exp(- pred0 - phi0 * b0[1]) * (1 - exp(- phi0* b0[2] * st0))) / (phi0 * b0[2])
} else{
psi<-t(weight)%*%exp(-phi0*(nodesr%*%b0))*st0*exp(-pred0)/2
}
logs<-log1p(-exp(gam0*log1mexp(-(rho0*psi)^kappa0)))
} else if(model=="ggamma"){
if(lme.model=="simple"){
psi = (exp(- pred0 - phi0 * b0[1]) * (1 - exp(- phi0* b0[2] * st0))) / (phi0 * b0[2])
} else{
psi<-t(weight)%*%exp(-phi0*(nodesr%*%b0))*st0*exp(-pred0)/2
}
logs<-pgamma((rho0*psi)^kappa0,shape=gam0,scale=1,lower.tail=FALSE,log.p=TRUE)
} else if(model=="weibullph"){
logs<- -c(weight%*%matrix(exp((kappa0-1)*log(nodes0)+phi0*(nodesr%*%b0)),
   nrow=length(weight)))*st0*kappa0*rho0^kappa0*exp(pred0)/2
} else if(model=="gllph"){
logs<- -c(weight%*%matrix(exp((kappa0-1)*log(nodes0)-log(1+(gam0*nodes0)^kappa0)+
   phi0*(nodesr%*%b0)),nrow=length(weight)))*st0*kappa0*rho0^kappa0*exp(pred0)/2
}
}
return(logs)
}
#--------------
# Posterior density for b given all other parameters
#--------------------
pdensb<-function(b0,y0,xalpha0,zlong0,Sigma0,sigma0,
  st0,pred0,xalpha_fixed0,alphatime0,
  xalpha_nodes,nodesr,nodes0,weight,
  phi0,rho0,kappa0,gam0,model,method,lme.model){
logs<-jm.surv0(b0=b0,st0=st0,pred0=pred0,
 xalpha_fixed0=xalpha_fixed0,alphatime0=alphatime0,
 xalpha_nodes=xalpha_nodes,nodesr=nodesr,nodes0=nodes0,weight=weight,
 phi0=phi0,rho0=rho0,kappa0=kappa0,gam0=gam0,
 model=model,method=method,lme.model=lme.model)
y.lden<-sum(dnorm(y0,mean=(xalpha0+zlong0%*%b0),sd=sigma0,log=TRUE))
b.lden<-den.mnorm(b0,mean=rep(0,length(b0)),sigma=Sigma0,log=TRUE)
lden<-y.lden+logs+b.lden
return(-lden)
}
########################################################################
#' Prediction of survival probabilities from a joint model fit
#' @description Dynamic predictions of survival probabilites as described by Rizopoulos (2016, 2020).
#' @keywords Joint moldeing, MCMC, Prediction, Survival probabilities
#' @param jmfit a \code{\link{jm.reg}} fit with \code{method = "MCMC"}.
#' @param newdata a data frame in which to look for variables with which to predict.
#'     The data frame must have columns that correspond by name to the
#'     original variables used to fit the joint model. 
#' @param n.sample an integer denoting how many MCMC samples to use. Default is \code{min(200, n.iter)}, where
#'      n.iter is the number of MCMC iterations (for all chains combined). 
#' @param st a numeric vector of times at which survival probabilities 
#'        P(\eqn{T > t | T > s}, history of the longitudinal response, covariate information) 
#'        are to be computed (see Details).
#'        If \code{NULL}, \eqn{s} is taken to be the last time the 
#'        subject provided a longitudinal measurement, and \code{t = seq(s, max.st, length = 15)}, where
#'        \code{max.st} is \code{floor(max(event times))}. A warning will be given If \eqn{s >} 
#'        the last time the individual provided a longitudinal measurement.
#' @param plot if \code{TRUE}, a plot of the survival curve is produced.
#' @param posterior.mean if \code{TRUE}, posterior means of the survival probabilities
#'          are used to produce the survival curve, otherwise posterior medians are used. 
#' @param include.y if \code{TRUE}, the fitted longitudinal profile is included 
#'          in the survival curve plot.
#' @param xlim x limits (a vector of length 2 of the form \code{c(minimum, maximum)}) of the plot (optional).
#' @param ylim limits for the left y axis (a vector of length 2 of the form \code{c(minimum, maximum)}) of the plot (optional).
#'        If \code{include.y = TRUE}, the left y axis is for the the fitted longitudinal profile, otherwise
#'        it is for predictions of survival probabilites.        
#' @param xlab a label (character expression) for the x axis (optional).
#' @param ylab.left a label (character expression) for the left y axis (optional). If \code{include.y = TRUE}, 
#'        the left y axis is for the the fitted longitudinal profile, otherwise
#'        it is for predictions of survival probabilites.
#' @param ylab.right a label (character expression) for the right y axis (optional). If \code{include.y = TRUE},
#'        the right y axis is for predictions of survival probabilites, otherwise
#'        it is ignored.
#' @param main main title (character expression) of the plot (optional).
#' @param adapt_delta the target average proposal acceptance probability during Stan’s adaptation period
#'        (see Details).
#' @param max_treedepth a positive integer specifying the maximum treedepth (see the Stan manual).
#' @param warmup a positive integer specifying the number of burnin iterations (see Details).
#' @param seed the seed for random number generation.
#' @details The algorithm proposed by Rizopoulos (2020) is used for dynamic predictions.
#'        The history of the longitudinal response up to time \eqn{s} along with covariate 
#'        information are taken into consideration for predictions.
#'        The algorithm is described as follows (see also Rizopoulos (2016, 2020)).
#'     \enumerate{ 
#'     \item Randomly take a realization from the posterior simulations of the joint model parameters. 
#'           Denote it by \eqn{\theta^{(1)}}. 
#'     \item Given \eqn{\theta^{(1)}}, the history of the longitudinal response, covariate information, 
#'           and survival up to time \eqn{s}, draw a realization of the random-effects vector \eqn{b} from its posterior distribution.
#'           For this step, the MCMC algorithm (implemented in Stan) is used with \eqn{n=1}, 
#'           \code{event time} = \eqn{s} and \code{status} = 0. The number of burnin iterations, adapt_delta
#'           and max_treedepth in Stan can be controlled using the arguments \code{warmup} (default is 200),
#'           \code{adapt_delta} (default is 0.90) and \code{max_treedepth} (default is 15).  
#'           Note that Rizopoulos (2020) used the Metropolis-Hastings algorithm with independent proposals (a 
#'           properly centered and scaled multivariate t distribution) to estimate \eqn{b}.
#'     \item Using \eqn{\theta^{(1)}} and \eqn{b}, compute P(\eqn{T > t | T > s}) = P(\eqn{T > t})/P(\eqn{T > s}).
#'     \item Repeat steps 1-3 \code{n.sample} times.   
#'          }
#'     Posterior summeries of the conditional probabilities are derived from the MCMC sample drawn 
#'     using the above algorithm.
#' @return Posterior summaries of the predictions of survival probabilites are returned.
#'         A plot of the survival curve is produced if \code{plot = TRUE}, with
#'         the predictions displayed by the solid line and the pointwise Bayesian intervals displayed by the dashed lines.
#' @author Shahedul Khan <khan@math.usask.ca>
#' @references Rizopoulos D, JMbayes: Joint Modeling of Longitudinal and Time-to-Event Data under a Bayesian Approach, 
#'           R package version 0.8-85, 2020, \url{https://cran.r-project.org/web/packages/JMbayes}.
#' @references Rizopoulos D, The R package JMbayes for fitting joint models for longitudinal and time-to-event data Using MCMC, 
#'           Journal of Statistical Software, 72(7): 1-45, 2016.
#' @seealso \code{\link{jm.lppredict}}, \code{\link{jm.lspredict}}, \code{\link{jm.reg}}
#' @examples
#'   # Example: AIDS data from package 'JM'
#'   library(JM)
#'   data(aids.id)
#'   data(aids)
#'   surv.fit <- coxph(Surv(Time, death) ~ drug + gender + prevOI + AZT,
#'           data = aids.id, x = TRUE)
#'   lme.fit <- lme(CD4 ~ obstime + obstime:drug + gender + prevOI + AZT,
#'       random =  ~ obstime | patient, data = aids)
#'   jmfit.ew0 <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "eweibull",  timevar = "obstime", form = "riz")
#'   newdata <- data.frame(aids[aids$patient == 409, ])
#'   jm.surv(jmfit.ew0, newdata = newdata, posterior.mean = FALSE)
#' @export

jm.surv<-function(jmfit,newdata,n.sample=NULL,
  st=NULL,plot=TRUE,posterior.mean=TRUE,include.y=TRUE,
  xlim=NULL,ylim=NULL,xlab=NULL,ylab.left=NULL,ylab.right=NULL,main=NULL,
  adapt_delta=0.9,max_treedepth=15,warmup=200,
  seed=sample.int(.Machine$integer.max, 1)){
start_time <- Sys.time()
if(class(jmfit)!="JMR"){
stop("jmfit must inherit from class JMR")
}
results<-jmfit$all.results
if(is.null(results)){
stop("\nThis function works only for method = 'MCMC'")
}
posterior<-as.matrix(results$stan.fit)
n.iter<-nrow(posterior)
if(is.null(n.sample)){
n.sample<-min(200,n.iter)
}
con<-list(adapt_delta=adapt_delta,max_treedepth=max_treedepth,warmup=warmup,
   n.sample=n.sample,qpoints=15)
n.sample<-con$n.sample
set.seed(seed)
# The new dataset must be a data frame
if(!is.data.frame(newdata)){
newdata<-data.frame(newdata)
}
# Extraact information
surv.model<-results$surv.model
method<-results$method
p1<-results$stan.data$p1
pp1<-results$stan.data$pp1
timevar<-results$timevar
surv.data<-results$surv.data
surv.fit<-results$surv.fit
lme.fit<-results$lme.fit
dat<-lme.fit$data
id.name<-colnames(lme.fit$groups)
if(length(unique(newdata[[id.name]]))>1)
stop(paste("\nThere are",length(unique(newdata[[id.name]])),"subjects in newdata.",
  "\npredictions are allowed for one subject/id at a time.",collapse="")) 
# Design matrix for the survival submodel
newsurvdata<-newdata[1,]
formula.surv <- formula(surv.fit)
survvar<-all.vars(formula.surv)[-(1:2)]
if(isFALSE(all(survvar %in% colnames(newsurvdata)))){
stop("\nVariable names in newdata must be
 the same as in the datasets that were used to fit the survival 
 and the longitudinal models. Not all variables for the 
 survival model are found in newdata.")
}
xsurv<-matrix(model.matrix(formula.surv, newsurvdata)[,-1],nrow=1)
# Design matrix for longitudinal submodel (fixed effects)
formula.fixed<-formula(lme.fit)
longvar<-all.vars(formula.fixed)
if(isFALSE(all(longvar %in% colnames(newdata)))){
stop("\nVariable names in newdata must be
 the same as in the datasets that were used to fit the survival 
 and the longitudinal models. Not all variables for the fixed 
 part of the longitudinal model are found in newdata.")
}
mframe.fdata<-model.frame(terms(formula.fixed),data=dat)
inf.fdata<-attr(mframe.fdata,"terms")
model.fdata<-model.frame(inf.fdata, data = newdata)
y<-as.vector(model.extract(model.fdata, "response"))
xlong<-matrix(model.matrix(formula.fixed,model.fdata),ncol=p1)
# Design matrix for longitudinal submodel (random effects)
formula.rand<-formula(lme.fit$modelStruct$reStruct[[1]])
randvar<-all.vars(formula.rand)
if(isFALSE(all(randvar %in% colnames(newdata)))){
stop("\nVariable names in newdata must be
 the same as in the datasets that were used to fit the survival 
 and the longitudinal models. Not all variables for the random 
 part of the longitudinal model are found in newdata.")
}
mframe.rdata<-model.frame(terms(formula.rand),data = dat)
inf.rdata<-attr(mframe.rdata,"terms")
model.rdata<-model.frame(inf.rdata, data = newdata)
zlong<-matrix(model.matrix(formula.rand,model.rdata),ncol=pp1)
ni<-nrow(zlong)
# st for prediction
if(!is.null(st)){
st<-sort(st)
if(st[1]>newdata[nrow(newdata),timevar])
warning(paste("min(st) = ", st[1], " > last time the individual 
     provided a longitudinal measurement",sep=""))
}
if(is.null(st)){
st0<-newdata[nrow(newdata),timevar]
max.st<-max(surv.data[,"st"])
st<-seq(st0,floor(max.st),length=15)
}
# location of the parameters in the MCMC samples
beta.loc<-results$beta.loc
phi.loc<-results$phi.loc
kappa.loc<-results$kappa.loc
rho.loc<-results$rho.loc
lcoef.loc<-results$alpha.loc
lcov.loc<-results$cov.loc
lres.loc<-results$sigma.loc
reffects.loc<-results$b.loc
if(surv.model=="eweibull" || surv.model=="ggamma" || surv.model=="gllph"){
gam.loc<-results$gam.loc
}
sam.dat<-posterior[sample(n.iter,n.sample),]
if(results$surv.model=="weibull" || results$surv.model=="llogistic" || 
 results$surv.model=="lnormal" || results$surv.model=="eweibull" || results$surv.model=="ggamma"){
lme.model<-ifelse(method=="riz",ifelse((results$fixed.model=="simple" && results$rand.model=="simple"),
 "simple","nsimple"),ifelse(results$rand.model=="simple","simple","nsimple"))
} else{
lme.model<-"nsimple"
}
# information for numerical integration
qpoints<-con$qpoints
nodes0<-GK.nodes(qpoints)$gk.n
weight<-GK.nodes(qpoints)$gk.w
nodes1<-sapply(st,function(st){0.5*(st*nodes0+st)})
newlongdata1<-newdata[1,]
newlongdata2<-newlongdata1[rep(seq_len(nrow(newlongdata1)), length(st)*qpoints),]
newlongdata2[[timevar]]<-c(nodes1)
newlongdata3<-model.frame(inf.rdata,data=newlongdata2)
newlongdata4<-model.frame(inf.fdata,data=newlongdata2)
nodes<-matrix(model.matrix(formula.rand,newlongdata3),ncol=pp1)
nodesf<-matrix(model.matrix(formula.fixed,newlongdata4),ncol=p1)
nodes.array<-list()
nodesf.array<-list()
for(i in 1:length(st)){
nodes.array[[i]]<-nodes[(1+(i-1)*qpoints):(i*qpoints),]
nodesf.array[[i]]<-nodesf[(1+(i-1)*qpoints):(i*qpoints),]
}
# create data for stan
beta<-sam.dat[,beta.loc]
phi<-sam.dat[,phi.loc]
rho<-sam.dat[,rho.loc]
kappa<-sam.dat[,kappa.loc]
if(surv.model=="eweibull" || surv.model=="ggamma" || surv.model=="gllph"){
gam<-sam.dat[,gam.loc]
} else{
gam<-NULL
}
alpha<-sam.dat[,lcoef.loc]
sigma<-sam.dat[,lres.loc]
pred0<-c(xsurv%*%t(beta))
if(ni>1){
y.mat<-t(replicate(n.sample,y))
} else{
y.mat<-matrix(rep(y,n.sample),ncol=1)
}
if(lme.model=="simple" && method=="riz"){
stan.dat<-results$stan.data
which_alphafixed<-stan.dat$which_alphafixed
which_alphatime<-stan.dat$which_alphatime
alphatime<-alpha[,which_alphatime]
alphafixed<-alpha[,which_alphafixed]
if(is.vector(alphafixed)){
alphafixed<-matrix(alphafixed,ncol=1)
}
xlongfixed<-xlong[,-which_alphatime]
if(is.vector(xlongfixed)){
xlongfixed<-matrix(xlongfixed,ncol=1)
}
xalpha_fixed<-c(matrix(xlongfixed[1,],nrow=1)%*%t(alphafixed))
} else{
alphatime<-alphafixed<-xalpha_fixed<-NULL
}
Sigma<-list()
Sigma.chol<-list()
xalpha<-matrix(0,nrow=n.sample,ncol=nrow(xlong))
xalpha_nodes<-matrix(0,nrow=n.sample,ncol=qpoints)
u<-matrix(0,nrow=sqrt(length(lcov.loc)),ncol=n.sample)
init.u<-rep(0,pp1)
for(i in 1:n.sample){
Sigma[[i]]<-matrix(sam.dat[i,lcov.loc],sqrt(length(lcov.loc)))
Sigma.chol[[i]]<-t(chol(Sigma[[i]]))
xalpha[i,]<-c(xlong%*%alpha[i,])
xalpha_nodes[i,]<-c(nodesf.array[[1]]%*%alpha[i,])
options(warn=-1)
fit.b<-tryCatch({nlminb(start=init.u,objective=pdensb,
 y0=y.mat[i,],xalpha0=xalpha[i,],zlong0=zlong,Sigma0=Sigma[[i]],sigma0=sigma[i],
  st0=st[1],pred0=pred0[i],xalpha_fixed0=xalpha_fixed[i],
  alphatime0=alphatime[i],xalpha_nodes=xalpha_nodes[i,],nodesr=nodes.array[[1]],
  nodes0=nodes1[,1],weight=weight,
  phi0=phi[i],rho0=rho[i],kappa0=kappa[i],gam0=gam[i],
  model=surv.model,method=method,lme.model=lme.model)},
  error = function(e){NULL})
options(warn=0)
if(!is.null(fit.b)){
if((fit.b$message=="relative convergence (4)" || 
  fit.b$message=="both X-convergence and relative convergence (5)") &&
  fit.b$convergence==0 && is.finite(fit.b$objective) && 
  fit.b$objective!=0 && !is.na(fit.b$objective)){
u[,i]<-solve(Sigma.chol[[i]])%*%fit.b$par
} else{
u[,i]<-init.u
}
}
if(is.null(fit.b)){
u[,i]<-init.u
}
}
inits<-list(list(b_unscaled=u))
par.stan<-c("b")
#--- setup stan model---------
if(method=="riz"){
if(surv.model=="weibull"){
if(lme.model=="simple"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   xalpha_fixed=xalpha_fixed,alphatime=alphatime,phi=phi,
   rho=rho,kappa=kappa,chol_cov_matrix=Sigma.chol)
#stanc("stan_weibullsurvrizsimple.stan")
#stan_model <- "stan_weibullsurvrizsimple.stan"
stan_model <- stanmodels$weibullsurvrizsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,xalpha1=xalpha_nodes,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_weibullsurvriz.stan")
#stan_model <- "stan_weibullsurvriz.stan"
stan_model <- stanmodels$weibullsurvriz
}
} else if(surv.model=="llogistic"){
if(lme.model=="simple"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   xalpha_fixed=xalpha_fixed,alphatime=alphatime,phi=phi,
   rho=rho,kappa=kappa,chol_cov_matrix=Sigma.chol)
#stanc("stan_llogisticsurvrizsimple.stan")
#stan_model <- "stan_llogisticsurvrizsimple.stan"
stan_model <- stanmodels$llogisticsurvrizsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,xalpha1=xalpha_nodes,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_llogisticsurvriz.stan")
#stan_model <- "stan_llogisticsurvriz.stan"
stan_model <- stanmodels$llogisticsurvriz
} 
} else if(surv.model=="lnormal"){
if(lme.model=="simple"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   xalpha_fixed=xalpha_fixed,alphatime=alphatime,phi=phi,
   rho=rho,kappa=kappa,chol_cov_matrix=Sigma.chol)
#stanc("stan_lnormalsurvrizsimple.stan")
#stan_model <- "stan_lnormalsurvrizsimple.stan"
stan_model <- stanmodels$lnormalsurvrizsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,xalpha1=xalpha_nodes,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_lnormalsurvriz.stan")
#stan_model <- "stan_lnormalsurvriz.stan"
stan_model <- stanmodels$lnormalsurvriz
}
} else if(surv.model=="eweibull"){
if(lme.model=="simple"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   xalpha_fixed=xalpha_fixed,alphatime=alphatime,phi=phi,
   rho=rho,kappa=kappa,gam=gam,chol_cov_matrix=Sigma.chol)
#stanc("stan_eweibullsurvrizsimple.stan")
#stan_model <- "stan_eweibullsurvrizsimple.stan"
stan_model <- stanmodels$eweibullsurvrizsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,xalpha1=xalpha_nodes,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,gam=gam,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_eweibullsurvriz.stan")
#stan_model <- "stan_eweibullsurvriz.stan"
stan_model <- stanmodels$eweibullsurvriz
}
} else if(surv.model=="ggamma"){
if(lme.model=="simple"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   xalpha_fixed=xalpha_fixed,alphatime=alphatime,phi=phi,
   rho=rho,kappa=kappa,gam=gam,chol_cov_matrix=Sigma.chol)
#stanc("stan_ggammasurvrizsimple.stan")
#stan_model <- "stan_ggammasurvrizsimple.stan"
stan_model <- stanmodels$ggammasurvrizsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,xalpha1=xalpha_nodes,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,gam=gam,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_ggammasurvriz.stan")
#stan_model <- "stan_ggammasurvriz.stan"
stan_model <- stanmodels$ggammasurvriz
}
} 
else if(surv.model=="weibullph"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,xalpha1=xalpha_nodes,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,weights=weight,
   nodes0=nodes1[,1],nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_weibullphsurvriz.stan")
#stan_model <- "stan_weibullphsurvriz.stan"
stan_model <- stanmodels$weibullphsurvriz
} 
else if(surv.model=="gllph"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,xalpha1=xalpha_nodes,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,gam=gam,weights=weight,
   nodes0=nodes1[,1],nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_gllphsurvriz.stan")
#stan_model <- "stan_gllphsurvriz.stan"
stan_model <- stanmodels$gllphsurvriz
}
} else{
if(surv.model=="weibull"){
if(lme.model=="sample"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   phi=phi,rho=rho,kappa=kappa,chol_cov_matrix=Sigma.chol)
#stanc("stan_weibullsurvsimple.stan")
#stan_model <- "stan_weibullsurvsimple.stan"
stan_model <- stanmodels$weibullsurvsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_weibullsurv.stan")
#stan_model <- "stan_weibullsurv.stan"
stan_model <- stanmodels$weibullsurv
}
} else if(surv.model=="llogistic"){
if(lme.model=="sample"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   phi=phi,rho=rho,kappa=kappa,chol_cov_matrix=Sigma.chol)
#stanc("stan_llogisticsurvsimple.stan")
#stan_model <- "stan_llogisticsurvsimple.stan"
stan_model <- stanmodels$llogisticsurvsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_llogisticsurv.stan")
#stan_model <- "stan_llogisticsurv.stan"
stan_model <- stanmodels$llogisticsurv
}
} else if(surv.model=="lnormal"){
if(lme.model=="sample"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   phi=phi,rho=rho,kappa=kappa,chol_cov_matrix=Sigma.chol)
#stanc("stan_lnormalsurvsimple.stan")
#stan_model <- "stan_lnormalsurvvsimple.stan"
stan_model <- stanmodels$lnormalsurvsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_lnormalsurv.stan")
#stan_model <- "stan_lnormalsurv.stan"
stan_model <- stanmodels$lnormalsurv
}
} else if(surv.model=="eweibull"){
if(lme.model=="sample"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   phi=phi,rho=rho,kappa=kappa,gam=gam,chol_cov_matrix=Sigma.chol)
#stanc("stan_eweibullsurvsimple.stan")
#stan_model <- "stan_eweibullsurvsimple.stan"
stan_model <- stanmodels$eweibullsurvsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,gam=gam,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_eweibullsurv.stan")
#stan_model <- "stan_eweibullsurv.stan"
stan_model <- stanmodels$eweibullsurv
}
} else if(surv.model=="ggamma"){
if(lme.model=="sample"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   phi=phi,rho=rho,kappa=kappa,gam=gam,chol_cov_matrix=Sigma.chol)
#stanc("stan_ggammasurvsimple.stan")
#stan_model <- "stan_ggammasurvsimple.stan"
stan_model <- stanmodels$ggammasurvsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,gam=gam,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_ggammasurv.stan")
#stan_model <- "stan_ggammasurv.stan"
stan_model <- stanmodels$ggammasurv
}
} else if(surv.model=="weibullph"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,weights=weight,
   nodes0=nodes1[,1],nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_weibullphsurv.stan")
#stan_model <- "stan_weibullphsurv.stan"
stan_model <- stanmodels$weibullphsurv
} else if(surv.model=="gllph"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=st[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,gam=gam,weights=weight,
   nodes0=nodes1[,1],nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_gllphsurv.stan")
#stan_model <- "stan_gllphsurv.stan"
stan_model <- stanmodels$gllphsurv
}
}
#---- Run stan -----------------
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
stan_fit <- quiet(suppressWarnings(rstan::sampling(object=stan_model, data = data.stan, 
  init = inits, pars = par.stan,
  warmup = con$warmup, iter = con$warmup+1, 
  chains = 1, thin = 1, seed=seed,
  control=list(adapt_delta=con$adapt_delta,max_treedepth=con$max_treedepth))))
posterior.b<-as.array(stan_fit)
b<-matrix(posterior.b[sample(con$warmup+1-con$warmup,1),1,][-dim(posterior.b)[3]],nrow=n.sample)
# -- Compute survival probabilities ------------
ssam<-matrix(0,nrow=n.sample,ncol=length(st))
ssam1<-matrix(0,nrow=n.sample,ncol=ncol(xalpha))
for(i in 1:n.sample){
s1<-rep(0,length(st))
for(j in 1:length(st)){
s1[j]<-jm.surv0(b0=b[i,],st0=st[j],pred0=pred0[i],
  xalpha_fixed0=xalpha_fixed[i],alphatime0=alphatime[i],
  xalpha_nodes=c(nodesf.array[[j]]%*%alpha[i,]),
  nodesr=nodes.array[[j]],nodes0=nodes1[,j],weight=weight,
  phi0=phi[i],rho0=rho[i],kappa0=kappa[i],gam0=gam[i],
  model=surv.model,method=method,lme.model=lme.model)
}
ssam[i,]<-exp(s1-s1[1])
ssam1[i,]<-xalpha[i,]+c(zlong%*%b[i,])
}
fresults<-cbind(t=st,mean=colMeans(ssam,na.rm=TRUE),sd=apply(ssam,2,stats:::sd,na.rm=TRUE),
      "2.5%"=apply(ssam,2,quantile,prob=0.025,na.rm=TRUE),
      "25%"=apply(ssam,2,quantile,prob=0.25,na.rm=TRUE),
      "50%"=apply(ssam,2,quantile,prob=0.5,na.rm=TRUE),
      "75%"=apply(ssam,2,quantile,prob=0.75,na.rm=TRUE),
      "97.5%"=apply(ssam,2,quantile,prob=0.975,na.rm=TRUE))
rownames(fresults)<-rep("",nrow(fresults))
resultname<-paste("Posterior summeries of P(T > t | T > ", round(st[1],digits=3),"):",collapse="")
ftime<-newdata[[timevar]]
if(plot){
if(posterior.mean){
pred.y<-colMeans(ssam1,na.rm=TRUE)
pred.s<-fresults[,2]
lower.s<-fresults[,4]
upper.s<-fresults[,8]
} else{
pred.y<-apply(ssam1,2,quantile,prob=0.5,na.rm=TRUE)
pred.s<-fresults[,6]
lower.s<-fresults[,4]
upper.s<-fresults[,8]
}
if(include.y){
if(is.null(ylim)){
ylim0<-c(min(getResponse(lme.fit)),max(getResponse(lme.fit)))
} else{
ylim0<-ylim
}
if(is.null(xlim)){
xlim0<-c(min(ftime),max(st))
} else{
xlim0<-xlim
}
if(is.null(xlab)){
xlab<-"Time"
} else{
xlab<-xlab
}
if(is.null(ylab.left)){
ylab.left<-formula.fixed[[2]]
} else{
ylab.left<-ylab.left
}
time1<-c(ftime,st)
y1<-c(y,rep(NA,length(st)))
pred.s1<-c(rep(NA,length(ftime)),pred.s)
lower.s1<-c(rep(NA,length(ftime)),lower.s)
upper.s1<-c(rep(NA,length(ftime)),upper.s)
par(mar=c(5, 4, 4, 6) + 0.1)
plot(time1,y1,pch=19,axes=FALSE,ylim=ylim0,xlim=xlim0,xlab="", ylab="",type="p")
lines(ftime,pred.y,lwd=1)
abline(v=max(ftime),lty=2,lwd=3)
axis(2, ylim=ylim0,las=1) 
mtext(ylab.left,side=2,line=2.75,cex=1.15)
box()
par(new=TRUE)
plot(time1,pred.s1,xlab="",ylab="",ylim=c(0,1),xlim=xlim0,axes=FALSE, type="l",lwd=1)
lines(time1,lower.s1,lty=3)
lines(time1,upper.s1,lty=3)
axis(4, ylim=c(0,1),las=1)
if(is.null(ylab.right)){
ylab.right<-list(bquote("P( T > t | T >" ~ .(st0)~")"))
mtext(do.call(expression, ylab.right),side=4,line=3.75,cex=1.15) 
} else{
mtext(ylab.right,side=4,line=3.75,cex=1.15) 
}
axis(1,pretty(range(time1),10))
mtext(xlab,side=1,col="black",line=3,cex=1.15)
if(is.null(main)){
subject.id<-newdata[[id.name]][1]
title<-list(bquote("Subject" ~ .(subject.id)))
mtext(do.call(expression, title),side=3,cex=1.5,line=0.5)
} else{
mtext(main,side=3,cex=1.5,line=0.5)
}
} else{
if(is.null(ylim)){
ylim0<-c(0,1)
} else{
ylim0<-ylim
}
if(is.null(xlim)){
xlim0<-c(min(st),max(st))
} else{
xlim0<-xlim
}
if(is.null(xlab)){
xlab<-"Time"
} else{
xlab<-xlab
}
par(mar=c(5, 4.5, 4, 3)+0.1)
plot(st,pred.s,type="l",axes=FALSE,ylab="",xlab=xlab,lwd=2,cex.lab=1.15,ylim=ylim0,xlim=xlim0)
lines(st,lower.s,lty=3)
lines(st,upper.s,lty=3)
if(is.null(ylab.left)){
ylab0<-list(bquote("P( T > t | T >" ~ .(st0)~")"))
mtext(do.call(expression, ylab0),side=2,line=3,cex=1.15) 
} else{
mtext(ylab.left,side=2,line=3,cex=1.15) 
}
axis(2, ylim=ylim0,las=1)
axis(1)
if(is.null(main)){
subject.id<-newdata[[id.name]][1]
title<-list(bquote("Subject" ~ .(subject.id)))
mtext(do.call(expression, title),side=3,cex=1.5,line=0.5)
} else{
mtext(main,side=3,cex=1.5,line=0.5)
}
box()
}
}
end_time <- Sys.time()
output<-list()
output[[resultname]]<-fresults
output[["time"]]<-noquote(paste(round(difftime(end_time,start_time,units="mins")[[1]],
  digits=2),"mins",collapse=""))
class(output)<-"JMR"
attr(output, "hidden")<-NULL
return(output)
}
########################################################################
#' Marginal predictions for the longitudinal outcome of a joint model
#' @description Given a \code{\link{jm.reg}} fit and a data frame, this function
#'      produces marginal predictions (population-level) for the longitudinal outcome. 
#' @keywords Joint modeling, MCMC, Prediction, Longitudinal
#' @param jmfit a \code{\link{jm.reg}} fit with \code{method = "MCMC"}.
#' @param newdata a data frame in which to look for variables with which to predict.
#'     The data frame must have columns that correspond by name to the
#'     original regressors used in the longitudinal part of the joint model fit.
#' @details Posterior summaries of X(s)\eqn{\alpha} (marginal predictions)
#'     are obtained from the MCMC samples, where X(s) is the fixed-effects
#'     design matrix and \eqn{\alpha} is the fixed-effects
#'     coefficient vector for the longitudinal model. 
#' @return This function returns posterior summaries of X(s)\eqn{\alpha} (marginal predictions).
#' @author Shahedul Khan <khan@math.usask.ca>
#' @seealso \code{\link{jm.lspredict}}, \code{\link{jm.surv}}, \code{\link{jm.reg}}
#' @examples 
#'   # Example: pbc data
#'   lme.fit <- lme(log(bilirubin) ~ drug + ns(futime, 2),
#'     data = pbc.long, random = ~ ns(futime, 2) | id)
#'   surv.fit <- coxph(Surv(st, status2) ~ drug * age, data = pbc.surv, x = TRUE)
#'   jmfit.gg2 <-  jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "ggamma", timevar = "futime", form = "riz")
#'   jm.summary(jmfit.gg2)
#'   newdata <-  expand.grid(futime = seq(min(pbc.long$futime), max(pbc.long$futime), 
#'       len = 25), drug = 0:1)
#'   pred <- jm.lppredict(jmfit = jmfit.gg2, newdata = newdata)
#'   pred
#'   # Plot of the predicted values
#'   futime <- pred[[1]][, 1]
#'   post.mean <- pred[[1]][, 2]
#'   lower <- pred[[1]][, 4]
#'   upper <- pred[[1]][, 8]
#'   drug <- factor(newdata$drug, levels = c(0, 1), labels = c("placebo", "D-penicillamine"))
#'   # Use the xyplot function from package 'lattice'
#'   library(lattice)
#'   xyplot(post.mean + lower + upper ~ futime | drug,
#'     type = "l", col = rep(1, 3), lty = c(1, 3, 3),
#'     lwd = 2,  ylab = "log(bilirubin)",  xlab = "futime")
#' @export

jm.lppredict<-function(jmfit,newdata){
if(class(jmfit)!="JMR"){
stop("jmfit must inherit from class JMR")
}
if(!is.data.frame(newdata)){
newdata<-data.frame(newdata)
}
results<-jmfit$all.results
p1<-results$stan.data$p1
timevar<-results$timevar
lme.fit<-results$lme.fit
dat<-lme.fit$data
posterior<-as.matrix(results$stan.fit)
#n.iter<-nrow(posterior)
#surv.model<-results$surv.model
#method<-results$method
#pp1<-results$stan.data$pp1
#surv.data<-results$surv.data
#surv.fit<-results$surv.fit
#id.name<-colnames(lme.fit$groups)
#if(isFALSE(id.name %in% colnames(newdata))){
#stop(paste("\nsubject identifier variable",id.name,"is not found in newdata",collapse=""))
#}
formula.fixed<-formula(lme.fit)
longvar<-all.vars(formula.fixed)
if(isFALSE(longvar[1] %in% colnames(newdata))){
newlongdataY<-data.frame(1,newdata)
colnames(newlongdataY)[1]<-longvar[1]
} else{
newlongdataY<-newdata
}
if(isFALSE(all(longvar %in% colnames(newlongdataY)))){
stop("\nVariable names in newdata must be
 the same as in the dataset that was used to fit the 
 longitudinal model. Not all variables for the 
 longitudinal model are found in newdata.")
}
mframe.fdata<-model.frame(terms(formula.fixed),data=dat)
inf.fdata<-attr(mframe.fdata,"terms")
model.fdata<-model.frame(inf.fdata, data = newlongdataY)
xlong<-matrix(model.matrix(formula.fixed,model.fdata),ncol=p1)
alpha<-posterior[,results$alpha.loc]
pop.pred<-apply(alpha,1,function(v) xlong%*% v)
quan<-apply(pop.pred,1,quantile,prob=c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)
pop.results<-cbind(t=newdata[[timevar]],
  mean=rowMeans(pop.pred,na.rm=TRUE),sd=apply(pop.pred,1,stats:::sd,na.rm=TRUE),
  "2.5%"=quan[1,],"25%"=quan[2,],"50%"=quan[3,],"75%"=quan[4,],"97.5%"=quan[5,])
#newlongdata.sp<-split(newlongdataY,newlongdataY[[id.name]])
#alpha<-posterior[,results$alpha.loc]
#pop.results<-NULL
#pop.results.list<-list()
#for(i in 1:length(newlongdata.sp)){
#split.dat<-newlongdata.sp[[i]]
#model.fdata<-model.frame(inf.fdata, data = split.dat)
#xlong<-matrix(model.matrix(formula.fixed,model.fdata),ncol=p1)
#pop.pred<-apply(alpha,1,function(v) xlong%*% v)
#pop.results0<-cbind(id=split.dat[[id.name]],t=split.dat[[timevar]],
#  mean=rowMeans(pop.pred,na.rm=TRUE),sd=apply(pop.pred,1,stats:::sd,na.rm=TRUE),
#  "2.5%"=apply(pop.pred,1,quantile,prob=0.025,na.rm=TRUE),
#  "25%"=apply(pop.pred,1,quantile,prob=0.25,na.rm=TRUE),
#  "50%"=apply(pop.pred,1,quantile,prob=0.5,na.rm=TRUE),
#  "75%"=apply(pop.pred,1,quantile,prob=0.75,na.rm=TRUE),
#  "97.5%"=apply(pop.pred,1,quantile,prob=0.975,na.rm=TRUE))
#pop.results<-rbind(pop.results,pop.results0)
#pop.results.list[[i]]<-pop.results0
#}
#colnames(pop.results)[1]<-id.name
#colnames(pop.results)[2]<-timevar
#if(plot){
#if(posterior.mean){
#pred1<-pop.results[,3]
#lower<-pop.results[,5]
#upper<-pop.results[,9]
#} else{
#pred1<-pop.results[,7]
#lower<-pop.results[,5]
#upper<-pop.results[,9]
#}
#ylim0<-c(min(c(pred1,lower,upper)),max(c(pred1,lower,upper)))
#xlim0<-c(min(pop.results[,2]),max(pop.results[,2]))
#ylim00<-c(min(pred1),max(pred1))
#subject<-factor(pop.results[,1],labels=sprintf("subject %d",unique(pop.results[,1])))
#if(length(pop.results.list)>1){
#if(!posterior.mean){
#par(mar=c(5, 4.5, 4, 8.5) + 0.1,xpd=TRUE)
#plot(pop.results.list[[1]][,2],pop.results.list[[1]][,7],type="l",
#   xlim=xlim0,ylim=ylim00, 
#   ylab=formula.fixed[[2]],xlab="Time",cex.lab=1.15,lty=1)
#for(i in 2:length(pop.results.list)){
#lines(pop.results.list[[i]][,2],pop.results.list[[i]][,7],lty=i)
#}
#legend("topright",inset=c(-0.35,0),lty=1:length(pop.results.list),
#  legend=sprintf("subject %d",unique(pop.results[,1])))
#} else{
#par(mar=c(5, 4.5, 4, 4) + 0.1)
#plot(pop.results.list[[1]][,2],pop.results.list[[1]][,3],type="l",
#   xlim=xlim0,ylim=ylim00, 
#   ylab=formula.fixed[[2]],xlab="Time",cex.lab=1.15,lty=1)
#for(i in 2:length(pop.results.list)){
#lines(pop.results.list[[i]][,2],pop.results.list[[i]][,3],lty=i)
#}
#legend("topleft",lty=1:length(pop.results.list),legend=sprintf("subject %d",unique(pop.results[,1])))
#}
#}
#windows()
#stripParams <- list(cex=1.25, lines=1.5)
#print(xyplot(pred1 + lower + upper ~ pop.results[,2] | subject,
#  xlim=xlim0,ylim=ylim0, 
#  type = "l", col = rep(1,3), lty = c(1,3,3), lwd = 2,
#  ylab = list(label=formula.fixed[[2]],fontsize=14),
#  xlab=list(label="Time",fontsize=14),
#  par.strip.text = stripParams,
#  par.settings = list(strip.background=list(col="gray"))))
#}
output<-list()
resultname<-paste("posterior summaries of the predicted values")
output[[resultname]]<-pop.results
class(output)<-"JMR"
attr(output, "hidden")<-NULL
return(output)
}
########################################################################
#' Subject-specific predictions for the longitudinal outcome of a joint model
#' @description Given a \code{\link{jm.reg}} fit and a data frame, this function
#'      produces subject-specific predictions for the longitudinal outcome.
#' @keywords Joint modeling, MCMC, Prediction, Longitudinal
#' @param jmfit a \code{\link{jm.reg}} fit with \code{method = "MCMC"}.
#' @param newdata a data frame in which to look for variables with which to predict.
#'     The data frame must have columns that correspond by name to the
#'     original variables used to fit the joint model. 
#' @param n.sample an integer denoting how many MCMC samples to use. Default is \code{min(200, n.iter)}, where
#'      n.iter is the number of MCMC iterations (for all chains combined). 
#' @param ptime a numeric vector of times at which predictions 
#'        are to be computed.
#'        If \code{NULL}, \code{ptime = seq(ptime0, floor(max.st), length=25)} is used,
#'        where \code{ptime0} is the last time the subject provided a longitudinal measurement, and
#'        \code{max.st} is \code{floor(max(event times))}.
#' @param plot if \code{TRUE}, a plot of the predictions is produced.
#' @param posterior.mean if \code{TRUE}, posterior means of the predictions are used in the plot,
#'          otherwise posterior medians are used. 
#' @param xlim x limits (a vector of length 2 of the form \code{c(minimum, maximum)}) of the plot (optional).
#' @param ylim y limits (a vector of length 2 of the form \code{c(minimum, maximum)}) of the plot (optional).
#' @param xlab a label (character expression) for the x axis (optional).
#' @param ylab a label (character expression) for the y axis (optional).
#' @param main main title (character expression) of the plot (optional).
#' @param adapt_delta the target average proposal acceptance probability during Stan’s adaptation period.
#' @param max_treedepth a positive integer specifying the maximum treedepth.
#' @param warmup a positive integer specifying the number of burnin iterations.
#' @param seed the seed for random number generation.
#' @details  This function computes subject-specific predictions for the longitudinal
#'       outcome based on the joint model. This is accomplished with a Monte Carlo simulation scheme,
#'       similar to the one described in \code{\link{jm.surv}}. Let \eqn{\theta^{(m)}} be a 
#'       realization taken randomly from the posterior simulations of the joint model parameters
#'       \eqn{(m = 1, 2, ...,} \code{n.sample)}. Also, let
#'       \eqn{b^{(m)}} be a realization of the random-effects vector b from its 
#'       posterior distribution given \eqn{\theta^{(m)}} (see the MCMC scheme described in \code{\link{jm.surv}}). 
#'       Posterior summaries of \eqn{x'_i(s)\alpha + w'_i(s)b_i} (subject-specific prediction)
#'       at time \eqn{s} are derived from the MCMC sample \eqn{\{(\theta^{(m)}, b^{(m)}),m = 1, 2, ...,} \code{n.sample\}}.
#' @return Posterior summaries of the subject-specific predictions.
#'         A plot of the predictions is produced if \code{plot = TRUE}, with
#'         the predictions displayed by the solid curve, the pointwise Bayesian intervals 
#'         displayed by the dashed curves, and \code{ptime[1]} displayed by the dashed vertical line.
#'         Note that with the default definition of \code{ptime}, \code{ptime[1]} is the last time the subject 
#'         provided a longitudinal measurement. 
#' @author Shahedul Khan <khan@math.usask.ca>
#' @references Rizopoulos D, The R package JMbayes for fitting joint models for longitudinal and time-to-event data Using MCMC, 
#'           Journal of Statistical Software, 72(7): 1-45, 2016.
#' @seealso \code{\link{jm.lppredict}}, \code{\link{jm.reg}}, \code{\link{jm.surv}}
#' @examples
#'   # Example: pbc data
#'   lme.fit <- lme(log(bilirubin) ~ drug + ns(futime, 2), data = pbc.long,
#'     random = ~ ns(futime, 2) | id)
#'   surv.fit <- coxph(Surv(st, status2) ~ drug * age, data = pbc.surv, x = TRUE)
#'   jmfit.wph2 <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "weibullph", timevar = "futime", form = "riz")
#'   jm.summary(jmfit.wph2)
#'   newdata2 <- data.frame(pbc.long[pbc.long$id == 2, ])
#'   newdata6 <- data.frame(pbc.long[pbc.long$id == 6, ])
#'   par(mfrow = c(1, 2))
#'   pred2 <- jm.lspredict(jmfit = jmfit.wph2, newdata = newdata2)
#'   pred6 <- jm.lspredict(jmfit = jmfit.wph2, newdata = newdata6)
#'   # For comparison, use the same xlim and ylim in the two plots.
#'   par(mfrow = c(1, 2))
#'   pred2 <-  jm.lspredict(jmfit = jmfit.wph2, newdata = newdata2,
#'       xlim = c(0, 14),  ylim = c(-2.25, 4))
#'   pred6 <- jm.lspredict(jmfit = jmfit.wph2, newdata = newdata6,
#'       xlim = c(0, 14), ylim = c(-2.25, 4))
#' @export

jm.lspredict<-function(jmfit,newdata,n.sample=NULL,
  ptime=NULL,plot=TRUE,posterior.mean=TRUE,
  xlim=NULL,ylim=NULL,xlab=NULL,ylab=NULL,main=NULL,
  adapt_delta=0.9,max_treedepth=15,warmup=200,
  seed=sample.int(.Machine$integer.max, 1)){
start_time <- Sys.time()
if(class(jmfit)!="JMR"){
stop("jmfit must inherit from class JMR")
}
results<-jmfit$all.results
if(is.null(results)){
stop("\nThis function works only for method = 'MCMC'")
}
posterior<-as.matrix(results$stan.fit)
n.iter<-nrow(posterior)
if(is.null(n.sample)){
n.sample<-min(200,n.iter)
}
con<-list(adapt_delta=adapt_delta,max_treedepth=max_treedepth,warmup=warmup,
   n.sample=n.sample,qpoints=15)
n.sample<-con$n.sample
set.seed(seed)
# The new dataset must be a data frame
if(!is.data.frame(newdata)){
newdata<-data.frame(newdata)
}
# Extraact information
surv.model<-results$surv.model
method<-results$method
p1<-results$stan.data$p1
pp1<-results$stan.data$pp1
timevar<-results$timevar
surv.data<-results$surv.data
surv.fit<-results$surv.fit
lme.fit<-results$lme.fit
dat<-lme.fit$data
id.name<-colnames(lme.fit$groups)
if(length(unique(newdata[[id.name]]))>1)
stop(paste("\nThere are",length(unique(newdata[[id.name]])),"subjects in newdata.",
  "\npredictions are allowed for one subject/id at a time.",collapse="")) 
# Design matrix for the survival submodel
newsurvdata<-newdata[1,]
formula.surv <- formula(surv.fit)
survvar<-all.vars(formula.surv)[-(1:2)]
if(isFALSE(all(survvar %in% colnames(newsurvdata)))){
stop("\nVariable names in newdata must be
 the same as in the datasets that were used to fit the survival 
 and the longitudinal models. Not all variables for the 
 survival model are found in newdata.")
}
xsurv<-matrix(model.matrix(formula.surv, newsurvdata)[,-1],nrow=1)
# Design matrix for longitudinal submodel (fixed effects)
formula.fixed<-formula(lme.fit)
longvar<-all.vars(formula.fixed)
if(isFALSE(all(longvar %in% colnames(newdata)))){
stop("\nVariable names in newdata must be
 the same as in the datasets that were used to fit the survival 
 and the longitudinal models. Not all variables for the fixed 
 part of the longitudinal model are found in newdata.")
}
mframe.fdata<-model.frame(terms(formula.fixed),data=dat)
inf.fdata<-attr(mframe.fdata,"terms")
model.fdata<-model.frame(inf.fdata, data = newdata)
y<-as.vector(model.extract(model.fdata, "response"))
xlong<-matrix(model.matrix(formula.fixed,model.fdata),ncol=p1)
# Design matrix for longitudinal submodel (random effects)
formula.rand<-formula(lme.fit$modelStruct$reStruct[[1]])
randvar<-all.vars(formula.rand)
if(isFALSE(all(randvar %in% colnames(newdata)))){
stop("\nVariable names in newdata must be
 the same as in the datasets that were used to fit the survival 
 and the longitudinal models. Not all variables for the random 
 part of the longitudinal model are found in newdata.")
}
mframe.rdata<-model.frame(terms(formula.rand),data = dat)
inf.rdata<-attr(mframe.rdata,"terms")
model.rdata<-model.frame(inf.rdata, data = newdata)
zlong<-matrix(model.matrix(formula.rand,model.rdata),ncol=pp1)
ni<-nrow(zlong)
# time for prediction
if(!is.null(ptime)){
ptime<-sort(ptime)
if(ptime[1]>newdata[nrow(newdata),timevar])
warning(paste("min(ptime) = ", ptime[1], " > last time the individual 
     provided a longitudinal measurement",sep=""))
}
if(is.null(ptime)){
ptime0<-newdata[nrow(newdata),timevar]
max.st<-max(surv.data[,"st"])
ptime<-seq(ptime0,floor(max.st),length=25)
}
# location of the parameters in the MCMC samples
beta.loc<-results$beta.loc
phi.loc<-results$phi.loc
kappa.loc<-results$kappa.loc
rho.loc<-results$rho.loc
lcoef.loc<-results$alpha.loc
lcov.loc<-results$cov.loc
lres.loc<-results$sigma.loc
reffects.loc<-results$b.loc
if(surv.model=="eweibull" || surv.model=="ggamma" || surv.model=="gllph"){
gam.loc<-results$gam.loc
}
sam.dat<-posterior[sample(n.iter,n.sample),]
if(results$surv.model=="weibull" || results$surv.model=="llogistic" || 
 results$surv.model=="lnormal" || results$surv.model=="eweibull" || results$surv.model=="ggamma"){
lme.model<-ifelse(method=="riz",ifelse((results$fixed.model=="simple" && results$rand.model=="simple"),
 "simple","nsimple"),ifelse(results$rand.model=="simple","simple","nsimple"))
} else{
lme.model<-"nsimple"
}
# information for numerical integration
qpoints<-con$qpoints
nodes0<-GK.nodes(qpoints)$gk.n
weight<-GK.nodes(qpoints)$gk.w
nodes1<-sapply(ptime[1],function(st){0.5*(st*nodes0+st)})
newlongdata1<-newdata[1,]
newlongdata2<-newlongdata1[rep(seq_len(nrow(newlongdata1)),qpoints),]
newlongdata2[[timevar]]<-c(nodes1)
newlongdata3<-model.frame(inf.rdata,data=newlongdata2)
newlongdata4<-model.frame(inf.fdata,data=newlongdata2)
nodes<-matrix(model.matrix(formula.rand,newlongdata3),ncol=pp1)
nodesf<-matrix(model.matrix(formula.fixed,newlongdata4),ncol=p1)
nodes.array<-list()
nodesf.array<-list()
for(i in 1:1){
nodes.array[[i]]<-nodes[(1+(i-1)*qpoints):(i*qpoints),]
nodesf.array[[i]]<-nodesf[(1+(i-1)*qpoints):(i*qpoints),]
}
# create data for stan
beta<-sam.dat[,beta.loc]
phi<-sam.dat[,phi.loc]
rho<-sam.dat[,rho.loc]
kappa<-sam.dat[,kappa.loc]
if(surv.model=="eweibull" || surv.model=="ggamma" || surv.model=="gllph"){
gam<-sam.dat[,gam.loc]
} else{
gam<-NULL
}
alpha<-sam.dat[,lcoef.loc]
sigma<-sam.dat[,lres.loc]
pred0<-c(xsurv%*%t(beta))
if(ni>1){
y.mat<-t(replicate(n.sample,y))
} else{
y.mat<-matrix(rep(y,n.sample),ncol=1)
}
if(lme.model=="simple" && method=="riz"){
stan.dat<-results$stan.data
which_alphafixed<-stan.dat$which_alphafixed
which_alphatime<-stan.dat$which_alphatime
alphatime<-alpha[,which_alphatime]
alphafixed<-alpha[,which_alphafixed]
if(is.vector(alphafixed)){
alphafixed<-matrix(alphafixed,ncol=1)
}
xlongfixed<-xlong[,-which_alphatime]
if(is.vector(xlongfixed)){
xlongfixed<-matrix(xlongfixed,ncol=1)
}
xalpha_fixed<-c(matrix(xlongfixed[1,],nrow=1)%*%t(alphafixed))
} else{
alphatime<-alphafixed<-xalpha_fixed<-NULL
}
Sigma<-list()
Sigma.chol<-list()
xalpha<-matrix(0,nrow=n.sample,ncol=nrow(xlong))
xalpha_nodes<-matrix(0,nrow=n.sample,ncol=qpoints)
u<-matrix(0,nrow=sqrt(length(lcov.loc)),ncol=n.sample)
init.u<-rep(0,pp1)
for(i in 1:n.sample){
Sigma[[i]]<-matrix(sam.dat[i,lcov.loc],sqrt(length(lcov.loc)))
Sigma.chol[[i]]<-t(chol(Sigma[[i]]))
xalpha[i,]<-c(xlong%*%alpha[i,])
xalpha_nodes[i,]<-c(nodesf.array[[1]]%*%alpha[i,])
options(warn=-1)
fit.b<-tryCatch({nlminb(start=init.u,objective=pdensb,
 y0=y.mat[i,],xalpha0=xalpha[i,],zlong0=zlong,Sigma0=Sigma[[i]],sigma0=sigma[i],
  st0=ptime[1],pred0=pred0[i],xalpha_fixed0=xalpha_fixed[i],
  alphatime0=alphatime[i],xalpha_nodes=xalpha_nodes[i,],nodesr=nodes.array[[1]],
  nodes0=nodes1[,1],weight=weight,
  phi0=phi[i],rho0=rho[i],kappa0=kappa[i],gam0=gam[i],
  model=surv.model,method=method,lme.model=lme.model)},
  error = function(e){NULL})
options(warn=0)
if(!is.null(fit.b)){
if((fit.b$message=="relative convergence (4)" || 
  fit.b$message=="both X-convergence and relative convergence (5)") &&
  fit.b$convergence==0 && is.finite(fit.b$objective) && 
  fit.b$objective!=0 && !is.na(fit.b$objective)){
u[,i]<-solve(Sigma.chol[[i]])%*%fit.b$par
} else{
u[,i]<-init.u
}
}
if(is.null(fit.b)){
u[,i]<-init.u
}
}
inits<-list(list(b_unscaled=u))
par.stan<-c("b")
#--- setup stan model---------
if(method=="riz"){
if(surv.model=="weibull"){
if(lme.model=="simple"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   xalpha_fixed=xalpha_fixed,alphatime=alphatime,phi=phi,
   rho=rho,kappa=kappa,chol_cov_matrix=Sigma.chol)
#stanc("stan_weibullsurvrizsimple.stan")
#stan_model <- "stan_weibullsurvrizsimple.stan"
stan_model <- stanmodels$weibullsurvrizsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,xalpha1=xalpha_nodes,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_weibullsurvriz.stan")
#stan_model <- "stan_weibullsurvriz.stan"
stan_model <- stanmodels$weibullsurvriz
}
} else if(surv.model=="llogistic"){
if(lme.model=="simple"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   xalpha_fixed=xalpha_fixed,alphatime=alphatime,phi=phi,
   rho=rho,kappa=kappa,chol_cov_matrix=Sigma.chol)
#stanc("stan_llogisticsurvrizsimple.stan")
#stan_model <- "stan_llogisticsurvrizsimple.stan"
stan_model <- stanmodels$llogisticsurvrizsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,xalpha1=xalpha_nodes,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_llogisticsurvriz.stan")
#stan_model <- "stan_llogisticsurvriz.stan"
stan_model <- stanmodels$llogisticsurvriz
} 
} else if(surv.model=="lnormal"){
if(lme.model=="simple"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   xalpha_fixed=xalpha_fixed,alphatime=alphatime,phi=phi,
   rho=rho,kappa=kappa,chol_cov_matrix=Sigma.chol)
#stanc("stan_lnormalsurvrizsimple.stan")
#stan_model <- "stan_lnormalsurvrizsimple.stan"
stan_model <- stanmodels$lnormalsurvrizsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,xalpha1=xalpha_nodes,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_lnormalsurvriz.stan")
#stan_model <- "stan_lnormalsurvriz.stan"
stan_model <- stanmodels$lnormalsurvriz
}
} else if(surv.model=="eweibull"){
if(lme.model=="simple"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   xalpha_fixed=xalpha_fixed,alphatime=alphatime,phi=phi,
   rho=rho,kappa=kappa,gam=gam,chol_cov_matrix=Sigma.chol)
#stanc("stan_eweibullsurvrizsimple.stan")
#stan_model <- "stan_eweibullsurvrizsimple.stan"
stan_model <- stanmodels$eweibullsurvrizsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,xalpha1=xalpha_nodes,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,gam=gam,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_eweibullsurvriz.stan")
#stan_model <- "stan_eweibullsurvriz.stan"
stan_model <- stanmodels$eweibullsurvriz
}
} else if(surv.model=="ggamma"){
if(lme.model=="simple"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   xalpha_fixed=xalpha_fixed,alphatime=alphatime,phi=phi,
   rho=rho,kappa=kappa,gam=gam,chol_cov_matrix=Sigma.chol)
#stanc("stan_ggammasurvrizsimple.stan")
#stan_model <- "stan_ggammasurvrizsimple.stan"
stan_model <- stanmodels$ggammasurvrizsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,xalpha1=xalpha_nodes,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,gam=gam,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_ggammasurvriz.stan")
#stan_model <- "stan_ggammasurvriz.stan"
stan_model <- stanmodels$ggammasurvriz
}
} 
else if(surv.model=="weibullph"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,xalpha1=xalpha_nodes,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,weights=weight,
   nodes0=nodes1[,1],nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_weibullphsurvriz.stan")
#stan_model <- "stan_weibullphsurvriz.stan"
stan_model <- stanmodels$weibullphsurvriz
} 
else if(surv.model=="gllph"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,xalpha1=xalpha_nodes,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,gam=gam,weights=weight,
   nodes0=nodes1[,1],nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_gllphsurvriz.stan")
#stan_model <- "stan_gllphsurvriz.stan"
stan_model <- stanmodels$gllphsurvriz
}
} else{
if(surv.model=="weibull"){
if(lme.model=="sample"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   phi=phi,rho=rho,kappa=kappa,chol_cov_matrix=Sigma.chol)
#stanc("stan_weibullsurvsimple.stan")
#stan_model <- "stan_weibullsurvsimple.stan"
stan_model <- stanmodels$weibullsurvsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_weibullsurv.stan")
#stan_model <- "stan_weibullsurv.stan"
stan_model <- stanmodels$weibullsurv
}
} else if(surv.model=="llogistic"){
if(lme.model=="sample"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   phi=phi,rho=rho,kappa=kappa,chol_cov_matrix=Sigma.chol)
#stanc("stan_llogisticsurvsimple.stan")
#stan_model <- "stan_llogisticsurvsimple.stan"
stan_model <- stanmodels$llogisticsurvsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_llogisticsurv.stan")
#stan_model <- "stan_llogisticsurv.stan"
stan_model <- stanmodels$llogisticsurv
}
} else if(surv.model=="lnormal"){
if(lme.model=="sample"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   phi=phi,rho=rho,kappa=kappa,chol_cov_matrix=Sigma.chol)
#stanc("stan_lnormalsurvsimple.stan")
#stan_model <- "stan_lnormalsurvvsimple.stan"
stan_model <- stanmodels$lnormalsurvsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_lnormalsurv.stan")
#stan_model <- "stan_lnormalsurv.stan"
stan_model <- stanmodels$lnormalsurv
}
} else if(surv.model=="eweibull"){
if(lme.model=="sample"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   phi=phi,rho=rho,kappa=kappa,gam=gam,chol_cov_matrix=Sigma.chol)
#stanc("stan_eweibullsurvsimple.stan")
#stan_model <- "stan_eweibullsurvsimple.stan"
stan_model <- stanmodels$eweibullsurvsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,gam=gam,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_eweibullsurv.stan")
#stan_model <- "stan_eweibullsurv.stan"
stan_model <- stanmodels$eweibullsurv
}
} else if(surv.model=="ggamma"){
if(lme.model=="sample"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,pred0=pred0,
   phi=phi,rho=rho,kappa=kappa,gam=gam,chol_cov_matrix=Sigma.chol)
#stanc("stan_ggammasurvsimple.stan")
#stan_model <- "stan_ggammasurvsimple.stan"
stan_model <- stanmodels$ggammasurvsimple
} else{
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,gam=gam,weights=weight,
   nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_ggammasurv.stan")
#stan_model <- "stan_ggammasurv.stan"
stan_model <- stanmodels$ggammasurv
}
} else if(surv.model=="weibullph"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,weights=weight,
   nodes0=nodes1[,1],nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_weibullphsurv.stan")
#stan_model <- "stan_weibullphsurv.stan"
stan_model <- stanmodels$weibullphsurv
} else if(surv.model=="gllph"){
data.stan<-list(N=n.sample,ni=ni,pp1=pp1,qpoints=qpoints,
   st=ptime[1],y=c(t(y.mat)),sigma=rep(sigma,each=ni),
   xalpha=xalpha,z_rand=zlong,
   pred0=pred0,phi=phi,rho=rho,kappa=kappa,gam=gam,weights=weight,
   nodes0=nodes1[,1],nodes=nodes.array[[1]],chol_cov_matrix=Sigma.chol)
#stanc("stan_gllphsurv.stan")
#stan_model <- "stan_gllphsurv.stan"
stan_model <- stanmodels$gllphsurv
}
}
#---- Run stan -----------------
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
stan_fit <- quiet(suppressWarnings(rstan::sampling(object=stan_model, data = data.stan, 
  init = inits, pars = par.stan,
  warmup = con$warmup, iter = con$warmup+1, 
  chains = 1, thin = 1, seed=seed,
  control=list(adapt_delta=con$adapt_delta,max_treedepth=con$max_treedepth))))
posterior.b<-as.array(stan_fit)
b<-matrix(posterior.b[sample(con$warmup+1-con$warmup,1),1,][-dim(posterior.b)[3]],nrow=n.sample)
# -- predictions ------------
pptime<-unique(c(seq(from=min(newdata[[timevar]]),to=max(newdata[[timevar]]),length.out=ni*10),ptime))
pdat0<-newdata[ni,][rep(seq_len(nrow(newdata[ni,])),length(pptime)),]
pdat0[[timevar]]<-pptime
fmodel.pdat0<-model.frame(inf.fdata, data = pdat0)
fpdat<-matrix(model.matrix(formula.fixed,fmodel.pdat0),ncol=p1)
rownames(fpdat)<-NULL
rmodel.pdat0<-model.frame(inf.rdata, data = pdat0)
rpdat<-matrix(model.matrix(formula.rand,rmodel.pdat0),ncol=pp1)
pred.mcmc<-apply(alpha,1,function(a) fpdat%*%a)+apply(b,1,function(a) rpdat%*%a)
fresults<-cbind(t=pdat0[[timevar]],mean=rowMeans(pred.mcmc,na.rm=TRUE),
       sd=apply(pred.mcmc,1,stats:::sd,na.rm=TRUE),
      "2.5%"=apply(pred.mcmc,1,quantile,prob=0.025,na.rm=TRUE),
      "25%"=apply(pred.mcmc,1,quantile,prob=0.25,na.rm=TRUE),
      "50%"=apply(pred.mcmc,1,quantile,prob=0.5,na.rm=TRUE),
      "75%"=apply(pred.mcmc,1,quantile,prob=0.75,na.rm=TRUE),
      "97.5%"=apply(pred.mcmc,1,quantile,prob=0.975,na.rm=TRUE))
colnames(fresults)[1]<-timevar
rownames(fresults)<-rep("",nrow(fresults))
resultname<-paste("predicted values for the longitudinal response (subject-specific)")
if(plot){
if(posterior.mean){
pred1<-fresults[,2]
lower<-fresults[,4]
upper<-fresults[,8]
} else{
pred1<-fresults[,6]
lower<-fresults[,4]
upper<-fresults[,8]
}
y1<-c(y,rep(NA,length(pptime)-length(y)))
t1<-c(newdata[[timevar]],rep(NA,length(pptime)-length(y)))
if(is.null(ylim)){
ylim0<-c(min(c(y,pred1,lower,upper)),max(c(y,pred1,lower,upper)))
} else{
ylim0<-ylim
}
if(is.null(xlim)){
xlim0<-c(min(fresults[,1]),max(fresults[,1]))
} else{
xlim0<-xlim
}
if(is.null(xlab)){
xlab0<-"Time"
} else{
xlab0<-xlab
}
if(is.null(ylab)){
ylab0<-formula.fixed[[2]]
} else{
ylab0<-ylab
}
par(mar=c(5, 4.25, 4, 2) + 0.1)
if(is.null(main)){
subject.id<-newdata[[id.name]][1]
plot(t1,y1,pch=19,type="p",ylim=ylim0,xlim=xlim0,
  xlab=xlab0,ylab=ylab0,cex.lab=1.15)
title<-list(bquote("Subject" ~ .(subject.id)))
mtext(do.call(expression, title),side=3,cex=1.5,line=0.5)
} else{
plot(t1,y1,pch=19,type="p",ylim=ylim0,xlim=xlim0,
  xlab=xlab0,ylab=ylab0,cex.lab=1.15,main=main)
}
abline(v=ptime[1],lty=2,lwd=3)
lines(fresults[,1],pred1,lwd=2)
lines(fresults[,1],lower,lty=3)
lines(fresults[,1],upper,lty=3)
}
end_time <- Sys.time()
output<-list()
output[[resultname]]<-fresults
output[["time"]]<-noquote(paste(round(difftime(end_time,start_time,units="mins")[[1]],
  digits=2),"mins",collapse=""))
class(output)<-"JMR"
attr(output, "hidden")<-"time"
return(output)
}



