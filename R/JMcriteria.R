#------------------------------------------
# JM log-likelihood given mcmc data (Weibull AFT)
#-------------------------------------------
jmllik.weibull<-function(mcmc.data,st,status,
  y,x,x_fixed=NULL,x_rand,z,z_rand,z_fixed,m,pp1,
  beta.loc,phi.loc,kappa.loc,gam.loc=NULL,rho.loc,
  alpha.loc,sigma.loc,b.loc,which_alphafixed=NULL,which_alphatime=NULL,
  id,id.length,nodes1=NULL,nodesr=NULL,nodesf=NULL,weight=NULL,qpoints=NULL,
  method,lme.model){
beta<-mcmc.data[beta.loc]
phi<-mcmc.data[phi.loc]
kappa<-mcmc.data[kappa.loc]
rho<-mcmc.data[rho.loc]
alpha<-mcmc.data[alpha.loc]
sigma<-mcmc.data[sigma.loc]
b<-matrix(mcmc.data[b.loc],ncol=pp1)
bb<-b[rep(1:m, times = id.length), ]
norm.lden<-dnorm(y,c(x%*%alpha) + rowSums(x_rand*bb),sigma,log=TRUE)
y.lden<-as.vector(tapply(norm.lden,id,sum))
pred0<-c(z%*%beta)
if(lme.model=="simple"){
alphafixed <- alpha[which_alphafixed]
alphatime <- alpha[which_alphatime]
if(method=="riz"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
psi = (exp(- pred0 - phi * (x_fixed %*% alphafixed + b[,1])) *
  (1 - exp(- phi * (alphatime + b[,2]) * st))) / (phi * (alphatime + b[,2]))
} else{
pred1<-pred0+phi*rowSums(z_rand*b)
psi = (exp(- pred0 - phi * b[,1]) * (1 - exp(- phi * b[,2] * st))) / (phi * b[,2])
}
} else{
bbb<-b[rep(1:m, times = rep(qpoints,m)), ]
if(method=="riz"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
psi<-c(t(weight)%*%matrix(exp(-phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*st*exp(-pred0)/2
} else{
pred1<-pred0+phi*rowSums(z_rand*b)
psi<-c(t(weight)%*%matrix(exp(-phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*exp(-pred0)/2
}
}
cum.haz<-(rho*psi)^kappa
logh <- log(kappa) + kappa * log(rho) + (kappa-1) * log(psi) - pred1
llik <- y.lden+status*logh-cum.haz
llik[!is.finite(llik)]<-NA
return(llik)
}
#------------------------------------------
# JM lok-likelihood given mcmc data (Log-logistic AFT)
#-------------------------------------------
jmllik.llogistic<-function(mcmc.data,st,status,
  y,x,x_fixed=NULL,x_rand,z,z_rand,z_fixed,m,pp1,
  beta.loc,phi.loc,kappa.loc,gam.loc=NULL,rho.loc,
  alpha.loc,sigma.loc,b.loc,which_alphafixed=NULL,which_alphatime=NULL,
  id,id.length,nodes1=NULL,nodesr=NULL,nodesf=NULL,weight=NULL,qpoints=NULL,
  method,lme.model){
beta<-mcmc.data[beta.loc]
phi<-mcmc.data[phi.loc]
kappa<-mcmc.data[kappa.loc]
rho<-mcmc.data[rho.loc]
alpha<-mcmc.data[alpha.loc]
sigma<-mcmc.data[sigma.loc]
b<-matrix(mcmc.data[b.loc],ncol=pp1)
bb<-b[rep(1:m, times = id.length), ]
norm.lden<-dnorm(y,c(x%*%alpha) + rowSums(x_rand*bb),sigma,log=TRUE)
y.lden<-as.vector(tapply(norm.lden,id,sum))
pred0<-c(z%*%beta)
if(lme.model=="simple"){
alphafixed <- alpha[which_alphafixed]
alphatime <- alpha[which_alphatime]
if(method=="riz"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
psi = (exp(- pred0 - phi * (x_fixed %*% alphafixed + b[,1])) *
  (1 - exp(- phi * (alphatime + b[,2]) * st))) / (phi * (alphatime + b[,2]))
} else{
pred1<-pred0+phi*rowSums(z_rand*b)
psi = (exp(- pred0 - phi * b[,1]) * (1 - exp(- phi * b[,2] * st))) / (phi * b[,2])
}
} else{
bbb<-b[rep(1:m, times = rep(qpoints,m)), ]
if(method=="riz"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
psi<-c(t(weight)%*%matrix(exp(-phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*st*exp(-pred0)/2
} else{
pred1<-pred0+phi*rowSums(z_rand*b)
psi<-c(t(weight)%*%matrix(exp(-phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*exp(-pred0)/2
}
}
cum.haz<- log(1+(rho*psi)^kappa)
logh <- log(kappa) + kappa * log(rho) + (kappa-1) * log(psi) - log(1+(rho*psi)^kappa) - pred1
llik <- y.lden+status*logh-cum.haz
llik[!is.finite(llik)]<-NA
return(llik)
}
#------------------------------------------
# JM lok-likelihood given mcmc data (Log-normal AFT)
#-------------------------------------------
jmllik.lnormal<-function(mcmc.data,st,status,
  y,x,x_fixed=NULL,x_rand,z,z_rand,z_fixed,m,pp1,
  beta.loc,phi.loc,kappa.loc,gam.loc=NULL,rho.loc,
  alpha.loc,sigma.loc,b.loc,which_alphafixed=NULL,which_alphatime=NULL,
  id,id.length,nodes1=NULL,nodesr=NULL,nodesf=NULL,weight=NULL,qpoints=NULL,
  method,lme.model){
beta<-mcmc.data[beta.loc]
phi<-mcmc.data[phi.loc]
kappa<-mcmc.data[kappa.loc]
rho<-mcmc.data[rho.loc]
alpha<-mcmc.data[alpha.loc]
sigma<-mcmc.data[sigma.loc]
b<-matrix(mcmc.data[b.loc],ncol=pp1)
bb<-b[rep(1:m, times = id.length), ]
norm.lden<-dnorm(y,c(x%*%alpha) + rowSums(x_rand*bb),sigma,log=TRUE)
y.lden<-as.vector(tapply(norm.lden,id,sum))
pred0<-c(z%*%beta)
if(lme.model=="simple"){
alphafixed <- alpha[which_alphafixed]
alphatime <- alpha[which_alphatime]
if(method=="riz"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
psi = (exp(- pred0 - phi * (x_fixed %*% alphafixed + b[,1])) *
  (1 - exp(- phi * (alphatime + b[,2]) * st))) / (phi * (alphatime + b[,2]))
} else{
pred1<-pred0+phi*rowSums(z_rand*b)
psi = (exp(- pred0 - phi * b[,1]) * (1 - exp(- phi * b[,2] * st))) / (phi * b[,2])
}
} else{
bbb<-b[rep(1:m, times = rep(qpoints,m)), ]
if(method=="riz"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
psi<-c(t(weight)%*%matrix(exp(-phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*st*exp(-pred0)/2
} else{
pred1<-pred0+phi*rowSums(z_rand*b)
psi<-c(t(weight)%*%matrix(exp(-phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*exp(-pred0)/2
}
}
cum.haz<- - plnorm(psi, meanlog = -log(rho), sdlog = 1/kappa, lower.tail = FALSE, log.p = TRUE)
logf <- dlnorm(psi, meanlog = -log(rho), sdlog = 1/kappa, log = TRUE)-pred1
llik <- y.lden+status*logf-(1-status)*cum.haz
llik[!is.finite(llik)]<-NA
return(llik)
}
#------------------------------------------
# JM lok-likelihood given mcmc data (Weibull PH)
#-------------------------------------------
jmllik.weibullph<-function(mcmc.data,st,status,
  y,x,x_fixed=NULL,x_rand,z,z_rand,z_fixed,m,pp1,
  beta.loc,phi.loc,kappa.loc,gam.loc=NULL,rho.loc,
  alpha.loc,sigma.loc,b.loc,which_alphafixed=NULL,which_alphatime=NULL,
  id,id.length,nodes1=NULL,nodesr=NULL,nodesf=NULL,weight=NULL,qpoints=NULL,
  method,lme.model){
beta<-mcmc.data[beta.loc]
phi<-mcmc.data[phi.loc]
kappa<-mcmc.data[kappa.loc]
rho<-mcmc.data[rho.loc]
alpha<-mcmc.data[alpha.loc]
sigma<-mcmc.data[sigma.loc]
b<-matrix(mcmc.data[b.loc],ncol=pp1)
bb<-b[rep(1:m, times = id.length), ]
norm.lden<-dnorm(y,c(x%*%alpha) + rowSums(x_rand*bb),sigma,log=TRUE)
y.lden<-as.vector(tapply(norm.lden,id,sum))
pred0<-c(z%*%beta)
bbb<-b[rep(1:m, times = rep(qpoints,m)), ]
if(method=="riz"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
cum.haz<-c(t(weight)%*%matrix(exp((kappa-1)*log(nodes1)+
  phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*
  st*kappa*rho^kappa*exp(pred0)/2
logh <- log(kappa) + kappa * log(rho) + (kappa-1) * log(st) + pred1
} else{
pred1<-pred0+phi*rowSums(z_rand*b)
cum.haz<-c(t(weight)%*%matrix(exp((kappa-1)*log(nodes1)+
  phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*kappa*rho^kappa*exp(pred0)/2
logh <- log(kappa) + kappa * log(rho) + (kappa-1) * log(st) + pred1
}
llik <- y.lden+status*logh-cum.haz
llik[!is.finite(llik)]<-NA
return(llik)
}
#------------------------------------------
# JM lok-likelihood given mcmc data (eweibull AFT)
#-------------------------------------------
jmllik.eweibull<-function(mcmc.data,st,status,
  y,x,x_fixed=NULL,x_rand,z,z_rand,z_fixed,m,pp1,
  beta.loc,phi.loc,kappa.loc,gam.loc=NULL,rho.loc,
  alpha.loc,sigma.loc,b.loc,which_alphafixed=NULL,which_alphatime=NULL,
  id,id.length,nodes1=NULL,nodesr=NULL,nodesf=NULL,weight=NULL,qpoints=NULL,
  method,lme.model){
beta<-mcmc.data[beta.loc]
phi<-mcmc.data[phi.loc]
kappa<-mcmc.data[kappa.loc]
gam<-mcmc.data[gam.loc]
rho<-mcmc.data[rho.loc]
alpha<-mcmc.data[alpha.loc]
sigma<-mcmc.data[sigma.loc]
b<-matrix(mcmc.data[b.loc],ncol=pp1)
bb<-b[rep(1:m, times = id.length), ]
norm.lden<-dnorm(y,c(x%*%alpha) + rowSums(x_rand*bb),sigma,log=TRUE)
y.lden<-as.vector(tapply(norm.lden,id,sum))
pred0<-c(z%*%beta)
if(lme.model=="simple"){
alphafixed <- alpha[which_alphafixed]
alphatime <- alpha[which_alphatime]
if(method=="riz"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
psi = (exp(- pred0 - phi * (x_fixed %*% alphafixed + b[,1])) *
  (1 - exp(- phi * (alphatime + b[,2]) * st))) / (phi * (alphatime + b[,2]))
} else{
pred1<-pred0+phi*rowSums(z_rand*b)
psi = (exp(- pred0 - phi * b[,1]) * (1 - exp(- phi * b[,2] * st))) / (phi * b[,2])
}
} else{
bbb<-b[rep(1:m, times = rep(qpoints,m)), ]
if(method=="riz"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
psi<-c(t(weight)%*%matrix(exp(-phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*st*exp(-pred0)/2
} else{
pred1<-pred0+phi*rowSums(z_rand*b)
psi<-c(t(weight)%*%matrix(exp(-phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*exp(-pred0)/2
}
}
logf<-log(kappa)+log(gam)+kappa*log(rho)+(kappa-1)*log(psi)+
   (gam-1)*log1mexp(-(rho*psi)^kappa)-(rho*psi)^kappa-pred1
logs<-log1p(-exp(gam*log1mexp(-(rho*psi)^kappa)))
if(any(!is.finite(logs))){
wh<-which(!is.finite(logs))
logs[wh]<-log(gam)-(rho*psi[wh])^kappa
}
cum.haz<- -logs
llik <- y.lden+status*logf-(1-status)*cum.haz
llik[!is.finite(llik)]<-NA
return(llik)
}
#------------------------------------------
# JM lok-likelihood given mcmc data (ggamma AFT)
#-------------------------------------------
jmllik.ggamma<-function(mcmc.data,st,status,
  y,x,x_fixed=NULL,x_rand,z,z_rand,z_fixed,m,pp1,
  beta.loc,phi.loc,kappa.loc,gam.loc=NULL,rho.loc,
  alpha.loc,sigma.loc,b.loc,which_alphafixed=NULL,which_alphatime=NULL,
  id,id.length,nodes1=NULL,nodesr=NULL,nodesf=NULL,weight=NULL,qpoints=NULL,
  method,lme.model){
beta<-mcmc.data[beta.loc]
phi<-mcmc.data[phi.loc]
kappa<-mcmc.data[kappa.loc]
gam<-mcmc.data[gam.loc]
rho<-mcmc.data[rho.loc]
alpha<-mcmc.data[alpha.loc]
sigma<-mcmc.data[sigma.loc]
b<-matrix(mcmc.data[b.loc],ncol=pp1)
bb<-b[rep(1:m, times = id.length), ]
norm.lden<-dnorm(y,c(x%*%alpha) + rowSums(x_rand*bb),sigma,log=TRUE)
y.lden<-as.vector(tapply(norm.lden,id,sum))
pred0<-c(z%*%beta)
if(lme.model=="simple"){
alphafixed <- alpha[which_alphafixed]
alphatime <- alpha[which_alphatime]
if(method=="riz"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
psi = (exp(- pred0 - phi * (x_fixed %*% alphafixed + b[,1])) *
  (1 - exp(- phi * (alphatime + b[,2]) * st))) / (phi * (alphatime + b[,2]))
} else{
pred1<-pred0+phi*rowSums(z_rand*b)
psi = (exp(- pred0 - phi * b[,1]) * (1 - exp(- phi * b[,2] * st))) / (phi * b[,2])
}
} else{
bbb<-b[rep(1:m, times = rep(qpoints,m)), ]
if(method=="riz"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
psi<-c(t(weight)%*%matrix(exp(-phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*st*exp(-pred0)/2
} else{
pred1<-pred0+phi*rowSums(z_rand*b)
psi<-c(t(weight)%*%matrix(exp(-phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*exp(-pred0)/2
}
}
logf<-log(kappa)+kappa*gam*log(rho)+(kappa*gam-1)*log(psi)-(rho*psi)^kappa-lgamma(gam)-pred1
cum.haz<--pgamma((rho*psi)^kappa,shape=gam,scale=1,lower.tail=FALSE,log.p=TRUE)
llik <- y.lden+status*logf-(1-status)*cum.haz
llik[!is.finite(llik)]<-NA
return(llik)
}
#------------------------------------------
# JM lok-likelihood given mcmc data (GLL PH)
#-------------------------------------------
jmllik.gllogistic<-function(mcmc.data,st,status,
  y,x,x_fixed=NULL,x_rand,z,z_rand,z_fixed,m,pp1,
  beta.loc,phi.loc,kappa.loc,gam.loc=NULL,rho.loc,
  alpha.loc,sigma.loc,b.loc,which_alphafixed=NULL,which_alphatime=NULL,
  id,id.length,nodes1=NULL,nodesr=NULL,nodesf=NULL,weight=NULL,qpoints=NULL,
  method,lme.model){
beta<-mcmc.data[beta.loc]
phi<-mcmc.data[phi.loc]
kappa<-mcmc.data[kappa.loc]
gam<-mcmc.data[gam.loc]
rho<-mcmc.data[rho.loc]
alpha<-mcmc.data[alpha.loc]
sigma<-mcmc.data[sigma.loc]
b<-matrix(mcmc.data[b.loc],ncol=pp1)
bb<-b[rep(1:m, times = id.length), ]
norm.lden<-dnorm(y,c(x%*%alpha) + rowSums(x_rand*bb),sigma,log=TRUE)
y.lden<-as.vector(tapply(norm.lden,id,sum))
pred0<-c(z%*%beta)
bbb<-b[rep(1:m, times = rep(qpoints,m)), ]
if(method=="riz"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
cum.haz<-c(t(weight)%*%matrix(exp((kappa-1)*log(nodes1)-log(1+(gam*nodes1)^kappa)+
   phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*
   st*kappa*rho^kappa*exp(pred0)/2
logh <- log(kappa) + kappa * log(rho) + (kappa-1) * log(st) - log(1+(gam*st)^kappa) + pred1
} else{
pred1<-pred0+phi*rowSums(z_rand*b)
cum.haz<-c(t(weight)%*%matrix(exp((kappa-1)*log(nodes1)-log(1+(gam*nodes1)^kappa)+
   phi*rowSums(nodesr*bbb)),nrow=qpoints))*
   st*kappa*rho^kappa*exp(pred0)/2
logh <- log(kappa) + kappa * log(rho) + (kappa-1) * log(st) - log(1+(gam*st)^kappa) + pred1
}
llik <- y.lden+status*logh-cum.haz
llik[!is.finite(llik)]<-NA
return(llik)
}
############################################
# Information criteria
#############################################
########################################################################
#' Information criteria (DIC and WAIC)
#' @description Computes DIC and WAIC for a Bayesian joint model fit. It uses a \code{\link{jm.reg}} fit as its argument.
#' @keywords Joint modeling, MCMC, DIC, WAIC
#' @param jmfit a \code{\link{jm.reg}} fit with \code{method = "MCMC"}.
#' @param posterior.mean if \code{TRUE} (default), the posterior mean is used to compute
#'        DIC, otherwise posterior median is used. Note that this argument
#'        is used only for DIC calculation.
#' @details This function uses the posterior draws from class "stanfit". See Gelman et al. (2014) 
#'       for a description of DIC and WAIC from a Bayesian perspective.
#' @return \code{DIC:} the values of log predictive density (lpd) or log-likelihood given a point estimate of the 
#'         fitted model, effective number of parameters (p), and DIC. If \code{posterior.mean = TRUE},
#'         the posterior mean is used as a point estimate, otherwise posterior median is used.
#' @return \code{WAIC:} the values of log pointwise predictive density (lppd) evaluated using the posterior 
#'         simulations, effective number of parameters (p), and WAIC.
#' @references Gelman A, Hwang J, and Vehtari A, Understanding predictive information criteria 
#'     for Bayesian models, Statistics and Computing, 24: 997-1016, 2014.  
#' @author Shahedul Khan <khan@math.usask.ca>
#' @seealso \code{\link{jm.plots}}, \code{\link{jm.reg}}, \code{\link{jm.resid}}
#' @examples 
#'   # Example: pbc data
#'   lme.fit <- lme(log(bilirubin) ~ drug + ns(futime, 2),
#'     data = pbc.long, random = ~ ns(futime, 2) | id)
#'   surv.fit <- coxph(Surv(st, status2) ~ drug * age, data = pbc.surv, x = TRUE)
#'   jmfit.ll2 <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "llogistic", timevar = "futime", form = "riz")
#'   jm.icriteria(jmfit.ll2)
#'   jmfit.wph2 <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "weibullph", timevar = "futime", form = "riz")
#'   jm.icriteria(jmfit.wph2)
#' @export


jm.icriteria<-function(jmfit,posterior.mean=TRUE){
if(class(jmfit)!="JMR"){
stop("jmfit must inherit from class JMR")
}
results<-jmfit$all.results
if(is.null(results)){
stop("\nThis function works only for method = 'MCMC'")
}
post.mean<-results$post.mean
post.med<-results$post.med
dat<-results$stan.data
st<-dat$st
y<-dat$y
status<-dat$status
x<-dat$x
x_rand<-dat$x_rand
z<-dat$z
z_rand<-dat$z_rand
z_fixed<-dat$z_fixed
m<-dat$m
pp1<-dat$pp1
id<-dat$id
beta.loc<-results$beta.loc
phi.loc<-results$phi.loc
kappa.loc<-results$kappa.loc
rho.loc<-results$rho.loc
alpha.loc<-results$alpha.loc
sigma.loc<-results$sigma.loc
b.loc<-results$b.loc
method<-results$method
if(results$surv.model=="weibull" || results$surv.model=="llogistic" || 
 results$surv.model=="lnormal" || results$surv.model=="eweibull" || results$surv.model=="ggamma"){
lme.model<-ifelse(method=="riz",ifelse((results$fixed.model=="simple" && results$rand.model=="simple"),
 "simple","nsimple"),ifelse(results$rand.model=="simple","simple","nsimple"))
} else{
lme.model<-"nsimple"
}
if(lme.model=="simple"){
nodesf<-nodesr<-nodes1<-weight<-qpoints<-NULL
x_fixed<-dat$x_fixed
which_alphafixed<-dat$which_alphafixed
which_alphatime<-dat$which_alphatime
} else{
nodesf<-dat$nodesf
nodesr<-dat$nodesr
nodes1<-c(dat$nodes1)
weight<-dat$weights
qpoints<-dat$qpoints
x_fixed<-which_alphafixed<-which_alphatime<-NULL
}
mcmc.sample<-as.matrix(results$stan.fit)
if(results$surv.model=="weibull"){
gam.loc<-NULL
jmllik<-jmllik.weibull
p.mean<-c(post.mean$beta,post.mean$phi,
  post.mean$rho,post.mean$kappa,
  post.mean$alpha,c(post.mean$cov_rand),post.mean$sigma,c(post.mean$b))
p.med<-c(post.med$beta,post.med$phi,
  post.med$rho,post.med$kappa,
  post.med$alpha,c(post.med$cov_rand),post.med$sigma,c(post.med$b))
} else if(results$surv.model=="llogistic"){
gam.loc<-NULL
jmllik<-jmllik.llogistic
p.mean<-c(post.mean$beta,post.mean$phi,
  post.mean$rho,post.mean$kappa,
  post.mean$alpha,c(post.mean$cov_rand),post.mean$sigma,c(post.mean$b))
p.med<-c(post.med$beta,post.med$phi,
  post.med$rho,post.med$kappa,
  post.med$alpha,c(post.med$cov_rand),post.med$sigma,c(post.med$b))
} else if(results$surv.model=="lnormal"){
gam.loc<-NULL
jmllik<-jmllik.lnormal
p.mean<-c(post.mean$beta,post.mean$phi,
  post.mean$rho,post.mean$kappa,
  post.mean$alpha,c(post.mean$cov_rand),post.mean$sigma,c(post.mean$b))
p.med<-c(post.med$beta,post.med$phi,
  post.med$rho,post.med$kappa,
  post.med$alpha,c(post.med$cov_rand),post.med$sigma,c(post.med$b))
} else if(results$surv.model=="weibullph"){
gam.loc<-NULL
jmllik<-jmllik.weibullph
p.mean<-c(post.mean$beta,post.mean$phi,
  post.mean$rho,post.mean$kappa,
  post.mean$alpha,c(post.mean$cov_rand),post.mean$sigma,c(post.mean$b))
p.med<-c(post.med$beta,post.med$phi,
  post.med$rho,post.med$kappa,
  post.med$alpha,c(post.med$cov_rand),post.med$sigma,c(post.med$b))
} else if(results$surv.model=="eweibull"){
gam.loc<-results$gam.loc
jmllik<-jmllik.eweibull
p.mean<-c(post.mean$beta,post.mean$phi,
  post.mean$rho,post.mean$kappa,post.mean$gam,
  post.mean$alpha,c(post.mean$cov_rand),post.mean$sigma,c(post.mean$b))
p.med<-c(post.med$beta,post.med$phi,
  post.med$rho,post.med$kappa,post.med$gam,
  post.med$alpha,c(post.med$cov_rand),post.med$sigma,c(post.med$b))
} else if(results$surv.model=="ggamma"){
gam.loc<-results$gam.loc
jmllik<-jmllik.ggamma
p.mean<-c(post.mean$beta,post.mean$phi,
  post.mean$rho,post.mean$kappa,post.mean$gam,
  post.mean$alpha,c(post.mean$cov_rand),post.mean$sigma,c(post.mean$b))
p.med<-c(post.med$beta,post.med$phi,
  post.med$rho,post.med$kappa,post.med$gam,
  post.med$alpha,c(post.med$cov_rand),post.med$sigma,c(post.med$b))
} else{
gam.loc<-results$gam.loc
jmllik<-jmllik.gllogistic
p.mean<-c(post.mean$beta,post.mean$phi,
  post.mean$rho,post.mean$kappa,post.mean$gam,
  post.mean$alpha,c(post.mean$cov_rand),post.mean$sigma,c(post.mean$b))
p.med<-c(post.med$beta,post.med$phi,
  post.med$rho,post.med$kappa,post.med$gam,
  post.med$alpha,c(post.med$cov_rand),post.med$sigma,c(post.med$b))
}
llik<-apply(mcmc.sample,1,jmllik,st=st,status=status,y=y,
  x=x,x_fixed=x_fixed,x_rand=x_rand,z=z,z_rand=z_rand,z_fixed=z_fixed,
  m=m,pp1=pp1,beta.loc=beta.loc,phi.loc=phi.loc,
  kappa.loc=kappa.loc,rho.loc=rho.loc,gam.loc=gam.loc,
  alpha.loc=alpha.loc,sigma.loc=sigma.loc,b.loc=b.loc,
  which_alphafixed=which_alphafixed,which_alphatime=which_alphatime,
  id=id,id.length=as.vector(table(id)),
  nodes1=nodes1,nodesr=nodesr,nodesf=nodesf,weight=weight,
  qpoints=qpoints,method=method,lme.model=lme.model)
if(posterior.mean){
llik0<-jmllik(p.mean,st=st,status=status,y=y,
  x=x,x_fixed=x_fixed,x_rand=x_rand,z=z,z_rand=z_rand,z_fixed=z_fixed,
  m=m,pp1=pp1,beta.loc=beta.loc,phi.loc=phi.loc,
  kappa.loc=kappa.loc,rho.loc=rho.loc,gam.loc=gam.loc,
  alpha.loc=alpha.loc,sigma.loc=sigma.loc,b.loc=b.loc,
  which_alphafixed=which_alphafixed,which_alphatime=which_alphatime,
  id=id,id.length=as.vector(table(id)),
  nodes1=nodes1,nodesr=nodesr,nodesf=nodesf,weight=weight,
  qpoints=qpoints,method=method,lme.model=lme.model)
} else{
llik0<-jmllik(p.med,st=st,status=status,y=y,
  x=x,x_fixed=x_fixed,x_rand=x_rand,z=z,z_rand=z_rand,z_fixed=z_fixed,
  m=m,pp1=pp1,beta.loc=beta.loc,phi.loc=phi.loc,
  kappa.loc=kappa.loc,rho.loc=rho.loc,gam.loc=gam.loc,
  alpha.loc=alpha.loc,sigma.loc=sigma.loc,b.loc=b.loc,
  which_alphafixed=which_alphafixed,which_alphatime=which_alphatime,
  id=id,id.length=as.vector(table(id)),
  nodes1=nodes1,nodesr=nodesr,nodesf=nodesf,weight=weight,
  qpoints=qpoints,method=method,lme.model=lme.model)
}
#--- Calculation of WAIC and DIC --------
pwaic<-sum(apply(llik,1,var,na.rm=TRUE),na.rm=TRUE)
lppd<-sum(log(rowSums(exp(llik),na.rm=TRUE)/dim(llik)[2]),na.rm=TRUE)
waic<--2*(lppd-pwaic)
waic.results<-cbind(lppd=lppd,p=pwaic,WAIC=waic)
rownames(waic.results)<-""
llik.bar<-mean(colSums(llik,na.rm=TRUE))
llik.hat<-sum(llik0,na.rm=TRUE)
p.dic<-2*(llik.hat-llik.bar)
dic<--2*(llik.hat-p.dic)
dic.results<-cbind(llik.hat,p.dic,dic)
if(posterior.mean){
colnames(dic.results)<-c("lpd(posterior mean)","p","DIC")
} else{
colnames(dic.results)<-c("lpd(posterior median)","p","DIC")
}
rownames(dic.results)<-""
cat(paste("  lpd = log predictive density or log-likelihood"),
    paste("        given a point estimate of the fitted model"),
    paste("  lppd = log pointwise predictive density"),
    paste("         evaluated using the posterior simulations"),
    paste("  p = effective number of parameters"),
    paste("  DIC = - 2 * (lpd - p)"), 
    paste("  WAIC = - 2 * (lppd - p)"),paste(" "),sep="\n")
output<-list(DIC=dic.results,WAIC=waic.results)
return(output)
}

