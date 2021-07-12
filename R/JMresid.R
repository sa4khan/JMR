#------------------------------------------
# JM cumulative hazard given mcmc data (Weibull AFT)
#-------------------------------------------
jmcumhaz.weibull<-function(mcmc.data,st,m,pp1,z,x_fixed=NULL,
  beta.loc,phi.loc,kappa.loc,gam.loc=NULL,rho.loc,
  alpha.loc,b.loc,which_alphafixed=NULL,which_alphatime=NULL,
  nodes1=NULL,nodesr=NULL,nodesf=NULL,weight=NULL,qpoints=NULL,
  method,lme.model){
beta<-mcmc.data[beta.loc]
phi<-mcmc.data[phi.loc]
kappa<-mcmc.data[kappa.loc]
rho<-mcmc.data[rho.loc]
alpha<-mcmc.data[alpha.loc]
b<-matrix(mcmc.data[b.loc],ncol=pp1)
pred0<-c(z%*%beta)
if(lme.model=="simple"){
alphafixed <- alpha[which_alphafixed]
alphatime <- alpha[which_alphatime]
if(method=="riz"){
psi = (exp(- pred0 - phi * (x_fixed %*% alphafixed + b[,1])) *
  (1 - exp(- phi * (alphatime + b[,2]) * st))) / (phi * (alphatime + b[,2]))
} else{
psi = (exp(- pred0 - phi * b[,1]) * (1 - exp(- phi * b[,2] * st))) / (phi * b[,2])
}
} else{
bbb<-b[rep(1:m, times = rep(qpoints,m)), ]
if(method=="riz"){
psi<-c(t(weight)%*%matrix(exp(-phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*st*exp(-pred0)/2
} else{
psi<-c(t(weight)%*%matrix(exp(-phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*exp(-pred0)/2
}
}
cum.haz<-(rho*psi)^kappa
return(cum.haz)
}
#------------------------------------------
# JM cumulative hazard given mcmc data (Log-logistic AFT)
#-------------------------------------------
jmcumhaz.llogistic<-function(mcmc.data,st,m,pp1,z,x_fixed=NULL,
  beta.loc,phi.loc,kappa.loc,gam.loc=NULL,rho.loc,
  alpha.loc,b.loc,which_alphafixed=NULL,which_alphatime=NULL,
  nodes1=NULL,nodesr=NULL,nodesf=NULL,weight=NULL,qpoints=NULL,
  method,lme.model){
beta<-mcmc.data[beta.loc]
phi<-mcmc.data[phi.loc]
kappa<-mcmc.data[kappa.loc]
rho<-mcmc.data[rho.loc]
alpha<-mcmc.data[alpha.loc]
b<-matrix(mcmc.data[b.loc],ncol=pp1)
pred0<-c(z%*%beta)
if(lme.model=="simple"){
alphafixed <- alpha[which_alphafixed]
alphatime <- alpha[which_alphatime]
if(method=="riz"){
psi = (exp(- pred0 - phi * (x_fixed %*% alphafixed + b[,1])) *
  (1 - exp(- phi * (alphatime + b[,2]) * st))) / (phi * (alphatime + b[,2]))
} else{
psi = (exp(- pred0 - phi * b[,1]) * (1 - exp(- phi * b[,2] * st))) / (phi * b[,2])
}
} else{
bbb<-b[rep(1:m, times = rep(qpoints,m)), ]
if(method=="riz"){
psi<-c(t(weight)%*%matrix(exp(-phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*st*exp(-pred0)/2
} else{
psi<-c(t(weight)%*%matrix(exp(-phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*exp(-pred0)/2
}
}
cum.haz<-log(1+(rho*psi)^kappa)
return(cum.haz)
}
#------------------------------------------
# JM cumulative hazard given mcmc data (Log-normal AFT)
#-------------------------------------------
jmcumhaz.lnormal<-function(mcmc.data,st,m,pp1,z,x_fixed=NULL,
  beta.loc,phi.loc,kappa.loc,gam.loc=NULL,rho.loc,
  alpha.loc,b.loc,which_alphafixed=NULL,which_alphatime=NULL,
  nodes1=NULL,nodesr=NULL,nodesf=NULL,weight=NULL,qpoints=NULL,
  method,lme.model){
beta<-mcmc.data[beta.loc]
phi<-mcmc.data[phi.loc]
kappa<-mcmc.data[kappa.loc]
rho<-mcmc.data[rho.loc]
alpha<-mcmc.data[alpha.loc]
b<-matrix(mcmc.data[b.loc],ncol=pp1)
pred0<-c(z%*%beta)
if(lme.model=="simple"){
alphafixed <- alpha[which_alphafixed]
alphatime <- alpha[which_alphatime]
if(method=="riz"){
psi = (exp(- pred0 - phi * (x_fixed %*% alphafixed + b[,1])) *
  (1 - exp(- phi * (alphatime + b[,2]) * st))) / (phi * (alphatime + b[,2]))
} else{
psi = (exp(- pred0 - phi * b[,1]) * (1 - exp(- phi * b[,2] * st))) / (phi * b[,2])
}
} else{
bbb<-b[rep(1:m, times = rep(qpoints,m)), ]
if(method=="riz"){
psi<-c(t(weight)%*%matrix(exp(-phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*st*exp(-pred0)/2
} else{
psi<-c(t(weight)%*%matrix(exp(-phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*exp(-pred0)/2
}
}
cum.haz<- - plnorm(psi, meanlog = -log(rho), sdlog = 1/kappa, lower.tail = FALSE, log.p = TRUE)
return(cum.haz)
}
#------------------------------------------
# JM cumulative hazard given mcmc data (Weibull PH)
#-------------------------------------------
jmcumhaz.weibullph<-function(mcmc.data,st,m,pp1,z,x_fixed=NULL,
  beta.loc,phi.loc,kappa.loc,gam.loc=NULL,rho.loc,
  alpha.loc,b.loc,which_alphafixed=NULL,which_alphatime=NULL,
  nodes1=NULL,nodesr=NULL,nodesf=NULL,weight=NULL,qpoints=NULL,
  method,lme.model){
beta<-mcmc.data[beta.loc]
phi<-mcmc.data[phi.loc]
kappa<-mcmc.data[kappa.loc]
rho<-mcmc.data[rho.loc]
alpha<-mcmc.data[alpha.loc]
b<-matrix(mcmc.data[b.loc],ncol=pp1)
pred0<-c(z%*%beta)
bbb<-b[rep(1:m, times = rep(qpoints,m)), ]
if(method=="riz"){
cum.haz<-c(t(weight)%*%matrix(exp((kappa-1)*log(nodes1)+
  phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*
  st*kappa*rho^kappa*exp(pred0)/2
} else{
cum.haz<-c(t(weight)%*%matrix(exp((kappa-1)*log(nodes1)+
  phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*kappa*rho^kappa*exp(pred0)/2
}
return(cum.haz)
}
#------------------------------------------
# JM cumulative hazard given mcmc data (eweibull AFT)
#-------------------------------------------
jmcumhaz.eweibull<-function(mcmc.data,st,m,pp1,z,x_fixed=NULL,
  beta.loc,phi.loc,kappa.loc,gam.loc=NULL,rho.loc,
  alpha.loc,b.loc,which_alphafixed=NULL,which_alphatime=NULL,
  nodes1=NULL,nodesr=NULL,nodesf=NULL,weight=NULL,qpoints=NULL,
  method,lme.model){
beta<-mcmc.data[beta.loc]
phi<-mcmc.data[phi.loc]
kappa<-mcmc.data[kappa.loc]
rho<-mcmc.data[rho.loc]
gam<-mcmc.data[gam.loc]
alpha<-mcmc.data[alpha.loc]
b<-matrix(mcmc.data[b.loc],ncol=pp1)
pred0<-c(z%*%beta)
if(lme.model=="simple"){
alphafixed <- alpha[which_alphafixed]
alphatime <- alpha[which_alphatime]
if(method=="riz"){
psi = (exp(- pred0 - phi * (x_fixed %*% alphafixed + b[,1])) *
  (1 - exp(- phi * (alphatime + b[,2]) * st))) / (phi * (alphatime + b[,2]))
} else{
psi = (exp(- pred0 - phi * b[,1]) * (1 - exp(- phi * b[,2] * st))) / (phi * b[,2])
}
} else{
bbb<-b[rep(1:m, times = rep(qpoints,m)), ]
if(method=="riz"){
psi<-c(t(weight)%*%matrix(exp(-phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*st*exp(-pred0)/2
} else{
psi<-c(t(weight)%*%matrix(exp(-phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*exp(-pred0)/2
}
}
logs<-log1p(-exp(gam*log1mexp(-(rho*psi)^kappa)))
if(any(!is.finite(logs))){
wh<-which(!is.finite(logs))
logs[wh]<-log(gam)-(rho*psi[wh])^kappa
}
cum.haz<- -logs
return(cum.haz)
}
#------------------------------------------
# JM cumulative hazard given mcmc data (ggamma AFT)
#-------------------------------------------
jmcumhaz.ggamma<-function(mcmc.data,st,m,pp1,z,x_fixed=NULL,
  beta.loc,phi.loc,kappa.loc,gam.loc=NULL,rho.loc,
  alpha.loc,b.loc,which_alphafixed=NULL,which_alphatime=NULL,
  nodes1=NULL,nodesr=NULL,nodesf=NULL,weight=NULL,qpoints=NULL,
  method,lme.model){
beta<-mcmc.data[beta.loc]
phi<-mcmc.data[phi.loc]
kappa<-mcmc.data[kappa.loc]
rho<-mcmc.data[rho.loc]
gam<-mcmc.data[gam.loc]
alpha<-mcmc.data[alpha.loc]
b<-matrix(mcmc.data[b.loc],ncol=pp1)
pred0<-c(z%*%beta)
if(lme.model=="simple"){
alphafixed <- alpha[which_alphafixed]
alphatime <- alpha[which_alphatime]
if(method=="riz"){
psi = (exp(- pred0 - phi * (x_fixed %*% alphafixed + b[,1])) *
  (1 - exp(- phi * (alphatime + b[,2]) * st))) / (phi * (alphatime + b[,2]))
} else{
psi = (exp(- pred0 - phi * b[,1]) * (1 - exp(- phi * b[,2] * st))) / (phi * b[,2])
}
} else{
bbb<-b[rep(1:m, times = rep(qpoints,m)), ]
if(method=="riz"){
psi<-c(t(weight)%*%matrix(exp(-phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*st*exp(-pred0)/2
} else{
psi<-c(t(weight)%*%matrix(exp(-phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*exp(-pred0)/2
}
}
cum.haz<--pgamma((rho*psi)^kappa,shape=gam,scale=1,lower.tail=FALSE,log.p=TRUE)
return(cum.haz)
}
#------------------------------------------
# JM cumulative hazard given mcmc data (GLL PH)
#-------------------------------------------
jmcumhaz.gllogistic<-function(mcmc.data,st,m,pp1,z,x_fixed=NULL,
  beta.loc,phi.loc,kappa.loc,gam.loc=NULL,rho.loc,
  alpha.loc,b.loc,which_alphafixed=NULL,which_alphatime=NULL,
  nodes1=NULL,nodesr=NULL,nodesf=NULL,weight=NULL,qpoints=NULL,
  method,lme.model){
beta<-mcmc.data[beta.loc]
phi<-mcmc.data[phi.loc]
kappa<-mcmc.data[kappa.loc]
rho<-mcmc.data[rho.loc]
gam<-mcmc.data[gam.loc]
alpha<-mcmc.data[alpha.loc]
b<-matrix(mcmc.data[b.loc],ncol=pp1)
pred0<-c(z%*%beta)
bbb<-b[rep(1:m, times = rep(qpoints,m)), ]
if(method=="riz"){
cum.haz<-c(t(weight)%*%matrix(exp((kappa-1)*log(nodes1)-log(1+(gam*nodes1)^kappa)+
   phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*
   st*kappa*rho^kappa*exp(pred0)/2
} else{
cum.haz<-c(t(weight)%*%matrix(exp((kappa-1)*log(nodes1)-log(1+(gam*nodes1)^kappa)+
   phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*
   st*kappa*rho^kappa*exp(pred0)/2
}
return(cum.haz)
}
########################################################################
#' Cox-Snell residuals for a jm.reg fit
#' @description Returns posterior summaries of the Cox-Snell residuals for a \code{\link{jm.reg}} fit.
#'   It also produces a residual plot.
#' @keywords Cox-Snell residual, Joint modeling, MCMC
#' @param jmfit a \code{\link{jm.reg}} fit with \code{method = "MCMC"}. 
#' @param posterior.mean returns posterior means of the residuals if \code{posterior.mean = TRUE}, otherwise
#'               returns posterior medians of the residuals.
#' @param n.sample an integer denoting how many MCMC samples to use. Default is \code{min(200, n.iter)}, where
#'      n.iter is the number of MCMC iterations (for all chains combined).  
#' @param xlim x limits (a vector of length 2 of the form \code{c(minimum, maximum)}) of the plot (optional).
#' @param ylim y limits (a vector of length 2 of the form \code{c(minimum, maximum)}) of the plot (optional).
#' @param xlab a label (character expression) for the x axis (optional).
#' @param ylab a label (character expression) for the y axis (optional).
#' @param main main title (character expression) of the plot (optional).
#' @param seed the seed for random number generation to draw \code{n.sample} samples.
#'      Default is \code{sample.int(.Machine$integer.max, 1)}.
#' @details A random sample of size \code{n.sample} is drawn from the posterior simulations, and the
#'   Cox-Snell residuals are computed for each sample. Following Rizopoulos and Ghosh (2011), the estimate of 
#'   the posterior expectation is obtained as the MCMC sample mean (or median) of the residuals 
#'   at each observed time point.
#'   A plot of residual vs. \eqn{-\textrm{log}\textrm{S}_{\textrm{KM}}}(redidual) is then produced (solid circles), 
#'   where \eqn{\textrm{S}_{\textrm{KM}}}(redidual) is the Kaplan-Meier estimate of the survivor function.
#'  Uncertainty in the plot can be assessed through the \code{n.sample} sets of residuals obtained
#'  from the posterior simulations (Zhou and Hanson, 2018). These are shown in the plot as grey lines.
#'  The plot also shows the unit slop line (in black). The plot of the Cox-Snell residuals should be roughly a straight line 
#'     with unit slope when the survival sub-model is adequate. 
#' @return \code{Cox-Snell Residuals:} posterior expectation of the residuals and \eqn{-\textrm{log}\textrm{S}_{\textrm{KM}}}(rediduals).
#' @return \code{plot:} see Details.
#' @references Rizopoulos D and Ghosh P, A Bayesian semiparametric multivariate joint model for multiple 
#'       longitudinal outcomes and a time‐to‐event, Statistics in Medicine, 30: 1366-1380, 2011.  
#' @references Zhou H and Hanson T, A unified framework for fitting Bayesian
#'       semiparametric models to arbitrarily censored survival data, including spatially referenced data, 
#'       Journal of the American Statistical Association, 113(522), 571-581, 2018.
#' @author Shahedul Khan <khan@math.usask.ca>
#' @seealso \code{\link{jm.icriteria}}, \code{\link{jm.plots}}, \code{\link{jm.reg}}
#' @examples 
#'   # Example: AIDS data from package 'JM'
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
#'   jm.resid(jmfit.wph0)
#' @export
############################################
# Cox-Snell residuals
#############################################
jm.resid<-function(jmfit,posterior.mean=TRUE,n.sample=NULL,
   xlim=NULL,ylim=NULL,xlab=NULL,ylab=NULL,main=NULL,
   seed=sample.int(.Machine$integer.max, 1)){
if(class(jmfit)!="JMR"){
stop("jmfit must inherit from class JMR")
}
results<-jmfit$all.results
if(is.null(results)){
stop("\nThis function works only for method = 'MCMC'")
}
dat<-results$stan.data
st<-dat$st
status<-dat$status
z<-dat$z
m<-dat$m
pp1<-dat$pp1
beta.loc<-results$beta.loc
phi.loc<-results$phi.loc
kappa.loc<-results$kappa.loc
rho.loc<-results$rho.loc
alpha.loc<-results$alpha.loc
sigma.loc<-results$sigma.loc
b.loc<-results$b.loc
method<-results$method
nodesf<-dat$nodesf
nodesr<-dat$nodesr
nodes1<-c(dat$nodes1)
weight<-dat$weights
qpoints<-dat$qpoints
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
if(results$surv.model=="weibull"){
gam.loc<-NULL
jmcumhaz<-jmcumhaz.weibull
} else if(results$surv.model=="llogistic"){
gam.loc<-NULL
jmcumhaz<-jmcumhaz.llogistic
} else if(results$surv.model=="lnormal"){
gam.loc<-NULL
jmcumhaz<-jmcumhaz.lnormal
} else if(results$surv.model=="weibullph"){
gam.loc<-NULL
jmcumhaz<-jmcumhaz.weibullph
} else if(results$surv.model=="eweibull"){
gam.loc<-results$gam.loc
jmcumhaz<-jmcumhaz.eweibull
} else if(results$surv.model=="ggamma"){
gam.loc<-results$gam.loc
jmcumhaz<-jmcumhaz.ggamma
} else{
gam.loc<-results$gam.loc
jmcumhaz<-jmcumhaz.gllogistic
}
mcmc.sample<-as.matrix(results$stan.fit)
n.iter<-nrow(mcmc.sample)
if(is.null(n.sample)){
n.sample<-min(200,n.iter)
}
set.seed(seed)
mcmc.sample1<-mcmc.sample[sample(n.iter,n.sample),]
#------- Residual Analysis -----------------
sam.dat<-apply(mcmc.sample1,1,jmcumhaz,st=st,m=m,pp1=pp1,z=z,x_fixed=x_fixed,
  beta.loc=beta.loc,phi.loc=phi.loc,kappa.loc=kappa.loc,
  gam.loc=gam.loc,rho.loc=rho.loc,alpha.loc=alpha.loc,b.loc=b.loc,
  which_alphafixed=which_alphafixed,which_alphatime=which_alphatime,
  nodes1=nodes1,nodesr=nodesr,nodesf=nodesf,weight=weight,
  qpoints=qpoints,method=method,lme.model=lme.model)
if(posterior.mean){
cumhaz<-apply(sam.dat,1,mean,na.rm=TRUE)
} else{
cumhaz<-apply(sam.dat,1,median,na.rm=TRUE)
}
H1<-list()
H2<-list()
for(i in 1:n.sample){
km.Rhat<-survfit(Surv(sam.dat[,i],status)~1,type="kaplan-meier")
H<-cbind(km.Rhat$time,-log(km.Rhat$surv))
H[which(!is.finite(H))]<-NA
H<-na.omit(H)
H1[[i]]<-H[,1]
H2[[i]]<-H[,2]
}
km.Rhat1<-survfit(Surv(cumhaz,status)~1,type="kaplan-meier")
H3<-cbind(km.Rhat1$time,-log(km.Rhat1$surv))
H3[which(!is.finite(H3))]<-NA
H3<-na.omit(H3)
if(is.null(xlim) || is.null(ylim)){
xlim<-ylim<-c(min(c(unlist(H1),unlist(H2))),max(c(unlist(H1),unlist(H2))))
}
if(is.null(xlab)){
xlab<-"Cox-Snell Residuals"
}
if(is.null(ylab)){
ylab<-as.expression(bquote('-log S'[KM]*'(Residuals)'))
}
if(is.null(main)){
par(mar = c(5,4.5,4,2) + 0.1)    
plot(H1[[1]],H2[[1]],type="l",lty=2,xlim=xlim,ylim=ylim,
  xlab=xlab,ylab=ylab,cex.lab=1.25,cex.main=1.5,col="grey85")
title<-list(bquote(atop("Cox-Snell Residual Plot for a Sample",
    "of" ~ .(n.sample)~"MCMC Simulations")))
mtext(do.call(expression, title),side=3,cex=1.5,line=0.5)
} else{
par(mar = c(5,4.5,4,2) + 0.1)    
plot(H1[[1]],H2[[1]],type="l",lty=2,xlim=xlim,ylim=ylim,
 xlab=xlab,ylab=ylab,cex.lab=1.25,cex.main=1.5,main=main,col="grey85")
}
for(i in 2:n.sample){
lines(H1[[i]],H2[[i]],type="l",lty=2,col="grey75")
}
abline(a=0,b=1,lwd=2)
lines(H3[,1],H3[,2],type="p",pch=20)
colnames(H3)<-c("residuals","-logS_KM(residuals)")
output<-list("Cox-Snell Residuals"=H3)
class(output)<-"JMR"
attr(output, "hidden") <-NULL
return(output)
}
