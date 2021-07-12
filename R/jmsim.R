#-------------
# Weibull PH and Generalized log-logistic PH
#-------------
jmsim.ph<-function(st,b,pred,alpha,data1,qpoints,nodes,weight,
   timevar,lme.formula,phi,rho,kappa,gam=NULL,logu,model,method){
nodes<-nodes
nodes0<-sapply(st,function(u){0.5*(u*nodes+u)})
newdat<-data1[rep(seq_len(nrow(data1)),qpoints),]
newdat[[timevar]]<-c(nodes0)
nodesf<-model.matrix(lme.formula,newdat)
nodesr<-cbind(1,nodes0)
if(method=="riz"){
if(model=="weibullph"){
logs<- -c(weight%*%matrix(exp((kappa-1)*log(nodes0)+phi*(nodesf %*% alpha + nodesr%*%b)),
   nrow=length(weight)))*st*kappa*rho^kappa*exp(pred)/2
} else{
logs<- -c(weight%*%matrix(exp((kappa-1)*log(nodes0)-log(1+(gam*nodes0)^kappa)+
   phi*(nodesf %*% alpha + nodesr%*%b)),nrow=length(weight)))*st*kappa*rho^kappa*exp(pred)/2
}
} else{
if(model=="weibullph"){
logs<- -c(weight%*%matrix(exp((kappa-1)*log(nodes0)+phi*(nodesr%*%b)),
   nrow=length(weight)))*st*kappa*rho^kappa*exp(pred)/2
} else{
logs<- -c(weight%*%matrix(exp((kappa-1)*log(nodes0)-log(1+(gam*nodes0)^kappa)+
   phi*(nodesr%*%b)),nrow=length(weight)))*st*kappa*rho^kappa*exp(pred)/2
}
}
return(logs-logu)
}
#-------------
# Log-normal AFT
#-------------
jmsim.lnormal<-function(st,beta,kappa,rho,alpha0,alphat,b,z,x,logu,form){
if(form=="riz"){
psi<-exp(-z%*%beta-phi*(x%*%alpha0+b[1]))*(1-exp(-phi*(alphat+b[2])*st))/(phi*(alphat+b[2]))
} else{
psi<-exp(-z%*%beta-phi*+b[1])*(1-exp(-phi*b[2]*st))/(phi*b[2])
}
logs<- plnorm(psi, meanlog = -log(rho), sdlog = 1/kappa, lower.tail = FALSE, log.p = TRUE)
return(logs-logu)
}
#-------------
# Generalized gamma AFT
#-------------
jmsim.ggamma<-function(st,beta,kappa,gam,rho,alpha0,alphat,b,z,x,logu,form){
if(form=="riz"){
psi<-exp(-z%*%beta-phi*(x%*%alpha0+b[1]))*(1-exp(-phi*(alphat+b[2])*st))/(phi*(alphat+b[2]))
} else{
psi<-exp(-z%*%beta-phi*+b[1])*(1-exp(-phi*b[2]*st))/(phi*b[2])
}
logs<-pgamma((rho*psi)^kappa,shape=gam,scale=1,lower.tail=FALSE,log.p=TRUE)
return(logs-logu)
}
#------------------------------------
# Main function
#-------------------------------------
#' Simulate from joint models
#' @description This function simulates longitudinal responses and event times from joint models.
#'     Available options for the time-to-event submodel are Weibull AFT, log-logistic AFT, log-normal AFT, 
#'     exponentiated Weibull AFT, generalized gamma AFT, Weibull PH and generalized log-logistic PH.
#'     For the longitudinal process, the linear mixed-effects model is assumed, 
#'     with \emph{a simple linear model for the random-effects component}
#'     (i.e., \eqn{b_0 + b_1} \code{times}, where \eqn{b_0} and \eqn{b_1} are the random intercept and random slope,
#'     respectively).
#' @keywords Joint model, Simulation
#' @param surv.formula a one-sided formula of the form \code{~ z1 + z2 + ... + zp}, where 
#'       \code{z1, z2, ..., zp} are baseline covariates for the time-to-event process. 
#' @param lme.formula a one-sided formula of the form \code{~ times} or \code{~ times + x2 + x3 + ... + xq}
#'       for the fixed-effects component of the longitudinal model, where 
#'       \code{times} is the time variable (time points at which longitudinal measurements are taken) and
#'       \code{x2, x3, ..., xq} are time-independent covariates.
#' @param surv.model the survival model to be used to describe the event process. Available options are given as follow.
#'        \itemize{ 
#'           \item AFT Models: \code{"weibull"}, \code{"lnormal"}, \code{"llogistic"}, \code{"eweibull"} and \code{"ggamma"}
#'                 for Weibull, log-normal, log-logistic, exponentiated Weibull and generalized gamma distributions, respectively. 
#'           \item PH Models: \code{"weibullph"} and \code{"gllph"} for Weibull and generalized log-logistic distributions, respectively.
#'         }
#'       The parameterization of these distributions is described in \code{\link{etime.reg}}.
#' @param form a character string to describe the formulation of the joint model. Available options are
#'      \code{form = "rizopoulos"} (can be abbreviated as "riz") for the formulation described by Rizopoulos (2012), and
#'      \code{form = "henderson"} (can be abbreviated as "hen") for the formulation proposed by Henderson et al. (2000). 
#'      See \code{\link{jm.reg}} for details.
#' @param surv.par a named list of the form \code{list(beta = a vector, phi = a scalar, logrho = a scalar,} 
#'      \code{logkappa = a scalar, loggamma = a scalar, status.rho = a scalar)} for the time-to-event model parameters, 
#'      where \code{beta} is the vector of regression coefficients for the baseline covariates as specified in
#'      \code{surv.formula} (without an intercept), \code{phi} is the association parameter, 
#'      \code{logrho} = log\eqn{\rho}, \code{logkappa} = log\eqn{\kappa}, \code{loggamma} = log\eqn{\gamma},
#'      and \code{status.rho} is the rate
#'      parameter of the exponential distribution from which censored data are generated.
#'      See \code{\link{etime.reg}} for parameterizations of the time-to-event distributions
#'      in terms of \eqn{\rho}, \eqn{\kappa} and \eqn{\gamma}. 
#' @param lme.par a named list of the form \code{list(alpha = a vector, sigma = a scalar, rand_cov = a 2 x 2 covariance matrix)} 
#'      for the longitudinal model parameters, 
#'      where \code{alpha} is the vector of fixed-effects regression coefficients as specified in
#'      \code{lme.formula} (includes an intercept), \code{sigma} is the 
#'      standard deviation of measurement errors (i.e., the value of \eqn{\sigma} for \eqn{\epsilon} ~ \eqn{N(0,\sigma^2)}),
#'      and \code{rand_cov} is the covariance matrix of the random effects \eqn{b_{0}} and \eqn{b_{1}}. 
#'      Note that only a simple linear model \eqn{b_{0} + b_{1}} \code{times} is allowed for random effects.
#' @param times a numeric vector for the time points at which longitudinal measurements are planned to be taken.
#' @param timevar the time variable (a character string) of the longitudinal model. For example,
#'        if \code{lme.formula = ~ times + x1}, then \code{timevar = "times"}.
#' @param Data a data frame containing the baseline covariates \code{z}'s for the time-to-event model and
#'        the time-independent covariates \code{x}'s for the longitudinal model. For \eqn{n} subjects, it must
#'        must have \eqn{n} rows.
#' @details This function simulates longitudinal responses and event times from joint models.
#'     The longitudinal model is of the form \eqn{y = x'(s)\alpha +w'(s)b + \epsilon},
#'     where \eqn{x'(s)\alpha=\alpha_0 + \alpha_1 s} or
#'     \eqn{\alpha_0 + \alpha_1 s + \alpha_2 x_2 + \alpha_3 x_3 + ... + \alpha_q x_q},
#'     \eqn{w'(s)b = b_0 + b_1 s}, \eqn{\epsilon} ~ \eqn{N(0,\sigma^2)} and \eqn{b} ~ \eqn{N(0, \Sigma_b)}.
#'     For an AFT event process, the survivor function is of the form
#'     \eqn{S_0(\int_0^t} exp\eqn{[-(z'\beta+\phi(x'(u)\alpha +w'(u)b)]du)} for \code{form = "rizopoulos"}, 
#'     and \eqn{S_0(\int_0^t} exp\eqn{[-(z'\beta+\phi w'(u)b)]du)} for \code{form = "henderson"}.
#'     For a PH event process, the survivor function is of the form exp[\eqn{-\Lambda(t)}], where
#'     \eqn{\Lambda(t)=\int_0^t \lambda_0(u)} exp\eqn{[z'\beta+\phi(x'(u)\alpha +w'(u)b)]du} for \code{form = "rizopoulos"},
#'     and \eqn{\Lambda(t)=\int_0^t \lambda_0(u)} exp\eqn{[z'\beta+\phi w'(u)b]du} for \code{form = "henderson"}.
#'     Event times are drawn by solving the equation \code{survivor function = U} for \eqn{t}, where
#'     \code{U} is a Uniform(0, 1) variable. Note that
#'     the amount of censoring can be controlled by changing the value of \code{status.rho}. 
#' @return A list with components:
#' @return \code{b:} random effects (a matrix).
#' @return \code{long.data:} simulated longitudinal data.
#' @return \code{surv.data:} simulated event time data.
#' @references Henderson R, Diggle P, and Dobson A, Joint modelling of longitudinal measurements and event time data, 
#'           Biostatistics, 1: 465-480, 2000.
#' @references Rizopoulos D, Joint Models for Longitudinal and Time-to-Event Data: 
#'          With Applications in R, Chapman and Hall/CRC, 2012.
#' @author Shahedul Khan <khan@math.usask.ca>
#' @seealso \code{\link{jm.reg}}, \code{\link{jm.summary}}
#' @examples 
#'   n <- 100
#'   surv.formula <- ~ z1 + z2
#'   lme.formula <- ~ times + x1 + x2
#'   cor.mat <- matrix(c(1, -0.2, -0.2, 1), ncol = 2, byrow = TRUE)
#'   sd <- c(0.5, 0.1)
#'   rand_cov <- (sd %*% t(sd)) * cor.mat
#'   lme.par <- list(alpha = c(10, -1, -1, 0.25), sigma = 0.5, rand_cov = rand_cov)
#'   surv.par <- list(beta = c(-0.5, 0.5), phi = -0.25, logrho = log(0.01),
#'     logkappa = log(2), status.rho = 0.003)
#'   z1 <- rbinom(n, 1, 0.5)
#'   z2 <- rnorm(n)
#'   # Data with all covariate information
#'   Data <- data.frame(z1 = z1, z2 = z2, x1 = z1, x2 = z2)
#'   # Time points at which longitudinal measurements are planned to be taken
#'   times <- c(0, 1, 2, 3, 4)
#'   timevar <- "times"
#'   jmsim.dat <- jm.sim(surv.formula = surv.formula, lme.formula = lme.formula,
#'     surv.model = "weibullph", form = "hen", surv.par = surv.par,
#'     lme.par = lme.par, times = times, timevar = timevar, Data = Data)
#'   surv.data <- jmsim.dat$surv.data
#'   long.data <- jmsim.dat$long.data
#'   surv.fit <- coxph(Surv(st, status) ~ z1 + z2, data = surv.data, x = TRUE)
#'   lme.fit <- lme(y ~ times + x1 + x2, random =  ~ times | id, data = long.data)
#'   swph <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "weibullph", timevar = "times", form = "hen")
#'   jm.summary(swph)
#' @export

jm.sim<-function(surv.formula,lme.formula,surv.model,form,
 surv.par,lme.par,times,timevar,Data){
if(isFALSE(surv.model=="weibull" || surv.model=="llogistic" || 
   surv.model=="eweibull" || surv.model=="ggamma" || 
   surv.model=="lnormal" || surv.model=="weibullph" || 
   surv.model=="gllph")){
stop("\nsurv.model must be 'weibull' or 'llogistic' or 'eweibull' or
 'ggamma' or 'lnormal' or 'weibullph' or 'gllph'")
}
survvar<-all.vars(surv.formula)
if(isFALSE(all(survvar %in% colnames(Data)))){
stop("\nvariable names in 'surv.formula' did not match
 with the variable names in 'Data'")
}
surv.par.names<-names(surv.par)
if(surv.model=="weibull" || surv.model=="llogistic" || 
  surv.model=="lnormal" || surv.model=="weibullph"){
survpar<-c("beta","phi","logrho","logkappa","status.rho")
if(isFALSE(all(survpar %in% surv.par.names))){
stop("use names to define 'surv.par' list; at least one
  element of 'surv.par' is missing or not recognized")
}
}
if(surv.model=="eweibull" || surv.model=="ggamma" || surv.model=="gllph"){
survpar<-c("beta","phi","logrho","logkappa","loggamma","status.rho")
if(isFALSE(all(survpar %in% surv.par.names))){
stop("use names to define 'surv.par' list; at least one
  element of 'surv.par' is missing or not recognized")
}
}  
if(length(surv.par$beta)!=length(survvar)){
stop("length of beta in 'surv.par' did not match with 'surv.formula'")
}
lmevar<-all.vars(lme.formula)
lme.par.names<-names(lme.par)
lmepar<-c("alpha","sigma","rand_cov")
if(isFALSE(all(lmepar %in% lme.par.names))){
stop("\nuse names to define 'lme.par' list; at least one
 element of 'lme.par' is missing or not recognized")
}
if(length(lme.par$alpha)!=(length(lmevar)+1)){
stop("\nlength of alpha in 'lme.par' did not match with 'lme.formula'")
}
if(lme.par$sigma<=0){
stop("\nsigma in 'lme.par' must be positive")
}
if(isFALSE(all(dim(lme.par$rand_cov)==c(2,2)))){
stop("\nrand_cov in 'lme.par' must be a 2 x 2 positive definite matrix")
}
if(any(eigen(lme.par$rand_cov)$values<=0) || lme.par$rand_cov[1,2] != lme.par$rand_cov[2,1] || 
  any(diag(lme.par$rand_cov)<=0)){
stop("\nrand_cov is not a covariance matrix")
}
if(isFALSE(all(lmevar %in% c(colnames(Data),timevar)))){
stop("\nvariable names in 'lme.formula' did not match
 with the variable names in data; check your variable
 names in 'Data', 'timevar' and 'lme.formula'")
}
beta<-surv.par$beta
phi<-surv.par$phi
rho<-exp(surv.par$logrho)
kappa<-exp(surv.par$logkappa)
if(surv.model=="eweibull" || surv.model=="ggamma" || surv.model=="gllph"){
gam<-exp(surv.par$loggamma)
} else{
gam<-NULL
}
alpha<-lme.par$alpha
rand_cov<-lme.par$rand_cov
alphat.loc<-which(c("int",lmevar)==timevar)
alphat<-alpha[alphat.loc]
alpha0<-alpha[-alphat.loc]
rform<-all.vars(lme.formula)[-which(lmevar==timevar)]
if(length(rform)>1){
rform1<-paste(rform,collapse= "+")
} 
if(length(rform)==1){
rform1<-rform
}
if(length(rform)>=1){
fixed.formula<-as.formula(paste("~",rform1,collapse=""))
}
if(length(rform)==0){
fixed.formula<-as.formula(paste("~","1",collapse=""))
}
fixed.dat<-matrix(model.matrix(fixed.formula,Data),ncol=(length(alpha)-1))
surv.mat<-matrix(model.matrix(surv.formula,data=Data)[,-1],ncol=length(survvar))
if(form=="riz" || form=="rizopoulos"){
method<-"riz"
} else{
method<-"hen"
}
nb<-nrow(rand_cov)
b.mat<-matrix(0,nrow=n,ncol=nb)
st<-NULL
#####
# Weibull AFT
#####
if(surv.model=="weibull"){
if(method=="riz"){
for(i in 1:n){
repeat{
u<-runif(1)
b<-sim.mnorm(1,rep(0,2),rand_cov)
din<-phi*(alphat+b[2])
t0<-c(exp(surv.mat[i,]%*%beta+phi*(fixed.dat[i,]%*%alpha0+b[1]))*phi*(alphat+b[2])*(-log(u))^(1/kappa)/rho)
if((1-t0)>0){
if((-log1p(-t0)/din)>0){
break
}
} 
}
st[i]<--log1p(-t0)/din
b.mat[i,]<-b
}
} else{
for(i in 1:n){
repeat{
u<-runif(1)
b<-sim.mnorm(1,rep(0,2),rand_cov)
din<-phi*b[2]
t0<-c(exp(surv.mat[i,]%*%beta+phi*b[1])*phi*b[2]*(-log(u))^(1/kappa)/rho)
if((1-t0)>0){
if((-log1p(-t0)/din)>0){
break
}
}
}
st[i]<--log1p(-t0)/din
b.mat[i,]<-b
}
}
}
#####
# Log-logistic AFT
#####
if(surv.model=="llogistic"){
if(method=="riz"){
for(i in 1:n){
repeat{
u<-runif(1)
b<-sim.mnorm(1,rep(0,2),rand_cov)
din<-phi*(alphat+b[2])
t0<-c(exp(surv.mat[i,]%*%beta+phi*(fixed.dat[i,]%*%alpha0+b[1]))*phi*(alphat+b[2])*((1-u)/u)^(1/kappa)/rho)
if((1-t0)>0){
if((-log1p(-t0)/din)>0){
break
}
}
}
st[i]<--log1p(-t0)/din
b.mat[i,]<-b
}
} else{
for(i in 1:n){
repeat{
u<-runif(1)
b<-sim.mnorm(1,rep(0,2),rand_cov)
din<-phi*b[2]
t0<-c(exp(surv.mat[i,]%*%beta+phi*b[1])*phi*b[2]*((1-u)/u)^(1/kappa)/rho)
if((1-t0)>0){
if((-log1p(-t0)/din)>0){
break
}
}
}
st[i]<--log1p(-t0)/din
b.mat[i,]<-b
}
}
}
#####
# Exponentiated Weibull AFT
#####
if(surv.model=="eweibull"){
if(method=="riz"){
for(i in 1:n){
repeat{
u<-runif(1)
b<-sim.mnorm(1,rep(0,2),rand_cov)
din<-phi*(alphat+b[2])
ter<-(-log1p(-exp(log1p(-u)/gam)))^(1/kappa)
t0<-c(exp(surv.mat[i,]%*%beta+phi*(fixed.dat[i,]%*%alpha0+b[1]))*phi*(alphat+b[2])*ter/rho)
if((1-t0)>0){
if((-log1p(-t0)/din)>0){
break
}
}
}
st[i]<--log1p(-t0)/din
b.mat[i,]<-b
}
} else{
for(i in 1:n){
repeat{
u<-runif(1)
b<-sim.mnorm(1,rep(0,2),rand_cov)
din<-phi*b[2]
ter<-(-log1p(-exp(log1p(-u)/gam)))^(1/kappa)
t0<-c(exp(surv.mat[i,]%*%beta+phi*b[1])*phi*b[2]*ter/rho)
if((1-t0)>0){
if((-log1p(-t0)/din)>0){
break
}
}
}
st[i]<--log1p(-t0)/din
b.mat[i,]<-b
}
}
}
#####
# Generalized gamma AFT
#####
if(surv.model=="ggamma"){
for(i in 1:n){
repeat{
u<-runif(1)
b<-sim.mnorm(1,rep(0,2),rand_cov)
options(warn=-1)
root<-tryCatch({uniroot(f=jmsim.ggamma,lower=0,upper=max(times)+100,
  beta=beta,kappa=kappa,gam=gam,rho=rho,alpha0=alpha0,alphat=alphat,
  b=b,z=surv.mat[i,],x=fixed.dat[i,],logu=log(u),form=method)},
  error = function(e){NULL})
options(warn=0)
if(!is.null(root)){
if(root$root>0){
break
}
}
}
st[i]<-root$root
b.mat[i,]<-b
}
}
#####
# Log-normal AFT
#####
if(surv.model=="lnormal"){
for(i in 1:n){
repeat{
u<-runif(1)
b<-sim.mnorm(1,rep(0,2),rand_cov)
options(warn=-1)
root<-tryCatch({uniroot(f=jmsim.lnormal,lower=0,upper=max(times)+100,
  beta=beta,kappa=kappa,rho=rho,alpha0=alpha0,alphat=alphat,
  b=b,z=surv.mat[i,],x=fixed.dat[i,],logu=log(u),form=method)},
  error = function(e){NULL})
options(warn=0)
if(!is.null(root)){
if(root$root>0){
break
}
}
}
st[i]<-root$root
b.mat[i,]<-b
}
}
#####
# Weibull and GLL PH
#####
if(surv.model=="weibullph" || surv.model=="gllph"){
pred<-c(surv.mat%*%beta)
qpoints<-21
nodes<-GK.nodes(qpoints)$gk.n
weight<-GK.nodes(qpoints)$gk.w
for(i in 1:n){
repeat{
u<-runif(1)
#repeat{
b<-sim.mnorm(1,rep(0,2),rand_cov)
#if(method="riz"){
#din<-phi*(alphat+b[2])
#} else{
#din<-phi*b[2]
#}
#if(din>0) break
#}
options(warn=-1)
root<-tryCatch({uniroot(f=jmsim.ph,lower=0,upper=max(times)+100,
   b=b,pred=pred[i],alpha=alpha,data1=Data[i,],
   qpoints=qpoints,nodes=nodes,weight=weight,
   timevar=timevar,lme.formula=lme.formula,
   phi=phi,rho=rho,kappa=kappa,gam=gam,logu=log(u),
   model=surv.model,method=method)},
   error = function(e){NULL})
options(warn=0)
if(!is.null(root)){
if(root$root>0){
break
}
}
}
st[i]<-root$root
b.mat[i,]<-b
}
}
#-----------
st1<-st
status.rho<-surv.par$status.rho
ctimes<-rexp(n,status.rho)
#ctimes<-runif(n,min(times),end.time)
st<-pmin(st1,ctimes)
status<-as.numeric(st1<=ctimes)
#-----------
length.t<-length(times)
rep.t<-rep(times,n)
id<-rep(1:n,each=length.t)
dat1<-data.frame(id=id,Data[rep(seq_len(nrow(Data)),each=length.t),],rep.t)
colnames(dat1)[ncol(dat1)]<-timevar
colnames(dat1)[2:(ncol(dat1)-1)]<-colnames(Data)
xlong<- model.matrix(lme.formula, data=dat1)
zlong<- cbind(1,rep.t)
bb<-b.mat[rep(seq_len(nrow(b.mat)),each=length.t),]
sigma<-lme.par$sigma
eps<-rnorm(length(id),mean=0,sd=sigma)
y<-c(xlong%*%alpha)+as.vector(rowSums(bb*zlong))+eps
dat2<-data.frame(dat1,y=y,st=rep(st,each=length.t),status=rep(status,each=length.t))
long.data<-dat2[dat2[,"st"]>=dat2[,timevar],]
rownames(long.data)<-NULL
surv.data<-data.frame(id=unique(id),Data,st=st,status=status)
output<-list(b=b.mat,long.data=long.data,surv.data=surv.data)
class(output)<-"JMR"
attr(output, "hidden") <-c("b")
return(output)
}
######################################



