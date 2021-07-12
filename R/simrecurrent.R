#########################################################
# Data Simulation: EW AFT
#########################################################
sim0.ewaft<-function(t0,xx,kappa,gam,rho,beta){
pred<-sum(xx*beta)
u<-rho*t0*exp(-pred)
repeat{
e<-rexp(1)
r0<-log1mexp(-u^kappa)
rr0<-log1p(-exp(gam*r0))
r1<-ifelse(is.finite(rr0),(rr0-e),(log(gam)-u^kappa-e))
s0<-log1mexp(r1)
ss0<-log1p(-exp(s0/gam))
s1<-ifelse(is.finite(ss0),(-ss0)^(1/kappa),(-(log(1/gam)+r1))^(1/kappa))
w<-as.numeric(s1*exp(pred)/rho-t0)
if(is.finite(w)) break
}
return(w)
}
####################################
#sim0.ewaft<-function(t0,xx,kappa,gam,rho,beta){
#pred<-sum(xx*beta)
#u<-rho*t0*exp(-pred)
#repeat{
#e<-rexp(1)
#r1<-log1p(-exp(gam*log1p(-exp(-u^kappa))))-e
#r2<-exp(log1p(-exp(r1))/gam)
#r3<-(-log1p(-r2))^(1/kappa)
#w<-as.numeric(r3*exp(pred)/rho-t0)
#if(is.finite(w)) break
#}
#return(w)
#}
#############################################
sim1.ewaft<-function(t0,xx,end.time,kappa,gam,rho,beta){
t00<-t0
t<-NULL
repeat{
tt0<-t0
t1<-t0+sim0.ewaft(t0=t0,xx=xx,kappa=kappa,gam=gam,rho=rho,beta=beta)
t<-c(t,t1)
t0<-t1
if(t0>end.time || tt0==t1) break
}
status<-rep(1,length(t))
status[which(t>end.time)]<-0
t[t>end.time]<-end.time
if(length(t)==1){
start<-t00
} else{
start<-c(t00,t[1:(length(t)-1)])
}
data.frame(start=start,stop=t,status=status,t(xx))
}
###################################################
sim.ewaft<-function(origin,end.time,design.mat,beta,kappa,gam,rho){
n<-length(end.time)
k<-0
sim.data<-NULL
repeat{
k<-k+1
dat0<-sim1.ewaft(t0=origin[k],xx=design.mat[k,],end.time=end.time[k],
   kappa=kappa,gam=gam,rho=rho,beta=beta)
sim.data<-rbind(sim.data,cbind(id=k,dat0))
if(k==n) break
}
sim.data<-sim.data[sim.data[,2]!=sim.data[,3],]
rownames(sim.data)<-NULL
return(data.frame(sim.data))
}
#########################################################
#########################################################
# Data Simulation: ggamma AFT
#########################################################
#########################################################
int.ggaft<-function(st,xx,rho,kappa,gam,beta){
pred<-sum(xx*beta)
u<-st*exp(-pred)
logf<-log(kappa)+kappa*gam*log(rho)+(kappa*gam-1)*log(u)-(rho*u)^kappa-lgamma(gam)-pred
logs<-pgamma((rho*u)^kappa,shape=gam,rate=1,lower.tail=FALSE,log=TRUE)
logh<-logf-logs
return(exp(logh))
}
################################
solve.ggaft<-function(w,t0,xx,rho,kappa,gam,beta,e){
t1<-t0+w
H<-integrate(f=int.ggaft,lower=t0,upper=t1,xx=xx,rho=rho,kappa=kappa,gam=gam,beta=beta)$value
equ<-log(H)-log(e)
#return(equ^2)
return(equ)
}
#########################################
sim0.ggaft<-function(t0,xx,kappa,gam,rho,beta,max.time){
#repeat{
#e<-rexp(1)
#init<-0.1
#root<-tryCatch({nlminb(start=init,objective=solve.ggaft,lower=0,
#  t0=t0,xx=xx,rho=rho,kappa=kappa,gam=gam,beta=beta,e=e)},
#  error = function(e){NULL})
#if(!is.null(root)) {
#if(is.finite(root$par) && is.finite(root$objective)){
#if(root$convergence==0 && root$par !=init && root$objective>0 &&
#   sqrt(root$objective)<1e-8 && root$iterations>1){
#w<-root$par
#break}}}
#}
repeat{
e<-rexp(1)
options(warn=-1)
root<-tryCatch({uniroot(f=solve.ggaft,lower=1e-8,upper=max.time+100,
  t0=t0,xx=xx,rho=rho,kappa=kappa,gam=gam,beta=beta,e=e)},
  error = function(u){NULL})
options(warn=0)
if(!is.null(root)){
if(root$root>0){
w<-root$root
break}}
}
return(w)
}
###########################################
sim1.ggaft<-function(t0,xx,end.time,max.time,kappa,gam,rho,beta){
t00<-t0
t<-NULL
repeat{
tt0<-t0
t1<-t0+sim0.ggaft(t0=t0,xx=xx,kappa=kappa,gam=gam,rho=rho,beta=beta,max.time=max.time)
t<-c(t,t1)
t0<-t1
if(t0>end.time || tt0==t1) break
}
status<-rep(1,length(t))
status[which(t>end.time)]<-0
t[t>end.time]<-end.time
if(length(t)==1){
start<-t00
} else{
start<-c(t00,t[1:(length(t)-1)])
}
data.frame(start=start,stop=t,status=status,t(xx))
}
###################################################
sim.ggaft<-function(origin,end.time,design.mat,beta,kappa,gam,rho){
n<-length(end.time)
max.time<-max(end.time)
k<-0
sim.data<-NULL
repeat{
k<-k+1
dat0<-sim1.ggaft(t0=origin[k],xx=design.mat[k,],end.time=end.time[k],max.time=max.time,
   kappa=kappa,gam=gam,rho=rho,beta=beta)
sim.data<-rbind(sim.data,cbind(id=k,dat0))
if(k==n) break
}
sim.data<-sim.data[sim.data[,2]!=sim.data[,3],]
rownames(sim.data)<-NULL
return(data.frame(sim.data))
}
#########################################################
#########################################################
# Data Simulation: lnormal AFT
#########################################################
#########################################################
int.lnaft<-function(st,xx,rho,kappa,beta){
pred<-sum(xx*beta)
u<-st*exp(-pred)
logf<-dlnorm(u,(-log(rho)),(1/kappa),log=TRUE)-pred
logs<-plnorm(u,(-log(rho)),(1/kappa),lower.tail = FALSE,log.p=TRUE)
logh<-logf-logs
return(exp(logh))
}
################################
solve.lnaft<-function(w,t0,xx,rho,kappa,beta,e){
t1<-t0+w
H<-integrate(f=int.lnaft,lower=t0,upper=t1,xx=xx,rho=rho,kappa=kappa,beta=beta)$value
equ<-log(H)-log(e)
#return(equ^2)
return(equ)
}
########################################
sim0.lnaft<-function(t0,xx,kappa,rho,beta,max.time){
#repeat{
#e<-rexp(1)
#init<-0.1
#root<-tryCatch({nlminb(start=init,objective=solve.lnaft,lower=0,
#  t0=t0,xx=xx,rho=rho,kappa=kappa,beta=beta,e=e)},
#  error = function(e){NULL})
#if(!is.null(root)) {
#if(is.finite(root$par) && is.finite(root$objective)){
#if(root$convergence==0 && root$par !=init && root$objective>0 &&
#   sqrt(root$objective)<1e-8 && root$iterations>1){
#w<-root$par
#break}}}
#}
repeat{
e<-rexp(1)
options(warn=-1)
root<-tryCatch({uniroot(f=solve.lnaft,lower=1e-8,upper=max.time+100,
  t0=t0,xx=xx,rho=rho,kappa=kappa,beta=beta,e=e)},
  error = function(u){NULL})
options(warn=0)
if(!is.null(root)){
if(root$root>0){
w<-root$root
break}}
}
return(w)
}
###########################################
sim1.lnaft<-function(t0,xx,end.time,max.time,kappa,rho,beta){
t00<-t0
t<-NULL
repeat{
tt0<-t0
t1<-t0+sim0.lnaft(t0=t0,xx=xx,kappa=kappa,rho=rho,beta=beta,max.time=max.time)
t<-c(t,t1)
t0<-t1
if(t0>end.time || tt0==t1) break
}
status<-rep(1,length(t))
status[which(t>end.time)]<-0
t[t>end.time]<-end.time
if(length(t)==1){
start<-t00
} else{
start<-c(t00,t[1:(length(t)-1)])
}
data.frame(start=start,stop=t,status=status,t(xx))
}
###################################################
sim.lnaft<-function(origin,end.time,design.mat,beta,kappa,rho){
n<-length(end.time)
max.time<-max(end.time)
k<-0
sim.data<-NULL
repeat{
k<-k+1
dat0<-sim1.lnaft(t0=origin[k],xx=design.mat[k,],end.time=end.time[k],max.time=max.time,
   kappa=kappa,rho=rho,beta=beta)
sim.data<-rbind(sim.data,cbind(id=k,dat0))
if(k==n) break
}
sim.data<-sim.data[sim.data[,2]!=sim.data[,3],]
rownames(sim.data)<-NULL
return(data.frame(sim.data))
}
########################################################
#########################################################
# Data Simulation: W AFT
#########################################################
#########################################################
sim0.waft<-function(t0,xx,kappa,rho,beta){
pred<-sum(xx*beta)
u<-rho*t0*exp(-pred)
repeat{
e<-rexp(1)
w<-(u^kappa+e)^(1/kappa)*exp(pred)/rho-t0
if(is.finite(w)) break
}
return(w)
}
###################################
sim1.waft<-function(t0,xx,end.time,kappa,rho,beta){
t00<-t0
t<-NULL
repeat{
tt0<-t0
t1<-t0+sim0.waft(t0=t0,xx=xx,kappa=kappa,rho=rho,beta=beta)
t<-c(t,t1)
t0<-t1
if(t0>end.time || tt0==t1) break
}
status<-rep(1,length(t))
status[which(t>end.time)]<-0
t[t>end.time]<-end.time
if(length(t)==1){
start<-t00
} else{
start<-c(t00,t[1:(length(t)-1)])
}
data.frame(start=start,stop=t,status=status,t(xx))
}
######################################
sim.waft<-function(origin,end.time,design.mat,beta,kappa,rho){
n<-length(end.time)
k<-0
sim.data<-NULL
repeat{
k<-k+1
dat0<-sim1.waft(t0=origin[k],xx=design.mat[k,],end.time=end.time[k],
   kappa=kappa,rho=rho,beta=beta)
sim.data<-rbind(sim.data,cbind(id=k,dat0))
if(k==n) break
}
sim.data<-sim.data[sim.data[,2]!=sim.data[,3],]
rownames(sim.data)<-NULL
return(data.frame(sim.data))
}
########################################################
#########################################################
# Data Simulation: LL AFT
#########################################################
##########################################################
sim0.llaft<-function(t0,xx,kappa,rho,beta){
pred<-sum(xx*beta)
u<-rho*t0*exp(-pred)
repeat{
e<-rexp(1)
w<-(expm1(log1p(u^kappa)+e))^(1/kappa)*exp(pred)/rho-t0
if(is.finite(w)) break
}
return(w)
}
###################################
sim1.llaft<-function(t0,xx,end.time,kappa,rho,beta){
t00<-t0
t<-NULL
repeat{
tt0<-t0
t1<-t0+sim0.llaft(t0=t0,xx=xx,kappa=kappa,rho=rho,beta=beta)
t<-c(t,t1)
t0<-t1
if(t0>end.time || tt0==t1) break
}
status<-rep(1,length(t))
status[which(t>end.time)]<-0
t[t>end.time]<-end.time
if(length(t)==1){
start<-t00
} else{
start<-c(t00,t[1:(length(t)-1)])
}
data.frame(start=start,stop=t,status=status,t(xx))
}
######################################
sim.llaft<-function(origin,end.time,design.mat,beta,kappa,rho){
n<-length(end.time)
k<-0
sim.data<-NULL
repeat{
k<-k+1
dat0<-sim1.llaft(t0=origin[k],xx=design.mat[k,],end.time=end.time[k],
   kappa=kappa,rho=rho,beta=beta)
sim.data<-rbind(sim.data,cbind(id=k,dat0))
if(k==n) break
}
sim.data<-sim.data[sim.data[,2]!=sim.data[,3],]
rownames(sim.data)<-NULL
return(data.frame(sim.data))
}
########################################################
#########################################################
# Data Simulation: W PH
#########################################################
#########################################################
sim0.wph<-function(t0,xx,kappa,rho,beta){
pred<-sum(xx*beta)
repeat{
e<-rexp(1)
w<-(((rho*t0)^kappa+e*exp(-pred))^(1/kappa))/rho-t0
if(is.finite(w)) break
}
return(w)
}
###################################
sim1.wph<-function(t0,xx,end.time,kappa,rho,beta){
t00<-t0
t<-NULL
repeat{
tt0<-t0
t1<-t0+sim0.wph(t0=t0,xx=xx,kappa=kappa,rho=rho,beta=beta)
t<-c(t,t1)
t0<-t1
if(t0>end.time || tt0==t1) break
}
status<-rep(1,length(t))
status[which(t>end.time)]<-0
t[t>end.time]<-end.time
if(length(t)==1){
start<-t00
} else{
start<-c(t00,t[1:(length(t)-1)])
}
data.frame(start=start,stop=t,status=status,t(xx))
}
######################################
sim.wph<-function(origin,end.time,design.mat,beta,kappa,rho){
n<-length(end.time)
k<-0
sim.data<-NULL
repeat{
k<-k+1
dat0<-sim1.wph(t0=origin[k],xx=design.mat[k,],end.time=end.time[k],
   kappa=kappa,rho=rho,beta=beta)
sim.data<-rbind(sim.data,cbind(id=k,dat0))
if(k==n) break
}
sim.data<-sim.data[sim.data[,2]!=sim.data[,3],]
rownames(sim.data)<-NULL
return(data.frame(sim.data))
}
#########################################################
#########################################################
# Data Simulation: GLL PH
#########################################################
#########################################################
sim0.gllph<-function(t0,xx,kappa,gam,rho,beta){
pred<-sum(xx*beta)
repeat{
e<-rexp(1)
r<-expm1((gam/rho)^kappa*e*exp(-pred)+log1p((gam*t0)^kappa))
w<-r^(1/kappa)/gam-t0
if(is.finite(w)) break
}
return(w)
}
#############################################
sim1.gllph<-function(t0,xx,end.time,kappa,gam,rho,beta){
t00<-t0
t<-NULL
repeat{
tt0<-t0
t1<-t0+sim0.gllph(t0=t0,xx=xx,kappa=kappa,gam=gam,rho=rho,beta=beta)
t<-c(t,t1)
t0<-t1
if(t0>end.time || tt0==t1) break
}
status<-rep(1,length(t))
status[which(t>end.time)]<-0
t[t>end.time]<-end.time
if(length(t)==1){
start<-t00
} else{
start<-c(t00,t[1:(length(t)-1)])
}
data.frame(start=start,stop=t,status=status,t(xx))
}
###################################################
sim.gllph<-function(origin,end.time,design.mat,beta,kappa,gam,rho){
n<-length(end.time)
k<-0
sim.data<-NULL
repeat{
k<-k+1
dat0<-sim1.gllph(t0=origin[k],xx=design.mat[k,],end.time=end.time[k],
   kappa=kappa,gam=gam,rho=rho,beta=beta)
sim.data<-rbind(sim.data,cbind(id=k,dat0))
if(k==n) break
}
sim.data<-sim.data[sim.data[,2]!=sim.data[,3],]
rownames(sim.data)<-NULL
return(data.frame(sim.data))
}
#########################
##########################
##########################
#' Simulate recurrent event data
#' @description Simulate recurrent event data for fixed covariates
#' @keywords AFT model, Recurrent event data, Simulation
#' @param surv.formula a one-sided formula of the form \code{~ z1 + z2 + ... + zp}, where 
#'       \code{z1, z2, ..., zp} are baseline covariates for the time-to-event process. 
#' @param surv.par a list of the form \code{list(beta = a vector, logrho = a scalar, logkappa = a scalar)}
#'      for the Weibull AFT, log-logistic AFT, log-normal AFT and Weibull PH models, and
#'      \code{list(beta = a vector, logrho = a scalar, logkappa = a scalar, loggamma = a scalar)}  
#'      for the exponentiated Weibull AFT, generalized gamma AFT and generalized log-logistic PH models,
#'      where \code{beta} is the vector of regression coefficients for the baseline covariates as specified in
#'      \code{surv.formula} (without an intercept), and \code{logrho}, \code{logkappa} and \code{loggamma} are the \code{log} of the
#'      parameters \eqn{\rho}, \eqn{\kappa} and \eqn{\gamma} (see \code{\link{etime.reg}} for 
#'      a description of the distribution parameterizations). Note that \code{surv.par}
#'      must be a named list as shown above. 
#' @param surv.model the model to be used to describe the recurrent event process. Available options are
#'       \code{"weibull"}, \code{"llogistic"}, \code{"lnormal"}, \code{"eweibull"} and \code{"ggamma"}
#'       for the AFT model, and \code{"weibullph"} and \code{"gllph"} for the PH model.
#' @param Data a data frame containing the covariates \code{z1, z2, ..., zp}  for the time-to-event model.
#'        For \eqn{n} subjects, it must must have \eqn{n} rows.
#' @param origin the time origin (an \eqn{n} x 1 vector for \eqn{n} subjects; typically, a vector of zeros).
#' @param end.time the end of follow-up time (an \eqn{n} x 1 vector for \eqn{n} subjects).
#' @details Simulate recurrent event data. See Section 2.4 of Cook and Lawless.
#' @return a data frame containing \code{id}, \code{start}, \code{stop}, \code{status} and covariates.
#' @references Cook RJ and Lawless J, The statistical analysis of recurrent events, Springer, 2007.
#' @author Shahedul Khan <khan@math.usask.ca>
#' @seealso \code{\link{etime.reg}}
#' @examples
#'   n <- 100 # number of subjects
#'   Data <- data.frame(z1 = rnorm(n), z2 = rbinom(n, 1, 0.5)) # data frame
#'   origin <- rep(0, n) # time origin for n subjects
#'   end.time <- runif(n, 10, 15) # end of follow-up time
#'   # Simulate data from the generalized gamma AFT model
#'   surv.par <- list(beta = c(0.5, -0.5), logrho = log(0.025),
#'       logkappa = log(2), loggamma = log(0.4))
#'   Data.gg <- sim.rec( ~ z1 + z2, surv.par, "ggamma", Data, origin, end.time)
#'   # Fit the generalized gamma AFT model
#'   fit.gg <- etime.reg(c(start, stop, status) ~ z1 + z2,
#'               surv.model = "ggamma", Data = Data.gg)
#'   fit.gg
#'   # Simulate data from the Weibull PH model
#'   surv.par <- list(beta = c(0.5, -0.5), logrho = log(0.1), logkappa = log(2))
#'   Data.wph <- sim.rec( ~ z1 + z2, surv.par, "weibullph", Data, origin, end.time)
#'   # Fit the Weibull PH model
#'   fit.wph <- etime.reg(c(start, stop, status) ~ z1 + z2,
#'               surv.model = "weibullph", Data = Data.wph)
#'   fit.wph
#' @export

sim.rec<-function(surv.formula,surv.par,surv.model,Data,origin,end.time){
if(isFALSE(surv.model=="weibull" || surv.model=="llogistic" || 
   surv.model=="eweibull" || surv.model=="lnormal" || surv.model=="ggamma" ||
   surv.model=="weibullph" || surv.model=="gllph")){
stop("surv.model must be one of 'weibull' or 'llogistic' or 'lnormal' or
     'eweibull' or 'ggamma' or 'weibullph' or 'gllph'")
}
survvar<-all.vars(surv.formula)
if(isFALSE(all(survvar %in% colnames(Data)))){
stop("variable names in 'surv.formula' did not match
     with the variable names in 'Data'")
}
surv.par.names<-names(surv.par)
if(surv.model=="weibull" || surv.model=="llogistic" || surv.model=="lnormal" || surv.model=="weibullph"){
survpar<-c("beta","logrho","logkappa")
if(isFALSE(all(survpar %in% surv.par.names))){
stop("use names to define 'surv.par' list; at least one
  element of 'surv.par' is missing or not recognized")
}
}
if(surv.model=="eweibull" || surv.model=="gllph" || surv.model=="ggamma"){
survpar<-c("beta","logrho","logkappa","loggamma")
if(isFALSE(all(survpar %in% surv.par.names))){
stop("use names to define 'surv.par' list; at least one
  element of 'surv.par' is missing or not recognized")
}
}  
if(length(origin)!=nrow(Data) || length(end.time)!=nrow(Data)){
stop("must be length(origin) = length(end.time) = nrow(Data)")
}
if(any(end.time<=origin)){
stop("end.time must be > origin")
}
if(length(surv.par$beta)!=length(survvar)){
stop("length of beta in 'surv.par' did not match with 'surv.formula'")
}
design.mat<-matrix(model.matrix(surv.formula,data=Data)[,-1],ncol=length(survvar))
colnames(design.mat)<-survvar
beta<-surv.par$beta
rho<-exp(surv.par$logrho)
kappa<-exp(surv.par$logkappa)
if(surv.model=="eweibull" || surv.model=="gllph" || surv.model=="ggamma"){
gamma<-exp(surv.par$loggamma)
}
if(surv.model=="weibull") {sim<-sim.waft(origin=origin,end.time=end.time,
   design.mat=design.mat,beta=beta,kappa=kappa,rho=rho)}
if(surv.model=="llogistic") {sim<-sim.llaft(origin=origin,end.time=end.time,
   design.mat=design.mat,beta=beta,kappa=kappa,rho=rho)}
if(surv.model=="lnormal") {sim<-sim.lnaft(origin=origin,end.time=end.time,
   design.mat=design.mat,beta=beta,kappa=kappa,rho=rho)}
if(surv.model=="weibullph") {sim<-sim.wph(origin=origin,end.time=end.time,
   design.mat=design.mat,beta=beta,kappa=kappa,rho=rho)}
if(surv.model=="eweibull") {sim<-sim.ewaft(origin=origin,end.time=end.time,
   design.mat=design.mat,beta=beta,kappa=kappa,gam=gamma,rho=rho)}
if(surv.model=="ggamma") {sim<-sim.ggaft(origin=origin,end.time=end.time,
   design.mat=design.mat,beta=beta,kappa=kappa,gam=gamma,rho=rho)}
if(surv.model=="gllph") {sim<-sim.gllph(origin=origin,end.time=end.time,
   design.mat=design.mat,beta=beta,kappa=kappa,gam=gamma,rho=rho)}
return(sim)
}
######################################
