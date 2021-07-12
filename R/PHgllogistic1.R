###################################################
# Negative log-likelihood of GLL PH
####################################################
Rgll.llik<-function(init,tstart,tstop,status,xx){
p<-ncol(xx)
beta<-init[1:p]
kappa<-exp(init[p+1])
gam<-exp(init[p+2])
rho<-exp(init[p+3])
pred<-as.matrix(xx)%*%beta
log.ht1<-log(kappa)+kappa*log(rho)+(kappa-1)*log(tstop)-log(1+(gam*tstop)^kappa)+pred
Ht<-(rho/gam)^kappa*(log(1+(gam*tstop)^kappa)-log(1+(gam*tstart)^kappa))*exp(pred)
lf<-status*log.ht1-Ht
return(-sum(lf))
}
#######################################################
# Fit of the GLL PH.
#################################################
Rgll.fit<-function(Formula,init=NULL,data,conf.int=0.95,iter.max=250){
variables<-all.vars(Formula)
data0<-data[,variables]
if(any(is.na(data0))){
stop("missing values in the data")
}
Formula1<-as.formula(paste("~",as.character(Formula)[3],collapse=""))
xx0<-model.matrix(Formula1,data)
cnames<-colnames(xx0)[-1]
xx<-as.matrix(xx0[,-1])
colnames(xx)<-cnames
tstart<-c(as.matrix(data[variables[1]]))
tstop<-c(as.matrix(data[variables[2]]))
if(any(tstop<=tstart)){
stop("stop time must be greater than the start time")
}
status<-c(as.matrix(data[variables[3]]))
p<-ncol(xx)
if(is.null(init)){
init<-rbind(c(rep(0,p),rep(-1,3)),c(rep(0,p),rep(-1.5,3)),c(rep(0,p),rep(-2,3)),
   c(rep(0,p),rep(0,3)),c(rep(0,p),rep(-3,3)),c(rep(0,p),rep(-4,3)),
   c(rep(0,p),rep(0.5,3)),c(rep(0,p),rep(1,3)))
}
code<-1
k<-0
repeat{
k<-k+1
options(warn=-1)
fit1<-tryCatch({nlminb(start=init[k,],objective=Rgll.llik,tstart=tstart,tstop=tstop,
  status=status,xx=xx,control=list(iter.max=iter.max))},error = function(e){NULL})
options(warn=0)
if(!is.null(fit1)){
if(fit1$message=="relative convergence (4)" & fit1$convergence==0 & is.finite(fit1$objective) & 
  fit1$objective!=0 & !is.na(fit1$objective)){
options(warn=-1)
hess<-fdHess(pars=fit1$par,fun=Rgll.llik,tstart=tstart,tstop=tstop,
  status=status,xx=xx)$Hessian
options(warn=0)
if(all(is.finite(hess))){
if(!any(eigen(hess)$values<=0)){
logL<--fit1$objective
code<-0
optimizer<-"nlminb"}}
if(code==1){
options(warn=-1)
hess<-hessian(func=Rgll.llik,x=fit1$par,tstart=tstart,tstop=tstop,
  status=status,xx=xx)
options(warn=0)
if(all(is.finite(hess))){
if(!any(eigen(hess)$values<=0)){
logL<--fit1$objective
code<-0
optimizer<-"nlminb"}}}
}} 
if(k==nrow(init) || code==0) break}
if(code==1){
k<-0
repeat{
k<-k+1
options(warn=-1)
fit1<-tryCatch({optim(par=init[k,],fn=Rgll.llik,tstart=tstart,tstop=tstop,
  status=status,xx=xx,method="BFGS",hessian=TRUE,control=list(maxit=iter.max))},
  error = function(e){NULL})
options(warn=0)
if(!is.null(fit1)){
hess<-fit1$hessian
if(all(is.finite(hess)) & fit1$convergence==0 & is.finite(fit1$value) & fit1$value!=0 & !is.na(fit1$value)){
if(!any(eigen(hess)$values<=0)){
code<-0
logL<--fit1$value
optimizer<-"optim"}}}
if(k==nrow(init) || code==0) break}}
if(code==1){
return(list(code=1,message="non-convergence; try with a different set of initial values"))
} else{
est<-fit1$par
cov.mat<-solve(hess)
rownames(cov.mat)<-colnames(cov.mat)<-c(colnames(xx),"log(kappa)","log(gamma)","log(rho)")
se<-sqrt(diag(cov.mat))
z<-est/se
ci.lo<-est-se*qnorm((1-conf.int)/2,lower.tail=FALSE)
ci.up<-est+se*qnorm((1-conf.int)/2,lower.tail=FALSE)
pval<-2*pnorm(-abs(z))
result1<-data.frame(est=est,se=se,z=z,
  "p value"=ifelse(pval<0.001,"<0.001",round(pval,digits=3)),lower.95=ci.lo,upper.95=ci.up)
colnames(result1)[5]<-paste("lower",sub("^(-?)0.", "\\1.", sprintf("%.2f",conf.int)),sep="")
colnames(result1)[6]<-paste("upper",sub("^(-?)0.", "\\1.", sprintf("%.2f",conf.int)),sep="")
rownames(result1)<-c(colnames(xx),"log(kappa)","log(gamma)","log(rho)")
result2<-cbind("exp(est)"=exp(est)[1:p],"lower.95"=exp(ci.lo)[1:p],upper.95=exp(ci.up)[1:p])
rownames(result2)<-colnames(xx)
colnames(result2)<-c(paste("exp(coef)"),
   paste("lower",sub("^(-?)0.", "\\1.", sprintf("%.2f",conf.int)),sep=""),
   paste("upper",sub("^(-?)0.", "\\1.", sprintf("%.2f",conf.int)),sep=""))
result3<-cbind("exp(-est)"=exp(-est)[1:p],"lower.95"=exp(-ci.up)[1:p],upper.95=exp(-ci.lo)[1:p])
rownames(result3)<-colnames(xx)
colnames(result3)<-c(paste("exp(-coef)"),
   paste("lower",sub("^(-?)0.", "\\1.", sprintf("%.2f",conf.int)),sep=""),
   paste("upper",sub("^(-?)0.", "\\1.", sprintf("%.2f",conf.int)),sep=""))
result0<-cbind(result2,result3)
dev<--2*logL
aic<-dev+2*length(est)
result4<-cbind(logL=logL,deviance=dev,AIC=aic)
rownames(result4)<-""
dat.sum<-cbind(n=length(tstop),events=sum(status),censored=sum(1-status),predictors=p)
rownames(dat.sum)<-""
model<-noquote(rbind("Generalized Log-Logistic PH",
  "rho, kappa and gamma are parameters", 
  "Includes as special cases the log-logistic (rho = gamma) and Weibull (gamma = 0) distributions"))
colnames(model)<-""
rownames(model)<-rep("",3)
results<-list(model=model,
  "data summary"=dat.sum,fit=result1,"exp(coef) and exp(-coef)"=result0,
  "fit criteria"=result4,optimizer=optimizer,cov=cov.mat,
  st=cbind(start=tstart,stop=tstop),status=status,design.mat=xx,
  surv.model="Rgllph",code=0)
}
class(results)<-"JMR"
attr(results, "hidden") <-c("optimizer","cov","st","status","design.mat","surv.model","code")
return(results)
}
##################################################
# Residual plot for gll
#################################################
Rgll.resid<-function(fit,plot=FALSE,conf.int=0.95,xlim=NULL,ylim=NULL,
    xlab=NULL,ylab=NULL,main=NULL){
tstart<-fit$st[,1]
tstop<-fit$st[,2]
status<-fit$status
xx<-as.matrix(fit$design.mat)
p<-ncol(xx)
beta<-fit$fit[1:p,1]
kappa<-exp(fit$fit[(p+1),1])
gam<-exp(fit$fit[(p+2),1])
rho<-exp(fit$fit[(p+3),1])
pred<-as.matrix(xx)%*%beta
Rhat1<-(rho/gam)^kappa*log(1+(gam*tstart)^kappa)*exp(pred)
Rhat2<-(rho/gam)^kappa*log(1+(gam*tstop)^kappa)*exp(pred)
km.Rhat<-survfit(Surv(Rhat1,Rhat2,status)~1, type="kaplan-meier",
  conf.int=conf.int,timefix=FALSE)
if(plot){
H<-cbind(km.Rhat$time,-log(km.Rhat$surv),-log(km.Rhat$upper),-log(km.Rhat$lower))
H[which(!is.finite(H))]<-NA
H<-na.omit(H)
if(is.null(xlim) || is.null(ylim)){
xlim<-ylim<-c(min(c(H)),max(c(H)))
}
if(is.null(xlab)){
xlab<-"Cox-Snell Residuals"
}
if(is.null(ylab)){
ylab<-"Estimated Cumulative Hazard"
}
if(is.null(main)){
main<-"Cox-Snell Residual Plot"
}
plot(H[,1],H[,2],type="p",pch=20,xlab="Residuals",
  ylab="Estimated Cumulative Hazard",
  xlim=xlim,ylim=ylim,main=main,
  cex.lab=1.25,cex.main=1.5)
lines(H[,1],H[,3],type="s",lty=2)
lines(H[,1],H[,4],type="s",lty=2)
abline(a=0,b=1,col="grey75",lwd=2)
}
output<-list(residuals=as.vector(km.Rhat$time))
class(output)<-"JMR"
attr(output, "hidden") <-NULL
return(output)
}
##################################################
