jmllik.2stage<-function(spar,st,status,z,z_rand,z_fixed,m,p2,
    alpha,b,nodes1,nodesr,nodesf,weight,qpoints,surv.model,method){
beta<-spar[1:p2]
phi<-spar[p2+1]
rho<-exp(spar[p2+2])
kappa<-exp(spar[p2+3])
if(surv.model=="eweibull" || surv.model=="ggamma" || surv.model=="gllph"){
gam<-exp(spar[p2+4])
} 
pred0<-c(z%*%beta)
bbb<-b[rep(1:m, times = rep(qpoints,m)), ]
if(method=="riz"){
if(surv.model=="weibull"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
psi<-c(t(weight)%*%matrix(exp(-phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*st*exp(-pred0)/2
cum.haz<-(rho*psi)^kappa
logh <- log(kappa) + kappa * log(rho) + (kappa-1) * log(psi) - pred1
llik <- status*logh-cum.haz
} else if(surv.model=="llogistic"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
psi<-c(t(weight)%*%matrix(exp(-phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*st*exp(-pred0)/2
cum.haz<- log(1+(rho*psi)^kappa)
logh <- log(kappa) + kappa * log(rho) + (kappa-1) * log(psi) - log(1+(rho*psi)^kappa) - pred1
llik <- status*logh-cum.haz
} else if(surv.model=="lnormal"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
psi<-c(t(weight)%*%matrix(exp(-phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*st*exp(-pred0)/2
cum.haz<- - plnorm(psi, meanlog = -log(rho), sdlog = 1/kappa, lower.tail = FALSE, log.p = TRUE)
logf <- dlnorm(psi, meanlog = -log(rho), sdlog = 1/kappa, log = TRUE)-pred1
llik <- status*logf-(1-status)*cum.haz
} else if(surv.model=="eweibull"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
psi<-c(t(weight)%*%matrix(exp(-phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*st*exp(-pred0)/2
logf<-log(kappa)+log(gam)+kappa*log(rho)+(kappa-1)*log(psi)+
   (gam-1)*log1mexp(-(rho*psi)^kappa)-(rho*psi)^kappa-pred1
logs<-log1p(-exp(gam*log1mexp(-(rho*psi)^kappa)))
if(any(!is.finite(logs))){
wh<-which(!is.finite(logs))
logs[wh]<-log(gam)-(rho*psi[wh])^kappa
}
cum.haz<- -logs
llik <- status*logf-(1-status)*cum.haz
} else if(surv.model=="ggamma"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
psi<-c(t(weight)%*%matrix(exp(-phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*st*exp(-pred0)/2
logf<-log(kappa)+kappa*gam*log(rho)+(kappa*gam-1)*log(psi)-(rho*psi)^kappa-lgamma(gam)-pred1
cum.haz<--pgamma((rho*psi)^kappa,shape=gam,scale=1,lower.tail=FALSE,log.p=TRUE)
llik <- status*logf-(1-status)*cum.haz
} else if(surv.model=="weibullph"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
cum.haz<-c(t(weight)%*%matrix(exp((kappa-1)*log(nodes1)+
  phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*
  st*kappa*rho^kappa*exp(pred0)/2
logh <- log(kappa) + kappa * log(rho) + (kappa-1) * log(st) + pred1
llik <- status*logh-cum.haz
} else if(surv.model=="gllph"){
pred1<-pred0+phi*(c(z_fixed%*%alpha)+rowSums(z_rand*b))
cum.haz<-c(t(weight)%*%matrix(exp((kappa-1)*log(nodes1)-log(1+(gam*nodes1)^kappa)+
   phi*(c(nodesf%*%alpha)+rowSums(nodesr*bbb))),nrow=qpoints))*
   st*kappa*rho^kappa*exp(pred0)/2
logh <- log(kappa) + kappa * log(rho) + (kappa-1) * log(st) - log(1+(gam*st)^kappa) + pred1
llik <- status*logh-cum.haz
}
} else{
if(surv.model=="weibull"){
pred1<-pred0+phi*rowSums(z_rand*b)
psi<-c(t(weight)%*%matrix(exp(-phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*exp(-pred0)/2
cum.haz<-(rho*psi)^kappa
logh <- log(kappa) + kappa * log(rho) + (kappa-1) * log(psi) - pred1
llik <- status*logh-cum.haz
} else if(surv.model=="llogistic"){
pred1<-pred0+phi*rowSums(z_rand*b)
psi<-c(t(weight)%*%matrix(exp(-phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*exp(-pred0)/2
cum.haz<- log(1+(rho*psi)^kappa)
logh <- log(kappa) + kappa * log(rho) + (kappa-1) * log(psi) - log(1+(rho*psi)^kappa) - pred1
llik <- status*logh-cum.haz
} else if(surv.model=="lnormal"){
pred1<-pred0+phi*rowSums(z_rand*b)
psi<-c(t(weight)%*%matrix(exp(-phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*exp(-pred0)/2
cum.haz<- - plnorm(psi, meanlog = -log(rho), sdlog = 1/kappa, lower.tail = FALSE, log.p = TRUE)
logf <- dlnorm(psi, meanlog = -log(rho), sdlog = 1/kappa, log = TRUE)-pred1
llik <- status*logf-(1-status)*cum.haz
} else if(surv.model=="eweibull"){
pred1<-pred0+phi*rowSums(z_rand*b)
psi<-c(t(weight)%*%matrix(exp(-phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*exp(-pred0)/2
logf<-log(kappa)+log(gam)+kappa*log(rho)+(kappa-1)*log(psi)+
   (gam-1)*log1mexp(-(rho*psi)^kappa)-(rho*psi)^kappa-pred1
logs<-log1p(-exp(gam*log1mexp(-(rho*psi)^kappa)))
if(any(!is.finite(logs))){
wh<-which(!is.finite(logs))
logs[wh]<-log(gam)-(rho*psi[wh])^kappa
}
cum.haz<- -logs
llik <- status*logf-(1-status)*cum.haz
} else if(surv.model=="ggamma"){
pred1<-pred0+phi*rowSums(z_rand*b)
psi<-c(t(weight)%*%matrix(exp(-phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*exp(-pred0)/2
logf<-log(kappa)+kappa*gam*log(rho)+(kappa*gam-1)*log(psi)-(rho*psi)^kappa-lgamma(gam)-pred1
cum.haz<--pgamma((rho*psi)^kappa,shape=gam,scale=1,lower.tail=FALSE,log.p=TRUE)
llik <- status*logf-(1-status)*cum.haz
} else if(surv.model=="weibullph"){
pred1<-pred0+phi*rowSums(z_rand*b)
cum.haz<-c(t(weight)%*%matrix(exp((kappa-1)*log(nodes1)+
  phi*rowSums(nodesr*bbb)),nrow=qpoints))*st*kappa*rho^kappa*exp(pred0)/2
logh <- log(kappa) + kappa * log(rho) + (kappa-1) * log(st) + pred1
llik <- status*logh-cum.haz
} else if(surv.model=="gllph"){
pred1<-pred0+phi*rowSums(z_rand*b)
cum.haz<-c(t(weight)%*%matrix(exp((kappa-1)*log(nodes1)-log(1+(gam*nodes1)^kappa)+
   phi*rowSums(nodesr*bbb)),nrow=qpoints))*
   st*kappa*rho^kappa*exp(pred0)/2
logh <- log(kappa) + kappa * log(rho) + (kappa-1) * log(st) - log(1+(gam*st)^kappa) + pred1
llik <- status*logh-cum.haz
}
}
return(-sum(llik,na.rm=TRUE))
}
#------------
fit.2stage<-function(st,status,z,z_rand,z_fixed,m,p2,
    alpha,b,nodes1,nodesr,nodesf,weight,qpoints,surv.model,method){
if(method=="riz" || method=="rizopoulos"){
method<-"riz"
} else{
method<-"hen"
}
if(surv.model=="eweibull" || surv.model=="ggamma" || surv.model=="gllph"){
p<-3
} else{
p<-2
}
init<-rbind(c(rep(0,p2+1),rep(-1,p)),c(rep(0,p2+1),rep(-1.5,p)),c(rep(0,p2+1),rep(-2,p)),
   c(rep(0,p2+1),rep(0,p)),c(rep(0,p2+1),rep(-3,p)),c(rep(0,p2+1),rep(0.5,p)),
   c(rep(0,p2+1),rep(1,p)),c(rep(0,p2+1),rep(1.5,p)))
code<-1
k<-0
repeat{
k<-k+1
options(warn=-1)
fit1<-tryCatch({nlminb(start=init[k,],objective=jmllik.2stage,st=st,status=status,z=z,
  z_rand=z_rand,z_fixed=z_fixed,m=m,p2=p2,
  alpha=alpha,b=b,nodes1=nodes1,
  nodesr=nodesr,nodesf=nodesf,weight=weight,
  qpoints=qpoints,surv.model=surv.model,method=method)},
  error = function(e){NULL})
options(warn=0)
if(!is.null(fit1)){
if(fit1$iterations>1 & fit1$convergence==0 & is.finite(fit1$objective) & 
  fit1$objective!=0 & !is.na(fit1$objective)){
code<-0
}}
if(k==nrow(init) || code==0) break
}
if(code==1){
return(list(fit="non-convergence",code=1))
} else{
result1<-fit1$par
if(surv.model=="weibull" || surv.model=="llogistic" || surv.model=="lnormal" || surv.model=="weibullph"){
names(result1)<-c(colnames(z),"association","log(rho)","log(kappa)")
} else{
names(result1)<-c(colnames(z),"association","log(rho)","log(kappa)","log(gamma)")
}
results<-list(fit=result1,code=0)
class(results)<-"JMR"
attr(results, "hidden") <-c("code")
return(results)
}
}

