jmfit.2par<-function(surv.fit,lme.fit,surv.model,
 fixed.model,rand.model,
 timevar,inits=NULL,form="rizopoulos",method="MCMC",
 warmup=1000,iter=3500,chains=2,thin=1,
 adapt_delta=0.8,max_treedepth=10,control=list()){
start_time <- Sys.time()
# Longitudinal data
id <- as.vector(unclass(lme.fit$groups[[1]]))
dat<-lme.fit$data
if(!any(colnames(dat)==timevar)){
stop("timevar is not found in the data")
}
formula.fixed <- formula(lme.fit)
model.fdata <- model.frame(terms(formula.fixed), data = dat)
inf.fdata<-attr(model.fdata, "terms")
x0<- model.matrix(formula.fixed, model.fdata)
p1<-ncol(x0)
x<-matrix(x0,ncol=p1)
formula.rand<- formula(lme.fit$modelStruct$reStruct[[1]])
model.rdata <- model.frame(terms(formula.rand), data = dat)
inf.rdata<-attr(model.rdata, "terms")
x_rand0<-model.matrix(formula.rand,model.rdata)
pp1<-ncol(x_rand0)
x_rand<-matrix(x_rand0,ncol=pp1)
y <- as.vector(model.response(model.fdata, "numeric"))
n<-nrow(x)
obs<-as.numeric(rownames(data.frame(x)))
n1<-as.numeric(tapply(obs,id,function(x) x[1]))
n2<-as.numeric(tapply(obs,id,function(x) x[(length(x))]))
id_freq<-as.vector(table(id))
# -------------------------------
# Survival data
m<-surv.fit$n
if(m!=length(unique(id))){
stop("number of subjects in the survival data 
  is not equal to the number of subjects in the 
  longitudinal data")
}
z<-matrix(surv.fit$x,nrow=m)
colnames(z)<-colnames(surv.fit$x)
p2<-ncol(z)
st0<-as.matrix(surv.fit$y)
st<-as.vector(st0[,1])
status<-as.vector(st0[,2])
# ------------------------------
# Design matrix of the link model (Survival process)
dat0 <- dat[!duplicated(id), ]
dat0[[timevar]] <-st
dat1 <- model.frame(inf.rdata, data = dat0)
dat11 <- model.frame(inf.fdata, data = dat0)
z_rand0<-model.matrix(formula.rand, dat1)
z_rand<-matrix(z_rand0,ncol=pp1)
z_fixed0<-model.matrix(formula.fixed, dat11)
z_fixed<-matrix(z_fixed0,ncol=p1)
# ---------------------------------
#-- Initial values for the longitudinal part --
alpha<-as.vector(lme.fit$coef$fixed)
b<-as.matrix(ranef(lme.fit))
colnames(b)<-NULL
sigma<-lme.fit$sigma
#-- Initial values for the survival part --
iqpoints<-31
inodes0<-GK.nodes(iqpoints)$gk.n
iweight<-GK.nodes(iqpoints)$gk.w
inodes1<-sapply(st,function(st){0.5*(st*inodes0+st)})
iid2 <- rep(1:length(st), each = iqpoints)
idat2 <- dat0[iid2, ]
idat2[[timevar]] <- c(inodes1)
idat3 <- model.frame(inf.rdata, data = idat2)
idat4 <- model.frame(inf.fdata, data = idat2)
inodesr<-matrix(model.matrix(formula.rand, idat3),ncol=pp1)
inodesf<-matrix(model.matrix(formula.fixed, idat4),ncol=p1)
ifit<-fit.2stage(st=st,status=status,z=z,z_rand=z_rand,z_fixed=z_fixed,m=m,p2=p2,
    alpha=alpha,b=b,nodes1=c(inodes1),nodesr=inodesr,nodesf=inodesf,weight=iweight,
    qpoints=iqpoints,surv.model=surv.model,method=form)
if(method=="2stage"){
lme.results<-c(fixed.effects(lme.fit),lme.fit$sigma)
names(lme.results)[length(lme.results)]<-"residual.sd"
results.2stage<-list("survival sub-model"=ifit$fit,"longitudinal sub-model"=lme.results)
class(results.2stage)<-"JMR"
attr(results.2stage, "hidden") <-NULL
return(results.2stage)
}
if(ifit$code==0){
beta<-beta_mu<-ifit$fit[1:p2]
phi<-phi_mu<-ifit$fit[p2+1]
beta_sd<-rep(1,p2)
phi_sd<-1
rho<-exp(ifit$fit[p2+2])
kappa<-exp(ifit$fit[p2+3])
} else{
surv.dat.frame<-data.frame(st=st,status=status,z)
cformula<-paste(colnames(surv.dat.frame)[-c(1,2)],collapse= "+")
Formula<-as.formula(paste("Surv(st,status)~",cformula,collapse=""))
#-- Weibull --
if(surv.model=="weibull"){
fit<-waft.fit(Formula,init=NULL,data=surv.dat.frame)
if(fit$code==0){
kappa<-exp(fit[[3]][(p2+1),1])
rho<-exp(fit[[3]][(p2+2),1])
beta<-c(fit[[3]][1:p2,1])
phi<-0.1
beta_mu<-rep(0,p2)
beta_sd<-rep(10,p2)
phi_mu<-0
phi_sd<-10
} 
if(fit$code!=0){
kappa<-1
rho<-1
beta<-rep(0,p2)
phi<-0.1
beta_mu<-rep(0,p2)
beta_sd<-rep(10,p2)
phi_mu<-0
phi_sd<-10
}
}
#-- Log-logistic --
if(surv.model=="llogistic"){
fit<-llaft.fit(Formula,init=NULL,data=surv.dat.frame)
if(fit$code==0){
kappa<-exp(fit[[3]][(p2+1),1])
rho<-exp(fit[[3]][(p2+2),1])
beta<-c(fit[[3]][1:p2,1])
phi<-0.1
beta_mu<-rep(0,p2)
beta_sd<-rep(10,p2)
phi_mu<-0
phi_sd<-10
} 
if(fit$code!=0){
kappa<-1
rho<-1
beta<-rep(0,p2)
phi<-0.1
beta_mu<-rep(0,p2)
beta_sd<-rep(10,p2)
phi_mu<-0
phi_sd<-10
}
}
#-- Log-normal --
if(surv.model=="lnormal"){
fit<-lnaft.fit(Formula,init=NULL,data=surv.dat.frame)
if(fit$code==0){
kappa<-exp(fit[[3]][(p2+1),1])
rho<-exp(fit[[3]][(p2+2),1])
beta<-c(fit[[3]][1:p2,1])
phi<-0.1
beta_mu<-rep(0,p2)
beta_sd<-rep(10,p2)
phi_mu<-0
phi_sd<-10
} 
if(fit$code!=0){
kappa<-1
rho<-1
beta<-rep(0,p2)
phi<-0.1
beta_mu<-rep(0,p2)
beta_sd<-rep(10,p2)
phi_mu<-0
phi_sd<-10
}
}
#-- Weibull PH --
if(surv.model=="weibullph"){
fit<-waft.fit(Formula,init=NULL,data=surv.dat.frame)
if(fit$code==0){
kappa<-exp(fit[[3]][(p2+1),1])
rho<-exp(fit[[3]][(p2+2),1])
beta<-c(-fit[[3]][1:p2,1]*exp(fit[[3]][(p2+1),1]))
phi<-0.1
beta_mu<-rep(0,p2)
beta_sd<-rep(10,p2)
phi_mu<-0
phi_sd<-10
} 
if(fit$code!=0){
kappa<-1
rho<-1
beta<-rep(0,p2)
phi<-0.1
beta_mu<-rep(0,p2)
beta_sd<-rep(10,p2)
phi_mu<-0
phi_sd<-10
}
}
}
# ------------------------------------
con<-list(alpha=as.array(c(alpha)),alpha_mu=as.array(c(alpha)),alpha_sd=as.array(c(rep(1,p1))),
    beta=as.array(c(beta)),beta_mu=as.array(c(beta_mu)),beta_sd=as.array(c(beta_sd)),
    phi=phi,phi_mu=phi_mu,phi_sd=phi_sd,
    b=b,rand_cov=matrix(getVarCov(lme.fit),ncol=pp1),nu=pp1+1,A=diag(pp1)*0.1,
    sigma=sigma,a0=0.1,a1=0.1,
    rho=rho,b0=1,b1=0.0005,
    kappa=kappa,c0=0.01,c1=0.01,
    qpoints=5,seed=sample.int(.Machine$integer.max, 1))
nmsC <- names(con)
con[(namc <- names(control))] <- control
if (length(noNms <- namc[!namc %in% nmsC])) 
  stop("unknown names in control: ", paste(noNms, collapse = ", "))
if(any(con$alpha_sd<=0))
stop("\nalpha_sd must be positive.")
if(any(con$beta_sd<=0))
stop("\nbeta_sd must be positive.")
if(con$phi_sd<=0)
stop("\nphi_sd must be positive.")
if(con$nu<pp1)
stop("\n inv-Wishart df must be >= ", paste(pp1, collapse = ", "))
if(any(eigen(con$A)$values<=0))
stop("\nA must be a symmetric positive definite matrix.")
if(!isSymmetric(con$A))
stop("\nA must be a symmetric positive definite matrix.")
if(con$sigma<0)
stop("\ninitial value for sigma must be positive.")
if(con$kappa<=0)
stop("\nkappa must be positive.")
if(con$rho<=0)
stop("\nrho must be positive.")
if(any(eigen(con$rand_cov)$values<=0))
stop("\nrand_cov must be a symmetric positive definite matrix.")
if(!isSymmetric(con$rand_cov))
stop("\nrand_cov must be a symmetric positive definite matrix.")
if(!(con$qpoints==3 || con$qpoints==5 || con$qpoints==7 || con$qpoints==15 || 
   con$qpoints==21 || con$qpoints==31 || con$qpoints==41))
stop("\nquad.points must be 3, 5, 7, 15, 21, 31 or 41.")
if(con$a0<=0 || con$a1<=0)
stop("\na0 and a1 must be positive.")
if(con$b0<=0 || con$b1<=0)
stop("\nb0 and b1 must be positive.")
if(con$c0<=0 || con$c1<=0)
stop("\nc0 and c1 must be positive.")
# ---------------------------------
# Data for numerical integration
qpoints<-con$qpoints
nodes0<-GK.nodes(qpoints)$gk.n
weight<-GK.nodes(qpoints)$gk.w
nodes1<-sapply(st,function(st){0.5*(st*nodes0+st)})
id2 <- rep(1:length(st), each = qpoints)
dat2 <- dat0[id2, ]
dat2[[timevar]] <- c(nodes1)
dat3 <- model.frame(inf.rdata, data = dat2)
dat4 <- model.frame(inf.fdata, data = dat2)
nodesr<-matrix(model.matrix(formula.rand, dat3),ncol=pp1)
nodesf<-matrix(model.matrix(formula.fixed, dat4),ncol=p1)
nodesr.array<-list()
nodesf.array<-list()
for(i in 1:m){
nodesr.array[[i]]<-nodesr[(1+(i-1)*qpoints):(i*qpoints),]
nodesf.array[[i]]<-nodesf[(1+(i-1)*qpoints):(i*qpoints),]
}
idd<-rep(1:m,times=as.vector(table(id)))
# ---------------------------------
# Initial values for stan, and parameters to be monitored
if(is.null(inits)){
inits<-rep(list(list(alpha=con[["alpha"]],beta=con[["beta"]],phi=con[["phi"]],
   squ_kappa=con[["kappa"]]^2,inv_rho2=1/con[["rho"]]^2,
   b_unscaled=solve(t(chol(getVarCov(lme.fit))))%*%t(con[["b"]]),
   rand_cov=con[["rand_cov"]],inv_sigma2=1/con[["sigma"]]^2)),chains)
}
par.stan<-c("beta","phi","rho","kappa","alpha","rand_cov","sigma","b")
# ---------------------------------
# Stan Model
if(form=="riz" || form=="rizopoulos"){
if(surv.model=="weibull"){
if(fixed.model=="simple" && rand.model=="simple"){
x_fixed0<-model.matrix(formula.fixed, dat[!duplicated(id), ])
which_alphatime<-which(colnames(x_fixed0)==timevar)
which_alphafixed<-(1:p1)[-which_alphatime]
x_fixed<-x_fixed0[,-which_alphatime]
if(is.vector(x_fixed)){
x_fixed<-matrix(x_fixed,ncol=1)
}
data.stan<-list(m=m,n=n,p1=p1,p2=p2,pp1=pp1,
  id_freq=id_freq,
  y=y,st=st,status=status,
  z=z,z_rand=z_rand,z_fixed=z_fixed,
  x=x,x_fixed=x_fixed,x_rand=x_rand,
  which_alphatime=which_alphatime,which_alphafixed=as.array(c(which_alphafixed)),
  nu=con[["nu"]],A=con[["A"]],
  alpha_mu=con[["alpha_mu"]],alpha_sd=con[["alpha_sd"]],
  beta_mu=con[["beta_mu"]],beta_sd=con[["beta_sd"]],
  phi_mu=con[["phi_mu"]],phi_sd=con[["phi_sd"]],
  a0=con[["a0"]],a1=con[["a1"]],b0=con[["b0"]],b1=con[["b1"]],
  c0=con[["c0"]],c1=con[["c1"]])
#stanc("stan_weibullrizsimple.stan")
#stan_model <- "stan_weibullrizsimple.stan"
stan_model<-stanmodels$weibullrizsimple
} else{
data.stan<-list(m=m,n=n,p1=p1,p2=p2,pp1=pp1,
    qpoints=con[["qpoints"]],id_freq=id_freq,
    y=y,st=st,status=status,
    weights=weight,nodesf=nodesf,nodesr=nodesr,
    z=z,z_rand=z_rand,z_fixed=z_fixed,
    x=x,x_rand=x_rand,
    nu=con[["nu"]],A=con[["A"]],
    alpha_mu=con[["alpha_mu"]],alpha_sd=con[["alpha_sd"]],
    beta_mu=con[["beta_mu"]],beta_sd=con[["beta_sd"]],
    phi_mu=con[["phi_mu"]],phi_sd=con[["phi_sd"]],
    a0=con[["a0"]],a1=con[["a1"]],b0=con[["b0"]],b1=con[["b1"]],
    c0=con[["c0"]],c1=con[["c1"]])
#stanc("stan_weibullriz.stan")
#stan_model <- "stan_weibullriz.stan"
stan_model<-stanmodels$weibullriz
}
method0<-rbind("Induced by the longitudinal value.",
  "Rizopoulos D (2012). Joint Models for Longitudinal and Time-to-Event", 
  "       Data With Applications in R, Chapman and Hall/CRC.",
  "Tseng YK, Hsieh F, and Wang JL (2005). Joint Modelling of Accelerated Failure",
  "       Time and Longitudinal Data, Biometrika 92: 587 - 603.")
colnames(method0)<-""
rownames(method0)<-c("","","","","")
method0<-noquote(method0)
}
if(surv.model=="llogistic"){
if(fixed.model=="simple" && rand.model=="simple"){
x_fixed0<-model.matrix(formula.fixed, dat[!duplicated(id), ])
which_alphatime<-which(colnames(x_fixed0)==timevar)
which_alphafixed<-(1:p1)[-which_alphatime]
x_fixed<-x_fixed0[,-which_alphatime]
if(is.vector(x_fixed)){
x_fixed<-matrix(x_fixed,ncol=1)
}
data.stan<-list(m=m,n=n,p1=p1,p2=p2,pp1=pp1,
  id_freq=id_freq,
  y=y,st=st,status=status,
  z=z,z_rand=z_rand,z_fixed=z_fixed,
  x=x,x_fixed=x_fixed,x_rand=x_rand,
  which_alphatime=which_alphatime,which_alphafixed=as.array(c(which_alphafixed)),
  nu=con[["nu"]],A=con[["A"]],
  alpha_mu=con[["alpha_mu"]],alpha_sd=con[["alpha_sd"]],
  beta_mu=con[["beta_mu"]],beta_sd=con[["beta_sd"]],
  phi_mu=con[["phi_mu"]],phi_sd=con[["phi_sd"]],
  a0=con[["a0"]],a1=con[["a1"]],b0=con[["b0"]],b1=con[["b1"]],
  c0=con[["c0"]],c1=con[["c1"]])
#stanc("stan_llogisticrizsimple.stan")
#stan_model <- "stan_llogisticrizsimple.stan"
stan_model<-stanmodels$llogisticrizsimple
} else{
data.stan<-list(m=m,n=n,p1=p1,p2=p2,pp1=pp1,
    qpoints=con[["qpoints"]],id_freq=id_freq,
    y=y,st=st,status=status,
    weights=weight,nodesf=nodesf,nodesr=nodesr,
    z=z,z_rand=z_rand,z_fixed=z_fixed,
    x=x,x_rand=x_rand,
    nu=con[["nu"]],A=con[["A"]],
    alpha_mu=con[["alpha_mu"]],alpha_sd=con[["alpha_sd"]],
    beta_mu=con[["beta_mu"]],beta_sd=con[["beta_sd"]],
    phi_mu=con[["phi_mu"]],phi_sd=con[["phi_sd"]],
    a0=con[["a0"]],a1=con[["a1"]],b0=con[["b0"]],b1=con[["b1"]],
    c0=con[["c0"]],c1=con[["c1"]])
#stanc("stan_llogisticriz.stan")
#stan_model <- "stan_llogisticriz.stan"
stan_model<-stanmodels$llogisticriz
}
method0<-rbind("Induced by the longitudinal value.",
  "Rizopoulos D (2012). Joint Models for Longitudinal and Time-to-Event", 
  "       Data With Applications in R, Chapman and Hall/CRC.",
  "Tseng YK, Hsieh F, and Wang JL (2005). Joint Modelling of Accelerated Failure",
  "       Time and Longitudinal Data, Biometrika 92: 587 - 603.")
colnames(method0)<-""
rownames(method0)<-c("","","","","")
method0<-noquote(method0)
}
if(surv.model=="lnormal"){
if(fixed.model=="simple" && rand.model=="simple"){
x_fixed0<-model.matrix(formula.fixed, dat[!duplicated(id), ])
which_alphatime<-which(colnames(x_fixed0)==timevar)
which_alphafixed<-(1:p1)[-which_alphatime]
x_fixed<-x_fixed0[,-which_alphatime]
if(is.vector(x_fixed)){
x_fixed<-matrix(x_fixed,ncol=1)
}
data.stan<-list(m=m,n=n,p1=p1,p2=p2,pp1=pp1,
  id_freq=id_freq,
  y=y,st=st,status=status,
  z=z,z_rand=z_rand,z_fixed=z_fixed,
  x=x,x_fixed=x_fixed,x_rand=x_rand,
  which_alphatime=which_alphatime,which_alphafixed=as.array(c(which_alphafixed)),
  nu=con[["nu"]],A=con[["A"]],
  alpha_mu=con[["alpha_mu"]],alpha_sd=con[["alpha_sd"]],
  beta_mu=con[["beta_mu"]],beta_sd=con[["beta_sd"]],
  phi_mu=con[["phi_mu"]],phi_sd=con[["phi_sd"]],
  a0=con[["a0"]],a1=con[["a1"]],b0=con[["b0"]],b1=con[["b1"]],
  c0=con[["c0"]],c1=con[["c1"]])
#stanc("stan_lnormalrizsimple.stan")
#stan_model <- "stan_lnormalrizsimple.stan"
stan_model<-stanmodels$lnormalrizsimple
} else{
data.stan<-list(m=m,n=n,p1=p1,p2=p2,pp1=pp1,
    qpoints=con[["qpoints"]],id_freq=id_freq,
    y=y,st=st,status=status,
    weights=weight,nodesf=nodesf,nodesr=nodesr,
    z=z,z_rand=z_rand,z_fixed=z_fixed,
    x=x,x_rand=x_rand,
    nu=con[["nu"]],A=con[["A"]],
    alpha_mu=con[["alpha_mu"]],alpha_sd=con[["alpha_sd"]],
    beta_mu=con[["beta_mu"]],beta_sd=con[["beta_sd"]],
    phi_mu=con[["phi_mu"]],phi_sd=con[["phi_sd"]],
    a0=con[["a0"]],a1=con[["a1"]],b0=con[["b0"]],b1=con[["b1"]],
    c0=con[["c0"]],c1=con[["c1"]])
#stanc("stan_lnormalriz.stan")
#stan_model <- "stan_lnormalriz.stan"
stan_model<-stanmodels$lnormalriz
}
method0<-rbind("Induced by the longitudinal value.",
  "Rizopoulos D (2012). Joint Models for Longitudinal and Time-to-Event", 
  "       Data With Applications in R, Chapman and Hall/CRC.",
  "Tseng YK, Hsieh F, and Wang JL (2005). Joint Modelling of Accelerated Failure",
  "       Time and Longitudinal Data, Biometrika 92: 587 - 603.")
colnames(method0)<-""
rownames(method0)<-c("","","","","")
method0<-noquote(method0)
}
if(surv.model=="weibullph"){
data.stan<-list(m=m,n=n,p1=p1,p2=p2,pp1=pp1,
    qpoints=con[["qpoints"]],id_freq=id_freq,
    y=y,st=st,status=status,
    weights=weight,nodesf=nodesf,nodesr=nodesr,nodes1=c(nodes1),
    z=z,z_rand=z_rand,z_fixed=z_fixed,
    x=x,x_rand=x_rand,
    nu=con[["nu"]],A=con[["A"]],
    alpha_mu=con[["alpha_mu"]],alpha_sd=con[["alpha_sd"]],
    beta_mu=con[["beta_mu"]],beta_sd=con[["beta_sd"]],
    phi_mu=con[["phi_mu"]],phi_sd=con[["phi_sd"]],
    a0=con[["a0"]],a1=con[["a1"]],b0=con[["b0"]],b1=con[["b1"]],
    c0=con[["c0"]],c1=con[["c1"]])
#stanc("stan_weibullphriz.stan")
#stan_model <- "stan_weibullphriz.stan"
stan_model<-stanmodels$weibullphriz
method0<-rbind("Induced by the longitudinal value.",
  "Rizopoulos D (2012). Joint Models for Longitudinal and Time-to-Event", 
  "       Data With Applications in R, Chapman and Hall/CRC.")
colnames(method0)<-""
rownames(method0)<-c("","","")
method0<-noquote(method0)
}
} else{
if(surv.model=="weibull"){
if(rand.model=="simple"){
data.stan<-list(m=m,n=n,p1=p1,p2=p2,pp1=pp1,
  id_freq=id_freq,
  y=y,st=st,status=status,
  z=z,z_rand=z_rand,
  x=x,x_rand=x_rand,
  nu=con[["nu"]],A=con[["A"]],
  alpha_mu=con[["alpha_mu"]],alpha_sd=con[["alpha_sd"]],
  beta_mu=con[["beta_mu"]],beta_sd=con[["beta_sd"]],
  phi_mu=con[["phi_mu"]],phi_sd=con[["phi_sd"]],
  a0=con[["a0"]],a1=con[["a1"]],b0=con[["b0"]],b1=con[["b1"]],
  c0=con[["c0"]],c1=con[["c1"]])
#stanc("stan_weibullsimple.stan")
#stan_model <- "stan_weibullsimple.stan"
stan_model<-stanmodels$weibullsimple
} else{
data.stan<-list(m=m,n=n,p1=p1,p2=p2,pp1=pp1,
    qpoints=con[["qpoints"]],id_freq=id_freq,
    y=y,st=st,status=status,
    weights=weight,nodesr=nodesr,
    z=z,z_rand=z_rand,
    x=x,x_rand=x_rand,
    nu=con[["nu"]],A=con[["A"]],
    alpha_mu=con[["alpha_mu"]],alpha_sd=con[["alpha_sd"]],
    beta_mu=con[["beta_mu"]],beta_sd=con[["beta_sd"]],
    phi_mu=con[["phi_mu"]],phi_sd=con[["phi_sd"]],
    a0=con[["a0"]],a1=con[["a1"]],b0=con[["b0"]],b1=con[["b1"]],
    c0=con[["c0"]],c1=con[["c1"]])
#stanc("stan_weibull.stan")
#stan_model <- "stan_weibull.stan"
stan_model<-stanmodels$weibull
}
method0<-rbind("Induced by the longitudinal value.",
  "Guo X and Carlin BP (2004). Separate and Joint Modeling of Longitudinal", 
  "       and Event Time Data Using Standard Computer Packages,",
  "       The American Statistician, 58: 1 - 9.",
  "Henderson R, Diggle P, and Dobson A (2000). Joint Modelling of Longitudinal", 
  "       Measurements and Event Time Data, Biostatistics, 1: 465 - 480.",
  "Tseng YK, Hsieh F, and Wang JL (2005). Joint Modelling of Accelerated Failure",
  "       Time and Longitudinal Data, Biometrika 92: 587 - 603.")
colnames(method0)<-""
rownames(method0)<-c("","","","","","","","")
method0<-noquote(method0)
}
if(surv.model=="llogistic"){
if(rand.model=="simple"){
data.stan<-list(m=m,n=n,p1=p1,p2=p2,pp1=pp1,
  id_freq=id_freq,
  y=y,st=st,status=status,
  z=z,z_rand=z_rand,
  x=x,x_rand=x_rand,
  nu=con[["nu"]],A=con[["A"]],
  alpha_mu=con[["alpha_mu"]],alpha_sd=con[["alpha_sd"]],
  beta_mu=con[["beta_mu"]],beta_sd=con[["beta_sd"]],
  phi_mu=con[["phi_mu"]],phi_sd=con[["phi_sd"]],
  a0=con[["a0"]],a1=con[["a1"]],b0=con[["b0"]],b1=con[["b1"]],
  c0=con[["c0"]],c1=con[["c1"]])
#stanc("stan_llogisticsimple.stan")
#stan_model <- "stan_llogisticsimple.stan"
stan_model<-stanmodels$llogisticsimple
} else{
data.stan<-list(m=m,n=n,p1=p1,p2=p2,pp1=pp1,
    qpoints=con[["qpoints"]],id_freq=id_freq,
    y=y,st=st,status=status,
    weights=weight,nodesr=nodesr,
    z=z,z_rand=z_rand,
    x=x,x_rand=x_rand,
    nu=con[["nu"]],A=con[["A"]],
    alpha_mu=con[["alpha_mu"]],alpha_sd=con[["alpha_sd"]],
    beta_mu=con[["beta_mu"]],beta_sd=con[["beta_sd"]],
    phi_mu=con[["phi_mu"]],phi_sd=con[["phi_sd"]],
    a0=con[["a0"]],a1=con[["a1"]],b0=con[["b0"]],b1=con[["b1"]],
    c0=con[["c0"]],c1=con[["c1"]])
#stanc("stan_llogistic.stan")
#stan_model <- "stan_llogistic.stan"
stan_model<-stanmodels$llogistic
}
method0<-rbind("Induced by the longitudinal value.",
  "Guo X and Carlin BP (2004). Separate and Joint Modeling of Longitudinal", 
  "       and Event Time Data Using Standard Computer Packages,",
  "       The American Statistician, 58: 1 - 9.",
  "Henderson R, Diggle P, and Dobson A (2000). Joint Modelling of Longitudinal", 
  "       Measurements and Event Time Data, Biostatistics, 1: 465 - 480.",
  "Tseng YK, Hsieh F, and Wang JL (2005). Joint Modelling of Accelerated Failure",
  "       Time and Longitudinal Data, Biometrika 92: 587 - 603.")
colnames(method0)<-""
rownames(method0)<-c("","","","","","","","")
method0<-noquote(method0)
}
if(surv.model=="lnormal"){
if(rand.model=="simple"){
data.stan<-list(m=m,n=n,p1=p1,p2=p2,pp1=pp1,
  id_freq=id_freq,
  y=y,st=st,status=status,
  z=z,z_rand=z_rand,
  x=x,x_rand=x_rand,
  nu=con[["nu"]],A=con[["A"]],
  alpha_mu=con[["alpha_mu"]],alpha_sd=con[["alpha_sd"]],
  beta_mu=con[["beta_mu"]],beta_sd=con[["beta_sd"]],
  phi_mu=con[["phi_mu"]],phi_sd=con[["phi_sd"]],
  a0=con[["a0"]],a1=con[["a1"]],b0=con[["b0"]],b1=con[["b1"]],
  c0=con[["c0"]],c1=con[["c1"]])
#stanc("stan_lnormalsimple.stan")
#stan_model <- "stan_lnormalsimple.stan"
stan_model<-stanmodels$lnormalsimple
} else{
data.stan<-list(m=m,n=n,p1=p1,p2=p2,pp1=pp1,
    qpoints=con[["qpoints"]],id_freq=id_freq,
    y=y,st=st,status=status,
    weights=weight,nodesr=nodesr,
    z=z,z_rand=z_rand,
    x=x,x_rand=x_rand,
    nu=con[["nu"]],A=con[["A"]],
    alpha_mu=con[["alpha_mu"]],alpha_sd=con[["alpha_sd"]],
    beta_mu=con[["beta_mu"]],beta_sd=con[["beta_sd"]],
    phi_mu=con[["phi_mu"]],phi_sd=con[["phi_sd"]],
    a0=con[["a0"]],a1=con[["a1"]],b0=con[["b0"]],b1=con[["b1"]],
    c0=con[["c0"]],c1=con[["c1"]])
#stanc("stan_lnormal.stan")
#stan_model <- "stan_lnormal.stan"
stan_model<-stanmodels$lnormal
}
method0<-rbind("Induced by the longitudinal value.",
  "Guo X and Carlin BP (2004). Separate and Joint Modeling of Longitudinal", 
  "       and Event Time Data Using Standard Computer Packages,",
  "       The American Statistician, 58: 1 - 9.",
  "Henderson R, Diggle P, and Dobson A (2000). Joint Modelling of Longitudinal", 
  "       Measurements and Event Time Data, Biostatistics, 1: 465 - 480.",
  "Tseng YK, Hsieh F, and Wang JL (2005). Joint Modelling of Accelerated Failure",
  "       Time and Longitudinal Data, Biometrika 92: 587 - 603.")
colnames(method0)<-""
rownames(method0)<-c("","","","","","","","")
method0<-noquote(method0)
}
if(surv.model=="weibullph"){
data.stan<-list(m=m,n=n,p1=p1,p2=p2,pp1=pp1,
    qpoints=con[["qpoints"]],id_freq=id_freq,
    y=y,st=st,status=status,
    weights=weight,nodesr=nodesr,nodes1=c(nodes1),
    z=z,z_rand=z_rand,
    x=x,x_rand=x_rand,
    nu=con[["nu"]],A=con[["A"]],
    alpha_mu=con[["alpha_mu"]],alpha_sd=con[["alpha_sd"]],
    beta_mu=con[["beta_mu"]],beta_sd=con[["beta_sd"]],
    phi_mu=con[["phi_mu"]],phi_sd=con[["phi_sd"]],
    a0=con[["a0"]],a1=con[["a1"]],b0=con[["b0"]],b1=con[["b1"]],
    c0=con[["c0"]],c1=con[["c1"]])
#stanc("stan_weibullph.stan")
#stan_model <- "stan_weibullph.stan"
stan_model<-stanmodels$weibullph
method0<-rbind("Induced by the longitudinal value.",
  "Guo X and Carlin BP (2004). Separate and Joint Modeling of Longitudinal", 
  "       and Event Time Data Using Standard Computer Packages,",
  "       The American Statistician, 58: 1 - 9.",
  "Henderson R, Diggle P, and Dobson A (2000). Joint Modelling of Longitudinal", 
  "       Measurements and Event Time Data, Biostatistics, 1: 465 - 480.")
colnames(method0)<-""
rownames(method0)<-c("","","","","","")
method0<-noquote(method0)
}
}
#--------------------------------
# Run stan
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7 -mtune=corei7')
stan_fit <- rstan::sampling(object=stan_model, data = data.stan, 
  init = inits, pars = par.stan, warmup = warmup, iter = iter, 
  chains = chains, thin = thin,seed=con$seed, 
  control=list(adapt_delta=adapt_delta,max_treedepth=max_treedepth))
if(!is.array(stan_fit) || length(dim(stan_fit))==0){
for(i in 1: length(inits)){
inits[[i]]$squ_kappa<-1
inits[[i]]$inv_rho2<-1
inits[[i]]$beta<-rep(0,p2)
}
stan_fit <- rstan::sampling(object=stan_model, data = data.stan, 
  init = inits, pars = par.stan, warmup = warmup, iter = iter, 
  chains = chains, thin = thin,seed=con$seed, 
  control=list(adapt_delta=adapt_delta,max_treedepth=max_treedepth))
}
#------------------------------------------------------
# Summarize
if(is.array(stan_fit) && length(dim(stan_fit))==3){
if(dim(stan_fit)[1]==iter-warmup && dim(stan_fit)[2]==chains && 
   dim(stan_fit)[3]==length(unlist(inits[[1]]))+1){
sum.stanfit<-rstan::summary(stan_fit,pars=c("beta","phi","rho","kappa","alpha","rand_cov","sigma"))$summary
sum.stanfit[,"Rhat"]<-round(sum.stanfit[,"Rhat"],digits=2)
sum.b<-rstan::summary(stan_fit,pars=c("b"))$summary
beta.loc<-1:p2
phi.loc<-p2+1
rho.loc<-p2+2
kappa.loc<-p2+3
alpha.loc<-(p2+4):(p2+3+p1)
cov.loc<-(p2+4+p1):(p2+3+p1+pp1^2)
sigma.loc<-p2+4+p1+pp1^2
b.loc<-(p2+5+p1+pp1^2):(p2+4+p1+pp1^2+pp1*m)
#---------------------------------------------------
# Save data
if(fixed.model=="simple" && rand.model=="simple" && (form=="riz" || form=="rizopoulos")){
stan.data<-list(m=m,n=n,p1=p1,p2=p2,pp1=pp1,id=idd,
    qpoints=con[["qpoints"]],id_freq=id_freq,
    y=y,st=st,status=status,
    weights=weight,nodesf=nodesf,nodesr=nodesr,nodes1=nodes1,
    z=z,z_rand=z_rand,z_fixed=z_fixed,
    x=x,x_rand=x_rand,x_fixed=x_fixed,
    which_alphatime=which_alphatime,which_alphafixed=as.array(c(which_alphafixed)),
    nu=con[["nu"]],A=con[["A"]],
    alpha=con[["alpha"]],alpha_mu=con[["alpha_mu"]],alpha_sd=con[["alpha_sd"]],
    beta=con[["beta"]],beta_mu=con[["beta_mu"]],beta_sd=con[["beta_sd"]],
    phi=con[["phi"]],phi_mu=con[["phi_mu"]],phi_sd=con[["phi_sd"]],
    b=con[["b"]],rand_cov=con[["rand_cov"]],nu=con[["nu"]],A=con[["A"]],
    sigma=con[["sigma"]],a0=con[["a0"]],a1=con[["a1"]],
    rho=con[["rho"]],b0=con[["b0"]],b1=con[["b1"]],
    kappa=con[["kappa"]],c0=con[["c0"]],c1=con[["c1"]])
} else{
stan.data<-list(m=m,n=n,p1=p1,p2=p2,pp1=pp1,id=idd,
    qpoints=con[["qpoints"]],id_freq=id_freq,
    y=y,st=st,status=status,
    weights=weight,nodesf=nodesf,nodesr=nodesr,nodes1=nodes1,
    z=z,z_rand=z_rand,z_fixed=z_fixed,
    x=x,x_rand=x_rand,
    nu=con[["nu"]],A=con[["A"]],
    alpha=con[["alpha"]],alpha_mu=con[["alpha_mu"]],alpha_sd=con[["alpha_sd"]],
    beta=con[["beta"]],beta_mu=con[["beta_mu"]],beta_sd=con[["beta_sd"]],
    phi=con[["phi"]],phi_mu=con[["phi_mu"]],phi_sd=con[["phi_sd"]],
    b=con[["b"]],rand_cov=con[["rand_cov"]],nu=con[["nu"]],A=con[["A"]],
    sigma=con[["sigma"]],a0=con[["a0"]],a1=con[["a1"]],
    rho=con[["rho"]],b0=con[["b0"]],b1=con[["b1"]],
    kappa=con[["kappa"]],c0=con[["c0"]],c1=con[["c1"]])
}
post.mean<-list(beta=sum.stanfit[beta.loc,1],phi=sum.stanfit[phi.loc,1],
  kappa=sum.stanfit[kappa.loc,1],rho=sum.stanfit[rho.loc,1],
  alpha=sum.stanfit[alpha.loc,1],cov_rand=matrix(sum.stanfit[cov.loc,1],ncol=pp1),
  sigma=sum.stanfit[sigma.loc,1],b=matrix(sum.b[,1],byrow=TRUE,ncol=pp1))
post.med<-list(beta=sum.stanfit[beta.loc,6],phi=sum.stanfit[phi.loc,6],
  kappa=sum.stanfit[kappa.loc,6],rho=sum.stanfit[rho.loc,6],
  alpha=sum.stanfit[alpha.loc,6],cov_rand=matrix(sum.stanfit[cov.loc,6],ncol=pp1),
  sigma=sum.stanfit[sigma.loc,6],b=matrix(sum.b[,6],byrow=TRUE,ncol=pp1))
surv.data<-cbind(id=unique(id),st=st,status=status,z)
colnames(surv.data)<-c("id","st","status",colnames(surv.fit$x))
long.data<-cbind(id,y,x)
colnames(long.data)<-c("id","y","Intercept",colnames(x0)[-1])
results<-list()
results[["stan.fit"]]<-stan_fit
results[["surv.data"]]<-surv.data
results[["long.data"]]<-long.data
results[["lme.fit"]]<-lme.fit
results[["surv.fit"]]<-surv.fit
results[["stan.data"]]<-stan.data
results[["post.mean"]]<-post.mean
results[["post.med"]]<-post.med
results[["beta.loc"]]<-beta.loc
results[["phi.loc"]]<-phi.loc
results[["kappa.loc"]]<-kappa.loc
results[["rho.loc"]]<-rho.loc
results[["alpha.loc"]]<-alpha.loc
results[["cov.loc"]]<-cov.loc
results[["sigma.loc"]]<-sigma.loc
results[["b.loc"]]<-b.loc
results[["nodes0"]]<-nodes0
results[["timevar"]]<-timevar
results[["fixed.model"]]<-fixed.model
results[["rand.model"]]<-rand.model
results[["ifit"]]<-ifit
results[["est.method"]]<-method
if(form=="riz" || form=="rizopoulos"){
results[["method"]]<-"riz"
} else{
results[["method"]]<-"hen"
}
l.model<-rbind("LME","alpha = fixed-effects coefficients", "b = random-effects coefficients",
  "rand_cov = covariance matrix for b","sigma = sd of random errors")
colnames(l.model)<-""
rownames(l.model)<-c("","","","","")
l.model<-noquote(l.model)
if(surv.model=="weibull"){
results[["surv.model"]]<-"weibull"
s.model<-rbind("Weibull AFT", "rho = rate parameter","kappa = shape parameter",
  "beta = regression coefficient vector",
  "phi = association parameter")
colnames(s.model)<-""
rownames(s.model)<-c("","","","","")
s.model<-noquote(s.model)
}
if(surv.model=="llogistic"){
results[["surv.model"]]<-"llogistic"
s.model<-rbind("Log-logistic AFT", "rho = rate parameter","kappa = shape parameter",
  "beta = regression coefficient vector",
  "phi = association parameter")
colnames(s.model)<-""
rownames(s.model)<-c("","","","","")
s.model<-noquote(s.model)
}
if(surv.model=="lnormal"){
results[["surv.model"]]<-"lnormal"
s.model<-rbind("Log-normal AFT", "rho = rate parameter","kappa = shape parameter",
  "beta = regression coefficient vector",
  "phi = association parameter")
colnames(s.model)<-""
rownames(s.model)<-c("","","","","")
s.model<-noquote(s.model)
}
if(surv.model=="weibullph"){
results[["surv.model"]]<-"weibullph"
s.model<-rbind("Weibull PH", "rho = rate parameter","kappa = shape parameter",
  "beta = regression coefficient vector",
  "phi = association parameter")
colnames(s.model)<-""
rownames(s.model)<-c("","","","","")
s.model<-noquote(s.model)
}
results[["method0"]]<-method0
results[["s.model"]]<-s.model
results[["l.model"]]<-l.model
end_time <- Sys.time()
if(qpoints==5 || qpoints==7){
nint<-noquote(paste("Gauss Legendre quadrature,", qpoints, "points",collapse=""))
} else{
nint<-noquote(paste("Gauss Kronrod quadrature,", qpoints, "points",collapse=""))
}
output<-list("Survival sub-model"=s.model,
   "Longitudinal sub-model"=l.model,
    "Association structure"=method0,
    stan_fit=sum.stanfit,
    time=noquote(paste(round(difftime(end_time,start_time,units="mins")[[1]],digits=2),"mins",collapse="")),
    integration=nint,all.results=results)
class(output)<-"JMR"
attr(output, "hidden") <-c("integration","all.results")
return(output)
} else{
stop("Stan can't start. Changing initial values may help")
}
} else{
stop("Stan can't start. Changing initial values may help")
}
}

