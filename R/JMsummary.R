########################################################################
#' Summary of a jm.reg fit
#' @description Produces a summary of the Bayesian joint model fit. It uses a \code{\link{jm.reg}} fit as its argument.
#' @keywords Joint modeling, MCMC
#' @param jmfit a \code{\link{jm.reg}} fit with \code{method = "MCMC"}.
#' @details This function uses the posterior draws and 
#'     the posterior summary from class "stanfit".
#' @return A summary of the survival data, longitudinal data, the joint model, number of MCMC chains,
#'    number of iterations per chain after warmup, and computaional time. It also gives
#'
#' @return \code{Posterior summary (survival sub-model):} posterior summary for the survival sub-model.
#' @return \code{Posterior summary (exp(coef) of the survival sub-model):} exp(coef) for the survival sub-model (for PH models only).
#' @return \code{Posterior summary (longitudinal sub-model):} posterior summary for the longitudinal sub-model.
#' @author Shahedul Khan <khan@math.usask.ca>
#' @seealso \code{\link{jm.icriteria}}, \code{\link{jm.plots}}, \code{\link{jm.reg}},
#'          \code{\link{jm.resid}}
#' @examples 
#'   # Example: pbc data
#'   lme.fit <- lme(log(bilirubin) ~ drug + ns(futime, 2),
#'     data = pbc.long, random = ~ ns(futime, 2) | id)
#'   surv.fit <- coxph(Surv(st, status2) ~ drug * age, data = pbc.surv, x = TRUE)
#'   jmfit.ll2 <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "llogistic", timevar = "futime", form = "riz")
#'   jm.summary(jmfit.ll2)
#'   jmfit.wph2 <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
#'       surv.model = "weibullph", timevar = "futime", form = "riz")
#'   jm.summary(jmfit.wph2)
#' @export

jm.summary<-function(jmfit){
if(class(jmfit)!="JMR"){
stop("jmfit must inherit from class JMR")
}
results<-jmfit$all.results
if(is.null(results)){
stop("\nThis function works only for method = 'MCMC'")
}
posterior <- as.array(results$stan.fit)
chains<-dim(posterior)[2]
surv.data<-results$surv.data
surv.covariates<-colnames(surv.data)[-(1:3)]
long.data<-results$long.data
long.covariates<-colnames(long.data)[-(1:2)]
p2<-results$stan.data$p2
#-----
beta.loc<-results$beta.loc
phi.loc<-results$phi.loc
kappa.loc<-results$kappa.loc
if(results$surv.model=="eweibull" || results$surv.model=="ggamma" || results$surv.model=="gll"){
gam.loc<-results$gam.loc
}
rho.loc<-results$rho.loc
alpha.loc<-results$alpha.loc
sigma.loc<-results$sigma.loc
cov.loc<-results$cov.loc
#----------
sd.cor<-function(vv,col2){
v.mat<-matrix(vv,ncol=sqrt(col2))
v.cor0<-cov2cor(v.mat)
v.cor0[lower.tri(v.cor0,diag=T)] <- NA
v.cor1<-c(na.omit(c(t(v.cor0))))
v.sd<-sqrt(diag(v.mat))
return(c(v.cor1,v.sd))
}
L.corr<-list()
L.sd<-list()
exp.beta<-list()
for(i in 1:chains){
if(results$surv.model=="weibullph" || results$surv.model=="gllph"){
exp.beta[[i]]<-exp(posterior[,i,][,c(beta.loc,phi.loc)])
}
L.sig1<-posterior[,i,][,cov.loc]
col2<-ncol(L.sig1)
L.cor.sd<-t(apply(L.sig1,1,sd.cor,col2=col2))
sd.ran<-L.cor.sd[,(ncol(L.cor.sd)-sqrt(col2)+1):ncol(L.cor.sd)]
L.cor2<-L.cor.sd[,1:(ncol(L.cor.sd)-sqrt(col2))]
if(is.vector(L.cor2)){L.cor2<-matrix(L.cor2,ncol=1)}
colnames(sd.ran)<-paste(rep("sd(b",ncol(sd.ran)),c(0:(ncol(sd.ran)-1)),rep(")",ncol(sd.ran)),sep="")
L.sd[[i]]<-sd.ran
vm<-sqrt(ncol(L.sig1))
mat <- matrix(0, vm, vm)
mat1<-which(mat ==0, arr.ind = T)
mat1[,1:2]<-mat1[,2:1]
mat2<-matrix(mat1[mat1[,"row"]<mat1[,"col"],],ncol=2)-1
gg<-NULL
for(k in 1:(vm*(vm-1)/2)){
gg<-c(gg,paste("cor(b",paste(mat2[k,1]),", b",paste(mat2[k,2]),")",sep = ''))
}
colnames(L.cor2)<-gg
L.corr[[i]]<-L.cor2
}
#--------
L.corr.mcmc<-as.mcmc.list(lapply(L.corr,mcmc))
L.corr.sum<-summary(L.corr.mcmc)
R.long.corr<-round(gelman.diag(L.corr.mcmc,multivariate=FALSE,
  autoburnin=FALSE)[[1]][,1],digits=2)
effn.corr<-effectiveSize(L.corr.mcmc)
#---------
L.sd.mcmc<-as.mcmc.list(lapply(L.sd,mcmc))
L.sd.sum<-summary(L.sd.mcmc)
R.long.sd<-round(gelman.diag(L.sd.mcmc,multivariate=FALSE,autoburnin=FALSE)[[1]][,1],digits=2)
effn.sd<-effectiveSize(L.sd.mcmc)
#----------
rfit<-jmfit$stan_fit
rownames(rfit)[beta.loc]<-surv.covariates
rownames(rfit)[phi.loc]<-"association"
#----------
if(results$surv.model=="weibullph" || results$surv.model=="gllph"){
exp.beta.mcmc<-as.mcmc.list(lapply(exp.beta,mcmc))
exp.beta.sum<-summary(exp.beta.mcmc)
R.exp.beta<-round(gelman.diag(exp.beta.mcmc,multivariate=FALSE,
  autoburnin=FALSE)[[1]][,1],digits=2)
effn.exp.beta<-effectiveSize(exp.beta.mcmc)
ebeta.sum<-cbind(exp.beta.sum[[1]][,1],exp.beta.sum[[1]][,4],
  exp.beta.sum[[1]][,2],exp.beta.sum[[2]],effn.exp.beta,R.exp.beta)
colnames(ebeta.sum)<-colnames(rfit)
exp.beta.rnames<-NULL
for(i in 1:length(surv.covariates)){
exp.beta.rnames[i]<-paste("exp(",surv.covariates[i],")",sep="")
}
rownames(ebeta.sum)[1:p2]<-exp.beta.rnames
rownames(ebeta.sum)[p2+1]<-"exp(association)"
}
#---------
rownames(rfit)[kappa.loc]<-"kappa"
rownames(rfit)[rho.loc]<-"rho"
if(results$surv.model=="eweibull" || results$surv.model=="ggamma"){
rownames(rfit)[gam.loc]<-"gamma"
}
rownames(rfit)[alpha.loc]<-long.covariates
rownames(rfit)[sigma.loc]<-"sd(resid)"
#-----------
surv.output<-rfit[1:(alpha.loc[1]-1),]
#-----------
sd.sum<-cbind(L.sd.sum[[1]][,1],L.sd.sum[[1]][,4],
  L.sd.sum[[1]][,2],L.sd.sum[[2]],effn.sd,R.long.sd)
colnames(sd.sum)<-colnames(rfit)
if(!is.matrix(L.corr.sum[[1]])){
corr.sum<-matrix(c(L.corr.sum[[1]][1],L.corr.sum[[1]][4],
  L.corr.sum[[1]][2],L.corr.sum[[2]],effn.corr,R.long.corr),nrow=1)
colnames(corr.sum)<-colnames(rfit)
rownames(corr.sum)<-gg
} else{
corr.sum<-cbind(L.corr.sum[[1]][,1],L.corr.sum[[1]][,4],
  L.corr.sum[[1]][,2],L.corr.sum[[2]],effn.corr,R.long.corr)
colnames(corr.sum)<-colnames(rfit)
rownames(corr.sum)<-gg
}
#-------
long.output0<-rfit[c(alpha.loc,sigma.loc),]
long.output<-rbind(long.output0,sd.sum,corr.sum)
#------
s.model<-results$s.model
lme.fit<-results$lme.fit
id <- as.vector(unclass(lme.fit$groups[[1]]))
surv.data<-results$surv.data
cens<-sum(1-surv.data[,"status"])
cens.p<-round(cens*100/nrow(surv.data),digits=1)
cat(paste("Survival data:"),
    paste("          ","number of observations = ",nrow(surv.data)),
    paste("          ","censoring = ",cens,"(",cens.p,"%)"),
    paste("Longitudinal data:"), 
    paste("          ","number of observations = ",length(id)),
    paste("          ","number of individuals = ",length(unique(id))),paste(" "),sep="\n")
if(results$surv.model=="gllph"){
cat(paste("Survival sub-model:",s.model[1,],"with",s.model[2,],sep=" "),
   paste("Longitudinal sub-model:","LME with b0, b1, ... = random-effects coefficients",sep=" "),
   paste("Association Structure: Induced by the longitudinal value"),paste(" "),sep="\n")
} else{
cat(paste("Survival sub-model:",s.model[1,],"with",s.model[2,],"and",s.model[3,],sep=" "),
   paste("Longitudinal sub-model:","LME with b0, b1, ... = random-effects coefficients",sep=" "),
   paste("Association Structure: Induced by the longitudinal value"),paste(" "),sep="\n")
}
cat(paste("Number of chains =",chains), 
   paste("Number of iterations per chain after warmup =",dim(posterior)[1]),
   paste("Time =",jmfit$time),paste(""),sep="\n")
if(results$surv.model=="weibullph" || results$surv.model=="gllph"){
output<-list("Posterior summary (survival sub-model)"=surv.output,
   "Posterior summary (exp(coef) of the survival sub-model)"=ebeta.sum,
   "Posterior summary (longitudinal sub-model)"=long.output)
} else{
output<-list("Posterior summary (survival sub-model)"=surv.output,
   "Posterior summary (longitudinal sub-model)"=long.output)
}
class(output)<-"JMR"
attr(output, "hidden")<-NULL
return(output)
}


