####################################################
# LR test for weibull as a submodel of eweibull or ggamma
###################################################
#' Likelihood ratio test
#' @description Likelihood ratio test for Weibull AFT as a submodel of the exponentiated Weibull AFT and
#'      generalized gamma AFT models. 
#' @keywords Exponentiated Weibull, Generalized gamma, Likelihood ratio test, Weibull
#' @param fit a list of two \code{etime.reg} fits: one must be the Weibull AFT fit and the other one is the exponentiated
#'      Weibull AFT fit or the generalized gamma AFT fit.
#' @details The log-likelihood values from the two fits are used for the likelihood ratio test.
#' @return \code{chi.sq:} chi-square test statistic.
#' @return \code{p value:} p value.
#' @author Shahedul Khan <khan@math.usask.ca>
#' @seealso \code{\link{etime.reg}}, \code{\link{etime.resid}}.
#' @examples 
#'   library(frailtypack)
#'   data(readmission)
#'   fit.gg <- etime.reg(c(t.start, t.stop, event) ~ sex + chemo + charlson,
#'       surv.model = "ggamma", Data = readmission)
#'   fit.ew <- etime.reg(c(t.start, t.stop, event) ~ sex + chemo + charlson,
#'       surv.model = "eweibull", Data = readmission)
#'   fit.w <- etime.reg(c(t.start, t.stop, event) ~ sex + chemo + charlson,
#'       surv.model = "weibull", Data = readmission)
#'   LR.test(list(fit.gg, fit.w))
#'   LR.test(list(fit.ew, fit.w))
#' @export

LR.test<-function(fit){
for(i in 1:length(fit)){
if(fit[[i]]$surv.model=="Rwaft" || fit[[i]]$surv.model=="waft") {L2<-fit[[i]][[5]][1]}
if(fit[[i]]$surv.model=="Rewaft" || fit[[i]]$surv.model=="ewaft" ||
   fit[[i]]$surv.model=="Rggaft" || fit[[i]]$surv.model=="ggaft") {L1<-fit[[i]][[5]][1]}
}
stat<-2*(L1-L2)
p.val<-pchisq(stat,1,lower.tail=FALSE)
results<-cbind(chi.sq=stat,"p value"=p.val)
rownames(results)<-""
return(results)
}
