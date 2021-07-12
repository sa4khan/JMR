#################################################
#' Cox-Snell residuals
#' @description Cox-Snell residuals of the \code{etime.reg} fit.
#' @keywords Cox-Snell residual
#' @param fit \code{etime.reg} fit.
#' @param plot returns residual plot if \code{TRUE} (default is \code{FALSE}).
#' @param conf.int confidence level for the pointwise confidence intervals of the residuals (if \code{plot = TRUE}).
#' @param xlim x limits (a vector of length 2 of the form \code{c(minimum,maximum)}) of the plot (optional).
#' @param ylim y limits (a vector of length 2 of the form \code{c(minimum,maximum)}) of the plot (optional).
#' @param main main title (character expression) of the plot (optional).
#' @details If \code{plot = TRUE}, this function returns both residuals and a residual plot. 
#'     The plot shows the unit slop line (in grey), residuals (solid circles) and pointwise confidence intervals (dashed lines).
#'     The plot of the Cox-Snell (hazard-based) residuals should be roughly a straight line 
#'     with unit slope when the model is adequate. 
#'     See Section 6.2 of Lawless and Section 3.7.3 of Cook and Lawless for details.
#' @return \code{residuals} Cox-Snell residuals.
#' @return \code{plot} residual plot if \code{plot=TRUE}.
#' @references Cook RJ and Lawless J, The statistical analysis of recurrent events, Springer, 2007.
#' @references Lawless J, Statistical models and methods for lifetime data, Wiley, 2003.
#' @author Shahedul Khan <khan@math.usask.ca>
#' @seealso \code{\link{etime.reg}}, \code{\link{LR.test}}.
#' @examples 
#' # Example 1: Recurrent Event Analysis
#'   library(frailtypack)
#'   data(readmission)
#'   fit.gg <- etime.reg(c(t.start, t.stop, event) ~ sex + chemo + charlson,
#'       surv.model = "ggamma", Data = readmission)
#'   fit.wph <- etime.reg(c(t.start, t.stop, event) ~ sex + chemo + charlson,
#'       surv.model = "weibullph", Data = readmission)
#'   etime.resid(fit.wph, plot = TRUE)
#'   par(mfrow = c(1, 2))
#'   etime.resid(fit.wph, plot = TRUE, xlim = c(0, 5.25),
#'     ylim = c(0, 5.25),  main = "Weibull PH")
#'   etime.resid(fit.gg, plot = TRUE, xlim = c(0, 5.25),
#'      ylim = c(0, 5.25),  main = "Generalized Gamma AFT")
#'
#' # Example 2: Non-recurrent Data (Right Censored)
#'   library(survival)
#'   fit.ln <- etime.reg(c(time, status) ~ karno + diagtime + age + prior + celltype + trt,
#'       Data = veteran, surv.model = "lnormal")
#'   etime.resid(fit.ln, plot = TRUE, main = "Log-normal AFT")
#' @export

etime.resid<-function(fit,plot=FALSE,conf.int=0.95,xlim=NULL,ylim=NULL,
    xlab=NULL,ylab=NULL,main=NULL){
if(class(fit)!="JMR"){
stop("'fit' must inherit from class JMR")
}
if(fit$surv.model=="ewaft") {resid<-ewaft.resid(fit,plot=plot,conf.int=conf.int,
  xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)}
if(fit$surv.model=="Rewaft") {resid<-Rewaft.resid(fit,plot=plot,conf.int=conf.int,
  xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)}
if(fit$surv.model=="waft") {resid<-waft.resid(fit,plot=plot,conf.int=conf.int,
  xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)}
if(fit$surv.model=="Rwaft") {resid<-Rwaft.resid(fit,plot=plot,conf.int=conf.int,
  xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)}
if(fit$surv.model=="ggaft") {resid<-ggaft.resid(fit,plot=plot,conf.int=conf.int,
  xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)}
if(fit$surv.model=="Rggaft") {resid<-Rggaft.resid(fit,plot=plot,conf.int=conf.int,
  xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)}
if(fit$surv.model=="llaft") {resid<-llaft.resid(fit,plot=plot,conf.int=conf.int,
  xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)}
if(fit$surv.model=="Rllaft") {resid<-Rllaft.resid(fit,plot=plot,conf.int=conf.int,
  xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)}
if(fit$surv.model=="lnaft") {resid<-lnaft.resid(fit,plot=plot,conf.int=conf.int,
  xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)}
if(fit$surv.model=="Rlnaft") {resid<-Rlnaft.resid(fit,plot=plot,conf.int=conf.int,
  xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)}
if(fit$surv.model=="wph") {resid<-wph.resid(fit,plot=plot,conf.int=conf.int,
  xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)}
if(fit$surv.model=="Rwph") {resid<-Rwph.resid(fit,plot=plot,conf.int=conf.int,
  xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)}
if(fit$surv.model=="gllph") {resid<-gll.resid(fit,plot=plot,conf.int=conf.int,
  xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)}
if(fit$surv.model=="Rgllph") {resid<-Rgll.resid(fit,plot=plot,conf.int=conf.int,
  xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main)}
return(resid)
}
