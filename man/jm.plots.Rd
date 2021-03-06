% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/JMplots.R
\name{jm.plots}
\alias{jm.plots}
\title{MCMC plots for posterior analysis}
\usage{
jm.plots(
  jmfit,
  pairs = TRUE,
  rhat = FALSE,
  neff = FALSE,
  trace = FALSE,
  density = FALSE,
  uncertainty.intervals = FALSE
)
}
\arguments{
\item{jmfit}{a \code{\link{jm.reg}} fit with \code{method = "MCMC"}.}

\item{pairs}{if \code{TRUE}, produce pair plots.}

\item{rhat}{if \code{TRUE}, produce a plot for rhats.}

\item{neff}{if \code{TRUE}, produce a plot for ratios of effective sample size to total sample size.}

\item{trace}{if \code{TRUE}, produce trace plots.}

\item{density}{if \code{TRUE}, produce density plots.}

\item{uncertainty.intervals}{if \code{TRUE}, produce uncertainty intervals computed from parameter draws.}
}
\value{
Plots for posterior analysis.
}
\description{
This function produces plots for posterior analysis. 
   The 'bayesplot' package is used for plotting. Available options are
   pair plots, rhat values as either points or a histogram, ratios of effective sample size to total sample size
    as either points or a histogram, trace plots, density plots, and plots of uncertainty intervals. See the
   'bayesplot' package for detail.
}
\details{
The 'bayesplot' package is used for plotting. See the
   'bayesplot' package for detail.
}
\examples{
# Example: AIDS data from package JM
  library(JM)
  data(aids.id)
  data(aids)
  surv.fit <- coxph(Surv(Time, death) ~ drug + gender + prevOI + AZT,
          data = aids.id, x = TRUE)
  lme.fit <- lme(CD4 ~ obstime + obstime:drug + gender + prevOI + AZT,
      random =  ~ obstime | patient, data = aids)
  # Weibull PH with form = "riz"
  jmfit.wph0 <- jm.reg(surv.fit = surv.fit, lme.fit = lme.fit,
      surv.model = "weibullph", timevar = "obstime", form = "riz")
  jm.plots(jmfit.wph0, pairs = TRUE, rhat = TRUE, neff = TRUE,
    trace = TRUE, density = TRUE, uncertainty.intervals = TRUE)
}
\references{
Gabry J and Mahr T, bayesplot: Plotting for Bayesian Models, R
   package version 1.7.2, \url{https://mc-stan.org/bayesplot}
}
\seealso{
\code{\link{jm.reg}}, \code{\link{jm.resid}}
}
\author{
Shahedul Khan <khan@math.usask.ca>
}
\keyword{Density}
\keyword{Effective}
\keyword{Joint}
\keyword{MCMC,}
\keyword{Pair}
\keyword{Rhat}
\keyword{Trace}
\keyword{Uncertainty}
\keyword{interval.}
\keyword{modeling,}
\keyword{plot,}
\keyword{plots,}
\keyword{sample}
\keyword{size,}
\keyword{values,}
