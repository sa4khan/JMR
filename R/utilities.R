#--------------------------
# log1mexp gives the value of log(1-exp(x))
#---------------------------
log1mexp <- function(x) {
ifelse(x <= log(2), log(-expm1(x)), log1p(-exp(x)))
}
#-------------------------
# log1pexp gives the value of log(1+exp(x))
#-------------------------
log1pexp <- function(x){
ifelse(x <= -37, exp(x),
ifelse(x <= 18, log1p(exp(x)),
 ifelse(x <= 33, x + exp(-x), x)))
}
#------------------------------
# Density of multivariate normal
#------------------------------
den.mnorm<-function (x, mean = rep(0, p), sigma = diag(p), log = FALSE) {
if (is.vector(x)) 
x <- matrix(x, ncol = length(x))
p <- ncol(x)
if (!missing(mean)) {
if (!is.null(dim(mean))) 
dim(mean) <- NULL
if (length(mean) != p) 
stop("mean and sigma have non-conforming size")
}
if (!missing(sigma)) {
if (p != ncol(sigma)) 
stop("x and sigma have non-conforming size")
if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
check.attributes = FALSE)) 
stop("sigma must be a symmetric matrix")
}
dec <- tryCatch(chol(sigma), error = function(e) e)
if (inherits(dec, "error")) {
x.is.mu <- colSums(t(x) != mean) == 0
logretval <- rep.int(-Inf, nrow(x))
logretval[x.is.mu] <- Inf
}
else {
tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
rss <- colSums(tmp^2)
logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
}
names(logretval) <- rownames(x)
if (log) 
logretval
else exp(logretval)
}
#-----------------
# Simulate from multivariate normal
# -----------------
sim.mnorm<-function (n=1,mu,Sigma,tol = 1e-06,empirical = FALSE,EISPACK = FALSE) {
p <- length(mu)
if (!all(dim(Sigma) == c(p, p))) 
stop("incompatible arguments")
if (EISPACK) 
stop("'EISPACK' is no longer supported by R", domain = NA)
eS <- eigen(Sigma, symmetric = TRUE)
ev <- eS$values
if (!all(ev >= -tol * abs(ev[1L]))) 
stop("'Sigma' is not positive definite")
X <- matrix(rnorm(p * n), n)
if (empirical) {
X <- scale(X, TRUE, FALSE)
X <- X %*% svd(X, nu = 0)$v
X <- scale(X, FALSE, TRUE)
}
X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
nm <- names(mu)
if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) 
nm <- dn[[1L]]
dimnames(X) <- list(nm, NULL)
if (n == 1) 
drop(X)
else t(X)
}
#------------------------
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
}
#--------------------
itexp <- function(u, m, t) { -log(1-u*(1-exp(-t*m)))/m }
rtexp <- function(n, m, t) { itexp(runif(n), m, t) }
