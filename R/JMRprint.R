#' Prints JMR objects
#' @param x A \code{JMR} object
#' @param digits minimal number of significant digits
#' @return No value is returned.
#' @author Shahedul Khan <khan@math.usask.ca>
#' @export

print.JMR<- function(x,digits=max(options()$digits - 4, 3)) {
hid <- attr(x, "hidden")
print(x[!names(x) %in% hid],digits=digits)
invisible(x)
}