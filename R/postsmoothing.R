#' Postsmoothing Equating Functions
#' 
#' These functions are used to smooth equating output.
#' 
#' Start with cubic spline postsmoothing.
#' 
#' @author Anthony Albano \email{tony.d.albano@@gmail.com}
#' @seealso \code{\link{presmoothing}}
#' @references Kolen, M. J., and Brennan, R. L. (2004). \emph{Test Equating,
#' Scaling, and Linking}. (2nd ed.), New York: Springer.
#' @examples
#'
#' # Create some examples here
#'
#' # From presmoothing examples 
#' # Bivariate example, preserving first 3 moments of total
#' # and anchor for x and y, and the covariance
#' # between anchor and total
#' # see equated scores in Wang (2009), Table 4
#' xvtab <- freqtab(KBneat$x, scales = list(0:36, 0:12))
#' yvtab <- freqtab(KBneat$y, scales = list(0:36, 0:12))
#' Y <- as.data.frame(yvtab)[, 1]
#' V <- as.data.frame(yvtab)[, 2]
#' scorefun <- cbind(Y, Y^2, Y^3, V, V^2, V^3, V*Y)
#' wang09 <- equate(xvtab, yvtab, type = "equip",
#'   method = "chained", smooth = "loglin",
#'   scorefun = scorefun)
#' 
#' @export
postsmoothing <- function(x, ...) UseMethod("postsmoothing")

#' @describeIn postsmoothing Default method for vector of equated scores.
#' @export
postsmoothing.default <- function(x, smoothmethod = c("none", "spline"),
  ...) {

	smoothmethod <- match.arg(smoothmethod)
	if (smoothmethod == "none")
		return(x)
	else if (smoothmethod == "spline")
		return(espline(x))
}

#----------------------------------------------------------------
# Exported internal function for spline smoothing

#' @rdname postsmoothing
#' @export
espline <- function() {
  
}
