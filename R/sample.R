#' Bootstrap Random Sampling from Frequency Tables
#' 
#' An extension of \code{\link{sample}} to objects of class
#' \dQuote{\code{freqtab}} for bootstrap sampling.
#' 
#' @param x object of class \dQuote{\code{freqtab}}, which is an array of
#' counts across one or more numeric dimensions.
#' @param size non-negative integer giving the sample size.
#' @param replace logical with default \code{TRUE} indicating whether
#' sampling should be with replacement.
#' @return A table array, as a \dQuote{\code{freqtab}} object, sampled
#' from the original \code{x}.
#' @author Anthony Albano \email{tony.d.albano@@gmail.com}
#' @seealso \code{\link{table}}, \code{\link{ftable}},
#' \code{\link{summary.freqtab}}, \code{\link{plot.freqtab}}
#' @keywords methods
#' @examples
#' 
#' # Sample with replacement from ACT math and compare results
#' set.seed(2021)
#' rx <- as.freqtab(ACTmath[, 1:2])
#' rxs <- sample.freqtab(rx)
#' summary(rx)
#' summary(rxs)
#' 
#' @export sample.freqtab
sample.freqtab <- function(x, size = sum(x), replace = TRUE) {

  xd <- as.data.frame(x)[x > 0, 1:margins(x), drop = FALSE]
  xi <- sample.int(nrow(xd), size, replace = replace,
    prob = x[x > 0]/sum(x))
  out <- freqtab(xd[xi, ], scales = scales(x, 1:margins(x)))
  return(out)
}
