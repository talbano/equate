\name{px}
\alias{px}
\alias{px.default}
\alias{px.freqtab}
\alias{fx}
\alias{fx.default}
\alias{fx.freqtab}
\title{Percentile Ranks and Cumulativie Frequencies}
\description{
These functions compute percentile ranks and cumulative frequency distributions for frequency tables.
}
\usage{
\method{px}{default}(x, y, ys, ...)

\method{px}{freqtab}(x, margin = 1, y, ymargin = 1, ...)

\method{fx}{default}(x, ...)

\method{fx}{freqtab}(x, margin = 1, ...)
}
\arguments{
  \item{x}{
either a vector of counts, or an object of class \dQuote{\code{freqtab}} from which counts will be taken.
}
  \item{y}{
an object of class \dQuote{\code{freqtab}} when \code{x} is as well, otherwise, a vector or \code{data.frame} of counts. See below for details.
}
  \item{ys}{
vector specifying the \code{y} score scale, when it is not contained in the first column of \code{y}. If \code{y} can be converted to a \code{data.frame}, it is assumed to be univariate with the first column containing the score scale and the second containing the counts.
}
  \item{margin, ymargin}{
integers specifying the margins for which frequencies or percentile ranks will be returned. \code{margin} applies to \code{x} and \code{ymargin} to \code{y}.
}
  \item{\dots}{
further arguments passed to or from other methods.
}
}
\details{
These functions compute percentile ranks and cumulative frequencies for a univariate distribution, and percentile ranks from one univariate distribution (\code{x}) corresponding to score values in another (\code{y}).
}
\value{
A vector is returned containing either percentile ranks or cumulative frequencies with length equal to \code{length(x)}.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{freqtab}}
}
\examples{
x <- as.freqtab(ACTmath[, 1:2], drop = TRUE)
y <- as.freqtab(ACTmath[, c(1, 3)], drop = TRUE)

# Percentile ranks for the x scale
round(px(x), 3)

# Percentile ranks in y for x each score
round(px(x, y = y), 3)

# Cumulative frequency distribution for x
round(fx(x), 3)
}
\keyword{univar}