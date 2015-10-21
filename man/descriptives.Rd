\name{descriptives}
\alias{summary.freqtab}
\alias{mean.freqtab}
\alias{sd.freqtab}
\alias{var.freqtab}
\alias{cov.freqtab}
\alias{cor.freqtab}
\alias{min.freqtab}
\alias{max.freqtab}
\alias{range.freqtab}
\alias{skew.freqtab}
\alias{kurt.freqtab}
\title{Descriptive Statistics for Frequency Tables}
\description{
These functions return descriptive statistics for a frequency table of class \dQuote{\code{\link{freqtab}}}.
}
\usage{
\method{summary}{freqtab}(object, margin = seq(margins(object)), ...)

\method{mean}{freqtab}(x, margin = 1, ...)

sd.freqtab(x, margin = 1)

var.freqtab(x, margin = 1)

cov.freqtab(x, margin = seq(margins(x)))

cor.freqtab(x, margin = seq(margins(x)))

\method{min}{freqtab}(x, margin = 1, ..., na.rm = FALSE)

\method{max}{freqtab}(x, margin = 1, ..., na.rm = FALSE)

\method{range}{freqtab}(x, margin = 1, ..., na.rm = FALSE)

skew.freqtab(x, margin = 1)

kurt.freqtab(x, margin = 1)

}
\arguments{
  \item{object, x}{
object of class \dQuote{\code{freqtab}}.
}
  \item{margin}{
integer vector specifying the margin(s) for which summary statistics will be returned. This defaults to \code{1} for univariate statistics, and \code{seq(margins(x))}, i.e., all the margins, for multivariate statistics (covariance and correlation).
}
  \item{\dots}{
further arguments passed to or from other methods.
}
  \item{na.rm}{
logical indicating whether missing values should be removed, currently ignored since frequency tables cannot contain missing values.
}
}
\details{
\code{mean}, \code{sd.freqtab}, \code{var.freqtab}, \code{skew.freqtab}, and \code{kurt.freqtab} return the mean, standard deviation, variance, skewness, and kurtosis. \code{min} and \code{max} return the minimum and maximum observed scores, and \code{range} returns both. \code{cov.freqtab} and \code{cor.freqtab} return the covariance and correlation matrices for one or more variables.

\code{summary} returns univariate statistics across one or more margins.
}
\value{
\code{summary} returns a data frame of summary statistics, including the mean, standard deviation, skewness, kurtosis, minimum, maximum, and number of observations for each variable in \code{margin}. Otherwise, a vector of length \code{length(margin)} is returned with the corresponding statistic for each variable.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{freqtab}}
}
\examples{
summary(as.freqtab(ACTmath[, 1:2]))

ny <- freqtab(KBneat$y, scales = list(0:36, 0:12))
summary(ny)
cov.freqtab(ny)
}
\keyword{methods}
