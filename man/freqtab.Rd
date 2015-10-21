\name{freqtab}
\alias{freqtab}
\alias{freqtab.default}
\alias{as.freqtab}
\alias{is.freqtab}
\alias{as.freqtab.table}
\alias{freqtab.table}
\alias{as.data.frame.freqtab}
\alias{head.freqtab}
\alias{tail.freqtab}
\alias{print.freqtab}
\alias{scales}
\alias{margin}
\alias{margins}
\alias{droplevels.freqtab}
\title{Frequency Distribution Tables}
\description{
Functions for creating and manipulating frequency tables of class \dQuote{\code{freqtab}}.
}
\usage{
\method{freqtab}{default}(x, scales, items, na.rm = TRUE, ...)

\method{freqtab}{table}(x, ...)

as.freqtab(x, scales, drop = FALSE, ...)

is.freqtab(x)

\method{as.data.frame}{freqtab}(x, row.names = NULL,
  optional = FALSE, drop = FALSE, ...)

scales(x, margin = 1)

margin(x, margin = 1)

margins(x)

\method{droplevels}{freqtab}(x, ...)

\method{head}{freqtab}(x, ...)

\method{tail}{freqtab}(x, ...)
}
\arguments{
  \item{x}{
either an object (vector or \code{data.frame}) containing total scores or item responses with which total scores will be calculated, or an object inheriting from class \dQuote{\code{freqtab}}. In the \code{freqtab} function, the first column in \code{x} must be the total test; any remaining columns may contain anchor scores. See below for details.
}
  \item{scales}{
list of vectors containing the score scales for each column in \code{x}. \code{length(scales)} should equal \code{ncol(x)}.
}
  \item{items}{
list of vectors of column indices (numbers or column names) with which total scores will be computed for \code{x} when it contains item responses.
}
  \item{na.rm}{
logical with default \code{TRUE} specifying whether or not missing item responses should be ignored, i.e., treated as 0, when calculating total scores.
}
  \item{row.names, optional}{
arguments passed to \code{as.data.frame}, currently ignored.
}
  \item{drop}{
logical, with default \code{FALSE}, indicating whether or not unused factor levels, or score values with zero counts, should be dropped. See below for details.
}
  \item{margin}{
integer vector specifying the margin(s) over which frequencies should be summed.
}
  \item{\dots}{
  further arguments passed to or from other functions.
}
}
\details{
\code{freqtab} creates a frequency table from a vector or \code{data.frame} of scores. When the \code{items} argument is included, scores are assumed to be item responses, which are summed to create total scores. The scores are tabulated and stored as an array, with dimensions for each variable. Note that in previous versions of the \dQuote{\code{freqtab}} class the frequency table was stored as a \code{data.frame}. This is no longer the case. Instead, the table is stored as an array and converted to a \code{data.frame} when printed or manipulated with the \code{head} and \code{tail} methods.

\code{as.data.frame} converts an object of class \dQuote{\code{freqtab}} to \dQuote{\code{data.frame}}. \code{droplevels} returns \code{x} with any unused factor levels, or levels with zero counts, removed.

When \code{x} is an object of class \dQuote{\code{table}}, \code{freqtab} simply modifies the attributes and converts to class \dQuote{\code{freqtab}}. In this case, \code{x} must already be structured similar to a \dQuote{\code{freqtab}} object, with the first dimension containing counts for total scores, and remaining dimensions containing counts for one or more anchor tests.

\code{as.freqtab} converts a 'flat' contingency table (see \code{\link{ftable}}) to class \dQuote{\code{freqtab}} with the appropriate attributes. A flat contingency table is the \code{data.frame} version of a \dQuote{\code{freqtab}} object, where the first column contains the total score scale, the last column contains counts, and the columns in between contain different anchor test score combinations. \code{is.freqtab} tests for class \dQuote{\code{freqtab}}.

\code{scales} extracts the measurement scales for the variables specified in \code{margin}, with \code{margin = 1} referring to the total score scale, and subsequent margins referring to anchor tests. \code{margin} is a wrapper for \code{\link{margin.table}}, which itself is a simple wrapper for summing over marginal counts, i.e., \code{apply(x, margin, sum)}. And \code{margins} returns the number of dimensions, i.e., score variables, in a frequency table.

The main difference between the \dQuote{\code{freqtab}} class and other tabulation classes, like \dQuote{\code{table}} and \dQuote{\code{ftable}}, is that the \code{dimnames}, i.e., the score scales, are required to be numeric. This facilitates plotting with \code{\link{plot.freqtab}}, equating with the \code{\link{equate}} function, and descriptive statistics with the \code{\link{summary.freqtab}} and other methods.
}
\value{
A table array with dimensions equal to the number of score scales. In most cases, this will be a univariate or bivariate distribution, but multivariate distributions are supported. \code{scales} and \code{margins} return numeric vectors.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{table}}, \code{\link{ftable}}, \code{\link{summary.freqtab}}, \code{\link{plot.freqtab}}
}
\examples{
# Univariate distribution with score scale
set.seed(2005)
x <- round(rnorm(1000, 100, 10))
head(freqtab(x, scales = 70:130))

# Existing frequency table converted to class "freqtab"
# The first score of zero, with zero counts, is dropped
head(as.freqtab(ACTmath[, 1:2], drop = TRUE))

# Bivariate distribution
# Reduce y to the anchor test margin (2)
ny <- freqtab(x = KBneat$y, scales = list(0:36, 0:12))
margin(ny, margin = 2)

# Summing scored item responses with the PISA data
attach(PISA)
r6items <- paste(items$itemid[items$clusterid == "r6"])
r7items <- paste(items$itemid[items$clusterid == "r7"])
pisa67 <- freqtab(students[students$book == 6, ],
	items = list(r6items, r7items),
	scales = list(0:16, 0:14))
detach(PISA)

# Scales for both margins
# Zero total score is unobserved
scales(pisa67, 1:2)
scales(droplevels(pisa67), 1:2)
}
\keyword{methods}
