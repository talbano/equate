% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.freqtab.R
\name{plot.freqtab}
\alias{plot.freqtab}
\alias{points.freqtab}
\alias{bfreqplot}
\alias{ufreqplot}
\title{Plotting Frequency Distributions}
\usage{
\method{plot}{freqtab}(
  x,
  y = NULL,
  xcol = 1,
  ycol,
  pch = 16,
  ylty = 1,
  xlab = names(dimnames(x))[1],
  addlegend = !missing(y),
  legendtext,
  ...
)

\method{points}{freqtab}(x, xcol = 1, pch = 16, ds = 50, dm = 100, ...)
}
\arguments{
\item{x}{univariate or bivariate score distribution of class
\dQuote{\code{\link{freqtab}}}.}

\item{y}{either an object of class \dQuote{\code{freqtab}}, where
frequencies will be extracted, or a vector or matrix of frequencies, to be
added to the plot of \code{x}. See below for details.}

\item{xcol, ycol}{colors used in plotting \code{x} and \code{y}.}

\item{pch}{plotting symbol used to plot bivariate points.}

\item{ylty}{line type used to plot frequencies in \code{y}.}

\item{xlab}{label for the x axis.}

\item{addlegend}{logical indicating whether or not a legend should be added.}

\item{legendtext}{character vector of text to be passed to the \code{legend}
argument of the \code{legend} function, defaulting to column names used in
\code{y}.}

\item{\dots}{further arguments passed to or from other methods, such as
graphical parameters besides \code{col}, \code{type}, and \code{pch}.}

\item{ds, dm}{integers for the scaling and center of the RGB density values,
with defaults of 50 and 100. These are used to convert the observed counts
in \code{x} to the [0, 255] range of RGB values.}
}
\value{
The univariate option produces a single line plot of \code{type =
"h"}. Frequencies from \code{y} are then superimposed. The bivariate option
produces a scatterplot with a marginal frequency plot for each distribution.
}
\description{
This function plots univariate and bivariate frequency tables of class
\dQuote{\code{\link{freqtab}}}.
}
\details{
For the points method, a scatterplot for \code{x} is added to the current
opened plot.

For the plot method, when \code{x} is univariate, i.e, having 2 columns, a
frequency plot is created for \code{x}. When \code{x} is bivariate, e.g.,
coming from a single group equating design or one form of a nonequivalent
groups design, a scatterplot is produced with frequency plots for the
marginal distributions.

\code{y} is used to superimpose lines, e.g., smoothed frequencies, over the
(marginal) frequencies of \code{x}.

Colors must be specified using \code{xcol} and \code{ycol}. When \code{ycol}
is missing, a vector of colors is created using \code{rainbow(ncol(y))}.
}
\examples{

x <- freqtab(KBneat$x, scales = list(0:36, 0:12))
plot(x)

xs <- loglinear(x, degrees = c(4, 1),
  stepup = TRUE, showWarnings = FALSE)
plot(x, xs, lwd = 2)

}
\seealso{
\code{\link{plot.table}}, \code{\link{plot.equate}},
\code{\link{lines}}, \code{\link{points}}
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\keyword{misc}
