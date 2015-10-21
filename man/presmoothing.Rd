\name{Presmoothing}
\alias{presmoothing}
\alias{presmoothing.default}
\alias{presmoothing.formula}
\alias{loglinear}
\alias{freqbump}
\alias{freqavg}
\title{Frequency Distribution Presmoothing}
\description{
These functions are used to smooth frequency distributions.
}
\usage{
\method{presmoothing}{default}(x, smoothmethod = c("none", "average", "bump",
  "loglinear"), jmin, asfreqtab = TRUE, ...)

\method{presmoothing}{formula}(x, data, ...)

loglinear(x, scorefun, degrees = list(4, 2, 2), grid,
  rmimpossible = FALSE, asfreqtab = TRUE, models,
  stepup = !missing(models),
  compare = FALSE, verbose = FALSE,
  showWarnings = TRUE, ...)
}
\arguments{
  \item{x}{
either an object of class \dQuote{\code{freqtab}} specifying a univariate or multivariate score distribution, or a \dQuote{\code{formula}} object.
}
  \item{smoothmethod}{
character string indicating the smoothing method to be used by \code{presmoothing}. \code{"none"} returns unsmoothed frequencies, \code{"bump"} adds a small frequency to each score value, \code{"average"} imputes small frequencies with average values, and \code{"loglinear"} fits loglinear models. See below for details.
}
  \item{jmin}{
for \code{smoothmethod = "average"}, the minimum frequency, as an integer, below which frequencies will be replaced (default is 1).
for \code{smoothmethod = "bump"}, the value to be added to each score point (as a probability, with default 1e-6).
}
  \item{asfreqtab}{
logical, with default \code{TRUE}, indicating whether or not a frequency table should be returned. For \code{smoothmethod = "average"} and \code{smoothmethod = "bump"}, the alternative is a vector of frequencies. For \code{loglinear}, there are other options.
}
  \item{data}{
an object of class \dQuote{\code{freqtab}}.
}
  \item{scorefun}{
matrix of score functions used in loglinear presmoothing, where each column includes a transformation of the score scale or interactions between score scales. If missing, \code{degrees} and \code{xdegree} will be used to construct polynomial score functions.
}
  \item{degrees}{
list of integer vectors, each one indicating the maximum polynomial score transformations to be computed for each variable at a given order of interactions. Defaults (\code{degrees = list(4, 2, 2)}) are provided for up to trivariate interactions. \code{degrees} are ignored if \code{scorefun} or \code{grid} are provided. See below for details.
}
  \item{grid}{
matrix with one column per margin in \code{x} and one row per term in the model. See below for details.
}
  \item{rmimpossible}{
logical, with default \code{FALSE}, indicating whether or not impossible scores should be removed before smoothing.
}
  \item{models}{
integer vector indicating which model terms should be grouped together when fitting multiple models.
}
  \item{stepup}{
logical, with default \code{FALSE}, indicating whether or not multiple (nested) models should be fit. If \code{TRUE} and \code{models} is missing, an attempt will be made to create it using \code{grid} and/or \code{degrees}. Otherwise, in the absence of \code{models}, each column in \code{scorefun} will define a new sequential model.
}
  \item{compare}{
logical, with default \code{FALSE}, indicating whether or not fit for nested models should be compared. If \code{TRUE}, \code{stepup} is also set to \code{TRUE} and only results from the model fit comparison are returned, that is, \code{verbose} is ignored.
}
  \item{verbose}{
logical, with default \code{FALSE}, indicating whether or not full \code{glm} output should be returned.
}
  \item{showWarnings}{
logical, with default \code{TRUE}, indicating whether or not warnings from \code{glm} should be shown. A warning is common when fitting sparse frequency distributions.
}
  \item{\dots}{
further arguments passed to other methods. For \code{presmoothing}, these are passed to \code{loglinear} and include those listed above.
}
}
\details{
Loglinear smoothing is a flexible procedure for reducing irregularities in a frequency distribution prior to equating, where the degree of each polynomial term determines the specific moment of the observed distribution that is preserved in the fitted distribution (see below for examples). The \code{loglinear} function is a wrapper for \code{\link{glm}}, and is used to simplify the creation of polynomial score functions and the fitting and comparing of multiple loglinear models.

\code{scorefun}, if supplied, must contain at least one score function of the scale score values. Specifying \code{degrees} is an alternative to supplying \code{scorefun}. Each list element in \code{degrees} should be a vector equal in length to the number of variables contained in \code{x}; there should also be one such vector for each possible level of interaction between the variables in \code{x}.

For example, the default \code{degrees = list(4, 2, 2)} is recycled to produce \code{list(c(4, 4, 4), c(2, 2, 2), c(2, 2, 2))}, resulting in polynomials to the fourth power for each univariate distribution, to the second power for each two-way interaction, and to the second power for the three-way interaction.

Terms can also be specified with \code{grid}, which is a matrix with each row containing integers specifying the powers for each variable at each interaction term, including main effects. For example, the main effect to the first power for the total score in a bivariate distribution would be \code{c(1, 0)}; the interaction to the second power would be \code{c(2, 2)}.

\code{stepup} is used to run nested models based on subsets of the columns in \code{scorefun}. Output will correspond to models based on columns 1 and 2, 1 through 3, 1 through 4, to 1 through \code{ncol(scorefun)}. This list of polynomial terms is then used to create a \code{grid} using \code{expand.grid}. The \code{grid} can also be supplied directly, in which case \code{degrees} will be ignored.

\code{compare} returns output as an \code{anova} table, comparing model fit for all the models run with \code{stepup = TRUE}.

The remaining smoothing methods make adjustments to scores with low or zero frequencies. \code{smoothmethod = "bump"} adds the proportion \code{jmin} to each score point and then adjusts the probabilities to sum to 1. \code{smoothmethod = "average"} replaces frequencies falling below the minimum \code{jmin} with averages of adjacent values.
}
\value{
When \code{smoothmethod = "average"} or \code{smoothmethod = "bump"}, either a smoothed frequency vector or table is returned. Otherwise, \code{loglinear} returns the following:
\item{}{when \code{compare = TRUE}, an anova table for model fit}
\item{}{when \code{asfreqtab = TRUE}, a smoothed frequency table}
\item{}{when \code{verbose = TRUE}, full \code{glm} output, for all nested models when \code{stepup = TRUE}}
\item{}{when \code{stepup = TRUE} and \code{verbose = FALSE}, a \code{data.frame} of fitted frequencies, with one column per model}
}
\references{
Holland, P. W., and Thayer, D. T. (1987). \emph{Notes on the use of log-linear models for fitting discrete probability distributions} (PSR Technical Rep. No. 87-79; ETS RR-87-31). Princeton, NJ: ETS.

Holland, P. W., and Thayer, D. T. (2000). Univariate and bivariate loglinear models for discrete test score distributions. \emph{Journal of Educational and Behavioral Statistics, 25}, 133--183.

Moses, T. P., and von Davier, A. A. (2006). \emph{A SAS macro for loglinear smoothing: Applications and implications} (ETS Research Rep. No. RR-06-05). Princeton, NJ: ETS.

Moses, T., and Holland, P. W. (2008). \emph{Notes on a general framework for observed score equating} (ETS Research Rep. No. RR-08-59). Princeton, NJ: ETS.

Wang, T. (2009). Standard errors of equating for the percentile rank-based equipercentile equating with log-linear presmoothing. \emph{Journal of Educational and Behavioral Statistics, 34}, 7--23.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{glm}}, \code{\link{loglin}}
}
\examples{
set.seed(2010)
x <- round(rnorm(1000, 100, 15))
xscale <- 50:150
xtab <- freqtab(x, scales = xscale)

# Adjust frequencies
plot(xtab, y = cbind(average = freqavg(xtab),
  bump = freqbump(xtab)))

# Smooth x preserving first 3 moments
xlog1 <- loglinear(xtab, degrees = 3,
  asfreqtab = FALSE)
plot(xtab, y = xlog1, legendtext = "degree = 3")

# Add "teeth" and "gaps" to x
# Smooth with formula interface
teeth <- c(.5, rep(c(1, 1, 1, 1, .5), 20))
xttab <- as.freqtab(cbind(xscale, c(xtab) * teeth))
xlog2 <- presmoothing(~ poly(total, 3, raw = TRUE),
  xttab, showWarnings = FALSE)

# Smooth xt using score functions that preserve 
# the teeth structure (also 3 moments)
teeth2 <- c(1, rep(c(0, 0, 0, 0, 1), 20))
xt.fun <- cbind(xscale, xscale^2, xscale^3)
xt.fun <- cbind(xt.fun, teeth2, xt.fun * teeth2)
xlog3 <- loglinear(xttab, xt.fun, showWarnings = FALSE)

# Plot to compare teeth versus no teeth
op <- par(no.readonly = TRUE)
par(mfrow = c(3, 1))
plot(xttab, main = "unsmoothed", ylim = c(0, 30))
plot(xlog2, main = "ignoring teeth", ylim = c(0, 30))
plot(xlog3, main = "preserving teeth", ylim = c(0, 30))
par(op)

# Bivariate example, preserving first 3 moments of total
# and anchor for x and y, and the covariance
# between anchor and total
# see equated scores in Wang (2009), Table 4
xvtab <- freqtab(KBneat$x, scales = list(0:36, 0:12))
yvtab <- freqtab(KBneat$y, scales = list(0:36, 0:12))
Y <- as.data.frame(yvtab)[, 1]
V <- as.data.frame(yvtab)[, 2]
scorefun <- cbind(Y, Y^2, Y^3, V, V^2, V^3, V*Y)
wang09 <- equate(xvtab, yvtab, type = "equip",
  method = "chained", smooth = "loglin",
  scorefun = scorefun)
wang09$concordance
}
\keyword{smooth}
\keyword{models}