\name{plot.equate}
\alias{plot.equate}
\alias{plot.equate.list}
\title{Plotting Equating Results}
\description{
This function plots equating functions from one or more objects of class \dQuote{\code{equate}}.
}
\usage{
\method{plot}{equate}(..., elist = NULL, add = FALSE, out = "eqs",
  xpoints, ypoints, addident = TRUE, identy, identcol = 1,
  rescale = c(0, 1), xlab = "Total Score", ylab,
  col = rainbow(length(x)), pch, lty = 1, lwd = 1, subset,
  morepars = NULL, addlegend = TRUE, legendtext,
  legendplace = "bottomright")
}
\arguments{
  \item{\dots}{
one or more equating objects, each containing results for equating the same two test forms.
}
  \item{elist}{
optional list of additional equatings to be plotted.
}
  \item{add}{
logical, with default \code{FALSE}, specifying whether to create a new plot or add to the current one.
}
  \item{out}{
character vector specifying the output to be plotted, either equating functions (\code{"eqs"}), standard errors (\code{"se"}), bias (\code{"bias"}), or RMSE (\code{"rmse"}).
}
  \item{xpoints, ypoints}{
optional vectors of the same length containing raw scores on forms X and Y, assuming a single group or equivalent groups design.
}
  \item{addident}{
logical, with default \code{TRUE}, for plotting the identity function. The result depends on \code{out}.
}
  \item{identy}{
vector of y coordinates for plotting the identity line. Defaults to the X scale when \code{out = "eqs"}, otherwise, a horizontal line with intercept 0.
}
  \item{identcol}{
color used for plotting the identity line.
}
  \item{rescale}{
intercept and slope, with default 0 and 1, used to rescale all lines before plotting.
}
  \item{xlab, ylab, col, pch, lty, lwd}{
graphical parameters passed to \code{par}, with \code{col}, \code{pch}, \code{lty}, and \code{lwd} recycled as necessary.
}
  \item{subset}{
vector for subsetting the output when multiple equating functions are included in \code{x}.
}
  \item{morepars}{
list of additional graphical parameters, excluding \code{xlab}, \code{ylab}, \code{col}, \code{pch}, \code{lty}, and \code{lwd}.
}
  \item{addlegend}{
logical, with default \code{TRUE}, indicating whether or not a legend should be added.
}
  \item{legendtext}{
character vector of text to be passed to the \code{legend} argument of the \code{legend} function, defaulting to a combination of the equating types and methods specified in each equating object.
}
  \item{legendplace}{
placement of the legend.
}
}
\details{
Equating functions (\code{out = "eqs"}) are plotted as lines based on the concordance table for each equating object that is supplied. Standard errors (\code{out = "se"}) default to bootstrap standard errors, if available, otherwise, analyitical standard errors are plotted. Bias (\code{out = "bias"}) and RMSE (\code{out = "rmse"}) are also taken from bootstrapping output.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{equate}}, \code{\link{plot.bootstrap}}
}
\examples{
x <- freqtab(KBneat$x, scales = list(0:36, 0:12))
y <- freqtab(KBneat$y, scales = list(0:36, 0:12))

set.seed(111213)
eq1 <- equate(x, y, type = "linear", method = "tuck",
  boot = TRUE, reps = 10)
eq2 <- equate(x, y, type = "equip", method = "freq",
  boot = TRUE, smooth = "log", reps = 10,
  showWarnings = FALSE)
eq3 <- equate(x, y, type = "linear", method = "braun",
  boot = TRUE, reps = 10, showWarnings = FALSE)

plot(eq1, eq2, eq3)
plot(eq1, eq2, eq3, out = "se", legendplace = "top",
  addident = FALSE)
}
\keyword{misc}
