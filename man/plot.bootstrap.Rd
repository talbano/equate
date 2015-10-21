\name{plot.bootstrap}
\alias{plot.bootstrap}
\title{Plotting Bootstrap Equating Results}
\description{
This function plots bootstrap equating results for objects of class \dQuote{\code{bootstrap}}.
}
\usage{
\method{plot}{bootstrap}(x, add = FALSE, out = "mean", xpoints,
  ypoints, addident = TRUE, identy, identcol = 1,
  rescale = c(0, 1), xlab = "Total Score", ylab,
  col = rainbow(length(x$args)), pch, lty = 1, subset,
  morepars = NULL, addlegend = TRUE, legendtext,
  legendplace = "bottomright", ...)
}
\arguments{
  \item{x}{
output from the \code{\link{bootstrap}} function.
}
  \item{add}{
logical, with default \code{FALSE}, specifying whether to create a new plot or add to the current one.
}
  \item{out}{
character vector specifying the output to be plotted, either the mean equated scores (\code{"mean"}), standard errors (\code{"se"}), bias (\code{"bias"}), or RMSE (\code{"rmse"}).
}
  \item{xpoints, ypoints}{
optional vectors of the same length containing raw scores on forms X and Y, assuming a single group or equivalent groups design.
}
  \item{addident}{
logical, with default \code{TRUE}, for plotting the identity function. The result depends on \code{out}.
}
  \item{identy}{
vector of y coordinates for plotting the identity line. Defaults to the identity function when \code{out = "eqs"}, otherwise, a horizontal line with intercept 0.
}
  \item{identcol}{
color used for plotting the identity line.
}
  \item{rescale}{
intercept and slope, with default 0 and 1, used to rescale all lines before plotting.
}
  \item{xlab, ylab, col, pch, lty}{
graphical parameters passed to \code{par}, with the lengths of col and lty recycled as necessary.
}
  \item{subset}{
vector for subsetting the output when multiple equating functions are included in \code{x}.
}
  \item{morepars}{
list of additional graphical parameters, excluding \code{xlab}, \code{ylab}, \code{col}, \code{pch}, \code{lty}.
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
  \item{\dots}{
further arguments passed to or from other methods, excluding graphical parameters.
}
}
\details{
Lines are plotted for the chosen output type, whether mean equated scores across replications (\code{out = "mean"}), standard errors (\code{out = "se"}), bias (\code{out = "bias"}) or RMSE (\code{out = "rmse"}). The result is similar to that of \code{\link{plot.equate}}.
}
\author{
Anthony Albano \email{tony.d.albano@gmail.com}
}
\seealso{
\code{\link{bootstrap}}, \code{\link{plot.equate}}
}
\examples{
set.seed(122713)
neat.x <- freqtab(KBneat$x, scales = list(0:36, 0:12))
neat.y <- freqtab(KBneat$y, scales = list(0:36, 0:12))
eqargs <- list(m.t = list(type = "mean", method = "t"),
	l.t = list(type = "lin", method = "t"),
	c.t = list(type = "circ", method = "t"))
bootout <- bootstrap(x = neat.x, y = neat.y, args = eqargs,
	reps = 20)
plot(bootout, out = "se", legendplace = "top")
}
\keyword{misc}
