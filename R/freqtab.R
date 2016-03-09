#----------------------------------------------------------------
# Functions for creating, manipulating, summarizing,
# and plotting frequency tables

#----------------------------------------------------------------
# Generic freqtab function

freqtab <- function(x, ...) UseMethod("freqtab")

#----------------------------------------------------------------
# Default data.frame method
# For data.frame of item responses
# or data.frame of total and anchor scores
# or vector of total scores

freqtab.default <- function(x, scales, items,
	na.rm = TRUE, ...) {

	x <- data.frame(x)
	if(!missing(items)) {
		if(!is.list(items))
			items <- list(items)
		x <- sapply(items, function(i) {
			if(is.factor(i))
				i <- as.character(i)
			rowSums(cbind(x[, i]), na.rm = na.rm)
			})
	}

	nx <- ncol(x)
	x <- as.data.frame(x,
		row.names = seq_len(nrow(x)))
	if(missing(scales))
		scales <- apply(x, 2, function(y)
			min(y):max(y))
	if(!is.list(scales))
		scales <- list(scales)
	for(i in 1:nx)
		x[, i] <- factor(as.character(x[, i]),
			levels = scales[[i]])
	
	out <- freqtab.table(table(x))

	return(out)
}

#----------------------------------------------------------------
# Assign freqtab class and attributes
# based on scales

# Note drop will drop levels with zero counts across all
# other factors - you can't drop if x is a vector
# adding levels with zero counts is achieved by
# including those values in scales and setting drop
# to FALSE
# You can't drop current zeros and then add new ones

as.freqtab <- function(x, scales, drop = FALSE, ...) {
	
	x <- as.data.frame(x)
	nx <- ncol(x) - 1

	# If x is a vector of counts, scales must be supplied
	if(nx == 0) {
		dm <- sapply(scales, length)
		if(nrow(x) != prod(dm))
			stop("'scales' do not match dimensions of 'x'")
		out <- array(unlist(x), dm, dimnames = scales)
	}
	else {
		if(drop)
			x <- x[x[, nx + 1] > 0, ]
		# Set scales
		if(missing(scales)) {
			if(nx == 1)
				scales <- list(unique(x[, 1]))
			else if(nx > 1)
				scales <- lapply(x[, -(nx + 1)], function(y)
					as.character(sort(unique(y))))
			else
				stop("'scales' cannot be inferred from 'x'")
		}
		else if(!is.list(scales))
			scales <- list(scales)
		# Refactor like in xtabs
		xf <- lapply(1:nx, function(i) {
			factor(as.character(x[, i]),
				levels = scales[[i]])[, drop = drop]
		})
    	out <- tapply(x[, nx + 1], xf, sum)
	    out[is.na(out)] <- 0
	}

	return(freqtab.table(as.table(out)))
}

#----------------------------------------------------------------
# Test for freqtab class

is.freqtab <- function(x) {
	
	return(class(x)[1] == "freqtab")
}

#----------------------------------------------------------------
# table method

freqtab.table <- function(x, ...) {
	
	# Set attributes, mainly dnn
	out <- c(x)
	attributes(out) <- attributes(x)
	nx <- margins(x)
	dnn <- "total"
	if(nx == 2)
		dnn[2] <- "anchor"
	else if(nx > 2)
		dnn[2:nx] <- paste("anchor",
			1:(nx - 1), sep = "")
	names(dimnames(out)) <- dnn
	class(out) <- c("freqtab", "table")
	
	return(out)
}

#----------------------------------------------------------------
# Function for converting to data.frame
# All factors are converted to numeric

as.data.frame.freqtab <- function(x, row.names = NULL,
	optional = FALSE, drop = FALSE, ...) {
	
	out <- as.data.frame.table(x, row.names = NULL,
		responseName = "count", stringsAsFactors = TRUE)
	out <- sapply(out, function(y)
		as.numeric(as.character(y)))
	if(drop)
		out <- out[out[, ncol(out)] > 0, ]

	return(as.data.frame(out))
}

#----------------------------------------------------------------
# head and tail methods

head.freqtab <- function(x, ...) {
	
	head(as.data.frame(x), ...)
}

tail.freqtab <- function(x, ...) {
	
	tail(as.data.frame(x), ...)
}

#----------------------------------------------------------------
# print method

print.freqtab <- function(x, ...) {
	
	print(as.data.frame(x))
}

#----------------------------------------------------------------
# Function for extracting scales

scales <- function(x, margin = 1) {
	
	if(length(margin) == 1)
		return(as.numeric(dimnames(x)[[margin]]))
	else if(length(margin) > 1)
		return(lapply(dimnames(x)[margin], as.numeric))
}

#----------------------------------------------------------------
# Wrapper for margin.table

margin <- function(x, margin = 1) {
	
	if(any(!margin %in% seq(margins(x))))
		stop("misspecified margins")

	return(margin.table(x, margin))
}

#----------------------------------------------------------------
# Function for getting the number of margins

margins <- function(x) {
	
	return(length(dim(x)))
}

#----------------------------------------------------------------
# Function for dropping unused levels

droplevels.freqtab <- function(x, ...) {
	
	xd <- as.data.frame(x)[x > 0, ]
	return(as.freqtab(droplevels(xd)))
}

#----------------------------------------------------------------
# Plot method
# For plotting univariate and bivariate frequency distributions

plot.freqtab <- function(x, y = NULL, xcol = 1,
	ycol, pch = 16, ylty = 1, xlab = "Total Test",
	addlegend = !missing(y), legendtext, ...) {
	
	nx <- margins(x)
	if(nx == 1)
		ufreqplot(x, y, xcol, ycol, ylty, xlab,
			horiz = FALSE, addlegend = addlegend,
			legendtext = legendtext, ...)
	else if(nx == 2)
		bfreqplot(x, y, xcol, ycol, pch, ylty, xlab,
			addlegend = addlegend,
			legendtext = legendtext, ...)
	else stop("'x' must be either univariate or bivariate")
}

#----------------------------------------------------------------
# Internal univariate plot

ufreqplot <- function(x, y = NULL, xcol = 1, ycol,
	ylty = 1, xlab = "Total Test", ylab = "Count",
	horiz = FALSE, addlegend = FALSE,
	legendtext, ...) {

	x <- as.data.frame(x)
	if(!is.null(y)) {
		if(is.freqtab(y))
			y <- cbind(y[, 2])
		else
			y <- cbind(y)
		if(missing(ycol))
			ycol <- rainbow(ncol(y))
	}
	
	if(horiz) {
		plot.default(round(range(0, x[, 2], y)),
			range(x[, 1]), type = "n", xlab = xlab,
			ylab = ylab, ...)
		segments(rep(0, nrow(x)), x[, 1], x[, 2],
			col = xcol)
		if(!is.null(y))
			matlines(y, x[, 1], col = ycol,
				lty = ylty, ...)
	}
	else {
		plot.default(range(x[, 1]),
			round(range(0, x[, 2], y)), type = "n",
			xlab = xlab, ylab = ylab, ...)
		segments(x[, 1], y0 = rep(0, nrow(x)),
			y1 = x[, 2], col = xcol)
		if(!is.null(y))
			matlines(x[, 1], y, col = ycol,
				lty = ylty, ...)
	}
	
	if(addlegend & !is.null(y)) {
		if(missing(legendtext))
			legendtext <- if(is.null(colnames(y)))
				1:ncol(y) else colnames(y)
		legend("topright", legend = legendtext,
			lty = ylty, col = ycol, bty = "n")
	}
}

#----------------------------------------------------------------
# Internal bivariate plot

bfreqplot <- function(x, y = NULL, xcol = 1,
	ycol, pch = 16, ylty = 1, xlab = "Total Test",
	ylab = "Anchor Test", addlegend = FALSE,
	legendtext, ...) {

	xtab <- margin(x)
	xvtab <- margin(x, 2)
	xd <- as.data.frame(x)

	if(!is.null(y)) {
		if(is.freqtab(y))
			y <- cbind(c(y))
		if(missing(ycol))
			ycol <- rainbow(ncol(y))
		ytab <- apply(y, 2, function(z)
			tapply(z, xd[, 1], sum))
		yvtab <- apply(y, 2, function(z)
			tapply(z, xd[, 2], sum))
	}
	else ytab <- yvtab <- NULL
	
	reset.par <- par(no.readonly = TRUE)
	nf <- layout(matrix(c(2, 4, 1, 3), 2, 2,
		byrow = TRUE), c(3, 1), c(1, 3), TRUE)
	par(mar = c(4, 4, 1, 1))
	plot(range(xtab), range(xvtab), type = "n",
		xlab = xlab, ylab = ylab, ...)
	points(x, xcol = xcol, pch = pch)

	par(mar = c(0, 4, 1, 1))
	ufreqplot(xtab, ytab, xcol, ycol, ylty,
		xlab = "", ylab = "", xaxt = "n", bty = "n")

	par(mar = c(4, 0, 1, 1))
	ufreqplot(xvtab, yvtab, xcol, ycol, ylty,
		xlab = "", ylab = "", yaxt = "n", bty = "n",
		horiz = TRUE)

	if(addlegend & !is.null(y)) {
		par(mar = c(0, 0, 0, 0))
		plot(0, 0, type = "n", bty = "n", xaxt = "n",
			yaxt = "n")
		if(missing(legendtext))
			legendtext <- if(is.null(colnames(y)))
				1:ncol(y) else colnames(y)
		legend("bottomleft", legend = legendtext,
			lty = ylty, col = ycol, bty = "n")
	}
	
	par(reset.par)
}

#----------------------------------------------------------------
# Points method

points.freqtab <- function(x, xcol = 1, pch = 16,
	ds = 50, dm = 100, ...) {
	
	x <- as.data.frame(x)
	if(ncol(x) != 3)
		stop("'x' must be a bivariate frequency table")
		
	index <- as.logical(x[, 3])
	xpoints <- x[index, 1]
	vpoints <- x[index, 2]
	if(sd(x[index, 3]) > 0)
		dens <- pmax(0, pmin(255,
			scale(x[index, 3]) * ds + dm))
	else dens <- rep(150, sum(index))
	rgbcol <- col2rgb(xcol)
	ptcol <- rgb(rgbcol[1], rgbcol[2], rgbcol[3],
		dens, maxColorValue = 255)
	points(xpoints, vpoints, col = ptcol,
		pch = pch, ...)	
}

#----------------------------------------------------------------
# Summary method

summary.freqtab <- function(object,
	margin = seq(margins(object)), ...) {
	
	out <- NULL
	for(i in margin) {
		xm <- margin(object, i)
		out <- rbind(out, data.frame(
			mean = mean(xm),
			sd = sd.freqtab(xm),
			skew = skew.freqtab(xm),
			kurt = kurt.freqtab(xm),
			min = min(xm),
			max = max(xm),
			n = sum(xm)))
	}
	rownames(out) <- names(dimnames(object))[margin]

	return(out)
}

#----------------------------------------------------------------
# Mean

mean.freqtab <- function(x, margin = 1, ...) {

	inmars <- margin %in% seq(margins(x))
	if(any(!inmars)) {
		margin <- margin[inmars]
		warning("misspecified margins ",
			"have been removed.")
	}
	out <- sapply(margin, function(i) {
		xm <- margin(x, i)
		sum(xm*scales(xm)/sum(xm))})

	return(out)
}

#----------------------------------------------------------------
# Standard deviation

sd.freqtab <- function(x, margin = 1) {

	return(sqrt(var.freqtab(x, margin)))
}

#----------------------------------------------------------------
# Variance

var.freqtab <- function(x, margin = 1) {
	
	inmars <- margin %in% seq(margins(x))
	if(any(!inmars)) {
		margin <- margin[inmars]
		warning("misspecified margins ",
			"have been removed.")
	}
	n <- sum(x)
	out <- sapply(margin, function(i) {
		xm <- margin(x, i)
		xsc <- scales(xm)
		(sum(xsc*xsc*xm) - (sum(xm*xsc)^2)/n)/(n - 1)})

	return(out)
}

#----------------------------------------------------------------
# Covariance

cov.freqtab <- function(x, margin = seq(margins(x))) {

	inmars <- margin %in% seq(margins(x))
	if(any(!inmars)) {
		margin <- margin[inmars]
		warning("misspecified margins ",
			"have been removed.")
	}
	n <- sum(x)
	nx <- length(margin)
	out <- matrix(nrow = nx, ncol = nx)
	for(i in 1:nx) {
		out[i, i:nx] <- out[i:nx, i] <-
			sapply(margin[i:nx], function(j) {
				xd <- as.data.frame(margin(x,
					unique(c(i, j))))
				nc <- ncol(xd)
				sum((xd[, 1] - mean(x, i))*
					(xd[, nc - 1] - mean(x, j))*
					xd[, nc])/(n - 1)})
	}

	attr(out, "dim") <- c(nx, nx)
	attr(out, "dimnames") <- list(names(dimnames(x))[margin],
		names(dimnames(x))[margin])
	return(out)
}

#----------------------------------------------------------------
# Correlation

cor.freqtab <- function(x, margin = seq(margins(x))) {

	sds <- 1/sd.freqtab(x, margin)
	covs <- cov.freqtab(x, margin)
	out <- diag(sds) %*% covs %*% diag(sds)
	attributes(out) <- attributes(covs)

	return(out)
}

#----------------------------------------------------------------
# Minimum

min.freqtab <- function(x, margin = 1, ..., na.rm = FALSE) {

	inmars <- margin %in% seq(margins(x))
	if(any(!inmars)) {
		margin <- margin[inmars]
		warning("misspecified margins ",
			"have been removed.")
	}
	out <- sapply(margin, function(i) {
		xm <- margin(x, i)
		min(scales(xm)[as.logical(xm)])})

	return(out)
}

#----------------------------------------------------------------
# Maximum

max.freqtab <- function(x, margin = 1, ..., na.rm = FALSE) {

	inmars <- margin %in% seq(margins(x))
	if(any(!inmars)) {
		margin <- margin[inmars]
		warning("misspecified margins ",
			"have been removed.")
	}
	out <- sapply(margin, function(i) {
		xm <- margin(x, i)
		max(scales(xm)[as.logical(xm)])})

	return(out)
}

#----------------------------------------------------------------
# Range

range.freqtab <- function(x, margin = 1, ...,
	na.rm = FALSE) {

	if(length(margin) == 1)
		out <- c(min(x), max(x))
	else if(length(margin) > 1)
		out <- lapply(margin, function(i)
			c(min(margin(x, i)), max(margin(x, i))))

	return(out)
}

#----------------------------------------------------------------
# Skewness

skew.freqtab <- function(x, margin = 1) {

	inmars <- margin %in% seq(margins(x))
	if(any(!inmars)) {
		margin <- margin[inmars]
		warning("misspecified margins ",
			"have been removed.")
	}
	n <- sum(x)
	out <- sapply(margin, function(i) {
		xm <- margin(x, i)
		xsc <- scales(xm)
		sum(((xsc - mean(xm))^3*xm))/(n)/
			(sum(((xsc - mean(xm))^2*xm))/(n - 1))^1.5})

	return(out)
}

#----------------------------------------------------------------
# Kurtosis

kurt.freqtab <- function(x, margin = 1) {

	inmars <- margin %in% seq(margins(x))
	if(any(!inmars)) {
		margin <- margin[inmars]
		warning("misspecified margins ",
			"have been removed.")
	}
	n <- sum(x)
	out <- sapply(margin, function(i) {
		xm <- margin(x, i)
		xsc <- scales(xm)
		sum(((xsc - mean(xm))^4*xm))/(n)/
			(sum(((xsc - mean(xm))^2*xm))/(n - 1))^2 - 3})

	return(out)
}

#----------------------------------------------------------------
# Cumulative frequency

fx <- function(x, ...) UseMethod("fx")

fx.default <- function(x, ...) {
	
	return(as.numeric(cumsum(x/sum(x))))
}

fx.freqtab <- function(x, margin = 1, ...) {
	
	if(!margin %in% seq(margins(x)))
		stop("'margin' not found in 'x'")
	
	return(fx.default(margin(x, margin)))
}

#----------------------------------------------------------------
# Percentile ranks with linear interpolation

px <- function(x, ...) UseMethod("px")

px.default <- function(x, y, ys, ...) {

	if(missing(y)) {
		x <- as.numeric(x/sum(x))
		p <- .5 * x[1]
		for(i in 2:length(x))
			p[i] <- sum(x[1:(i - 1)]) + .5 * x[i]
	}
	else {
		y <- as.data.frame(y)
		if(ncol(y) == 2) {
			ys <- y[, 1]
			y <- y[, 2]
		}
		xs <- floor(x + .5)
		yn <- sum(y)
		f <- sapply(xs, function(xi)
			sum(y[ys <= xi])/yn)
		flow <- sapply(xs, function(xi)
			sum(y[ys <= xi - 1])/yn)
		p <- flow + (x - (xs - .5)) * (f - flow)
	}
	return(p)
}

px.freqtab <- function(x, margin = 1,
	y, ymargin = 1, ...) {

	if(!margin %in% seq(margins(x)))
		stop("'margin' not found in 'x'")
	if(missing(y))
		p <- px.default(margin(x, margin))
	else {
		if(!ymargin %in% seq(margins(y)))
			stop("'ymargin' not found in 'y'")
		p <- px.default(scales(margin(x, margin)),
			as.data.frame(margin(y, ymargin)))
	}

	return(p)
}

#----------------------------------------------------------------
