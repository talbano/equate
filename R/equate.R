#----------------------------------------------------------------
# Equating functions for the equate class including
# a main function, conversion function, utility
# functions, summary, plot, and print methods

#----------------------------------------------------------------
# Generic equating function

equate <- function(x, ...) UseMethod("equate")

#----------------------------------------------------------------
# Equating a list of two or more data sets with one
# or more equatings

equate.list <- function(x, args, ...) {
	
	if(length(x) < 2)
		stop("'x' must contain at least 2 data sets")
	if(length(x) == 2 & (is.null(names(x)) |
		any(names(x) == "")))
		names(x) <- c("x", "y")
	if(is.null(names(x)) | any(names(x) == ""))
		stop("all elements of 'x' must be named")
	if(!is.list(args))
		stop("'args' must be a list")
	if(!all(sapply(args, is.list)))
		stop("each element in 'args' must be a list")
	if(length(x) > 2) {
		if(!all(unlist(lapply(args, function(z)
			"x" %in% names(z)))))
			stop("'args' must include 'x' for each equating")
	}
	dots <- list(...)
	neq <- length(args)
	
	out <- vector("list", length = neq)
	for(i in 1:neq) {
		eqargs <- args[[i]]
		eqargs[names(dots)] <- dots
		if(is.null(eqargs$x))
			eqargs$x <- names(x)[1]
		if(is.null(eqargs$y))
			eqargs$y <- names(x)[2]
		eqargs$x <- x[[eqargs$x]]
		eqargs$y <- x[[eqargs$y]]
		out[[i]] <- do.call(equate.freqtab, eqargs)
	}
	names(out) <- names(args)

	return(as.equate.list(out))
}

#----------------------------------------------------------------
# Assign equate.list class

as.equate.list <- function(x) {
	
	class(x) <- "equate.list"
	return(x)
}

#----------------------------------------------------------------
# Test for equate.list class

is.equate.list <- function(x) {
	
	return(class(x)[1] == "equate.list")
}

#----------------------------------------------------------------
# Extract method for equate.list

assign("[.equate.list", function(x, i){

	out <- as.equate.list(NextMethod("["))
	return(out)
})

#----------------------------------------------------------------
# Equating a frequency table

equate.freqtab <- function(x, y, type = c("identity",
	"mean", "linear", "general linear", "circle-arc",
	"equipercentile"), method = c("none", "nominal weights",
	"tucker", "levine", "frequency estimation", "chained",
	"braun/holland"), name, lowp, highp, boot = FALSE,
	verbose = TRUE, ...) {
	
	if(missing(y)) 
        y <- margin(x, 2)
	if(missing(lowp))
		lowp <- c(min(scales(x)), min(scales(y)))
	else if(pmatch(lowp[1], "obs", 0))
		lowp <- c(min(x), min(y))
	else if(length(lowp) == 1)
		lowp[2] <- lowp[1]
	if(missing(highp))
		highp <- c(max(scales(x)), max(scales(y)))
	else if(pmatch(highp[1], "obs", 0))
		highp <- c(max(x), max(y))
	else if(length(highp) == 1)
		highp[2] <- highp[1]
	
	type <- match.arg(type)
	method <- match.arg(method)
	if(type %in% c("identity", "mean", "linear", "general linear"))
		eqout <- linear(x, y, type = type, method = method,
			lowp = lowp, highp = highp, verbose = verbose, ...)
	else if(type == "equipercentile")
		eqout <- equipercentile(x, y, type = type, method = method,
			lowp = lowp, highp = highp, verbose = verbose, ...)
	else if(type == "circle-arc")
		eqout <- circlearc(x, y, type = type, method = method,
			lowp = lowp, highp = highp, verbose = verbose, ...)

	if(verbose) {
		if(missing(name)) {
			name <- ifelse(method == "none", type,
				paste(method, type))
			xname <- ifelse(exists(deparse(substitute(x)), 1,
				inherits = FALSE), deparse(substitute(x)), "x")
			if(margins(y) == margins(x))
				yname <- ifelse(exists(deparse(substitute(y)), 1,
					inherits = FALSE), deparse(substitute(y)), "y")
			else yname <- xname
			name <- paste(gsub("\\b(\\w)", "\\U\\1", name,
				perl = TRUE), "Equating:", xname, "to", yname)
		}
		design <- if(margins(y) < margins(x)) "single group"
			else if(margins(x) == 1) "equivalent groups"
			else "nonequivalent groups"
		out <- c(list(name = name, type = type,
			method = method, design = design), eqout)
		out <- as.equate(out)
		if(boot)
			out$bootstraps <- bootstrap.equate(out, ...)
	}
	else out <- eqout

	return(out)
}

#----------------------------------------------------------------
# Equating a vector of scores

equate.default <- function(x, y, ...) {

	if(is.equate(y))
		yx <- convert(x, y, ...)
	else if(is.equate.list(y)) {
		if(is.composite(y[[1]])) {
			yx <- sapply(y[-1], function(z) convert(x, z))
			wcs <- if(!is.null(y[[1]]$wcs)) y[[1]]$wcs
				else y[[1]]$wc
			yx <- c(yx %*% wcs)
		}
		else yx <- sapply(y, function(z) convert(x, z))
	}
	else
		stop("'y' must be an 'equate' or 'equate.list' object")
		
	return(yx)
}

#----------------------------------------------------------------
# Score conversion function

convert <- function(x, y, ...) {
	
	if(!is.equate(y))
		stop("'y' must be an 'equate' object")
		
	if(y$type == "equipercentile") {
		if(y$method == "none") {
			if(y$smoothmethod == "none") {
				xtab <- y$x
				ytab <- y$y
			}
			else {
				xtab <- y$xsmooth
				ytab <- y$ysmooth
			}
			p <- px(x, xtab)
		}
		else if(y$method == "frequency estimation") {
			xtab <- margin(y$xsynthetic)
			ytab <- margin(y$ysynthetic)
			p <- px(x, xtab)
		}
		else {
			if(y$smoothmethod == "none") {
				xtab <- margin(y$x)
				xvtab <- margin(y$x, 2)
				ytab <- margin(y$y)
				yvtab <- margin(y$y, 2)
			}
			else {
				xtab <- margin(y$xsmooth)
				xvtab <- margin(y$xsmooth, 2)
				ytab <- margin(y$ysmooth)
				yvtab <- margin(y$ysmooth, 2)
			}
			vx <- equip(px(x, xtab), xvtab)$yx
			p <- px(vx, yvtab)
		}
		yx <- equip(p, ytab)$yx
	}
	else if(y$type == "circle-arc") {
		yx <- circ(x, y$points[1, ], y$points[2, ],
			y$points[3, ], y$coef[1], y$coef[2], y$coef[3:4],
			y$coef[5], simple = y$simple)
	}
	else {
		yx <- lin(x, y$coef[1], y$coef[2])
	}
	names(yx) <- NULL

	return(yx)
}

#----------------------------------------------------------------
# Assign equate class

as.equate <- function(x) {
	
	class(x) <- "equate"
	return(x)
}

#----------------------------------------------------------------
# Test for equate class

is.equate <- function(x) {
	
	return(class(x)[1] == "equate")
}

#----------------------------------------------------------------
# Summary method

summary.equate <- function(object, ...) {

	xtab <- as.data.frame(margin(object$x))
	ytab <- as.data.frame(margin(object$y))
	if(object$type == "equipercentile" &&
		object$smoothmethod != "none") {
		xtab <- data.frame(xtab,
			smooth = c(margin(object$xsmooth)))
		ytab <- data.frame(ytab,
			smooth = c(margin(object$ysmooth)))
	}
	xvtab <- yvtab <- NULL
	if(object$design == "nonequivalent groups") {
		xvtab <- as.data.frame(margin(object$x, 2))
		yvtab <- as.data.frame(margin(object$y, 2))
		if(object$type == "equipercentile" &&
			object$smoothmethod != "none") {
			xvtab <- data.frame(xvtab,
				smooth = c(margin(object$xsmooth, 2)))
			yvtab <- data.frame(yvtab,
				smooth = c(margin(object$ysmooth, 2)))
		}
		if(object$type != "composite" && (object$method ==
			"frequency estimation" | object$method ==
			"braun/holland")) {
			xtab <- data.frame(xtab,
				synthetic = c(margin(object$xsynthetic)))
			ytab <- data.frame(ytab,
				synthetic = c(margin(object$ysynthetic)))
			xvtab <- data.frame(xvtab,
				synthetic = c(margin(object$xsynthetic, 2)))
			yvtab <- data.frame(yvtab,
				synthetic = c(margin(object$ysynthetic, 2)))
		}
	}
	yxtab <- as.freqtab(cbind(object$concordance[, 2], xtab$count))
	yxstats <- summary(yxtab)
	rownames(yxstats) <- "observed"
	
	xstats <- ystats <- NULL
	for(i in 2:ncol(xtab)) {
		xstats <- rbind(xstats,
			summary(as.freqtab(xtab[, c(1, i)])))
		ystats <- rbind(ystats,
			summary(as.freqtab(ytab[, c(1, i)])))
	}
	rownames(xstats) <- rownames(ystats) <-
		colnames(xtab)[-1]
	xvstats <- yvstats <- NULL
	if(object$design == "nonequivalent groups") {
		for(i in 2:ncol(xvtab)) {
			xvstats <- rbind(xvstats,
				summary(as.freqtab(xvtab[, c(1, i)])))
			yvstats <- rbind(yvstats,
				summary(as.freqtab(yvtab[, c(1, i)])))
		}
		rownames(xvstats) <- rownames(yvstats) <-
			colnames(xvtab)[-1]
	}

	out <- c(object[names(object) %in% c("name", "type", "method",
		"design", "ident", "ws", "smoothmethod")],
		freqtab = list(list(x = xtab, y = ytab,
			yx = yxtab)),
		stats = list(list(x = xstats, y = ystats,
			yx = yxstats)))
	if(!is.null(object$bootstraps)) {
		xw <- xtab[, 2]/sum(xtab[, 2])
		out$error <- data.frame(
			bias = mean(object$bootstraps$bias),
			abs.bias = mean(abs(object$bootstraps$bias)),
			w.bias = mean(object$bootstraps$bias*xw),
			wabs.bias = mean(abs(object$bootstraps$bias*xw)),
			se = mean(object$bootstraps$se),
			w.se = mean(object$bootstraps$se*xw),
			rmse = mean(object$bootstraps$rmse))
	}
	if(object$design == "nonequivalent groups") {
		out$freqtab <- c(out$freqtab, list(xv = xvtab,
			yv = yvtab))
		out$stats <- c(out$stats, list(xv = xvstats,
			yv = yvstats))
	}	
	class(out) <- "summary.equate"

	return(out)
}

#----------------------------------------------------------------
# Summary method for equate.list

summary.equate.list <- function(object, ...) {

	lapply(object, summary.equate)
}

#----------------------------------------------------------------
# Summary print method

print.summary.equate <- function(x, ...) {

	cat("\n")
	cat(x$name, "\n\n")
	cat("Design:", x$design, "\n\n")
	if(x$type == "equipercentile") {
		cat("Smoothing Method: ")
		cat(switch(x$smoothmethod,
			bump = "adjusted frequency presmoothing",
			average = "average frequency presmoothing",
			loglinear = "loglinear presmoothing", "none"),
			"\n\n")
	}
	if(x$design == "nonequivalent groups" &&
		x$type != "composite" &&
		x$method != "chained")
		cat("Synthetic Weighting for x:",
			x$ws, "\n\n")
	lstats <- lapply(x$stats, function(y)
		t(round(y, 3)))
	out <- data.frame(matrix(unlist(lstats),
		ncol = ncol(x$stats[[1]]), byrow = TRUE))
	colnames(out) <- colnames(x$stats[[1]])
	scales <- names(unlist(lapply(lstats, colnames)))
	scales <- gsub("1|2|3|4", "", scales)
	dists <- unlist(lapply(lstats, colnames))
	dists <- gsub("(etic)|(erved)", "", dists)
	rownames(out) <- paste(scales, dists, sep = ".")
	cat("Summary Statistics:\n")
	print(out)
	cat("\n")
	
	if(!is.null(x$error)) {
		cat("Mean Error:\n")
		print(round(x$error, 4), row.names = FALSE)
		cat("\n")
	}
	
	invisible(x)
}

#----------------------------------------------------------------
# Plot method

plot.equate <- function(..., elist = NULL, add = FALSE,
	out = "eqs", xpoints, ypoints, addident = TRUE,
	identy, identcol = 1, rescale = c(0, 1),
	xlab = "Total Score", ylab, col = rainbow(length(x)),
	pch, lty = 1, lwd = 1, subset, morepars = NULL, addlegend = TRUE,
	legendtext, legendplace = "bottomright") {
	
	x <- c(list(...), elist)
	if(missing(subset)) subset <- 1:length(x)
	x <- x[subset]
	nx <- length(x)
	xscale <- scales(x[[1]]$x)

	out <- match.arg(tolower(out),
		c("se", "bias", "eqs", "rmse"))
	if(out == "se") {
		y <- lapply(x, function(z) {
				if(is.null(z$bootstraps)) z$se
				else z$bootstraps$se
			})
	}
	else if(out == "bias")
		y <- lapply(x, function(z) z$bootstraps$bias)
	else if(out == "rmse")
		y <- lapply(x, function(z) z$bootstraps$rmse)
	else if(out == "eqs")
		y <- lapply(x, function(z) z$concordance[, 2])
	else
		stop("'out' must be one of 'eqs', 'se', 'bias' ",
			"or rmse")
	if(any(unlist(lapply(y, is.null))))
		stop("one or more equatings does not contain ", out)
	y <- lapply(y, function(z) z*rescale[2] + rescale[1])
	
	if(missing(ylab))
		ylab <- switch(out, eqs = "Equated Score",
			se = "Standard Error", bias = "Bias",
			rmse = "RMSE")
	if(!is.null(morepars)) {
		nopars <- c("xlab", "ylab", "col", "lty", "pch", "lwd")
		noparsl <- nopars %in% names(morepars)
		if(any(noparsl)) {
			warning("the following graphical parameter(s)",
				" must be specified outside of 'morepars': ",
				paste(nopars[noparsl], collapse = ", "))
			morepars <- morepars[!names(morepars) %in%
				nopars]
		}
	}
	
	if(!add) {
		do.call(plot, c(list(x = range(xscale),
			y = range(y), xlab = xlab, ylab = ylab,
			type = "n"), morepars))
		if(!missing(xpoints) && is.freqtab(xpoints))
			do.call(points.freqtab, c(list(x = xpoints,
				xcol = "lightgray"), morepars))
		else if(!missing(xpoints) & !missing(ypoints))
			do.call(points, c(list(x = xpoints,
				y = ypoints, col = "lightgray"), morepars))
	}
	
	if(addident) {
		if(missing(identy))
			identy <- switch(out, eqs = xscale,
				rep(0, length(xscale)))
		lines(xscale, identy*rescale[2] + rescale[1],
			col = identcol)
	}
	
	col <- rep(col, length = nx)
	lty <- rep(lty, length = nx)
	lwd <- rep(lwd, length = nx)
	for(i in 1:nx)
		lines(xscale, y[[i]], col = col[i],
			lty = lty[i], lwd = lwd[i])
	if(!missing(pch)) {
		pch <- rep(pch, length = nx)
		for(i in 1:nx)
			points(xscale, y[[i]], col = col[i],
				pch = pch[i])
	}

	if(addlegend) {
		if(missing(legendtext)) {
			legendtext <- c("identity", "mean", "linear",
				"general", "circle", "equip", "composite")
			legendtext <- lapply(x, function(z)
				legendtext[charmatch(substr(z$type, 1, 2),
				legendtext)])
			if(x[[1]]$design == "nonequivalent groups") {
				methods <- c("nW", "chain", "b/H", "tucker",
					"levine", "fE")
				methods <- lapply(x, function(z)
					methods[charmatch(substr(z$method, 1, 1),
					methods)])
				legendtext <- paste(legendtext,
					methods, sep = ": ")
				legendtext[grep("ident", legendtext)] <-
					"identity"
				legendtext[grep("comp", legendtext)] <-
					"composite"
			}
			legendtext <- gsub("\\b(\\w)", "\\U\\1",
				legendtext, perl = TRUE)
		}
		if(addident) {
			legendtext <- c("Identity", legendtext)
			lty = c(1, lty)
			col = c(identcol, col)
		}
		legend(legendplace, legend = legendtext,
			lty = lty, col = col, bty = "n")
	}
}

#----------------------------------------------------------------
# Plot method for equate.list

plot.equate.list <- function(x, ...) {

	plot.equate(elist = x, ...)
}

#----------------------------------------------------------------
# Print method

print.equate <- function(x, ...) {

	cat("\n")
	cat(x$name, "\n\n")
	cat("Design:", x$design, "\n\n")
	if(x$type == "equipercentile") {
		cat("Smoothing Method: ")
		cat(switch(x$smoothmethod,
			bump = "adjusted frequency presmoothing",
			average = "average frequency presmoothing",
			loglinear = "loglinear presmoothing", "none"),
			"\n\n")
	}
	if(!is.null(x$ws))
		cat("Synthetic Weighting for x:",
			x$ws, "\n\n")
	xm <- margin(x$x)
	stats <- rbind(x = summary(xm),
		y = summary(margin(x$y)),
		yx = summary(as.freqtab(cbind(x$concordance[, 2],
			c(xm)))))
	cat("Summary Statistics:\n")
		print(round(stats, 2))
		cat("\n")

	if(!is.null(x$coef)) {
		cat("Coefficients:\n")
		print(round(x$coef, 4))
		cat("\n")
	}

	invisible(x)
}

#----------------------------------------------------------------