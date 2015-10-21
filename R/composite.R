#----------------------------------------------------------------
# Composite equating functions including
# a main function, conversion function, utility
# functions, summary, plot, and print methods

#----------------------------------------------------------------
# Generic composite function

composite <- function(x, ...) UseMethod("composite")

#----------------------------------------------------------------
# equate.list method

composite.equate.list <- function(x, wc, name, symmetric = FALSE,
	p = 1, verbose = TRUE, ...) {
	
	if(missing(wc))
		wc <- rep(1/length(x), length(x))
	if(symmetric) {
		if(!all(sapply(x, function(z) z$type) %in%
			c("identity", "mean", "linear", "general"))) {
				warning("all components must be linear to create ",
					"symmetric weights")
				wcs <- wc
				symmetric <- FALSE
			}
		else {
			slopes <- sapply(x, function(z) z$coef[2])
			wcs <- (wc*(1 + slopes^p)^-(1/p))/
					sum(wc*(1 + slopes^p)^-(1/p))
		}
	}
	else wcs <- wc
		
	yx <- composite.default(sapply(x, function(z) z$conc$yx),
		wcs)

	if(verbose) {
		if(missing(name))
			name <- paste("Composite Equating:",
				strsplit(x[[1]]$name, ": ")[[1]][2])
		out <- list(name = name, type = "composite",
			design = x[[1]]$design, x = x[[1]]$x, y = x[[1]]$y,
			concordance = data.frame(scale = x[[1]]$conc$scale,
				yx = yx), wc = wc, symmetric = symmetric)
		if(symmetric) out$wcs <- wcs
		out <- as.composite(out)
		out <- as.equate.list(c(list(out), x))
	}
	else out <- yx
	
	return(out)
}

#----------------------------------------------------------------
# Default method

composite.default <- function(x, wc, ...) {
	
	yx <- x %*% wc
	return(yx)
}

#----------------------------------------------------------------
# list method

composite.list <- function(x, wc, name,	symmetric = FALSE,
	p = 1, verbose = TRUE, ...) {
	
	if(!all(sapply(x, function(z) is.equate(z))))
		stop("all elements of 'x' must be class 'equate'")
	
	return(composite(as.equate.list(x), wc, name,
		symmetric, p, verbose, ...))
}

#----------------------------------------------------------------
# Convert to composite class

as.composite <- function(x) {
	
	class(x) <- c("composite", "equate")
	return(x)
}

#----------------------------------------------------------------
# Check for composite class

is.composite <- function(x) {
	
	return(class(x)[1] == "composite")
}

#----------------------------------------------------------------
# Print method

print.composite <- function(x, ...) {

	cat("\n")
	cat(x$name, "\n\n")
	cat("Design:", x$design, "\n\n")

	stats <- rbind(x = summary(margin(x$x)),
		y = summary(margin(x$y)),
		yx = summary(as.freqtab(cbind(x$concordance[, 2],
			c(margin(x$x))))))
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
