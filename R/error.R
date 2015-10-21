#----------------------------------------------------------------
# Bias and standard error of equating

#----------------------------------------------------------------
# Generic bootstrapping method

bootstrap <- function(x, ...) UseMethod("bootstrap")

#----------------------------------------------------------------
# Default method

bootstrap.default <- function(x, y, ...) {

	if(!is.freqtab(x) | is.freqtab(y))
		stop("'x' and 'y' must be frequency tables")
	else do.call(bootstrap.freqtab, c(list(x = x, y = y),
		list(...)))
}
		
#----------------------------------------------------------------
# Method for equate class

bootstrap.equate <- function(x, xp = x$x, yp = x$y, ...) {
	
	dots <- list(...)
	if(is.character(xp))
		xp <- x[[xp]]
	if(is.character(yp))
		yp <- x[[yp]]
	rmnames <- c("x", "y", "yx", "concordance",
		"bootstraps", "coefficients", "synthstats",
		"xsynthetic", "ysynthetic", "xsmooth", "ysmooth",
		"points")
	args <- x[-pmatch(rmnames, names(x), nomatch = 0)]
	dots[pmatch(rmnames, names(dots), nomatch = 0)] <- NULL
	mi <- pmatch(names(dots), names(args), nomatch = 0)
	args[mi] <- dots[as.logical(mi)]
	dots <- dots[!as.logical(mi)]
	do.call(bootstrap.freqtab, c(list("x" = xp, "y" = yp),
		args, dots))
}

#----------------------------------------------------------------
# Method for freqtab class

bootstrap.freqtab <- function(x, y, xn = sum(x),
	yn = sum(y), reps = 100, crit, args,
	eqs = FALSE, ...) {
	
	dots <- list(...)[names(list(...) != "")]
	if(missing(args)) {
		args <- list(dots)
		neq <- 1
		args[[1]]["verbose"] <- FALSE
	}
	else {
		neq <- length(args)
		for(i in 1:neq) {
			args[[i]][names(dots)] <- dots
			args[[i]]["verbose"] <- FALSE
		}
	}
	if(missing(y)) {
		yn <- xn
		y <- NULL
		xs <- scales(x, 1)
		ys <- scales(x, 2)
		xd <- as.data.frame(as.data.frame(x)[x > 0, 1:2])
		xp <- x[x > 0]/sum(x)
		xni <- nrow(xd)
		eqmats <- lapply(rep(NA, neq), matrix,
			nrow = length(xs), ncol = reps)
		for(i in 1:reps) {
			xi <- sample.int(xni, xn, replace = TRUE, prob = xp)
			xtemp <- freqtab(xd[xi, ], scales = list(xs, ys))
			for(j in 1:neq)
				eqmats[[j]][, i] <- do.call("equate",
					c(list(x = xtemp), args[[j]]))
		}
	}
	else {
		nx <- margins(x)
		ny <- margins(y)
		xs <- scales(x, 1:nx)
		ys <- scales(y, 1:ny)
		xd <- as.data.frame(as.data.frame(x)[x > 0, 1:nx])
		yd <- as.data.frame(as.data.frame(y)[y > 0, 1:ny])
		xp <- x[x > 0]/sum(x)
		yp <- y[y > 0]/sum(y)
		xni <- nrow(xd)
		yni <- nrow(yd)
		eqmats <- lapply(rep(NA, neq), matrix,
			nrow = length(scales(x, 1)), ncol = reps)
		for(i in 1:reps) {
			xi <- sample.int(xni, xn, replace = TRUE, prob = xp)
			xtemp <- freqtab(xd[xi, ], scales = xs)
			yi <- sample.int(yni, yn, replace = TRUE, prob = yp)
			ytemp <- freqtab(yd[yi, ], scales = ys)
			for(j in 1:neq)
				eqmats[[j]][, i] <- do.call("equate",
					c(list(x = xtemp, y = ytemp), args[[j]]))
		}
	}
	names(eqmats) <- names(args)
	out <- list(x = x, y = y, reps = reps, xn = xn, yn = yn,
		args = args, mean = sapply(eqmats, apply, 1, mean),
		se = sapply(eqmats, apply, 1, sd))
	if(!missing(crit)) {
		out$bias <- sapply(eqmats, apply, 1, mean) - crit
		out$rmse <- sqrt(out$bias^2 + out$se^2)
	}
	if(neq == 1)
		out[-(1:6)] <- lapply(out[-(1:6)], c)
	if(eqs)
		out$eqs <- if(neq == 1) eqmats[[1]] else eqmats
	out <- as.bootstrap(out)

	return(out)
}

#----------------------------------------------------------------
# Assign bootstrap class

as.bootstrap <- function(x) {
	
	class(x) <- "bootstrap"
	return(x)
}

#----------------------------------------------------------------
# Test for bootstrap class

is.bootstrap <- function(x) {
	
	return(class(x)[1] == "bootstrap")
}

#----------------------------------------------------------------
# Print method

print.bootstrap <- function(x, ...) {
	
	nf <- length(x$args)
	cat("\nBootstrap Equating Error\n\n")
	cat("Design:", if(is.null(x$y)) "single group"
			else if(margins(x$x) == 1) "equivalent groups"
			else "nonequivalent groups", "\n\n")
	cat("Replications:", x$reps, "\n\n")
	cat("Sample Sizes: x =", paste(x$xn, "; y =", sep = ""),
		x$yn, "\n\n")
}

#----------------------------------------------------------------
# Summary method

summary.bootstrap <- function(object, weights,
	subset, ...) {

	if(missing(subset))
		subset <- 1:length(scales(object$x))
	if(missing(weights))
		weights <- c(margin(object$x))[subset]/
			sum(margin(object$x)[subset])
	out <- data.frame(se = apply(cbind(object$se)[subset, ],
		2, mean), w.se = apply(cbind(object$se)[subset, ] *
		weights, 2, mean))
	if(!is.null(object$bias)) {
		tempbias <- cbind(object$bias)[subset, ]
		out$bias <- apply(tempbias, 2, mean)
		out$a.bias <- apply(abs(tempbias), 2, mean)
		out$w.bias <- apply(tempbias * weights, 2, mean)
		out$wa.bias <- apply(abs(tempbias * weights), 2, mean)
		out$rmse <- apply(cbind(object$rmse)[subset, ],
			2, mean)
		out$w.rmse <- apply(cbind(object$rmse)[subset, ] *
			weights, 2, mean)
	}
	class(out) <- c("summary.bootstrap", "data.frame")

	return(out)
}

#----------------------------------------------------------------
# Plot method

plot.bootstrap <- function(x, add = FALSE, out = "mean",
	xpoints, ypoints, addident = TRUE, identy,
	identcol = 1, rescale = c(0, 1), xlab = "Total Score",
	ylab, col = rainbow(length(x$args)), pch, lty = 1,
	subset, morepars = NULL, addlegend = TRUE,
	legendtext, legendplace = "bottomright", ...) {
	
	if(missing(subset)) subset <- 1:length(x$args)
	x$args <- x$args[subset]
	nx <- length(subset)
	xscale <- scales(x$x)

	out <- match.arg(tolower(out),
		c("se", "bias", "mean", "rmse"))
	if(!out %in% names(x))
		stop(paste("'x' does not contain", out))
	y <- cbind(cbind(x[[out]])[, subset])
	y <- apply(y, 2, function(z) z*rescale[2] + rescale[1])
	
	if(missing(ylab))
		ylab <- switch(out, mean = "Mean Equated Score",
			se = "Standard Error", bias = "Bias",
			rmse = "RMSE")
	if(!is.null(morepars)) {
		nopars <- c("xlab", "ylab", "col", "lty", "pch")
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
		if(!missing(xpoints) & !missing(ypoints))
			do.call(points, c(list(x = xpoints,
				y = ypoints, col = "lightgray"), morepars))
	}
	
	if(addident) {
		if(missing(identy))
			identy <- switch(out, mean = xscale,
				rep(0, length(xscale)))
		lines(xscale, identy*rescale[2] + rescale[1],
			col = identcol)
	}

	lty <- rep(lty, length = nx)
	col <- rep(col, length = nx)
	for(i in 1:nx)
		lines(xscale, y[, i], col = col[i],
			lty = lty[i])
	if(!missing(pch)) {
		pch <- rep(pch, length = nx)
		for(i in 1:nx)
			points(xscale, y[, i], col = col[i],
				pch = pch[i])
	}

	if(addlegend) {
		if(missing(legendtext)) {
			legendtext <- c("identity", "mean", "linear",
				"general", "circle", "equip", "composite")
			legendtext <- lapply(x$args, function(z)
				legendtext[charmatch(substr(z$type, 1, 2),
				legendtext)])
			if(margins(x$x) == 2 & !is.null(x$y)) {
				methods <- c("nW", "chain", "b/H", "tucker",
					"levine", "fE")
				methods <- lapply(x$args, function(z)
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
# Linear standard error under random groups

lnse <- function(x, y) {

	nx <- sum(x[, 2])
	ny <- sum(y[, 2])
	mx <- mean(x)
	sdx <- sd.freqtab(x)
	vary <- cov.freqtab(y)
	skewterm <- skew.freqtab(x)/nx + skew.freqtab(y)/ny
	kurtterm <- (kurt.freqtab(x) - 1)/(4*nx) +
		(kurt.freqtab(y) - 1)/(4*ny)
	se <- vector()
	for(i in 1:nrow(x)) {
		xterm <- (x[i, 1] - mx)/sdx
		se[i] <- vary*(1/nx + 1/ny + skewterm*xterm +
			kurtterm*xterm^2)
	}

	return(se)
}

#----------------------------------------------------------------
# Equipercentile standard error under random groups

eqse <- function(q, g0, gm, xn, yn) {

	return(((1 - q)*q/(xn*g0^2) + (1/(yn*g0^2))*
		(gm - q^2 + ((q - gm)^2)/g0))^.5)
}

#----------------------------------------------------------------
