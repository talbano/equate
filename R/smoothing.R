#----------------------------------------------------------------
# Methods for presmoothing frequency distributions

#----------------------------------------------------------------
# Generic smoothing function

presmoothing <- function(x, ...) UseMethod("presmoothing")

#----------------------------------------------------------------
# Main presmoothing function

presmoothing.default <- function(x, smoothmethod = c("none",
	"average", "bump", "loglinear"), jmin,
	asfreqtab = TRUE, ...) {

	smoothmethod <- match.arg(smoothmethod)
	if(smoothmethod == "none")
		return(x)
	else if(smoothmethod == "average")
		return(freqavg(x, jmin = jmin,
			asfreqtab = asfreqtab))
	else if(smoothmethod == "bump")
		return(freqbump(x, jmin = jmin,
			asfreqtab = asfreqtab))
	else if(smoothmethod == "loglinear")
		return(loglinear(x, asfreqtab = asfreqtab, ...))
}

#----------------------------------------------------------------
# Formula method

presmoothing.formula <- function(x, data, ...) {
	
	formula <- terms(as.formula(x))
	if(attr(formula, "response") == 0) {
		formula <- terms(reformulate(attr(formula,
			"term.labels"), "count"))
	}
	if(attr(formula, "intercept") == 0) {
		attributes(formula)$intercept <- 1
		warning("an intercept has been added to the model")
	}

	scorefun <- as.data.frame(model.matrix(formula,
		data = as.data.frame(data))[, -1])
	loglinear(data, scorefun, ...)
}


#----------------------------------------------------------------
# Internal function for loglinear smoothing

loglinear <- function(x, scorefun, degrees = list(4, 2, 2), grid,
	rmimpossible = FALSE, asfreqtab = TRUE, models,
	stepup = !missing(models), compare = FALSE,
	verbose = FALSE, showWarnings = TRUE, ...) {

	# Powers in higher order interactions should never 
	# be larger than in lower - they will be ignored
		
	xd <- as.data.frame(x)
	nx <- ncol(xd) - 1
	
	if(rmimpossible & nx > 1) {
		keepi <- apply(xd, 1, function(y)
			all(y[-1] <= y[1]))
		xd <- xd[keepi, ]
	}
	else
		keepi <- rep(TRUE, nrow(xd))
	
	if(missing(scorefun)) {
		if(missing(grid)) {
			# Create a grid
			if(length(degrees) < nx) # must be at least 0 for higher orders
				degrees[(length(degrees) + 1):nx] <- 0
			degrees <- lapply(degrees, function(y)
				rep(y, nx)[1:nx])
			# Start grid without intercept
			grid <- cbind(expand.grid(lapply(degrees[[1]],
				function(y) 0:y))[-1, ])
			# Remove higher order interactions as necessary
			if(nx > 1) {
				for(i in 2:nx) {
					# Make sure higher orders don't contain larger powers
					# They're already excluded in grid
					#degrees[[i]] <- pmin(degrees[[i - 1]],
					#	degrees[[i]])
					rm1 <- apply(grid, 1, function(y)
						sum(y == 0) == (nx - i))
					rm2 <- apply(sapply(1:nx, function(j)
						grid[, j] > degrees[[i]][j]), 1, any)
					grid <- grid[!(rm1 & rm2), ]
				}
				# Sort grid by orders
				grid <- cbind(grid[order(apply(grid, 1, function(y)
					sum(y == 0)), decreasing = T), ])
				os <- factor(nx - apply(grid, 1, function(y)
					sum(y == 0)))
				grid <- do.call("rbind", by(grid, os,
					function(y) y[order(apply(y, 1, max)), ]))
			}
			else
				os <- factor(nx - apply(grid, 1, function(y)
					sum(y == 0)))
		}
		if(stepup | compare) {
			# Create model index
			models <- as.numeric(factor(paste(os,
				apply(grid, 1, max), sep = ".")))
			mnames <- unique(paste(os, apply(grid, 1, max),
				sep = "."))
		}
		scorefun <- NULL
		for(j in 1:nrow(grid)) {
			tempfun <- sapply(1:nx, function(k)
				xd[, k]^grid[j, k])
			scorefun <- cbind(scorefun,
				apply(tempfun, 1, prod))
		}
		colnames(scorefun) <- apply(grid, 1, paste,
			collapse = ".")
	}
	else {
		if(stepup | compare) {
			if(missing(models))
				models <- 1:ncol(scorefun)
			mnames <- unique(models)
		}
		scorefun <- scorefun[keepi, ]
	}
	
	if(nrow(scorefun) != nrow(xd))
		stop("'scorefun' must contain the same ",
			"number of rows as 'x'")
	scorefun <- data.frame(f = xd[, nx + 1],
		scorefun, check.names = FALSE)
	if(stepup | compare) {
		if(ncol(scorefun) < 3)
			stop(paste("cannot run multiple models with only",
				ncol(scorefun) - 1, "model terms"))
		snames <- colnames(scorefun)[-1]
		if(showWarnings)
			out <- lapply(unique(models), function(i)
				glm(scorefun[, c("f", snames[models <= i])],
					family = poisson))
		else
			suppressWarnings(out <- lapply(unique(models), function(i)
				glm(scorefun[, c("f", snames[models <= i])],
					family = poisson)))
		names(out) <- mnames
	}
	else if(showWarnings)
		out <- glm(scorefun, family = poisson)
	else
		suppressWarnings(out <- glm(scorefun, family = poisson))

	if(compare) {
		nm <- length(out)
		resdf <- as.numeric(lapply(out, function(y) y$df.residual))
		resdev <- as.numeric(lapply(out, function(y) y$deviance))
		aic <- as.numeric(lapply(out, AIC))
		bic <- as.numeric(lapply(out, BIC))
		tab <- data.frame(resdf, resdev, aic, bic,
			c(NA, -diff(resdf)), c(NA, 
			-diff(resdev)))
		vars <- lapply(out, function(y) paste(deparse(formula(y)), 
			collapse = "\n"))
		dimnames(tab) <- list(1:nm, c("Resid. Df", "Resid. Dev", 
			"AIC", "BIC", "Df", "Deviance"))
		tab <- stat.anova(tab, test = "Chisq", scale = 1, 
			df.scale = Inf,
			n = length(out[[order(resdf)[1]]]$residuals))
		return(structure(tab,
			heading = c("Analysis of Deviance Table\n",
			paste("Model ", format(1:nm), ": ", vars, 
			sep = "", collapse = "\n")), class = c("anova", 
			"data.frame")))
	}
	else if(verbose)
		return(out)
	else if(stepup)
		return(data.frame(lapply(out, fitted),
			check.names = FALSE))
	else if(asfreqtab)
		return(as.freqtab(cbind(xd[, 1:nx],
			out$fitted), scales = scales(x, 1:nx)))
	else
		return(out$fitted)
}

#----------------------------------------------------------------
# Frequency adjustment
# Bump frequencies upward by a small amount

freqbump <- function(x, jmin = 1e-6, asfreqtab = FALSE, ...) {

	x <- as.data.frame(x)
	nc <- ncol(x)
	f <- x[, nc]/sum(x[, nc])
	out <- (f + jmin)/(1 + (max(x[, 1]) + 1)*jmin)
	out <- ((f + jmin)/sum(f + jmin))*sum(x[, nc])
	
	if(asfreqtab)
		return(as.freqtab(cbind(x[, -nc], out)))
	else
		return(out)
}

#----------------------------------------------------------------
# Frequency averaging

freqavg <- function(x, jmin = 1, asfreqtab = FALSE, ...) {
	
	xtab <- x <- as.data.frame(x)
	nc <- ncol(xtab)
	if(nc > 2)
		stop("frequency averaging only supported ",
			"for univariate 'x'")
		
	x <- cbind(x, 0, 0, 0, 0, 0)
	ks <- 1
	while(ks <= nrow(x) & x[ks, 2] < jmin)
		ks <- ks + 1
	x[1:ks, 3] <- 1
	x[1:ks, 4] <- ks

	lls <- nrow(x)
	while(lls >= 0 & x[lls, 2] < jmin)
		lls <- lls - 1
	x[lls:nrow(x), 3] <- lls
	x[lls:nrow(x), 4] <- nrow(x)

	ss <- ks + 1
	tts <- lls - 1
	for(j in ss:tts) {
		if(x[j, 2] < jmin) {
			lls <- j
			ks <- j
			while(lls >= 1 & x[lls, 2] < jmin)
				lls <- lls - 1
			while(ks <= nrow(x) & x[ks, 2] < jmin)
				ks <- ks + 1
			x[lls:ks, 3] <- lls
			x[lls:ks, 4] <- ks
		}
	}

	for(p in 1:(nrow(x) - 1)) {
		if(x[p, 4] > 0 & x[p, 4] == x[p + 1, 3]) {
			if(x[p, 3] > 0)
				lls <- x[p, 3]
			if(x[p + 1, 4] > 0)
				ks <- x[p + 1, 4]
			x[lls:ks, 3] <- lls
			x[lls:ks, 4] <- ks
		}
	}

	for(j in 1:nrow(x)) {
		lls <- x[j, 3]
		if(lls == 0)
			lls <- j
		ks <- x[j, 4]
		if(ks == 0)
			ks <- j
		sumit <- 0
		sumit <- sumit + sum(x[lls:ks, 2])
		for(i in lls:ks) {
			x[i, 5] <- sumit
			x[i, 6] <- x[j, 4] - x[j, 3] + 1
			x[i, 7] <- x[i, 5]/x[i, 6]
		}
		j <- j + x[j, 4] - x[j, 3]
	}

	colnames(x)[c(1, 2, 6, 7)] <-
		c("score", "count", "b", "acount")

	if(asfreqtab)
		return(as.freqtab(cbind(xtab[, -nc],
			x[, 7])))
	else
		return(x[, 7])
}

#----------------------------------------------------------------
