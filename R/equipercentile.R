#----------------------------------------------------------------
# Equipercentile equating functions

#----------------------------------------------------------------
# The main function

equipercentile <- function(x, y, method = "none",
	ws = -1, smoothmethod = c("none", "bump",
	"average", "loglinear"), lowp = c(min(scales(x)),
	min(scales(y))), highp = c(max(scales(x)), max(scales(y))),
	verbose = FALSE, ...) {

	smoothmethod <- match.arg(smoothmethod)
	if (missing(y))
		y <- margin(x, 2)
	if (margins(y) < margins(x)) {
		xtab <- presmoothing(x, smoothmethod, ...)
		ytab <- margin(xtab, 2)
		xtab <- margin(xtab)
		x <- margin(x)
	} else {
		xtab <- presmoothing(x, smoothmethod, ...)
		ytab <- presmoothing(y, smoothmethod, ...)
	}

	if (method == "none")
		yxs <- equip(xtab, ytab, ly = lowp[2], ky = highp[2])
	else if (method == "frequency estimation") {
		stabs <- synthetic(xtab, ytab, ws, method)
		yxs <- equip(margin(stabs$xsynthetic),
			margin(stabs$ysynthetic), ly = lowp[2], ky = highp[2])
	} else if (method == "chained") {
		vx <- equip(margin(xtab),
			margin(xtab, 2))$yx
		pvyx <- px(vx, margin(ytab, 2))
		yxs <- equip(pvyx, margin(ytab))
	}
	yx <- yxs$yx
	
	if (!verbose)
		out <- yx
	else {
		out <- c(list(x = x, y = y, concordance =
				data.frame(scale = scales(x), yx = yx),
			points = data.frame(low = lowp, high = highp,
				row.names = c("x", "y")),
			smoothmethod = smoothmethod), list(...)[names(list(...))
				%in% c("jmin", "degree", "xdegree", "scorefun")])
		if (method == "frequency estimation")
			out <- c(out, stabs)
		if (smoothmethod != "none") {
			out$xsmooth <- xtab
			out$ysmooth <- ytab
		}
		if (method == "none" & is.design(x, "eg") & is.design(y, "eg"))
		  out$concordance$se <- yxs$se
	}
	return(out)
}

#----------------------------------------------------------------
# Equipercentile conversion of x to y
# Does this work for scales with negatives? Yes
# Does it work for shrunken or stretched scales? Yes, now.
# Previously, it is assumed that scores were in 1 point increments

equip <- function(x, y, ly = min(scales(y)), ky = max(scales(y))) {

	yscale <- scales(y)
	yn <- sum(y)
	if (!is.freqtab(x)) {
		prank <- sort(unique(x))
		xscale <- yscale
		xn <- 0
	} else {
		prank <- round(px(x), 10)
		xscale <- scales(x)
		xn <- sum(x)
		se <- rep(0, length = length(prank))
	}
	sn <- length(yscale)
	yinc <- round(diff(yscale), 8)
	yincl <- c(yinc[1]/2, yinc/2)
	yinch <- c(yinc/2, yinc[sn - 1]/2)
	yx <- numeric(length = length(prank))
	fy <- round(fx(y), 10)
	xnone <- prank == 0
	xone <- prank == 1
	xbot <- sum(xnone) + 1
	xtop <- sum(!xone)
	yxi <- xbot:xtop
	xyone <- which(xscale[xone] > (ky + yinch[sn])) + xtop

	yx[xnone] <- ly - yincl[1]
	yx[xone] <- ky + yinch[sn]
	yx[xyone] <- xscale[xyone]
	if (any(yx == 0)) {
		yu <- sapply(yxi, function(i) sum(fy <= prank[i]) + 1)
		yu2 <- yu - 1
		yu2[yu2 > 0] <- fy[yu2]
		g0 <- fy[yu] - yu2
		yx[yxi] <- yscale[yu] - yincl[yu] + ((prank[yxi] - yu2)/g0) *
		  (yinch + yincl)[yu]
		# standard errors
		if (xn)
			se[yxi] <- egese(prank[yxi], g0, yu2, xn, yn)
		if (any(y == 0)) {
			yxil <- (xbot + sum(prank[!xnone] <= min(fy))):xtop
			yl <- sapply(yxil, function(i) sum(fy < prank[i]))
			yl2 <- fy[yl + 1]
			yxtemp <- yscale[yl] + yincl[yl] +
				((prank[yxil] - fy[yl])/(yl2 - fy[yl])) *
				(yinch + yincl)[yl]
			yx[yxil] <- (yx[yxil] + yxtemp)/2
		}
	}

	if (xn) {
		out <- list(yx = yx)
		out$se <- se
	} else
		out <- list(yx = yx[match(x, prank)])

	return(out)
}

#----------------------------------------------------------------
# Standard errors

egese <- function(q, g0, gm, xn, yn) {
  return(((1 - q) * q/(xn * g0^2) + (1/(yn * g0^2)) *
  (gm - q^2 + ((q - gm)^2)/g0))^.5)
}

#----------------------------------------------------------------
