#----------------------------------------------------------------
# Linear equating functions

#----------------------------------------------------------------
# The main function

linear <- function(x, y, type = "linear", method = "none",
	ws = -1, wax, way, wbx, wby, lowp = c(min(scales(x)),
	min(scales(y))), midp = c(median(c(lowp[1], highp[1])),
	median(c(lowp[2], highp[2]))), highp = c(max(scales(x)),
	max(scales(y))), sx = highp[1] - lowp[1], sy = highp[2] -
	lowp[2], cx = midp[1], cy = midp[2], internal = TRUE,
	lts = FALSE, verbose = FALSE, ...) {

	if (missing(y)) {
	  y <- margin(x, 2)
	  xsg <- x
	}
  else xsg <- NULL
	if (margins(y) < margins(x))
		x <- margin(x, 1)
	xscale <- scales(x)
	
	if (method == "chained") {
		if (type == "mean") {
			slope1 <- slope2 <- 1
		} else {
			slope1 <- sd.freqtab(y)/sd.freqtab(y, 2)
			slope2 <- sd.freqtab(x, 2)/sd.freqtab(x)
		}
		intercept <- slope1 * (mean(x, 2) -
			slope2 * mean(x) - mean(y, 2)) + mean(y)
		slope <- slope1 * slope2
		yx <- lin(xscale, intercept, slope)
	} else {
		if (method == "none") {
		  if (is.design(x, "cb") & is.design(y, "cb")) {
		    sigmax <- sqrt((var.freqtab(x, 1) + var.freqtab(y, 1))/2)
		    sigmay <- sqrt((var.freqtab(x, 2) + var.freqtab(y, 2))/2)
		    mux <- (mean(x, 1) + mean(y, 1))/2
		    muy <- (mean(x, 2) + mean(y, 2))/2
		  } else {
			  sigmax <- sd.freqtab(x)
			  sigmay <- sd.freqtab(y)
			  mux <- mean(x)
			  muy <- mean(y)
		  }
		} else {
			synth <- synthetic(x, y, ws, method, internal, lts)
			stats <- synth$synthstats
			if (!lts) {
				sigmax <- stats$sd[1]
				sigmay <- stats$sd[2]
				mux <- stats$m[1]
				muy <- stats$m[2]
			} else {
				sigmax <- synth$gamma[1]
				sigmay <- synth$gamma[2]
				mux <- mean(x)
				muy <- mean(y) + synth$gamma[2] *
					(mean(x, 2) - mean(y, 2))
			}
		}
		if (type %in% c("identity", "general linear")) {
			if (missing(wax)) wax <- 0
			if (missing(way)) way <- 0
			if (missing(wbx)) wbx <- 0
			if (missing(wby)) wby <- 0
		} else if (type == "mean") {
			if (missing(wax)) wax <- 0
			if (missing(way)) way <- 0
			if (missing(wbx)) wbx <- 1
			if (missing(wby)) wby <- 1
		} else if (type == "linear") {
			if (missing(wax)) wax <- 1
			if (missing(way)) way <- 1
			if (missing(wbx)) wbx <- 1
			if (missing(wby)) wby <- 1
		}
		if (cx == "cy") cx <- muy
		if (cy == "cx") cy <- mux
		if (sx == "sy") sx <- sigmay
		if (sy == "sx") sy <- sigmax
		# Put sigma on the scale of s
		# sigmayr <- sigmay*sx/sigmax
		# Then, sigmaxr <- sx
		slope <- (way * sigmay + (1 - way) * sy)/
			(wax * sigmax + (1 - wax) * sx)
		if (is.na(slope)) slope <- 1
		intercept <- (wby * muy + (1 - wby) * cy) -
			slope * (wbx * mux + (1 - wbx) * cx)
		yx <- lin(xscale, intercept, slope)
	}

	if (verbose) {
		out <- list(x = x, y = y, concordance = data.frame(scale =
			xscale, yx = yx), internal = internal, lts = lts,
			coefficients = c(intercept = intercept, slope = slope),
			points = data.frame(low = lowp, mid = midp,
				high = highp, row.names = c("x", "y")))
		if (method != "chained") {
			out$coefficients[c("cx", "cy", "sx", "sy")] <-
				c(cx, cy, sx, sy)
			out$weights <- c(wax = wax, way = way,
				wbx = wbx, wby = wby)
		}
		if (!method %in% c("none", "chained") & !lts)
			out <- c(out, synth)
		if(!lts)
		  out$concordance$se <- lse(out, type, method, xsg)
	}
	else out <- yx

	return(out)
}

#----------------------------------------------------------------
# Linear conversion of x to y

lin <- function(x, intercept, slope) {

	out <- x * slope + intercept
	names(out) <- NULL
	
	return(out)
}

#----------------------------------------------------------------
# Standard errors
# The functions dlin, ds, dcl, dp, and omeg are based on code
# by Zu (2012), copyright Educational Testing Service (www.ets.org)

lse <- function(x, type, method, xsg) {
  if (margins(x$x) > 2 | margins(x$y) > 2) {
    out <- NULL
  } else if (type == "identity") {
    out <- rep(0, length(scales(x$x)))
  } else if (type == "mean") {
    if (is.design(x$x, "sg"))
      out <- sgmse(x$x, x$y, xsg)
    else if (is.design(x$x, "eg") && is.design(x$y, "eg"))
      out <- egmse(x$x, x$y)
    else out <- NULL
  } else if (type == "linear") {
    if (method == "none") {
      if (is.design(x$x, "sg"))
        out <- sglse(x$x, x$y, xsg)
      else if (is.design(x$x, "eg") && is.design(x$y, "eg"))
        out <- eglse(x$x, x$y)
      else if (is.design(x$x, "cb") && is.design(x$y, "cb")) {
        out <- cblse(x$x, x$y)
      } else out <- rep(0, length(scales(x$x)))
    } else if (method %in% c("chained", "tucker", "levine")) {
      ox <- omeg(x$x)
      oy <- omeg(x$y)
      omegan <- omegag <- matrix(0, 10, 10)
      omegan[1:5, 1:5] <- ox$n
      omegan[6:10, 6:10] <- oy$n
      omegag[1:5, 1:5] <- ox$g
      omegag[6:10, 6:10] <- oy$g
    
      mx <- mean(x$x)
      mxv <- mean(x$x, 2)
      myv <- mean(x$y, 2)
      sxmat <- cov.freqtab(x$x) * (sum(x$x) - 1)/sum(x$x)
      symat <- cov.freqtab(x$y) * (sum(x$y) - 1)/sum(x$y)
      sx <- sxmat[1]
      sxv <- sxmat[2, 2]
      sy <- symat[1]
      syv <- symat[2, 2]
      xscale <- x$concordance$scale
      if (method == "chained") {
        gd <- dcl(mx, mxv, myv, sqrt(sx), sqrt(sxv),
          sqrt(sy), sqrt(syv), xscale)
      } else if (method %in% c("tucker", "levine")) {
        cxv <- sxmat[1, 2]
        cyv <- symat[1, 2]
        if (method == "tucker") {
          d1x <- 0
          d1xv1 <- 1/sxv
          d1v1 <- -cxv/sxv^2
          d2y <- 0
          d2yv2 <- 1/syv
          d2v2 <- -cyv/syv^2
        } else if (method == "levine" & x$internal) {
          d1x <- 1/cxv
          d1xv1 <- -sx/cxv^2
          d1v1 <- 0
          d2y <- 1/cyv
          d2yv2 <- -sy/cyv^2
          d2v2 <- 0
        } else if (method == "levine" & !x$internal) {
          d1x <- 1/(sxv + cxv)
          d1xv1 <- (sxv - sx)/(sxv + cxv)^2
          d1v1 <- -(sx + cxv)/(sxv + cxv)^2
          d2y <- 1/(syv + cyv)
          d2yv2 <- (syv - sy)/(syv + cyv)^2
          d2v2 <- -(sy + cyv)/(syv + cyv)^2
        }
        gld <- dlin(x, xscale)
        gsd <- ds(mxv, myv, sxv, syv, x$ws,
          x$gamma[1], x$gamma[2],
          d1x, d1xv1, d1v1, d2y, d2yv2, d2v2)
        gd <- gld %*% gsd
      }
      out <- data.frame(n = sqrt(diag(gd %*% omegan %*% t(gd))),
        g = sqrt(diag(gd %*% omegag %*% t(gd))))
    } else out <- NULL
  } else out <- NULL
  return(out)
}

# Single group
sgmse <- function(x, y, xsg) {
  return(NULL)
}

# Mean equivalent groups
egmse <- function(x, y) {
  return(NULL)
}

# Linear equivalent groups
# Braun and Holland (1982) Observed-score test equating: A mathematical
# analysis of some ETS equating procedures
eglse <- function(x, y) {
  nx <- sum(x)
  ny <- sum(y)
  skewterm <- skew.freqtab(x)/nx + skew.freqtab(y)/ny
  kurtterm <- (kurt.freqtab(x) - 1)/(4 * nx) +
    (kurt.freqtab(y) - 1)/(4 * ny)
  xz <- (scales(x) - mean(x))/sd.freqtab(x)
  out <- sqrt(var.freqtab(y) * (1/nx + 1/ny + skewterm * xz +
    kurtterm * xz^2))
  return(out)
}

# Linear single group
sglse <- function(x, y, xsg) {
  return(NULL)
}

# Linear counterbalanced
cblse <- function(x, y) {
  xscale <- scales(x, 1)
  sigmax <- sqrt((var.freqtab(x, 1) + var.freqtab(y, 1))/2)
  sigmay <- (var.freqtab(x, 2) + var.freqtab(y, 2))/2
  mux <- (mean(x, 1) + mean(y, 1))/2
  muy <- (mean(x, 2) + mean(y, 2))/2
  z <- (xscale - mux)/sigmax
  ns <- sum(x) + sum(y)
  rho <- cor.freqtab(x)[2]
  n <- sigmay^2 * (1 - rho) * ((z^2 * (1 + rho) + 2)/ns)
  return(n)
}

# Linear derivatives
dlin <- function(x, xscale) {
  msx <- x$synth$mean[1]
  ssx <- x$synth$sd[1]
  ssy <- x$synth$sd[2]
  return(cbind(-ssy/ssx, 1, -0.5 * ssy * (xscale - msx)/ssx^3,
    0.5 * (xscale - msx)/ssx/ssy))
}

# Synthetic derivatives
# mx, mxv, myv, sx, sxv, sy, syv, xscale
ds <- function(mxv, myv, sxv, syv, wx,
  g1, g2, d1x, d1xv1, d1v1, d2y, d2yv2, d2v2) {
  
  wy <- 1 - wx
  temp1 <- 2 * g1 * (-wy * (sxv - syv) + wx * wy *
      (mxv - myv)^2)
  temp2 <- 2 * g2 * (wx * (sxv - syv) + wx * wy *
      (mxv - myv)^2)
  out <- rbind(c(1, -wy * g1, -wy * (mxv - myv) * d1x,
      -wy * (mxv - myv) * d1xv1, -wy * (mxv - myv) * d1v1,
      0, wy * g1, 0, 0, 0),
    c(0, wx * g2, 0, 0, 0, 1, -wx * g2,
      wx * (mxv - myv) * d2y, wx * (mxv - myv) * d2yv2,
      wx * (mxv - myv) * d2v2),
    c(0, 2 * wx * wy * g1^2 * (mxv - myv),
      1 + temp1 * d1x, temp1 * d1xv1, -wy * g1^2 + temp1 * d1v1,
      0, -2 * wx * wy * g1^2 * (mxv - myv), 0, 0, wy * g1^2),
    c(0, 2 * wx * wy * g2^2 * (mxv - myv),
      0, 0, wx * g2^2, 0, -2 * wx * wy * g2^2 * (mxv - myv),
      1 + temp2 * d2y, temp2 * d2yv2, -wx * g2^2 + temp2 * d2v2))
  
  return(out)
}

# Chained derivatives
dcl <- function(mx, mxv, myv, sx, sxv, sy, syv, xscale) {
  temp <- sxv/sx * (xscale - mx) + mxv - myv
  return(cbind(-sxv * sy/sx/syv, sy/syv,
    -0.5 * sxv * sy/sx^3/syv * (xscale - mx), 0,
    0.5 * sy/sx/syv/sxv * (xscale - mx), 1, -sy/syv,
    0.5/syv/sy * temp, 0, -0.5 * sy/syv^3 * temp))
}

# Duplication matrix
dp <- function(p) {
  out <- matrix(0, p * p, p * (p + 1)/2)
  count <- 0
  for (j in 1:p) {
    for (i in j:p) {
      count <- count + 1
      if (i == j) {
        out[(j - 1) * p + j, count] <- 1
      } else {
        out[(j - 1) * p + i, count] <- 1
        out[(i - 1) * p + j, count] <- 1
      }
    }
  }
  return(out)
}

# Function for getting omegas
omeg <- function(x) {
  xn <- sum(x)
  xm <- margins(x)
  ps <- xm * (xm + 1)/2 
  
  mx <- mean(x, 1:xm)
  smat <- cov.freqtab(x) * (xn - 1)/xn
  xr <- x[x > 0]
  xd <- as.matrix(as.data.frame(x)[x > 0, ])
  
  # omegan
  dup <- dp(xm)
  ismat <- solve(smat)
  w <- .5 * t(dup) %*% (ismat %x% ismat) %*% dup
  omegan <- rbind(cbind(smat, matrix(0, xm, ps)),
    cbind(matrix(0, ps, xm), solve(w)))/xn
  
  # omegag
  s12 <- matrix(0, xm, ps)
  s22 <- matrix(0, ps, ps)
  zi0 <- sweep(xd[, 1:xm], 2, mx)
  for (i in seq_along(xr)) {
    difi <- zi0[i, ] %*% t(zi0[i, ]) - smat
    vdifi <- difi[!upper.tri(difi)]
    s12 <- s12 + (zi0[i, ] %*% t(vdifi)) * xr[i]
    s22 <- s22 + (vdifi %*% t(vdifi)) * xr[i]
  }
  omegag <- rbind(cbind(smat, s12/xn),
    cbind(t(s12)/xn, s22/xn))/xn
  
  return(list(n = omegan, g = omegag))
}

#----------------------------------------------------------------
