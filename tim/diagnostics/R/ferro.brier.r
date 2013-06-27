# Point and Interval Estimates for the Brier Score
#
# Description:
#
#      Returns point estimates, standard errors and bootstrapped
#      pointwise and simultaneous confidence intervals for both
#      conditional and unconditional Brier scores of ensemble-based
#      probabilistic forecasts of binary events.
#
# Usage:
#
#      brier.score(x, y = NULL, th, size, above = TRUE, con = TRUE,
#                  unc = TRUE, cover = NULL, boot = 1000, simul = FALSE)
#
# Arguments:
#
#      x: Matrix with observations in the first row and members of
#         corresponding ensemble forecasts in subsequent rows.
#
#      y: Optional matrix with second set of observations and forecasts
#         in the same form as `x'.
#
#     th: Vector of thresholds defining different events. The same
#         thresholds are used for the observations and forecasts.
#
#   size: Vector of ensemble sizes for which to estimate the expected
#         Brier scores.
#
#  above: Logical: if `TRUE' (the default), events are defined by
#         threshold exceedance; if `FALSE', events are defined by
#         threshold deficit.
#
#    con: Logical: if `TRUE' (the default), computes standard errors for
#         the true conditional Brier scores.
#
#    unc: Logical: if `TRUE' (the default), computes standard errors for
#         the true unconditional Brier scores.
#
#  cover: Coverage probability for confidence intervals; if `NULL', no
#         confidence intervals are computed.
#
#   boot: Number of bootstrap samples from which to compute confidence
#         intervals.
#
#  simul: Logical: if `FALSE' (the default), only pointwise confidence
#         intervals are computed; if `TRUE', simultaneous confidence
#         intervals for all thresholds are also computed.
#
# Details:
#
#      See Ferro (2007).
#
# Output:
#
#      A list with three components (`x', `y' and `difference'
#      corresponding to the scores for the x- and y-forecasts and the
#      differences between the scores for the x- and y-forecasts) is
#      returned invisibly. Each component is itself a list containing
#      the following elements:
#
#      members: Vector of ensemble sizes `size'.
#
#   thresholds: Vector of thresholds `th'.
#
#       scores: Matrix of estimates of the Brier score for each
#                 threshold (rows) and ensemble size (columns).
#
#  reliability: Matrix of the reliability components for the estimates
#               in `scores'.
#
#   resolution: Matrix of the resolution components for the estimates
#               in `scores'.
#
#    certainty: Matrix of the uncertainty components for the estimates
#               in `scores'.
#
#   adjustment: Matrix of the adjustment components for the estimates
#               in `scores'.
#
#          con: List containing the following components:
#
#               std.err: Matrix of conditional standard errors for the
#                        estimates in `scores'.
#
#                 point: List containing lower and upper limits of
#                        conditional pointwise confidence intervals for
#                        the estimates in `scores'.
#
#                 simul: List containing lower and upper limits of
#                        conditional simultaneous confidence intervals
#                        for the estimates in `scores', plus estimates of
#                        the global and local coverages of the intervals.
#
#          unc: List with the same form as `con' but for unconditional
#               standard errors and intervals.
#
# Author:
#
#      Chris Ferro <c.a.t.ferro@reading.ac.uk> 30 January 2007
#
# References:
#
#      Ferro CAT (2007) Comparing probabilistic forecasting systems with
#      the Brier score. Weather and Forecasting. In press. Available at
#
#      http://www.met.rdg.ac.uk/~sws02caf/Publications/brier.pdf
#
#      $Id$
#
# Examples:
#
#      obs <- rnorm(30)
#      mod <- matrix(NA, 5, 30)
#      for(i in 1:5) mod[i, ] <- rnorm(30, 0.5 * obs, sqrt(1 - 0.5^2))
#      data <- rbind(obs, mod)
#      out <- brier.score(data,, c(0, 2), c(5, Inf))

brier.score <- function(x, y = NULL, th, size, above = TRUE, con = TRUE, unc = TRUE, cover = NULL, boot = 1000, simul = FALSE) {
  # check arguments
  if(!is.null(cover) && !con && !unc)
    stop("One of `con' and `unc' must equal `TRUE'")
  # negate data if below threshold
  if(!above) {
    x <- -x
    if(!is.null(y)) y <- -y
    th <- -th
  }
  # data
  xobs <- x[1, ]
  xens <- x[-1,, drop = FALSE]
  
  if(!is.null(y)) {
    yobs <- y[1, ]
    yens <- y[-1,, drop = FALSE]
  }
  if(is.null(cover))
    alpha <- NULL
  else
    alpha <- (1 - cover) / 2
  # dimensions
#  m <- nrow(xens) # not needed? need mx and my?
  n <- length(xobs)
  nt <- length(th)
  ns <- length(size)
  # point estimates and standard errors
  temp <- brier.point(xobs, xens, th, size, con, unc)
  xb <- temp$b
  xb0 <- temp$b0
  xs.con <- temp$s.con
  xs.unc <- temp$s.unc
  xout <- list(members = size, thresholds = th, scores = xb, reliability = temp$relia, resolution = temp$resol, certainty = temp$uncer, adjustment = temp$adjus)
  if(con) xout$con$std.err <- xs.con
  if(unc) xout$unc$std.err <- xs.unc
  if(!is.null(y)) {
    temp <- brier.point(yobs, yens, th, size, con, unc)
    yb <- temp$b
    yb0 <- temp$b0
    ys.con <- temp$s.con
    ys.unc <- temp$s.unc
    yout <- list(members = size, thresholds = th, scores = yb, reliability = temp$relia, resolution = temp$resol, certainty = temp$uncer, adjustment = temp$adjus)
    if(con) yout$con$std.err <- ys.con
    if(unc) yout$unc$std.err <- ys.unc
    d <- xb - yb
    covar <- brier.cov(xobs, xens, yobs, yens, th, size)
    ds.con <- sqrt(xs.con^2 + ys.con^2)
    ds.unc <- sqrt(xs.unc^2 + ys.unc^2 - 2 * covar / n)
    dout <- list(members = size, thresholds = th, scores = d)
    if(con) dout$con$std.err <- ds.con
    if(unc) dout$unc$std.err <- ds.unc
  }
  # return point estimates
  if(is.null(alpha)) {
    if(!above) th <- -th
    xout$thresholds <- th
    if(is.null(y)) {
      out <- xout
    } else {
      yout$thresholds <- th
      dout$thresholds <- th
      out <- list(x = xout, y = yout, difference = dout)
    }
    return(invisible(out))
  }
  # pointwise confidence intervals
  if(con) {
    # conditional
    print("Conditional bootstrap (x)...")
    temp <- brier.boot(xobs, xens,,, th, size, boot, con = TRUE)
    xbb.con <- temp$x$score
    xsb.con <- temp$x$std.err
    temp <- brier.int(xbb.con, xsb.con, xb0, xs.con, xb, alpha, con = TRUE)
    xrb.con <- temp$r
    xtb.con <- temp$t
    xlo.con <- temp$lo
    xup.con <- temp$up
    xout$con$point$lower <- xlo.con
    xout$con$point$upper <- xup.con
    if(!is.null(y)) {
      print("Conditional bootstrap (y)...")
      temp <- brier.boot(yobs, yens,,, th, size, boot, con = TRUE)
      ybb.con <- temp$x$score
      ysb.con <- temp$x$std.err
      temp <- brier.int(ybb.con, ysb.con, yb0, ys.con, yb, alpha, con = TRUE)
      yrb.con <- temp$r
      ytb.con <- temp$t
      ylo.con <- temp$lo
      yup.con <- temp$up
      yout$con$point$lower <- ylo.con
      yout$con$point$upper <- yup.con
      dsb.con <- sqrt(xsb.con^2 + ysb.con^2)
      dtb.con <- (xrb.con - yrb.con) / dsb.con
      zb <- aperm(apply(dtb.con, c(1, 2), quantile, probs = c(1 - alpha, alpha), na.rm = TRUE), c(2, 3, 1))
      dlo.con <- d - ds.con * zb[, , 1]
      dup.con <- d - ds.con * zb[, , 2]
      dout$con$point$lower <- dlo.con
      dout$con$point$upper <- dup.con
    }
  }
  if(unc) {
    # unconditional
    if(is.null(y)) {
      print("Unconditional bootstrap (x)...")
      temp <- brier.boot(xobs, xens,,, th, size, boot, con = FALSE)
      xbb.unc <- temp$x$score
      xsb.unc <- temp$x$std.err
      temp <- brier.int(xbb.unc, xsb.unc, xb0, xs.unc, xb, alpha, con = FALSE)
      xrb.unc <- temp$r
      xtb.unc <- temp$t
      xlo.unc <- temp$lo
      xup.unc <- temp$up
      xout$unc$point$lower <- xlo.unc
      xout$unc$point$upper <- xup.unc
    } else {
      print("Unconditional bootstrap (x and y)...")
      temp <- brier.boot(xobs, xens, yobs, yens, th, size, boot, con = FALSE)
      xbb.unc <- temp$x$score
      xsb.unc <- temp$x$std.err
      ybb.unc <- temp$y$score
      ysb.unc <- temp$y$std.err
      covarb <- temp$cov
      temp <- brier.int(xbb.unc, xsb.unc, xb0, xs.unc, xb, alpha, con = FALSE)
      xrb.unc <- temp$r
      xtb.unc <- temp$t
      xlo.unc <- temp$lo
      xup.unc <- temp$up
      xout$unc$point$lower <- xlo.unc
      xout$unc$point$upper <- xup.unc
      temp <- brier.int(ybb.unc, ysb.unc, yb0, ys.unc, yb, alpha, con = FALSE)
      yrb.unc <- temp$r
      ytb.unc <- temp$t
      ylo.unc <- temp$lo
      yup.unc <- temp$up
      yout$unc$point$lower <- ylo.unc
      yout$unc$point$upper <- yup.unc
      dsb.unc <- sqrt(xsb.unc^2 + ysb.unc^2 - 2 * covarb / n)
      dtb.unc <- (xrb.unc - yrb.unc) / dsb.unc
      zb <- aperm(apply(dtb.unc, c(1, 2), quantile, probs = c(1 - alpha, alpha), na.rm = TRUE), c(2, 3, 1))
      dlo.unc <- d - ds.unc * zb[, , 1]
      dup.unc <- d - ds.unc * zb[, , 2]
      dout$unc$point$lower <- dlo.unc
      dout$unc$point$upper <- dup.unc
    }
  }
  # return pointwise intervals
  if(!simul) {
    if(!above) th <- -th
    xout$thresholds <- th
    if(is.null(y)) {
      out <- xout
    } else {
      yout$thresholds <- th
      dout$thresholds <- th
      out <- list(x = xout, y = yout, difference = dout)
    }
    return(invisible(out))
  }
  # simultaneous confidence intervals
  if(con) {
    # conditional
    temp <- brier.sim(xtb.con, xb, xs.con, alpha)
    xlo.con.sim <- pmax(temp$lo, 0)
    xup.con.sim <- pmin(temp$up, 1)
    xcov.con.sim <- temp$sim
    xcov.con.loc <- temp$loc
    xout$con$simul$lower <- xlo.con.sim
    xout$con$simul$upper <- xup.con.sim
    xout$con$simul$cover <- rbind(global = xcov.con.sim, local = xcov.con.loc)
    if(!is.null(y)) {
      temp <- brier.sim(ytb.con, yb, ys.con, alpha)
      ylo.con.sim <- pmax(temp$lo, 0)
      yup.con.sim <- pmin(temp$up, 1)
      ycov.con.sim <- temp$sim
      ycov.con.loc <- temp$loc
      yout$con$simul$lower <- ylo.con.sim
      yout$con$simul$upper <- yup.con.sim
      yout$con$simul$cover <- rbind(global = ycov.con.sim, local = ycov.con.loc)
      temp <- brier.sim(dtb.con, d, ds.con, alpha)
      dlo.con.sim <- temp$lo
      dup.con.sim <- temp$up
      dcov.con.sim <- temp$sim
      dcov.con.loc <- temp$loc
      dout$con$simul$lower <- dlo.con.sim
      dout$con$simul$upper <- dup.con.sim
      dout$con$simul$cover <- rbind(global = dcov.con.sim, local = dcov.con.loc)
    }
  }
  if(unc) {
    # unconditional
    temp <- brier.sim(xtb.unc, xb, xs.unc, alpha)
    xlo.unc.sim <- pmax(temp$lo, 0)
    xup.unc.sim <- pmin(temp$up, 1)
    xcov.unc.sim <- temp$sim
    xcov.unc.loc <- temp$loc
    xout$unc$simul$lower <- xlo.unc.sim
    xout$unc$simul$upper <- xup.unc.sim
    xout$unc$simul$cover <- rbind(global = xcov.unc.sim, local = xcov.unc.loc)
    if(!is.null(y)) {
      temp <- brier.sim(ytb.unc, yb, ys.unc, alpha)
      ylo.unc.sim <- pmax(temp$lo, 0)
      yup.unc.sim <- pmin(temp$up, 1)
      ycov.unc.sim <- temp$sim
      ycov.unc.loc <- temp$loc
      yout$unc$simul$lower <- ylo.unc.sim
      yout$unc$simul$upper <- yup.unc.sim
      yout$unc$simul$cover <- rbind(global = ycov.unc.sim, local = ycov.unc.loc)
      temp <- brier.sim(dtb.unc, d, ds.unc, alpha)
      dlo.unc.sim <- temp$lo
      dup.unc.sim <- temp$up
      dcov.unc.sim <- temp$sim
      dcov.unc.loc <- temp$loc
      dout$unc$simul$lower <- dlo.unc.sim
      dout$unc$simul$upper <- dup.unc.sim
      dout$unc$simul$cover <- rbind(global = dcov.unc.sim, local = dcov.unc.loc)
    }
  }
  # output
  if(!above) th <- -th
  xout$thresholds <- th
  if(is.null(y)) {
    out <- xout
  } else {
    yout$thresholds <- th
    dout$thresholds <- th
    out <- list(x = xout, y = yout, difference = dout)
  }
  invisible(out)
}


brier.point <- function(obs, sim, th, size, con, unc) {
  m <- nrow(sim)
  n <- length(obs)
  ns <- length(size)
  nt <- length(th)
  b <- b0 <- matrix(NA, nt, ns)
  s.con <- s.unc <- matrix(NA, nt, ns)
  relia <- resol <- adjus <- uncer <- matrix(NA, nt, ns)
  for(i in 1:nt) {
    I <- (obs > th[i])
    Q <- apply(sim > th[i], 2, mean)
    Ibar <- mean(I)
    Qbar <- table(I, Q)[, as.character(Q)]
    Qbar <- Qbar["TRUE", ] / apply(Qbar, 2, sum)
    relia[i, ] <- mean((Q - Qbar)^2)
    resol[i, ] <- mean((Qbar - Ibar)^2)
    uncer[i, ] <- Ibar * (1 - Ibar)
    sharp <- mean((Q - 0.5)^2)
    adjus[i, ] <- (1 - m / size) * (0.25 - sharp) / (m - 1)
    for(j in 1:ns) {
      b[i, j] <- brier1(I, Q, m, size[j])
      b0[i, j] <- brier.true(I, Q, size[j])
      if(con) s.con[i, j] <- brier.sd(I, Q, m, n, size[j])
      if(unc) s.unc[i, j] <- brier.sd(I, Q, m, n, size[j], TRUE)
    }
  }
  list(b = b, b0 = b0, s.con = s.con, s.unc = s.unc, relia = relia, resol = resol, uncer = uncer, adjus = adjus)
}


brier.cov <- function(xobs, xsim, yobs, ysim, th, size) {
  mx <- nrow(xsim)
  my <- nrow(ysim)
  nt <- length(th)
  ns <- length(size)
  covar <- matrix(NA, nt, ns)
  for(i in 1:nt) {
    Ix <- (xobs > th[i])
    Qx <- apply(xsim > th[i], 2, mean)
    Iy <- (yobs > th[i])
    Qy <- apply(ysim > th[i], 2, mean)
    for(j in 1:ns) {
      M <- size[j]
      Wx <- (Qx - Ix)^2 - (1 - mx / M) * Qx * (1 - Qx) / (mx - 1)
      Wy <- (Qy - Iy)^2 - (1 - my / M) * Qy * (1 - Qy) / (my - 1)
      covar[i, j] <- mean(Wx * Wy) - mean(Wx) * mean(Wy)
    }
  }
  covar
}


brier1 <- function(I, Q, m, M) {
  if(m == 1 && M == 1) return(mean((Q - I)^2))
  mean((Q - I)^2 - (1 - m / M) * Q * (1 - Q) / (m - 1))
}


brier.sd <- function(I, Q, m, n, M, sample = FALSE) {
  if(sample) {
    return(sd((Q - I)^2 - (1 - m / M) * Q * (1 - Q) / (m - 1)) / sqrt(n))
  }
  if(m < 4) {
    Q2 <- Q * Q
    Q3 <- Q * Q2
    Q4 <- Q * Q3
  } else {
    Q2 <- Q * (m * Q - 1) / (m - 1)
    Q3 <- Q2 * (m * Q - 2) / (m - 2)
    Q4 <- Q3 * (m * Q - 3) / (m - 3)
  }
  sqrt((1 - 1 / M) * max(0, 2 * (3 - 2 * m) * (1 - 1 / M) * mean(Q4) / (m * (m - 1)) + 4 * (m - 2 + (3 - 2 * m) / M) * mean(Q3) / (m * (m - 1)) + 2 * (1 + (2 * m - 3) / (M - 1) + (7 - 5 * m) / (2 * M * (M - 1))) * mean(Q2) / (m * (m - 1)) + mean(Q) / (m * (m - 1) * M) + 8 * mean(I * Q3) / m - 12 * mean(I * Q2) / m + 4 * mean(I * Q) / m)) / sqrt(n)
}


brier.true <- function(I, Q, M) {
  (1 - 1 / M) * mean(Q^2) + mean(Q) / M - 2 * mean(I * Q) + mean(I)
}


brier.boot <- function(xobs, xsim, yobs = NULL, ysim = NULL, th, size, boot, con) {
  mx <- nrow(xsim)
  my <- nrow(ysim)
  n <- length(xobs)
  ns <- length(size)
  nt <- length(th)
  covb <- array(NA, c(nt, ns, boot))
  xbb <- xsb <- array(NA, c(nt, ns, boot))
  ybb <- ysb <- array(NA, c(nt, ns, boot))
  if(con) {
    # conditional x
    for(k in 1:boot) {
      xsimb <- apply(xsim, 2, sample, replace = TRUE)
      for(i in 1:nt) {
        I <- (xobs > th[i])
        Q <- apply(xsimb > th[i], 2, mean)
        for(j in 1:ns) {
          xbb[i, j, k] <- brier1(I, Q, mx, size[j])
          xsb[i, j, k] <- brier.sd(I, Q, mx, n, size[j])
        }
      }
    }
    out <- list(x = list(score = xbb, std.err = xsb))
    # conditional y
    if(!is.null(yobs)) {
      for(k in 1:boot) {
        ysimb <- apply(ysim, 2, sample, replace = TRUE)
        for(i in 1:nt) {
          I <- (yobs > th[i])
          Q <- apply(ysimb > th[i], 2, mean)
          for(j in 1:ns) {
            ybb[i, j, k] <- brier1(I, Q, my, size[j])
            ysb[i, j, k] <- brier.sd(I, Q, my, n, size[j])
          }
        }
      }
      out$y <- list(score = ybb, std.err = ysb)
    }
  } else {
    # unconditional x
    if(is.null(yobs)) {
      for(k in 1:boot) {
        index <- sample(n, replace = TRUE)
        obsb <- xobs[index]
        simb <- xsim[, index]
        for(i in 1:nt) {
          I <- (obsb > th[i])
          Q <- apply(simb > th[i], 2, mean)
          for(j in 1:ns) {
            xbb[i, j, k] <- brier1(I, Q, mx, size[j])
            xsb[i, j, k] <- brier.sd(I, Q, mx, n, size[j], TRUE)
          }
        }
      }
      out <- list(x = list(score = xbb, std.err = xsb))
    } else {
      # unconditional x and y
      for(k in 1:boot) {
        index <- sample(n, replace = TRUE)
        xobsb <- xobs[index]
        yobsb <- yobs[index]
        xsimb <- xsim[, index]
        ysimb <- ysim[, index]
        for(i in 1:nt) {
          Ix <- (xobsb > th[i])
          Qx <- apply(xsimb > th[i], 2, mean)
          Iy <- (yobsb > th[i])
          Qy <- apply(ysimb > th[i], 2, mean)
          for(j in 1:ns) {
            xbb[i, j, k] <- brier1(Ix, Qx, mx, size[j])
            xsb[i, j, k] <- brier.sd(Ix, Qx, mx, n, size[j], TRUE)
            ybb[i, j, k] <- brier1(Iy, Qy, my, size[j])
            ysb[i, j, k] <- brier.sd(Iy, Qy, my, n, size[j], TRUE)
          }
        }
        covb[,, k] <- brier.cov(xobsb, xsimb, yobsb, ysimb, th, size)
      }
      out <- list(x = list(score = xbb, std.err = xsb))
      out$y <- list(score = ybb, std.err = ysb)
      out$cov <- covb
    }
  }
  out
}


brier.int <- function(bb, sb, b0, s, b, alpha, con) {
  nt <- dim(bb)[1]
  ns <- dim(bb)[2]
  boot <- dim(bb)[3]
  lo <- up <- matrix(NA, nt, ns)
  rb <- tb <- array(NA, c(nt, ns, boot))
  if(con) {
    for(i in 1:nt) {
      for(j in 1:ns) {
        rb[i, j, ] <- bb[i, j, ] - b0[i, j]
        tb[i, j, ] <- rb[i, j, ] / sb[i, j, ]
        zb <- quantile(tb[i, j, ], c(1 - alpha, alpha), na.rm = TRUE)
        if(s[i, j] == 0) {
          lo[i, j] <- b[i, j]
          up[i, j] <- b[i, j]
        } else {
          lo[i, j] <- b[i, j] - s[i, j] * zb[1]
          up[i, j] <- b[i, j] - s[i, j] * zb[2]
        }
      }
    }
  } else {
    for(i in 1:nt) {
      for(j in 1:ns) {
        rb[i, j, ] <- bb[i, j, ] - b[i, j]
        tb[i, j, ] <- rb[i, j, ] / sb[i, j, ]
        zb <- quantile(tb[i, j, ], c(1 - alpha, alpha), na.rm = TRUE)
        if(s[i, j] == 0) {
          lo[i, j] <- b[i, j]
          up[i, j] <- b[i, j]
        } else {
          lo[i, j] <- b[i, j] - s[i, j] * zb[1]
          up[i, j] <- b[i, j] - s[i, j] * zb[2]
        }
      }
    }
  }
  lo[lo < 0] <- 0
  up[up > 1] <- 1
  list(r = rb, t = tb, lo = lo, up = up)
}


brier.sim <- function(tb, b, s, alpha) {
  nt <- dim(tb)[1]
  ns <- dim(tb)[2]
  boot <- dim(tb)[3]
  sim <- loc <- numeric(ns)
  lo <- up <- matrix(NA, nt, ns)
  for(j in 1:ns) {
    k <- 1
    tr <- t(apply(tb[, j, ], 1, rank))
    tbs <- t(apply(tb[, j, ], 1, sort, na.last = TRUE))
    out <- apply((tr <= k) | (tr >= boot + 1 - k), 2, any)
    p.global <- sum(out) / boot
    while((p.global < 2 * alpha) && (k < boot)) {
      k <- k + 1
      out <- apply((tr <= k) | (tr >= boot + 1 - k), 2, any)
      p.global <- sum(out) / boot
    }
    ok <- apply(!is.na(tbs), 1, sum)
    for(i in 1:nt) {
      if(s[i, j] == 0) {
        lo[i, j] <- b[i, j]
        up[i, j] <- b[i, j]
      } else {
        lo[i, j] <- b[i, j] - s[i, j] * tbs[i, ok[i] + 1 - k]
        up[i, j] <- b[i, j] - s[i, j] * tbs[i, k]
      }
    }
    sim[j] <- 1 - p.global
    loc[j] <- 1 - 2 * k / boot
  }
  list(sim = sim, loc = loc, lo = lo, up = up)
}
