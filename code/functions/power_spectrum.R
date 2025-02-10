# Computes the power spectrum from a time series
# (code adapted from Rogers et al (2023) Ecol Lett)

# Arguments:
# ts: time series to use (time in 1st column, state variable in 2nd column)
# lam: regularization coefficient
# scale: whether to scale to zero mean and unit standard deviation
# trim: whether to trim NAs at the time series edges
# plot: whether to plot results

# Output:
# results: list containing (1) data frame with all frequencies, periods, and power; and
# (2) vector with frequency and period of maximum power

power_spectrum <- function(ts, lam, scale, trim, plot) {
  # select species time series
  xts <- ts[ , 2]
  # standardize data
  if (sd(xts) > 0 & scale == TRUE) {
    xtss <- (xts - mean(xts, na.rm = TRUE)) / sd(xts, na.rm = TRUE)
  } else {
    xtss <- xts - mean(xts, na.rm = TRUE)
  }
  # get spectrum
  if (trim) {  # ts length, exclude NAs on ends
    tsl <- length(which(!is.na(xtss))[1]:which(!is.na(xtss))[length(which(!is.na(xtss)))])
  } else { # do not trim NAs on ends
    tsl <- length(xtss)
  }
  tsle <- floor(tsl/2) * 2 #ts length rounded down to even number (if odd)
  fi <- (1:(tsle/2)) / tsle # frequencies (cycles per timestep)
  per <- 1 / fi # periods (timesteps per cycle)
  wi <- 2 * pi * fi # frequencies in radians per timestep
  times <- 1:length(xtss)
  cosbf <- cos(outer(times, wi))
  sinbf <- sin(outer(times,wi))
  allbf <- cbind(cosbf,sinbf) # all basis functions
  y <- xtss[complete.cases(xtss)] # remove missing timepoints
  X <- allbf[complete.cases(xtss), ] # remove missing timepoints
  coefsdir <- solve(t(X) %*% X + lam*diag(ncol(X))) %*% t(X) %*% y
  lmr <- sqrt(coefsdir[1:(length(coefsdir)/2), ]^2 + coefsdir[(length(coefsdir)/2+1):length(coefsdir), ]^2)
  results_df <- data.frame(frequency = fi, period = per, power = lmr)
  freq_max <- results_df$frequency[which.max(results_df$power)] # frequency of maximum power
  period_max <- results_df$period[which.max(results_df$power)] # period of maximum power
  results <- list(results_df, c(freq_max, period_max))
  # plot power as a function of period
  if (plot) {
    plot(results_df$frequency, results_df$power, type = "l")
  }
  return(results)
}
