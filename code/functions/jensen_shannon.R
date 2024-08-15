# Function that computes the Jensen-Shannon divergence for two samples 
# from continuous probability distributions

# Arguments:
# x: sample from first distribution
# y: sample from second distribution
# start_bin: starting point of first bin
# end_bin: end point of last bin
# bin_number: number of bins to use

# Output:
# jsd: value of the Jensen-Shannon divergence

jensen_shannon <- function(x, y, start_bin, end_bin, bin_number) {
  # sorting x and y
  x_sorted <- sort(x)
  y_sorted <- sort(y)
  # establishing bins
  bin_cutoff <- seq(start_bin, end_bin, length = bin_number)
  bin_lower <- c(bin_cutoff[-length(bin_cutoff)])
  bin_upper <- c(bin_cutoff[-1])
  # compute frequencies of x and y inside each bin
  freq_x <- c()
  freq_y <- c()
  for (i in 1:length(bin_lower)) {
    freq_x[i] <- length(x[x >= bin_lower[i] & x < bin_upper[i]])
    freq_y[i] <- length(y[y >= bin_lower[i] & y < bin_upper[i]])
  }
  freq_x <- freq_x / sum(freq_x)
  freq_y <- freq_y / sum(freq_y)
  # compute frequencies for mixture distribution
  freq_z <- 0.5 * (freq_x + freq_y)
  # compute Kullback–Leibler divergences
  kl_x <- sum(freq_x * log2(freq_x/freq_z), na.rm = TRUE)
  kl_y <- sum(freq_y * log2(freq_y/freq_z), na.rm = TRUE)
  # compute Jensen–Shannon divergence
  jsd <- 0.5 * kl_x + 0.5 * kl_y
  # return value
  return(jsd)
}
