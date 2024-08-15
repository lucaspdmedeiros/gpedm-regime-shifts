# Function that iterates the single-species logistic map with noise

# Arguments:
# t_max: number of time steps
# x: vector of initial abundances
# p: vector of parameters 
# control_change: rate of change of control parameter
# control noise: noise level in control values

# Output:
# data frame with abundances over time

logistic_map <- function(t_max, x, p, control_change, control_noise) {
  # initial condition
  x1 <- x[1]
  # iterate map
  for (t in 1:(t_max-1)) {
    z <- rnorm(1, p[2], p[3])
    x1[t+1] <- p[1] * x1[t] * (1 - x1[t]) * exp(z)
    if (x1[t+1] >= 1) {
      x1[t+1] <- 0.9999
    }
    if (x1[t+1] <= 0) {
      x1[t+1] <- 0.0001
    }
    p[1] <- p[1] + control_change
  }
  # return results data frame
  results_df <- data.frame(time = 1:t_max, x1 = x1)
  return(results_df)
}
