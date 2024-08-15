# Function that iterates a 2-species competition harvesting map with noise

# Arguments:
# t_max: number of time steps
# x: vector of initial abundances
# p: vector of parameters 
# control_change: rate of change of control parameter
# control noise: noise level in control values

# Output:
# data frame with abundances over time

harvesting_2sp_map <- function(t_max, x, p, control_change, control_noise) {
  # initial condition
  x1 <- x[1]
  x2 <- x[2]
  # control level
  control <- p[4] + rnorm(1, 0, control_noise)
  # iterate map
  for (t in 1:(t_max-1)) {
    z <- rnorm(2, p[7], p[8])
    x1[t+1] <- x1[t] * exp(p[1] * (1 - x1[t] / p[2] - p[3] * x2[t]) - (control[t] * x1[t]) / (p[5]^2 + x1[t]^2)) * exp(z[1])
    x2[t+1] <- x2[t] * exp(p[1] * (1 - x2[t] / p[2] - p[6] * x1[t])) * exp(z[2])
    if (x1[t+1] <= 0) {
      x1[t+1] <- 0.0001
    }
    if (x2[t+1] <= 0) {
      x2[t+1] <- 0.0001
    }
    p[4] <- p[4] + control_change
    control <- c(control, p[4] + rnorm(1, 0, control_noise))
  }
  # return results data frame
  results_df <- data.frame(time = 1:t_max, x1 = x1, x2 = x2, control = control)
  return(results_df)
}
