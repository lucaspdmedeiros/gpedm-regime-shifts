# Function that iterates the single-species harvesting map with noise

# Arguments:
# t_max: number of time steps
# x: vector of initial abundances
# p: vector of parameters 
# control_change: rate of change of control parameter

# Output:
# data frame with abundances over time

harvesting_1sp_map <- function(t_max, x, p, control_change, control_noise) {
  # initial condition
  x1 <- x[1]
  # iterate map
  for (t in 1:(t_max-1)) {
    z <- rnorm(1, p[5], p[6])
    x1[t+1] <- x1[t] * exp(p[1] * (1 - x1[t] / p[2]) - (p[3] * x1[t]) / (p[4]^2 + x1[t]^2)) * exp(z)
    p[3] <- p[3] + control_change
  }
  # return results data frame
  results_df <- data.frame(time = 1:t_max, x1 = x1)
  return(results_df)
}
