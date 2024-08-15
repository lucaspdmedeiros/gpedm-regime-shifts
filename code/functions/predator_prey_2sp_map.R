# Function that iterates a 2-species predator-prey map with noise

# Arguments:
# t_max: number of time steps
# x: vector of initial abundances
# p: vector of parameters 
# control_change: rate of change of control parameter

# Output:
# data frame with abundances over time

predator_prey_2sp_map <- function(t_max, x, p, control_change, control_noise) {
  # initial condition
  x1 <- x[1]
  x2 <- x[2]
  # iterate map
  for (t in 1:(t_max-1)) {
    z <- rnorm(2, p[7], p[8])
    x1[t+1] <- x1[t] * exp((p[1] * (1 - x1[t] / p[2]) - (p[3] * x2[t]) / (1 + p[3] * p[4] * x1[t]))) * exp(z[1])
    x2[t+1] <- x2[t] * exp(((p[5] * p[3] * x1[t]) / (1 + p[3] * p[4] * x1[t]) - p[6])) * exp(z[2])
    p[3] <- p[3] + control_change
  }
  # return results data frame
  results_df <- data.frame(time = 1:t_max, x1 = x1, x2 = x2)
  return(results_df)
}
