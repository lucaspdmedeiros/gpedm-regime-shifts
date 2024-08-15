# Function that iterates a 3-species competition map with noise

# Arguments:
# t_max: number of time steps
# x: vector of initial abundances
# p: vector of parameters 
# control_change: rate of change of control parameter

# Output:
# data frame with abundances over time

competition_3sp_map <- function(t_max, x, p, control_change, control_noise) {
  # initial condition
  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  # iterate map
  for (t in 1:(t_max-1)) {
    z <- rnorm(3, p[8], p[9])
    x1[t+1] <- x1[t] * exp(p[1] * (1 - x1[t] - p[2] * x2[t] - p[3] * x3[t])) * exp(z[1])
    x2[t+1] <- x2[t] * exp(p[1] * (1 - p[4] * x1[t] - x2[t] - p[5] * x3[t])) * exp(z[2])
    x3[t+1] <- x3[t] * exp(p[1] * (1 - p[6] * x1[t] - p[7] * x2[t] - x3[t])) * exp(z[3])
    p[1] <- p[1] + control_change
  }
  # return results data frame
  results_df <- data.frame(time = 1:t_max, x1 = x1, x2 = x2, x3 = x3)
  return(results_df)
}
