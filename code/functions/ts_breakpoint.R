# Function that performs the change-point analysis by finding the breakpoint(s) 
# of a time series that minimize the sum of squares computed for each of
# two (or three) time series segments.

# Arguments:
# ts: data frame with time in first column and variable in second column
# n_breaks: number of breakpoints to try
# n_segments: number of time series segments to try (1 or 2)

# Outputs, list with elements:
# t_break: time point that best breaks time series into two different segments
# ss: sum of squares of best outcome

ts_breakpoint <- function(ts, n_breaks, n_segments) {
  if (n_segments == 1) {
    # determine breakpoints to test
    middlepoint <- nrow(ts) / 2
    time_points <- seq(floor(middlepoint - n_breaks/2), floor(middlepoint + n_breaks/2), by = 1)
    # to store results
    rel_ss <- c()
    for (i in 1:length(time_points)) {
      # sum of squares for first segment
      ts1 <- ts[1:time_points[i], ]
      x1 <- as.numeric(unlist(ts1[ , 2]))
      ss1 <- sum((x1 - mean(x1, na.rm = TRUE))^2, na.rm = TRUE)
      # sum of squares for second segment
      ts2 <- ts[(time_points[i]+1):nrow(ts), ]
      x2 <- as.numeric(unlist(ts2[ , 2]))
      ss2 <- sum((x2 - mean(x2, na.rm = TRUE))^2, na.rm = TRUE)
      # total sum of squares
      x <- as.numeric(unlist(ts[ , 2]))
      ss_total <- sum((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE)
      # compute relative sum of squares
      rel_ss[i] <- (ss1 + ss2) / ss_total
    }
    # most likely breakpoint
    t_break <- time_points[which.min(rel_ss)]
    # sum of squares using that breakpoint
    ss_breakpoint <- min(rel_ss)
    # results list
    results_list <- list(t_break, ss_breakpoint)
  }
  if (n_segments == 2) {
    # determine breakpoints to test
    middlepoint <- nrow(ts) / 2
    time_points <- seq(floor(middlepoint - n_breaks/2), floor(middlepoint + n_breaks/2), by = 1)
    time_points_df <- expand.grid(time_points, time_points)
    time_points_df <- time_points_df[time_points_df$Var1 < time_points_df$Var2, ]
    time_points_df <- time_points_df[abs(time_points_df$Var1 - time_points_df$Var2) > floor(nrow(ts) / 10), ]
    # to store results
    rel_ss <- c()
    for (i in 1:nrow(time_points_df)) {
      # sum of squares for first segment
      ts1 <- ts[1:time_points_df$Var1[i], ]
      x1 <- as.numeric(unlist(ts1[ , 2]))
      ss1 <- sum((x1 - mean(x1, na.rm = TRUE))^2, na.rm = TRUE)
      # sum of squares for second segment
      ts2 <- ts[(time_points_df$Var1[i]+1):(time_points_df$Var2[i]), ]
      x2 <- as.numeric(unlist(ts2[ , 2]))
      ss2 <- sum((x2 - mean(x2, na.rm = TRUE))^2, na.rm = TRUE)
      # sum of squares for third segment
      ts3 <- ts[(time_points_df$Var2[i]+1):nrow(ts), ]
      x3 <- as.numeric(unlist(ts3[ , 2]))
      ss3 <- sum((x3 - mean(x3, na.rm = TRUE))^2, na.rm = TRUE)
      # total sum of squares
      x <- as.numeric(unlist(ts[ , 2]))
      ss_total <- sum((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE)
      # compute relative sum of squares
      rel_ss[i] <- (ss1 + ss2 + ss3) / ss_total
    }
    # most likely breakpoints
    t_break <- as.numeric(time_points_df[which.min(rel_ss), ])
    # sum of squares using those breakpoints
    ss_breakpoint <- min(rel_ss)
    # results list
    results_list <- list(t_break, ss_breakpoint)
  }
  return(list(t_break, ss_breakpoint))
}
