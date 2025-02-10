# Code for Fig 3. For the 2-species competition model with harvesting, generate 
# time series while gradually increasing harvest rate until population collapse. 
# Then, train a GP-EDM for successively larger training data sets and extrapolate 
# to unseen values of the control parameter to reconstruct the bifurcation diagram.

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/logistic_1sp_map.R")
source("code/functions/competition_3sp_map.R")
source("code/functions/predator_prey_2sp_map.R")
source("code/functions/harvesting_2sp_map.R")
source("code/functions/ts_breakpoint.R")
if (!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if (!require(gridExtra)) {install.packages("gridExtra"); library(gridExtra)}
if (!require(viridis)) {install.packages("viridis"); library(viridis)}
if (!require(GPEDM)) {install.packages("GPEDM"); library(GPEDM)}
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}

# set up for EWS analysis ------------------------------ 
# set seed to make results reproducible
set.seed(6)
# levels of control parameter to use as time series endpoints
control_gp <- c(0.456, 0.4624, 0.464, 0.472)
# alternative scenario 1
#set.seed(2)
#control_gp <- c(0.4608, 0.4672, 0.4688, 0.4768)
# alternative scenario 2
#set.seed(4)
#control_gp <- c(0.4544, 0.4608, 0.4624, 0.4704)
# whether to save plot
save_plot <- TRUE
# model to use
func_name <- "harvesting_2sp_map"
# whether to use native or delay coordinates
coordinates <- "delay"
# whether to use lags for control
control_lags <- FALSE
# rate of change of control parameter
control_change <- 0
# noise level in control values
control_noise <- 0
# to save results
results_no_noise <- data.frame()
# model settings to generate time series without control parameter change
n_sp <- 2
func <- harvesting_2sp_map
x_init <- c(x1 = 0.5, x2 = 0.5)
control_true <- seq(0.42, 0.48, by = 0.0002)
parms <- list(r = 2.7, k = 1, c = 0.2, p = NA, a = 0.1, b = 0.1, 
              s_mean = NA, s_sd = NA)
log_abund <- FALSE
# full time series length
t_max_full <- 300
# final time series length
t_max <- 100
# generate time series for different control values
for (i in 1:length(control_true)) {
  # parameters without noise
  parms$p <- control_true[i]
  parms$s_mean <- 0
  parms$s_sd <- 0
  # simulate model without noise
  curr_results <- func(t_max = t_max_full, x = x_init, p = unlist(parms), 
                       control_change = 0, control_noise = 0)
  # save results without noise
  curr_results$control <- control_true[i]
  results_no_noise <- rbind(results_no_noise, tail(curr_results, t_max))
}
# change abundances to log scale
if (log_abund) {
  results_no_noise[ , -c(1, ncol(results_no_noise))] <- 
    log(results_no_noise[ , -c(1, ncol(results_no_noise))])
}
# model settings to generate time series with gradual control parameter change
if (func_name == "harvesting_2sp_map") {
  n_sp <- 2
  func <- harvesting_2sp_map
  x_init <- c(x1 = 0.6722805, x2 = 1.8093019)
  control_true <- seq(0.42, 0.48, by = 0.0002)
  s <- 0.02
  parms <- list(r = 2.7, k = 1, c = 0.2, p = 0.4, a = 0.1, b = 0.1, s_mean = -(s^2)/2, s_sd = s)
  control_change <- 0.0004
  control_noise <- 0
  control_name <- "Harvest rate"
  log_abund <- FALSE
}
# grid of embedding dimension values and time lag value
E_values <- 1:9
tau <- 1
# defining coordinates to use
if (coordinates == "native") {
  E <- 1
  sp_plot <- 1
  if (n_sp == 1) {
    sp_input <- 1
  }
  if (n_sp == 2) {
    sp_input <- c(1, 2)
  }
  if (n_sp == 3) {
    sp_input <- c(1, 2, 3)
  }
}
if (coordinates == "delay") {
  sp_input <- 1
}
# full time series length
t_max_full <- 500
# final time series length
t_max <- 100
# fixed value for the inverse length scale of the control parameter
# (set to NA for phi to be estimated from data)
fixed_phi_control <- NA

# generate time series while gradually increasing control ------------------------------ 
# simulate model
ts <- func(t_max = t_max_full, x = x_init, p = unlist(parms), control_change = control_change,
           control_noise = control_noise)
ts <- ts[ts$control > 0.42 & ts$control <= 0.48, ]
ts$time <- 1:nrow(ts)
# determine tipping point through change-point analysis
sp_breakpoint <- c()
ss_breakpoint <- c()
breakpoint_results <- ts_breakpoint(ts = ts[ , c("time", "x1")], 
                                    n_breaks = 100,
                                    n_segments = 1)
sp_breakpoint <- breakpoint_results[[1]]
ss_breakpoint <- breakpoint_results[[2]]
# compute rolling variance (Early Warning Signal)
window_times <- ts$time[!is.na(match(round(ts$control, 4), round(control_gp, 4)))]
control_noiseless <- seq(parms$p, length = nrow(ts), by = control_change)
window_size <- 30
var_ts <- c()
window_end <- c()
for (i in 1:(nrow(ts)-window_size)) {
  ts_curr <- ts[i:(i+window_size), ]
  window_end[i] <- tail(ts_curr$time, 1)
  var_ts[i] <- var(ts_curr$x1)
}
plot_df <- ts
plot_df$var <- c(rep(NA, window_end[1]-1), var_ts)
plot_df$var[plot_df$control > tail(control_gp, 1)] <- NA
# plot time series and rolling variance
fig <- ggplot(data = plot_df) +
  geom_line(aes(x = control, y = x1), size = 0.9, color = "#919191") +
  geom_line(aes(x = control, y = var * 30), size = 0.9, color = "#2171B5") +
  geom_vline(xintercept = control_gp, linetype = "dashed", size = 0.8) +
  geom_vline(xintercept = ts$control[sp_breakpoint], size = 0.8) +
  xlab(label = control_name) +
  scale_x_continuous(limits = c(min(control_true), max(control_true)),
                     breaks = seq(control_true[1], control_true[length(control_true)], by = 0.01)) +
  scale_y_continuous(name = "Abundance (species 1)", 
                     sec.axis = sec_axis(~./30, name = "Rolling variance")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y.right = element_text(color = "#2171B5"),
        axis.text.y.right = element_text(color = "#2171B5"),
        axis.line.y.right = element_line(color = "#2171B5"), 
        axis.ticks.y.right = element_line(color = "#2171B5"))
if (save_plot) {
  ggsave(paste("figs/fig3/fig3_", func_name, "_ts.pdf", sep = ""),
         fig, width = 32, height = 8, units = "cm")
}
# plot control parameter over time
fig <- ggplot(data = ts, aes(x = time, y = control)) +
  geom_line(size = 0.8, color = "#919191") +
  ylab(label = control_name) +
  xlab(label = "Time") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

# fit GP-EDM to data to determine best embedding dimension ------------------------------ 
# create lagged data frame for GP-EDM
training_df <- ts[ts$control <= control_gp[1], ]
# log abundances
if (log_abund) {
  training_df$x1 <- log(training_df$x1)
  training_df$x2 <- log(training_df$x2)
}
# for each value of E, fit GP-EDM and evaluate prediction accuracy
if (coordinates == "delay") {
  R2 <- c()
  for (i in 1:length(E_values)) {
    # input names
    inputs <- paste(names(training_df)[sp_input + 1], E_values[1:i], sep = "_")
    if (control_lags) {
      inputs <- c(inputs, paste("control", E_values[1:i], sep = "_"))
    } else {
      inputs <- c(inputs, "control_1")
    }
    # create lags
    training_lags <- as.data.frame(makelags(data = training_df, y = colnames(training_df)[-1], 
                                            time = "time", E = max(E_values), tau = tau))
    full_training_lags <- cbind(training_df, training_lags)
    # fit GP-EDM on training data
    if (sp_input == 1) {
      GP_fit_sp1 <- fitGP(data = full_training_lags, y = "x1", 
                          x = inputs, scaling = "local",
                          fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_control, NA, NA),
                          predictmethod = "loo")
    }
    if (sp_input == 2) {
      GP_fit_sp2 <- fitGP(data = full_training_lags, y = "x2", 
                          x = inputs, scaling = "local",
                          fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_control, NA, NA),
                          predictmethod = "loo")
    }
    # compute R2
    if (sp_input == 1) {
      full_training_lags$obs <- GP_fit_sp1$outsampresults$obs
      full_training_lags$predmean <- GP_fit_sp1$outsampresults$predmean
    }
    if (sp_input == 2) {
      full_training_lags$obs <- GP_fit_sp2$outsampresults$obs
      full_training_lags$predmean <- GP_fit_sp2$outsampresults$predmean
    }
    R2_df <- full_training_lags
    R2_df <- na.omit(R2_df)
    if (log_abund) {
      R2_df$obs <- exp(R2_df$obs)
      R2_df$predmean <- exp(R2_df$predmean)
    }
    R2[i] <- 1 - (sum((R2_df$obs - R2_df$predmean)^2) / sum((R2_df$obs - mean(R2_df$obs))^2))
  }
  # select best embedding dimension
  E <- which.max(R2)
}

# reconstruct bifurcation diagrams by sequentially increasing time series ------------------------------ 
# input names
if (coordinates == "delay") {
  inputs <- paste(names(training_df)[sp_input + 1], 1:E, sep = "_")
  if (control_lags) {
    inputs <- c(inputs, paste("control", 1:E, sep = "_"))
  } else {
    inputs <- c(inputs, "control_1")
  }
}
if (coordinates == "native") {
  inputs <- paste(rep(names(training_df)[sp_input + 1], each = E), 1:E, sep = "_")
  inputs <- c(inputs, "control_1")
}
# to store plots
plot_list <- list()
# loop over successively larger training data sets
for (i in 1:length(control_gp)) {
  print(i)
  # create lagged data frame for GP-EDM
  training_df <- ts[ts$control <= control_gp[i], ]
  # log abundances
  if (log_abund) {
    training_df$x1 <- log(training_df$x1)
    training_df$x2 <- log(training_df$x2)
  }
  # create lags
  training_lags <- as.data.frame(makelags(data = training_df, y = colnames(training_df)[-1], 
                                          time = "time", E = E, tau = tau))
  full_training_lags <- cbind(training_df, training_lags)
  # fit GP-EDM on training data
  if (sp_input == 1) {
    GP_fit_sp1 <- fitGP(data = full_training_lags, y = "x1", 
                        x = inputs, scaling = "local",
                        fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_control, NA, NA))
  }
  if (sp_input == 2) {
    GP_fit_sp2 <- fitGP(data = full_training_lags, y = "x2", 
                        x = inputs, scaling = "local",
                        fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_control, NA, NA))
  }
  
  # perform forecasts to reconstruct bifurcation diagram ------------------------------ 
  full_test_set <- data.frame()
  # starting abundance values
  x_init <- head(training_df, E + 1)[ , -c(1, ncol(training_df))]
  # loop over all control parameter values
  for (j in 1:length(control_true)) {
    print(j)
    # build test set to make predictions
    if (n_sp == 2) {
      test_df <- data.frame(time = 1:(E + 1),
                            x1 = as.numeric(x_init[ , 1]),
                            x2 = as.numeric(x_init[ , 2]),
                            control = control_true[j])
    }
    # iterate forecasts for t_max_full time steps
    for (k in 1:t_max_full) {
      # build data frame to make forecast
      test_lags <- makelags(data = test_df, y = colnames(training_df)[-1],
                            time = "time", E = E, tau = tau, forecast = TRUE)
      # make forecast
      if (coordinates == "delay") {
        if (sp_input == 1) {
          forecast_sp1 <- predict(GP_fit_sp1, newdata = test_lags)$outsampresults$predmean
          test_df[k + E + 1, "x1"] <- forecast_sp1
        }
        if (sp_input == 2) {
          forecast_sp2 <- predict(GP_fit_sp2, newdata = test_lags)$outsampresults$predmean
          test_df[k + E + 1, "x2"] <- forecast_sp2
        }
      }
      if (coordinates == "native") {
        if (n_sp == 2) {
          forecast_sp1 <- predict(GP_fit_sp1, newdata = test_lags)$outsampresults$predmean
          test_df[k + E + 1, "x1"] <- forecast_sp1
          forecast_sp2 <- predict(GP_fit_sp2, newdata = test_lags)$outsampresults$predmean
          test_df[k + E + 1, "x2"] <- forecast_sp2
        }
      }
      # update time and control level in test set
      test_df[k + E + 1, "time"] <- test_df[k + E, "time"] + 1
      test_df[k + E + 1, "control"] <- control_true[j]
    }
    # plot time series
    fig <- ggplot(data = test_df, aes(x = time, y = x1)) +
      geom_line(size = 0.8, color = "#919191") +
      xlab(label = "Time") +
      ylab(label = "Abundance (species 1)") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(size = 1),
            axis.title = element_text(size = 18),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14))
    # update starting abundance values
    x_init <- tail(test_df, E + 1)[ , -c(1, ncol(test_df))]
    # extract tail from test set
    test_df <- tail(test_df, t_max)
    # merge with full data frame
    full_test_set <- rbind(full_test_set, test_df)
  }
  # exponentiate abundances
  if (log_abund) {
    full_test_set$x1 <- exp(full_test_set$x1)
    full_test_set$x2 <- exp(full_test_set$x2)
  }
  # create results data frame
  results_no_noise$type <- "true"
  full_test_set$type <- "predicted"
  plot_df <- full_test_set
  plot_ts <- ts
  plot_ts$type <- "unobserved"
  plot_ts$type[plot_ts$control <= control_gp[i]] <- "observed"
  # plot predicted bifurcation diagram together with data points
  fig <- ggplot() +
    geom_point(data = subset(plot_ts, type == "unobserved"), 
               aes(x = control, y = x1), 
               size = 2.5, alpha = 0.8, color = "#DEDEDE") +
    geom_point(data = subset(plot_ts, type == "observed"), 
               aes(x = control, y = x1), 
               size = 2.5, alpha = 0.8, color = "#454545") +
    geom_point(data = plot_df, aes(x = control, y = x1),
               size = 0.9, color = "#CB181D") +
    geom_vline(xintercept = control_gp[i], linetype = "dashed", size = 1) +
    geom_vline(xintercept = plot_ts$control[sp_breakpoint], size = 1) +
    xlab(label = control_name) +
    ylab(label = "Abundance (species 1)") +
    theme_bw() +
    scale_x_continuous(limits = c(min(control_true), max(control_true)),
                       breaks = seq(control_true[1], control_true[length(control_true)], by = 0.02)) +
    scale_y_continuous(limits = c(0, 0.7),
                       labels = scales::number_format(accuracy = 0.1)) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          legend.position = "none")
  # save plot
  if (save_plot) {
    ggsave(paste("figs/fig3/fig3_", func_name, "_predicted_bif_diag_", paste(sp_input, collapse = "_"), "_E_", E, 
                 "_control_level_", i, ".pdf", sep = ""), 
           fig, width = 13, height = 10, units = "cm")
  }
}
