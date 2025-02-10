# For the two-species predator-prey model, generate time
# series for a few values of the control parameter. Then, train a GP-EDM
# model and extrapolate to unseen values of the control parameter to
# reconstruct the bifurcation diagram. Then, plot the power spectrum 
# and Lyapunov Exponent of true and predicted time series for certain 
# control values.

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/predator_prey_2sp_map.R")
source("code/functions/power_spectrum.R")
source("code/functions/le_qr_decomp.R")
if (!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if (!require(ggExtra)) {install.packages("ggExtra"); library(ggExtra)}
if (!require(viridis)) {install.packages("viridis"); library(viridis)}
if (!require(GPEDM)) {install.packages("GPEDM"); library(GPEDM)}
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}

# set up ------------------------------ 
# whether to save plots
save_plots <- TRUE
# model to use
func_name <- "predator_prey_2sp_map"
# whether to use native or delay coordinates
coordinates <- "delay"
# rate of change of control parameter
control_change <- 0
# noise level in control values
control_noise <- 0
# control parameter, initial conditions, and GP-EDM settings for selected model
if (func_name == "predator_prey_2sp_map") {
  set.seed(10)
  n_sp <- 2
  func <- predator_prey_2sp_map
  x_init <- c(x1 = 1, x2 = 1)
  control_true <- seq(2.4, 3.4, by = 0.005)
  control_gp <- c(2.4, 2.6, 3.15, 3.4)
  control_test <- seq(2.7, 3.4, by = 0.1)
  parms <- list(r = 2.5, k = 4, p = NA, h = 0.1, e = 0.5, m = 2, 
                s_mean = NA, s_sd = NA)
  s <- 0.02
  control_name <- "Attack rate"
  log_abund <- FALSE
}
# to name files
if (log_abund) {
  log_abund_name <- "yes_log_abund"
} else {
  log_abund_name <- "no_log_abund"
}
# grid of embedding dimension values and time lag value
E_values <- 1:9
tau <- 1
# define species to use in plots and analyses
if (coordinates == "native") {
  E <- 1
  sp_input <- 1:n_sp
  sp_plot <- 1
}
if (coordinates == "delay") {
  sp_input <- 1
  sp_plot <- sp_input
}
# full time series length
t_max_full <- 500
# final time series length
t_max <- 100
# time series length to train GP-EDM
t_max_train <- 50
# fixed value for the inverse length scale of the control parameter
# (set to NA for phi to be estimated from data)
fixed_phi_control <- NA
# data frames to save results
results_no_noise <- data.frame()
results_noise <- data.frame()

# build time series from model ------------------------------ 
# data frame to store power spectrum of each time series
power_spectrum_df <- data.frame()
for (i in 1:length(control_true)) {
  # parameters without noise
  parms$p <- control_true[i]
  parms$s_mean <- 0
  parms$s_sd <- 0
  # simulate model without noise
  curr_results <- func(t_max = t_max_full, x = x_init, p = unlist(parms), 
                       control_change = control_change, control_noise = control_noise)
  # compute time series power spectrum
  ts_ps <- tail(curr_results, t_max)[ , c(1, sp_plot+1)]
  power_spectrum_curr <- power_spectrum(ts = ts_ps, lam = 0.01, 
                                        scale = TRUE, trim = TRUE, plot = FALSE)[[1]]
  power_spectrum_curr$control <- control_true[i]
  power_spectrum_curr$type <- "true"
  power_spectrum_df <- rbind(power_spectrum_df, power_spectrum_curr)
  # save results without noise
  curr_results$control <- control_true[i]
  results_no_noise <- rbind(results_no_noise, tail(curr_results, t_max))
  # parameters with noise
  parms$p <- control_true[i]
  parms$s_mean <- -(s^2)/2
  parms$s_sd <- s
  # simulate model with noise
  curr_results <- func(t_max = t_max_full, x = x_init, p = unlist(parms), 
                       control_change = control_change, control_noise = control_noise)
  # save results with noise
  curr_results$control <- control_true[i]
  results_noise <- rbind(results_noise, tail(curr_results, t_max_train))
}
# change abundances to log scale
if (log_abund) {
  results_noise[ , -c(1, ncol(results_noise))] <- 
    log(results_noise[ , -c(1, ncol(results_noise))])
  results_no_noise[ , -c(1, ncol(results_no_noise))] <- 
    log(results_no_noise[ , -c(1, ncol(results_no_noise))])
}

# plot time series used to train the GP-EDM model ------------------------------ 
# extract data frames
tol <- 10^(-10)
training_df <- subset(results_noise, (abs(results_noise$control - control_gp[1]) < tol) | 
                        (abs(results_noise$control - control_gp[2]) < tol) |
                        (abs(results_noise$control - control_gp[3]) < tol) |
                        (abs(results_noise$control - control_gp[4]) < tol))
# create data frame with all lags
training_list <- split(training_df, training_df$control)
full_training_lags <- data.frame()
for (i in 1:length(training_list)) {
  training_lags <- as.data.frame(makelags(data = training_list[[i]], y = colnames(training_df)[sp_input+1], 
                                          time = "time", E = max(E_values), tau = tau))
  training_lags$control_1 <- c(NA, training_list[[i]]$control[2:nrow(training_list[[i]])])
  training_lags_full <- cbind(training_list[[i]], training_lags)
  full_training_lags <- rbind(full_training_lags, training_lags_full)
}

# fit GP-EDM to data to determine best embedding dimension ------------------------------ 
# for each value of E, fit GP-EDM and evaluate prediction accuracy
if (coordinates == "delay") {
  R2 <- c()
  phi_control <- c()
  for (i in 1:length(E_values)) {
    # input names
    inputs <- paste(names(training_df)[sp_input + 1], E_values[1:i], sep = "_")
    inputs <- c(inputs, "control_1")
    # fit GP-EDM on training data
    if (sp_input == 1) {
      GP_fit_sp1 <- fitGP(data = full_training_lags, y = "x1", 
                          x = inputs, scaling = "local",
                          fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_control, NA, NA),
                          predictmethod = "loo")
      phi_control[i] <- GP_fit_sp1$pars[E_values[i]+1]
    }
    if (sp_input == 2) {
      GP_fit_sp2 <- fitGP(data = full_training_lags, y = "x2", 
                          x = inputs, scaling = "local",
                          fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_control, NA, NA),
                          predictmethod = "loo")
      phi_control[i] <- GP_fit_sp2$pars[E_values[i]+1]
    }
    if (sp_input == 3) {
      GP_fit_sp3 <- fitGP(data = full_training_lags, y = "x3", 
                          x = inputs, scaling = "local",
                          fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_control, NA, NA),
                          predictmethod = "loo")
      phi_control[i] <- GP_fit_sp3$pars[E_values[i]+1]
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
    if (sp_input == 3) {
      full_training_lags$obs <- GP_fit_sp3$outsampresults$obs
      full_training_lags$predmean <- GP_fit_sp3$outsampresults$predmean
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

# fit GP-EDM to data using best embedding dimension ------------------------------ 
# input names
if (coordinates == "delay") {
  inputs <- paste(names(training_df)[sp_input + 1], 1:E, sep = "_")
  inputs <- c(inputs, "control_1")
}
if (coordinates == "native") {
  inputs <- paste(rep(names(training_df)[sp_input + 1], each = E), 1:E, sep = "_")
  inputs <- c(inputs, "control_1")
}
# fit GP-EDM in delay coordinates
if (coordinates == "delay") {
  if (sp_input == 1) {
    GP_fit_sp1 <- fitGP(data = full_training_lags, y = "x1", 
                        x = inputs, scaling = "local",
                        fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_control, NA, NA),
                        predictmethod = "loo")
    full_training_lags$obs <- GP_fit_sp1$outsampresults$obs
    full_training_lags$predmean <- GP_fit_sp1$outsampresults$predmean
  }
  if (sp_input == 2) {
    GP_fit_sp2 <- fitGP(data = full_training_lags, y = "x2", 
                        x = inputs, scaling = "local",
                        fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_control, NA, NA),
                        predictmethod = "loo")
    full_training_lags$obs <- GP_fit_sp2$outsampresults$obs
    full_training_lags$predmean <- GP_fit_sp2$outsampresults$predmean
  }
  if (sp_input == 3) {
    GP_fit_sp3 <- fitGP(data = full_training_lags, y = "x3", 
                        x = inputs, scaling = "local",
                        fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_control, NA, NA),
                        predictmethod = "loo")
    full_training_lags$obs <- GP_fit_sp3$outsampresults$obs
    full_training_lags$predmean <- GP_fit_sp3$outsampresults$predmean
  }
  plot_df <- full_training_lags
  if (log_abund) {
    plot_df$obs <- exp(plot_df$obs)
    plot_df$predmean <- exp(plot_df$predmean)
  }
}
# fit GP-EDM in native coordinates
if (coordinates == "native") {
  if (n_sp == 1) {
    GP_fit_sp1 <- fitGP(data = full_training_lags, y = "x1", 
                        x = inputs, scaling = "local",
                        fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_control, NA, NA),
                        predictmethod = "loo")
  }
  if (n_sp == 2) {
    GP_fit_sp1 <- fitGP(data = full_training_lags, y = "x1", 
                        x = inputs, scaling = "local",
                        fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_control, NA, NA),
                        predictmethod = "loo")
    GP_fit_sp2 <- fitGP(data = full_training_lags, y = "x2", 
                        x = inputs, scaling = "local",
                        fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_control, NA, NA),
                        predictmethod = "loo")
  }
  if (n_sp == 3) {
    GP_fit_sp1 <- fitGP(data = full_training_lags, y = "x1", 
                        x = inputs, scaling = "local",
                        fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_control, NA, NA),
                        predictmethod = "loo")
    GP_fit_sp2 <- fitGP(data = full_training_lags, y = "x2", 
                        x = inputs, scaling = "local",
                        fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_control, NA, NA),
                        predictmethod = "loo")
    GP_fit_sp3 <- fitGP(data = full_training_lags, y = "x3", 
                        x = inputs, scaling = "local",
                        fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_control, NA, NA),
                        predictmethod = "loo")
  }
  if (sp_plot == 1) {
    full_training_lags$obs <- GP_fit_sp1$outsampresults$obs
    full_training_lags$predmean <- GP_fit_sp1$outsampresults$predmean
  }
  if (sp_plot == 2) {
    full_training_lags$obs <- GP_fit_sp2$outsampresults$obs
    full_training_lags$predmean <- GP_fit_sp2$outsampresults$predmean
  }
  if (sp_plot == 3) {
    full_training_lags$obs <- GP_fit_sp3$outsampresults$obs
    full_training_lags$predmean <- GP_fit_sp3$outsampresults$predmean
  }
  plot_df <- full_training_lags
  if (log_abund) {
    plot_df$obs <- exp(plot_df$obs)
    plot_df$predmean <- exp(plot_df$predmean)
  }
  # compute R2
  R2_df <- full_training_lags
  R2_df <- na.omit(R2_df)
  if (log_abund) {
    R2_df$obs <- exp(R2_df$obs)
    R2_df$predmean <- exp(R2_df$predmean)
  }
  R2 <- 1 - (sum((R2_df$obs - R2_df$predmean)^2) / sum((R2_df$obs - mean(R2_df$obs))^2))
}

# make predictions for the whole range of control parameter values ------------------------------ 
full_test_set <- data.frame()
# loop over all control parameter values
for (i in 1:length(control_true)) {
  print(i)
  # choosing starting abundance values
  curr_control <- control_gp[which.min(sqrt((control_true[i] - control_gp)^2))]
  x_init <- tail(subset(training_df, abs(control - curr_control) < tol), E + 1)[ , -c(1, ncol(training_df))]
  # build test set to make predictions
  if (n_sp == 1) {
    test_df <- data.frame(time = 1:(E + 1),
                          x1 = as.numeric(x_init),
                          control = control_true[i])
  }
  if (n_sp == 2) {
    test_df <- data.frame(time = 1:(E + 1),
                          x1 = as.numeric(x_init[ , 1]),
                          x2 = as.numeric(x_init[ , 2]),
                          control = control_true[i])
  }
  if (n_sp == 3) {
    test_df <- data.frame(time = 1:(E + 1),
                          x1 = as.numeric(x_init[ , 1]),
                          x2 = as.numeric(x_init[ , 2]),
                          x3 = as.numeric(x_init[ , 3]),
                          control = control_true[i])
  }
  # iterate forecasts for t_max_full time steps
  for (j in 1:t_max_full) {
    # build data frame to make forecast
    test_lags <- makelags(data = test_df, y = colnames(training_df)[sp_input+1],
                          time = "time", E = E, tau = tau, forecast = TRUE)
    test_lags$control_1 <- control_true[i]
    # make forecast
    if (coordinates == "delay") {
      if (sp_input == 1) {
        forecast_sp1 <- predict(GP_fit_sp1, newdata = test_lags)$outsampresults$predmean
        test_df[j + E + 1, "x1"] <- forecast_sp1
      }
      if (sp_input == 2) {
        forecast_sp2 <- predict(GP_fit_sp2, newdata = test_lags)$outsampresults$predmean
        test_df[j + E + 1, "x2"] <- forecast_sp2
      }
      if (sp_input == 3) {
        forecast_sp3 <- predict(GP_fit_sp3, newdata = test_lags)$outsampresults$predmean
        test_df[j + E + 1, "x3"] <- forecast_sp3
      }
    }
    if (coordinates == "native") {
      if (n_sp == 1) {
        forecast_sp1 <- predict(GP_fit_sp1, newdata = test_lags)$outsampresults$predmean
        test_df[j + E + 1, "x1"] <- forecast_sp1
      }
      if (n_sp == 2) {
        forecast_sp1 <- predict(GP_fit_sp1, newdata = test_lags)$outsampresults$predmean
        test_df[j + E + 1, "x1"] <- forecast_sp1
        forecast_sp2 <- predict(GP_fit_sp2, newdata = test_lags)$outsampresults$predmean
        test_df[j + E + 1, "x2"] <- forecast_sp2
      }
      if (n_sp == 3) {
        forecast_sp1 <- predict(GP_fit_sp1, newdata = test_lags)$outsampresults$predmean
        test_df[j + E + 1, "x1"] <- forecast_sp1
        forecast_sp2 <- predict(GP_fit_sp2, newdata = test_lags)$outsampresults$predmean
        test_df[j + E + 1, "x2"] <- forecast_sp2
        forecast_sp3 <- predict(GP_fit_sp3, newdata = test_lags)$outsampresults$predmean
        test_df[j + E + 1, "x3"] <- forecast_sp3
      }
    }
    # update time and control level in test set
    test_df[j + E + 1, "time"] <- test_df[j + E, "time"] + 1
    test_df[j + E + 1, "control"] <- control_true[i]
  }
  # extract tail from test set
  test_df <- tail(test_df, t_max)
  # compute time series power spectrum
  ts_ps <- test_df[ , c(1, sp_plot+1)]
  power_spectrum_curr <- power_spectrum(ts = ts_ps, lam = 0.01, 
                                        scale = TRUE, trim = TRUE, plot = FALSE)[[1]]
  power_spectrum_curr$control <- control_true[i]
  power_spectrum_curr$type <- "predicted"
  power_spectrum_df <- rbind(power_spectrum_df, power_spectrum_curr)
  # merge with full data frame
  full_test_set <- rbind(full_test_set, test_df)
}

# plot true and predicted power spectrum for certain control levels ------------------------------ 
# data frame with chosen control levels
power_spectrum_plot_df <- subset(power_spectrum_df, (abs(power_spectrum_df$control - control_test[1]) < tol) | 
                                   (abs(power_spectrum_df$control - control_test[2]) < tol) |
                                   (abs(power_spectrum_df$control - control_test[3]) < tol) |
                                   (abs(power_spectrum_df$control - control_test[4]) < tol) |
                                   (abs(power_spectrum_df$control - control_test[5]) < tol) |
                                   (abs(power_spectrum_df$control - control_test[6]) < tol) |
                                   (abs(power_spectrum_df$control - control_test[7]) < tol) |
                                   (abs(power_spectrum_df$control - control_test[8]) < tol))
power_spectrum_plot_df$p <- power_spectrum_plot_df$control
# plot power spectrum 
fig <- ggplot(data = power_spectrum_plot_df, aes(x = frequency, y = power, color = type)) +
  geom_line(size = 1) +
  facet_wrap(~p, labeller = "label_both", nrow = length(control_test)/2) +
  xlab(label = "Frequency") +
  ylab(label = "Power") +
  scale_color_manual(values = c("#CB181D", "#919191")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white", size = 1),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'))
if (save_plots) {
  ggsave(paste("figs/fig_", func_name, "_true_vs_predicted_power_spectrum.pdf", sep = ""), 
         fig, width = 16, height = 16, units = "cm")
}

# compute Lyapunov exponents for true dynamics for certain control levels ------------------------------ 
# merge data frames with true and predicted bifurcation diagrams
results_no_noise$time <- rep(1:t_max, times = length(unique(control_true)))
results_no_noise$type <- "true"
full_test_set$time <- rep(1:t_max, times = length(unique(control_true)))
full_test_set$type <- "predicted"
plot_df <- rbind(results_no_noise, full_test_set)
if (log_abund) {
  plot_df[ , -c(1, ncol(plot_df)-1, ncol(plot_df))] <- 
    exp(plot_df[ , -c(1, ncol(plot_df)-1, ncol(plot_df))])
}
# take subset of data frame to plot
plot_df <- subset(plot_df, (abs(plot_df$control - control_test[1]) < tol) | 
                    (abs(plot_df$control - control_test[2]) < tol) |
                    (abs(plot_df$control - control_test[3]) < tol) |
                    (abs(plot_df$control - control_test[4]) < tol) |
                    (abs(plot_df$control - control_test[5]) < tol) |
                    (abs(plot_df$control - control_test[6]) < tol) |
                    (abs(plot_df$control - control_test[7]) < tol) |
                    (abs(plot_df$control - control_test[8]) < tol))
# plot true and predicted together
plot_df$x <- plot_df[ , sp_plot + 1]
plot_df$p <- as.factor(plot_df$control)
plot_df$type <- factor(plot_df$type, levels = c("predicted", "true"))
# function that creates Jacobian matrix
Jacobian <- function(x1, x2, r, k, p, h, e, m) {
  dx1_dx1 = exp(r * (1 - x1 / k) - p * x2 / (1 + p * h * x1)) - 
    x1 * (exp(r * (1 - x1 / k) - p * x2 / (1 + p * h * x1)) * (r * (1 / k) - p * x2 * (p * h) / (1 + p * h * x1)^2))
  dx1_dx2 = -(x1 * (exp(r * (1 - x1/k) - p * x2 / (1 + p * h * x1)) * (p / (1 + p * h * x1))))
  dx2_dx1 = x2 * (exp(e * p * x1 / (1 + p * h * x1) - m) * (e * p / (1 + p * h * x1) - e * p * x1 * p * h / (1 + p * h * x1)^2))
  dx2_dx2 = exp(e * p * x1/(1 + p * h * x1) - m)
  J <- matrix(c(dx1_dx1, dx1_dx2, dx2_dx1, dx2_dx2), nrow = 2, ncol = 2, byrow = TRUE)
  return(J)
}
# compute Lyapunov exponent for each value of the control parameter for true dynamics
max_le_true <- c()
for (i in 1:length(control_test)) {
  # current data frame
  curr_plot_df <- subset(plot_df, (abs(plot_df$control - control_test[i]) < tol) & 
                           type == "true")
  # Jacobian matrix for each state 
  J <- list()
  for (j in 1:nrow(curr_plot_df)) {
    J[[j]] <- Jacobian(x1 = as.numeric(curr_plot_df[j, "x1"]), x2 = as.numeric(curr_plot_df[j, "x2"]), 
                       r = 2.5, k = 4, p = as.numeric(curr_plot_df[j, "control"]), h = 0.1, e = 0.5, m = 2)
  }
  # compute largest Lyapunov exponent from sequence of Jacobian matrices
  max_le_true[i] <- max(le_qr_decomp(J))
}

# compute Lyapunov exponents for predicted dynamics for certain control levels ------------------------------ 
max_le_predicted <- c()
for (i in 1:length(control_test)) {
  # current data frame
  curr_plot_df <- subset(plot_df, (abs(plot_df$control - control_test[i]) < tol) & 
                           type == "predicted")
  # obtain gradient by fitting a GPEDM model to the predicted time series
  curr_plot_df_lags <- as.data.frame(makelags(data = curr_plot_df, y = "x1", 
                                              time = "time", E = E, tau = tau))
  curr_plot_df_lags <- cbind(curr_plot_df, curr_plot_df_lags)
  GP_fit_sp1_predicted <- fitGP(data = curr_plot_df_lags, y = "x1", 
                                x = inputs[-length(inputs)], scaling = "local",
                                fixedpars = GP_fit_sp1$pars[-c(E+1, length(GP_fit_sp1$pars))],
                                predictmethod = "loo", returnGPgrad = TRUE)
  grad_df <- GP_fit_sp1_predicted$GPgrad
  # remove NAs
  grad_df <- na.omit(grad_df)
  # create dummy matrix
  shifted_I <- cbind(diag(E-1), rep(0, E-1))
  # Jacobian matrix for each state 
  J <- list()
  for (j in 1:nrow(grad_df)) {
    # current GPEDM model gradient
    grad <- as.matrix(grad_df[j, 1:E])
    J[[j]] <- rbind(grad, shifted_I)
  }
  # compute largest Lyapunov exponent from sequence of Jacobian matrices
  max_le_predicted[i] <- max(le_qr_decomp(J))
}

# plot true and predicted Lyapunov exponents for certain control levels ------------------------------ 
# create data frame
le_plot_df <- data.frame(max_le_true = max_le_true, 
                         max_le_predicted = max_le_predicted,
                         p = control_test)
# plot
fig <- ggplot(data = le_plot_df, aes(x = max_le_true, y = max_le_predicted, color = p)) +
  geom_point(size = 3.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", size = 1) +
  geom_vline(xintercept = 0, size = 0.5) +
  geom_hline(yintercept = 0, size = 0.5) +
  scale_color_viridis() +
  xlab(label = "Largest Lyapunov exponent\ncomputed from true dynamics") +
  ylab(label = "Largest Lyapunov exponent\ncomputed from predicted dynamics") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white", size = 1),
        legend.position = "top",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 14),
        legend.key.height = unit(0.6, 'cm'),
        legend.key.width = unit(0.8, 'cm'))
if (save_plots) {
  ggsave(paste("figs/fig_", func_name, "_true_vs_predicted_lyapunov_exponent.pdf", sep = ""), 
         fig, width = 16, height = 16, units = "cm")
}

