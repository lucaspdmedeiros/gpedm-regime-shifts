# For the continuous-time Hastings-Powell food chain model, generate time
# series for a few values of the control parameter. Then, train a GP-EDM
# model and extrapolate to unseen values of the control parameter to
# reconstruct the bifurcation diagram.

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/jensen_shannon.R")
if (!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if (!require(viridis)) {install.packages("viridis"); library(viridis)}
if (!require(GPEDM)) {install.packages("GPEDM"); library(GPEDM)}
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(demodelr)) {install.packages("demodelr"); library(demodelr)}

# set up ------------------------------ 
# set seed to make results reproducible
set.seed(1)
# whether to save plots
save_plots <- TRUE
# whether to use native or delay coordinates
coordinates <- "delay"
# control parameter, initial conditions, and GP-EDM settings for selected model
n_sp <- 3
deterministic_part <- c(dx1 ~ x1 * (1 - x1) - ((a1 * x1 * x2) / (1 + p * x1)), 
                        dx2 ~ ((a1 * x1 * x2) / (1 + p * x1)) - ((a2 * x2 * x3) / (1 + b2 * x2)) - d1 * x2,
                        dx3 ~ ((a2 * x2 * x3) / (1 + b2 * x2)) - d2 * x3)
stochastic_part <-  c(dx1 ~ 1, dx2 ~ 1, dx3 ~ 1)
x_init <- c(x1 = 0.4386305, x2 = 0.2785410, x3 = 9.520890)
parms <- list(a1 = 5.0, p = NA, a2 = 0.1, 
              b2 = 2.0, d1 = 0.4, d2 = 0.01)
s <- 0.000001
control_true <- seq(2, 3, by = 0.005)
control_gp <- c(2, 2.3, 2.6, 3)
control_name <- "Half-saturation coefficient of the consumer"
log_abund <- FALSE
time_step <- 0.01 
n_points <- 100000
sampling_freq <- 500
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
t_max_train <- (n_points / 2) / sampling_freq
# fixed value for the inverse length scale of the control parameter
# (set to NA for phi to be estimated from data)
fixed_phi_control <- NA

# build time series from model ------------------------------ 
# data frames to save results
results_no_noise <- data.frame()
results_noise <- data.frame()
results_max_min <- data.frame()
for (i in 1:length(control_true)) {
  print(i)
  # parameters without noise
  parms$p <- control_true[i]
  # simulate model without noise
  curr_results <- as.data.frame(euler_stochastic(deterministic_rate = deterministic_part,
                                                 stochastic_rate = stochastic_part,
                                                 initial_condition = x_init,
                                                 parameters = unlist(parms),
                                                 deltaT = time_step,
                                                 n_steps = n_points,
                                                 D = 0))
  # remove transient
  curr_results <- tail(curr_results, nrow(curr_results) / 2)
  # obtain local maxima and minima for chosen species
  x <- curr_results[ , sp_plot + 1]
  max_min_locations <- which(abs(diff(sign(diff(x)))) > 0) + 1
  results_max_min_curr <- curr_results[max_min_locations, ]
  results_max_min_curr$control <- control_true[i]
  results_max_min <- rbind(results_max_min, results_max_min_curr)
  # sample points according to sampling frequency
  curr_results <- curr_results[seq(1, nrow(curr_results), by = sampling_freq), ]
  # save results without noise
  curr_results$control <- control_true[i]
  results_no_noise <- rbind(results_no_noise, curr_results)
  # simulate model with noise
  curr_results <- as.data.frame(euler_stochastic(deterministic_rate = deterministic_part,
                                                 stochastic_rate = stochastic_part,
                                                 initial_condition = x_init,
                                                 parameters = unlist(parms),
                                                 deltaT = time_step,
                                                 n_steps = n_points,
                                                 D = s))
  # remove transient
  curr_results <- tail(curr_results, nrow(curr_results) / 2)
  # sample points according to sampling frequency
  curr_results <- curr_results[seq(1, nrow(curr_results), by = sampling_freq), ]
  # save results with noise
  curr_results$control <- control_true[i]
  results_noise <- rbind(results_noise, curr_results)
}
# change abundances to log scale
if (log_abund) {
  results_noise[ , -c(1, ncol(results_noise))] <- 
    log(results_noise[ , -c(1, ncol(results_noise))])
  results_no_noise[ , -c(1, ncol(results_no_noise))] <- 
    log(results_no_noise[ , -c(1, ncol(results_no_noise))])
}
# plot true bifurcation diagram
results_max_min$abundance <- results_max_min[ , sp_plot + 1]
fig <- ggplot(data = results_max_min, aes(x = control, y = abundance)) +
  geom_point(size = 1, color = "#929292") +
  geom_vline(xintercept = control_gp, linetype = "dashed", size = 1) +
  xlab(label = control_name) +
  ylab(label = paste("Abundance (species ", sp_plot, ")", sep = "")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))
# change data frame column names
names(results_no_noise)[1] <- "time"
names(results_noise)[1] <- "time"
names(results_max_min)[1] <- "time"

# plot time series used to train the GP-EDM model ------------------------------ 
# extract data frames
tol <- 10^(-10)
training_df <- subset(results_noise, (abs(results_noise$control - control_gp[1]) < tol) | 
                        (abs(results_noise$control - control_gp[2]) < tol) |
                        (abs(results_noise$control - control_gp[3]) < tol) |
                        (abs(results_noise$control - control_gp[4]) < tol))
# plot time series used to train the GP-EDM model
if (n_sp == 1) {
  plot_df <- gather(training_df, "species", "density", x1)
}
if (n_sp == 2) {
  plot_df <- gather(training_df, "species", "density", x1:x2)
}
if (n_sp == 3) {
  plot_df <- gather(training_df, "species", "density", x1:x3)
}
fig <- ggplot(data = plot_df, aes(x = time, y = density, group = species, color = species)) +
  geom_line(size = 1) +
  facet_wrap(~control, labeller = "label_both", nrow = length(control_gp)) +
  xlab(label = "Time") +
  ylab(label = "Population abundance") +
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
    }
    if (sp_input == 2) {
      GP_fit_sp2 <- fitGP(data = full_training_lags, y = "x2", 
                          x = inputs, scaling = "local",
                          fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_control, NA, NA),
                          predictmethod = "loo")
    }
    if (sp_input == 3) {
      GP_fit_sp3 <- fitGP(data = full_training_lags, y = "x3", 
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
  # compute prediction errors
  plot_df$residuals <- plot_df$obs - plot_df$predmean
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
  # compute prediction errors
  plot_df$residuals <- plot_df$obs - plot_df$predmean
  # compute R2
  R2_df <- full_training_lags
  R2_df <- na.omit(R2_df)
  if (log_abund) {
    R2_df$obs <- exp(R2_df$obs)
    R2_df$predmean <- exp(R2_df$predmean)
  }
  R2 <- 1 - (sum((R2_df$obs - R2_df$predmean)^2) / sum((R2_df$obs - mean(R2_df$obs))^2))
}
# plot observations and predictions
plot_df$time <- rep(1:t_max_train, times = length(control_gp))
plot_df$p <- plot_df$control
fig <- ggplot(data = plot_df) +
  geom_line(aes(x = time, y = obs), size = 1, color = "black") +
  geom_line(aes(x = time, y = predmean), size = 1, color = "#CB181D") +
  xlab(label = "Time") +
  ylab(label = paste("Abundance (species ", sp_plot, ")", sep = "")) +
  facet_wrap(~p,
             labeller = "label_both", nrow = 5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.2),
        strip.text = element_text(size = 16),
        strip.background = element_rect(fill = "white", size = 1.2),
        title = element_text(size = 14),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14))

# residual plots ------------------------------ 
if (coordinates == "delay") {
  # plot residual autocorrelation
  residuals_df <- as.data.frame(makelags(data = plot_df, y = "residuals", 
                                         time = "time", E = max(E_values), tau = tau))
  auto_cor <- apply(X = residuals_df, MARGIN = 2, FUN = function(x, y) cor(x, y, use = "complete.obs"), 
                    y = plot_df$residuals)
  auto_cor_df <- data.frame(lag = 1:max(E_values), auto_cor = auto_cor)
  fig <- ggplot(data = auto_cor_df) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
    geom_bar(aes(x = lag, y = auto_cor), stat="identity", size = 1, width = 0.5,
             color = "black", fill = "gray70") +
    scale_x_continuous(n.breaks = max(E_values)) +
    xlab(label = "Lag") +
    ylab(label = "Residual\nautocorrelation") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.2),
          strip.text = element_text(size = 16),
          strip.background = element_rect(fill = "white", size = 1.2),
          title = element_text(size = 14),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18))
  # plot residual distribution
  fig <- ggplot(data = plot_df) +
    geom_density(aes(x = residuals), size = 1, color = "black", fill = "gray70") +
    xlab(label = "Residuals") +
    ylab(label = "Density") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.2),
          strip.text = element_text(size = 16),
          strip.background = element_rect(fill = "white", size = 1.2),
          title = element_text(size = 14),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18))
  # plot residuals vs state variables
  residuals_df <- plot_df[ , c(inputs[1:E], "residuals")]
  names(residuals_df) <- c(c(paste("lag", 1:E, sep = " "), "residuals"))
  residuals_df <- gather(residuals_df, "lag", "value", names(residuals_df)[-ncol(residuals_df)])
  fig <- ggplot(data = residuals_df) +
    geom_point(aes(x = residuals, y = value), size = 2) +
    facet_wrap(~lag, ncol = ceiling(E/2)) +
    xlab(label = "Residuals") +
    ylab(label = "Abundance") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.2),
          strip.text = element_text(size = 18),
          strip.background = element_rect(fill = "white", size = 1.2),
          title = element_text(size = 14),
          axis.title = element_text(size = 22),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 14),
          legend.key.size = unit(0.6, 'cm'))
}

# make predictions for the whole range of control parameter values ------------------------------ 
full_test_set <- data.frame()
results_max_min_predictions <- data.frame()
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
  # obtain local maxima and minima for chosen species
  x <- test_df[ , sp_plot + 1]
  max_min_locations <- which(abs(diff(sign(diff(x)))) > 0) + 1
  results_max_min_curr <- test_df[max_min_locations, ]
  results_max_min_predictions <- rbind(results_max_min_predictions, results_max_min_curr)
  # merge with full data frame
  full_test_set <- rbind(full_test_set, test_df)
}

# compute Jensen-Shannon divergence for each value of control parameter ------------------------------ 
# merge data frames with true and predicted bifurcation diagrams
results_max_min$type <- "true"
results_max_min_predictions$type <- "predicted"
# set negative predictions to zero
results_max_min_predictions[, sp_plot + 1][results_max_min_predictions[, sp_plot + 1] < 0] <- 0
# start and end points for bins
start_bin <- min(c(results_max_min[ , sp_plot + 1], results_max_min_predictions[ , sp_plot + 1]))
end_bin <- max(c(results_max_min[ , sp_plot + 1], results_max_min_predictions[ , sp_plot + 1]))
bin_number <- list(10, 11, 12, 13, 14, 15)
# loop over control values and compute Jensen-Shannon divergence
jsd <- c()
for (i in 1:length(control_true)) {
  # samples from probability distributions
  x <- results_max_min[results_max_min$control == control_true[i], sp_plot + 1]
  y <- results_max_min_predictions[results_max_min_predictions$control == control_true[i], sp_plot + 1]
  # establish bin size according to the Freedman–Diaconis rule
  bin_size <- 2 * (quantile(c(x, y))[4] - quantile(c(x, y))[2]) / 
    length(c(x, y))^(1/3)
  # compute average Jensen-Shannon divergence for a range of bin numbers
  jsd[i] <- mean(sapply(bin_number, jensen_shannon, 
                        x = x, y = y, start_bin = start_bin, end_bin = end_bin))
}
# plot of Jensen-Shannon divergence for different values of control parameter
plot_df <- data.frame(control = control_true,
                      jsd = jsd)
fig <- ggplot(data = plot_df, aes(x = control, y = jsd)) +
  geom_line(size = 0.8) +
  geom_vline(xintercept = control_gp, linetype = "dashed", size = 0.8) +
  scale_y_continuous(limits = c(0, 1)) +
  xlab(label = control_name) +
  ylab(label = "Jensen-Shannon divergence") +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'))
if (save_plots) {
  ggsave(paste("figs/fig2/fig2_hastings_powell_jsd_", coordinates, 
               "_sp_input_", paste(sp_input, collapse = "_"), "_control_phi_", fixed_phi_control,
               "_E_", E, "_", log_abund_name, ".pdf", sep = ""), 
         fig, width = 14, height = 6, units = "cm")
  
}

# plot true and predicted bifurcation diagrams ------------------------------ 
results_max_min_predictions$abundance <- results_max_min_predictions[ , sp_plot + 1]
plot_df <- rbind(results_max_min, results_max_min_predictions)
if (log_abund) {
  plot_df[ , -c(1, ncol(plot_df)-1, ncol(plot_df))] <- 
    exp(plot_df[ , -c(1, ncol(plot_df)-1, ncol(plot_df))])
}
# plot true and predicted together
fig <- ggplot(data = plot_df, aes(x = control, y = abundance, color = type)) +
  geom_point(size = 0.8) +
  geom_vline(xintercept = control_gp, linetype = "dashed", size = 0.8) +
  scale_color_manual(values = c("#CB181D", "#919191")) +
  scale_x_continuous(limits = c(min(control_true), max(control_true))) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  xlab(label = control_name) +
  ylab(label = paste("Abundance (species ", sp_plot, ")", sep = "")) +
  theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'))
if (save_plots) {
  ggsave(paste("figs/fig2/fig2_hastings_powell_true_vs_predicted_bif_diag_", coordinates, 
               "_sp_input_", paste(sp_input, collapse = "_"), "_control_phi_", fixed_phi_control,
               "_E_", E, "_", log_abund_name, ".pdf", sep = ""), 
         fig, width = 14, height = 10, units = "cm")
}
