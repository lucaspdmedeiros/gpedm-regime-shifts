# Code for Fig 1. Illustration of how GP-EDM works to predict dynamical 
# regimes for unobserved levels of the control parameter using the 2-species 
# predator-prey model.

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/logistic_1sp_map.R")
source("code/functions/competition_3sp_map.R")
source("code/functions/predator_prey_2sp_map.R")
source("code/functions/harvesting_1sp_map.R")
source("code/functions/harvesting_2sp_map.R")
if (!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if (!require(viridis)) {install.packages("viridis"); library(viridis)}
if (!require(GPEDM)) {install.packages("GPEDM"); library(GPEDM)}
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(latex2exp)) {install.packages("latex2exp"); library(latex2exp)}

# set up ------------------------------ 
# set seed to make results reproducible
set.seed(9)
# whether to save plots
save_plots <- TRUE
# model to use
func_name <- "predator_prey_2sp_map"
# rate of change of control parameter
control_change <- 0
# noise level in control values
control_noise <- 0
# parameter values
r <- 2.5
k <- 4
p <- c(2.4, 2.6, 3.15, 3.4)
h <- 0.1
e <- 0.5
m <- 2
# range of values for species abundances
x <- seq(1, 4, by = 0.05)
y <- 0.55
# create grid with control parameter and abundances
var_grid <- expand_grid(p, x, y)

# plot analytical conditional responses ------------------------------ 
# x1 and x2 at next time step
x_next <- var_grid$x * exp((r * (1 - var_grid$x / k) - (var_grid$p * var_grid$y) / (1 + var_grid$p * h * var_grid$x)))
y_next <- var_grid$y * exp(((e * var_grid$p * var_grid$x) / (1 + var_grid$p * h * var_grid$x) - m))
# put everything in a single data frame
plot_df <- cbind(var_grid, x_next, y_next)
plot_df$p <- as.factor(plot_df$p)
plot_df$y <- as.factor(plot_df$y)
# plot
fig <- ggplot() +
  geom_line(data = plot_df, aes(x = x, y = x_next, group = p, color = p), size = 1) + 
  scale_color_viridis(name = "Attack\nrate (p)", discrete = TRUE) +
  xlab(label = TeX("Current prey population ($x_{t}$)")) +
  ylab(label = TeX("Future prey population ($x_{t+1}$)")) +
  scale_y_continuous(breaks = c(2.0, 2.5)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.2),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "white", size = 1.2),
        title = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.7, 'cm'))
if (save_plots) {
  ggsave(paste("figs/fig1/fig1_", func_name, "_analytical_cond_resp.pdf", sep = ""), 
         fig, width = 9, height = 7, units = "cm")
}

# plot conditional responses from GP-EDM ------------------------------ 
# control parameter, initial conditions, and GP-EDM settings for selected model
if (func_name == "predator_prey_2sp_map") {
  n_sp <- 2
  func <- predator_prey_2sp_map
  x_init <- c(x1 = 1, x2 = 1)
  control_true <- seq(2.4, 3.4, by = 0.005)
  control_gp <- c(2.4, 2.6, 3.15, 3.4)
  control_test <- c(2.5, 2.8, 3, 3.3)
  parms <- list(r = 2.5, k = 4, p = NA, h = 0.1, e = 0.5, m = 2, 
                s_mean = NA, s_sd = NA)
  E <- 5
  sp_input <- 1
  tau <- 1
  s <- 0.02
  control_name <- "Attack rate"
  log_abund <- FALSE
}
# full time series length
t_max_full <- 500
# final time series length
t_max <- 100
# time series length to train GP-EDM
t_max_train <- 50
# data frames to save results
results_no_noise <- data.frame()
results_noise <- data.frame()

# build time series from model ------------------------------ 
for (i in 1:length(control_true)) {
  # parameters without noise
  parms$p <- control_true[i]
  parms$s_mean <- 0
  parms$s_sd <- 0
  # simulate model without noise
  curr_results <- func(t_max = t_max_full, x = x_init, p = unlist(parms), 
                       control_change = control_change, control_noise = control_noise)
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

# fit GP-EDM to training data with certain control parameter values ------------------------------ 
# extract data frames
tol <- 10^(-10)
training_df <- subset(results_noise, (abs(results_noise$control - control_gp[1]) < tol) | 
                        (abs(results_noise$control - control_gp[2]) < tol) |
                        (abs(results_noise$control - control_gp[3]) < tol) |
                        (abs(results_noise$control - control_gp[4]) < tol))
training_df$time <- rep(1:t_max_train, 4)
# plot time series used to train the GP-EDM model
plot_df <- gather(training_df, species, density, x1:x2)
names(plot_df) <- c("time", "p", "species", "density")
plot_df$p <- as.factor(plot_df$p)
fig <- ggplot(data = plot_df, aes(x = time, y = density, color = p)) +
  geom_line(aes(linetype = species), size = 1) +
  scale_color_viridis(name = "Attack\nrate (p)", discrete = TRUE) +
  facet_wrap(~p, labeller = "label_both", ncol = length(control_gp), scales = "free_x") +
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
        legend.key.size = unit(0.4, 'cm'))
if (save_plots) {
  ggsave(paste("figs/fig1/fig1_", func_name, "_ts_different_driver_values.pdf", sep = ""), 
         fig, width = 32, height = 6, units = "cm")
}
# create lagged data frame for GP-EDM
training_list <- split(training_df, training_df$control)
full_training_lags <- data.frame()
for (i in 1:length(training_list)) {
  training_lags <- as.data.frame(makelags(data = training_list[[i]], y = colnames(training_df)[sp_input+1], 
                                          time = "time", E = E, tau = tau))
  training_lags$control_1 <- c(NA, training_list[[i]]$control[2:nrow(training_list[[i]])])
  training_lags_full <- cbind(training_list[[i]], training_lags)
  full_training_lags <- rbind(full_training_lags, training_lags_full)
}
# fit GP-EDM on training data
GP_fit_sp1 <- fitGP(data = full_training_lags, y = "x1", 
                    x = colnames(training_lags), scaling = "local",
                    fixedpars = c(rep(NA, length(colnames(training_lags))-1), NA, NA, NA))
GP_fit_sp2 <- fitGP(data = full_training_lags, y = "x2", 
                    x = colnames(training_lags), scaling = "local",
                    fixedpars = c(rep(NA, length(colnames(training_lags))-1), NA, NA, NA))

# examine conditional responses of fitted model ------------------------------ 
# extract conditional responses
x1_1 <- seq(1, 4, by = 0.05)
x1_2 <- mean(full_training_lags$x1_2, na.rm = TRUE)
x1_3 <- mean(full_training_lags$x1_3, na.rm = TRUE)
x1_4 <- mean(full_training_lags$x1_4, na.rm = TRUE)
x1_5 <- mean(full_training_lags$x1_5, na.rm = TRUE)
# create grid with control parameter and abundances
cond_resp_df <- expand_grid(p, x1_1, x1_2, x1_3, x1_4, x1_5)
names(cond_resp_df) <- c("control_1", paste("x1_", 1:E, sep = ""))
cond_resp_df$x1 <- predict(GP_fit_sp1, newdata = cond_resp_df)$outsampresults$predmean
names(cond_resp_df) <- c("p", paste("x1_", 1:E, sep = ""), "x1")
cond_resp_df$p <- as.factor(cond_resp_df$p)
plot_training_df <- full_training_lags[ , c("x1", paste("x1_", 1:E, sep = ""), "control_1")]
names(plot_training_df) <- c("x1", paste("x1_", 1:E, sep = ""), "p")
plot_training_df$p <- as.factor(plot_training_df$p)
# plot
fig <- ggplot() +
  geom_line(data = cond_resp_df, aes(x = x1_1, y = x1, group = p, color = p), size = 1) + 
  scale_color_viridis(name = "Attack\nrate (p)", discrete = TRUE) +
  xlab(label = TeX("Current prey population ($x_{t}$)")) +
  ylab(label = TeX("Future prey population ($x_{t+1}$)")) +
  scale_y_continuous(breaks = c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.2),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = "white", size = 1.2),
        title = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.7, 'cm'))
if (save_plots) {
  ggsave(paste("figs/fig1/fig1_", func_name, "_GPEDM_cond_resp.pdf", sep = ""), 
         fig, width = 9, height = 7, units = "cm")
}

# extrapolate GP-EDM to unseen control parameter values ------------------------------ 
# extract data frames
tol <- 10^(-10)
training_df_test <- subset(results_no_noise, (abs(results_no_noise$control - control_test[1]) < tol) | 
                             (abs(results_no_noise$control - control_test[2]) < tol) |
                             (abs(results_no_noise$control - control_test[3]) < tol) |
                             (abs(results_no_noise$control - control_test[4]) < tol))
training_df_test$time <- rep(1:t_max, 4)
# to save results
full_test_set <- data.frame()
# loop over control parameter values
for (i in 1:length(control_test)) {
  print(i)
  # choosing starting abundance values
  curr_control <- control_gp[which.min(sqrt((control_test[i] - control_gp)^2))]
  x_init <- tail(subset(training_df, abs(control - curr_control) < tol), E + 1)[ , -c(1, ncol(training_df))]
  # build test set to make predictions
  test_df <- data.frame(time = 1:(E + 1),
                        x1 = as.numeric(x_init[ , 1]),
                        x2 = as.numeric(x_init[ , 2]),
                        control = control_test[i])
  # iterate forecasts for t_max_full time steps
  for (j in 1:t_max_full) {
    # build data frame to make forecast
    test_lags <- makelags(data = test_df, y = colnames(training_df)[sp_input+1],
                          time = "time", E = E, tau = tau, forecast = TRUE)
    test_lags$control_1 <- control_test[i]
    # make forecasts
    forecast_sp1 <- predict(GP_fit_sp1, newdata = test_lags)$outsampresults$predmean
    test_df[j + E + 1, "x1"] <- forecast_sp1
    forecast_sp2 <- predict(GP_fit_sp2, newdata = test_lags)$outsampresults$predmean
    test_df[j + E + 1, "x2"] <- forecast_sp2
    # update time and control level in test set
    test_df[j + E + 1, "time"] <- test_df[j + E, "time"] + 1
    test_df[j + E + 1, "control"] <- control_test[i]
  }
  # extract tail from test set
  test_df <- tail(test_df, t_max)
  # merge with full data frame
  full_test_set <- rbind(full_test_set, test_df)
}
full_test_set$time <- rep(1:t_max, 4)

# plotting observed and predicted dynamics for unobserved driver levels ------------------------------ 
# merge data frames with true and predicted bifurcation diagrams
results_no_noise$type <- "true"
full_test_set$type <- "predicted"
plot_df <- rbind(results_no_noise, full_test_set)
# plot bifurcation diagram
fig <- ggplot(data = plot_df, aes(x = control, y = x1, color = type)) +
  geom_point(size = 1.7) +
  scale_color_manual(values = c("#CB181D", "#919191")) +
  xlab(label = "Attack rate (p)") +
  ylab(label = "Prey population abundance (x)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1),
        axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.position = "none")
if (save_plots) {
  ggsave(paste("figs/fig1/fig1_", func_name, "_bif_diag.pdf", sep = ""), 
         fig, width = 18, height = 12, units = "cm")
}
# plot true and predicted time series
training_df_test$type <- "True"
full_test_set$type <- "Predicted"
plot_df <- rbind(training_df_test, full_test_set)
names(plot_df) <- c("time", "x", "y", "p", "type")
plot_df$p <- as.factor(plot_df$p)
plot_df$type <- factor(plot_df$type, levels = c("Predicted", "True"))
fig <- ggplot(data = plot_df, aes(x = time, y = x, group = type, color = type)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("#CB181D", "#919191")) +
  facet_grid(type~p, labeller = "label_both", scales = "free_x") +
  xlab(label = "Time") +
  ylab(label = "Prey population abundance (x)") +
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
        legend.key.size = unit(0.4, 'cm'))
if (save_plots) {
  ggsave(paste("figs/fig1/fig1_", func_name, "_ts_different_driver_values_predicted.pdf", sep = ""), 
         fig, width = 30, height = 10, units = "cm")
}
