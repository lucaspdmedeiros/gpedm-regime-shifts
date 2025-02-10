# Code for Fig 4. For the experimental microbial ecosystem from Becks et al 2005,
# train a GP-EDM model using data from different dilution rate treatments and then 
# extrapolate to unseen values of dilution rate to reconstruct the bifurcation diagram.

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
if (!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if (!require(viridis)) {install.packages("viridis"); library(viridis)}
if (!require(GPEDM)) {install.packages("GPEDM"); library(GPEDM)}
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(latex2exp)) {install.packages("latex2exp"); library(latex2exp)}
if (!require(scales)) {install.packages("scales"); library(scales)}

# set up ------------------------------ 
# whether to save plots
save_plots <- TRUE
# whether to use native or delay coordinates
coordinates <- "delay"
# whether to log species abundances
log_abund <- TRUE
# number of species
n_sp <- 3
# grid of embedding dimension values and time lag value
E_values <- 1:9
tau <- 1
# defining coordinates to use
if (coordinates == "native") {
  E <- 1
  sp_input <- c(1, 2, 3)
  sp_plot <- 1
}
if (coordinates == "delay") {
  sp_input <- 1
  sp_plot <- sp_input
}
# whether to scale data within each treatment
treatment_scaling <- TRUE
# whether to remove transient dynamics
remove_transient <- TRUE
transient_length <- 5
# values of dilution rate to extrapolate GP-EDM
dilution_grid <- seq(from = 0.45, to = 0.9, by = 0.001)
# full time series length
t_max_full <- 1000
# final time series length
t_max <- 100
# fixed value for the inverse length scale of dilution rate 
# (set to NA for phi to be estimated from data)
fixed_phi_dilution <- 1

# load and format data ------------------------------ 
# load data 
data <- read.csv("data/microbial_experiment/becks_etal_2005.csv")
names(data) <- c("time", "predator", "prey1", "prey2", "dilution", "replicate")
# data to use for analyses
training_df <- data
# change abundances to log scale
if (log_abund) {
  training_df[ , c("predator", "prey1", "prey2")] <- 
    log(training_df[ , c("predator", "prey1", "prey2")])
}
# extract dilution/replicate combinations
dilutions <- c(0.45, 0.45, 0.5, 0.5, 0.5, 0.5, 0.5, 0.75, 0.9)
replicates <- c(1, 2, 1, 2, 3, 4, 5, 1, 1)
combinations <- data.frame(dilutions = dilutions, replicates = replicates)
combination_name <- apply(combinations, 1, paste, collapse = "_")
# add variable combining dilution rate and replicate
training_df$combination <- apply(training_df[ , c("dilution", "replicate")], 1, paste, collapse = "_")
# save indexes of transient dynamics
transient_indexes <- c()
for (i in 1:nrow(combinations)) {
  dilution_curr <- combinations$dilutions[i]
  replicate_curr <- combinations$replicates[i]
  if (!((dilution_curr == 0.5 & replicate_curr == 4) | 
        (dilution_curr == 0.5 & replicate_curr == 5))) {
    transient_indexes_curr <- which(training_df$dilution == dilution_curr & 
                                      training_df$replicate == replicate_curr)[1:transient_length]
    transient_indexes <- c(transient_indexes, transient_indexes_curr)
  }
}
# remove transient dynamics
if (remove_transient) {
  training_df <- training_df[-transient_indexes, ]
}
# scale to zero mean within dilution treatments
if (treatment_scaling) {
  mean_sp1 <- c()
  mean_sp2 <- c()
  mean_sp3 <- c()
  for (i in 1:length(unique(dilutions))) {
    dilution_curr <- unique(dilutions)[i]
    training_df_curr <- subset(training_df, dilution == dilution_curr)
    curr_mean_sp1 <- rep(mean(training_df_curr$predator, na.rm = TRUE), nrow(training_df_curr))
    curr_mean_sp2 <- rep(mean(training_df_curr$prey1, na.rm = TRUE), nrow(training_df_curr))
    curr_mean_sp3 <- rep(mean(training_df_curr$prey2, na.rm = TRUE), nrow(training_df_curr))
    training_df_curr$predator <- (training_df_curr$predator - curr_mean_sp1)
    training_df_curr$prey1 <- (training_df_curr$prey1 - curr_mean_sp2)
    training_df_curr$prey2 <- (training_df_curr$prey2 - curr_mean_sp3)
    training_df[training_df$dilution == dilution_curr, ] <- training_df_curr
    mean_sp1 <- c(mean_sp1, curr_mean_sp1)
    mean_sp2 <- c(mean_sp2, curr_mean_sp2)
    mean_sp3 <- c(mean_sp3, curr_mean_sp3)
  }
}

# plot time series ------------------------------ 
# scale time series for plot
plot_df <- training_df
plot_df[, c("predator", "prey1", "prey2")] <- scale(plot_df[, c("predator", "prey1", "prey2")])
# create long format data frame
plot_df <- gather(plot_df, "species", "abundance", "predator":"prey2")
# plot time series for each dilution/replicate combination
fig <- ggplot(data = plot_df, aes(x = time, y = abundance, 
                                  group = species, color = species)) +
  geom_line(size = 1) +
  facet_wrap(dilution~replicate, labeller = "label_both", nrow = 5) +
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
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.size = unit(0.6, 'cm'))

# fit GP-EDM to data to determine best embedding dimension ------------------------------ 
# create list of each subset of data set
training_list <- list()
for (i in 1:length(unique(training_df$combination))) {
  comb_curr <- unique(training_df$combination)[i]
  training_list[[i]] <- subset(training_df, combination == comb_curr)
}
# create data frame with all lags
full_training_lags <- data.frame()
for (i in 1:length(training_list)) {
  training_lags <- as.data.frame(makelags(data = training_list[[i]], y = colnames(training_df)[sp_input+1], 
                                          time = "time", E = max(E_values), tau = tau))
  training_lags$dilution_1 <- c(NA, training_list[[i]]$dilution[2:nrow(training_list[[i]])])
  training_lags_full <- cbind(training_list[[i]], training_lags)
  full_training_lags <- rbind(full_training_lags, training_lags_full)
}
# for each value of E, create lagged data frame for GP-EDM and evaluate prediction accuracy
if (coordinates == "delay") {
  R2 <- c()
  phi_dilution <- c()
  for (i in 1:length(E_values)) {
    print(i)
    # input names
    inputs <- paste(names(training_df)[sp_input + 1], E_values[1:i], sep = "_")
    inputs <- c(inputs, "dilution_1")
    # fit GP-EDM to training data and compute prediction accuracy
    if (sp_input == 1) {
      GP_fit_sp1 <- fitGP(data = full_training_lags, y = "predator", 
                          x = inputs, scaling = "local",
                          fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_dilution, NA, NA),
                          predictmethod = "loo")
      if (treatment_scaling) {
        full_training_lags$obs <- GP_fit_sp1$outsampresults$obs + mean_sp1
        full_training_lags$predmean <- GP_fit_sp1$outsampresults$predmean + mean_sp1
      } else {
        full_training_lags$obs <- GP_fit_sp1$outsampresults$obs
        full_training_lags$predmean <- GP_fit_sp1$outsampresults$predmean
      }
      R2_df <- full_training_lags
      phi_dilution[i] <- GP_fit_sp1$pars[E_values[i]+1]
    }
    if (sp_input == 2) {
      GP_fit_sp2 <- fitGP(data = full_training_lags, y = "prey1", 
                          x = inputs, scaling = "local",
                          fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_dilution, NA, NA),
                          predictmethod = "loo")
      if (treatment_scaling) {
        full_training_lags$obs <- GP_fit_sp2$outsampresults$obs + mean_sp2
        full_training_lags$predmean <- GP_fit_sp2$outsampresults$predmean + mean_sp2
      } else {
        full_training_lags$obs <- GP_fit_sp2$outsampresults$obs
        full_training_lags$predmean <- GP_fit_sp2$outsampresults$predmean
      }
      R2_df <- full_training_lags
      phi_dilution[i] <- GP_fit_sp2$pars[E_values[i]+1]
    }
    if (sp_input == 3) {
      GP_fit_sp3 <- fitGP(data = full_training_lags, y = "prey2", 
                          x = inputs, scaling = "local",
                          fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_dilution, NA, NA),
                          predictmethod = "loo")
      if (treatment_scaling) {
        full_training_lags$obs <- GP_fit_sp3$outsampresults$obs + mean_sp3
        full_training_lags$predmean <- GP_fit_sp3$outsampresults$predmean + mean_sp3
      } else {
        full_training_lags$obs <- GP_fit_sp3$outsampresults$obs
        full_training_lags$predmean <- GP_fit_sp3$outsampresults$predmean
      }
      R2_df <- full_training_lags
      phi_dilution[i] <- GP_fit_sp3$pars[E_values[i]+1]
    }
    # compute R2
    R2_df <- na.omit(R2_df)
    if (log_abund) {
      R2_df$obs <- exp(R2_df$obs)
      R2_df$predmean <- exp(R2_df$predmean)
    }
    R2[i] <- 1 - (sum((R2_df$obs - R2_df$predmean)^2) / sum((R2_df$obs - mean(R2_df$obs))^2))
  }
  # select best embedding dimension
  if (is.na(fixed_phi_dilution)) {
    E <- order(R2, decreasing = TRUE)[which(phi_dilution[order(R2, decreasing = TRUE)] > 0.5)[1]]
  } else {
    E <- which.max(R2)
  }
}

# fit GP-EDM to data using best embedding dimension ------------------------------ 
# fit GP-EDM in delay coordinates
if (coordinates == "delay") {
  # input names
  inputs <- paste(names(training_df)[sp_input + 1], 1:E, sep = "_")
  inputs <- c(inputs, "dilution_1")
  # fit GP-EDM
  if (sp_input == 1) {
    GP_fit_sp1 <- fitGP(data = full_training_lags, y = "predator", 
                        x = inputs, scaling = "local",
                        fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_dilution, NA, NA),
                        predictmethod = "loo")
    if (treatment_scaling) {
      full_training_lags$obs <- GP_fit_sp1$outsampresults$obs + mean_sp1
      full_training_lags$predmean <- GP_fit_sp1$outsampresults$predmean + mean_sp1
    } else {
      full_training_lags$obs <- GP_fit_sp1$outsampresults$obs
      full_training_lags$predmean <- GP_fit_sp1$outsampresults$predmean
    }
  }
  if (sp_input == 2) {
    GP_fit_sp2 <- fitGP(data = full_training_lags, y = "prey1", 
                        x = inputs, scaling = "local",
                        fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_dilution, NA, NA),
                        predictmethod = "loo")
    if (treatment_scaling) {
      full_training_lags$obs <- GP_fit_sp2$outsampresults$obs + mean_sp2
      full_training_lags$predmean <- GP_fit_sp2$outsampresults$predmean + mean_sp2
    } else {
      full_training_lags$obs <- GP_fit_sp2$outsampresults$obs
      full_training_lags$predmean <- GP_fit_sp2$outsampresults$predmean
    }
  }
  if (sp_input == 3) {
    GP_fit_sp3 <- fitGP(data = full_training_lags, y = "prey2", 
                        x = inputs, scaling = "local",
                        fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_dilution, NA, NA),
                        predictmethod = "loo")
    if (treatment_scaling) {
      full_training_lags$obs <- GP_fit_sp3$outsampresults$obs + mean_sp3
      full_training_lags$predmean <- GP_fit_sp3$outsampresults$predmean + mean_sp3
    } else {
      full_training_lags$obs <- GP_fit_sp3$outsampresults$obs
      full_training_lags$predmean <- GP_fit_sp3$outsampresults$predmean
    }
  }
  plot_df <- full_training_lags
  # compute prediction errors
  plot_df$residuals <- plot_df$obs - plot_df$predmean
}
# fit GP-EDM in native coordinates
if (coordinates == "native") {
  # input names
  inputs <- paste(rep(names(training_df)[sp_input + 1], each = E), 1:E, sep = "_")
  inputs <- c(inputs, "dilution_1")
  # fit GP-EDM
  GP_fit_sp1 <- fitGP(data = full_training_lags, y = "predator", 
                      x = inputs, scaling = "local",
                      fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_dilution, NA, NA),
                      predictmethod = "loo")
  GP_fit_sp2 <- fitGP(data = full_training_lags, y = "prey1", 
                      x = inputs, scaling = "local",
                      fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_dilution, NA, NA),
                      predictmethod = "loo")
  GP_fit_sp3 <- fitGP(data = full_training_lags, y = "prey2", 
                      x = inputs, scaling = "local",
                      fixedpars = c(rep(NA, length(inputs)-1), fixed_phi_dilution, NA, NA),
                      predictmethod = "loo")
  # compute R2
  if (treatment_scaling) {
    full_training_lags$obs_sp1 <- GP_fit_sp1$outsampresults$obs + mean_sp1
    full_training_lags$predmean_sp1 <- GP_fit_sp1$outsampresults$predmean + mean_sp1
    full_training_lags$obs_sp2 <- GP_fit_sp2$outsampresults$obs + mean_sp2
    full_training_lags$predmean_sp2 <- GP_fit_sp2$outsampresults$predmean + mean_sp2
    full_training_lags$obs_sp3 <- GP_fit_sp3$outsampresults$obs + mean_sp3
    full_training_lags$predmean_sp3 <- GP_fit_sp3$outsampresults$predmean + mean_sp3
  } else {
    full_training_lags$obs_sp1 <- GP_fit_sp1$outsampresults$obs
    full_training_lags$predmean_sp1 <- GP_fit_sp1$outsampresults$predmean
    full_training_lags$obs_sp2 <- GP_fit_sp2$outsampresults$obs
    full_training_lags$predmean_sp2 <- GP_fit_sp2$outsampresults$predmean
    full_training_lags$obs_sp3 <- GP_fit_sp3$outsampresults$obs
    full_training_lags$predmean_sp3 <- GP_fit_sp3$outsampresults$predmean
  }
  R2_df <- full_training_lags
  R2_df <- na.omit(R2_df)
  if (log_abund) {
    R2_df$obs_sp1 <- exp(R2_df$obs_sp1)
    R2_df$predmean_sp1 <- exp(R2_df$predmean_sp1)
    R2_df$obs_sp2 <- exp(R2_df$obs_sp2)
    R2_df$predmean_sp2 <- exp(R2_df$predmean_sp2)
    R2_df$obs_sp3 <- exp(R2_df$obs_sp3)
    R2_df$predmean_sp3 <- exp(R2_df$predmean_sp3)
  }
  R2_sp1 <- 1 - (sum((R2_df$obs_sp1 - R2_df$predmean_sp1)^2) / sum((R2_df$obs_sp1 - mean(R2_df$obs_sp1))^2))
  R2_sp2 <- 1 - (sum((R2_df$obs_sp2 - R2_df$predmean_sp2)^2) / sum((R2_df$obs_sp2 - mean(R2_df$obs_sp2))^2))
  R2_sp3 <- 1 - (sum((R2_df$obs_sp3 - R2_df$predmean_sp3)^2) / sum((R2_df$obs_sp3 - mean(R2_df$obs_sp3))^2))
  # data frame for plotting
  if (sp_plot == 1) {
    plot_df <- full_training_lags
    plot_df$obs <- full_training_lags$obs_sp1
    plot_df$predmean <- full_training_lags$predmean_sp1
  }
  if (sp_plot == 2) {
    plot_df <- full_training_lags
    plot_df$obs <- full_training_lags$obs_sp2
    plot_df$predmean <- full_training_lags$predmean_sp2
  }
  if (sp_plot == 3) {
    plot_df <- full_training_lags
    plot_df$obs <- full_training_lags$obs_sp3
    plot_df$predmean <- full_training_lags$predmean_sp3
  }
}
# plot observations and predictions
if (log_abund) {
  plot_df$obs <- exp(plot_df$obs)
  plot_df$predmean <- exp(plot_df$predmean)
}
# compute prediction errors
plot_df$residuals <- plot_df$obs - plot_df$predmean
if (sp_input == 1) {
  sp_name <- "Predator abundance (ind/ml)"
}
if (sp_input == 2) {
  sp_name <- "Preferred prey abundance (cells/ml)"
}
if (sp_input == 3) {
  sp_name <- "Less-preferred prey abundance (cells/ml)"
}
# plot observations and predictions
fig <- ggplot(data = plot_df) +
  geom_line(aes(x = time, y = obs), size = 1, color = "black") +
  geom_line(aes(x = time, y = predmean), size = 1, color = "#CB181D") +
  xlab(label = "Time") +
  ylab(label = sp_name) +
  facet_wrap(dilution~replicate,
             labeller = "label_both", nrow = 3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.2),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white", size = 1.2),
        title = element_text(size = 14),
        axis.title = element_text(size = 17),
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

# make predictions for a range of dilution rate values ------------------------------ 
# observed dilution rate values
dilution_obs <- unique(data$dilution)
full_test_set <- data.frame()
# loop over all dilution rate values
for (i in 1:length(dilution_grid)) {
  print(i)
  # choosing starting abundance values
  tol <- 10^(-10)
  curr_dilution <- dilution_obs[which.min(sqrt((dilution_grid[i] - dilution_obs)^2))]
  abund_init <- tail(subset(training_df, abs(dilution - curr_dilution) < tol), E + 1)[ , -c(1, ncol(training_df))]
  # build test set to make predictions
  test_df <- data.frame(time = 1:(E + 1),
                        predator = as.numeric(abund_init[ , 1]),
                        prey1 = as.numeric(abund_init[ , 2]),
                        prey2 = as.numeric(abund_init[ , 3]),
                        dilution = dilution_grid[i])
  # iterate forecasts for t_max_full time steps
  for (j in 1:t_max_full) {
    # build data frame to make forecast
    test_lags <- makelags(data = test_df, y = colnames(training_df)[sp_input+1],
                          time = "time", E = E, tau = tau, forecast = TRUE)
    test_lags$dilution_1 <- dilution_grid[i]
    # make forecast
    if (coordinates == "delay") {
      if (sp_input == 1) {
        forecast_sp1 <- predict(GP_fit_sp1, newdata = test_lags)$outsampresults$predmean
        test_df[j + E + 1, "predator"] <- forecast_sp1
      }
      if (sp_input == 2) {
        forecast_sp2 <- predict(GP_fit_sp2, newdata = test_lags)$outsampresults$predmean
        test_df[j + E + 1, "prey1"] <- forecast_sp2
      }
      if (sp_input == 3) {
        forecast_sp3 <- predict(GP_fit_sp3, newdata = test_lags)$outsampresults$predmean
        test_df[j + E + 1, "prey2"] <- forecast_sp3
      }
    }
    if (coordinates == "native") {
      forecast_sp1 <- predict(GP_fit_sp1, newdata = test_lags)$outsampresults$predmean
      test_df[j + E + 1, "predator"] <- forecast_sp1
      forecast_sp2 <- predict(GP_fit_sp2, newdata = test_lags)$outsampresults$predmean
      test_df[j + E + 1, "prey1"] <- forecast_sp2
      forecast_sp3 <- predict(GP_fit_sp3, newdata = test_lags)$outsampresults$predmean
      test_df[j + E + 1, "prey2"] <- forecast_sp3
    }
    # update time and control level in test set
    test_df[j + E + 1, "time"] <- test_df[j + E, "time"] + 1
    test_df[j + E + 1, "dilution"] <- dilution_grid[i]
  }
  # extract tail from test set
  test_df <- tail(test_df, t_max)
  # merge with full data frame
  full_test_set <- rbind(full_test_set, test_df)
}

# plot predicted bifurcation diagrams ------------------------------ 
# data frame to plot
plot_df <- full_test_set
# backtransform logged abundances
if (log_abund) {
  plot_df[ , -c(1, ncol(plot_df))] <- 
    exp(plot_df[ , -c(1, ncol(plot_df))])
  training_df[ , 2:(n_sp+1)] <- 
    exp(training_df[ , 2:(n_sp+1)])
}
# plot settings
if (coordinates == "delay") {
  R2_chosen <- R2[E]
  if (sp_input == 1) {
    sp_col_data <- "#DEDEDE"
    sp_col_pred <- "#FB6A4A"
    phi_values <- GP_fit_sp1$pars[1:(E+1)]
    names(phi_values) <- GP_fit_sp1$inputs$x_names
    max_abund <- max(plot_df$predator)
    min_abund <- min(plot_df$predator)
    plot_df$abundance <- plot_df$predator
    training_df$abundance <- training_df$predator
  }
  if (sp_input == 2) {
    sp_col_data <- "#919191"
    sp_col_pred <- "#CB181D"
    phi_values <- GP_fit_sp2$pars[1:(E+1)]
    names(phi_values) <- GP_fit_sp2$inputs$x_names
    max_abund <- max(plot_df$prey1)
    min_abund <- min(plot_df$prey1)
    plot_df$abundance <- plot_df$prey1
    training_df$abundance <- training_df$prey1
  }
  if (sp_input == 3) {
    sp_col_data <- "#454545"
    sp_col_pred <- "#67000D"
    phi_values <- GP_fit_sp3$pars[1:(E+1)]
    names(phi_values) <- GP_fit_sp3$inputs$x_names
    max_abund <- max(plot_df$prey2)
    min_abund <- min(plot_df$prey2)
    plot_df$abundance <- plot_df$prey2
    training_df$abundance <- training_df$prey2
  }
}
if (coordinates == "native") {
  if (sp_plot == 1) {
    R2_chosen <- R2_sp1
    phi_values <- GP_fit_sp1$pars[1:length(GP_fit_sp1$inputs$x_names)]
    names(phi_values) <- GP_fit_sp1$inputs$x_names
    max_abund <- max(plot_df$predator)
    min_abund <- min(plot_df$predator)
    plot_df$abundance <- plot_df$predator
    training_df$abundance <- training_df$predator
  }
  if (sp_plot == 2) {
    R2_chosen <- R2_sp2
    phi_values <- GP_fit_sp2$pars[1:length(GP_fit_sp2$inputs$x_names)]
    names(phi_values) <- GP_fit_sp2$inputs$x_names
    max_abund <- max(plot_df$prey1)
    min_abund <- min(plot_df$prey1)
    plot_df$abundance <- plot_df$prey1
    training_df$abundance <- training_df$prey2
  }
  if (sp_plot == 3) {
    R2_chosen <- R2_sp3
    phi_values <- GP_fit_sp3$pars[1:length(GP_fit_sp3$inputs$x_names)]
    names(phi_values) <- GP_fit_sp3$inputs$x_names
    max_abund <- max(plot_df$prey2)
    min_abund <- min(plot_df$prey2)
    plot_df$abundance <- plot_df$prey2
    training_df$abundance <- training_df$prey2
  }
}
# plot predicted bifurcation diagram
fig <- ggplot() +
  geom_vline(xintercept = dilution_obs, linetype = "dashed", size = 0.7) +
  geom_point(data = plot_df, aes(x = dilution, y = abundance),
             size = 0.7, color = sp_col_pred) +
  geom_point(data = training_df, aes(x = dilution, y = abundance),
             shape = 21, size = 3.5, fill = sp_col_data) +
  scale_x_continuous(limits = c(min(dilution_grid), max(dilution_grid))) +
  scale_y_continuous(labels = label_number(accuracy = 0.1)) +
  xlab(label = "Dilution rate") +
  ylab(label = "Population abundance") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1),
        plot.title = element_text(size = 22),
        plot.subtitle = element_text(size = 18),
        axis.title = element_text(size = 22),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.position = "none",
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.key.size = unit(0.6, 'cm'))

# save plot ------------------------------ 
if (save_plots) {
  # to name files
  if (treatment_scaling) {
    treatment_scaling_name <- "yes_treat_scaling"
  } else {
    treatment_scaling_name <- "no_treat_scaling"
  }
  if (log_abund) {
    log_abund_name <- "yes_log_abund"
  } else {
    log_abund_name <- "no_log_abund"
  }
  if (remove_transient) {
    remove_transient_name <- "no_transient"
  } else {
    remove_transient_name <- "yes_transient"
  }
  # save plot
  ggsave(paste("figs/fig4/fig4_becks_etal_2005_predicted_bif_diag_", coordinates, 
               "_sp_input_", sp_input, "_dilution_phi_", fixed_phi_dilution,
               "_E_", E, "_", log_abund_name, "_", treatment_scaling_name,
               "_", remove_transient_name, ".pdf", sep = ""), 
         fig, width = 26, height = 12, units = "cm")
}
