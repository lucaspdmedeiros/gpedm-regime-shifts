# Code for Fig 5. For the Lake Zurich plankton ecosystem from Merz et al 2023,
# train a GP-EDM model using successively larger training data sets from large green 
# algae and gradually decreasing phosphate concentration. Then, extrapolate the model
# to unseen values of phosphate concentration to reconstruct the bifurcation diagram.

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/ts_breakpoint.R")
source("code/functions/power_spectrum.R")
source("code/functions/le_qr_decomp.R")
if (!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if (!require(gridExtra)) {install.packages("gridExtra"); library(gridExtra)}
if (!require(viridis)) {install.packages("viridis"); library(viridis)}
if (!require(GPEDM)) {install.packages("GPEDM"); library(GPEDM)}
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(MARSS)) {install.packages("MARSS"); library(MARSS)}
if (!require(zoo)) {install.packages("zoo"); library(zoo)}
if (!require(plyr)) {install.packages("plyr"); library(plyr)}
if (!require(tidyverse)) {install.packages("tidyverse"); library(tidyverse)}
if (!require(doBy)) {install.packages("doBy"); library(doBy)}
if (!require(stats)) {install.packages("stats"); library(stats)}
if (!require(zoo)) {install.packages("zoo"); library(zoo)}

# set up and load data ------------------------------ 
# whether to save plots
save_plot <- TRUE
# whether to filter phosphate time series
filter_phosphate <- TRUE
# whether to interpolate missing data points
interpolate <- FALSE
# whether to log transform abundances
log_abund <- FALSE
# maximum lag for species abundance
max_E_sp <- 9
# maximum lag for phosphate
max_E_phosphate <- 1
# maximum lag for temperature (if NA, does not include temperature)
max_E_temperature <- 2
# fixed value for the inverse length scale of the control parameter
# (set to NA for phi to be estimated from data)
fixed_phi_control <- NA
# time lag
tau <- 1
# full time series length
t_max_full <- 1000
# final time series length
t_max <- 100
# temporal resolution
temporal_resolution <- "month"
# lake to use
lake <- "zurich"
# species, phosphate levels, and end points to use
sp_input <- "Green algae (large)"
phosphate_grid <- seq(0.02, 0.07, length = 301)
end_points <- c(106, 118, 166, 226)
# phosphate levels to use for conditional responses, power spectrum, and Lyapunov exponent
phosphate_levels <- seq(0.02, 0.07, by = 0.01)
# load data
data <- read.csv(paste("data/lake_plankton/lake_", lake,
                       ".csv", sep = ""))
# change species names
names(data)[names(data) == "Ci"] <- "Ciliates"
names(data)[names(data) == "Cr1"] <- "Cryptophytes (small)"
names(data)[names(data) == "Cr2"] <- "Cryptophytes (large)"
names(data)[names(data) == "Cy"] <- "Cyanobacteria"
names(data)[names(data) == "Di1"] <- "Diatoms (small)"
names(data)[names(data) == "Di2"] <- "Diatoms (large)"
names(data)[names(data) == "Go1"] <- "Gold algae (small)"
names(data)[names(data) == "Go2"] <- "Gold algae (large)"
names(data)[names(data) == "Gr1"] <- "Green algae (small)"
names(data)[names(data) == "Gr2"] <- "Green algae (large)"
names(data)[names(data) == "He"] <- "Large herbivores"
names(data)[names(data) == "Mi"] <- "Mixotrophic flagellates"
names(data)[names(data) == "Om"] <- "Omnivores"
names(data)[names(data) == "Pr"] <- "Invertebrate predators"
names(data)[names(data) == "Ro"] <- "Herbivore rotifers"
names(data)[names(data) == "temperature"] <- "Temperature"
names(data)[names(data) == "phosphorus"] <- "Phosphate"
month <- month.abb[data$month]
data$month_year <- as.yearmon(paste(month, data$year))
data$quarter_year <- as.yearqtr(data$month_year)
# monthly or quarterly averages
if (temporal_resolution == "month") {
  summ_df <- data %>% group_by(month_year) %>% summarise_all(mean, na.rm = TRUE)
  summ_df$Time <- summ_df$month_year
  summ_df$Year <- summ_df$year
  summ_df$Month <- summ_df$month
}
if (temporal_resolution == "quarter") {
  summ_df <- data %>% group_by(quarter_year) %>% summarise_all(mean, na.rm = TRUE)
  summ_df$Time <- summ_df$quarter_year
  summ_df$Year <- summ_df$year
  summ_df$Month <- summ_df$month
}

# reconstruct bifurcation diagrams for each species ------------------------------ 
# current data set
training_df <- summ_df
# interpolate missing data
if (interpolate) {
  # specify time points for spline to go through (same number of points as original time series)
  xout <- seq(1, nrow(training_df), by = 1)
  # fit cubic hermite interpolation for each species
  for (j in 1:length(sp_input)) {
    # fit spline
    sp_spline <- spline(x = xout, y = unlist(training_df[ , sp_input[j]]), method = "fmm", xout = xout)
    sp_spline$y[sp_spline$y < 0] <- 0
    # add interpolated data to data frame
    training_df[ , sp_input[j]] <- sp_spline$y
  }
}
# filter phosphate time series
window_size <- 25
mean_phosphate <- c()
for (j in 1:(nrow(training_df)-window_size+1)) {
  mean_phosphate[j] <- mean(training_df$Phosphate[j:(j+window_size-1)], na.rm = TRUE)
}
length_diff <- nrow(training_df) - length(mean_phosphate)
filtered_phosphate <- c(rep(NA, length_diff/2), mean_phosphate, rep(NA, length_diff/2))
points_to_keep <- !is.na(filtered_phosphate)
if (filter_phosphate) {
  training_df$Phosphate <- filtered_phosphate
}
training_df <- training_df[points_to_keep, ]
if (window_size == 13) {
  training_df <- training_df[7:nrow(training_df), ]
}
# current data frame
curr_training_df <- training_df[ , c("Time", "Month", "Phosphate", "Temperature", sp_input)]
# compute rolling variance
window_size <- 48
var_ts <- c()
window_end <- c()
for (j in 1:(nrow(curr_training_df)-window_size)) {
  ts_curr <- curr_training_df[j:(j+window_size), ]
  window_end[j] <- as.character(tail(ts_curr$Time, 1))
  var_ts[j] <- var(ts_curr[ , sp_input], na.rm = TRUE)
}
# plot abundance, variance, and phosphate time series ------------------------------ 
print(curr_training_df[end_points[1], c("Time", "Phosphate")])
print(curr_training_df[end_points[2], c("Time", "Phosphate")])
print(curr_training_df[end_points[3], c("Time", "Phosphate")])
print(curr_training_df[end_points[4], c("Time", "Phosphate")])
subset_training_df <- curr_training_df[ , c("Time", "Phosphate", sp_input)]
subset_training_df$Variance <- c(rep(NA, window_size), var_ts)
subset_training_df$Variance[subset_training_df$Time > as.yearmon(curr_training_df[end_points[4], "Time"])] <- NA
# plot with abundance and phosphate
fig <- ggplot(data = subset_training_df) +
  geom_line(aes(x = Time, y = `Green algae (large)`), color = "#919191", size = 0.9) +
  geom_line(aes(x = Time, y = Phosphate * 1000000), color = "#35B779", size = 0.9) +
  geom_vline(xintercept = curr_training_df$Time[end_points], linetype = "dashed", size = 0.8) +
  geom_vline(xintercept = curr_training_df$Time[end_points[2]+1], size = 0.8) +
  xlab(label = "Month") +
  scale_x_yearmon(n = 10) + 
  scale_y_continuous(name = "Large green algae\nabundance", 
                     sec.axis = sec_axis(~./1000000, name = "Phosphate concentration\n(smoothed)")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1),
        axis.title = element_text(size = 17),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y.right = element_text(color = "#35B779"),
        axis.text.y.right = element_text(color = "#35B779"),
        axis.line.y.right = element_line(color = "#35B779"), 
        axis.ticks.y.right = element_line(color = "#35B779"))
if (save_plot) {
  ggsave(paste("figs/fig5/fig5_", lake, "_phosphate_ts.pdf", sep = ""),
         fig, width = 32, height = 8, units = "cm")
}
# plot with abundance and variance
fig <- ggplot(data = subset_training_df) +
  geom_line(aes(x = Time, y = `Green algae (large)`), color = "#919191", size = 0.9) +
  geom_line(aes(x = Time, y = Variance / 10000), color = "#2171B5", size = 0.9) +
  geom_vline(xintercept = curr_training_df$Time[end_points], linetype = "dashed", size = 0.8) +
  geom_vline(xintercept = curr_training_df$Time[end_points[2]+1], size = 0.8) +
  xlab(label = "Month") +
  scale_x_yearmon(n = 10) + 
  scale_y_continuous(name = "Large green algae\nabundance", 
                     sec.axis = sec_axis(~.*10000, name = "Rolling variance")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1),
        axis.title = element_text(size = 17),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y.right = element_text(color = "#2171B5"),
        axis.text.y.right = element_text(color = "#2171B5"),
        axis.line.y.right = element_line(color = "#2171B5"), 
        axis.ticks.y.right = element_line(color = "#2171B5"))
if (save_plot) {
  ggsave(paste("figs/fig5/fig5_", lake, "_variance_ts_window_", window_size,
               ".pdf", sep = ""),
         fig, width = 26, height = 8, units = "cm")
}
# determine tipping point through change-point analysis for each species
sp_breakpoints <- c()
ss_breakpoints <- c()
for (j in 1:length(sp_input)) {
  breakpoint_results <- ts_breakpoint(ts = training_df[ , c("Time", sp_input[j])], 
                                      n_breaks = 300,
                                      n_segments = 1)
  sp_breakpoints[j] <- breakpoint_results[[1]]
  ss_breakpoints[j] <- breakpoint_results[[2]]
}
breakpoint_df <- data.frame(variable = sp_input,
                            breakpoint = sp_breakpoints,
                            breakpoint_date = curr_training_df$Time[sp_breakpoints],
                            breakpoint_ss = ss_breakpoints,
                            phosphate = curr_training_df$Phosphate[sp_breakpoints])
# create time series of monthly temperature averages
temp_month_avg <- tapply(curr_training_df$Temperature, curr_training_df$Month, mean)
# log transform abundances
if (log_abund) {
  for (j in 1:length(sp_input)) {
    curr_training_df[ , sp_input[j]][curr_training_df[ , sp_input[j]] == 0] <- 1
  }
  curr_training_df[ , sp_input] <- 
    log(curr_training_df[ , sp_input])
}
# loop over species
for (j in 1:length(sp_input)) {
  
  # determine best embedding dimension ------------------------------ 
  sp_input_curr <- sp_input[j]
  # select data for the current species
  curr_training_df_sp <- curr_training_df[1:end_points[1], ]
  # create list of inputs to try
  input_list <- list()
  list_length <- 1
  if (!is.na(max_E_temperature)) {
    for (k in 1:max_E_sp) {
      sp_lags <- paste(sp_input_curr, 1:k, sep = "_")
      for (l in 1:max_E_phosphate) {
        phosphate_lags <- paste("Phosphate", 1:l, sep = "_")
        for (m in 1:max_E_temperature) {
          temperature_lags <- paste("Temperature", 1:m, sep = "_")
          all_inputs <- c(sp_lags, phosphate_lags, temperature_lags)
          if (length(all_inputs) <= max_E_sp) {
            input_list[[list_length]] <- all_inputs
            list_length <- list_length + 1
          }
        }
      }
    }
  } else {
    for (k in 1:max_E_sp) {
      sp_lags <- paste(sp_input_curr, 1:k, sep = "_")
      for (l in 1:max_E_phosphate) {
        phosphate_lags <- paste("Phosphate", 1:l, sep = "_")
        all_inputs <- c(sp_lags, phosphate_lags)
        if (length(all_inputs) <= max_E_sp) {
          input_list[[list_length]] <- all_inputs
          list_length <- list_length + 1
        }
      }
    }
  }
  # compute R2 for all sets of inputs 
  R2 <- c()
  for (k in 1:length(input_list)) {
    # input names
    inputs <- input_list[[k]]
    # inverse length scale values
    E_sp <- sum(grepl(pattern = sp_input_curr, x = inputs, fixed = TRUE))
    E_phosphate <- sum(grepl(pattern = "Phosphate", x = inputs, fixed = TRUE))
    E_temperature <- sum(grepl(pattern = "Temperature", x = inputs, fixed = TRUE))
    fixed_pars <- c(rep(NA, E_sp), rep(fixed_phi_control, E_phosphate), 
                    rep(NA, E_temperature), NA, NA)
    # create lags
    training_lags <- as.data.frame(makelags(data = curr_training_df_sp, y = colnames(curr_training_df_sp)[-1], 
                                            time = "Time", E = max_E_sp, tau = tau))
    full_training_lags <- cbind(curr_training_df_sp, training_lags)
    # fit GP-EDM on training data
    GP_fit_sp <- fitGP(data = full_training_lags, y = sp_input_curr, 
                       x = inputs, scaling = "local",
                       fixedpars = fixed_pars,
                       predictmethod = "loo")
    # compute R2
    full_training_lags$obs <- GP_fit_sp$outsampresults$obs
    full_training_lags$predmean <- GP_fit_sp$outsampresults$predmean
    R2_df <- full_training_lags
    R2_df <- na.omit(R2_df)
    if (log_abund) {
      R2_df$obs <- exp(R2_df$obs)
      R2_df$predmean <- exp(R2_df$predmean)
    }
    R2[k] <- 1 - (sum((R2_df$obs - R2_df$predmean)^2) / sum((R2_df$obs - mean(R2_df$obs))^2))
  }
  # select best embedding dimension
  best_inputs <- input_list[[which.max(R2)]]
  print(best_inputs)
  print(max(R2))
  
  # perform forecasts to reconstruct bifurcation diagram ------------------------------ 
  # best input names
  inputs <- best_inputs
  # inverse length scale values
  E_sp <- sum(grepl(pattern = sp_input_curr, x = inputs, fixed = TRUE))
  E_phosphate <- sum(grepl(pattern = "Phosphate", x = inputs, fixed = TRUE))
  E_temperature <- sum(grepl(pattern = "Temperature", x = inputs, fixed = TRUE))
  fixed_pars <- c(rep(NA, E_sp), rep(fixed_phi_control, E_phosphate), 
                  rep(NA, E_temperature), NA, NA)
  # starting abundance values
  max_E_curr <- max(c(E_sp, E_phosphate, E_temperature))
  x_init <- as.numeric(unlist(c(tail(curr_training_df_sp[ , sp_input_curr], max_E_curr + 1))))
  # loop over phosphate levels
  for (k in 1:length(end_points)) {
    print(k)
    # create lagged data frame for GP-EDM
    curr_training_df_sp <- curr_training_df[1:end_points[k], ]
    # create temperature values
    start_month <- unlist(c(tail(curr_training_df_sp[ , "Month"], max_E_curr + 1)))[1]
    month_order <- c(names(temp_month_avg)[start_month:length(temp_month_avg)],
                     names(temp_month_avg)[1:(start_month-1)])
    temp_month_avg_ordered <- temp_month_avg[as.numeric(month_order)]
    temp_month_avg_ordered <- rep(temp_month_avg_ordered, length = t_max_full + max_E_curr)
    # log abundances
    if (log_abund) {
      curr_training_df_sp[ , sp_input[j]][curr_training_df_sp[ , sp_input[j]] == 0] <- 1
      curr_training_df_sp[ , sp_input[j]] <- log(curr_training_df_sp[ , sp_input[j]])
    }
    # create lags
    training_lags <- as.data.frame(makelags(data = curr_training_df_sp, y = colnames(curr_training_df_sp)[-1], 
                                            time = "Time", E = max_E_sp, tau = tau))
    full_training_lags <- cbind(curr_training_df_sp, training_lags)
    # fit GP-EDM on training data
    GP_fit_sp <- fitGP(data = full_training_lags, y = sp_input_curr, 
                       x = inputs, scaling = "local",
                       fixedpars = fixed_pars, predictmethod = "loo")
    print(GP_fit_sp$pars)
    
    # examine conditional responses of fitted model ------------------------------ 
    if (sp_input_curr == "Green algae (large)") {
      # extract conditional responses
      abund_1 <- seq(500, 80000, by = 500)
      abund_2 <- mean(full_training_lags$`Green algae (large)_2`, na.rm = TRUE)
      abund_3 <- mean(full_training_lags$`Green algae (large)_3`, na.rm = TRUE)
      abund_4 <- mean(full_training_lags$`Green algae (large)_4`, na.rm = TRUE)
      temp_1 <- mean(full_training_lags$Temperature_1, na.rm = TRUE)
      temp_2 <- mean(full_training_lags$Temperature_2, na.rm = TRUE)
      # create grid with control parameter and abundances
      cond_resp_df <- expand_grid(abund_1, abund_2, abund_3, abund_4, phosphate_levels, temp_1, temp_2)
      names(cond_resp_df) <- inputs
      cond_resp_df$abund <- predict(GP_fit_sp, newdata = cond_resp_df)$outsampresults$predmean
      names(cond_resp_df) <- c(inputs, "Green algae (large)")
      cond_resp_df$Phosphate_1 <- as.factor(cond_resp_df$Phosphate_1)
      # plot
      fig <- ggplot() +
        geom_line(data = cond_resp_df, aes(x = `Green algae (large)_1`, y = `Green algae (large)`, 
                                           group = `Phosphate_1`, color = `Phosphate_1`), size = 1) + 
        scale_color_viridis(name = "Phosphate\nconcentration", discrete = TRUE) +
        xlab(label = "Current abundance") +
        ylab(label = "Future abundance") +
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
      if (save_plot) {
        ggsave(paste("figs/fig_", lake, "_", sp_input[j], "_cond_resp_", temporal_resolution, 
                     "_end_point_", end_points[k], "_phi_control_", fixed_phi_control, ".pdf", sep = ""), 
               fig, width = 14, height = 9, units = "cm")
      }
    }
    
    # residual plots ------------------------------ 
    plot_df <- full_training_lags
    # compute prediction errors
    plot_df$obs <- GP_fit_sp$outsampresults$obs
    plot_df$predmean <- GP_fit_sp$outsampresults$predmean
    plot_df$residuals <- plot_df$obs - plot_df$predmean
    # plot residual autocorrelation
    residuals_df <- as.data.frame(makelags(data = plot_df, y = "residuals", 
                                           time = "Time", E = max_E_sp, tau = tau))
    auto_cor <- apply(X = residuals_df, MARGIN = 2, FUN = function(x, y) cor(x, y, use = "complete.obs"), 
                      y = plot_df$residuals)
    auto_cor_df <- data.frame(lag = 1:max_E_sp, auto_cor = auto_cor)
    fig <- ggplot(data = auto_cor_df) +
      geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
      geom_bar(aes(x = lag, y = auto_cor), stat="identity", size = 1, width = 0.5,
               color = "black", fill = "gray70") +
      scale_x_continuous(n.breaks = max_E_sp) +
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
    residuals_df <- plot_df[ , c(best_inputs, "residuals")]
    residuals_df <- gather(residuals_df, "lag", "value", names(residuals_df)[-ncol(residuals_df)])
    fig <- ggplot(data = residuals_df) +
      geom_point(aes(x = residuals, y = value), size = 2) +
      facet_wrap(~lag) +
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
    
    # perform forecasts for each phosphate level ------------------------------ 
    full_test_set <- data.frame()
    for (l in 1:length(phosphate_grid)) {
      # build test set to make predictions
      test_df <- data.frame(time = 1:(max_E_curr+1),
                            x = x_init,
                            phosphate = phosphate_grid[l],
                            temperature = temp_month_avg_ordered[1:(max_E_curr+1)])
      names(test_df) <- c("Time", sp_input_curr, "Phosphate", "Temperature")
      # iterate forecasts for t_max_full time steps
      for (m in 1:t_max_full) {
        # build data frame to make forecast
        test_lags <- makelags(data = test_df, y = colnames(test_df)[-1],
                              time = "Time", E = max_E_curr, tau = tau, forecast = TRUE)
        # make forecast
        forecast_sp <- predict(GP_fit_sp, newdata = test_lags)$outsampresults$predmean
        if (forecast_sp < 0) {
          forecast_sp <- 0
        }
        test_df[m + max_E_curr + 1, sp_input_curr] <- forecast_sp
        # update time and control level in test set
        test_df[m + max_E_curr + 1, "Time"] <- test_df[m + max_E_curr, "Time"] + 1
        test_df[m + max_E_curr + 1, "Phosphate"] <- phosphate_grid[l]
        test_df[m + max_E_curr + 1, "Temperature"] <- temp_month_avg_ordered[m + max_E_curr + 1]
      }
      # update starting abundance values
      x_init <- as.numeric(unlist(c(tail(test_df[ , sp_input_curr], max_E_curr + 1))))
      # extract tail from test set
      test_df <- tail(test_df, t_max)
      # merge with full data frame
      full_test_set <- rbind(full_test_set, test_df)
    }
    # exponentiate abundances
    if (log_abund) {
      full_test_set[ , sp_input_curr] <- exp(full_test_set[ , sp_input_curr])
    }
    # create results data frame
    plot_df <- full_test_set
    plot_df$abundance <- plot_df[ , sp_input_curr]
    # defining plot axis
    if (sp_input_curr == "Green algae (large)") {
      max_abund <- 100000
    }
    # plot predicted bifurcation diagram together with data points
    plot_training_df <- training_df
    plot_training_df$type <- "unobserved"
    plot_training_df$type[1:end_points[k]] <- "observed"
    fig <- ggplot() +
      geom_point(data = subset(plot_training_df, type == "unobserved"), 
                 aes(x = Phosphate, y = `Green algae (large)`), 
                 size = 2.5, alpha = 0.8, color = "#DEDEDE") +
      geom_point(data = subset(plot_training_df, type == "observed"), 
                 aes(x = Phosphate, y = `Green algae (large)`), 
                 size = 2.5, alpha = 0.8, color = "#454545") +
      geom_point(data = plot_df, aes(x = Phosphate, y = abundance), size = 0.9, color = "#CB181D") +
      geom_vline(xintercept = tail(curr_training_df_sp$Phosphate, 1), linetype = "dashed", size = 1) +
      geom_vline(xintercept = breakpoint_df$phosphate, size = 1) +
      scale_x_reverse() +
      scale_y_continuous(limits = c(0, max_abund)) +
      xlab(label = "Phosphate concentration") +
      ylab(label = "Population abundance") +
      theme_bw() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(size = 1),
            title = element_text(size = 22),
            axis.title = element_text(size = 22),
            axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18),
            legend.position = "none")
    # plot info
    if (filter_phosphate) {
      filter_phosphate_name <- "yes_filter_phosphate"
    } else {
      filter_phosphate_name <- "no_filter_phosphate"
    }
    if (interpolate) {
      interpolate_name <- "yes_interpolate"
    } else {
      interpolate_name <- "no_interpolate"
    }
    if (log_abund) {
      log_abund_name <- "yes_log_abund"
    } else {
      log_abund_name <- "no_log_abund"
    }
    if (!is.na(max_E_temperature)) {
      temperature_lag <- paste("E_temperature", E_temperature, sep = "_")
    } else {
      temperature_lag <- "E_temperature_NA"
    }
    sp_lag <- paste("E_sp", E_sp, sep = "_")
    phosphate_lag <- paste("E_phosphate", E_phosphate, sep = "_")
    # save plot
    if (save_plot) {
      ggsave(paste("figs/fig5/fig5_", lake, "_", sp_input[j], "_predicted_bif_diag_", temporal_resolution, 
                   "_end_point_", end_points[k], "_phi_control_", fixed_phi_control, "_",
                   filter_phosphate_name, "_", interpolate_name, "_", log_abund_name, "_",
                   sp_lag, "_", temperature_lag, "_", phosphate_lag, "_smaller_window.pdf", sep = ""), 
             fig, width = 15, height = 10, units = "cm")
    }
    
    # compute power spectrum and Lyapunov exponent at certain phosphate levels ------------------------------ 
    if (sp_input_curr == "Green algae (large)") {
      # tolerance to subset data frame
      tol <- 10^(-10)
      # take subset of data frame
      sub_plot_df <- subset(plot_df, (abs(plot_df$Phosphate - phosphate_levels[1]) < tol) | 
                              (abs(plot_df$Phosphate - phosphate_levels[2]) < tol) |
                              (abs(plot_df$Phosphate - phosphate_levels[3]) < tol) |
                              (abs(plot_df$Phosphate - phosphate_levels[4]) < tol) |
                              (abs(plot_df$Phosphate - phosphate_levels[5]) < tol) |
                              (abs(plot_df$Phosphate - phosphate_levels[6]) < tol))
      # data frame to store power spectrum
      power_spectrum_df <- data.frame()
      # data frame to store largest Lyapunov exponent
      max_le_predicted_df <- c()
      # loop over chosen phosphate levels
      for (n in 1:length(phosphate_levels)) {
        # current data frame with lags
        curr_plot_df <- subset(sub_plot_df, (abs(sub_plot_df$Phosphate - phosphate_levels[n]) < tol))
        curr_plot_df_lags <- as.data.frame(makelags(y = curr_plot_df[ , c("Green algae (large)", "Temperature")], 
                                                    time = "time", E = E_sp, tau = tau))
        curr_plot_df_lags <- cbind(curr_plot_df, curr_plot_df_lags)
        # compute power spectrum from predicted time series
        power_spectrum_curr <- power_spectrum(ts = curr_plot_df[ , c("Time", "Green algae (large)")], lam = 0.01, 
                                              scale = TRUE, trim = TRUE, plot = FALSE)[[1]]
        power_spectrum_curr$phosphate <- phosphate_levels[n]
        power_spectrum_df <- rbind(power_spectrum_df, power_spectrum_curr)
        # obtain gradient by fitting a GPEDM model to the predicted time series
        which_phosphate <- which(grepl("Phosphate", inputs))
        GP_fit_sp_predicted <- fitGP(data = curr_plot_df_lags, y = sp_input_curr, 
                                     x = inputs[-which_phosphate], scaling = "local",
                                     fixedpars = GP_fit_sp$pars[-c(which_phosphate, length(GP_fit_sp$pars))],
                                     predictmethod = "loo", returnGPgrad = TRUE)
        grad_df <- GP_fit_sp_predicted$GPgrad
        # remove NAs
        grad_df <- na.omit(grad_df)
        # create dummy matrix
        shifted_I <- cbind(diag(E_sp-1), rep(0, E_sp-1))
        # Jacobian matrix for each state 
        J <- list()
        for (p in 1:nrow(grad_df)) {
          # current GPEDM model gradient
          grad <- as.matrix(grad_df[p, 1:E_sp])
          J[[p]] <- rbind(grad, shifted_I)
        }
        # compute largest Lyapunov exponent from sequence of Jacobian matrices
        max_le_predicted <- max(le_qr_decomp(J))
        max_le_predicted_curr_df <- data.frame(max_le_predicted = max_le_predicted,
                                               phosphate = phosphate_levels[n])
        max_le_predicted_df <- rbind(max_le_predicted_df, max_le_predicted_curr_df)
      }
      # plot power spectrum
      fig <- ggplot(data = power_spectrum_df, aes(x = frequency, y = power, color = phosphate)) +
        geom_line(size = 1) +
        facet_wrap(~phosphate, labeller = "label_both", nrow = length(phosphate_levels)/2) +
        xlab(label = "Frequency") +
        ylab(label = "Power") +
        scale_color_viridis() +
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
              legend.key.height = unit(0.5, 'cm'),
              legend.key.width = unit(1.2, 'cm'))
      if (save_plot) {
        ggsave(paste("figs/fig_", lake, "_", sp_input[j], "_power_spectrum_", temporal_resolution, 
                     "_end_point_", end_points[k], "_phi_control_", fixed_phi_control, "_",
                     filter_phosphate_name, "_", interpolate_name, "_", log_abund_name, "_",
                     sp_lag, "_", temperature_lag, "_", phosphate_lag, "_smaller_window.pdf", sep = ""), 
               fig, width = 14, height = 14, units = "cm")
      }
      # plot largest Lyapunov exponent
      fig <- ggplot(data = max_le_predicted_df, aes(x = phosphate, y = max_le_predicted, color = phosphate)) +
        geom_point(size = 5) +
        geom_hline(yintercept = 0, size = 0.5) +
        scale_color_viridis() +
        xlab(label = "Phosphate concentration") +
        ylab(label = "Largest Lyapunov exponent") +
        scale_x_reverse() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_rect(size = 1),
              axis.title = element_text(size = 18),
              axis.text.x = element_text(size = 14),
              axis.text.y = element_text(size = 14),
              strip.text = element_text(size = 14),
              strip.background = element_rect(fill = "white", size = 1),
              legend.position = "none")
      if (save_plot) {
        ggsave(paste("figs/fig_", lake, "_", sp_input[j], "_lyapunov_exponent_", temporal_resolution, 
                     "_end_point_", end_points[k], "_phi_control_", fixed_phi_control, "_",
                     filter_phosphate_name, "_", interpolate_name, "_", log_abund_name, "_",
                     sp_lag, "_", temperature_lag, "_", phosphate_lag, "_smaller_window.pdf", sep = ""), 
               fig, width = 12, height = 10, units = "cm")
      }
    }
  }
}