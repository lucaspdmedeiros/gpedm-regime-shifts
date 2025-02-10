# Code that explores the Swiss lakes plankton data sets to:
# 1) quantify the trend in phosphate over time
# 2) plot time series for species and environmental variables
# 3) perform GP-EDM predictions to check which species are more predictable

# cleaning wd, loading functions and packages ------------------------------ 
rm(list = ls(all = TRUE))
source("code/functions/ts_breakpoint.R")
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

# for each data set compute trend in phosphate and plot time series ------------------------------ 
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
# temporal resolution
temporal_resolution <- "month"
# lake to use
lake <- "zurich"
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
  summ_df_long <- gather(data = summ_df, key = "variable", value = "value", -Time)
  summ_df$Year <- summ_df$year
  summ_df$Month <- summ_df$month
}
# compute trend in phosphate
phosphate <- as.numeric(scale(summ_df$Phosphate))
time <- 1:nrow(summ_df)
mod <- lm(phosphate ~ time)
print(mod$coefficients)
print(summary(mod)$adj.r.squared)

# perform leave-one-out predictions for each species ------------------------------
# data frame to store R2 results
results_df <- data.frame()
# species names to use
sp_input <- c("Cyanobacteria", "Green algae (small)", "Green algae (large)", "Diatoms (small)", "Diatoms (large)",
              "Gold algae (small)", "Gold algae (large)", "Cryptophytes (small)", "Cryptophytes (large)", 
              "Mixotrophic flagellates", "Large herbivores", "Omnivores", "Invertebrate predators")
# loop over lakes
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
                            breakpoint_date = training_df$Time[sp_breakpoints],
                            breakpoint_ss = ss_breakpoints,
                            phosphate = training_df$Phosphate[sp_breakpoints])
# plot time series
summ_df_long <- gather(data = subset(training_df, select = -c(Year, Month, year, month, month_year, quarter_year)), 
                       key = "variable", value = "value", -Time)
fig <- ggplot() +
  geom_line(data = summ_df_long, aes(x = Time, y = value, 
                                     group = variable, color = variable), size = 0.7) +
  geom_vline(data = breakpoint_df, aes(xintercept = breakpoint_date, 
                                       group = variable), linetype = "dashed", size = 1) +
  xlab(label = "Time") +
  ylab(label = "Value") +
  facet_wrap(~variable, ncol = 3, scales = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1.2),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = "white", size = 1.2),
        axis.title = element_text(size = 17),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        legend.position = "none",
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "in"))
# save plot
if (save_plot) {
  if (filter_phosphate) {
    ggsave(paste("figs/fig5/fig_lake_plankton_", lake, "_ts_", temporal_resolution, "_filter_phosphate.pdf", sep = ""), 
           fig, width = 16, height = 10, units = "in")
  } else {
    ggsave(paste("figs/fig5/fig_lake_plankton_", lake, "_ts_", temporal_resolution, ".pdf", sep = ""), 
           fig, width = 16, height = 10, units = "in")
  }
}
# log transform abundances
if (log_abund) {
  for (j in 1:length(sp_input)) {
    training_df[ , sp_input[j]][training_df[ , sp_input[j]] == 0] <- 1
  }
  training_df[ , sp_input] <- 
    log(training_df[ , sp_input])
}
# loop over species
for (j in 1:length(sp_input)) {
  # take only a subset of all data for that species
  sp_input_curr <- sp_input[j]
  breakpoint <- 106
  curr_training_df_sp <- training_df[1:breakpoint, c("Time", "Month", "Phosphate", "Temperature", sp_input_curr)]
  
  # determine best embedding dimension ------------------------------ 
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
  if (is.na(fixed_phi_control)) {
    best_inputs <- input_list[[which.max(R2)]]
  } else {
    best_inputs <- input_list[[which.max(R2)]]
  }
  print(best_inputs)
  print(max(R2))
  
  # plot observations and predictions for best inputs ------------------------------ 
  # best input names
  inputs <- best_inputs
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
                     fixedpars = fixed_pars, predictmethod = "loo")
  # add results to data frame
  phi_sp <- GP_fit_sp$pars[1:E_sp]
  if (length(phi_sp) > 1) {
    phi_sp <- paste(round(phi_sp, 3), collapse = ", ")
  } else {
    phi_sp <- round(phi_sp, 3)
  }
  phi_phosphate <- GP_fit_sp$pars[(E_sp+1):(E_sp+E_phosphate)]
  if (length(phi_phosphate) > 1) {
    phi_phosphate <- paste(round(phi_phosphate, 3), collapse = ", ")
  } else {
    phi_phosphate <- round(phi_phosphate, 3)
  }
  phi_temperature <- GP_fit_sp$pars[(E_sp+E_phosphate+1):
                                      (E_sp+E_phosphate+E_temperature)]
  if (length(phi_temperature) > 1) {
    phi_temperature <- paste(round(phi_temperature, 3), collapse = ", ")
  } else {
    phi_temperature <- round(phi_temperature, 3)
  }
  curr_df <- data.frame(lake = lake,
                        species = sp_input[j],
                        E_sp = E_sp,
                        phi_sp = phi_sp,
                        E_phosphate = E_phosphate,
                        phi_phosphate = phi_phosphate,
                        E_temperature = E_temperature,
                        phi_temperature = phi_temperature,
                        R2 = max(R2))
  results_df <- rbind(results_df, curr_df)
  # plot true and predicted abundances
  full_training_lags$obs <- GP_fit_sp$outsampresults$obs
  full_training_lags$predmean <- GP_fit_sp$outsampresults$predmean
  plot_df <- full_training_lags
  plot_df$time <- training_df[1:breakpoint, ]$Time
  if (log_abund) {
    plot_df$obs <- exp(plot_df$obs)
    plot_df$predmean <- exp(plot_df$predmean)
  }
  fig <- ggplot(data = plot_df) +
    geom_line(aes(x = time, y = obs), size = 1, color = "black") +
    geom_line(aes(x = time, y = predmean), size = 1, color = "#E41A1C") +
    xlab(label = "Time") +
    ylab(label = "Population abundance") +
    ggtitle(label = paste("Lake ", lake, ", ", sp_input[j], ", R2 = ", round(max(R2), 2), sep = ""),
            subtitle = paste("Inputs = ", paste(inputs, collapse = ", "), sep = "")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1.2),
          strip.text = element_text(size = 16),
          strip.background = element_rect(fill = "white", size = 1.2),
          plot.title = element_text(size = 14),
          plot.subtitle = element_text(size = 10),
          axis.title = element_text(size = 18),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14))
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
    ggsave(paste("figs/fig5/fig_lake_plankton_", lake, "_", sp_input[j], "_loo_predictions_", temporal_resolution, 
                 "_end_point_", breakpoint, "_phi_control_", fixed_phi_control, "_",
                 filter_phosphate_name, "_", interpolate_name, "_", log_abund_name, ".pdf", sep = ""), 
           fig, width = 7, height = 3, units = "in")
  }
}
# save tables with results
write.csv(x = breakpoint_df, file = paste("figs/fig5/breakpoints_lake_plankton_", temporal_resolution, "_",
                                          interpolate_name, "_", log_abund_name, 
                                          ".csv", sep = ""), row.names = FALSE)
write.csv(x = results_df, file = paste("figs/fig5/loo_predictions_lake_plankton_", temporal_resolution, 
                                       "_end_point_", breakpoint, "_phi_control_", fixed_phi_control, "_",
                                       filter_phosphate_name, "_", interpolate_name, "_", log_abund_name, 
                                       "_max_E_sp_", max_E_sp, "_max_E_phosphate_", max_E_phosphate,
                                       "_max_E_temperature_", max_E_temperature, ".csv", sep = ""), row.names = FALSE)
