##################################################################################################
orderly2::orderly_strict_mode()
orderly2::orderly_parameters(n_samples = 100, n_particles = 8, n_chains = 2, step_size_var = 0.03, dt = 1,
                             restart = FALSE, rerun_every = 100,
                             include_NAs = TRUE)

params <- c("alpha", "beta", "gamma", "sigma", "asc_rate", "dispersion")
orderly2::orderly_artefact(description = "The posterior samples", "fitting_samples.rds")
orderly2::orderly_artefact(description = "Log_likelihood of chains", "ll_plot.png")
orderly2::orderly_artefact(description = "The posterior density of all chains", "Posterior_densities//00_all_plot.png")
orderly2::orderly_artefact(description = "The posterior sampling trajectories of all chains", "Posterior_chains//00_all_plot.png")
orderly2::orderly_artefact(description = "The parameters used", "parameters_used.txt")
for(i in 1:length(params)){
  orderly2::orderly_artefact(description = "Posterior density", sprintf("Posterior_densities//%s_plot.png", params[i]))
  orderly2::orderly_artefact(description = "Posterior trajectories", sprintf("Posterior_chains//%s_plot.png", params[i]))
}
orderly2::orderly_artefact(description = "MCMC diagnostics", "mcmc_diagnostics.txt")

orderly2::orderly_artefact(description = "Fits to data", "Fit_to_data/Infected_Herds_1.png")
orderly2::orderly_artefact(description = "Fits to data", "Fit_to_data/Infected_Proportion_Herds_1.png")
orderly2::orderly_artefact(description = "Fits to data", "Fit_to_data/Infected_Herds_2.png")
orderly2::orderly_artefact(description = "Fits to data", "Fit_to_data/Infected_Proportion_Herds_2.png")

orderly2::orderly_artefact(description = "Fits to data", "Fit_to_data/Declared_Infected_Herds_1.png")
orderly2::orderly_artefact(description = "Fits to data", "Fit_to_data/Declared_Infected_Proportion_Herds_1.png")
orderly2::orderly_artefact(description = "Fits to data", "Fit_to_data/Declared_Infected_Herds_2.png")
orderly2::orderly_artefact(description = "Fits to data", "Fit_to_data/Declared_Infected_Proportion_Herds_2.png")
##################################################################################################
library(cowflu)
library(dust2)
library(monty)
library(ggplot2)
library(gridExtra)
library(tictoc)
library(coda)
library(dplyr)
library(tidyr)
tic()
##################################################################################################
## Set parameters
pars <- cowflu:::cowflu_inputs(
  alpha = 0.05,
  beta = 1.45,
  gamma = 1.4,
  sigma = 1.5,
  asc_rate = 0.6,
  dispersion = 1,
  cowflu:::cowflu_fixed_inputs(
    n_herds_per_region = cowflu:::usda_data$n_herds_per_region,
    p_region_export = cowflu:::movement$p_region_export,
    p_cow_export = cowflu:::movement$p_cow_export,
    n_cows_per_herd = cowflu:::usda_data$n_cows_per_herd,
    movement_matrix = cowflu:::movement$movement_matrix,
    time_test = 19, #Day 136
    start_herd = 26940, #26804 is where Texas starts
    start_count = 5,
    condition_on_export = TRUE,
    likelihood_choice = "incidence"))

## Set priors
# prior <- monty::monty_dsl({
#   alpha ~ Beta(a = 1, b = 25)
#   beta ~ Uniform(min = 0.05, max = 4) #maybe 1 and 2.5
#   gamma ~ Uniform(min = 0.05, max = 4) #0 and 1
#   sigma ~ Uniform(min = 0.05, max = 5) #0 and 2
#   asc_rate ~ Beta(a = 5, b = 1)
#   #dispersion ~ Exponential(mean = 1)
# })

prior <- monty::monty_dsl({
  alpha ~ Uniform(min = 0, max = 0.1)
  beta ~ Uniform(min = 0.05, max = 3) #maybe 1 and 2.5
  gamma ~ Uniform(min = 0.05, max = 2) #0 and 1
  sigma ~ Uniform(min = 0.05, max = 2) #0 and 2
  asc_rate ~ Beta(a = 1, b = 1)
  #dispersion ~ Exponential(mean = 1)
})

## Pack the priors
pars_fixed <- pars[-(16:21)]
prior_packer <- monty::monty_packer(c("alpha", "beta", "gamma", "sigma", "asc_rate", "dispersion"), fixed = pars_fixed)

## With this packer we can convert from a list of name-value pairs suitable for
## initialising a dust2 system into a vector of parameters suitable for use with monty:
prior_packer$pack(pars)

## Load data to fit to
data_outbreaks <- cowflu:::process_data_incidence(cowflu:::outbreaks_data$weekly_outbreaks_data)
## Add extra "NA" weeks to output the model fit to these.
if(include_NAs){
## Generate rows for weeks 1-13
weeks_1_to_13 <- data.frame(
  day = seq(7, 7 * 13, by = 7),   # Days: 7, 14, 21, ..., 91
  week = 1:13,                    # Weeks 1 to 13
  positive_tests = I(lapply(1:13, function(x) rep(NA_real_, 48)))  # 48 NAs in each row
)

## Generate rows to fill to week 50
last_day <- max(data_outbreaks$day)  # Find the last day in the current data
last_week <- max(data_outbreaks$week)  # Find the last week in the current data

## Create the data frame to week 50
weeks_to_50 <- data.frame(
  day = seq(last_day + 7, last_day + 7 * (50 - last_week), by = 7),
  week = (last_week + 1):50,
  positive_tests = I(lapply((last_week + 1):50, function(x) rep(NA_real_, 48)))
)

## Append the new rows to the original data
data_outbreaks <- rbind(weeks_1_to_13, data_outbreaks, weeks_to_50)
}

set.seed(1)

## Build a particle filter
data_week <- dust2::dust_filter_data(data_outbreaks, time = "week")
filter <- dust2::dust_filter_create(cowflu:::cows(), 0, #0 is "time_start"
                                    data_week, n_particles = n_particles, n_threads = 32,
                                    dt=dt)
## Build a likelihood
## History save options are:
##  "S_herd"  "S_region" "E_herd" "E_region"  "I_herd"  "I_region"  "R_herd"  "R_region"
##  "outbreak_herd" "outbreak_region"
likelihood <- dust2::dust_likelihood_monty(filter, prior_packer,
                                           save_state = FALSE,
                                           save_history = c("outbreak_region", "infected_herds_region"))

## We can combine the prior and the likelihood to create a posterior:
posterior <- prior + likelihood

##Build the sampler
if(restart){
  sampler <- monty::monty_sampler_random_walk(diag(length(params) ) * step_size_var,
                                              rerun_every = rerun_every,
                                              rerun_random = TRUE)  #0.02 was baseline
}else{
  sampler <- monty::monty_sampler_random_walk(diag(length(params) ) * step_size_var)  #0.02 was baseline
}


## Run the samples
samples <- monty::monty_sample(posterior, sampler, n_samples, n_chains = n_chains,
                               initial = prior_packer$pack(pars))
saveRDS(samples, "fitting_samples.rds")
##################################################################################################
## Fit to data plotting
dir.create("Fit_to_data")
plot_indices <- list(first = c(1,25), second = c(26,48))

## First, the TRUE number of infected herds per state:
for(i in 1:length(plot_indices)){
  regions_to_plot <- plot_indices[[i]]
  regions_to_plot <- regions_to_plot[1]:regions_to_plot[2]
  n_states <- length(regions_to_plot)

  ## Extract the  data across all chains and iterations
  I_data <- samples$observations$history[49:96,,,]
  n_timepoints <- dim(I_data)[2]
  ## Reshape the array by collapsing dimensions 3 and 4
  I_data <- array(I_data, dim = c(48, n_timepoints, n_samples * n_chains))
  I_data <- I_data[regions_to_plot, , ]
  ## Initialize an empty list to store data frames for each state
  data_list <- list()
  ## Loop over each state and calculate summary statistics
  for (state in 1:n_states) {
    ## Extract data for the current state
    state_data <- I_data[state, , ]
    ## Calculate mean, lower, and upper CI across particles for each time point
    state_summary <- apply(state_data, 1, function(x) {
      mean_val <- mean(x)
      ci_low <- quantile(x, 0.025)
      ci_high <- quantile(x, 0.975)
      return(c(mean = mean_val, lower_ci = ci_low, upper_ci = ci_high))
    })
    ## Convert to a data frame
    state_df <- as.data.frame(t(state_summary))
    ## Add time and state information
    state_df$Time <- data_week$week
    state_df$US_state <- state
    ## Append to the list
    data_list[[state]] <- state_df
  }

  ## Combine all state data frames into one data frame
  result_df <- dplyr::bind_rows(data_list)

  ## Rename columns
  colnames(result_df) <- c("mean_infected", "lower_ci_infected", "upper_ci_infected", "Time", "US_state")
  ## Add number of herds per state:
  result_df$total_herds <- cowflu:::usda_data$n_herds_per_region[result_df$US_state]
  ## Replace the US_state numbers with actual names
  result_df$US_state <- factor(result_df$US_state, levels = 1:n_states, labels = cowflu:::usda_data$US_States[regions_to_plot])

  ## Create the plot
  ggplot(result_df, aes(x = Time, y = mean_infected, group = US_state)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower_ci_infected, ymax = upper_ci_infected), alpha = 0.2) +
    facet_wrap(~ US_state, ncol = 5, scales = "free_y") +
    theme_minimal() +
    labs(x = "Time (Weeks)", y = "Infected Herds", title = "Infected Herds per state") +
    theme(strip.text = element_text(size = 8)) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) -> my_plot

  ## Create the proportion plot
  ggplot(result_df, aes(x = Time, y = mean_infected/total_herds, group = US_state)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower_ci_infected/total_herds, ymax = upper_ci_infected/total_herds), alpha = 0.2) +
    facet_wrap(~ US_state, ncol = 5, scales = "free_y") +
    theme_minimal() +
    labs(x = "Time (Weeks)", y = "Infected Proportion of Herds", title = "Infected Proportion of Herds by state") +
    theme(strip.text = element_text(size = 8)) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) -> my_plot_prop

  file_path <- sprintf("Fit_to_data/Infected_Herds_%s.png", i)
  ggsave(filename = file_path, plot = my_plot, width = 1280/96, height = 720/96, dpi = 200)
  file_path <- sprintf("Fit_to_data/Infected_Proportion_Herds_%s.png", i)
  ggsave(filename = file_path, plot = my_plot_prop, width = 1280/96, height = 720/96, dpi = 200)
}

######################
### Now do the same, but for outbreaks declared, and show the fit to data.
## We will need to reorder the data into a data frame:
real_outbreaks_data <- data_week %>%
  tidyr::unnest_longer(positive_tests) %>%
  mutate(state = rep(cowflu:::usda_data$US_States, times = nrow(data_week)),
         total_herds = rep(cowflu:::usda_data$n_herds_per_region, times = nrow(data_week)))
real_outbreaks_data <- real_outbreaks_data[,c(2,3,4,5)]
colnames(real_outbreaks_data) <- c("Time", "positive_tests", "US_state", "total_herds")

plot_indices <- list(first = c(1,25), second = c(26,48))
for(i in 1:length(plot_indices)){
  regions_to_plot <- plot_indices[[i]]
  regions_to_plot <- regions_to_plot[1]:regions_to_plot[2]
  n_states <- length(regions_to_plot)

  ## Extract the  data across all chains and iterations
  I_data <- samples$observations$history[1:48,,,]
  n_timepoints <- dim(I_data)[2]
  ## Reshape the array by collapsing dimensions 3 and 4
  I_data <- array(I_data, dim = c(48, n_timepoints, n_samples * n_chains))
  I_data <- I_data[regions_to_plot, , ]
  ## Initialize an empty list to store data frames for each state
  data_list <- list()
  ## Loop over each state and calculate summary statistics
  for (state in 1:n_states) {
    ## Extract data for the current state
    state_data <- I_data[state, , ]
    ## Calculate mean, lower, and upper CI across particles for each time point
    state_summary <- apply(state_data, 1, function(x) {
      mean_val <- mean(x)
      ci_low <- quantile(x, 0.025)
      ci_high <- quantile(x, 0.975)
      return(c(mean = mean_val, lower_ci = ci_low, upper_ci = ci_high))
    })
    ## Convert to a data frame
    state_df <- as.data.frame(t(state_summary))
    ## Add time and state information
    state_df$Time <- data_week$week
    state_df$US_state <- state
    ## Append to the list
    data_list[[state]] <- state_df
  }

  ## Combine all state data frames into one data frame
  result_df <- dplyr::bind_rows(data_list)

  ## Rename columns
  colnames(result_df) <- c("mean_infected", "lower_ci_infected", "upper_ci_infected", "Time", "US_state")
  ## Add number of herds per state:
  result_df$total_herds <- cowflu:::usda_data$n_herds_per_region[result_df$US_state]
  ## Replace the US_state numbers with actual names
  result_df$US_state <- factor(result_df$US_state, levels = 1:n_states, labels = cowflu:::usda_data$US_States[regions_to_plot])

  ## Filter the data we fit to:
  specific_plot_outbreaks_data <- filter(real_outbreaks_data, US_state %in% cowflu:::usda_data$US_States[regions_to_plot])

  ## Create the plot
  ggplot(result_df, aes(x = Time, y = mean_infected, group = US_state)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower_ci_infected, ymax = upper_ci_infected), alpha = 0.2) +
    geom_point(data = specific_plot_outbreaks_data,
               aes(x = Time, y = positive_tests, group = US_state), col = "red") +
    facet_wrap(~ US_state, ncol = 5, scales = "free_y") +
    theme_minimal() +
    labs(x = "Time (Weeks)", y = "Declared Infected Herds", title = "Declared Infected Herds per state") +
    theme(strip.text = element_text(size = 8)) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) -> my_plot

  ## Create the proportion plot
  ggplot(result_df, aes(x = Time, y = mean_infected/total_herds, group = US_state)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower_ci_infected/total_herds, ymax = upper_ci_infected/total_herds), alpha = 0.2) +
    geom_point(data = filter(specific_plot_outbreaks_data,
                             !is.na(positive_tests)),
               aes(x = Time, y = positive_tests/total_herds, group = US_state), col = "red") +
    facet_wrap(~ US_state, ncol = 5, scales = "free_y") +
    theme_minimal() +
    labs(x = "Time (Weeks)", y = "Declared Infected Proportion of Herds", title = "Declared Infected Proportion of Herds by state") +
    theme(strip.text = element_text(size = 8)) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) -> my_plot_prop

  file_path <- sprintf("Fit_to_data/Declared_Infected_Herds_%s.png", i)
  ggsave(filename = file_path, plot = my_plot, width = 1280/96, height = 720/96, dpi = 200)
  file_path <- sprintf("Fit_to_data/Declared_Infected_Proportion_Herds_%s.png", i)
  ggsave(filename = file_path, plot = my_plot_prop, width = 1280/96, height = 720/96, dpi = 200)
}
##################################################################################################
## Diagnostic Plotting

dir.create("Posterior_densities")
dir.create("Posterior_chains")

ll_plot <- cowflu:::plot_chains_ll(samples)
ggsave(filename = "ll_plot.png", plot = ll_plot, width = 8, height = 6, units = "in", dpi = 300)
posterior_plot <- cowflu:::plot_param_posterior(samples, one_panel = TRUE)
ggsave(filename = "Posterior_densities//00_all_plot.png", plot = posterior_plot, width = 6*1.25, height = 8*1.25, units = "in", dpi = 300)
posterior_plot <- cowflu:::plot_param_posterior(samples)
chains_plot <- cowflu:::plot_param_traj(samples, one_panel = TRUE)
ggsave(filename = "Posterior_chains//00_all_plot.png", plot = chains_plot, width = 6*1.25, height = 8*1.25, units = "in", dpi = 300)
chains_plot <- cowflu:::plot_param_traj(samples)

for(i in 1:length(params)){
  ggsave(filename = sprintf("Posterior_densities//%s_plot.png", params[i]), plot = posterior_plot[[i]], width = 8, height = 6, units = "in", dpi = 300)
  ggsave(filename = sprintf("Posterior_chains//%s_plot.png", params[i]), plot = chains_plot[[i]], width = 8, height = 6, units = "in", dpi = 300)
}

dev.off()

## Convert to coda object:
samples_coda <- coda::as.mcmc.list(samples)
#coda::gelman.plot(samples_coda)
## Calculate Effective Sample Sizes (ESS)
ess <- coda::effectiveSize(samples_coda)
results <- data.frame(ESS = ess)
## Calculate PSRF (Potential Scale Reduction Factor)
## This can fail for very small trial runs:
tryCatch({
  # Code that may throw an error
  psrf <- coda::gelman.diag(samples_coda)$psrf
  results <- data.frame(ESS = ess, PSRF = psrf[,1], PSRF_95CI = psrf[,2])
}, error = function(e) {
  # Code to run if there's an error
  cat("An error occurred, the leading minor of order 3 is not positive.\n")
})

## Write the results to a .txt file
write.table(results, file = "mcmc_diagnostics.txt", sep = "\t", col.names = TRUE, row.names = TRUE)

##################################################################################################
duration <- toc()

## Print an output .txt of the parameters used:
param_string <- sprintf("duration ran: %s mins\n
  n_samples: %s \n
  n_particles: %s \n
  n_chains: %s \n
  dt: %s \n
  step_size_var: %s \n
  restarts:  %s \n
  restart_every: %s ",   (duration$toc - duration$tic)/60,
                        n_samples,
                        n_particles, n_chains,
                        dt, step_size_var,
                        restart, rerun_every)

fileConn<-file("parameters_used.txt")
writeLines(param_string, fileConn)
close(fileConn)
