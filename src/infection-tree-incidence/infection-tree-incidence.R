##################################################################################################
orderly2::orderly_strict_mode()
orderly2::orderly_resource("Full_samples_03_incidence.rds")
orderly2::orderly_resource("Samples03_info.txt")
# orderly2::orderly_parameters(n_samples = 100, n_particles = 8, n_chains = 2, step_size_var = 0.03, dt = 1,
#                              restart = FALSE, rerun_every = 100)
#
# orderly2::orderly_dependency("fitting-survival",
#                              "latest(parameter:n_samples == this:n_samples &&
#                                     parameter:n_particles == this:n_particles &&
#                                     parameter:n_chains == this:n_chains &&
#                                     parameter:step_size_var == this:step_size_var &&
#                                     parameter:dt == this:dt &&
#                                     parameter:restart == this:restart &&
#                                     parameter:rerun_every == this:rerun_every)",
#                              "fitting_samples.rds")

params <- c("alpha", "beta", "gamma", "sigma", "asc_rate", "dispersion")

orderly2::orderly_artefact(description = "MCMC diagnostics", "mcmc_diagnostics.txt")

orderly2::orderly_artefact(description = "Fits to data", "Fit_to_data/Infected_Herds_1.png")
orderly2::orderly_artefact(description = "Fits to data", "Fit_to_data/Infected_Proportion_Herds_1.png")
orderly2::orderly_artefact(description = "Fits to data", "Fit_to_data/Infected_Herds_2.png")
orderly2::orderly_artefact(description = "Fits to data", "Fit_to_data/Infected_Proportion_Herds_2.png")
orderly2::orderly_artefact(description = "Fits to data", "Fit_to_data/Declared_Infected_Herds_1.png")
orderly2::orderly_artefact(description = "Fits to data", "Fit_to_data/Declared_Infected_Herds_2.png")

orderly2::orderly_artefact(description = "Prior/Posterior comparison", "Prior_posterior/Alpha.png")
orderly2::orderly_artefact(description = "Prior/Posterior comparison", "Prior_posterior/Beta.png")
orderly2::orderly_artefact(description = "Prior/Posterior comparison", "Prior_posterior/Gamma.png")
orderly2::orderly_artefact(description = "Prior/Posterior comparison", "Prior_posterior/Sigma.png")
orderly2::orderly_artefact(description = "Prior/Posterior comparison", "Prior_posterior/asc_rate.png")
orderly2::orderly_artefact(description = "Prior/Posterior comparison", "Prior_posterior/Dispersion.png")

orderly2::orderly_artefact(description = "Histogram of generation time", "Generation_time.png")

for(j in 1:sum(cowflu:::outbreaks_data$weekly_outbreaks_data$Positive_Tests > 0) ){
  orderly2::orderly_artefact(description = "Infection trees", sprintf("Infection_trees/state_weeks/Outbreak_%s.png", j))
  orderly2::orderly_artefact(description = "Infection trees", sprintf("Infection_trees/state/Outbreak_%s.png", j))
}
##################################################################################################
library(ggplot2)
library(dplyr)
library(cowplot)
library(maps)
library(ggrepel)
library(sf)
library(stringr)
library(coda)
library(monty)
## Load the samples
samples <- readRDS("Full_samples_03_incidence.rds")
##################################################################################################
## MCMC diagnostics:
ess <- coda::effectiveSize(coda::as.mcmc.list(samples))
results <- data.frame(ESS = ess)
## Calculate PSRF (Potential Scale Reduction Factor)
## This can fail for very small trial runs:
tryCatch({
  # Code that may throw an error
  psrf <- coda::gelman.diag(coda::as.mcmc.list(samples))$psrf
  results <- data.frame(ESS = ess, PSRF = psrf[,1], PSRF_95CI = psrf[,2])
}, error = function(e) {
  # Code to run if there's an error
  cat("An error occurred, the leading minor of order 3 is not positive.\n")
  # Optional: You can log the error or take alternative actions here
})
## Write the results to a .txt file
write.table(results, file = "mcmc_diagnostics.txt", sep = "\t", col.names = TRUE, row.names = TRUE)
##################################################################################################
##################################################################################################
## Fit to data plots:
total_steps <- dim(samples$observations$history)[3]
burnin_steps <- 401

## Load data:
data_outbreaks <- cowflu:::process_data_incidence(cowflu:::outbreaks_data$weekly_outbreaks_data)
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

data_week <- dust2::dust_filter_data(data_outbreaks, time = "week")
## Here we plot the assoc. trajectories too.
## Fit to data plotting
dir.create("Fit_to_data")
plot_indices <- list(first = c(1,25), second = c(26,48))

## First, the TRUE number of infected herds per state:
for(i in 1:length(plot_indices)){
  regions_to_plot <- plot_indices[[i]]
  regions_to_plot <- regions_to_plot[1]:regions_to_plot[2]
  n_states <- length(regions_to_plot)

  ## Extract the  data across all chains and iterations
  I_data <- samples$observations$history[49:96,,burnin_steps:total_steps,]
  n_timepoints <- dim(I_data)[2]
  n_samples <- dim(I_data)[3]
  n_chains <- dim(I_data)[4]
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
    ylim(c(0,1.01)) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) -> my_plot_prop

  file_path <- sprintf("Fit_to_data/Infected_Herds_%s.png", i)
  ggsave(filename = file_path, plot = my_plot, width = 1280/96, height = 720/96, dpi = 200)
  file_path <- sprintf("Fit_to_data/Infected_Proportion_Herds_%s.png", i)
  ggsave(filename = file_path, plot = my_plot_prop, width = 1280/96, height = 720/96, dpi = 200)

}
##################################################################################################
## And now with the declared outbreaks that we actually fit to:
## We will need to reorder the data into a data frame:
real_outbreaks_data <- data_week %>%
  tidyr::unnest_longer(positive_tests) %>%
  mutate(state = rep(cowflu:::usda_data$US_States, times = nrow(data_week)))
real_outbreaks_data <- real_outbreaks_data[,c(2,3,4)]
colnames(real_outbreaks_data) <- c("Time", "outbreak_detected", "US_state")

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
  ## Replace the US_state numbers with actual names
  result_df$US_state <- factor(result_df$US_state, levels = 1:n_states, labels = cowflu:::usda_data$US_States[regions_to_plot])

  ## Filter the data we fit to:
  specific_plot_outbreaks_data <- filter(real_outbreaks_data, US_state %in% cowflu:::usda_data$US_States[regions_to_plot])

  ## Create the plot
  ggplot(result_df, aes(x = Time, y = mean_infected, group = US_state)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower_ci_infected, ymax = upper_ci_infected), alpha = 0.2) +
    geom_point(data = specific_plot_outbreaks_data,
               aes(x = Time, y = outbreak_detected, group = US_state), col = "red") +
    facet_wrap(~ US_state, ncol = 5, scales = "free_y") +
    theme_minimal() +
    labs(x = "Time (Weeks)", y = "Declared Infected Herds", title = "Declared Infected Herds per state") +
    theme(strip.text = element_text(size = 8)) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    ) -> my_plot

  file_path <- sprintf("Fit_to_data/Declared_Infected_Herds_%s.png", i)
  ggsave(filename = file_path, plot = my_plot, width = 1280/96, height = 720/96, dpi = 200)
}
##################################################################################################
## Tidy up:
rm(state_summary, ess, file_path, i, I_data, n_states, n_timepoints,
   regions_to_plot, state, data_list, data_outbreaks, my_plot, my_plot_prop, plot_indices, psrf,
   real_outbreaks_data, result_df, results, specific_plot_outbreaks_data, state_data, state_df)
##################################################################################################
##################################################################################################
## Plot prior and posterior distributions for each parameter:
dir.create("Prior_posterior")
## ALPHA

parameter_data <- data.frame( value = rbeta(100000, shape1 = 1, shape2 = 40), distr = rep("prior", 100000))
parameter_data <- rbind(parameter_data, data.frame( value = as.numeric(samples$pars[1, burnin_steps:total_steps,]),
                                                    distr = rep("posterior", (total_steps - burnin_steps + 1)*n_chains)))
ggplot(parameter_data) +
  geom_density(aes(x = value, fill = distr, color = distr), alpha = 0.5) +
  theme_classic() +
  labs(
    title = "Alpha Prior - Beta(1,40)",
    x = "Alpha value",
    y = "Density",
    color = "Distribution",
    fill = "Distribution") -> plot

file_path <- "Prior_posterior/Alpha.png"
ggsave(filename = file_path, plot = plot, width = 1280/192, height = 720/192, dpi = 200)

##Beta
parameter_data <- data.frame( value = runif(100000, 0.05, 3), distr = rep("prior", 100000))
parameter_data <- rbind(parameter_data, data.frame( value = as.numeric(samples$pars[2, burnin_steps:total_steps,]),
                                                    distr = rep("posterior", (total_steps - burnin_steps + 1)*n_chains)))
ggplot(parameter_data) +
  geom_density(aes(x = value, fill = distr, color = distr), alpha = 0.5) +
  theme_classic() +
  labs(
    title = "Beta Prior - Uniform(0.05,3)",
    x = "Beta value",
    y = "Density",
    color = "Distribution",
    fill = "Distribution") -> plot

file_path <- "Prior_posterior/Beta.png"
ggsave(filename = file_path, plot = plot, width = 1280/192, height = 720/192, dpi = 200)

## Gamma
parameter_data <- data.frame( value = runif(100000, 0.05, 3), distr = rep("prior", 100000))
parameter_data <- rbind(parameter_data, data.frame( value = as.numeric(samples$pars[3, burnin_steps:total_steps,]),
                                                    distr = rep("posterior", (total_steps - burnin_steps + 1)*n_chains)))
ggplot(parameter_data) +
  geom_density(aes(x = value, fill = distr, color = distr), alpha = 0.5) +
  theme_classic() +
  labs(
    title = "Gamma Prior - Uniform(0.05,3)",
    x = "Gamma value",
    y = "Density",
    color = "Distribution",
    fill = "Distribution") -> plot

file_path <- "Prior_posterior/Gamma.png"
ggsave(filename = file_path, plot = plot, width = 1280/192, height = 720/192, dpi = 200)

## Sigma
parameter_data <- data.frame( value = runif(100000, 0.05, 5), distr = rep("prior", 100000))
parameter_data <- rbind(parameter_data, data.frame( value = as.numeric(samples$pars[4, burnin_steps:total_steps,]),
                                                    distr = rep("posterior", (total_steps - burnin_steps + 1)*n_chains)))
ggplot(parameter_data) +
  geom_density(aes(x = value, fill = distr, color = distr), alpha = 0.5) +
  theme_classic() +
  labs(
    title = "Sigma Prior - Uniform(0.05,5)",
    x = "Sigma value",
    y = "Density",
    color = "Distribution",
    fill = "Distribution") -> plot

file_path <- "Prior_posterior/Sigma.png"
ggsave(filename = file_path, plot = plot, width = 1280/192, height = 720/192, dpi = 200)

## Asc_rate
parameter_data <- data.frame( value = rbeta(100000, shape1 = 5, shape2 = 1), distr = rep("prior", 100000))
parameter_data <- rbind(parameter_data, data.frame( value = as.numeric(samples$pars[5, burnin_steps:total_steps,]),
                                                    distr = rep("posterior", (total_steps - burnin_steps + 1)*n_chains)))
ggplot(parameter_data) +
  geom_density(aes(x = value, fill = distr, color = distr), alpha = 0.5) +
  theme_classic() +
  labs(
    title = "Ascertainment Rate Prior - Beta(5,1)",
    x = "asc_rate value",
    y = "Density",
    color = "Distribution",
    fill = "Distribution") -> plot

file_path <- "Prior_posterior/asc_rate.png"
ggsave(filename = file_path, plot = plot, width = 1280/192, height = 720/192, dpi = 200)

## Dispersion
parameter_data <- data.frame( value = rexp(100000, rate = 1/2), distr = rep("prior", 100000))
parameter_data <- rbind(parameter_data, data.frame( value = as.numeric(samples$pars[6, burnin_steps:total_steps,]),
                                                    distr = rep("posterior", (total_steps - burnin_steps + 1)*n_chains)))
ggplot(parameter_data) +
  geom_density(aes(x = value, fill = distr, color = distr), alpha = 0.5) +
  theme_classic() +
  labs(
    title = "Dispersion Prior - Exp(mean = 2)",
    x = "dispersion value",
    y = "Density",
    color = "Distribution",
    fill = "Distribution") -> plot

file_path <- "Prior_posterior/Dispersion.png"
ggsave(filename = file_path, plot = plot, width = 1280/192, height = 720/192, dpi = 200)

##############################################################################################
##############################################################################################
## Tidy up:
rm(parameter_data, plot, file_path)
##############################################################################################
##############################################################################################
## Now, for the actual infection trees.
## For each outbreak, we'll calculate the probability of source for each state/week
outbreaks_data <- cowflu:::outbreaks_data$weekly_outbreaks_data
outbreaks_data <- filter(outbreaks_data, Positive_Tests > 0)
start_date <- as.Date("2023-12-18")
outbreaks_data$day <- as.numeric(outbreaks_data$Week_Beginning - start_date)
outbreaks_data$Week <- outbreaks_data$day/7

n_extra_samples <- 100
## Build a kde density for the generation time, from sigma and gamma;
N_gen_time_samples <- (dim(samples$pars)[2] - burnin_steps + 1) * dim(samples$pars)[3] * n_extra_samples
gen_time_samples <- rep(NA, N_gen_time_samples)

for(i in 1:(dim(samples$pars)[2] - burnin_steps + 1)){
  for(j in 1:dim(samples$pars)[3]){
    T_E <- rexp(100, samples$pars["sigma",i + burnin_steps - 1, j])
    T_I <- rexp(100, samples$pars["gamma",i + burnin_steps - 1, j])
    gen_time_samples[((i-1)*dim(samples$pars)[3] + j - 1)*n_extra_samples + 1:n_extra_samples] <- T_E + T_I
  }
}

png("Generation_time.png", width = 800, height = 600)  # Opens the PNG device with specified dimensions
par(cex = 1.5)
hist(gen_time_samples, breaks = 800, xlim = c(0, 10), main = "Generation time distribution", xlab = "Generation time (weeks)")
dev.off()
##Tidy up
gen_time_kde <- density(gen_time_samples)
rm(N_gen_time_samples, gen_time_samples)

##Produce a plot for each outbreak:
dir.create("Infection_trees")
dir.create("Infection_trees/state_weeks")
dir.create("Infection_trees/state")
for(j in 1:length(outbreaks_data$State)){

  Outbreak <- outbreaks_data[j,]
  Weeks_to_consider <- Outbreak$Week - 1:8
  Probabilities_of_Source <- data.frame(State = rep(tolower(cowflu:::usda_data$US_States), length(Weeks_to_consider )),
                                        Week = sort(rep(Weeks_to_consider, length(cowflu:::usda_data$US_States))),
                                        Probability_of_source = rep(NA, length(cowflu:::usda_data$US_States)*length(Weeks_to_consider)) ,
                                        Currently_Infected_Herds = rep(NA, length(cowflu:::usda_data$US_States)*length(Weeks_to_consider)) )
  ## Populate the "currently infected" column:
  for(i in 1:length(Probabilities_of_Source$State)){
    n_state <- 48 + which(tolower(cowflu:::usda_data$US_States) == Probabilities_of_Source$State[i])
    n_week <- Probabilities_of_Source$Week[i]
    Probabilities_of_Source$Currently_Infected_Herds[i] <- mean(samples$observations$history[n_state, n_week,201:4000,])
    Probabilities_of_Source$prop_infected[i] <- Probabilities_of_Source$Currently_Infected_Herds[i]/cowflu:::usda_data$n_herds_per_region[n_state-48]
  }

  ## Populate the probabilities:
  for(i in 1:length(Probabilities_of_Source$State)){
    gen_time_prob <- approx(gen_time_kde$x, gen_time_kde$y, Outbreak$Week - Probabilities_of_Source$Week[i])$y
    #Probability that the source exports, multiplied by prob it's to the outbreak destination state
    # Note this currently doesn't capture the question about whether or not an infected cow actually got through the border
    # Technically coming from within SHOULD be more likely...
    export_prob <- cowflu:::movement$p_region_export[which(tolower(cowflu:::usda_data$US_States) == Probabilities_of_Source$State[i])] * (Probabilities_of_Source$Currently_Infected_Herds[i]/cowflu:::usda_data$n_herds_per_region[which(tolower(cowflu:::usda_data$US_States) == Probabilities_of_Source$State[i])])
    Probabilities_of_Source$Probability_of_source[i] <- gen_time_prob * export_prob
  }

  ## Scale Probability_of_source column so it sums to 1:
  Probabilities_of_Source$Probability_of_source <- Probabilities_of_Source$Probability_of_source/sum(Probabilities_of_Source$Probability_of_source)

  ## Order the data frame so highest probabilities are at the top:
  Probabilities_of_Source <- Probabilities_of_Source[order(Probabilities_of_Source$Probability_of_source, decreasing = TRUE),]

  ##Now we plot the top three probabilities:
  # 1. Get map data for US without Alaska and Hawaii
  us_states <- map_data("state")

  # 2. Prepare your outbreak data (assuming single state for outbreak)
  outbreak_state <- Outbreak$State[1]
  outbreak_week <- Outbreak$Week[1]

  # Get the coordinates of the outbreak state
  outbreak_coords <- us_states %>%
    filter(region == outbreak_state) %>%
    summarize(lat = mean(lat), long = mean(long))
  ##Jitter them...
  outbreak_coords$lat <- outbreak_coords$lat + 0.001
  outbreak_coords$long <- outbreak_coords$long + 0.001

  # 3. Find top 3 most probable source states for the outbreak week
  top_sources <- Probabilities_of_Source %>%
    slice(1:3)

  # 4. Get the coordinates of the top source states
  source_coords <- us_states %>%
    filter(region %in% top_sources$State) %>%
    group_by(region) %>%
    summarize(lat = mean(lat), long = mean(long))
  # Merge with probabilities
  source_coords <- source_coords %>%
    left_join(top_sources, by = c("region" = "State"))
  # Rank the sources by probability and add a curvature based on rank
  source_coords <- source_coords %>%
    mutate(rank = row_number(),  # Rank the sources by their row number (1, 2, 3)
           curvature = case_when(
             rank == 1 ~ 0.2,
             rank == 2 ~ -0.2,
             rank == 3 ~ 0.4
           ))
  # Extract the prop_infected for the most probable week
  infected_to_fill <- filter(Probabilities_of_Source, Week == top_sources$Week[1])
  us_states <- us_states %>%
    left_join(infected_to_fill[,c("State", "prop_infected")], by = c("region" = "State"))

  # 5. Base map
  p <- ggplot() +
    geom_polygon(data = us_states, aes(x = long, y = lat, group = group,
                                       fill = prop_infected),
                 color = "black") +
    coord_fixed(1.3) +  # Ensure map has the correct proportions
    theme_minimal() +
    scale_fill_gradient(low = "white", high = "red", limits = c(0,1))

  # 6. Mark the outbreak state with a star
  p <- p +
    geom_point(data = outbreak_coords, aes(x = long, y = lat),
               shape = 8, size = 5, color = "red") +  # Star for outbreak state
    geom_text_repel(data = outbreak_coords, aes(x = long, y = lat, label = str_to_title(outbreak_state)),
                    nudge_y = 1, size = 4)


  # 7. Draw arrows from sources to outbreak, with color intensity by probability
  # For rank 1 (first most probable source)
  p <- p +
    geom_curve(data = source_coords %>% filter(rank == 1),
               aes(x = long, y = lat, xend = outbreak_coords$long, yend = outbreak_coords$lat,
                   color = Probability_of_source),
               arrow = arrow(length = unit(0.3, "cm")),
               curvature = 0.9,  # Curve for first rank
               size = 1.5)

  # For rank 2 (second most probable source)
  p <- p +
    geom_curve(data = source_coords %>% filter(rank == 2),
               aes(x = long, y = lat, xend = outbreak_coords$long, yend = outbreak_coords$lat,
                   color = Probability_of_source),
               arrow = arrow(length = unit(0.3, "cm")),
               curvature = -0.2,  # Curve for second rank
               size = 1.5)

  # For rank 3 (third most probable source)
  p <- p +
    geom_curve(data = source_coords %>% filter(rank == 3),
               aes(x = long, y = lat, xend = outbreak_coords$long, yend = outbreak_coords$lat,
                   color = Probability_of_source),
               arrow = arrow(length = unit(0.3, "cm")),
               curvature = 0.4,  # Larger curve for third rank
               size = 1.5)
  # 8. Add points and labels for the source states
  p <- p +
    geom_point(data = unique(source_coords[,c("region", "lat", "long")]), aes(x = long, y = lat), color = "blue", size = 3) +
    geom_text_repel(data = unique(source_coords[,c("region", "lat", "long")]), aes(x = long, y = lat, label = str_to_title(region)), size = 3)

  # 9. Customize color gradient and labels
  p <- p +
    theme_void() +
    scale_color_gradient(low = "lightblue", high = "darkblue") +
    labs(title = sprintf("    %s - Week %s (%s)", str_to_title(Outbreak$State), Outbreak$Week, Outbreak$Week_Beginning),
         color = "Probability of Source",
         fill = "Proportion infected   ") +
    theme(legend.position = "bottom",
          legend.box = "vertical") +
    guides(color = guide_colourbar(order = 1),  # Color legend on top (order 1)
           fill = guide_colourbar(order = 2))

  # 10. Save the plot
  file_path <- sprintf("Infection_trees/state_weeks/Outbreak_%s.png", j)
  ggsave(filename = file_path, plot = p, bg = "white", width = 1280/96, height = 720/96, dpi = 200)

  ## Repeat this plot, but we want to sum together all weeks for states, to see the three most plausible states

  ## Sum the Probabilities of source data_frame over week:
  Probabilities_of_Source <- Probabilities_of_Source %>%
    group_by(State) %>%
    summarize(
      # Sum Probability_of_source over all weeks for each state
      Probability_of_source = sum(Probability_of_source),

      # Get the row with the highest Week for Currently_Infected_Herds and prop_infected
      Currently_Infected_Herds = Currently_Infected_Herds[which.max(Week)],
      prop_infected = prop_infected[which.max(Week)]
    )

  ## Order the data frame so highest probabilities are at the top:
  Probabilities_of_Source <- Probabilities_of_Source[order(Probabilities_of_Source$Probability_of_source, decreasing = TRUE),]

  ##Now we plot the top three probabilities:
  # 1. Get map data for US without Alaska and Hawaii
  us_states <- map_data("state")

  # 2. Prepare your outbreak data (assuming single state for outbreak)
  outbreak_state <- Outbreak$State[1]
  outbreak_week <- Outbreak$Week[1]

  # Get the coordinates of the outbreak state
  outbreak_coords <- us_states %>%
    filter(region == outbreak_state) %>%
    summarize(lat = mean(lat), long = mean(long))
  ##Jitter them...
  outbreak_coords$lat <- outbreak_coords$lat + 0.001
  outbreak_coords$long <- outbreak_coords$long + 0.001

  # 3. Find top 3 most probable source states for the outbreak week
  top_sources <- Probabilities_of_Source %>%
    slice(1:3)

  # 4. Get the coordinates of the top source states
  source_coords <- us_states %>%
    filter(region %in% top_sources$State) %>%
    group_by(region) %>%
    summarize(lat = mean(lat), long = mean(long))
  # Merge with probabilities
  source_coords <- source_coords %>%
    left_join(top_sources, by = c("region" = "State"))
  # Rank the sources by probability and add a curvature based on rank
  source_coords <- source_coords %>%
    mutate(rank = row_number(),  # Rank the sources by their row number (1, 2, 3)
           curvature = case_when(
             rank == 1 ~ 0.2,
             rank == 2 ~ -0.2,
             rank == 3 ~ 0.4
           ))
  # Extract the prop_infected for the last week
  infected_to_fill <- Probabilities_of_Source
  us_states <- us_states %>%
    left_join(infected_to_fill[,c("State", "prop_infected")], by = c("region" = "State"))

  # 5. Base map
  p <- ggplot() +
    geom_polygon(data = us_states, aes(x = long, y = lat, group = group,
                                       fill = prop_infected),
                 color = "black") +
    coord_fixed(1.3) +  # Ensure map has the correct proportions
    theme_minimal() +
    scale_fill_gradient(low = "white", high = "red", limits = c(0,1))

  # 6. Mark the outbreak state with a star
  p <- p +
    geom_point(data = outbreak_coords, aes(x = long, y = lat),
               shape = 8, size = 5, color = "red") +  # Star for outbreak state
    geom_text_repel(data = outbreak_coords, aes(x = long, y = lat, label = str_to_title(outbreak_state)),
                    nudge_y = 1, size = 4)


  # 7. Draw arrows from sources to outbreak, with color intensity by probability
  # For rank 1 (first most probable source)
  p <- p +
    geom_curve(data = source_coords %>% filter(rank == 1),
               aes(x = long, y = lat, xend = outbreak_coords$long, yend = outbreak_coords$lat,
                   color = Probability_of_source),
               arrow = arrow(length = unit(0.3, "cm")),
               curvature = 0.9,  # Curve for first rank
               size = 1.5)

  # For rank 2 (second most probable source)
  p <- p +
    geom_curve(data = source_coords %>% filter(rank == 2),
               aes(x = long, y = lat, xend = outbreak_coords$long, yend = outbreak_coords$lat,
                   color = Probability_of_source),
               arrow = arrow(length = unit(0.3, "cm")),
               curvature = -0.2,  # Curve for second rank
               size = 1.5)

  # For rank 3 (third most probable source)
  p <- p +
    geom_curve(data = source_coords %>% filter(rank == 3),
               aes(x = long, y = lat, xend = outbreak_coords$long, yend = outbreak_coords$lat,
                   color = Probability_of_source),
               arrow = arrow(length = unit(0.3, "cm")),
               curvature = 0.4,  # Larger curve for third rank
               size = 1.5)
  # 8. Add points and labels for the source states
  p <- p +
    geom_point(data = unique(source_coords[,c("region", "lat", "long")]), aes(x = long, y = lat), color = "blue", size = 3) +
    geom_text_repel(data = unique(source_coords[,c("region", "lat", "long")]), aes(x = long, y = lat, label = str_to_title(region)), size = 3)

  # 9. Customize color gradient and labels
  p <- p +
    theme_void() +
    scale_color_gradient(low = "lightblue", high = "darkblue") +
    labs(title = sprintf("    %s - Week %s (%s)", str_to_title(Outbreak$State), Outbreak$Week, Outbreak$Week_Beginning),
         color = "Probability of Source",
         fill = "Proportion infected   ") +
    theme(legend.position = "bottom",
          legend.box = "vertical") +
    guides(color = guide_colourbar(order = 1),  # Color legend on top (order 1)
           fill = guide_colourbar(order = 2))

  # 10. Save the plot
  file_path <- sprintf("Infection_trees/state/Outbreak_%s.png", j)
  ggsave(filename = file_path, plot = p, bg = "white", width = 1280/96, height = 720/96, dpi = 200)


}
