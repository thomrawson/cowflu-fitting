##################################################################################################
orderly2::orderly_strict_mode()
orderly2::orderly_resource("Full_samples_02.rds")
orderly2::orderly_resource("Samples02_info.txt")
orderly2::orderly_parameters(n_particles = 2, total_samples = 1000)
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

params <- c("alpha", "beta", "gamma", "sigma", "asc_rate")

orderly2::orderly_artefact(description = "Plots", "Trajectory_plots/True_infected_herds.png")
orderly2::orderly_artefact(description = "Plots", "Trajectory_plots/Declared_outbreaks.png")
orderly2::orderly_artefact(description = "Plots", "Trajectory_plots/Currently_infected.png")
orderly2::orderly_artefact(description = "Plots", "Trajectory_plots/Total_infected.png")
orderly2::orderly_artefact(description = "Plots data", "Trajectory_plots/Plot_data.rds")

orderly2::orderly_artefact(description = "Counterfactual 95% CIs", "Counterfactual_differences.csv")
orderly2::orderly_artefact(description = "The parameters used", "parameters_used.txt")

##################################################################################################
library(ggplot2)
library(dplyr)
library(cowplot)
library(cowflu)
library(monty)
library(tictoc)
tic()
## Load the samples
samples <- readRDS("Full_samples_02.rds")
samples <- samples$pars
##################################################################################################
## We want to build samples of simulations for this parameter set.
burn_in_steps <- 401
total_steps <- dim(samples)[2]
samples <- samples[, burn_in_steps:total_steps, ]
## Collapse chains together:
dim(samples) <- c(5, dim(samples)[2]*dim(samples)[3])
set.seed(1)

## We can't run for all, filter to how many we want
samples <- samples[, round(seq(from = 1, to = 360000, length.out = total_samples)) ]

## Pre-allocate arrays:
Realtime_infected1 <- array(dim = c(dim(samples)[2], 51))
Confirmed_outbreaks1 <- array(dim = c(dim(samples)[2], 51))
True_infected_herds1 <- array(dim = c(dim(samples)[2], 51))
Total_infected1 <- array(dim = c(dim(samples)[2], 51))

## Prep progress bar
pb <- txtProgressBar(min = 0, max = dim(samples)[2], style = 3)

for(i in 1:dim(samples)[2]){
## Sample the model:
pars <- cowflu:::cowflu_inputs(alpha = samples[1,i], beta = samples[2,i], gamma = samples[3,i],
                      sigma = samples[4,i], asc_rate = samples[5,i],
                      dispersion = 1,
                      cowflu:::cowflu_fixed_inputs(p_region_export = cowflu:::movement$p_region_export,
                                          p_cow_export = cowflu:::movement$p_cow_export,
                                          movement_matrix = cowflu:::movement$movement_matrix,
                                          time_test = 19,
                                          n_herds_per_region = cowflu:::usda_data$n_herds_per_region,
                                          n_cows_per_herd = cowflu:::usda_data$n_cows_per_herd,
                                          start_herd = 26940, #26804 is where Texas starts.
                                          start_count = 5,
                                          condition_on_export = TRUE
                      ))


sys <- dust2::dust_system_create(cowflu:::cows(), pars, n_particles = n_particles, dt = 1)
dust2::dust_system_set_state_initial(sys)

index <- seq.int(pars$n_herds + 1, length.out = pars$n_regions)
index <- c(outer(index, (pars$n_herds + pars$n_regions) * (0:4), "+")) #*0:3 will get the SEIR, but 0:4 will also include the number of outbreaks
##Also want the "true" # of infected herds, so add a number of regions on the end again.
index <- c(index, index[length(index)] + 1:48 )

s <- dust2::dust_system_simulate(sys, 0:50, index)
s2 <- array(s, c(pars$n_regions, 6, n_particles, 51))
## Sum over dim 1, we only want national totals:
s2 <- apply(s2, c(2, 3, 4), sum)
##Fill the vectors:
Realtime_infected1[i,] <-  apply(s2[3, ,], c(2), mean)
Confirmed_outbreaks1[i,] <-  apply(s2[5, ,], c(2), mean)
True_infected_herds1[i,] <-  apply(s2[6, ,], c(2), mean)
## To calculate cumulative infections, we aggregate cumulatively, and subtract the number of recovered
Total_infected1[i,] <- cumsum(apply(s2[3, ,], c(2), mean)) - apply(s2[4, ,], c(2), mean)
## Update the progress bar:
setTxtProgressBar(pb, i)
}
## Close the progress bar
close(pb)


##Save the plot quantiles from these vectors:
Weeks <- (1:(dim(Total_infected1)[2])) - 1
Factual_data1 <- data.frame(metric = rep(NA, length(Weeks)),
                           Week = rep(NA, length(Weeks)),
                           mean_value = rep(NA, length(Weeks)),
                           lower_value = rep(NA, length(Weeks)),
                           upper_value = rep(NA, length(Weeks)))
Factual_data2 <- data.frame(metric = rep(NA, length(Weeks)),
                            Week = rep(NA, length(Weeks)),
                            mean_value = rep(NA, length(Weeks)),
                            lower_value = rep(NA, length(Weeks)),
                            upper_value = rep(NA, length(Weeks)))
Factual_data3 <- data.frame(metric = rep(NA, length(Weeks)),
                            Week = rep(NA, length(Weeks)),
                            mean_value = rep(NA, length(Weeks)),
                            lower_value = rep(NA, length(Weeks)),
                            upper_value = rep(NA, length(Weeks)))
Factual_data4 <- data.frame(metric = rep(NA, length(Weeks)),
                            Week = rep(NA, length(Weeks)),
                            mean_value = rep(NA, length(Weeks)),
                            lower_value = rep(NA, length(Weeks)),
                            upper_value = rep(NA, length(Weeks)))
## Populate the data frame:
for(i in 1:length(Weeks)){
Factual_data1$metric[i] <- "Currently_infected"
Factual_data1$Week[i] <- Weeks[i]
Factual_data1$mean_value[i] <- mean(Realtime_infected1[,i])
Factual_data1$lower_value[i] <- quantile(Realtime_infected1[,i], 0.025)
Factual_data1$upper_value[i] <- quantile(Realtime_infected1[,i], 0.975)

Factual_data2$metric[i] <- "Total_infected"
Factual_data2$Week[i] <- Weeks[i]
Factual_data2$mean_value[i] <- mean(Total_infected1[,i])
Factual_data2$lower_value[i] <- quantile(Total_infected1[,i], 0.025)
Factual_data2$upper_value[i] <- quantile(Total_infected1[,i], 0.975)

Factual_data3$metric[i] <- "Declared_outbreaks"
Factual_data3$Week[i] <- Weeks[i]
Factual_data3$mean_value[i] <- mean(Confirmed_outbreaks1[,i])
Factual_data3$lower_value[i] <- quantile(Confirmed_outbreaks1[,i], 0.025)
Factual_data3$upper_value[i] <- quantile(Confirmed_outbreaks1[,i], 0.975)

Factual_data4$metric[i] <- "True_infected_herds"
Factual_data4$Week[i] <- Weeks[i]
Factual_data4$mean_value[i] <- mean(True_infected_herds1[,i])
Factual_data4$lower_value[i] <- quantile(True_infected_herds1[,i], 0.025)
Factual_data4$upper_value[i] <- quantile(True_infected_herds1[,i], 0.975)

}
Factual_data <- rbind(Factual_data1, Factual_data2, Factual_data3, Factual_data4)

##################################################################################################
#Do it again for the counterfactual, where interventions don't happen
#Pre-allocate arrays:
Realtime_infected2 <- array(dim = c(dim(samples)[2], 51))
Confirmed_outbreaks2 <- array(dim = c(dim(samples)[2], 51))
True_infected_herds2 <- array(dim = c(dim(samples)[2], 51))
Total_infected2 <- array(dim = c(dim(samples)[2], 51))

## Prep progress bar
pb <- txtProgressBar(min = 0, max = dim(samples)[2], style = 3)

for(i in 1:dim(samples)[2]){
  ## Sample the model:
  pars <- cowflu:::cowflu_inputs(alpha = samples[1,i], beta = samples[2,i], gamma = samples[3,i],
                                 sigma = samples[4,i], asc_rate = samples[5,i],
                                 dispersion = 1,
                                 cowflu:::cowflu_fixed_inputs(p_region_export = cowflu:::movement$p_region_export,
                                                              p_cow_export = cowflu:::movement$p_cow_export,
                                                              movement_matrix = cowflu:::movement$movement_matrix,
                                                              time_test = 190000,
                                                              n_herds_per_region = cowflu:::usda_data$n_herds_per_region,
                                                              n_cows_per_herd = cowflu:::usda_data$n_cows_per_herd,
                                                              start_herd = 26940, #26804 is where Texas starts.
                                                              start_count = 5,
                                                              condition_on_export = TRUE
                                 ))


  sys <- dust2::dust_system_create(cowflu:::cows(), pars, n_particles = n_particles, dt = 1)
  dust2::dust_system_set_state_initial(sys)

  index <- seq.int(pars$n_herds + 1, length.out = pars$n_regions)
  index <- c(outer(index, (pars$n_herds + pars$n_regions) * (0:4), "+")) #*0:3 will get the SEIR, but 0:4 will also include the number of outbreaks
  ##Also want the "true" # of infected herds, so add a number of regions on the end again.
  index <- c(index, index[length(index)] + 1:48 )

  s <- dust2::dust_system_simulate(sys, 0:50, index)
  s2 <- array(s, c(pars$n_regions, 6, n_particles, 51))
  ## Sum over dim 1, we only want national totals:
  s2 <- apply(s2, c(2, 3, 4), sum)
  ##Fill the vectors:
  Realtime_infected2[i,] <-  apply(s2[3, ,], c(2), mean)
  Confirmed_outbreaks2[i,] <-  apply(s2[5, ,], c(2), mean)
  True_infected_herds2[i,] <-  apply(s2[6, ,], c(2), mean)
  ##To calculate cumulative infections, we aggregate cumulatively, and subtract the number of recovered
  Total_infected2[i,] <- cumsum(apply(s2[3, ,], c(2), mean)) - apply(s2[4, ,], c(2), mean)
  ## Update the progress bar:
  setTxtProgressBar(pb, i)
}

## Close the progress bar
close(pb)

##Save the plot quantiles from these vectors:
Weeks <- (1:(dim(Total_infected2)[2])) - 1
Counterfactual_data1 <- data.frame(metric = rep(NA, length(Weeks)),
                            Week = rep(NA, length(Weeks)),
                            mean_value = rep(NA, length(Weeks)),
                            lower_value = rep(NA, length(Weeks)),
                            upper_value = rep(NA, length(Weeks)))
Counterfactual_data2 <- data.frame(metric = rep(NA, length(Weeks)),
                            Week = rep(NA, length(Weeks)),
                            mean_value = rep(NA, length(Weeks)),
                            lower_value = rep(NA, length(Weeks)),
                            upper_value = rep(NA, length(Weeks)))
Counterfactual_data3 <- data.frame(metric = rep(NA, length(Weeks)),
                            Week = rep(NA, length(Weeks)),
                            mean_value = rep(NA, length(Weeks)),
                            lower_value = rep(NA, length(Weeks)),
                            upper_value = rep(NA, length(Weeks)))
Counterfactual_data4 <- data.frame(metric = rep(NA, length(Weeks)),
                            Week = rep(NA, length(Weeks)),
                            mean_value = rep(NA, length(Weeks)),
                            lower_value = rep(NA, length(Weeks)),
                            upper_value = rep(NA, length(Weeks)))
## Populate the data frame:
for(i in 1:length(Weeks)){
  Counterfactual_data1$metric[i] <- "Currently_infected"
  Counterfactual_data1$Week[i] <- Weeks[i]
  Counterfactual_data1$mean_value[i] <- mean(Realtime_infected2[,i])
  Counterfactual_data1$lower_value[i] <- quantile(Realtime_infected2[,i], 0.025)
  Counterfactual_data1$upper_value[i] <- quantile(Realtime_infected2[,i], 0.975)

  Counterfactual_data2$metric[i] <- "Total_infected"
  Counterfactual_data2$Week[i] <- Weeks[i]
  Counterfactual_data2$mean_value[i] <- mean(Total_infected2[,i])
  Counterfactual_data2$lower_value[i] <- quantile(Total_infected2[,i], 0.025)
  Counterfactual_data2$upper_value[i] <- quantile(Total_infected2[,i], 0.975)

  Counterfactual_data3$metric[i] <- "Declared_outbreaks"
  Counterfactual_data3$Week[i] <- Weeks[i]
  Counterfactual_data3$mean_value[i] <- mean(Confirmed_outbreaks2[,i])
  Counterfactual_data3$lower_value[i] <- quantile(Confirmed_outbreaks2[,i], 0.025)
  Counterfactual_data3$upper_value[i] <- quantile(Confirmed_outbreaks2[,i], 0.975)

  Counterfactual_data4$metric[i] <- "True_infected_herds"
  Counterfactual_data4$Week[i] <- Weeks[i]
  Counterfactual_data4$mean_value[i] <- mean(True_infected_herds2[,i])
  Counterfactual_data4$lower_value[i] <- quantile(True_infected_herds2[,i], 0.025)
  Counterfactual_data4$upper_value[i] <- quantile(True_infected_herds2[,i], 0.975)

}
Counterfactual_data <- rbind(Counterfactual_data1, Counterfactual_data2, Counterfactual_data3, Counterfactual_data4)

## Stick the two together
Factual_data$factual <- "True measures"
Counterfactual_data$factual <- "Counterfactual -\n No measures"

Plot_data <- rbind(Factual_data, Counterfactual_data)

## And what about more stringent methods? Say, test up to 100 cows, and start it week 15 (a month earlier)
## Pre-allocate arrays:
Realtime_infected3 <- array(dim = c(dim(samples)[2], 51))
Confirmed_outbreaks3 <- array(dim = c(dim(samples)[2], 51))
True_infected_herds3 <- array(dim = c(dim(samples)[2], 51))
Total_infected3 <- array(dim = c(dim(samples)[2], 51))

## Prep progress bar
pb <- txtProgressBar(min = 0, max = dim(samples)[2], style = 3)

for(i in 1:dim(samples)[2]){
  ## Sample the model:
  pars <- cowflu:::cowflu_inputs(alpha = samples[1,i], beta = samples[2,i], gamma = samples[3,i],
                                 sigma = samples[4,i], asc_rate = samples[5,i],
                                 dispersion = 1,
                                 cowflu:::cowflu_fixed_inputs(p_region_export = cowflu:::movement$p_region_export,
                                                              p_cow_export = cowflu:::movement$p_cow_export,
                                                              movement_matrix = cowflu:::movement$movement_matrix,
                                                              time_test = 15,
                                                              n_test = 100,
                                                              n_herds_per_region = cowflu:::usda_data$n_herds_per_region,
                                                              n_cows_per_herd = cowflu:::usda_data$n_cows_per_herd,
                                                              start_herd = 26940, #26804 is where Texas starts.
                                                              start_count = 5,
                                                              condition_on_export = TRUE
                                 ))


  sys <- dust2::dust_system_create(cowflu:::cows(), pars, n_particles = n_particles, dt = 1)
  dust2::dust_system_set_state_initial(sys)

  index <- seq.int(pars$n_herds + 1, length.out = pars$n_regions)
  index <- c(outer(index, (pars$n_herds + pars$n_regions) * (0:4), "+")) #*0:3 will get the SEIR, but 0:4 will also include the number of outbreaks
  ##Also want the "true" # of infected herds, so add a number of regions on the end again.
  index <- c(index, index[length(index)] + 1:48 )

  s <- dust2::dust_system_simulate(sys, 0:50, index)
  s2 <- array(s, c(pars$n_regions, 6, n_particles, 51))
  ## Sum over dim 1, we only want national totals:
  s2 <- apply(s2, c(2, 3, 4), sum)
  ##Fill the vectors:
  Realtime_infected3[i,] <-  apply(s2[3, ,], c(2), mean)
  Confirmed_outbreaks3[i,] <-  apply(s2[5, ,], c(2), mean)
  True_infected_herds3[i,] <-  apply(s2[6, ,], c(2), mean)
  ##To calculate cumulative infections, we aggregate cumulatively, and subtract the number of recovered
  Total_infected3[i,] <- cumsum(apply(s2[3, ,], c(2), mean)) - apply(s2[4, ,], c(2), mean)
  ## Update the progress bar:
  setTxtProgressBar(pb, i)
}

## Close the progress bar
close(pb)

##Save the plot quantiles from these vectors:
Weeks <- (1:(dim(Total_infected3)[2])) - 1
Counterfactual_data1 <- data.frame(metric = rep(NA, length(Weeks)),
                                   Week = rep(NA, length(Weeks)),
                                   mean_value = rep(NA, length(Weeks)),
                                   lower_value = rep(NA, length(Weeks)),
                                   upper_value = rep(NA, length(Weeks)))
Counterfactual_data2 <- data.frame(metric = rep(NA, length(Weeks)),
                                   Week = rep(NA, length(Weeks)),
                                   mean_value = rep(NA, length(Weeks)),
                                   lower_value = rep(NA, length(Weeks)),
                                   upper_value = rep(NA, length(Weeks)))
Counterfactual_data3 <- data.frame(metric = rep(NA, length(Weeks)),
                                   Week = rep(NA, length(Weeks)),
                                   mean_value = rep(NA, length(Weeks)),
                                   lower_value = rep(NA, length(Weeks)),
                                   upper_value = rep(NA, length(Weeks)))
Counterfactual_data4 <- data.frame(metric = rep(NA, length(Weeks)),
                                   Week = rep(NA, length(Weeks)),
                                   mean_value = rep(NA, length(Weeks)),
                                   lower_value = rep(NA, length(Weeks)),
                                   upper_value = rep(NA, length(Weeks)))
## Populate the data frame:
for(i in 1:length(Weeks)){
  Counterfactual_data1$metric[i] <- "Currently_infected"
  Counterfactual_data1$Week[i] <- Weeks[i]
  Counterfactual_data1$mean_value[i] <- mean(Realtime_infected3[,i])
  Counterfactual_data1$lower_value[i] <- quantile(Realtime_infected3[,i], 0.025)
  Counterfactual_data1$upper_value[i] <- quantile(Realtime_infected3[,i], 0.975)

  Counterfactual_data2$metric[i] <- "Total_infected"
  Counterfactual_data2$Week[i] <- Weeks[i]
  Counterfactual_data2$mean_value[i] <- mean(Total_infected3[,i])
  Counterfactual_data2$lower_value[i] <- quantile(Total_infected3[,i], 0.025)
  Counterfactual_data2$upper_value[i] <- quantile(Total_infected3[,i], 0.975)

  Counterfactual_data3$metric[i] <- "Declared_outbreaks"
  Counterfactual_data3$Week[i] <- Weeks[i]
  Counterfactual_data3$mean_value[i] <- mean(Confirmed_outbreaks3[,i])
  Counterfactual_data3$lower_value[i] <- quantile(Confirmed_outbreaks3[,i], 0.025)
  Counterfactual_data3$upper_value[i] <- quantile(Confirmed_outbreaks3[,i], 0.975)

  Counterfactual_data4$metric[i] <- "True_infected_herds"
  Counterfactual_data4$Week[i] <- Weeks[i]
  Counterfactual_data4$mean_value[i] <- mean(True_infected_herds3[,i])
  Counterfactual_data4$lower_value[i] <- quantile(True_infected_herds3[,i], 0.025)
  Counterfactual_data4$upper_value[i] <- quantile(True_infected_herds3[,i], 0.975)

}
Counterfactual_data <- rbind(Counterfactual_data1, Counterfactual_data2, Counterfactual_data3, Counterfactual_data4)

Counterfactual_data$factual <- "Counterfactual -\n Stronger measures"

Plot_data <- rbind(Plot_data, Counterfactual_data)

##################################################################################################
##################################################################################################
## Plotting
dir.create("Trajectory_plots")
## Save the plot data
saveRDS(Plot_data, "Trajectory_plots/Plot_data.rds")

## Total infected
ggplot(filter(Plot_data, metric == "Total_infected")) +
  geom_line(aes(x = Week, y = mean_value, colour = factual)) +
  geom_ribbon(aes(x= Week, ymin = lower_value, ymax = upper_value, fill = factual),
              alpha = 0.2) +
  theme_minimal() +
  labs(x = "Time (Weeks)", y = "Infected Dairy Cows", title = "Total cumulative infected dairy cattle",
       color = "Scenario", fill = "Scenario") +
  theme(strip.text = element_text(size = 8)) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 20),    # Title font size
    axis.title.x = element_text(size = 16),  # X-axis title font size
    axis.title.y = element_text(size = 16),  # Y-axis title font size
    axis.text.x = element_text(size = 12),   # X-axis labels font size
    axis.text.y = element_text(size = 12),   # Y-axis labels font size
    legend.text = element_text(size = 14),   # Legend text font size
    legend.title = element_text(size = 16),  # Legend title font size
    legend.key.spacing.y = unit(0.5, 'cm')
  )  -> p
ggsave(filename = "Trajectory_plots/Total_infected.png",
       plot = p, bg = "white", width = 1280/96, height = 720/96, dpi = 200)


## Currently infected
ggplot(filter(Plot_data, metric == "Currently_infected")) +
  geom_line(aes(x = Week, y = mean_value, colour = factual)) +
  geom_ribbon(aes(x= Week, ymin = lower_value, ymax = upper_value, fill = factual),
              alpha = 0.2) +
  theme_minimal() +
  labs(x = "Time (Weeks)", y = "Infected Dairy Cows", title = "Weekly infected dairy cattle") +
  theme(strip.text = element_text(size = 8)) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 20),    # Title font size
    axis.title.x = element_text(size = 16),  # X-axis title font size
    axis.title.y = element_text(size = 16),  # Y-axis title font size
    axis.text.x = element_text(size = 12),   # X-axis labels font size
    axis.text.y = element_text(size = 12),   # Y-axis labels font size
    legend.text = element_text(size = 14),   # Legend text font size
    legend.title = element_text(size = 16),  # Legend title font size
    legend.key.spacing.y = unit(0.5, 'cm')
  )  -> p
ggsave(filename = "Trajectory_plots/Currently_infected.png",
       plot = p, bg = "white", width = 1280/96, height = 720/96, dpi = 200)

## Declared outbreaks
ggplot(filter(Plot_data, metric == "Declared_outbreaks")) +
  geom_line(aes(x = Week, y = mean_value, colour = factual)) +
  geom_ribbon(aes(x= Week, ymin = lower_value, ymax = upper_value, fill = factual),
              alpha = 0.2) +
  theme_minimal() +
  labs(x = "Time (Weeks)", y = "Declared outbreaks", title = "Total declared outbreaks") +
  theme(strip.text = element_text(size = 8)) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 20),    # Title font size
    axis.title.x = element_text(size = 16),  # X-axis title font size
    axis.title.y = element_text(size = 16),  # Y-axis title font size
    axis.text.x = element_text(size = 12),   # X-axis labels font size
    axis.text.y = element_text(size = 12),   # Y-axis labels font size
    legend.text = element_text(size = 14),   # Legend text font size
    legend.title = element_text(size = 16),  # Legend title font size
    legend.key.spacing.y = unit(0.5, 'cm')
  )  -> p
ggsave(filename = "Trajectory_plots/Declared_outbreaks.png",
       plot = p, bg = "white", width = 1280/96, height = 720/96, dpi = 200)

## True infected herds
ggplot(filter(Plot_data, metric == "True_infected_herds")) +
  geom_line(aes(x = Week, y = mean_value, colour = factual)) +
  geom_ribbon(aes(x= Week, ymin = lower_value, ymax = upper_value, fill = factual),
              alpha = 0.2) +
  theme_minimal() +
  labs(x = "Time (Weeks)", y = "Infected Dairy Herds", title = "Total infected herds") +
  theme(strip.text = element_text(size = 8)) +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(size = 20),    # Title font size
    axis.title.x = element_text(size = 16),  # X-axis title font size
    axis.title.y = element_text(size = 16),  # Y-axis title font size
    axis.text.x = element_text(size = 12),   # X-axis labels font size
    axis.text.y = element_text(size = 12),   # Y-axis labels font size
    legend.text = element_text(size = 14),   # Legend text font size
    legend.title = element_text(size = 16),  # Legend title font size
    legend.key.spacing.y = unit(0.5, 'cm')
  )  -> p
ggsave(filename = "Trajectory_plots/True_infected_herds.png",
       plot = p, bg = "white", width = 1280/96, height = 720/96, dpi = 200)

##################################################################################################
##################################################################################################
## Calculate the 95% CIs for the differences of the counterfactuals
No_measures_Realtime_infected <- Realtime_infected2 - Realtime_infected1
No_measures_Confirmed_outbreaks <- Confirmed_outbreaks2 - Confirmed_outbreaks1
No_measures_True_infected_herds <- True_infected_herds2 - True_infected_herds1
No_measures_Total_infected <- Total_infected2 - Total_infected1

Stronger_measures_Realtime_infected <- Realtime_infected3 - Realtime_infected1
Stronger_measures_Confirmed_outbreaks <- Confirmed_outbreaks3 - Confirmed_outbreaks1
Stronger_measures_True_infected_herds <- True_infected_herds3 - True_infected_herds1
Stronger_measures_Total_infected <- Total_infected3 - Total_infected1

## Export in a data frame
Counterfactual_differences <- data.frame( measures = rep(c("No measures", "Stronger measures"),4),
                                          metric = c(rep("Realtime_infected",2),
                                                     rep("Confirmed_outbreaks",2),
                                                     rep("Declared_infected_herds",2),
                                                     rep("Total_infected",2)),
                                          mean_value = c(mean(No_measures_Realtime_infected),
                                                         mean(Stronger_measures_Realtime_infected),
                                                         mean(No_measures_Confirmed_outbreaks),
                                                         mean(Stronger_measures_Confirmed_outbreaks),
                                                         mean(No_measures_True_infected_herds),
                                                         mean(Stronger_measures_True_infected_herds),
                                                         mean(No_measures_Total_infected),
                                                         mean(Stronger_measures_Total_infected)),
                                          lower_value = c(quantile(No_measures_Realtime_infected, 0.025),
                                                          quantile(Stronger_measures_Realtime_infected, 0.025),
                                                          quantile(No_measures_Confirmed_outbreaks, 0.025),
                                                          quantile(Stronger_measures_Confirmed_outbreaks, 0.025),
                                                          quantile(No_measures_True_infected_herds, 0.025),
                                                          quantile(Stronger_measures_True_infected_herds, 0.025),
                                                          quantile(No_measures_Total_infected, 0.025),
                                                          quantile(Stronger_measures_Total_infected, 0.025)),
                                          upper_value = c(quantile(No_measures_Realtime_infected, 0.975),
                                                          quantile(Stronger_measures_Realtime_infected, 0.975),
                                                          quantile(No_measures_Confirmed_outbreaks, 0.975),
                                                          quantile(Stronger_measures_Confirmed_outbreaks, 0.975),
                                                          quantile(No_measures_True_infected_herds, 0.975),
                                                          quantile(Stronger_measures_True_infected_herds, 0.975),
                                                          quantile(No_measures_Total_infected, 0.975),
                                                          quantile(Stronger_measures_Total_infected, 0.975))
                                          )
write.csv(Counterfactual_differences, "Counterfactual_differences.csv", row.names = FALSE)

duration <- toc()


## Print an output .txt of the parameters used
param_string <- sprintf("duration ran: %s mins\n
  total_samples: %s \n
  n_particles: %s \n
  ",   (duration$toc - duration$tic)/60,
                        total_samples,
                        n_particles)

fileConn<-file("parameters_used.txt")
writeLines(param_string, fileConn)
close(fileConn)
