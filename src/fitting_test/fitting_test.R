##################################################################################################
orderly2::orderly_strict_mode()
orderly2::orderly_parameters(n_samples = 100, n_particles = 8, n_chains = 1, step_size_var = 0.02, dt = 1)

params <- c("alpha", "beta", "gamma", "sigma", "asc_rate")
orderly2::orderly_artefact(description = "The posterior samples", "fitting_samples.rds")
orderly2::orderly_artefact(description = "Log_likelihood of chains", "ll_plot.png")
orderly2::orderly_artefact(description = "The posterior density of all chains", "Posterior_densities//00_all_plot.png")
orderly2::orderly_artefact(description = "The posterior sampling trajectories of all chains", "Posterior_chains//00_all_plot.png")
orderly2::orderly_artefact(description = "The parameters used", "parameters_used.txt")
for(i in 1:length(params)){
  orderly2::orderly_artefact(description = "Posterior density", sprintf("Posterior_densities//%s_plot.png", params[i]))
  orderly2::orderly_artefact(description = "Posterior trajectories", sprintf("Posterior_chains//%s_plot.png", params[i]))
}
##################################################################################################
library(cowflu)
library(dust2)
library(monty)
library(ggplot2)
library(gridExtra)
library(tictoc)
tic()
##################################################################################################
## Set parameters
pars <- cowflu:::cowflu_inputs(
  alpha = 0.2,
  beta = 1.5,
  gamma = 0.1,
  sigma = 0.125,
  asc_rate = 1,
  dispersion = 1,
  cowflu:::cowflu_fixed_inputs(
    n_herds_per_region = c(3, 7, 11),
    p_region_export = c(.5, .5, .5),
    p_cow_export = c(0.2, 0.2, 0.2),
    n_cows_per_herd = c(rep(200, 3), rep(1000, 7), rep(3000, 11)),
    movement_matrix = cbind(c(.6, .2, .2), c(.2, .6, .2), c(.2, .2, .6)),
    time_test = 1,
    start_herd = 4,
    start_count = 5,
    likelihood_choice = "survival"))

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
  alpha ~ Beta(a = 1, b = 25)
  beta ~ Uniform(min = 1, max = 2.5) #maybe 1 and 2.5
  gamma ~ Uniform(min = 0.05, max = 1) #0 and 1
  sigma ~ Uniform(min = 0.05, max = 2) #0 and 2
  asc_rate ~ Beta(a = 5, b = 1)
  #dispersion ~ Exponential(mean = 1)
})

## Pack the priors
pars_fixed <- pars[-(15:19)]
prior_packer <- monty::monty_packer(c("alpha", "beta", "gamma", "sigma", "asc_rate"), fixed = pars_fixed)

## With this packer we can convert from a list of name-value pairs suitable for
## initialising a dust2 system into a vector of parameters suitable for use with monty:
prior_packer$pack(pars)

## Create data to fit to
positive_tests <- replicate(75, list(c(1, 1, 1)))
data_outbreaks <- data.frame(week = 1:75, outbreak_detected = I(positive_tests))

set.seed(1)

## Build a particle filter
data_week <- dust2::dust_filter_data(data_outbreaks, time = "week")
filter <- dust2::dust_filter_create(cowflu:::cows(), 0, #0 is "time_start"
                                    data_week, n_particles = n_particles, n_threads = 32,
                                    dt=dt)
## Build a likelihood
likelihood <- dust2::dust_likelihood_monty(filter, prior_packer)
## We can combine the prior and the likelihood to create a posterior:
posterior <- prior + likelihood

##Build the sampler
sampler <- monty::monty_sampler_random_walk(diag(5) * step_size_var)  #0.02 was baseline

## Run the samples
samples <- monty::monty_sample(posterior, sampler, n_samples, n_chains = n_chains,
                               initial = prior_packer$pack(pars))
saveRDS(samples, "fitting_samples.rds")

##################################################################################################
## Plotting

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

duration <- toc()
##################################################################################################
## Print an output .txt of the parameters used:
param_string <- sprintf("duration ran: %s mins\n
  n_samples: %s \n
  n_particles: %s \n
  n_chains: %s \n
  dt: %s \n
  step_size_var: %s",   (duration$toc - duration$tic)/60,
                        n_samples,
                        n_particles, n_chains,
                        dt, step_size_var)

fileConn<-file("parameters_used.txt")
writeLines(param_string, fileConn)
close(fileConn)
