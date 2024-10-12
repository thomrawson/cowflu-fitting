orderly2::orderly_parameters(n_samples = 40, n_particles = 2, n_chains = 2, step_size_var = 0.03, dt = 1, restart = TRUE, rerun_every = 100)
params <- c("alpha", "beta", "gamma", "sigma", "asc_rate")
##################################################################################################
##Assume these package versions are installed:
#remotes::install_github("mrc-ide/dust2@mrc-5713-2", upgrade = FALSE)
#remotes::install_github("mrc-ide/monty@rerun-and-observe", upgrade = FALSE)
#remotes::install_github("mrc-ide/cowflu@new-dust", upgrade = FALSE)

library(cowflu)
library(dust2)
library(monty)
##################################################################################################
## Set parameters
pars <- cowflu:::cowflu_inputs(
  alpha = 0.05,
  beta = 2.45,
  gamma = 1.4,
  sigma = 3.5,
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
    likelihood_choice = "survival"))

## Set priors
prior <- monty::monty_dsl({
  alpha ~ Beta(a = 1, b = 25)
  beta ~ Uniform(min = 1, max = 3) #maybe 1 and 2.5
  gamma ~ Uniform(min = 0.05, max = 2) #0 and 1
  sigma ~ Uniform(min = 0.05, max = 4) #0 and 2
  asc_rate ~ Beta(a = 5, b = 1)
})

## Pack the priors
pars_fixed <- pars[-(15:19)]
prior_packer <- monty::monty_packer(c("alpha", "beta", "gamma", "sigma", "asc_rate"), fixed = pars_fixed)

## With this packer we can convert from a list of name-value pairs suitable for
## initialising a dust2 system into a vector of parameters suitable for use with monty:
prior_packer$pack(pars)

## Load data to fit to
data_outbreaks <- cowflu:::process_data_outbreak(cowflu:::outbreaks_data$weekly_outbreaks_data)
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
                                           save_history = c("outbreak_region"))

## We can combine the prior and the likelihood to create a posterior:
posterior <- prior + likelihood

##Build the sampler
if(restart){
  sampler <- monty::monty_sampler_random_walk(diag(5) * step_size_var,
                                              rerun_every = rerun_every,
                                              rerun_random = TRUE)  #0.02 was baseline
}else{
  sampler <- monty::monty_sampler_random_walk(diag(5) * step_size_var)  #0.02 was baseline
}


## Run the samples
samples <- monty::monty_sample(posterior, sampler, n_samples, n_chains = n_chains,
                               initial = prior_packer$pack(pars) )

## Here currently throws an error due to the dimension mis-match:
unpacked_state <- filter$packer_history$unpack(samples$observations$history)
