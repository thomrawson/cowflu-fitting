##################################################################################################
orderly2::orderly_strict_mode()
orderly2::orderly_parameters(n_particles = 160, n_iterations = 10, dt = 1)

params <- c("alpha", "beta", "gamma", "sigma", "asc_rate")
orderly2::orderly_artefact(description = "The dataframe", "n_particles_density_df.rds")

orderly2::orderly_artefact(description = "The parameters used", "parameters_used.txt")

##################################################################################################
library(cowflu)
library(dust2)
library(monty)
library(tictoc)
tic()
##################################################################################################
n_particles_density_df <- data.frame(n_particles = rep(n_particles, n_iterations),
                                     density = rep(NA, n_iterations),
                                     dt = rep(dt, n_iterations))

pars <- cowflu:::cowflu_inputs(alpha = 0.05, beta = 2.45, gamma = 1.4, sigma = 3.5, asc_rate = 1,
                               dispersion = 1,
                               cowflu:::cowflu_fixed_inputs(p_region_export = cowflu:::movement$p_region_export,
                                                            p_cow_export = cowflu:::movement$p_cow_export,
                                                            movement_matrix = cowflu:::movement$movement_matrix,
                                                            time_test = 19,  #Was day 136
                                                            n_herds_per_region = cowflu:::usda_data$n_herds_per_region,
                                                            n_cows_per_herd = cowflu:::usda_data$n_cows_per_herd,
                                                            start_herd = 26940, #26804 is where Texas starts.
                                                            start_count = 5,
                                                            condition_on_export = TRUE,
                                                            likelihood_choice = "survival"
                               ))

data <- cowflu:::process_data_outbreak(cowflu:::outbreaks_data$weekly_outbreaks_data)
data_week <- dust2::dust_filter_data(data, time = "week")

for(i in 1:n_iterations){
  if( (i %% 10) == 0){
    print(i)
  }

  filter <- dust2::dust_filter_create(cowflu:::cows(), 0, data_week, n_particles = n_particles, n_threads = 8,
                                      dt=dt)

  hold_density <- dust2::dust_likelihood_run(filter, pars)

  n_particles_density_df$density[i] <- hold_density

}

saveRDS(n_particles_density_df, "n_particles_density_df.rds")

duration <- toc()
##################################################################################################
## Print an output .txt of the parameters used:
param_string <- sprintf("duration ran: %s mins\n
  n_iterations: %s \n
  n_particles: %s \n
  dt: %s ",   (duration$toc - duration$tic)/60,
                        n_iterations,
                        n_particles,
                        dt)

fileConn<-file("parameters_used.txt")
writeLines(param_string, fileConn)
close(fileConn)

