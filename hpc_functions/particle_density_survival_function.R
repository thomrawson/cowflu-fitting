#Actual Model - WEEKLY FITTING

particle_density_survival_comparison <- function( dt = 1, n_iterations = 40, max_particles = 320){
  
  library(cowflu)
  
  ## Identifying optimal n_particles:
  
  ## We're going to measure how variance in density estimates decreases as n_particles goes up.
  
  ## Data frame with a col for n_particles (make it be a multiple of.... 8? 32?)
  ## and cols for mean, sd, and var of density
  dt_to_run <- dt
  particle_numbers_to_run <- n_iterations
  
  n_particles_density_df <- data.frame(n_particles = seq(from = 8, to = max_particles, by = max_particles/particle_numbers_to_run),
                                       mean_density = rep(NA,particle_numbers_to_run),
                                       mean_density_lower_ci = rep(NA,particle_numbers_to_run),
                                       mean_density_upper_ci = rep(NA,particle_numbers_to_run),
                                       sd_density = rep(NA,particle_numbers_to_run),
                                       var_density = rep(NA,particle_numbers_to_run),
                                       dt = rep(dt_to_run, particle_numbers_to_run))
  
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
  
  for(i in 1:particle_numbers_to_run){
    print(i)
    n_particles_to_run <- n_particles_density_df$n_particles[i]
    
    filter <- dust2::dust_filter_create(cowflu:::cows(), 0, data_week, n_particles = n_particles_to_run, n_threads = 8,
                                        dt=dt_to_run)
    
    hold_densities <- rep(NA, 500)
    
    for(j in 1:length(hold_densities)){
      hold_densities[j] <- dust2::dust_likelihood_run(filter, pars)
    }
    n_particles_density_df$mean_density[i] <- mean(hold_densities)
    n_particles_density_df$mean_density_lower_ci[i] <- quantile(hold_densities, probs = c(0.025))
    n_particles_density_df$mean_density_upper_ci[i] <- quantile(hold_densities, probs = c(0.975))
    n_particles_density_df$sd_density[i] <- sd(hold_densities)
    n_particles_density_df$var_density[i] <- var(hold_densities)
    
  }
  
  library(ggplot2)
  ggplot(n_particles_density_df) +
    geom_line(aes(x = n_particles, y = var_density), size = 1) +
    theme_classic() +
    xlab("Number of particles") +
    ylab("Variance in density") +
    ggtitle("Impact of particles on density variance for Survival model") -> var_plot
  
  ggsave("var_plot_survival.png", var_plot, width = 20, height = 20, units = "cm")
  
  ggplot(n_particles_density_df) +
    geom_line(aes(x = log(n_particles), y = log(var_density)), size = 1) +
    theme_classic() +
    xlab("log( Number of particles )") +
    ylab("log( Variance ) in density") +
    ggtitle("Impact of particles on density variance for Survival model") -> log_var_plot
  
  ggsave("log_var_plot_survival.png", log_var_plot, width = 20, height = 20, units = "cm")
  
  ggplot(n_particles_density_df) +
    geom_line(aes(x = n_particles, y = mean_density), size = 1) +
    geom_ribbon(aes(x = n_particles, ymin = mean_density_lower_ci, ymax = mean_density_upper_ci),
                alpha = 0.2) +
    theme_classic() +
    xlab("Number of particles") +
    ylab("Mean density") +
    ggtitle("Impact of particles on mean for Survival model") -> mean_plot
  
  ggsave("mean_plot_survival.png", mean_plot, width = 20, height = 20, units = "cm")
  
  outputs_list <- list(data = n_particles_density_df, var_plot = var_plot, log_var_plot = log_var_plot, mean_plot = mean_plot)
  
  return(outputs_list)
  

  
}
