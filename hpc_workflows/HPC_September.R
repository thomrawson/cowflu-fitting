library(hipercow)
hipercow_init(driver = "windows")
hipercow_configuration()

id <- task_create_expr(sessionInfo())
task_status(id)
task_result(id)
#task_cancel(id)

hipercow_provision()
windows_authenticate()

hipercow_environment_create(sources = "hpc_functions/particle_density_survival_function.R")

resources <- hipercow_resources(cores = 32)


## 
test <- task_create_expr(particle_density_survival_comparison(n_iterations = 4, 
                                                                      max_particles = 64), 
                         resources = resources)

task_status(test)
task_log_show(test)
task_info(test)$times
## Let's look at the log likelihood:
outputs <- task_result(test)

## Here I ACTUALLY run it
particles_comparison <- task_create_expr(particle_density_survival_comparison(n_iterations = 40, 
                                                              max_particles = 320), 
                         resources = resources)

task_status(particles_comparison)
task_log_show(particles_comparison)
task_info(particles_comparison)$times
## Let's look at the log likelihood:
particles_comparison_outputs <- task_result(particles_comparison)

#####################################################
# FITTING

#15mins
First_test <- task_create_expr(orderly2::orderly_run("fitting-survival", 
                                                     parameters = list(n_samples = 200, 
                                                                       n_particles = 8, 
                                                                       n_chains = 4)), 
                               resources = resources)
task_status(First_test)
task_log_show(First_test)
task_info(First_test)$times

###################################
## These are shorter to assess the impact of the cvc matrix variance (step size)
# 20240926-221923-d9f76503
# 9.6hrs
low_var <- task_create_expr(orderly2::orderly_run("fitting-survival", 
                                                     parameters = list(n_samples = 2000, 
                                                                       n_particles = 160, 
                                                                       n_chains = 4,
                                                                       step_size_var = 0.01)), 
                               resources = resources)
task_status(low_var)
task_log_show(low_var)
task_info(low_var)$times
############
#20240926-221938-3b5ddb50
#9.6hrs
mid_var <- task_create_expr(orderly2::orderly_run("fitting-survival", 
                                                  parameters = list(n_samples = 2000, 
                                                                    n_particles = 160, 
                                                                    n_chains = 4,
                                                                    step_size_var = 0.03)), 
                            resources = resources)
task_status(mid_var)
task_log_show(mid_var)
task_info(mid_var)$times
############
#20240926-221944-8d6f9d9a
#9.7hrs
high_var <- task_create_expr(orderly2::orderly_run("fitting-survival", 
                                                  parameters = list(n_samples = 2000, 
                                                                    n_particles = 160, 
                                                                    n_chains = 4,
                                                                    step_size_var = 0.06)), 
                            resources = resources)
task_status(high_var)
task_log_show(high_var)
task_info(high_var)$times

##################################
#These are actually meant to be useful
First_Fit <- task_create_expr(orderly2::orderly_run("fitting-survival", 
                                                   parameters = list(n_samples = 10000, 
                                                                     n_particles = 160, 
                                                                     n_chains = 4,
                                                                     step_size_var = 0.03)), 
                             resources = resources)
task_status(First_Fit)
task_log_show(First_Fit)
task_info(First_Fit)$times

###########################
Second_Fit_lower_dt <- task_create_expr(orderly2::orderly_run("fitting-survival", 
                                                    parameters = list(n_samples = 5000, 
                                                                      n_particles = 160, 
                                                                      n_chains = 4,
                                                                      dt = 0.25,
                                                                      step_size_var = 0.03)), 
                              resources = resources)
task_status(Second_Fit_lower_dt)
task_log_show(Second_Fit_lower_dt)
task_info(Second_Fit_lower_dt)$times
##############################
Third_Fit_BIG <- task_create_expr(orderly2::orderly_run("fitting-survival", 
                                                        parameters = list(n_samples = 40000, 
                                                                          n_particles = 240, 
                                                                          n_chains = 4,
                                                                          step_size_var = 0.03)), 
                                  resources = resources)
task_status(Third_Fit_BIG)
task_log_show(Third_Fit_BIG)
task_info(Third_Fit_BIG)$times

###################################
###################################
## These are to assess the impact of dt
# 20240927-113914-9af605e2
# 
lowest_dt <- task_create_expr(orderly2::orderly_run("fitting-survival", 
                                                  parameters = list(n_samples = 2000, 
                                                                    n_particles = 160, 
                                                                    n_chains = 4,
                                                                    dt = 1/14,
                                                                    step_size_var = 0.03)), 
                            resources = resources)
task_status(lowest_dt)
task_log_show(lowest_dt)
task_info(lowest_dt)$times
############
# 20240927-114009-da52444d
#
low_dt <- task_create_expr(orderly2::orderly_run("fitting-survival", 
                                                    parameters = list(n_samples = 2000, 
                                                                      n_particles = 160, 
                                                                      n_chains = 4,
                                                                      dt = 1/7,
                                                                      step_size_var = 0.03)), 
                              resources = resources)
task_status(low_dt)
task_log_show(low_dt)
task_info(low_dt)$times
############
#20240927-114136-8344eaf5
#
mid_dt <- task_create_expr(orderly2::orderly_run("fitting-survival", 
                                                 parameters = list(n_samples = 2000, 
                                                                   n_particles = 160, 
                                                                   n_chains = 4,
                                                                   dt = 1/4,
                                                                   step_size_var = 0.03)), 
                           resources = resources)
task_status(mid_dt)
task_log_show(mid_dt)
task_info(mid_dt)$times

##################################
#20240927-114214-4e73bd81
#
high_dt <- task_create_expr(orderly2::orderly_run("fitting-survival", 
                                                 parameters = list(n_samples = 2000, 
                                                                   n_particles = 160, 
                                                                   n_chains = 4,
                                                                   dt = 1/2,
                                                                   step_size_var = 0.03)), 
                           resources = resources)
task_status(high_dt)
task_log_show(high_dt)
task_info(high_dt)$times

##########################################
##########################################
#This is me doing some big runs of incidence, varying step-size again
#
#
First_Incidence <- task_create_expr(orderly2::orderly_run("fitting-incidence", 
                                                    parameters = list(n_samples = 5000, 
                                                                      n_particles = 160, 
                                                                      n_chains = 4,
                                                                      step_size_var = 0.01)), 
                              resources = resources)
task_status(First_Incidence)
task_log_show(First_Incidence)
task_info(First_Incidence)$times
############
#
#
Second_Incidence <- task_create_expr(orderly2::orderly_run("fitting-incidence", 
                                                          parameters = list(n_samples = 5000, 
                                                                            n_particles = 160, 
                                                                            n_chains = 4,
                                                                            step_size_var = 0.01)), 
                                    resources = resources)
task_status(Second_Incidence)
task_log_show(Second_Incidence)
task_info(Second_Incidence)$times
#############
#
#
Third_Incidence <- task_create_expr(orderly2::orderly_run("fitting-incidence", 
                                                          parameters = list(n_samples = 5000, 
                                                                            n_particles = 160, 
                                                                            n_chains = 4,
                                                                            step_size_var = 0.01)), 
                                    resources = resources)
task_status(Third_Incidence)
task_log_show(Third_Incidence)
task_info(Third_Incidence)$times
