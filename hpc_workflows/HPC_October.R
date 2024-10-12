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

##################################################################################################
#Latin Hypercube Sampling
#########################
#########################
# (fixed) 20241003-122422-3e17c00e
# 11.7hrs
# New priors: 20241004-102019-1ae5b975
small_LHS <- task_create_expr(orderly2::orderly_run("LHC-survival", 
                                                  parameters = list(LHS_samples = 10000, 
                                                                     n_particles = 160, 
                                                                    n_iterations = 1,
                                                                    dt = 1)), 
                            resources = resources)
task_status(small_LHS)
task_log_show(small_LHS)
task_info(small_LHS)$times
#########################
# (fixed) 20241003-141306-12f4f5e2
# 32hrs
big_LHS <- task_create_expr(orderly2::orderly_run("LHC-survival", 
                                                    parameters = list(LHS_samples = 40000, 
                                                                      n_particles = 160, 
                                                                      n_iterations = 1,
                                                                      dt = 1)), 
                              resources = resources)
task_status(big_LHS)
task_log_show(big_LHS)
task_info(big_LHS)$times
#########################
# (Oct6) 
# 
biggest_LHS <- task_create_expr(orderly2::orderly_run("LHC-survival", 
                                                  parameters = list(LHS_samples = 40000, 
                                                                    n_particles = 160, 
                                                                    n_iterations = 4,
                                                                    dt = 1)), 
                            resources = resources)
task_status(biggest_LHS)
task_log_show(biggest_LHS)
task_info(biggest_LHS)$times
#########################
# Try 1 particle
# (Oct6) 
#
humour_me <- task_create_expr(orderly2::orderly_run("LHC-survival", 
                                                      parameters = list(LHS_samples = 40000, 
                                                                        n_particles = 1, 
                                                                        n_iterations = 100,
                                                                        dt = 1)), 
                                resources = resources)
task_status(humour_me)
#task_cancel(humour_me)
task_log_show(humour_me)
task_info(humour_me)$times

#######################
# AND FOR INCIDENCE
#########################
# (fixed) 20241003-135029-4ee17fd8
#
# New Priors: 20241004-102058-d76f3e3d
small_LHS_inc <- task_create_expr(orderly2::orderly_run("LHC-incidence", 
                                                    parameters = list(LHS_samples = 10000, 
                                                                      n_particles = 160, 
                                                                      n_iterations = 1,
                                                                      dt = 1)), 
                              resources = resources)
task_status(small_LHS_inc)
task_log_show(small_LHS_inc)
task_info(small_LHS_inc)$times
#########################
# (fixed) 20241003-152153-446153e3
# 
big_LHS_inc <- task_create_expr(orderly2::orderly_run("LHC-incidence", 
                                                  parameters = list(LHS_samples = 40000, 
                                                                    n_particles = 160, 
                                                                    n_iterations = 1,
                                                                    dt = 1)), 
                            resources = resources)
task_status(big_LHS_inc)
task_log_show(big_LHS_inc)
task_info(big_LHS_inc)$times
#########################
# (Oct6) 
# 
biggest_LHS_inc <- task_create_expr(orderly2::orderly_run("LHC-incidence", 
                                                      parameters = list(LHS_samples = 40000, 
                                                                        n_particles = 160, 
                                                                        n_iterations = 4,
                                                                        dt = 1)), 
                                resources = resources)
task_status(biggest_LHS_inc)
task_log_show(biggest_LHS_inc)
task_info(biggest_LHS_inc)$times
#########################
# Try 1 particle
# (Oct6) 
#
humour_me_inc <- task_create_expr(orderly2::orderly_run("LHC-incidence", 
                                                    parameters = list(LHS_samples = 40000, 
                                                                      n_particles = 1, 
                                                                      n_iterations = 100,
                                                                      dt = 1)), 
                              resources = resources)
task_status(humour_me_inc)
task_log_show(humour_me_inc)
task_info(humour_me_inc)$times

#####################################################
#####################################################
#####################################################
# FITTING
# (Oct6) 20241007-141740-405288b0
#
First_test <- task_create_expr(orderly2::orderly_run("fitting-survival", 
                                                     parameters = list(n_samples = 200, 
                                                                       n_particles = 8, 
                                                                       n_chains = 4,
                                                                       restart = TRUE,
                                                                       rerun_every = 50)), 
                               resources = resources)
task_status(First_test)
task_log_show(First_test)
task_info(First_test)$times

##################################
#These are actually meant to be useful
# (Oct6) 20241007-141801-ffbae18a
# 42hrs
First_Fit <- task_create_expr(orderly2::orderly_run("fitting-survival", 
                                                    parameters = list(n_samples = 10000, 
                                                                      n_particles = 160, 
                                                                      n_chains = 4,
                                                                      step_size_var = 0.03,
                                                                      restart = TRUE,
                                                                      rerun_every = 50)), 
                              resources = resources)
task_status(First_Fit)
task_log_show(First_Fit)
task_info(First_Fit)$times

###########################
# (Oct6) 20241007-141830-a270c157
# 
Second_Fit_lower_dt <- task_create_expr(orderly2::orderly_run("fitting-survival", 
                                                              parameters = list(n_samples = 5000, 
                                                                                n_particles = 160, 
                                                                                n_chains = 4,
                                                                                dt = 0.25,
                                                                                step_size_var = 0.03,
                                                                                restart = TRUE,
                                                                                rerun_every = 50)), 
                                        resources = resources)
task_status(Second_Fit_lower_dt)
task_log_show(Second_Fit_lower_dt)
task_info(Second_Fit_lower_dt)$times
##############################
#(Oct6)
# 
Third_Fit_BIG <- task_create_expr(orderly2::orderly_run("fitting-survival", 
                                                        parameters = list(n_samples = 40000, 
                                                                          n_particles = 200, 
                                                                          n_chains = 4,
                                                                          step_size_var = 0.03,
                                                                          restart = TRUE,
                                                                          rerun_every = 50)), 
                                  resources = resources)
task_status(Third_Fit_BIG)
task_log_show(Third_Fit_BIG)
task_info(Third_Fit_BIG)$times

###################################
###################################
##########################################
##########################################
#This is me doing some big runs of incidence, 
# (Oct6) 20241007-141927-a1d2f463
#
First_Incidence <- task_create_expr(orderly2::orderly_run("fitting-incidence", 
                                                          parameters = list(n_samples = 10000, 
                                                                            n_particles = 240, 
                                                                            n_chains = 4,
                                                                            restart = TRUE,
                                                                            rerun_every = 50)), 
                                    resources = resources)
task_status(First_Incidence)
task_log_show(First_Incidence)
task_info(First_Incidence)$times
############
#(Oct6)
#
Second_Incidence <- task_create_expr(orderly2::orderly_run("fitting-incidence", 
                                                           parameters = list(n_samples = 20000, 
                                                                             n_particles = 240, 
                                                                             n_chains = 4,
                                                                             restart = TRUE,
                                                                             rerun_every = 50)), 
                                     resources = resources)
task_status(Second_Incidence)
task_log_show(Second_Incidence)
task_info(Second_Incidence)$times
#############
#(Oct6)
#
Third_Incidence <- task_create_expr(orderly2::orderly_run("fitting-incidence", 
                                                          parameters = list(n_samples = 30000, 
                                                                            n_particles = 240, 
                                                                            n_chains = 4,
                                                                            dt = 0.5,
                                                                            restart = TRUE,
                                                                            rerun_every = 50)), 
                                    resources = resources)
task_status(Third_Incidence)
task_log_show(Third_Incidence)
task_info(Third_Incidence)$times

#############
#(Oct6)
#
Fourth_Incidence <- task_create_expr(orderly2::orderly_run("fitting-incidence", 
                                                          parameters = list(n_samples = 40000, 
                                                                            n_particles = 240, 
                                                                            n_chains = 4,
                                                                            restart = TRUE,
                                                                            rerun_every = 50)), 
                                    resources = resources)
task_status(Fourth_Incidence)
task_log_show(Fourth_Incidence)
task_info(Fourth_Incidence)$times

######################################