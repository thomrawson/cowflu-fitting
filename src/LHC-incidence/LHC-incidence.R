##################################################################################################
orderly2::orderly_strict_mode()
orderly2::orderly_parameters(LHS_samples = 10000, n_particles = 160, n_iterations = 1, dt = 1)

params <- c("alpha", "beta", "gamma", "sigma", "asc_rate", "dispersion")
orderly2::orderly_artefact(description = "The dataframe", "Parameters_df.rds")
orderly2::orderly_artefact(description = "plot", "Alpha.png")
orderly2::orderly_artefact(description = "plot", "Beta.png")
orderly2::orderly_artefact(description = "plot", "Gamma.png")
orderly2::orderly_artefact(description = "plot", "Sigma.png")
orderly2::orderly_artefact(description = "plot", "asc_rate.png")
orderly2::orderly_artefact(description = "plot", "dispersion.png")

orderly2::orderly_artefact(description = "PCA cumulative proportions", "pca_importance.txt")

orderly2::orderly_artefact(description = "The parameters used", "parameters_used.txt")

##################################################################################################
library(cowflu)
library(dust2)
library(monty)
library(tictoc)
library(lhs)
library(ggplot2)
tic()
##################################################################################################
## Generate the parameter values via Latin Hypercube (LHC) sampling:
## randomLHS(samples, parameters)

Parameters_df <- as.data.frame(randomLHS(LHS_samples, length(params)) )
colnames(Parameters_df) <- params

## Make an empty column for the log-likelihood values:
Parameters_df$density <- rep(NA, LHS_samples)

## Now we scale those LHC values by the intended parameter ranges:
# alpha ~ Beta(a = 1, b = 25)
# beta ~ Uniform(min = 1, max = 3) #maybe 1 and 2.5
# gamma ~ Uniform(min = 0.05, max = 2) #0 and 1
# sigma ~ Uniform(min = 0.05, max = 4) #0 and 2
# asc_rate ~ Beta(a = 5, b = 1)
# dispersion ~ Exponential(mean = 0.5)

## For the Beta distr. we use the inverse cumulative distribution function (CDF)
## (also known as the quantile function) of the Beta distribution.
## Using the qbeta() function, which returns the quantile corresponding to a given
## probability (in this case, the uniform values)
Parameters_df$alpha <- qbeta(Parameters_df$alpha , shape1 = 1, shape2 = 7)
Parameters_df$beta <- 1 + (4.5 - 1) * Parameters_df$beta
Parameters_df$gamma <- 0 + (3.5 - 0) * Parameters_df$gamma
Parameters_df$sigma <- 0 + (6.5 - 0) * Parameters_df$sigma
Parameters_df$asc_rate <- qbeta(Parameters_df$asc_rate , shape1 = 1, shape2 = 1)
Parameters_df$dispersion <- 1 + qexp(Parameters_df$dispersion, rate = 1/2)

## Prepare the data
data <- cowflu:::process_data_incidence(cowflu:::outbreaks_data$weekly_outbreaks_data)
data_week <- dust2::dust_filter_data(data, time = "week")

##Prep progress bar
pb <- txtProgressBar(min = 0, max = LHS_samples, style = 3)

filter <- dust2::dust_filter_create(cowflu:::cows(), 0, data_week, n_particles = n_particles, n_threads = 32,
                                    dt=dt)

## Now we can calculate the log-likelihood for each set of parameters:
for( i in 1:LHS_samples){
  pars <- cowflu:::cowflu_inputs(alpha = Parameters_df$alpha[i],
                                 beta = Parameters_df$beta[i],
                                 gamma = Parameters_df$gamma[i],
                                 sigma = Parameters_df$sigma[i],
                                 asc_rate = Parameters_df$asc_rate[i],
                                 dispersion = Parameters_df$dispersion[i],
                                 cowflu:::cowflu_fixed_inputs(p_region_export = cowflu:::movement$p_region_export,
                                                              p_cow_export = cowflu:::movement$p_cow_export,
                                                              movement_matrix = cowflu:::movement$movement_matrix,
                                                              time_test = 19,  #Was day 136
                                                              n_herds_per_region = cowflu:::usda_data$n_herds_per_region,
                                                              n_cows_per_herd = cowflu:::usda_data$n_cows_per_herd,
                                                              start_herd = 26940, #26804 is where Texas starts.
                                                              start_count = 5,
                                                              condition_on_export = TRUE,
                                                              likelihood_choice = "incidence"
                                 ))

  hold_density <- 0
  for(j in 1:n_iterations){
    hold_density <- hold_density + (dust2::dust_likelihood_run(filter, pars)/n_iterations)
  }

  ## Save the log-likelihood value:
  Parameters_df$density[i] <- hold_density

  ## Update the progress bar:
  setTxtProgressBar(pb, i)
}

## Close the progress bar
close(pb)

saveRDS(Parameters_df, "Parameters_df.rds")
#######################################
## Now make plots with the results:
######################################

## Perform PCA
pca_result <- prcomp(Parameters_df[, 1:6], scale. = TRUE)  # Replace with actual parameter columns
#  pca_df <- as.data.frame(pca_result$x)
# #summary(pca_result)
# # Add the log-likelihood column back
# pca_df$log_likelihood <- Parameters_df$density
#
# library(ggplot2)
# # PCA plot
# ggplot(pca_df, aes(x = PC1, y = PC2, color = log_likelihood)) +
#   geom_point() +
#   scale_color_viridis_c() +
#   labs(x = "PC1", y = "PC2", color = "Log Likelihood") +
#   ggtitle("PCA of Parameter Space Colored by Log Likelihood") +
#   theme_minimal()


# ALPHA
# Fit a quadratic model: y ~ x + x^2
quadraticModel <- lm(density ~ alpha + I(alpha^2), data=Parameters_df)
# Generate a sequence of x values for a smooth curve
x_vals <- seq(min(Parameters_df$alpha), max(Parameters_df$alpha), length.out = 100)
# Predict y values using the quadratic model
yPredict <- predict(quadraticModel, list(alpha = x_vals ))

ggplot() +
  geom_point(data = Parameters_df, aes(x = alpha, y = density), alpha = 0.5) +
  geom_line(aes(x = x_vals, y = yPredict), color = "red", size = 1, alpha = 0.4) +
  labs(title = "Alpha",
       x = "Alpha",
       y = "Density") +
  theme_classic() -> plot

ggsave(filename = "Alpha.png", plot = plot, width = 8, height = 6, units = "in", dpi = 300)
########
# BETA
quadraticModel <- lm(density ~ beta + I(beta^2), data=Parameters_df)
x_vals <- seq(min(Parameters_df$beta), max(Parameters_df$beta), length.out = 100)
yPredict <- predict(quadraticModel, list(beta = x_vals ))

ggplot() +
  geom_point(data = Parameters_df, aes(x = beta, y = density), alpha = 0.5) +
  geom_line(aes(x = x_vals, y = yPredict), color = "red", size = 1, alpha = 0.4) +
  labs(title = "Beta",
       x = "Beta",
       y = "Density") +
  theme_classic() -> plot

ggsave(filename = "Beta.png", plot = plot, width = 8, height = 6, units = "in", dpi = 300)
##########
# GAMMA
quadraticModel <- lm(density ~ gamma + I(gamma^2), data=Parameters_df)
x_vals <- seq(min(Parameters_df$gamma), max(Parameters_df$gamma), length.out = 100)
yPredict <- predict(quadraticModel, list(gamma = x_vals ))

ggplot() +
  geom_point(data = Parameters_df, aes(x = gamma, y = density), alpha = 0.5) +
  geom_line(aes(x = x_vals, y = yPredict), color = "red", size = 1, alpha = 0.4) +
  labs(title = "Gamma",
       x = "Gamma",
       y = "Density") +
  theme_classic() -> plot

ggsave(filename = "Gamma.png", plot = plot, width = 8, height = 6, units = "in", dpi = 300)
##########
# SIGMA
quadraticModel <- lm(density ~ sigma + I(sigma^2), data=Parameters_df)
x_vals <- seq(min(Parameters_df$sigma), max(Parameters_df$sigma), length.out = 100)
yPredict <- predict(quadraticModel, list(sigma = x_vals ))

ggplot() +
  geom_point(data = Parameters_df, aes(x = sigma, y = density), alpha = 0.5) +
  geom_line(aes(x = x_vals, y = yPredict), color = "red", size = 1, alpha = 0.4) +
  labs(title = "Sigma",
       x = "Sigma",
       y = "Density") +
  theme_classic() -> plot

ggsave(filename = "Sigma.png", plot = plot, width = 8, height = 6, units = "in", dpi = 300)
##########
# ASC_RATE
quadraticModel <- lm(density ~ asc_rate + I(asc_rate^2), data=Parameters_df)
x_vals <- seq(min(Parameters_df$asc_rate), max(Parameters_df$asc_rate), length.out = 100)
yPredict <- predict(quadraticModel, list(asc_rate = x_vals ))

ggplot() +
  geom_point(data = Parameters_df, aes(x = asc_rate, y = density), alpha = 0.5) +
  geom_line(aes(x = x_vals, y = yPredict), color = "red", size = 1, alpha = 0.4) +
  labs(title = "Ascertainment Rate",
       x = "asc_rate",
       y = "Density") +
  theme_classic() -> plot

ggsave(filename = "asc_rate.png", plot = plot, width = 8, height = 6, units = "in", dpi = 300)
##########
# DISPERSION
quadraticModel <- lm(density ~ dispersion + I(dispersion^2), data=Parameters_df)
x_vals <- seq(min(Parameters_df$dispersion), max(Parameters_df$dispersion), length.out = 100)
yPredict <- predict(quadraticModel, list(dispersion = x_vals ))

ggplot() +
  geom_point(data = Parameters_df, aes(x = dispersion, y = density), alpha = 0.5) +
  geom_line(aes(x = x_vals, y = yPredict), color = "red", size = 1, alpha = 0.4) +
  labs(title = "Dispersion hyperparameter",
       x = "dispersion",
       y = "Density") +
  theme_classic() -> plot

ggsave(filename = "dispersion.png", plot = plot, width = 8, height = 6, units = "in", dpi = 300)
##########


duration <- toc()
##################################################################################################
## Print an output .txt of the parameters used:
param_string <- sprintf("duration ran: %s mins\n
  LHS_samples: %s \n
  n_iterations: %s \n
  n_particles: %s \n
  dt: %s ",   (duration$toc - duration$tic)/60,
                        LHS_samples,
                        n_iterations,
                        n_particles,
                        dt)

fileConn<-file("parameters_used.txt")
writeLines(param_string, fileConn)
close(fileConn)

pca_result_text <- summary(pca_result)$importance
write.table(pca_result_text, file = "pca_importance.txt", sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


