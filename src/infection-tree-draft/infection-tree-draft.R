## We want to calculate, for each confirmed outbreak, what was the most likely SOURCE of that initial outbreak?
## There's three parts to this:
## 1) Which of our "TRUE" modelled outbreaks COULD have been the source? i.e. how many weeks back do we want to go?
## 2) What is the probability of that "shipment" coming? This is based on the movement matrix
## 3) What is the probability of an actual infected herd being that source herd?
## 4) What is the probability that it was a "spillover" from the environment, and not trade-related? (This will be the hardest part)

## I want to take some samples object as the input, and extract the "TRUE INFECTED HERDS" data from it.
## For now, I've just copied the samples from: 20241012-144052-ac490bf5
## n_samples: 10000
## n_particles: 160
## n_chains: 4
## dt: 1
## step_size_var: 0.03
## restarts:  TRUE
## restart_every: 50

fitting_samples <- readRDS("fitting_samples.rds")
## First, let's extract the mean and 95%CIs for the model parameters:
## Easiest with coda:
library(coda)
samples_coda <- coda::as.mcmc.list(fitting_samples)
summary_stats <- summary(samples_coda)

## Now let's load in the actual data
outbreaks_data <- cowflu:::outbreaks_data$weekly_outbreaks_data
## Filter out all zeros:
library(dplyr)
outbreaks_data <- filter(outbreaks_data, Positive_Tests > 0)
start_date <- as.Date("2023-12-18")
outbreaks_data$day <- as.numeric(outbreaks_data$Week_Beginning - start_date)
outbreaks_data$Week <- outbreaks_data$day/7

##So let's take one as an example.
## Kansas, has 1 outbreak on week 17 day 119

## First question, who COULD have been infectious at the right time to infect?
## Let's assume it can't be the current week.
## So just the last week?
## Sigma is the incubation period
## Gamma is the recovery rate.
## Sigma has a mean of 1.5/week, i.e. 0.214 chance of recovery per day (sort of...)
## Gamma has a mean 1.4/week

