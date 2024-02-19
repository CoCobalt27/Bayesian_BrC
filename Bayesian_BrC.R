# Bayesian inference algorithm to resolve BrC conc. and its AAE
# Parallel computation for speed boost
# Notes: `future` and `furrr` are packages for parallel computing

# 0. Packages ####
library(tidyverse)
library(readxl)
library(LaplacesDemon)
library(lubridate)
library(future)
library(furrr)

# 1. Data import ####
load("./measurements.RData")

# To show mass absorption efficiency (MAE) data of BC
mae

# The prior distribution parameters of MAE data of BrC
# Note: MAE = K * wavelength^(-AAE)
BrC_K <- 290000 
BrC_AAE <- 1.84

# 2. Model-I ####
# Model input: AE-33 + TCA data
# Model assumption: MAEs of BC are fixed and from reference samples
# Model output: for each observation (K_BrC, AAE_BrC, [BC], [BrC], [WtC])

## 2.1 Functions ####
# Likelihood function
likelihood <- function(i, y){
  # Check on inputs
  if(length(y) != 5) stop("Length of input is wrong!")
  
  # BrC MAE
  names(y) <- c("K", "AAE", "BC", "BrC", "WtC")
  # Predicted absorption
  mae_BrC <- y["K"] * (mae$wavelength_nm ^ -y["AAE"])
  abs_pred <- mae$BC_abs * y["BC"] + mae_BrC * y["BrC"]
  abs_meas <- unlist(dt[i, 3:9])
  abs_unc  <- unlist(unc[i, 3:9])
  dmvn(abs_meas, abs_pred, diag(abs_unc^2))
}

# Prior distribution
prior <- function(y){
  # Check on inputs
  if(length(y) != 5) stop("Length of input is wrong!")
  # BrC MAE
  names(y) <- c("K", "AAE", "BC", "BrC", "WtC")
  
  # Parameters for prior distribution
  K_avg <- BrC_K
  K_sd <- K_avg * 0.05
  AAE_expected <- BrC_AAE
  
  dnorm(y["K"], K_avg, K_sd) * dexp(y["AAE"], 1/(AAE_expected - 1))
}

# Reference distribution
q0 <- function(i){
  runif(1, BrC_K * 0.95,  BrC_K * 1.05) -> K
  runif(1, 1, 5) -> AAE
  as.vector(rdirichlet(1, c(1, 1, 1))) * dt$TC[i] -> conc
  c(K, AAE, conc)
}

# MCMC sampling starts here:
MCMC <- function(j){
  if(j%%5 == 0) message(j)
  length_of_mc <- 5000
  
  y0 <- q0(j)
  mc_samples <- matrix(rep(y0, length_of_mc), ncol = 5, byrow = TRUE)
  
  for(i in 1:(length_of_mc - 1)){
    y_o <- unlist(mc_samples[i, ])
    y_n <- q0(j)
    
    threshold <- min(c(1, prior(y_n)*likelihood(j, y_n)/(prior(y_o)*likelihood(j, y_o))))
    
    rand <- runif(1)
    
    judge <- (rand < threshold)
    
    ifelse(is.na(judge), FALSE, judge) -> judge
    
    if(judge){
      mc_samples[i + 1, ] <- y_n
    } else {
      mc_samples[i + 1, ] <- y_o
    }
  }
  
  mc_samples %>%
    as.data.frame() %>%
    rename(K = 1, AAE = 2, BC = 3, BrC = 4, WtC = 5) %>%
    slice(seq(round(length_of_mc/5), length_of_mc, 25)) -> mc_samples
  
  rbind(summarise_all(mc_samples, median),
        summarise_all(mc_samples, mad)) %>%
  mutate(key = c("avg", "sd"), index = j)
}

## 2.2 Parallel computation ####
# The `workers` parameter should be less than the number of cores in your PC
plan(multisession, workers = 2)

system.time(
            do.call(rbind,
                    future_imap(1:nrow(dt), ~ MCMC(.x),
                                .options = furrr_options(seed = 2))) -> mod_result
            )

mod_result %>%
  mutate(Time = rep(dt$Time, each = 2)) %>%
  select(Time, key, everything()) %>%
  select(-index) -> mod_result

#### END ####