# Load packages
library(haven)
library(lavaan)
library(dplyr)
library(tidyr)
library(semtree)
library(ggplot2)
library(reshape2)
library(mvtnorm)
library(devtools)
library(OpenMx)

# Load of versions of the semtree package
if ("semtree" %in% (.packages())) {
  detach("package:semtree", unload = TRUE)
}
# Load the modded version
load_all("C:/Users/amori/Desktop/Masterarbeit/semtree_Update/Simulation/semtree_mod")
# Load the old version
load_all("C:/Users/amori/Desktop/Masterarbeit/semtree_Update/Simulation/semtree_old")

# Seed
set.seed(1)

# Define parameters
n <- 500
B <- 5000


# Set up result matrix
res <- matrix(NA, nrow = B, ncol = 2)
colnames(res) <- c("Fix", "Old")
res <- as.data.frame(res)

# Start loop
for (i in 1:B) {
  
  # Simulate data
  long_data <- rmvnorm(n = n, mean = c(0, 0),
                       sigma = matrix(c(1, 0.2, 0.2, 1), nrow = 2, ncol = 2))
  colnames(long_data) <- c("VAEsdq", "aaktiv7")
  long_data <- as.data.frame(long_data)
  
  # Simulate random covaraite -> only 1 continuous
  long_data$z <- rnorm(n)
  
  # Create MCAR pattern
  # Randomly choose ids of missing rows
  miss_idx <- sample.int(n, n/2)
  long_data$VAEsdq[miss_idx] <- NA
  
  
  ######### SEM ##########
  # Define model
  v1 <- var(long_data$VAEsdq,  na.rm = TRUE)
  v2 <- var(long_data$aaktiv7, na.rm = TRUE)
  
  model_mx <- mxModel(
    "covariance_model",
    type = "RAM",
    manifestVars = c("VAEsdq","aaktiv7"),
    # variances: reasonable starts, >0 bounds
    mxPath(from = c("VAEsdq","aaktiv7"), arrows = 2,
           values = c(v1, v2), lbound = 1e-6,
           labels = c("covariance_model_S11", "covariance_model_S22")),
    # covariance: start near 0, keep away from ±1
    mxPath(from = "VAEsdq", to = "aaktiv7", arrows = 2,
           values = 0, lbound = -0.99, ubound = 0.99,
           labels = "covariance_model_S12"),
    # means
    mxPath(from = "one", to = c("VAEsdq","aaktiv7"), arrows = 1,
           labels = c("covariance_model_M11", "covariance_model_M12")),
    mxData(observed = long_data, type = "raw")
  ) 
  
  # Run model
  fit_mx <- mxTryHard(model_mx)
  
  
  
  
  
  ########## SEM Tree modded version ############
  # Define control settings
  control_mod <- semtree_mod::semtree.control(
    method = "score",
    min.N = 50,
    alpha = 0.05,
    max.depth = 3
  )
  
  # Calculate SEM Tree
  sem_tree_mod <- semtree_mod::semtree(
    model = fit_mx,            # SEM model fitted with openMX
    data = long_data,  # Data used for the tree
    control = control_mod, # Control settings
    predictors = c("z")
  )
  
  # Result
  res[i, "Fix"] <- sem_tree_mod$caption == "TERMINAL"
  
  
  
  
  
  ########## SEM Tree old version ############
  # Define control settings
  control_old <- semtree_old::semtree.control(
    method = "score",
    min.N = 50,
    alpha = 0.05,
    max.depth = 3
  )
  
  # Calculate SEM Tree
  sem_tree_old <- semtree_old::semtree(
    model = fit_mx,            # SEM model fitted with openMX
    data = long_data,  # Data used for the tree
    control = control_old, # Control settings
    predictors = c("z")
  )
  
  # Result
  res[i, "Old"] <- sem_tree_old$caption == "TERMINAL"
  
}




# Calculate the Type 1 Error Rate
split_rates <- 1- colMeans(res, na.rm = TRUE)

#  split_rates
# Fix    Old 
# 0.0566 0.0470 

se <- sqrt(split_rates * (1 - split_rates) / B)
ci <- cbind(lo = pmax(0, split_rates - 1.96 * se),
            hi = pmin(1, split_rates + 1.96 * se))
ci
