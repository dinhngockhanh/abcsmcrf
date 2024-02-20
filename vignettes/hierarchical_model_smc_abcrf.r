# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/R"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)

library(ggplot2)
library(gridExtra)
library(grid)
library(invgamma)

setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)







# ======================================================================
set.seed(1)
#---Model for the hierarchical model example
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters) {
    theta1 <- parameters$theta1
    theta2 <- parameters$theta2
    n_samples_per_parameter_set <- 20 # Number of samples per parameter set
    nNoise <- 500 # Number of noise variables
    #   Make simulations
    y <- matrix(NA, length(theta1), n_samples_per_parameter_set)
    for (i in 1:length(theta1)) {
        y[i, ] <- rnorm(n_samples_per_parameter_set, theta1[i], sqrt(theta2[i]))
    }
    #   Compute some summary statistics (mean, variance, MAD)
    summary <- matrix(NA, length(theta1), 3)
    for (i in 1:length(theta1)) {
        summary[i, ] <- c(mean(y[i, ]), var(y[i, ]), mad(y[i, ]))
    }
    data <- cbind(parameters, summary)
    #   Add some others summary statistics
    x <- data[, -c(1:2)]
    x <- cbind(
        x[, 1] + x[, 2], x[, 1] + x[, 3], x[, 2] + x[, 3], x[, 1] + x[, 2] + x[, 3],
        x[, 1] * x[, 2], x[, 1] * x[, 3], x[, 2] * x[, 3], x[, 1] * x[, 2] * x[, 3]
    )
    data <- cbind(data, x)
    #   Add noise statistics
    noise <- matrix(runif(nrow(parameters) * nNoise), nrow(parameters), nNoise)
    data <- cbind(data, noise)
    #   Add column names
    data <- data.frame(data)
    colnames(data) <- c(
        "theta1", "theta2",
        "expectation", "variance", "mad",
        "sum_esp_var",
        "sum_esp_mad",
        "sum_var_mad",
        "sum_esp_var_mad",
        "prod_esp_var",
        "prod_esp_mad",
        "prod_var_mad",
        "prod_esp_var_mad",
        paste0("noise_", c(1:nNoise))
    )
    return(data)
}
#---Model for parameter perturbation
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters, after perturbation
perturb <- function(parameters) {
    parameters$theta1 <- parameters$theta1 + runif(nrow(parameters), min = -1, max = 1)
    parameters$theta2 <- pmax(parameters$theta2 + runif(nrow(parameters), min = -1, max = 1), 0)
    return(parameters)
}
#---Target statistics
theta2 <- 1 / rgamma(1, shape = 4, rate = 3)
theta1 <- rnorm(1, 0, sqrt(theta2))
target <- model(data.frame(theta1 = theta1, theta2 = theta2))[-c(1:2)]
#---Initial guesses for parameters (sampled from prior distributions)
parameters_initial <- data.frame(
    theta1 = runif(1000, 0, 10),
    theta2 = runif(1000, 0, 10)
)
#---Run SMC-ABCRF
smcabcrf_test(
    target = target,
    model = model,
    perturb = perturb,
    parameters_initial = parameters_initial,
    nIter = 7, # Number of iterations
    nKeep = rep(1000, 7), # Number of particles for each iteration
    # ntree = 2000,
    parallel = T
)
# smcabcrf(
#     target = target,
#     model = model,
#     perturb = perturb,
#     parameters_initial = parameters_initial,
#     nIter = 7, # Number of iterations
#     nKeep = rep(1000, 7), # Number of particles for each iteration
#     # ntree = 2000,
#     parallel = T
# )
