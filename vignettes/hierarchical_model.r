# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/0224_test/hierarchical"
R_libPaths <- ""
R_libPaths_extra <- "/Users/xiangzijin/SMC-RF/R"
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
# =============================================SET UP INITIAL PARAMETERS
N <- 1000
n_samples_per_parameter_set <- 20
nNoise <- 50
alpha <- 4
beta <- 3
# ======================================================================
set.seed(1)
#---Model for the hierarchical model example
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters,
                  n_samples_per_parameter_set,
                  nNoise) {
    theta1 <- parameters$theta1
    theta2 <- parameters$theta2
    n_samples_per_parameter_set <- n_samples_per_parameter_set # Number of samples per parameter set
    nNoise <- nNoise # Number of noise variables
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
# ======================================Model for parameter perturbation
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters, after perturbation
# perturb <- function(parameters) {
#     parameters$theta1 <- parameters$theta1 + runif(nrow(parameters), min = -1, max = 1)
#     parameters$theta2 <- parameters$theta2 + runif(nrow(parameters), min = -1, max = 1)
#     # parameters$theta2 <- pmax(parameters$theta2 + runif(nrow(parameters), min = -0.5, max = 0.5), 0)
#     return(parameters)
# }
perturb <- function(parameter) {
    # parameter <- parameter + runif(1, min = -0.5, max = 0.5)
    parameter <- parameter + runif(1, min = -1, max = 1)
    # parameters$theta2 <- pmax(parameters$theta2 + runif(nrow(parameters), min = -0.5, max = 0.5), 0)
    return(parameter)
}
range <- data.frame(
    parameter = c("theta1", "theta2"),
    min = c(-Inf, 0),
    max = c(Inf, Inf)
)
theta1_min <- range$min[which(range$parameter == "theta1")]
theta1_max <- range$max[which(range$parameter == "theta1")]

# =====================================================Target statistics
theta2 <- 1 / rgamma(1, shape = alpha, rate = beta)
theta1 <- rnorm(1, 0, sqrt(theta2))
parameters_target <- data.frame(
    theta1 = theta1,
    theta2 = theta2
)
target <- model(
    parameters = parameters_target,
    n_samples_per_parameter_set = n_samples_per_parameter_set,
    nNoise = nNoise
)[-c(1:ncol(parameters_target))]

# ========================================Initial guesses for parameters
# ====================================(sampled from prior distributions)
theta2 <- 1 / rgamma(N, shape = alpha, rate = beta)
theta1 <- rnorm(N, 0, sqrt(theta2))
parameters_initial <- data.frame(
    theta1 = theta1,
    theta2 = theta2
)
# ===============================True marginal posteriors for parameters
# ================(sampled from prior distributions)-hierarchical normal
s_2 <- target[, "variance"] * (n_samples_per_parameter_set - 1)
ybar <- target[, "expectation"]
mean_theta1 <- n_samples_per_parameter_set / (n_samples_per_parameter_set + 1) * ybar
scale_theta1 <- sqrt((2 * (3 + s_2 / 2 + n_samples_per_parameter_set * ybar^2 / (2 * n_samples_per_parameter_set + 2))) / ((n_samples_per_parameter_set + 1) * (n_samples_per_parameter_set + 8)))
t_deviate <- rt(N * 10, df = n_samples_per_parameter_set + 8)
theta1_true <- mean_theta1 + t_deviate * scale_theta1
shape_theta2 <- n_samples_per_parameter_set / 2 + 4
scale_theta2 <- 0.5 * (s_2 + 6 + n_samples_per_parameter_set * ybar^2 / (n_samples_per_parameter_set + 1))
theta2_true <- rinvgamma(N * 10, shape = shape_theta2, scale = 1 / scale_theta2)
parameters_truth <- data.frame(
    theta1 = theta1_true,
    theta2 = theta2_true
)
# ========================================================Run SMC-ABCDRF
drf_output <- smcdrf(
    target = target,
    model = model,
    n_samples_per_parameter_set = n_samples_per_parameter_set,
    nNoise = nNoise,
    perturb = perturb,
    range = range,
    parameters_initial = parameters_initial,
    nIter = 7, # Number of iterations
    nParticles = rep(1000, 7), # Number of particles for each iteration
    # ntree = 2000,
    parallel = T,
)
filename <- "ABCSMC_DRF_output.rda"
save(drf_output, file = filename)
# =========================================================Run SMC-ABCRF
rf_output <- smcabcrf(
    target = target,
    model = model,
    n_samples_per_parameter_set = n_samples_per_parameter_set,
    nNoise = nNoise,
    perturb = perturb,
    parameters_initial = parameters_initial,
    nIter = 7, # Number of iterations
    range = range,
    nParticles = rep(N, 7), # Number of particles for each iteration
    # ntree = 2000,
    parallel = T
)
filename <- "ABCSMC_RF_output.rda"
save(rf_output, file = filename)
# =========================================Plot SMC-ABCRF for Parameters
plotting_smcrf(
    parameters_truth = parameters_truth,
    parameters_initial = parameters_initial,
    parameters_id = colnames(parameters_initial),
    outputdata = drf_output
)
plotting_smcrf(
    parameters_truth = parameters_truth,
    parameters_initial = parameters_initial,
    parameters_id = colnames(parameters_initial),
    outputdata = rf_output
)
# ==============================================Plot SMC-ABCRF for Stats
plotting_smcrf(
    parameters_id = names(target),
    outputdata = drf_output,
    Plot_stats = TRUE
)
plotting_smcrf(
    parameters_id = names(target),
    outputdata = rf_output,
    Plot_stats = TRUE
)
