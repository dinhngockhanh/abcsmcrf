# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/adaptive"
R_libPaths <- ""
R_libPaths_extra <- "/Users/xiangzijin/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zhihan - Macbook
# R_workplace <- "/Users/lexie/Documents/DNA/SMC-RF/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/lexie/Documents/DNA/SMC-RF/R"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(ggplot2)
library(gridExtra)
library(grid)
library(invgamma)
library(SimBIID)
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)
set.seed(1)
# ==============================Model for the hierarchical model example
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters) {
    theta1 <- parameters$theta1
    theta2 <- parameters$theta2
    nSamples <- 10
    nNoise <- 50
    #   Make simulations
    y <- matrix(NA, length(theta1), nSamples)
    for (i in 1:length(theta1)) {
        y[i, ] <- rnorm(nSamples, theta1[i], sqrt(theta2[i]))
    }
    #   Compute some summary statistics (mean, variance, MAD)
    summary <- matrix(NA, length(
        theta1
    ), 3)
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
# =====================================================Target statistics
alpha <- 4
beta <- 3
theta2 <- 1 / rgamma(1, shape = alpha, rate = beta)
theta1 <- rnorm(1, 0, sqrt(theta2))
parameters_target <- data.frame(
    theta1 = theta1,
    theta2 = theta2
)
statistics_target <- model(parameters = parameters_target)[-c(1:ncol(parameters_target))]
# ===============================True marginal posteriors for parameters
nSamples <- 10
s_2 <- statistics_target[, "variance"] * (nSamples - 1)
ybar <- statistics_target[, "expectation"]
shape_theta2 <- nSamples / 2 + alpha
scale_theta2 <- 0.5 * (s_2 + 2 * beta + nSamples * ybar^2 / (nSamples + 1))
theta2_true <- rinvgamma(10000, shape = shape_theta2, scale = 1 / scale_theta2)
theta1_true <- rnorm(10000, mean = nSamples * ybar / (nSamples + 1), sd = sqrt(2 * theta2_true / (nSamples + 1))) # removed in sd a 2 here
parameters_truth <- data.frame(
    theta1 = theta1_true,
    theta2 = theta2_true
)
# ======================================Model for parameter perturbation
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters, after perturbation
perturb <- function(parameters) {
    for (i in 1:ncol(parameters)) parameters[[i]] <- parameters[[i]] + runif(nrow(parameters), min = -0.1, max = 0.1)
    return(parameters)
}
# ======================================Define ranges for the parameters
range <- data.frame(
    parameter = c("theta1", "theta2"),
    min = c(-Inf, 0),
    max = c(Inf, Inf)
)
# ========================================Initial guesses for parameters
# ====================================(sampled from prior distributions)
theta2 <- 1 / rgamma(10000, shape = alpha, rate = beta)
theta1 <- rnorm(10000, 0, sqrt(theta2))
parameters_initial <- data.frame(
    theta1 = theta1,
    theta2 = theta2
)
# ====================================Labels for parameters in the plots
parameters_labels <- data.frame(
    parameter = c("theta1", "theta2"),
    label = c(deparse(expression(theta[1])), deparse(expression(theta[2])))
)
# ===================================================================DRF
#---Run ABC-DRF
# abcrf_results <- smcrf(
#     method = "smcrf-multi-param",
#     statistics_target = statistics_target,
#     parameters_initial = parameters_initial,
#     model = model,
#     perturb = perturb,
#     range = range,
#     nParticles = rep(10000, 1),
#     num.trees = 2500,
#     parallel = TRUE
# )
# save(abcrf_results, file = "drf.rda")
# drf_results <- load("drf.rda")
#---Define limits for the parameters
limits <- data.frame(
    parameter = c("theta1", "theta2"),
    min = c(0, 0),
    max = c(3, 3)
)
#---Plot posterior joint distributions against other methods
# plots_joint <- plot_compare_joint(
#     abc_results = abcrf_results,
#     parameters_labels = parameters_labels,
#     parameters_truth = parameters_truth,
#     lims = limits
# )
# #---Plot marginal distributions compare
# plots_marginal <- plot_compare_marginal(
#     # plots = plots_marginal,
#     abc_results = abcrf_results,
#     parameters_truth = parameters_truth,
#     parameters_labels = parameters_labels,
#     plot_hist = TRUE
# )

# ========================================SMC-RF for multiple parameters
#---Run SMC-RF for multiple parameters
smcrf_results_multi_param <- smcrf(
    method = "smcrf-multi-param",
    statistics_target = statistics_target,
    parameters_initial = parameters_initial,
    model = model,
    # perturb = perturb,
    perturb = "Beaumont",
    range = range,
    # num.trees = 2500,
    nParticles = rep(500, 5),
    parallel = TRUE
)
save(smcrf_results_multi_param, file = "smc-drf.rda")
load("smc-drf.rda")
#---Plot joint distributions
plot_smcrf_joint(
    smcrf_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels,
    lims = limits
)
#---Plot posterior joint distributions against other methods
plots_joint <- plot_compare_joint(
    # plots = plots_joint,
    abc_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels,
    lims = limits
)
#---Plot marginal distributions compare
plots_marginal <- plot_compare_marginal(
    # plots = plots_marginal,
    abc_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels,
    parameters_truth = parameters_truth,
    plot_hist = TRUE
)
plot_smcrf_marginal(
    smcrf_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels,
    parameters_truth = parameters_truth,
    plot_hist = TRUE
)

# ========================================SMC-RF for multiple parameters
#---Run SMC-RF for multiple parameters
smcrf_results_single_param <- smcrf(
    method = "smcrf-single-param",
    statistics_target = statistics_target,
    parameters_initial = parameters_initial,
    model = model,
    # perturb = perturb,
    perturb = "Beaumont",
    range = range,
    # num.trees = 2500,
    nParticles = rep(1000, 5),
    parallel = TRUE
)
save(smcrf_results_single_param, file = "smc-rf.rda")
load("smc-rf.rda")
#---Plot marginal distributions compare
plots_marginal <- plot_compare_marginal(
    plots = plots_marginal,
    abc_results = smcrf_results_single_param,
    parameters_labels = parameters_labels,
    parameters_truth = parameters_truth,
    plot_hist = TRUE
)
plot_smcrf_marginal(
    smcrf_results = smcrf_results_single_param,
    parameters_labels = parameters_labels,
    parameters_truth = parameters_truth,
    plot_hist = TRUE
)
