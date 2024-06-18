# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/hierarchical/new_extreme_cases/0.5;1"
R_libPaths <- ""
R_libPaths_extra <- "/Users/xiangzijin/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zhihan - Macbook
# R_workplace <- "/Users/lexie/Documents/DNA/hierarchical/0.5;3"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/lexie/Documents/DNA/SMC-RF/R"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(ggplot2)
library(gridExtra)
library(grid)
library(invgamma)
library(patchwork)
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)
set.seed(1)
# ===================================Function for plotting extrame cases
plot_hierarchical_extreme <- function(
    smcrf_results,
    drf_results,
    abcrf_results,
    parameters_truth,
    limits) {
    library(ggplot2)
    library(patchwork)
    # ========================================Plot multidimensional contours
    nIterations <- smcrf_results[["nIterations"]]
    Y_drf <- drf_results[["Iteration_2"]]$parameters
    Y_abcrf <- abcrf_results[["Iteration_2"]]$parameters
    Y_smcrf <- smcrf_results[[paste0("Iteration_", nIterations + 1)]]$parameters
    Y_truth <- parameters_truth
    # ==========================TOP plot = marginal distribution for theta_1
    p_top <- ggplot() +
        geom_histogram(
            data = Y_drf,
            aes(x = theta1, y = ..density..),
            fill = "royalblue2", color = "royalblue2", alpha = 0.5
        ) +
        geom_histogram(
            data = Y_abcrf, aes(x = theta1, y = ..density..),
            fill = "magenta4", color = "magenta4", alpha = 0.5
        ) +
        geom_histogram(
            data = Y_smcrf, aes(x = theta1, y = ..density..),
            fill = "salmon", color = "salmon", alpha = 0.5
        ) +
        geom_histogram(
            data = Y_truth, aes(x = theta1, y = ..density..),
            fill = "goldenrod2", color = "goldenrod2", alpha = 0.5
        ) +
        theme_void() +
        xlim(c(limits$min[which(limits$parameter == "theta1")], limits$max[which(limits$parameter == "theta1")]))

    # ========================Right plot = marginal distribution for theta_2
    p_right <- ggplot() +
        geom_histogram(
            data = Y_drf,
            aes(x = theta2, y = ..density..),
            fill = "royalblue2", color = "royalblue2", alpha = 0.5
        ) +
        geom_histogram(
            data = Y_abcrf, aes(x = theta2, y = ..density..),
            fill = "magenta4", color = "magenta4", alpha = 0.5
        ) +
        geom_histogram(
            data = Y_smcrf, aes(x = theta2, y = ..density..),
            fill = "salmon", color = "salmon", alpha = 0.5
        ) +
        geom_histogram(
            data = Y_truth, aes(x = theta2, y = ..density..),
            fill = "goldenrod2", color = "goldenrod2", alpha = 0.5
        ) +
        theme_void() +
        coord_flip() +
        xlim(c(limits$min[which(limits$parameter == "theta2")], limits$max[which(limits$parameter == "theta2")]))

    # ===============Main plot = joint distribution from DRF, TRUE and SMCRF
    p_main <- ggplot() +
        geom_density_2d_filled(data = Y_truth, aes(x = theta1, y = theta2)) +
        geom_density_2d(data = Y_drf, aes(x = theta1, y = theta2), colour = "cyan1", linewidth = 2, bins = 4) +
        geom_density_2d(data = Y_smcrf, aes(x = theta1, y = theta2), colour = "firebrick2", linewidth = 2, bins = 4) +
        theme_minimal() +
        xlim(c(limits$min[which(limits$parameter == "theta1")], limits$max[which(limits$parameter == "theta1")])) +
        ylim(c(limits$min[which(limits$parameter == "theta2")], limits$max[which(limits$parameter == "theta2")])) +
        labs(x = expression(theta[1]), y = expression(theta[2])) +
        theme(
            text = element_text(size = 30),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none"
        )

    # Combine the plots using patchwork
    combined_plot <- p_top + plot_spacer() + p_main + p_right +
        plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
    # Save the plot to a file
    file_name <- paste0("NORMAL_drf_contour_with_densities_theta1_theta2.png")
    ggsave(file_name, plot = combined_plot, width = 17, height = 17, units = "in", dpi = 150)
}
# ==============================Model for the hierarchical model example
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters) {
    theta1 <- parameters$theta1
    theta2 <- parameters$theta2
    nSamples <- 10
    nNoise <- 0
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
        "prod_esp_var_mad"
        # paste0("noise_", c(1:nNoise))
    )
    return(data)
}

# =====================================================Target statistics
alpha <- 1
beta <- 1
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
# perturb <- function(parameters) {
#     for (i in 1:ncol(parameters)) parameters[[i]] <- parameters[[i]] + runif(nrow(parameters), min = -0.1, max = 0.1)
#     return(parameters)
# }
perturb <- function(parameters) {
    if (any(grepl("theta1", colnames(parameters)))) {
        parameters[["theta1"]] <- parameters[["theta1"]] + runif(nrow(parameters), min = -0.1, max = 0.1)
    } else if (any(grepl("theta2", colnames(parameters)))) {
        parameters[["theta2"]] <- parameters[["theta2"]] + runif(nrow(parameters), min = -0.1, max = 0.1)
    }
    return(parameters)
}

# ======================================Define ranges for the parameters
range <- data.frame(
    parameter = c("theta1", "theta2"),
    min = c(-Inf, 0),
    max = c(Inf, Inf)
)
#--limits for 0.5;3
# para_limits <- data.frame(
#     parameter = c("theta1", "theta2"),
#     min = c(-5, 0),
#     max = c(10, 50)
# )
#--limits for 1.5;3
# para_limits <- data.frame(
#     parameter = c("theta1", "theta2"),
#     min = c(-3, 0),
#     max = c(10, 100)
# )
#--limits for 1;1
para_limits <- data.frame(
    parameter = c("theta1", "theta2"),
    min = c(-3, 0),
    max = c(10, 40)
)
# #--limits for 0.5;1
# para_limits <- data.frame(
#     parameter = c("theta1", "theta2"),
#     min = c(-5, 0),
#     max = c(10, 30)
# )
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
# drf_results <- smcrf(
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
# save(drf_results, file = "drf.rda")
load("drf.rda")
#---Plot posterior joint distributions against other methods
plots_joint <- plot_compare_joint(
    abc_results = drf_results,
    parameters_labels = parameters_labels,
    # lims = limits,
    lims = para_limits,
    parameters_truth = parameters_truth
)
#---Plot marginal distributions compare
plots_marginal <- plot_compare_marginal(
    # plots = plots_marginal,
    abc_results = drf_results,
    xlimit = para_limits,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)

# ========================================SMC-RF for multiple parameters
#---Run SMC-RF for multiple parameters
# smcrf_results_multi_param <- smcrf(
#     method = "smcrf-multi-param",
#     statistics_target = statistics_target,
#     parameters_initial = parameters_initial,
#     model = model,
#     perturb = perturb,
#     range = range,
#     num.trees = 2500,
#     nParticles = rep(5000, 2),
#     parallel = TRUE
# )
# save(smcrf_results_multi_param, file = "smc-drf.rda")
load("smc-drf.rda")
# #---Plot joint distributions
plot_smcrf_joint(
    smcrf_results = smcrf_results_multi_param,
    parameters_truth = parameters_truth,
    lims = para_limits,
    parameters_labels = parameters_labels
)
#---Plot posterior joint distributions against other methods
plots_joint <- plot_compare_joint(
    plots = plots_joint,
    abc_results = smcrf_results_multi_param,
    lims = para_limits,
    parameters_labels = parameters_labels
)
#---Plot marginal distributions compare
plots_marginal <- plot_compare_marginal(
    plots = plots_marginal,
    abc_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels,
    xlimit = para_limits,
    parameters_truth = parameters_truth,
    plot_hist = TRUE
)
plot_smcrf_marginal(
    smcrf_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels,
    xlimit = para_limits,
    parameters_truth = parameters_truth,
    plot_hist = TRUE
)
# ========================================ABC-RF for single parameters
#---Run SMC-RF for single parameters
# abcrf_results <- smcrf(
#     method = "smcrf-single-param",
#     statistics_target = statistics_target,
#     parameters_initial = parameters_initial,
#     model = model,
#     perturb = perturb,
#     range = range,
#     ntree = 2500,
#     nParticles = rep(10000, 1),
#     parallel = TRUE
# )
# save(abcrf_results, file = "abcrf.rda")
load("abcrf.rda")
#---Plot marginal distributions compare
plots_marginal <- plot_compare_marginal(
    plots = plots_marginal,
    abc_results = abcrf_results,
    parameters_labels = parameters_labels,
    xlimit = para_limits,
    parameters_truth = parameters_truth,
    plot_hist = TRUE
)
# =============================================Plot marginal joint plots
plot_hierarchical_extreme(
    smcrf_results = smcrf_results_multi_param,
    drf_results = drf_results,
    abcrf_results = abcrf_results,
    limits = para_limits,
    parameters_truth = parameters_truth
)
library(abcrf)
for (para in colnames(parameters_initial)) {
    #   Train the ABCrf model
    RFmodel <- abcrf_results$Iteration_1$rf_model[[para]]
    #   Plot the out-of-bag estimates (equivalent to cross-validation)
    png(paste0("NORMAL_abcrf_", para, "_out_of_bag.png"))
    plot(abcrf_results$Iteration_1$reference[, para],
        RFmodel$model.rf$predictions,
        xlab = "True value",
        ylab = "Out-of-bag estimate"
    ) + abline(a = 0, b = 1, col = "red")
    dev.off()
    #   Can the error be lowered by increasing the number of trees?
    ref_table <- cbind(abcrf_results$Iteration_1$reference[, para], abcrf_results$Iteration_1$reference[-c(1:ncol(parameters_initial))])
    colnames(ref_table)[1] <- "para"
    oob_error <- err.regAbcrf(RFmodel, training = ref_table, paral = TRUE)
    png(paste0("NORMAL_abcrf_", para, "_error_by_ntree.png"))
    plot(oob_error[, "ntree"], oob_error[, "oob_mse"], type = "l", xlab = "Number of trees", ylab = "Out-of-bag MSE")
    dev.off()
    #   Variable Importance of each statistic in inferring theta
    png(paste0("NORMAL_abcrf_", para, "_variable_importance.png"))
    plot(x = RFmodel)
    dev.off()
}
