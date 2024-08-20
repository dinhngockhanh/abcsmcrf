library(abcsmcrf)
library(invgamma)
# ===================================Function for plotting extrame cases
plot_hierarchical_extreme <- function(
    drf_results,
    parameters_truth,
    limits) {
    library(ggplot2)
    library(patchwork)
    # ========================================Plot multidimensional contours
    Y_drf <- drf_results[["Iteration_2"]]$parameters
    indx.sim <- sample(1:nrow(Y_drf), size = 500, replace = T)
    Y_drf <- Y_drf[indx.sim, ]
    Y_truth <- parameters_truth
    # ==========================TOP plot = marginal distribution for theta_1
    p_top <- ggplot() +
        geom_histogram(
            data = Y_drf,
            aes(x = theta1, y = ..density..),
            fill = "royalblue2", color = "royalblue2", alpha = 0.5
        ) +
        geom_histogram(
            data = Y_truth, aes(x = theta1, y = ..density..),
            fill = "firebrick1", color = "firebrick1", alpha = 0.5
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
            data = Y_truth, aes(x = theta2, y = ..density..),
            fill = "firebrick1", color = "firebrick1", alpha = 0.5
        ) +
        theme_void() +
        coord_flip() +
        xlim(c(limits$min[which(limits$parameter == "theta2")], limits$max[which(limits$parameter == "theta2")]))

    # ===============Main plot = joint distribution from DRF, TRUE and SMCRF
    p_main <- ggplot() +
        geom_density_2d_filled(data = Y_drf, aes(x = theta1, y = theta2)) +
        geom_density_2d(data = Y_truth, aes(x = theta1, y = theta2), colour = "firebrick1", linewidth = 2, bins = 4) +
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
        "mean", "variance", "mad",
        "sum_mean_var",
        "sum_mean_mad",
        "sum_var_mad",
        "sum_mean_var_mad",
        "prod_mean_var",
        "prod_mean_mad",
        "prod_var_mad",
        "prod_mean_var_mad",
        paste0("noise_", c(1:nNoise))
    )
    return(data)
}
# =====================================================Target statistics
alpha <- 4
beta <- 5
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
ybar <- statistics_target[, "mean"]
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
bounds <- data.frame(
    parameter = c("theta1", "theta2"),
    min = c(-Inf, 0),
    max = c(Inf, Inf)
)
# =================================Define limits for plotting parameters
limits <- data.frame(
    parameter = c("theta1", "theta2"),
    min = c(0, 0),
    max = c(5, 8)
)
# ========================================Initial guesses for parameters
# ====================================(sampled from prior distributions)
theta2 <- 1 / rgamma(100000, shape = alpha, rate = beta)
theta1 <- rnorm(100000, 0, sqrt(theta2))
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
drf_results <- smcrf(
    method = "smcrf-multi-param",
    statistics_target = statistics_target,
    parameters_initial = parameters_initial,
    model = model,
    perturb = perturb,
    bounds = bounds,
    nParticles = rep(20000, 1),
    num.trees = 2500,
    # compute.variable.importance = TRUE,
    parallel = TRUE
)
# =============================================Plot marginal joint plots
plot_hierarchical_extreme(
    drf_results = drf_results,
    limits = limits,
    parameters_truth = parameters_truth
)
# # ==========================================Plot the Variable Importance
# png(paste0("NORMAL_drf_CART_variable_importance.png"), width = 800, height = 800, res = 150)
# n.var <- 30
# variable_importance <- drf_results[["Iteration_1"]]$rf_model$variable.importance
# names(imp) <- c(
#     "mean", "variance", "mad",
#     "mean+var",
#     "mean+mad",
#     "var+mad",
#     "mean+var+mad",
#     "mean*var",
#     "mean*mad",
#     "var*mad",
#     "mean*var*mad",
#     paste0("noise_", c(1:50))
# )
# ord <- rev(order(variable_importance, decreasing = TRUE)[1:n.var])
# xlim <- c(0, max(variable_importance) + 1)
# dotchart(variable_importance[ord], pch = 19, xlab = "Variable Importance", ylab = "", xlim = xlim, main = NULL, bg = "white", cex = 0.7)
# dev.off()
