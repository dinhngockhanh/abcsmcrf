# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
# R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/0224_test/hierarchical"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/xiangzijin/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zhihan - Macbook
R_workplace <- "/Users/lexie/Documents/DNA/vignettes/birth-death-mcmc"
R_libPaths <- ""
R_libPaths_extra <- "/Users/lexie/Documents/DNA/SMC-RF/R"
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



set.seed(1)



nNoise <- 20
N <- 1000
# =====================================Model for the birth-death process
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters) {
    library(parallel)
    library(pbapply)
    library(data.table)
    nTimes <- 10
    nNoise <- 0
    #   Make simulations & compute summary statistics (population sizes at each time point)
    cl <- makePSOCKcluster(detectCores() - 1)
    clusterExport(cl, varlist = c("BD_model"))
    stats <- pblapply(
        cl = cl, X = 1:nrow(parameters),
        FUN = function(i) BD_model(parameters$lambda[i], parameters$mu[i], nTimes)
    )
    stopCluster(cl)
    stats <- rbindlist(stats)
    class(stats) <- "data.frame"
    #   Add noise statistics
    noise <- matrix(runif(nrow(parameters) * nNoise), nrow(parameters), nNoise)
    data <- cbind(stats, noise)
    #   Add column names
    data <- data.frame(data)
    if (nNoise > 0) {
        colnames(data) <- c(
            colnames(stats),
            paste0("noise_", c(1:nNoise))
        )
    }
    return(data)
}
BD_model <- function(lambda, mu, nTimes) {
    #   This simulates a birth-death process
    #   Input:  lambda = birth rate
    #           mu = death rate
    times <- seq(1, nTimes) / nTimes
    npt <- length(times)
    alpha <- rep(0, npt)
    beta <- rep(0, npt)
    tdiff <- c(times[1], diff(times))
    Fsim <- rep(0, npt)
    Zsim <- rep(0, npt)
    # compute alpha and beta vectors
    if (lambda == mu) {
        alpha <- lambda * tdiff / (1 + lambda * tdiff)
        beta <- alpha
    } else {
        elm <- exp((lambda - mu) * tdiff)
        alpha <- mu * (elm - 1) / (lambda * elm - mu)
        beta <- lambda * alpha / mu
    }
    # Now generate observations on F and Z. We assume F > 0 at each stage,
    # so we are conditioning on Z > 0 at the end
    zv <- 1 # start from single individual
    Fsim[1] <- sample(1:zv, size = 1, prob = dbinom(1:zv, size = zv, 1 - alpha[1]))
    Zsim[1] <- rnbinom(1, size = Fsim[1], prob = 1 - beta[1]) + Fsim[1]
    for (j in seq(2, npt)) {
        zv <- Zsim[j - 1]
        Fsim[j] <- sample(1:zv, size = 1, prob = dbinom(1:zv, size = zv, 1 - alpha[j]))
        Zsim[j] <- rnbinom(1, size = Fsim[j], prob = 1 - beta[j]) + Fsim[j]
    }
    stats <- data.frame(matrix(c(lambda, mu, Zsim[1:npt]), nrow = 1))
    colnames(stats) <- c("lambda", "mu", paste0("Z_", 1:npt))
    # stats <- data.frame(matrix(c(lambda, mu, Fsim[1:npt], Zsim[1:npt]), nrow = 1))
    # colnames(stats) <- c("lambda", "mu", paste0("F_", 1:npt), paste0("Z_", 1:npt))
    return(stats)
}
# =====================================================Target statistics
parameters_truth <- data.frame(
    lambda = 10,
    mu = 2
)
statistics_target <- model(parameters = parameters_truth)[-c(1:ncol(parameters_truth))]
# ======================================Model for parameter perturbation
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters, after perturbation
perturb <- function(parameters) {
    for (i in 1:ncol(parameters)) parameters[[i]] <- parameters[[i]] + runif(nrow(parameters), min = -0.5, max = 0.5)
    return(parameters)
}
# ======================================Define ranges for the parameters
range <- data.frame(
    parameter = c("lambda", "mu"),
    min = c(0, 0),
    max = c(15, 15)
)
# ========================================Initial guesses for parameters
# ====================================(sampled from prior distributions)
lambda <- runif(10000, 0, 15)
mu <- runif(10000, 0, lambda)
parameters_initial <- data.frame(
    lambda = lambda,
    mu = mu
)
# ====================================Labels for parameters in the plots
parameters_labels <- data.frame(
    parameter = c("lambda", "mu"),
    label = c(deparse(expression(lambda)), deparse(expression(mu)))
)
# ==========================================SMC-RF for single parameters
#---Run SMC-RF for single parameters
smcrf_results_single_param <- smcrf(
    method = "smcrf-single-param",
    statistics_target = statistics_target,
    parameters_initial = parameters_initial,
    model = model,
    perturb = perturb,
    range = range,
    nParticles = rep(10000, 7),
    parallel = TRUE
)
#---Plot marginal distributions
plot_smcrf_marginal(
    smcrf_results = smcrf_results_single_param,
    parameters_labels = parameters_labels,
    plot_statistics = FALSE
)
#---Plot joint distributions
plot_smcrf_joint(
    smcrf_results = smcrf_results_single_param,
    parameters_labels = parameters_labels
)
# ========================================SMC-RF for multiple parameters
#---Run SMC-RF for multiple parameters
smcrf_results_multi_param <- smcrf(
    method = "smcrf-multi-param",
    statistics_target = statistics_target,
    parameters_initial = parameters_initial,
    model = model,
    perturb = perturb,
    range = range,
    nParticles = rep(10000, 7),
    parallel = TRUE
)
#---Plot marginal distributions
plot_smcrf_marginal(
    smcrf_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels,
    plot_statistics = FALSE
)
#---Plot joint distributions
plot_smcrf_joint(
    smcrf_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels
)


# =========Compare final posterior distributions among different methods
plots <- plot_compare_marginal(
    abc_results = smcrf_results_single_param,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels,
    plot_statistics = TRUE
)
plots <- plot_compare_marginal(
    plots = plots,
    abc_results = smcrf_results_multi_param,
    plot_statistics = TRUE
)





# =======================================================ABC-MCMC
library(EasyABC)
toy_prior = list(c("unif", 0, 15), c("unif", 0, 15))
toy_model <- function(parameters) {
    library(parallel)
    library(pbapply)
    library(data.table)
    nTimes <- 10
    nNoise <- 0
    stats <- BD_model(parameters[1], parameters[2], nTimes)
    # stats <- rbindlist(stats)
    class(stats) <- "data.frame"
    #   Add noise statistics
    noise <- matrix(runif(nNoise), 1, nNoise)
    data <- cbind(stats, noise)
    #   Add column names
    data <- data.frame(data)
    if (nNoise > 0) {
        colnames(data) <- c(
            colnames(stats),
            paste0("noise_", c(1:nNoise))
        )
    }
    data <- data[-c(1:length(parameters))]
    data <- unlist(data)
    return(data)
}
# ABC_rejection(model=model,prior=prior,nb_simul=3, prior_test="X1 > X2")
sum_stat_obs <- unlist(statistics_target)
ABC_Marjoram<-ABC_mcmc(
    method="Marjoram", model=toy_model, prior=toy_prior,
    prior_test="X1 > X2",
    n_rec = 70000,
    use_seed = FALSE, 
    verbose = FALSE,
    progress_bar = TRUE,
    summary_stat_target=sum_stat_obs)

#---Set up color scheme for plotting
color_scheme <- c(
    "True Posterior" = "black",
    "SMC-RF for single parameters" = "salmon",
    "SMC-RF for multiple parameters" = "royalblue4",
    "ABC-MCMC" = "#74403b"
)
#---Set up legend order for plotting
legend_order <- c(
    "True Posterior",
    "SMC-RF for single parameters",
    "SMC-RF for multiple parameters",
    "ABC-MCMC"
)
alpha = 0.3
for (i in 1:dim(ABC_Marjoram$stats)[2]) {
    posterior_df <- data.frame(
        value = ABC_Marjoram$stats[,i],
        legend = "ABC-MCMC"
        )
    p <- ggplot() + 
        geom_density(data = posterior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2) + 
        scale_fill_manual(values = color_scheme, name = "", breaks = legend_order) +
        scale_color_manual(values = color_scheme, name = "", breaks = legend_order) +
        guides(fill = guide_legend(nrow = 1, keywidth = 2.5, keyheight = 1)) +
        theme(
            text = element_text(size = 50),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5)
        )
    file_name <- paste0("marginal-Z_", i, ".png")
    png(file_name, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
    print(p)
    dev.off()
}


for (i in 1:length(parameters_labels$parameter)) {
    parameter_id <- parameters_labels$parameter[i]
    #   Begin plot
    p <- ggplot()
    #   Plot true posterior distribution (if provided)
    if (!is.null(parameters_truth)) {
        true_posterior_df <- data.frame(value = parameters_truth[[parameter_id]], legend = "True Posterior")
        p <- p + geom_density(data = true_posterior_df, aes(x = value, fill = legend, color = legend), linewidth = 2)
    }
    # #   Plot prior distribution
    # prior_df <- data.frame(value = smcrf_results[["Iteration_1"]]$parameters[[parameter_id]], legend = "Prior Distribution")
    # p <- p + geom_density(data = prior_df, aes(x = value, fill = legend, color = legend), alpha = 0.2, linewidth = 2)
    #   Plot posterior distribution for each iteration
    posterior_df <- data.frame(
        value = ABC_Marjoram$param[,i],
        legend = "ABC-MCMC"
    )
    p <- p + geom_density(data = posterior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
    #   Add label for parameter
    if ("label" %in% colnames(parameters_labels)) {
        p <- p + labs(x = eval(parse(text = parameters_labels$label[i])))
    } else {
        p <- p + labs(x = parameter_id)
    }
    #   Beautify plot
    p <- p +
        scale_fill_manual(values = color_scheme, name = "", breaks = legend_order) +
        scale_color_manual(values = color_scheme, name = "", breaks = legend_order) +
        guides(fill = guide_legend(nrow = 1, keywidth = 2.5, keyheight = 1)) +
        theme(
            text = element_text(size = 50),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.position = "top",
            legend.justification = c(0, 0.5)
        )
    #   Print marginal distribution plot
    file_name <- paste0("MCMC-marginal-", parameter_id, ".png")
    png(file_name, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
    print(p)
    dev.off()
}



posterior_df <- data.frame(
    x = ABC_Marjoram$param[,1],
    y = ABC_Marjoram$param[,2],
    legend = "ABC-MCMC"
)
p <- ggplot() + geom_density_2d_filled(data = posterior_df, aes(x = x, y = y), show.legend = TRUE)
if ("label" %in% colnames(parameters_labels)) {
        p <- p + labs(x = eval(parse(text = parameters_labels$label[1])), y = eval(parse(text = parameters_labels$label[2])))
    } else {
        p <- p + labs(x = parameter_id[1], y = parameter_id[2])
    }
p <- p +
    guides(color = guide_legend(nrow = 1, keywidth = 5, override.aes = list(lwd = 10))) +
    theme(
        text = element_text(size = 50),
        panel.background = element_rect(fill = "white", colour = "white"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"),
        legend.position = "top",
        legend.justification = c(0, 0.5)
    )
#---Print joint distribution plot
file_name <- paste0("mcmc-joint-", parameters_labels$parameter[1], "-", parameters_labels$parameter[2], ".png")
png(file_name, res = 150, width = 30, height = 31, units = "in", pointsize = 12)
print(p)