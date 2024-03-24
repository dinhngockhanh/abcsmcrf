# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/vignettes"
R_libPaths <- ""
R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
# R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/0224_test/hierarchical"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/xiangzijin/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zhihan - Macbook
# R_workplace <- "/Users/lexie/Documents/DNA/SMC-RF/vignettes/birth_death"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/lexie/Documents/DNA/SMC-RF/R"
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



# =====================================Model for the birth-death process
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters, parallel = TRUE) {
    nTimes <- 10
    nNoise <- 0
    if (exists("nSimulations")) nSimulations <<- nSimulations + nrow(parameters)
    #   Make simulations & compute summary statistics (population sizes at each time point)
    if (parallel) {
        library(parallel)
        library(pbapply)
        library(data.table)
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cl, varlist = c("BD_model"))
        stats <- pblapply(
            cl = cl, X = 1:nrow(parameters),
            FUN = function(i) BD_model(parameters$lambda[i], parameters$mu[i], nTimes)
        )
        stopCluster(cl)
        stats <- rbindlist(stats)
        class(stats) <- "data.frame"
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            stats <- rbind(stats, BD_model(parameters$lambda[i], parameters$mu[i], nTimes))
        }
    }
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
statistics_target <- model(parameters = parameters_truth, parallel = FALSE)[-c(1:ncol(parameters_truth))]
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
# ==============================================================ABC-MCMC
#---Run ABC-MCMC
abc_mcmc_results <- abc_mcmc(
    statistics_target = statistics_target,
    model = model,
    parameters_labels = parameters_labels,
    prior_distributions = list(c("unif", 0, 15), c("unif", 0, 15)),
    prior_test = "X1 > X2",
    nParticles = 1000, method = "Marjoram_original", progress_bar = TRUE
)
#---Plot posterior joint distributions against other methods
plots <- plot_compare_joint(
    abc_results = abc_mcmc_results,
    parameters_labels = parameters_labels
)
# ===================================================================DRF
#---Run ABC-RF
abcrf_results <- smcrf(
    method = "smcrf-multi-param",
    statistics_target = statistics_target,
    parameters_initial = parameters_initial,
    model = model,
    perturb = perturb,
    range = range,
    nParticles = rep(10000, 1),
    parallel = TRUE
)
#---Plot posterior joint distributions against other methods
plots <- plot_compare_joint(
    plots = plots,
    abc_results = abcrf_results,
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
#---Plot joint distributions
plot_smcrf_joint(
    smcrf_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels
)
#---Plot posterior joint distributions against other methods
plots <- plot_compare_joint(
    plots = plots,
    abc_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels
)
