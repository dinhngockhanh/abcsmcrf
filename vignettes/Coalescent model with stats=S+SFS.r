# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/0329_sfs_for_paper/coala_npop=1000_nsim=10000/new_results/10000sim;npop=1000;abcrf&abc-rej;onlyS"
# R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/0329_sfs_for_paper/coala_npop=1000_nsim=10000/new_results/10000sim npop=1000;abcrf&abc-rej"
# R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/0329_sfs_for_paper/coala_npop=1000_nsim=10000/new_results/10000sim;npop=1000;abcrf;onlySFS"
R_libPaths <- ""
R_libPaths_extra <- "/Users/xiangzijin/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macmini
# R_workplace <- "/Users/khanhngocdinh/Documents/Zijin/0328_sfs_coala"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/khanhngocdinh/Documents/Zijin/SMC-RF/R"
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
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)
# =========================Model for the Site Frequency Spectrum (SFS)
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters, parallel = TRUE) {
    nSamples <- 1000
    nNoise <- 0
    if (exists("nSimulations")) nSimulations <<- nSimulations + nrow(parameters)
    #   Make simulations & compute summary statistics (allele count)
    if (parallel) {
        library(parallel)
        library(pbapply)
        library(data.table)
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cl, varlist = c("SFS_model"))
        stats <- pblapply(
            cl = cl, X = 1:nrow(parameters),
            FUN = function(i) {
                SFS_model(theta = parameters$theta[i], n = nSamples)
            }
        )
        stopCluster(cl)
        stats <- rbindlist(stats)
        class(stats) <- "data.frame"
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            stats <- rbind(stats, SFS_model(theta = parameters$theta[i], n = nSamples))
        }
    }
    #   Add noise statistics
    noise <- matrix(runif(nrow(parameters) * nNoise), nrow(parameters), nNoise)
    #   Add column names
    data <- data.frame(cbind(stats, noise))
    if (nNoise > 0) {
        colnames(data) <- c(
            colnames(stats),
            paste0("noise_", c(1:nNoise))
        )
    }
    return(data)
}

SFS_model <- function(theta, n) {
    library(coala)
    model <- coal_model(n, 0) +
        locus_single(1) +
        feat_mutation(par_const(theta)) +
        sumstat_sfs()
    sim_data <- simulate(model, nsim = 1)
    sfs <- create_abc_sumstat(sim_data, model)
    indices <- seq_along(sfs)
    weights <- 1 / (match(sfs, sfs)^2)
    weighted_sfs <- sfs * weights
    colnames(weighted_sfs) <- colnames(sfs)
    sval <- sum(sfs)
    sval_weighted <- sum(weighted_sfs)
    mean_sfs <- mean(sfs)
    lvec <- floor(sqrt(n))
    stats <- data.frame(matrix(c(theta, sval, sfs[1:lvec]), nrow = 1))
    colnames(stats) <- c("theta", "Mutation_count_S", paste0("SFS_", 1:lvec))
    return(stats)
}
# =====================================================Target statistics
set.seed(1)
theta <- runif(1, 1, 20)
parameters_target <- data.frame(
    theta = theta
)
statistics_target <- model(parameters = parameters_target, parallel = FALSE)[-c(1:ncol(parameters_target))]
# ======================================Model for parameter perturbation
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters, after perturbation
perturb <- function(parameters) {
    for (i in 1:ncol(parameters)) parameters[[i]] <- parameters[[i]] + runif(nrow(parameters), min = -1, max = 1)
    return(parameters)
}
# ======================================Define ranges for the parameters
range <- data.frame(
    parameter = c("theta"),
    min = c(1),
    max = c(20)
)
# ========================================Initial guesses for parameters
# ====================================(sampled from prior distributions)
set.seed(1)
theta <- runif(10000, 1, 20)
parameters_initial <- data.frame(
    theta = theta
)
# ====================================Labels for parameters in the plots
parameters_labels <- data.frame(
    parameter = c("theta"),
    label = c(deparse(expression(theta)))
)
# ========================================================True posterior
theta_single <- function(sample_size, a, b, s, theta) {
    n <- sample_size
    ln <- sum(1 / 1:(n - 1))
    cons <- (pgamma(a * ln, s + 1, 1, lower = FALSE) - pgamma(b * ln, s + 1, 1, lower = FALSE)) / ln # could use log(n) for ln
    print(cons)
    dens <- dpois(s, theta * ln) / cons
    integrand <- function(x, s, n) {
        dpois(s, x * sum(1 / 1:(n - 1)))
    }
    val <- integrate(integrand, a, b, s, n)
    print(val)
    return(dens)
}
density <- c()
test_theta <- sample(parameters_initial$theta, 1000)
for (theta in test_theta) {
    density <- c(density, theta_single(1000, 0, 20, statistics_target$Mutation_count_S, theta))
}
parameters_truth <- data.frame(
    theta = test_theta,
    density = density
)
# ================================================================ABC-RF
#---Run ABC-RF
abcrf_results <- smcrf(
    method = "smcrf-single-param",
    statistics_target = statistics_target,
    parameters_initial = parameters_initial,
    model = model,
    perturb = perturb,
    range = range,
    nParticles = rep(10000, 1),
    parallel = TRUE
)
#---Plot posterior marginal distributions against other methods
plots <- plot_compare_marginal(
    parameters_truth = parameters_truth,
    abc_results = abcrf_results,
    parameters_labels = parameters_labels,
    plot_statistics = TRUE,
    xlimit = range,
    plot_hist = TRUE,
    plot_prior = TRUE
)
