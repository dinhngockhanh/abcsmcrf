# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/0329_sfs_for_paper/coala_npop=1000_nsim=10000/new_results/10000sim npop=1000;abcrf&abc-rej"
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
    # stats <- data.frame(matrix(c(theta, sval, sfs[1:lvec]), nrow = 1))
    # colnames(stats) <- c("theta", "Mutation_count_S", paste0("SFS_", 1:lvec))
    stats <- data.frame(matrix(c(theta, sval, weighted_sfs[1:lvec]), nrow = 1))
    colnames(stats) <- c("theta", "Mutation_count_S", paste0("SFS_", 1:lvec))
    return(stats)
}
# =====================================================Target statistics
set.seed(1)
theta <- runif(1, 1, 20)
parameters_ground_truth <- data.frame(
    theta = theta
)
statistics_target <- model(parameters = parameters_ground_truth, parallel = FALSE)[-c(1:ncol(parameters_ground_truth))]
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
save(abcrf_results, file = "abcrf_results.rda")

load("/Users/xiangzijin/Documents/ABC_SMCRF/0329_sfs_for_paper/coala_npop=1000_nsim=10000/new_results/10000sim npop=1000;abcrf&abc-rej/abcrf_results.rda")
#---Plot posterior marginal distributions against other methods
plots <- plot_compare_marginal(
    # plots = plots,
    parameters_truth = parameters_ground_truth,
    abc_results = abcrf_results,
    parameters_labels = parameters_labels,
    plot_statistics = TRUE,
    xlimit = c(0, 20),
    plot_hist = TRUE,
    plot_prior = TRUE
)
# ========================================
#   Plot the out-of-bag estimates (equivalent to cross-validation)
png(paste0("NEUTRAL_abcrf_theta_out_of_bag.png"))
plot(abcrf_results$Iteration_1$parameters$theta,
    abcrf_results$Iteration_1$rf_model$model.rf$predictions,
    xlab = "True value",
    ylab = "Out-of-bag estimate"
) + abline(a = 0, b = 1, col = "red")
dev.off()
#   Can the error be lowered by increasing the number of trees?
library(abcrf)
oob_error <- err.regAbcrf(abcrf_results$Iteration_1$rf_model, training = abcrf_results$Iteration_1$reference, paral = T)
png(paste0("NEUTRAL_abcrf_theta_error_by_ntree.png"))
plot(oob_error[, "ntree"], oob_error[, "oob_mse"], type = "l", xlab = "Number of trees", ylab = "Out-of-bag MSE")
dev.off()
#   Variance Importance of each statistic in inferring gamma
png(paste0("NEUTRAL_abcrf_theta_variable_importance.png"), width = 1500, height = 800, res = 150)
n.var <- min(30, length(abcrf_results$Iteration_1$rf_model$model.rf$variable.importance))
imp <- abcrf_results$Iteration_1$rf_model$model.rf$variable.importance
names(imp) <- colnames(statistics_target)
ord <- rev(order(imp, decreasing = TRUE)[1:n.var])
xmin <- 0
xlim <- c(xmin, max(imp) + 1)
dotchart(imp[ord], pch = 19, xlab = "Variable Importance", ylab = "", xlim = xlim, main = NULL, bg = "white", cex = 0.7)
dev.off()
# =========================================================ABC-Rejection
#---Run ABC
abc_rej_results <- abc_rejection(
    statistics_target = statistics_target,
    model = model,
    parameters_labels = parameters_labels,
    prior_distributions = list(c("unif", 0, 20)),
    tolerance_quantile = 0.05,
    nParticles = 10000, progress_bar = TRUE
)


save(abc_rej_results, file = "abcrej_results.rda")
#---Plot marginal distributions compare
plots_marginal <- plot_compare_marginal(
    plots = plots,
    abc_results = abc_rej_results,
    parameters_labels = parameters_labels,
    parameters_truth = parameters_ground_truth,
    plot_statistics = TRUE,
    # xlimit = c(0, 20),
    plot_hist = TRUE,
    plot_prior = FALSE
)


# # ==========================================SMC-RF for single parameters
# #---Run SMC-RF for single parameters
# smcrf_results <- smcrf(
#     method = "smcrf-single-param",
#     statistics_target = statistics_target,
#     parameters_initial = parameters_initial,
#     model = model,
#     perturb = perturb,
#     range = range,
#     nParticles = rep(2000, 5),
#     parallel = TRUE
# )
# #---Plot marginal distributions
# plot_smcrf_marginal(
#     smcrf_results = smcrf_results,
#     parameters_truth = parameters_ground_truth,
#     parameters_labels = parameters_labels,
#     plot_statistics = TRUE,
#     plot_hist = TRUE,
#     bin_counts = 55
# )
# #---Plot posterior marginal distributions against other methods
# plots <- plot_compare_marginal(
#     plots = plots,
#     abc_results = smcrf_results,
#     parameters_labels = parameters_labels,
#     plot_statistics = TRUE,
#     xlimit = c(0, 20),
#     plot_hist = TRUE,
#     bin_counts = 55,
#     plot_prior = TRUE
# )

# # =============================================ABC-package-ABC-Rejection
# #---Run ABC
# load("/Users/xiangzijin/Documents/ABC_SMCRF/0329_sfs_for_paper/coala_npop=1000_nsim=10000/new_results/10000sim npop=1000;abcrf&abc-rej/abcrf_results.rda")
# abc_params <- abcrf_results[["Iteration_1"]][["reference"]][, 1]
# abc_stats <- abcrf_results[["Iteration_1"]][["reference"]][, c(3:12)]
# library(abc)
# posterior <- abc(statistics_target[2:11], abc_params, abc_stats, 0.05, method = "rejection")
# png(paste0("ABC_package_test_10_SFS.png"))
# hist(posterior, breaks = 20)
# dev.off()

# library(coala)
sfs_test <- as.matrix(c(112, 57, 24, 34, 16, 29, 8, 10, 15), nrow = 1)
sfs_test
new_sfs_test <- sapply(sfs_test, function(x) x / (which(sfs_test == x))^2)
new_sfs_test
# test_model <- coal_model(10, 50) +
#     feat_mutation(par_prior("theta", runif(1, 1, 5))) +
#     sumstat_sfs()
# sim_data <- simulate(test_model, nsim = 2000, seed = 17)
# sim_param <- create_abc_param(sim_data, test_model)
# sim_sumstat <- create_abc_sumstat(sim_data, test_model)
# posterior <- abc(sfs_test, sim_param, sim_sumstat, 0.05, method = "rejection")
# hist(posterior, breaks = 20)
