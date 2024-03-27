# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/R"
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/0326_sfs"
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
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)



# =========================Model for the Allele Frequency Spectrum (AFS)
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters, parallel = TRUE) {
    nSamples <- 100
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
                SFS_model(theta = parameters$theta[i], beta = 0, model_type = 1, n = nSamples)
            }
        )
        stopCluster(cl)
        stats <- rbindlist(stats)
        class(stats) <- "data.frame"
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            stats <- rbind(stats, SFS_model(theta = parameters$theta[i], beta = 0, model_type = 1, n = nSamples))
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
SFS_model <- function(theta, beta, model_type, n) {
    #   This simulates haplotype and site information
    #   Model_type = 1 for regular coalescent (beta = 0)
    #              = 2 for exponential growth coalescent model
    #              = 3 for pure death process with rate 1
    #---Set up arrays and initialisation
    hist <- list() # list of sets of common ancestor labels
    muts <- list() # list of mutation lists
    muts[[2 * n - 1]] <- 0
    mutindx <- 0 # label of last mutation to arise
    #---Generate coalescence times and related vectors
    ctimes <- rep(0, n) # coalescence times set to 0
    stimes <- rep(0, n) # jump times set to 0
    if (model_type != 3) {
        #   Generate ordinary coalescence times
        for (j in n:2) {
            ctimes[j] <- rexp(1, rate = j * (j - 1) / 2)
        }
        if (model_type == 2) {
            #   Calculate times from GT eqn (2.7)
            #   First, compute standard coalescent jump times
            stimes[n] <- ctimes[n]
            for (l in seq(n - 1, 2)) {
                stimes[l] <- ctimes[l] + stimes[l + 1]
            }
            #   Now compute exp growth jump times from GT (2.7)
            for (j in seq(2, n)) {
                stimes[j] <- log(1 + beta * stimes[j]) / beta
            }
            #   Then compute exp growth coal times by differencing
            ctimes[n] <- stimes[n]
            for (l in seq(n - 1, 2)) {
                ctimes[l] <- stimes[l] - stimes[l + 1]
            }
        }
    } else {
        for (j in n:2) {
            ctimes[j] <- rexp(1, rate = j)
        }
    }
    #---Simulate tree topology
    hist[[1]] <- 1:n
    for (j in n:2) { # number of things to choose from
        current <- hist[[n + 1 - j]] # current ancestor labels
        indx <- sample(1:j, 2, replace = F) # which two to merge
        current[indx[1]] <- 2 * n + 1 - j
        current[indx[2]] <- 2 * n + 1 - j
        hist[[n + 2 - j]] <- unique(current)
    }
    #---Pour mutations down the tree
    for (j in n:2) {
        parent <- setdiff(hist[[j]], hist[[j - 1]]) # parent label
        children <- setdiff(hist[[j - 1]], hist[[j]]) # children's labels
        for (l in children) {
            s <- rpois(1, theta * ctimes[n + 2 - j] / 2)
            if (s == 0) {
                muts[[l]] <- muts[[parent]]
            }
            if (s > 0) {
                mutstoadd <- seq(mutindx + 1, mutindx + s)
                mutindx <- mutindx + s # new last mutation label
                muts[[l]] <- union(muts[[parent]], mutstoadd)
            }
        }
        notmerged <- setdiff(hist[[j - 1]], children) # Set of labels w/ no merges at this level
        for (l in notmerged) {
            s <- rpois(1, theta * ctimes[n + 2 - j] / 2) # generate mutations on branch label l
            if (s > 0) {
                mutstoadd <- seq(mutindx + 1, mutindx + s)
                mutindx <- mutindx + s
                muts[[l]] <- union(muts[[l]], mutstoadd)
            }
        }
    }
    # #---Find the haplotype distribution
    # alleles <- table(sapply(muts[1:n], paste, collapse = " "))
    # dimnames(alleles) <- NULL
    # #---Find allele frequency spectrum
    # afs <- tabulate(alleles, nbins = n)
    # nalleles <- sum(afs)
    # ncnt <- sum(1:n * afs)
    #---Find site frequency spectrum
    b1 <- unlist(muts[1:n])
    b2 <- tabulate(b1)
    sfs <- tabulate(b2, nbins = n - 1)
    sval <- sum(sfs)
    # #---Combine table of statistics:
    # #---number of alleles, the homozygoity statistic and the first sqrt(n) frequencies in the ESF
    # ss <- 0
    # for (i in 1:n) ss <- ss + afs[i] * (i / n)^2
    lvec <- floor(sqrt(n))
    if (model_type == 1) {
        # stats <- data.frame(matrix(c(theta, nalleles), nrow = 1))
        # colnames(stats) <- c("theta", "Allele_count_K")
        stats <- data.frame(matrix(c(theta, sval, sfs[1:lvec]), nrow = 1))
        colnames(stats) <- c("theta", "Mutation_count_S", paste0("SFS_", 1:lvec))
        # stats <- data.frame(matrix(c(theta, nalleles, sval, ss, afs[1:lvec]), nrow = 1))
        # colnames(stats) <- c("theta", "Allele_count_K", "Mutation_count_S", "Homozygosity_statistic_F", paste0("AFS_", 1:lvec))
    } else {
        stats <- data.frame(matrix(c(theta, beta, nalleles), nrow = 1))
        colnames(stats) <- c("theta", "beta", "Allele_count_K")
    }
    return(stats)
}
# =====================================================Target statistics
theta <- 10
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
theta <- runif(1000000, 1, 20)
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
    nParticles = rep(200000, 1),
    parallel = TRUE
)
#---Plot posterior marginal distributions against other methods
plots <- plot_compare_marginal(
    # plots = plots,
    parameters_truth = parameters_ground_truth,
    abc_results = abcrf_results,
    parameters_labels = parameters_labels,
    plot_statistics = TRUE
)
# ========================================
#   Plot the out-of-bag estimates (equivalent to cross-validation)
abcrf_results[["Iteration_1"]][["rf_model"]][["model.rf"]]$predictions
abcrf_results$Iteration_1$rf_model$model.rf$predictions
abcrf_results$Iteration_1$parameters$theta
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
png(paste0("NEUTRAL_abcrf_theta_variance_importance.png"), width = 1500, height = 800, res = 150)
n.var <- min(30, length(abcrf_results$Iteration_1$rf_model$model.rf$variable.importance))
imp <- abcrf_results$Iteration_1$rf_model$model.rf$variable.importance
names(imp) <- colnames(statistics_target)
ord <- rev(order(imp, decreasing = TRUE)[1:n.var])
xmin <- 0
xlim <- c(xmin, max(imp) + 1)
dotchart(imp[ord], pch = 19, xlab = "Variable Importance", ylab = "", xlim = xlim, main = NULL, bg = "white", cex = 0.7)
dev.off()
# # =========================================================ABC-rejection
# #---Run ABC-rejection
# abc_rej_results <- abc_rejection(
#     statistics_target = statistics_target,
#     model = model,
#     parameters_labels = parameters_labels,
#     prior_distributions = list(c("unif", 1, 20)),
#     nParticles = 10000, tolerance_quantile = 0.1, progress_bar = TRUE
# )
# #---Plot posterior marginal distributions against other methods
# plots <- plot_compare_marginal(
#     abc_results = abc_rej_results,
#     parameters_labels = parameters_labels,
#     parameters_truth = parameters_ground_truth,
#     plot_statistics = TRUE
# )
# # ==============================================================ABC-MCMC
# #---Run ABC-rejection
# abc_mcmc_results <- abc_mcmc(
#     statistics_target = statistics_target,
#     model = model,
#     parameters_labels = parameters_labels,
#     prior_distributions = list(c("unif", 1, 20)),
#     nParticles = 1000, method = "Marjoram_original", progress_bar = TRUE
# )
# #---Plot posterior marginal distributions against other methods
# plots <- plot_compare_marginal(
#     plots = plots,
#     abc_results = abc_mcmc_results,
#     parameters_labels = parameters_labels,
#     plot_statistics = TRUE
# )
# # ===============================================================ABC-SMC
# #---Find minimum tolerance compatible with noisiness in model
# parameters_test <- do.call(rbind, replicate(1000, parameters_ground_truth, simplify = FALSE))
# statistics_test <- model(parameters = parameters_test)[-c(1:ncol(parameters_test))]
# distance_matrix <- as.matrix(dist(statistics_test, method = "euclidean"))^2
# tolerance_min <- sum(distance_matrix) / (nrow(distance_matrix) * (ncol(distance_matrix) - 1))
# #---Run ABC-SMC
# abc_smc_results <- abc_smc(
#     statistics_target = statistics_target,
#     model = model,
#     parameters_labels = parameters_labels,
#     prior_distributions = list(c("unif", 1, 20)),
#     nParticles = 1000, method = "Beaumont", progress_bar = TRUE,
#     tolerance = c(2 * tolerance_min, 1.5 * tolerance_min, tolerance_min),
#     dist_weights = rep(1, ncol(statistics_target))
# )
# #---Plot posterior marginal distributions against other methods
# plots <- plot_compare_marginal(
#     plots = plots,
#     abc_results = abc_smc_results,
#     parameters_labels = parameters_labels,
#     plot_statistics = TRUE
# )

# # ==========================================SMC-RF for single parameters
# #---Run SMC-RF for single parameters
# smcrf_results <- smcrf(
#     method = "smcrf-single-param",
#     statistics_target = statistics_target,
#     parameters_initial = parameters_initial,
#     model = model,
#     perturb = perturb,
#     range = range,
#     nParticles = rep(1000, 10),
#     parallel = TRUE
# )
# #---Plot marginal distributions
# plot_smcrf_marginal(
#     smcrf_results = smcrf_results,
#     parameters_truth = parameters_ground_truth,
#     parameters_labels = parameters_labels,
#     plot_statistics = TRUE
# )
# #---Plot posterior marginal distributions against other methods
# plots <- plot_compare_marginal(
#     plots = plots,
#     abc_results = smcrf_results,
#     parameters_labels = parameters_labels,
#     plot_statistics = TRUE
# )
