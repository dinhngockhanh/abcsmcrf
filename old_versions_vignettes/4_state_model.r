# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ZIJIN - Macbook
R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/CME/4-state/11params/RF_nRuns=10;nsims=3000"
R_libPaths <- ""
R_libPaths_extra <- "/Users/xiangzijin/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - HPC
# R_workplace <- getwd()
# R_libPaths <- "/burg/iicd/users/zx2406/rpackages"
# R_libPaths_extra <- "/burg/iicd/users/zx2406/R_smcrf"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(ggplot2)
library(gridExtra)
library(grid)
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)
# ======================================================== Gillespie SSA
SSA <- function(initial_state, parameters, reaction_propensities, reaction_stoichiometries, time_points) {
    stopifnot(all(diff(time_points) > 0))
    state <- initial_state
    time <- 0
    state_out <- matrix(0, nrow = length(time_points), ncol = length(initial_state))
    next_time_point_index <- 1
    tmp <- 0
    while (time < time_points[length(time_points)]) {
        tmp <- tmp + 1
        #   Compute reaction propensities
        props <- reaction_propensities(state, parameters, time)
        #   Compute total propensity
        total_prop <- sum(props)
        if (total_prop == 0) {
            delta_time <- Inf
        } else {
            #   Compute time to next reaction & reaction to occur
            delta_time <- rexp(1, rate = total_prop)
            reaction <- sample.int(length(props), size = 1, prob = props / total_prop)
            #   Update state
            state <- mapply(function(x, y) x + y, state, reaction_stoichiometries[[reaction]], SIMPLIFY = FALSE)
        }
        #   Update time
        time_next <- time + delta_time
        if (time < time_points[next_time_point_index] && time_next >= time_points[next_time_point_index]) {
            state_out[next_time_point_index, ] <- unlist(state)
            next_time_point_index <- next_time_point_index + 1
        }
        time <- time_next
    }
    return(state_out)
}
# =======================================Model for the chemical reaction
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model_target_statistics <- function(parameters, parallel = TRUE) {
    nNoise <- 0
    if (exists("nSimulations")) {
        nSimulations <<- nSimulations + nrow(parameters)
    }
    if (parallel) {
        library(parallel)
        library(pbapply)
        library(data.table)
        library(ggplot2)
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cl, varlist = c("SSA", "mRNA_model", "mRNA_model_multi_runs", "Hog1p_values", "Hog1p_time_points_sec"))
        stats <- pblapply(
            cl = cl, X = 1:nrow(parameters),
            FUN = function(i) {
                # mRNA_model(
                #     10^(parameters$log_k12[i]),
                #     10^(parameters$log_k21a[i]),
                #     10^(parameters$log_k21b[i]),
                #     10^(parameters$log_k23[i]),
                #     10^(parameters$log_k32[i]),
                #     10^(parameters$log_k34[i]),
                #     10^(parameters$log_k43[i]),
                #     10^(parameters$log_kr2[i]),
                #     10^(parameters$log_kr3[i]),
                #     10^(parameters$log_kr4[i]),
                #     10^(parameters$log_gamma[i])
                # )
                mRNA_model_multi_runs(
                    10^(parameters$log_k12[i]),
                    10^(parameters$log_k21a[i]),
                    10^(parameters$log_k21b[i]),
                    10^(parameters$log_k23[i]),
                    10^(parameters$log_k32[i]),
                    10^(parameters$log_k34[i]),
                    10^(parameters$log_k43[i]),
                    10^(parameters$log_kr2[i]),
                    10^(parameters$log_kr3[i]),
                    10^(parameters$log_kr4[i]),
                    10^(parameters$log_gamma[i]),
                    nRuns = 1000
                )
            }
        )
        stopCluster(cl)
        stats <- rbindlist(stats)
        class(stats) <- "data.frame"
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            stats <- rbind(stats, mRNA_model_multi_runs(
                10^(parameters$log_k12[i]),
                10^(parameters$log_k21a[i]),
                10^(parameters$log_k21b[i]),
                10^(parameters$log_k23[i]),
                10^(parameters$log_k32[i]),
                10^(parameters$log_k34[i]),
                10^(parameters$log_k43[i]),
                10^(parameters$log_kr2[i]),
                10^(parameters$log_kr3[i]),
                10^(parameters$log_kr4[i]),
                10^(parameters$log_gamma[i]),
                # nRuns = 1000
                nRuns = 1000
            ))
        }
    }
    new_stats <- data.frame(
        log_k12 = log10(stats$k12),
        log_k21a = log10(stats$k21a),
        log_k21b = log10(stats$k21b),
        log_k23 = log10(stats$k23),
        log_k32 = log10(stats$k32),
        log_k34 = log10(stats$k34),
        log_k43 = log10(stats$k43),
        log_kr2 = log10(stats$kr2),
        log_kr3 = log10(stats$kr3),
        log_kr4 = log10(stats$kr4),
        log_gamma = log10(stats$gamma)
    )
    stats <- cbind(new_stats, stats[, -c(1:11)])
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
model_inference <- function(parameters, parallel = TRUE) {
    nNoise <- 0
    if (exists("nSimulations")) {
        nSimulations <<- nSimulations + nrow(parameters)
    }
    if (parallel) {
        library(parallel)
        library(pbapply)
        library(data.table)
        library(ggplot2)
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cl, varlist = c("SSA", "mRNA_model", "mRNA_model_multi_runs", "Hog1p_values", "Hog1p_time_points_sec"))
        stats <- pblapply(
            cl = cl, X = 1:nrow(parameters),
            FUN = function(i) {
                mRNA_model_multi_runs(
                    10^(parameters$log_k12[i]),
                    10^(parameters$log_k21a[i]),
                    10^(parameters$log_k21b[i]),
                    10^(parameters$log_k23[i]),
                    10^(parameters$log_k32[i]),
                    10^(parameters$log_k34[i]),
                    10^(parameters$log_k43[i]),
                    10^(parameters$log_kr2[i]),
                    10^(parameters$log_kr3[i]),
                    10^(parameters$log_kr4[i]),
                    10^(parameters$log_gamma[i]),
                    nRuns = 100
                )
            }
        )
        stopCluster(cl)
        stats <- rbindlist(stats)
        class(stats) <- "data.frame"
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            stats <- rbind(stats, mRNA_model_multi_runs(
                10^(parameters$log_k12[i]),
                10^(parameters$log_k21a[i]),
                10^(parameters$log_k21b[i]),
                10^(parameters$log_k23[i]),
                10^(parameters$log_k32[i]),
                10^(parameters$log_k34[i]),
                10^(parameters$log_k43[i]),
                10^(parameters$log_kr2[i]),
                10^(parameters$log_kr3[i]),
                10^(parameters$log_kr4[i]),
                10^(parameters$log_gamma[i]),
                nRuns = 100
            ))
        }
    }
    new_stats <- data.frame(
        log_k12 = log10(stats$k12),
        log_k21a = log10(stats$k21a),
        log_k21b = log10(stats$k21b),
        log_k23 = log10(stats$k23),
        log_k32 = log10(stats$k32),
        log_k34 = log10(stats$k34),
        log_k43 = log10(stats$k43),
        log_kr2 = log10(stats$kr2),
        log_kr3 = log10(stats$kr3),
        log_kr4 = log10(stats$kr4),
        log_gamma = log10(stats$gamma)
    )
    stats <- cbind(new_stats, stats[, -c(1:11)])
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
model <- function(parameters, parallel = TRUE) {
    nNoise <- 0
    if (exists("nSimulations")) {
        nSimulations <<- nSimulations + nrow(parameters)
    }
    if (parallel) {
        library(parallel)
        library(pbapply)
        library(data.table)
        library(ggplot2)
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cl, varlist = c("SSA", "mRNA_model", "mRNA_model_multi_runs", "Hog1p_values", "Hog1p_time_points_sec"))
        stats <- pblapply(
            cl = cl, X = 1:nrow(parameters),
            FUN = function(i) {
                mRNA_model(
                    10^(parameters$log_k12[i]),
                    10^(parameters$log_k21a[i]),
                    10^(parameters$log_k21b[i]),
                    10^(parameters$log_k23[i]),
                    10^(parameters$log_k32[i]),
                    10^(parameters$log_k34[i]),
                    10^(parameters$log_k43[i]),
                    10^(parameters$log_kr2[i]),
                    10^(parameters$log_kr3[i]),
                    10^(parameters$log_kr4[i]),
                    10^(parameters$log_gamma[i])
                )
            }
        )
        stopCluster(cl)
        stats <- rbindlist(stats)
        class(stats) <- "data.frame"
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            stats <- rbind(stats, mRNA_model(
                10^(parameters$log_k12[i]),
                10^(parameters$log_k21a[i]),
                10^(parameters$log_k21b[i]),
                10^(parameters$log_k23[i]),
                10^(parameters$log_k32[i]),
                10^(parameters$log_k34[i]),
                10^(parameters$log_k43[i]),
                10^(parameters$log_kr2[i]),
                10^(parameters$log_kr3[i]),
                10^(parameters$log_kr4[i]),
                10^(parameters$log_gamma[i])
            ))
        }
    }
    new_stats <- data.frame(
        log_k12 = log10(stats$k12),
        log_k21a = log10(stats$k21a),
        log_k21b = log10(stats$k21b),
        log_k23 = log10(stats$k23),
        log_k32 = log10(stats$k32),
        log_k34 = log10(stats$k34),
        log_k43 = log10(stats$k43),
        log_kr2 = log10(stats$kr2),
        log_kr3 = log10(stats$kr3),
        log_kr4 = log10(stats$kr4),
        log_gamma = log10(stats$gamma)
    )
    stats <- cbind(new_stats, stats[, -c(1:11)])
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
mRNA_model <- function(k12, k21_a, k21_b, k23, k32, k34, k43, kr2, kr3, kr4, gamma) {
    time_points_min <- c(4, 6, 10, 20, 30, 40)
    time_points_sec <- time_points_min * 60
    initial_state <- list(
        OFF = 1,
        S2 = 0,
        S3 = 0,
        S4 = 0,
        mRNA = 0
    )
    parameters <- list(
        k12 = k12,
        k21_a = k21_a,
        k21_b = k21_b,
        k23 = k23,
        k32 = k32,
        k34 = k34,
        k43 = k43,
        kr2 = kr2,
        kr3 = kr3,
        kr4 = kr4,
        gamma = gamma
    )
    reaction_propensities <- function(state, parameters, time) {
        c(
            parameters$k12 * state$OFF,
            max(0, parameters$k21_a - parameters$k21_b * Hog1p_values[which.min(abs(Hog1p_time_points_sec - time))]) * state$S2,
            parameters$k23 * state$S2,
            parameters$k32 * state$S3,
            parameters$k34 * state$S3,
            parameters$k43 * state$S4,
            parameters$kr2 * state$S2,
            parameters$kr3 * state$S3,
            parameters$kr4 * state$S4,
            parameters$gamma * state$mRNA
        )
    }
    reaction_stoichiometries <- list(
        list(OFF = -1, S2 = 1, S3 = 0, S4 = 0, mRNA = 0),
        list(OFF = 1, S2 = -1, S3 = 0, S4 = 0, mRNA = 0),
        list(OFF = 0, S2 = -1, S3 = 1, S4 = 0, mRNA = 0),
        list(OFF = 0, S2 = 1, S3 = -1, S4 = 0, mRNA = 0),
        list(OFF = 0, S2 = 0, S3 = -1, S4 = 1, mRNA = 0),
        list(OFF = 0, S2 = 0, S3 = 1, S4 = -1, mRNA = 0),
        list(OFF = 0, S2 = 0, S3 = 0, S4 = 0, mRNA = 1),
        list(OFF = 0, S2 = 0, S3 = 0, S4 = 0, mRNA = 1),
        list(OFF = 0, S2 = 0, S3 = 0, S4 = 0, mRNA = 1),
        list(OFF = 0, S2 = 0, S3 = 0, S4 = 0, mRNA = -1)
    )
    tmp <- SSA(
        initial_state = initial_state,
        parameters = parameters,
        reaction_propensities = reaction_propensities,
        reaction_stoichiometries = reaction_stoichiometries,
        time_points = time_points_sec
    )
    stats <- data.frame(matrix(c(k12, k21_a, k21_b, k23, k32, k34, k43, kr2, kr3, kr4, gamma, tmp[, 5]), nrow = 1))
    # colnames(stats) <- c("k12", "k21a", "k21b", "k23", "k32", "k34", "k43", "kr2", "kr3", "kr4", "gamma", paste0("mRNA_", time_points_min, "min"))
    colnames(stats) <- c("log_k12", "log_k21a", "log_k21b", "log_k23", "log_k32", "log_k34", "log_k43", "log_kr2", "log_kr3", "log_kr4", "log_gamma", paste0("mRNA_", time_points_min, "min"))
    return(stats)
}
mRNA_model_multi_runs <- function(k12, k21_a, k21_b, k23, k32, k34, k43, kr2, kr3, kr4, gamma, nRuns) {
    time_points_min <- c(4, 6, 10, 20, 30, 40)
    time_points_sec <- time_points_min * 60
    initial_state <- list(
        OFF = 1,
        S2 = 0,
        S3 = 0,
        S4 = 0,
        mRNA = 0
    )
    parameters <- list(
        k12 = k12,
        k21_a = k21_a,
        k21_b = k21_b,
        k23 = k23,
        k32 = k32,
        k34 = k34,
        k43 = k43,
        kr2 = kr2,
        kr3 = kr3,
        kr4 = kr4,
        gamma = gamma
    )
    reaction_propensities <- function(state, parameters, time) {
        c(
            parameters$k12 * state$OFF,
            max(0, parameters$k21_a - parameters$k21_b * Hog1p_values[which.min(abs(Hog1p_time_points_sec - time))]) * state$S2,
            parameters$k23 * state$S2,
            parameters$k32 * state$S3,
            parameters$k34 * state$S3,
            parameters$k43 * state$S4,
            parameters$kr2 * state$S2,
            parameters$kr3 * state$S3,
            parameters$kr4 * state$S4,
            parameters$gamma * state$mRNA
        )
    }
    reaction_stoichiometries <- list(
        list(OFF = -1, S2 = 1, S3 = 0, S4 = 0, mRNA = 0),
        list(OFF = 1, S2 = -1, S3 = 0, S4 = 0, mRNA = 0),
        list(OFF = 0, S2 = -1, S3 = 1, S4 = 0, mRNA = 0),
        list(OFF = 0, S2 = 1, S3 = -1, S4 = 0, mRNA = 0),
        list(OFF = 0, S2 = 0, S3 = -1, S4 = 1, mRNA = 0),
        list(OFF = 0, S2 = 0, S3 = 1, S4 = -1, mRNA = 0),
        list(OFF = 0, S2 = 0, S3 = 0, S4 = 0, mRNA = 1),
        list(OFF = 0, S2 = 0, S3 = 0, S4 = 0, mRNA = 1),
        list(OFF = 0, S2 = 0, S3 = 0, S4 = 0, mRNA = 1),
        list(OFF = 0, S2 = 0, S3 = 0, S4 = 0, mRNA = -1)
    )
    results <- matrix(NA, nrow = nRuns, ncol = length(time_points_min))
    for (i in 1:nRuns) {
        tmp <- SSA(
            initial_state = initial_state,
            parameters = parameters,
            reaction_propensities = reaction_propensities,
            reaction_stoichiometries = reaction_stoichiometries,
            time_points = time_points_sec
        )
        results[i, ] <- tmp[, 5]
    }
    # stats <- data.frame(matrix(c(k12, k21_a, k21_b, k23, k32, k34, k43, kr2, kr3, kr4, gamma, c(apply(results, 2, mean))), nrow = 1))
    # colnames(stats) <- c("k12", "k21a", "k21b", "k23", "k32", "k34", "k43", "kr2", "kr3", "kr4", "gamma", paste0("mRNA_", time_points_min, "min_mean"))
    stats <- data.frame(matrix(c(k12, k21_a, k21_b, k23, k32, k34, k43, kr2, kr3, kr4, gamma, c(apply(results, 2, mean)), c(apply(results, 2, var))), nrow = 1))
    colnames(stats) <- c("k12", "k21a", "k21b", "k23", "k32", "k34", "k43", "kr2", "kr3", "kr4", "gamma", paste0("mRNA_", time_points_min, "min_mean"), paste0("mRNA_", time_points_min, "min_var"))
    return(stats)
}

input_signal <- function(t,
                         r1,
                         r2,
                         eta,
                         A,
                         M) {
    part1 <- (1 - exp(-r1 * t)) * exp(-r2 * t)
    part2 <- 1 + ((1 - exp(-r1 * t)) * exp(-r2 * t)) / M
    hog1p_star <- A * (part1 / part2)^eta
    return(hog1p_star)
}
r1 <- 6.9e-5
r2 <- 3.6e-3
eta <- 3.1
A <- 9.3e9
M <- 6.4e-4
Hog1p_time_points_min <- seq(0, 40, 0.001)
Hog1p_time_points_sec <- Hog1p_time_points_min * 60
Hog1p_values <- input_signal(t = Hog1p_time_points_sec, r1 = r1, r2 = r2, eta = eta, A = A, M = M)
# =====================================================Target statistics
parameters_truth <- data.frame(
    log_k12 = log10(1.29),
    log_k21a = log10(3200),
    log_k21b = log10(7710),
    log_k23 = log10(0.0191),
    log_k32 = log10(0.0175),
    log_k34 = log10(0.133),
    log_k43 = log10(0.0083),
    log_kr2 = log10(0.0098),
    log_kr3 = log10(1.01),
    log_kr4 = log10(0.0016),
    log_gamma = log10(0.0020)
)
statistics_target <- model_target_statistics(parameters = parameters_truth, parallel = FALSE)[-c(1:ncol(parameters_truth))]
print(statistics_target)
save(statistics_target, file = "statistics_target.rda")
load("statistics_target.rda")
# ======================================Model for parameter perturbation
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters, after perturbation
perturb <- function(parameters) {
    parameters[["log_k12"]] <- parameters[["log_k12"]] + runif(nrow(parameters), min = -1.0, max = 1.0)
    parameters[["log_k21a"]] <- parameters[["log_k21a"]] + runif(nrow(parameters), min = -1.0, max = 1.0)
    parameters[["log_k21b"]] <- parameters[["log_k21b"]] + runif(nrow(parameters), min = -1.0, max = 1.0)
    parameters[["log_k23"]] <- parameters[["log_k23"]] + runif(nrow(parameters), min = -1.0, max = 1.0)
    parameters[["log_k32"]] <- parameters[["log_k32"]] + runif(nrow(parameters), min = -1.0, max = 1.0)
    parameters[["log_k34"]] <- parameters[["log_k34"]] + runif(nrow(parameters), min = -1.0, max = 1.0)
    parameters[["log_k43"]] <- parameters[["log_k43"]] + runif(nrow(parameters), min = -1.0, max = 1.0)
    parameters[["log_kr2"]] <- parameters[["log_kr2"]] + runif(nrow(parameters), min = -1.0, max = 1.0)
    parameters[["log_kr3"]] <- parameters[["log_kr3"]] + runif(nrow(parameters), min = -1.0, max = 1.0)
    parameters[["log_kr4"]] <- parameters[["log_kr4"]] + runif(nrow(parameters), min = -1.0, max = 1.0)
    parameters[["log_gamma"]] <- parameters[["log_gamma"]] + runif(nrow(parameters), min = -1.0, max = 1.0)
    return(parameters)
}
# ======================================Define ranges for the parameters
range <- data.frame(
    parameter = c("log_k12", "log_k21a", "log_k21b", "log_k23", "log_k32", "log_k34", "log_k43", "log_kr2", "log_kr3", "log_kr4", "log_gamma"),
    min = c(-4, 0, 0, -4, -4, -4, -4, -4, -4, -4, -4),
    max = c(2, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2)
)
# ========================================Initial guesses for parameters
# ====================================(sampled from prior distributions)
log_k12 <- runif(100000, -4, 2)
log_k21a <- runif(100000, 0, 4)
log_k21b <- runif(100000, 0, 4)
log_k23 <- runif(100000, -4, 2)
log_k32 <- runif(100000, -4, 2)
log_k34 <- runif(100000, -4, 2)
log_k43 <- runif(100000, -4, 2)
log_kr2 <- runif(100000, -4, 2)
log_kr3 <- runif(100000, -4, 2)
log_kr4 <- runif(100000, -4, 2)
log_gamma <- runif(100000, -4, 2)
parameters_initial <- data.frame(
    log_k12 = log_k12,
    log_k21a = log_k21a,
    log_k21b = log_k21b,
    log_k23 = log_k23,
    log_k32 = log_k32,
    log_k34 = log_k34,
    log_k43 = log_k43,
    log_kr2 = log_kr2,
    log_kr3 = log_kr3,
    log_kr4 = log_kr4,
    log_gamma = log_gamma
)
# ====================================Labels for parameters in the plots
parameters_labels <- data.frame(
    parameter = c("log_k12", "log_k21a", "log_k21b", "log_k23", "log_k32", "log_k34", "log_k43", "log_kr2", "log_kr3", "log_kr4", "log_gamma"),
    label = c(
        "expression(k[12])", "expression(k[21*a])", "expression(k[21*b])", "expression(k[23])", "expression(k[32])",
        "expression(k[34])", "expression(k[43])", "expression(k[r2])", "expression(k[r3])", "expression(k[r4])", "expression(gamma)"
    )
)
# ===================================================================DRF
#---Run DRF
abcrf_results <- smcrf(
    method = "smcrf-multi-param",
    statistics_target = statistics_target,
    parameters_initial = parameters_initial,
    model = model_inference,
    perturb = perturb,
    range = range,
    # num.trees = 2500,
    nParticles = rep(1000, 1), # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # nParticles = rep(20, 1), # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    parallel = TRUE
)
save(abcrf_results, file = "munsky_drf.rda")
drf_results <- load("munsky_drf.rda")
#---Plot marginal distributions compare
plots_marginal <- plot_compare_marginal(
    # plots = plots_marginal,
    xlimit = range,
    abc_results = abcrf_results,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)
# # ========================================SMC-RF for multiple parameters
# #---Run SMC-RF for multiple parameters
# smcrf_results_multi_param <- smcrf(
#     method = "smcrf-multi-param",
#     statistics_target = statistics_target,
#     parameters_initial = parameters_initial,
#     model = model,
#     perturb = perturb,
#     range = range,
#     nParticles = rep(20, 5),
#     num.trees = 2500,
#     parallel = TRUE
# )
# save(smcrf_results_multi_param, file = "smc-drf.rda")
# load("smc-drf.rda")
# #---Plot marginal distributions compare
# plots_marginal <- plot_compare_marginal(
#     plots = plots_marginal,
#     xlimit = range,
#     abc_results = smcrf_results_multi_param,
#     parameters_truth = parameters_truth,
#     parameters_labels = parameters_labels,
#     plot_hist = TRUE
# )
# plot_smcrf_marginal(
#     smcrf_results = smcrf_results_multi_param,
#     parameters_labels = parameters_labels,
#     plot_hist = TRUE
# )
# =================================================TRAJECTORIES BOX PLOT
# plot_boxplot <- function(drf_results, smcrf_results, statistics_c3_small, statistics_c3_large) {
#     library(ggplot2)
#     library(reshape2) # Ensure reshape2 is loaded for melting data frames
#     # stats_id <- c("P", "E", "S", "ES")
#     stats_id <- c("mRNA")
#     nIterations <- drf_results[["nIterations"]]
#     plots <- list()
#     color_scheme <- c(
#         "Prior Distribution" = "gray",
#         "True Posterior Distribution" = "black",
#         "ABC-rejection" = "forestgreen",
#         "ABC-RF" = "magenta4",
#         "DRF" = "royalblue2",
#         "MCMC" = "goldenrod2",
#         "ABC-MCMC" = "goldenrod2",
#         "ABC-SMC" = "goldenrod2",
#         "SMC-RF for single parameters" = "salmon",
#         "SMC-RF for multiple parameters" = "salmon"
#     )
#     #---Set up legend order for plotting
#     legend_order <- c(
#         "Prior Distribution",
#         "True Posterior Distribution",
#         "ABC-rejection",
#         "ABC-MCMC",
#         "ABC-SMC",
#         "ABC-RF",
#         "DRF",
#         "MCMC",
#         "SMC-RF for single parameters",
#         "SMC-RF for multiple parameters"
#     )
#     # Factor levels for time
#     time_levels <- as.character(c(4, 6, 10, 20, 30, 40))

#     # for (stat_id in stats_id) {
#     #     # Extract data for the current statistic ID
#     #     posterior_data_smcrf <- smcrf_results[[paste0("Iteration_", nIterations + 1)]][["statistics"]][
#     #         1:1000,
#     #         grepl(paste0("^", stat_id, "_"), colnames(smcrf_results[[paste0("Iteration_", nIterations + 1)]][["statistics"]]))
#     #     ]
#     #     posterior_data_drf <- drf_results[[paste0("Iteration_", 1 + 1)]][["statistics"]][
#     #         1:1000,
#     #         grepl(paste0("^", stat_id, "_"), colnames(drf_results[[paste0("Iteration_", 1 + 1)]][["statistics"]]))
#     #     ]
#     #     posterior_data_drf_long <- melt(posterior_data_drf, variable.name = "Stat_type", value.name = "Value")
#     #     posterior_data_drf_long$time <- paste0(1 * as.integer(gsub("[^0-9]", "", posterior_data_drf_long$Stat_type)))
#     #     posterior_data_drf_long$legend <- "DRF"
#     #     print(posterior_data_drf_long)
#     #     posterior_data_smcrf_long <- melt(posterior_data_smcrf, variable.name = "Stat_type", value.name = "Value")
#     #     posterior_data_smcrf_long$time <- paste0(1 * as.integer(gsub("[^0-9]", "", posterior_data_smcrf_long$Stat_type)))
#     #     posterior_data_smcrf_long$legend <- "SMC-RF for multiple parameters"
#     #     print(posterior_data_smcrf_long)
#     #     total_posterior <- rbind(posterior_data_drf_long, posterior_data_smcrf_long)
#     #     target <- smcrf_results[["statistics_target"]][, grepl(paste0("^", stat_id, "_"), colnames(smcrf_results[["statistics_target"]]))]
#     #     target_long <- reshape2::melt(target, variable.name = "Stat_type", value.name = "Value")
#     #     target_long$time <- paste0(1 * as.integer(gsub("[^0-9]", "", target_long$Stat_type)))
#     #     # target_c3_small <- statistics_c3_small[, grepl(paste0("^", stat_id, "_"), colnames(statistics_c3_small))]
#     #     # target_c3_small_long <- reshape2::melt(target_c3_small, variable.name = "Stat_type", value.name = "Value")
#     #     # target_c3_small_long$time <- paste0(1 * as.integer(gsub("[^0-9]", "", target_c3_small_long$Stat_type)))
#     #     # target_c3_large <- statistics_c3_large[, grepl(paste0("^", stat_id, "_"), colnames(statistics_c3_large))]
#     #     # target_c3_large_long <- reshape2::melt(target_c3_large, variable.name = "Stat_type", value.name = "Value")
#     #     # target_c3_large_long$time <- paste0(1 * as.integer(gsub("[^0-9]", "", target_c3_large_long$Stat_type)))
#     #     # Convert time to factor with explicit ordering
#     #     posterior_data_drf_long$time <- factor(posterior_data_drf_long$time, levels = time_levels)
#     #     posterior_data_smcrf_long$time <- factor(posterior_data_smcrf_long$time, levels = time_levels)
#     #     target_long$time <- factor(target_long$time, levels = time_levels)
#     #     # target_c3_small_long$time <- factor(target_c3_small_long$time, levels = time_levels)
#     #     # target_c3_large_long$time <- factor(target_c3_large_long$time, levels = time_levels)
#     #     plots[[stat_id]] <- ggplot() +
#     #         geom_boxplot(data = total_posterior, aes(x = time, y = Value, color = legend, fill = legend), alpha = 0.8) +
#     #         geom_line(data = target_long, aes(x = time, y = Value, group = 1), color = "black", size = 1.5, show.legend = FALSE) +
#     #         geom_point(data = target_long, aes(x = time, y = Value, group = 1), color = "black", size = 10, show.legend = FALSE) +
#     #         # geom_line(data = target_c3_small_long, aes(x = time, y = Value, group = 1), color = "orange", size = 1.5, show.legend = FALSE) +
#     #         # geom_point(data = target_c3_small_long, aes(x = time, y = Value, group = 1), color = "orange", size = 10, show.legend = FALSE) +
#     #         # geom_line(data = target_c3_large_long, aes(x = time, y = Value, group = 1), color = "forestgreen", size = 1.5, show.legend = FALSE) +
#     #         # geom_point(data = target_c3_large_long, aes(x = time, y = Value, group = 1), color = "forestgreen", size = 10, show.legend = FALSE) +
#     #         labs(x = "Time", y = stat_id) +
#     #         ylim(0, 100) +
#     #         scale_fill_manual(values = color_scheme, name = "") +
#     #         scale_color_manual(values = color_scheme, name = "") +
#     #         theme(
#     #             text = element_text(size = 50),
#     #             panel.background = element_rect(fill = "white", colour = "white"),
#     #             panel.grid.major = element_line(colour = "white"),
#     #             panel.grid.minor = element_line(colour = "white"),
#     #             plot.margin = unit(c(1, 1, 1, 1), "cm"),
#     #             legend.position = "top",
#     #             legend.justification = c(0, 0.5)
#     #         )
#     #     # Output to PNG file
#     #     file_name <- paste0("boxplot-statistics=", stat_id, ".png")
#     #     png(file_name, res = 150, width = 3000, height = 1500, units = "px")
#     #     print(plots[[stat_id]])
#     #     dev.off()
#     # }
#     # ==============================The one for average from 10 sims
#     # for (stat_id in stats_id) {
#     #     # Extract data for the current statistic ID
#     #     # posterior_data_smcrf <- smcrf_results[[paste0("Iteration_", nIterations + 1)]][["statistics"]][
#     #     #     1:1000,
#     #     #     grepl(paste0("^", stat_id, "_", "[0-9]+", "_mean"), colnames(smcrf_results[[paste0("Iteration_", nIterations + 1)]][["statistics"]]))
#     #     # ]
#     #     posterior_data_drf <- drf_results[[paste0("Iteration_", 1 + 1)]][["statistics"]][
#     #         1:1000,
#     #         grepl(paste0("^", stat_id, "_", "[0-9]+", "min_mean"), colnames(drf_results[[paste0("Iteration_", 1 + 1)]][["statistics"]]))
#     #     ]

#     #     posterior_data_drf_long <- melt(posterior_data_drf, variable.name = "Stat_type", value.name = "Value")
#     #     posterior_data_drf_long$time <- paste0(1 * as.integer(gsub("[^0-9]", "", posterior_data_drf_long$Stat_type)))
#     #     posterior_data_drf_long$legend <- "DRF"
#     #     print(posterior_data_drf_long)
#     #     # posterior_data_smcrf_long <- melt(posterior_data_smcrf, variable.name = "Stat_type", value.name = "Value")
#     #     # posterior_data_smcrf_long$time <- paste0(1 * as.integer(gsub("[^0-9]", "", posterior_data_smcrf_long$Stat_type)))
#     #     # posterior_data_smcrf_long$legend <- "SMC-RF for multiple parameters"
#     #     # print(posterior_data_smcrf_long)
#     #     # total_posterior <- rbind(posterior_data_drf_long, posterior_data_smcrf_long)
#     #     target <- drf_results[["statistics_target"]][, grepl(paste0("^", stat_id, "_", "[0-9]+", "min_mean"), colnames(drf_results[["statistics_target"]]))]
#     #     target_long <- reshape2::melt(target, variable.name = "Stat_type", value.name = "Value")
#     #     target_long$time <- paste0(1 * as.integer(gsub("[^0-9]", "", target_long$Stat_type)))
#     #     # Convert time to factor with explicit ordering
#     #     posterior_data_drf_long$time <- factor(posterior_data_drf_long$time, levels = time_levels)
#     #     # posterior_data_smcrf_long$time <- factor(posterior_data_smcrf_long$time, levels = time_levels)
#     #     target_long$time <- factor(target_long$time, levels = time_levels)

#     #     plots[[stat_id]] <- ggplot() +
#     #         geom_boxplot(data = posterior_data_drf_long, aes(x = time, y = Value, color = legend, fill = legend), alpha = 0.8) +
#     #         geom_line(data = target_long, aes(x = time, y = Value, group = 1), color = "black", size = 1.5, show.legend = FALSE) +
#     #         geom_point(data = target_long, aes(x = time, y = Value, group = 1), color = "black", size = 10, show.legend = FALSE) +
#     #         labs(x = "Time", y = stat_id) +
#     #         ylim(0, 1000) +
#     #         scale_fill_manual(values = color_scheme, name = "") +
#     #         scale_color_manual(values = color_scheme, name = "") +
#     #         theme(
#     #             text = element_text(size = 50),
#     #             panel.background = element_rect(fill = "white", colour = "white"),
#     #             panel.grid.major = element_line(colour = "white"),
#     #             panel.grid.minor = element_line(colour = "white"),
#     #             plot.margin = unit(c(1, 1, 1, 1), "cm"),
#     #             legend.position = "top",
#     #             legend.justification = c(0, 0.5)
#     #         )

#     #     # Output to PNG file
#     #     file_name <- paste0("boxplot-statistics-mean=", stat_id, ".png")
#     #     png(file_name, res = 150, width = 3000, height = 1500, units = "px")
#     #     print(plots[[stat_id]])
#     #     dev.off()
#     # }
#     # ==============================The one for var from 10 sims
#     for (stat_id in stats_id) {
#         # Extract data for the current statistic ID
#         # posterior_data_smcrf <- smcrf_results[[paste0("Iteration_", nIterations + 1)]][["statistics"]][
#         #     1:1000,
#         #     grepl(paste0("^", stat_id, "_", "[0-9]+", "_mean"), colnames(smcrf_results[[paste0("Iteration_", nIterations + 1)]][["statistics"]]))
#         # ]
#         posterior_data_drf <- drf_results[[paste0("Iteration_", 1 + 1)]][["statistics"]][
#             1:1000,
#             grepl(paste0("^", stat_id, "_", "[0-9]+", "min_var"), colnames(drf_results[[paste0("Iteration_", 1 + 1)]][["statistics"]]))
#         ]
#         posterior_data_drf_long <- melt(posterior_data_drf, variable.name = "Stat_type", value.name = "Value")
#         posterior_data_drf_long$time <- paste0(1 * as.integer(gsub("[^0-9]", "", posterior_data_drf_long$Stat_type)))
#         posterior_data_drf_long$legend <- "DRF"
#         print(posterior_data_drf_long)
#         # posterior_data_smcrf_long <- melt(posterior_data_smcrf, variable.name = "Stat_type", value.name = "Value")
#         # posterior_data_smcrf_long$time <- paste0(1 * as.integer(gsub("[^0-9]", "", posterior_data_smcrf_long$Stat_type)))
#         # posterior_data_smcrf_long$legend <- "SMC-RF for multiple parameters"
#         # print(posterior_data_smcrf_long)
#         # total_posterior <- rbind(posterior_data_drf_long, posterior_data_smcrf_long)
#         target <- drf_results[["statistics_target"]][, grepl(paste0("^", stat_id, "_", "[0-9]+", "min_var"), colnames(drf_results[["statistics_target"]]))]
#         target_long <- reshape2::melt(target, variable.name = "Stat_type", value.name = "Value")
#         target_long$time <- paste0(1 * as.integer(gsub("[^0-9]", "", target_long$Stat_type)))
#         # Convert time to factor with explicit ordering
#         posterior_data_drf_long$time <- factor(posterior_data_drf_long$time, levels = time_levels)
#         # posterior_data_smcrf_long$time <- factor(posterior_data_smcrf_long$time, levels = time_levels)
#         target_long$time <- factor(target_long$time, levels = time_levels)
#         plots[[stat_id]] <- ggplot() +
#             geom_boxplot(data = posterior_data_drf_long, aes(x = time, y = Value, color = legend, fill = legend), alpha = 0.8) +
#             geom_line(data = target_long, aes(x = time, y = Value, group = 1), color = "black", size = 1.5, show.legend = FALSE) +
#             geom_point(data = target_long, aes(x = time, y = Value, group = 1), color = "black", size = 10, show.legend = FALSE) +
#             labs(x = "Time", y = stat_id) +
#             ylim(0, 1000) +
#             scale_fill_manual(values = color_scheme, name = "") +
#             scale_color_manual(values = color_scheme, name = "") +
#             theme(
#                 text = element_text(size = 50),
#                 panel.background = element_rect(fill = "white", colour = "white"),
#                 panel.grid.major = element_line(colour = "white"),
#                 panel.grid.minor = element_line(colour = "white"),
#                 plot.margin = unit(c(1, 1, 1, 1), "cm"),
#                 legend.position = "top",
#                 legend.justification = c(0, 0.5)
#             )
#         # Output to PNG file
#         file_name <- paste0("boxplot-statistics-var=", stat_id, ".png")
#         png(file_name, res = 150, width = 3000, height = 1500, units = "px")
#         print(plots[[stat_id]])
#         dev.off()
#     }
# }
# setwd("/Users/xiangzijin/Documents/ABC_SMCRF/CME/4-state/0523_results_nRuns=1;nSims=3000")
# load("munsky_drf.rda")
# load("smc-drf.rda")
# plot_boxplot(
#     drf_results = abcrf_results,
#     smcrf_results = smcrf_results_multi_param
# )
# setwd("/Users/xiangzijin/Documents/ABC_SMCRF/CME/4-state/0524_results_nRuns=10;nsims=3000")
# setwd("/Users/xiangzijin/Documents/ABC_SMCRF/CME/4-state/0529_results_nRuns=100;nsims=1000")
# load("munsky_drf.rda")
# load("smc-drf.rda")
# plot_boxplot(
#     drf_results = abcrf_results
# )
# ============   Variance Importance of each statistic in inferring gamma
load("abc-rf.rda")
statistics_target <- abcrf_results$statistics_target
png(paste0("NEUTRAL_abcrf_theta_variable_importance_ver.png"), width = 3000, height = 6000, res = 150, pointsize = 80)
# png(paste0("NEUTRAL_abcrf_theta_variable_importance_ver.png"), width = 30, height = 45, units = "in", res = 150)
n.var <- min(30, length(abcrf_results$Iteration_1$rf_model$model.rf$variable.importance))
imp <- abcrf_results$Iteration_1$rf_model$model.rf$variable.importance
names(imp) <- colnames(statistics_target)
ord <- rev(order(imp, decreasing = TRUE)[1:n.var])
xmin <- 0
xlim <- c(xmin, max(imp) + 1)
dotchart(imp[ord], pch = 19, xlab = "Variable Importance", ylab = "", xlim = xlim, main = NULL, bg = "white", cex = 0.5)
dev.off()
