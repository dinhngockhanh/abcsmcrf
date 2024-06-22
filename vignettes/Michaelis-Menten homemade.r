# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/R"
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ZIJIN - Macbook
R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/CME/m-w/home_made/old_example/timepoints/1-10"
# R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/adaptive"
R_libPaths <- ""
R_libPaths_extra <- "/Users/xiangzijin/SMC-RF/R"
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
        props <- reaction_propensities(state, parameters)
        #   Compute total propensity
        total_prop <- sum(props)
        if (total_prop == 0) {
            delta_time <- Inf
        } else {
            #   Compute time to next reaction & reaction to occur
            time_begin <- Sys.time()
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
model <- function(parameters, parallel = TRUE) {
    nNoise <- 0
    if (exists("nSimulations")) {
        nSimulations <<- nSimulations + nrow(parameters)
    }
    if (parallel) {
        library(parallel)
        library(pbapply)
        library(data.table)
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cl, varlist = c("SSA", "MW_model"))
        stats <- pblapply(
            cl = cl, X = 1:nrow(parameters),
            FUN = function(i) {
                MW_model(parameters$c1[i], parameters$c2[i], parameters$c3[i])
            }
        )
        stopCluster(cl)
        stats <- rbindlist(stats)
        class(stats) <- "data.frame"
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            stats <- rbind(stats, MW_model(parameters$c1[i], parameters$c2[i], parameters$c3[i]))
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
MW_model <- function(c1_input, c2_input, c3_input, parallel = FALSE) {
    time_points <- seq(1, 10, by = 1)
    nA <- 6.023e23
    vol <- 1e-15
    initial_state <- list(
        P = 0,
        E = round(2e-7 * nA * vol),
        S = round(5e-7 * nA * vol),
        ES = 0
    )
    parameters <- list(
        c1 = c1_input,
        c2 = (10^c2_input) / (nA * vol),
        c3 = (10^c3_input)
    )
    reaction_propensities <- function(state, parameters) {
        c(
            parameters$c1 * state$ES,
            parameters$c2 * state$E * state$S,
            parameters$c3 * state$ES
        )
    }
    reaction_stoichiometries <- list(
        list(P = 1, E = 1, S = 0, ES = -1),
        list(P = 0, E = -1, S = -1, ES = 1),
        list(P = 0, E = 1, S = 1, ES = -1)
    )
    time_begin <- Sys.time()
    tmp <- SSA(
        initial_state = initial_state,
        parameters = parameters,
        reaction_propensities = reaction_propensities,
        reaction_stoichiometries = reaction_stoichiometries,
        time_points = time_points
    )
    time_end <- Sys.time()
    print(time_end - time_begin)
    stats <- data.frame(matrix(c(c1_input, c2_input, c3_input, c(apply(tmp, 2, c))), nrow = 1))
    colnames(stats) <- c("c1", "c2", "c3", paste0("E_", 1:length(time_points)), paste0("S_", 1:length(time_points)), paste0("ES_", 1:length(time_points)), paste0("P_", 1:length(time_points)))
    return(stats)
}
# ========================The one with average
# MW_model <- function(c1_input, c2_input, c3_input, parallel = FALSE) {
#     nRuns <- 10
#     time_points <- seq(1, 10, by = 1)
#     nA <- 6.023e23
#     vol <- 1e-15
#     initial_state <- list(
#         P = 0,
#         E = round(2e-7 * nA * vol),
#         S = round(5e-7 * nA * vol),
#         ES = 0
#     )
#     parameters <- list(
#         c1 = c1_input,
#         c2 = (10^c2_input) / (nA * vol),
#         c3 = (10^c3_input)
#     )
#     reaction_propensities <- function(state, parameters) {
#         c(
#             parameters$c1 * state$ES,
#             parameters$c2 * state$E * state$S,
#             parameters$c3 * state$ES
#         )
#     }
#     reaction_stoichiometries <- list(
#         list(P = 1, E = 1, S = 0, ES = -1),
#         list(P = 0, E = -1, S = -1, ES = 1),
#         list(P = 0, E = 1, S = 1, ES = -1)
#     )
#     time_begin <- Sys.time()
#     #------------
#     MMs <- matrix(0, nrow = nRuns, ncol = length(time_points) * 4)
#     for (i in 1:nRuns) {
#         tmp <- SSA(
#             initial_state = initial_state,
#             parameters = parameters,
#             reaction_propensities = reaction_propensities,
#             reaction_stoichiometries = reaction_stoichiometries,
#             time_points = time_points
#         )
#         MMs[i, ] <- c(apply(tmp, 2, c))
#     }
#     time_end <- Sys.time()
#     print(time_end - time_begin)
#     stats <- data.frame(matrix(c(c1_input, c2_input, c3_input, c(apply(MMs, 2, mean)), c(apply(MMs, 2, var))), nrow = 1))
#     colnames(stats) <- c("c1", "c2", "c3", paste0("E_", 1:length(time_points), "_mean"), paste0("S_", 1:length(time_points), "_mean"), paste0("ES_", 1:length(time_points), "_mean"), paste0("P_", 1:length(time_points), "_mean"), paste0("E_", 1:length(time_points), "_var"), paste0("S_", 1:length(time_points), "_var"), paste0("ES_", 1:length(time_points), "_var"), paste0("P_", 1:length(time_points), "_var"))
#     return(stats)
# }
# =====================================================Target statistics
parameters_truth <- data.frame(
    c1 = 0.1,
    c2 = 6,
    c3 = -4
)
statistics_target <- model(parameters = parameters_truth, parallel = FALSE)[-c(1:ncol(parameters_truth))]
statistics_target
# ===========================================Target statistics with c3=-5
parameters_c3_small <- data.frame(
    c1 = 0.1,
    c2 = 6,
    c3 = -5
)
statistics_c3_small <- model(parameters = parameters_c3_small, parallel = FALSE)[-c(1:ncol(parameters_c3_small))]
statistics_c3_small
# ===========================================Target statistics with c3=-5
parameters_c3_large <- data.frame(
    c1 = 0.1,
    c2 = 6,
    c3 = -3
)
statistics_c3_large <- model(parameters = parameters_c3_large, parallel = FALSE)[-c(1:ncol(parameters_c3_large))]
statistics_c3_large
# ======================================Model for parameter perturbation
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters, after perturbation
perturb <- function(parameters) {
    if (any(grepl("c1", colnames(parameters)))) {
        parameters[["c1"]] <- parameters[["c1"]] + runif(nrow(parameters), min = -0.05, max = 0.05)
    } else if (any(grepl("c2", colnames(parameters)))) {
        parameters[["c2"]] <- parameters[["c2"]] + runif(nrow(parameters), min = -0.1, max = 0.1)
    } else if (any(grepl("c3", colnames(parameters)))) {
        parameters[["c3"]] <- parameters[["c3"]] + runif(nrow(parameters), min = -0.1, max = 0.1)
    }
    return(parameters)
}
# ======================================Define ranges for the parameters
range <- data.frame(
    parameter = c("c1", "c2", "c3"),
    min = c(0, 5, -5),
    max = c(1, 7, -3)
)
# ========================================Initial guesses for parameters
# ====================================(sampled from prior distributions)
c1 <- runif(100000, 0, 1)
c2 <- runif(100000, 5, 7)
c3 <- runif(100000, -5, -3)
parameters_initial <- data.frame(
    c1 = c1,
    c2 = c2,
    c3 = c3
)
# ====================================Labels for parameters in the plots
parameters_labels <- data.frame(
    parameter = c("c1", "c2", "c3"),
    label = c("expression(c[1])", "expression(c[2])", "expression(c[3])")
)
# # ========================================DRF
# #---Run DRF
# abcrf_results <- smcrf(
#     method = "smcrf-multi-param",
#     statistics_target = statistics_target,
#     parameters_initial = parameters_initial,
#     model = model,
#     perturb = perturb,
#     range = range,
#     nParticles = rep(20000, 1),
#     num.trees = 2500,
#     parallel = TRUE
# )
# save(abcrf_results, file = "drf.rda")
drf_results <- load("drf.rda")
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
#     nParticles = rep(4000, 5),
#     num.trees = 2500,
#     parallel = TRUE
# )
# save(smcrf_results_multi_param, file = "smc-drf.rda")
load("smc-drf.rda")
# #---Plot marginal distributions compare
plots_marginal <- plot_compare_marginal(
    plots = plots_marginal,
    xlimit = range,
    abc_results = smcrf_results_multi_param,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)
plot_smcrf_marginal(
    smcrf_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)
# # ===================================================================RF
# #---Run ABC-RF
# abcrf_results <- smcrf(
#     method = "smcrf-single-param",
#     statistics_target = statistics_target,
#     parameters_initial = parameters_initial,
#     model = model,
#     # perturb = perturb,
#     perturb = "Beaumont",
#     range = range,
#     ntree = 2500,
#     save_model = FALSE,
#     nParticles = rep(20000, 1),
#     parallel = TRUE
# )
# save(abcrf_results, file = "abc-rf.rda")
# load("abc-rf.rda")
# plots_marginal <- plot_compare_marginal(
#     # plots = plots_marginal,
#     xlimit = range,
#     abc_results = abcrf_results,
#     parameters_truth = parameters_truth,
#     parameters_labels = parameters_labels,
#     plot_hist = TRUE
# )
# # ========================================SMC-RF for single parameters
# #---Run SMC-RF for single parameters
# smcrf_results_single_param <- smcrf(
#     method = "smcrf-single-param",
#     statistics_target = statistics_target,
#     parameters_initial = parameters_initial,
#     model = model,
#     # perturb = perturb,
#     perturb = "Beaumont",
#     range = range,
#     nParticles = rep(4000, 5),
#     ntree = 2500,
#     save_model = FALSE,
#     parallel = TRUE
# )
# save(smcrf_results_single_param, file = "smc-rf.rda")
# load("smc-rf.rda")
# #---Plot marginal distributions compare
# plots_marginal <- plot_compare_marginal(
#     plots = plots_marginal,
#     xlimit = range,
#     abc_results = smcrf_results_single_param,
#     parameters_truth = parameters_truth,
#     parameters_labels = parameters_labels,
#     plot_hist = TRUE
# )
# plot_smcrf_marginal(
#     smcrf_results = smcrf_results_single_param,
#     parameters_labels = parameters_labels,
#     plot_hist = TRUE
# )
# ============================Plot joint distributions compare c1 and c2
# smcrf_results_multi_param
# nIterations <- smcrf_results_multi_param[["nIterations"]]
# posterior_df <- data.frame(
#     x = smcrf_results_multi_param[[paste0("Iteration_", nIterations + 1)]]$parameters[[parameters_labels$parameter[1]]],
#     y = smcrf_results_multi_param[[paste0("Iteration_", nIterations + 1)]]$parameters[[parameters_labels$parameter[2]]]
# )
# plots <- ggplot() + geom_density_2d_filled(data = posterior_df, aes(x = x, y = y), show.legend = FALSE)
# plots <- plots +
#     theme(
#         text = element_text(size = 50),
#         panel.background = element_rect(fill = "white", colour = "white"),
#         panel.grid.major = element_line(colour = "white"),
#         panel.grid.minor = element_line(colour = "white"),
#         legend.position = "top",
#         legend.justification = c(0, 0.5)
#     )
# #---Print joint distribution plot
# file_name <- paste0("Joint-parameters=", parameters_labels$parameter[1], "-vs-", parameters_labels$parameter[2], ".png")
# png(file_name, res = 150, width = 30, height = 31, units = "in", pointsize = 12)
# print(plots)
# dev.off()
# return(plots)

# =================================================TRAJECTORIES BOX PLOT
plot_boxplot <- function(drf_results, smcrf_results) {
    library(ggplot2)
    library(reshape2) # Ensure reshape2 is loaded for melting data frames
    stats_id <- c("P", "E", "S", "ES")
    nIterations <- smcrf_results[["nIterations"]]
    plots <- list()
    color_scheme <- c(
        "Prior Distribution" = "gray",
        "True Posterior Distribution" = "black",
        "ABC-rejection" = "forestgreen",
        "ABC-RF" = "magenta4",
        "ABC-DRF" = "royalblue2",
        "MCMC" = "goldenrod2",
        "ABC-MCMC" = "goldenrod2",
        "ABC-SMC" = "goldenrod2",
        "ABC-SMC-RF" = "salmon",
        "ABC-SMC-DRF" = "salmon"
    )
    #---Set up legend order for plotting
    legend_order <- c(
        "Prior Distribution",
        "True Posterior Distribution",
        "ABC-rejection",
        "ABC-MCMC",
        "ABC-SMC",
        "ABC-RF",
        "ABC-DRF",
        "MCMC",
        "ABC-SMC-RF",
        "ABC-SMC-DRF"
    )
    # Factor levels for time
    time_levels <- as.character(seq(1, 10, by = 1))

    for (stat_id in stats_id) {
        # Extract data for the current statistic ID
        posterior_data_smcrf <- smcrf_results[[paste0("Iteration_", nIterations + 1)]][["statistics"]][
            1:1000,
            grepl(paste0("^", stat_id, "_"), colnames(smcrf_results[[paste0("Iteration_", nIterations + 1)]][["statistics"]]))
        ]
        posterior_data_drf <- drf_results[[paste0("Iteration_", 1 + 1)]][["statistics"]][
            1:1000,
            grepl(paste0("^", stat_id, "_"), colnames(drf_results[[paste0("Iteration_", 1 + 1)]][["statistics"]]))
        ]
        posterior_data_drf_long <- melt(posterior_data_drf, variable.name = "Stat_type", value.name = "Value")
        posterior_data_drf_long$time <- paste0(1 * as.integer(gsub("[^0-9]", "", posterior_data_drf_long$Stat_type)))
        posterior_data_drf_long$legend <- "ABC-DRF"
        posterior_data_smcrf_long <- melt(posterior_data_smcrf, variable.name = "Stat_type", value.name = "Value")
        posterior_data_smcrf_long$time <- paste0(1 * as.integer(gsub("[^0-9]", "", posterior_data_smcrf_long$Stat_type)))
        posterior_data_smcrf_long$legend <- "ABC-SMC-DRF"
        total_posterior <- rbind(posterior_data_drf_long, posterior_data_smcrf_long)
        target <- smcrf_results[["statistics_target"]][, grepl(paste0("^", stat_id, "_"), colnames(smcrf_results[["statistics_target"]]))]
        target_long <- reshape2::melt(target, variable.name = "Stat_type", value.name = "Value")
        target_long$time <- paste0(1 * as.integer(gsub("[^0-9]", "", target_long$Stat_type)))
        # Convert time to factor with explicit ordering
        posterior_data_drf_long$time <- factor(posterior_data_drf_long$time, levels = time_levels)
        posterior_data_smcrf_long$time <- factor(posterior_data_smcrf_long$time, levels = time_levels)
        target_long$time <- factor(target_long$time, levels = time_levels)
        plots[[stat_id]] <- ggplot() +
            geom_boxplot(data = total_posterior, aes(x = time, y = Value, color = legend, fill = legend), alpha = 0.8) +
            geom_line(data = target_long, aes(x = time, y = Value, group = 1), color = "black", size = 1.5, show.legend = FALSE) +
            geom_point(data = target_long, aes(x = time, y = Value, group = 1), color = "black", size = 10, show.legend = FALSE) +
            labs(x = "Time", y = stat_id) +
            scale_fill_manual(values = color_scheme, name = "") +
            scale_color_manual(values = color_scheme, name = "") +
            theme(
                text = element_text(size = 50),
                panel.background = element_rect(fill = "white", colour = "white"),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white"),
                plot.margin = unit(c(1, 1, 1, 1), "cm"),
                legend.position = "top",
                legend.justification = c(0, 0.5)
            )
        # Output to PNG file
        file_name <- paste0("boxplot-statistics=", stat_id, ".png")
        png(file_name, res = 150, width = 3000, height = 1500, units = "px")
        print(plots[[stat_id]])
        dev.off()
    }
}
load("smc-drf.rda")
plot_boxplot(
    drf_results = abcrf_results,
    smcrf_results = smcrf_results_multi_param
)
