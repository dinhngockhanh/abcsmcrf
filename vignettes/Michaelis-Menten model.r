library(abcsmcrf)
# =====================================Function to plot the trajectories
plot_boxplot <- function(drf_results, smcrf_results) {
    library(ggplot2)
    library(reshape2)
    stats_id <- c("P", "E", "S", "ES")
    nIterations <- smcrf_results[["nIterations"]]
    plots <- list()
    color_scheme <- c(
        "ABC-DRF" = "#0072B5",
        "ABC-SMC-DRF" = "#BC3C29"
    )
    #---Set up legend order for plotting
    legend_order <- c(
        "ABC-DRF",
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
        stats <- parLapply(
            cl = cl, 1:nrow(parameters),
            function(i) MW_model(parameters$c1[i], parameters$c2[i], parameters$c3[i])
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
# =====================================================Target statistics
parameters_truth <- data.frame(
    c1 = 0.1,
    c2 = 6,
    c3 = -4
)
statistics_target <- model(parameters = parameters_truth, parallel = FALSE)[-c(1:ncol(parameters_truth))]
# ====================================================Prior distribution
dprior <- function(parameters, parameter_id = "all") {
    probs <- rep(1, nrow(parameters))
    if (parameter_id %in% c("all", "c1")) {
        probs <- probs * dunif(parameters[["c1"]], min = 0, max = 1)
    }
    if (parameter_id %in% c("all", "c2")) {
        probs <- probs * dunif(parameters[["c2"]], min = 5, max = 7)
    }
    if (parameter_id %in% c("all", "c3")) {
        probs <- probs * dunif(parameters[["c3"]], min = -5, max = -3)
    }
    return(probs)
}
rprior <- function(Nparameters) {
    parameters <- data.frame(
        c1 = runif(Nparameters, min = 0, max = 1),
        c2 = runif(Nparameters, min = 5, max = 7),
        c3 = runif(Nparameters, min = -5, max = -3)
    )
    return(parameters)
}
# =============================================Perturbation distribution
dperturb <- function(parameters, parameters_previous, parameters_previous_sampled, iteration, parameter_id = "all") {
    probs <- rep(1, nrow(parameters))
    if (parameter_id %in% c("all", "c1")) {
        probs <- probs * dunif(parameters[["c1"]], min = pmax(0, parameters_previous[["c1"]] - 0.05), max = pmin(1, parameters_previous[["c1"]] + 0.05))
    }
    if (parameter_id %in% c("all", "c2")) {
        probs <- probs * dunif(parameters[["c2"]], min = pmax(5, parameters_previous[["c2"]] - 0.1), max = pmin(7, parameters_previous[["c2"]] + 0.1))
    }
    if (parameter_id %in% c("all", "c3")) {
        probs <- probs * dunif(parameters[["c3"]], min = pmax(-5, parameters_previous[["c3"]] - 0.1), max = pmin(-3, parameters_previous[["c3"]] + 0.1))
    }
    return(probs)
}
rperturb <- function(parameters_unperturbed, parameters_previous_sampled, iteration) {
    parameters_perturbed <- parameters_unperturbed
    parameters_perturbed[["c1"]] <- runif(
        n = nrow(parameters_perturbed),
        min = pmax(0, parameters_unperturbed[["c1"]] - 0.05),
        max = pmin(1, parameters_unperturbed[["c1"]] + 0.05)
    )
    parameters_perturbed[["c2"]] <- runif(
        n = nrow(parameters_perturbed),
        min = pmax(5, parameters_unperturbed[["c2"]] - 0.1),
        max = pmin(7, parameters_unperturbed[["c2"]] + 0.1)
    )
    parameters_perturbed[["c3"]] <- runif(
        n = nrow(parameters_perturbed),
        min = pmax(-5, parameters_unperturbed[["c3"]] - 0.1),
        max = pmin(-3, parameters_unperturbed[["c3"]] + 0.1)
    )
    return(parameters_perturbed)
}
# ====================================Labels for parameters in the plots
parameters_labels <- data.frame(
    parameter = c("c1", "c2", "c3"),
    label = c("expression(c[1])", "expression(c[2])", "expression(c[3])")
)
# ===================================================================DRF
#---Run DRF
drf_results <- smcrf(
    method = "smcrf-multi-param",
    statistics_target = statistics_target,
    model = model,
    rprior = rprior,
    dprior = dprior,
    nParticles = rep(20000, 1),
    num.trees = 2500,
    parallel = TRUE
)
#---Plot marginal distributions compare
plots_marginal <- plot_compare_marginal(
    xlimit = data.frame(parameter = c("c1", "c2", "c3"), min = c(0, 5, -5), max = c(1, 7, -3)),
    abc_results = drf_results,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)
# ========================================SMC-RF for multiple parameters
#---Run SMC-RF for multiple parameters
smcrf_results_multi_param <- smcrf(
    method = "smcrf-multi-param",
    statistics_target = statistics_target,
    model = model,
    rprior = rprior,
    dprior = dprior,
    nParticles = rep(4000, 5),
    num.trees = 2500,
    rperturb = rperturb,
    dperturb = dperturb,
    parallel = TRUE
)
#---Plot marginal distributions compare
plots_marginal <- plot_compare_marginal(
    plots = plots_marginal,
    xlimit = data.frame(parameter = c("c1", "c2", "c3"), min = c(0, 5, -5), max = c(1, 7, -3)),
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
#---Plot box plot showing the trajectory
plot_boxplot(
    drf_results = drf_results,
    smcrf_results = smcrf_results_multi_param
)
# ===============================================MORRIS SENSITIVITY TEST
# ==========================================Model for Morris sensitivity
library(sensitivity)
morris_model <- function(parameters, parallel = TRUE) {
    parameters <- data.frame(parameters)
    colnames(parameters) <- c("c1", "c2", "c3")
    data <- t(as.matrix(model(parameters, parallel = TRUE)[-c(1:3)]))
    return(data)
}
x <- morrisMultOut(
    model = morris_model, factors = 3, r = 1000,
    design = list(type = "oat", levels = 5, grid.jump = 3), binf = c(0, 5, -5), bsup = c(1, 7, -3), scale = FALSE
)
#---Function to plot the morris output
plot_morris <- function(morris_output) {
    sensitivity_df <- data.frame(
        mu_star = apply(morris_output$ee, 2, function(morris_output) mean(abs(morris_output))),
        sigma = apply(morris_output$ee, 2, sd),
        label = c("c[1]", "c[2]", "c[3]")
    )
    plot_output <- ggplot() +
        geom_point(data = sensitivity_df, aes(x = mu_star, y = sigma), size = 10) +
        geom_text(data = sensitivity_df, aes(x = mu_star, y = sigma, label = label), vjust = 1, hjust = -0.5, parse = TRUE, size = 20) +
        xlim(0, 0.02) +
        ylim(0, 0.02) +
        theme(
            text = element_text(size = 50),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white"),
            legend.key.size = unit(10, "points"),
            legend.position = "top",
            legend.justification = c(0, 0.5)
        ) +
        labs(x = expression(mu^"*"), y = expression(sigma))
    file_name <- paste0("sensitivity-plot-morris.png")
    png(file_name, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
    print(plot_output)
    dev.off()
}
plot_morris(x)
