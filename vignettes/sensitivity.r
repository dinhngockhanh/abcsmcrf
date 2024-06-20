library(sensitivity)
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ZIJIN - Macbook
# R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/CME/m-w/home_made/old_example/timepoints/1-10;sensitivity"
R_workplace <- "/Users/xiangzijin/Downloads"
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
        clusterExport(cl, varlist = c("SSA", "MW_model", "MW_model_ses"))
        stats <- pblapply(
            cl = cl, X = 1:nrow(parameters),
            FUN = function(i) {
                MW_model_ses(parameters[i, 1], parameters[i, 2], parameters[i, 3])[-c(1:3)]
            }
        )
        stopCluster(cl)
        stats <- as.matrix(t(rbindlist(stats)))
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            stats <- cbind(stats, MW_model_ses(parameters[i, 1], parameters[i, 2], parameters[i, 3])[-c(1:3)])
        }
    }
    #   Add column names
    data <- data.frame(stats)
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
    stats <- data.frame(matrix(c(c(apply(tmp, 2, c))), nrow = 1))
    colnames(stats) <- NULL
    # colnames(stats) <- c(paste0("E_", 1:length(time_points)), paste0("S_", 1:length(time_points)), paste0("ES_", 1:length(time_points)), paste0("P_", 1:length(time_points)))
    return(stats)
}
MW_model_ses <- function(c1_input, c2_input, c3_input, parallel = FALSE) {
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
x <- morrisMultOut(
    model = model, factors = 3, r = 1000,
    design = list(type = "oat", levels = 5, grid.jump = 3), binf = c(0, 5, -5), bsup = c(1, 7, -3), scale = FALSE
)
save(x, file = "morris_sensitivity.rda")
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
load("morris_sensitivity.rda")
plot_morris(x)
