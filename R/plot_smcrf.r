plot_smcrf_marginal <- function(smcrf_results,
                                parameters_truth = NULL,
                                parameters_labels = NULL,
                                statistics_labels = NULL,
                                plot_statistics = FALSE,
                                alpha = 0.3) {
    nIterations <- smcrf_results[["nIterations"]]
    if (is.null(parameters_labels)) parameters_labels <- smcrf_results[["parameters_labels"]]
    if (is.null(statistics_labels)) statistics_labels <- smcrf_results[["statistics_labels"]]
    #---Set up color scheme for plotting
    color_scheme <- c(
        "True Posterior" = "black", "Data statistic" = "black",
        "Prior Distribution" = "gray",
        setNames(rev(rainbow(nIterations)), paste0("Iter. ", 1:nIterations))
    )
    #---Set up legend order for plotting
    legend_order <- c("True Posterior", "Data statistic", "Prior Distribution", paste0("Iter. ", 1:nIterations))
    #---Plot marginal distributions for each parameter
    for (i in 1:length(parameters_labels$parameter)) {
        parameter_id <- parameters_labels$parameter[i]
        #   Begin plot
        p <- ggplot()
        #   Plot true posterior distribution (if provided)
        if (!is.null(parameters_truth)) {
            true_posterior_df <- data.frame(value = parameters_truth[[parameter_id]], legend = "True Posterior")
            if (nrow(parameters_truth) == 1) {
                p <- p + geom_vline(data = true_posterior_df, aes(xintercept = value), linetype = "solid", linewidth = 5)
            } else {
                p <- p + geom_density(data = true_posterior_df, aes(x = value, fill = legend, color = legend), alpha = 1, linewidth = 2)
            }
        }
        #   Plot prior distribution
        prior_df <- data.frame(value = smcrf_results[["Iteration_1"]]$parameters[[parameter_id]], legend = "Prior Distribution")
        p <- p + geom_density(data = prior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
        #   Plot posterior distribution for each iteration
        for (iteration in 1:nIterations) {
            posterior_df <- data.frame(
                value = smcrf_results[[paste0("Iteration_", iteration + 1)]]$parameters[[parameter_id]],
                legend = paste0("Iter. ", iteration)
            )
            p <- p + geom_density(data = posterior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
        }
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
        file_name <- paste0(smcrf_results[["method"]], "-marginal-parameter=", parameter_id, ".png")
        png(file_name, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
        print(p)
        dev.off()
    }
    #---Plot marginal distributions for each statistic
    if (plot_statistics) {
        for (i in 1:length(statistics_labels$ID)) {
            statistic_id <- statistics_labels$ID[i]
            p <- ggplot()
            #   Plot statistic value from target data
            data_df <- data.frame(value = smcrf_results$statistics_target[[statistic_id]], legend = "Data statistic")
            p <- p + geom_vline(data = data_df, aes(xintercept = value), linetype = "solid", linewidth = 5)
            #   Plot prior distribution
            prior_df <- data.frame(value = smcrf_results[["Iteration_1"]]$statistics[[statistic_id]], legend = "Prior Distribution")
            p <- p + geom_density(data = prior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
            #   Plot posterior distribution for each iteration
            for (iteration in 1:nIterations) {
                posterior_df <- data.frame(
                    value = smcrf_results[[paste0("Iteration_", iteration + 1)]]$statistics[[statistic_id]],
                    legend = paste0("Iter. ", iteration)
                )
                p <- p + geom_density(data = posterior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
            }
            #   Add label for parameter
            if ("label" %in% colnames(statistics_labels)) {
                p <- p + labs(x = eval(parse(text = statistics_labels$label[i])))
            } else {
                p <- p + labs(x = statistic_id)
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
            file_name <- paste0(smcrf_results[["method"]], "-marginal-statistic=", statistic_id, ".png")
            png(file_name, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
            print(p)
            dev.off()
        }
    }
}

plot_compare_marginal <- function(plots = NULL,
                                  abc_results,
                                  parameters_truth = NULL,
                                  parameters_labels = NULL,
                                  statistics_labels = NULL,
                                  xlimit = NULL,
                                  sample_num = NULL,
                                  plot_statistics = FALSE,
                                  alpha = 0.3) {
    if (is.null(parameters_labels)) parameters_labels <- abc_results[["parameters_labels"]]
    if (is.null(statistics_labels)) statistics_labels <- abc_results[["statistics_labels"]]
    method <- abc_results[["method"]]
    #---Set up color scheme for plotting
    color_scheme <- c(
        "True Posterior" = "black",
        "ABC-rejection" = "forestgreen",
        "ABC-RF" = "royalblue4",
        "DRF" = "royalblue4",
        "ABC-MCMC" = "goldenrod2",
        "ABC-SMC" = "magenta4",
        "SMC-RF for single parameters" = "salmon",
        "SMC-RF for multiple parameters" = "salmon"
    )
    #---Set up legend order for plotting
    legend_order <- c(
        "True Posterior",
        "ABC-rejection",
        "ABC-MCMC",
        "ABC-SMC",
        "ABC-RF",
        "DRF",
        "SMC-RF for single parameters",
        "SMC-RF for multiple parameters"
    )
    #---Begin plots
    if (is.null(plots)) {
        plots <- list()
        for (parameter_id in parameters_labels$parameter) {
            plots$parameters[[parameter_id]] <- ggplot()
        }
        if (plot_statistics) {
            for (statistic_id in statistics_labels$ID) {
                plots$statistics[[statistic_id]] <- ggplot()
            }
        }
        new_plot <- TRUE
    } else {
        new_plot <- FALSE
    }
    #---Plot true posterior parameter distributions (if provided)
    if (!is.null(parameters_truth)) {
        for (parameter_id in parameters_labels$parameter) {
            true_posterior_df <- data.frame(value = parameters_truth[[parameter_id]], legend = "True Posterior")
            if (nrow(parameters_truth) == 1) {
                plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                    geom_vline(data = true_posterior_df, aes(xintercept = value), linetype = "solid", linewidth = 5)
            } else {
                plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                    geom_density(data = true_posterior_df, aes(x = value, fill = legend, color = legend), alpha = 0.8, linewidth = 2)
            }
        }
    }
    #---Plot statistic value from target data
    if (plot_statistics & new_plot) {
        for (statistic_id in statistics_labels$ID) {
            data_df <- data.frame(value = abc_results$statistics_target[[statistic_id]], legend = "Data statistic")
            plots$statistics[[statistic_id]] <- plots$statistics[[statistic_id]] +
                geom_vline(data = data_df, aes(xintercept = value), linetype = "solid", linewidth = 5)
        }
    }
    #---Extract final posterior distributions
    if (method == "smcrf-single-param") {
        nIterations <- abc_results[["nIterations"]]
        if (nIterations == 1) {
            legend_label <- "ABC-RF"
        } else {
            legend_label <- "SMC-RF for single parameters"
        }
        parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters
        if (!is.null(sample_num)) {
            parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters[1:sample_num, , drop = FALSE]
        }
        if (plot_statistics) statistics_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$statistics
    } else if (method == "smcrf-multi-param") {
        nIterations <- abc_results[["nIterations"]]
        if (nIterations == 1) {
            legend_label <- "DRF"
        } else {
            legend_label <- "SMC-RF for multiple parameters"
        }
        legend_label <- "SMC-RF for multiple parameters"
        parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters
        if (plot_statistics) statistics_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$statistics
    } else if (method == "abc-rejection") {
        legend_label <- "ABC-rejection"
        parameters_values <- abc_results[["Iteration_2"]]$parameters
        if (plot_statistics) statistics_values <- abc_results[["Iteration_2"]]$statistics
    } else if (method == "abc-smc") {
        legend_label <- "ABC-SMC"
        parameters_values <- abc_results[["Iteration_2"]]$parameters
        if (plot_statistics) statistics_values <- abc_results[["Iteration_2"]]$statistics
    } else if (method == "abc-mcmc") {
        legend_label <- "ABC-MCMC"
        parameters_values <- abc_results[["Iteration_2"]]$parameters
        if (plot_statistics) statistics_values <- abc_results[["Iteration_2"]]$statistics
    }
    #---Plot marginal distributions for each parameter
    for (parameter_id in parameters_labels$parameter) {
        posterior_df <- data.frame(
            value = parameters_values[[parameter_id]],
            legend = legend_label
        )
        plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
            geom_density(data = posterior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
        if (!is.null(xlimit)) {
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                xlim(xlimit[1], xlimit[2])
        }
    }
    if (plot_statistics) {
        for (statistic_id in statistics_labels$ID) {
            posterior_df <- data.frame(
                value = statistics_values[[statistic_id]],
                legend = legend_label
            )
            plots$statistics[[statistic_id]] <- plots$statistics[[statistic_id]] +
                geom_density(data = posterior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
        }
    }
    #---Add label for parameters and statistics
    if (new_plot == TRUE) {
        if ("label" %in% colnames(parameters_labels)) {
            for (parameter_id in parameters_labels$parameter) {
                plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                    labs(x = eval(parse(text = parameters_labels$label[which(parameters_labels$parameter == parameter_id)])))
            }
        } else {
            for (parameter_id in parameters_labels$parameter) {
                plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                    labs(x = parameter_id)
            }
        }
        if (plot_statistics) {
            if ("label" %in% colnames(statistics_labels)) {
                for (statistic_id in statistics_labels$ID) {
                    plots$statistics[[statistic_id]] <- plots$statistics[[statistic_id]] +
                        labs(x = eval(parse(text = statistics_labels$label[which(statistics_labels$ID == statistic_id)])))
                }
            } else {
                for (statistic_id in statistics_labels$ID) {
                    plots$statistics[[statistic_id]] <- plots$statistics[[statistic_id]] +
                        labs(x = statistic_id)
                }
            }
        }
    }
    #---Beautify plots
    if (new_plot) {
        for (parameter_id in parameters_labels$parameter) {
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
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
        }
        if (plot_statistics) {
            for (statistic_id in statistics_labels$ID) {
                plots$statistics[[statistic_id]] <- plots$statistics[[statistic_id]] +
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
            }
        }
    }
    #---Print marginal distribution plots
    for (parameter_id in parameters_labels$parameter) {
        file_name <- paste0("comparison-marginal-parameter=", parameter_id, ".png")
        png(file_name, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
        print(plots$parameters[[parameter_id]])
        dev.off()
    }
    if (plot_statistics) {
        for (statistic_id in statistics_labels$ID) {
            file_name <- paste0("comparison-marginal-statistic=", statistic_id, ".png")
            png(file_name, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
            print(plots$statistics[[statistic_id]])
            dev.off()
        }
    }
    return(plots)
}


plot_smcrf_joint <- function(smcrf_results,
                             parameters_truth = NULL,
                             parameters_labels = NULL,
                             lims = NULL,
                             nBins = 5) {
    library(MASS)
    library(dplyr)
    nIterations <- smcrf_results[["nIterations"]]
    method <- smcrf_results[["method"]]
    if (is.null(parameters_labels)) parameters_labels <- smcrf_results[["parameters_labels"]]
    if (nrow(parameters_labels) != 2) stop("ERROR: plot_smcrf_joint only works for two parameters. Please check parameters_labels.")
    #---Set up color scheme for plotting
    color_scheme <- c(
        "Prior Distribution" = "gray",
        setNames(rev(rainbow(nIterations)), paste0("Iter. ", 1:nIterations))
    )
    #---Set up legend order for plotting
    legend_order <- c("True Posterior", "Prior Distribution", paste0("Iter. ", 1:nIterations))
    #---Function to apply limits to data (if provided)
    apply_lims <- function(df) {
        df <- df %>%
            filter(
                x >= lims$min[which(lims$ID == parameters_labels$parameter[1])],
                x <= lims$max[which(lims$ID == parameters_labels$parameter[1])],
                y >= lims$min[which(lims$ID == parameters_labels$parameter[2])],
                y <= lims$max[which(lims$ID == parameters_labels$parameter[2])]
            )
    }
    #---Begin plot
    p <- ggplot()
    #---Plot true posterior distribution (if provided)
    if (!is.null(parameters_truth)) {
        true_posterior_df <- data.frame(
            x = parameters_truth[[parameters_labels$parameter[1]]],
            y = parameters_truth[[parameters_labels$parameter[2]]]
        )
        if (!is.null(lims)) true_posterior_df <- apply_lims(true_posterior_df)
        p <- p + geom_density_2d_filled(data = true_posterior_df, aes(x = x, y = y), show.legend = FALSE)
    }
    #---Plot prior distribution
    prior_df <- data.frame(
        x = smcrf_results[["Iteration_1"]]$parameters[[parameters_labels$parameter[1]]],
        y = smcrf_results[["Iteration_1"]]$parameters[[parameters_labels$parameter[2]]],
        legend = "Prior Distribution"
    )
    if (!is.null(lims)) prior_df <- apply_lims(prior_df)
    p <- p + geom_density_2d(data = prior_df, aes(x = x, y = y, color = legend), linewidth = 3, bins = nBins)
    #---Plot posterior distribution for each iteration
    for (iteration in 1:nIterations) {
        posterior_df <- data.frame(
            x = smcrf_results[[paste0("Iteration_", iteration + 1)]]$parameters[[parameters_labels$parameter[1]]],
            y = smcrf_results[[paste0("Iteration_", iteration + 1)]]$parameters[[parameters_labels$parameter[2]]],
            legend = paste0("Iter. ", iteration)
        )
        if (!is.null(lims)) posterior_df <- apply_lims(posterior_df)
        p <- p + geom_density_2d(data = posterior_df, aes(x = x, y = y, color = legend), linewidth = 3, bins = nBins)
    }
    #---Add label for parameter
    if ("label" %in% colnames(parameters_labels)) {
        p <- p + labs(x = eval(parse(text = parameters_labels$label[1])), y = eval(parse(text = parameters_labels$label[2])))
    } else {
        p <- p + labs(x = parameter_id[1], y = parameter_id[2])
    }
    #---Beautify plot
    p <- p +
        scale_color_manual(values = color_scheme, name = "", breaks = legend_order) +
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
    file_name <- paste0(smcrf_results[["method"]], "-joint-parameters=", parameters_labels$parameter[1], "-vs-", parameters_labels$parameter[2], ".png")
    png(file_name, res = 150, width = 30, height = 31, units = "in", pointsize = 12)
    print(p)
    dev.off()
}

plot_compare_joint <- function(plots = NULL,
                               abc_results,
                               parameters_truth = NULL,
                               parameters_labels = NULL,
                               lims = NULL,
                               nBins = 5) {
    library(MASS)
    library(dplyr)
    method <- abc_results[["method"]]
    if (is.null(parameters_labels)) parameters_labels <- abc_results[["parameters_labels"]]
    if (nrow(parameters_labels) != 2) stop("ERROR: plot_compare_joint only works for two parameters. Please check parameters_labels.")
    #---Set up color scheme for plotting
    color_scheme <- c(
        "True Posterior" = "black",
        "ABC-rejection" = "forestgreen",
        "ABC-RF" = "royalblue4",
        "DRF" = "royalblue4",
        "ABC-MCMC" = "goldenrod2",
        "ABC-SMC" = "magenta4",
        "SMC-RF for single parameters" = "salmon",
        "SMC-RF for multiple parameters" = "salmon"
    )
    #---Set up legend order for plotting
    legend_order <- c(
        "True Posterior",
        "ABC-rejection",
        "ABC-MCMC",
        "ABC-SMC",
        "ABC-RF",
        "DRF",
        "SMC-RF for single parameters",
        "SMC-RF for multiple parameters"
    )
    #---Function to apply limits to data (if provided)
    apply_lims <- function(df) {
        df <- df %>%
            filter(
                x >= lims$min[which(lims$ID == parameters_labels$parameter[1])],
                x <= lims$max[which(lims$ID == parameters_labels$parameter[1])],
                y >= lims$min[which(lims$ID == parameters_labels$parameter[2])],
                y <= lims$max[which(lims$ID == parameters_labels$parameter[2])]
            )
    }
    #---Begin plots
    if (is.null(plots)) {
        plots <- ggplot()
        new_plot <- TRUE
    } else {
        new_plot <- FALSE
    }
    #---Plot true posterior parameter distributions (if provided)
    if (!is.null(parameters_truth)) {
        true_posterior_df <- data.frame(
            x = parameters_truth[[parameters_labels$parameter[1]]],
            y = parameters_truth[[parameters_labels$parameter[2]]]
        )
        if (!is.null(lims)) true_posterior_df <- apply_lims(true_posterior_df)
        plots <- plots + geom_density_2d_filled(data = true_posterior_df, aes(x = x, y = y), show.legend = FALSE)
    }
    #---Extract final posterior distributions
    if (method == "smcrf-single-param") {
        nIterations <- abc_results[["nIterations"]]
        if (nIterations == 1) {
            legend_label <- "ABC-RF"
        } else {
            legend_label <- "SMC-RF for single parameters"
        }
        posterior_df <- data.frame(
            x = abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters[[parameters_labels$parameter[1]]],
            y = abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    } else if (method == "smcrf-multi-param") {
        nIterations <- abc_results[["nIterations"]]
        if (nIterations == 1) {
            legend_label <- "DRF"
        } else {
            legend_label <- "SMC-RF for multiple parameters"
        }
        posterior_df <- data.frame(
            x = abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters[[parameters_labels$parameter[1]]],
            y = abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    } else if (method == "abc-rejection") {
        legend_label <- "ABC-rejection"
        posterior_df <- data.frame(
            x = abc_results[["Iteration_2"]]$parameters[[parameters_labels$parameter[1]]],
            y = abc_results[["Iteration_2"]]$parameters[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    } else if (method == "abc-smc") {
        legend_label <- "ABC-SMC"
        posterior_df <- data.frame(
            x = abc_results[["Iteration_2"]]$parameters[[parameters_labels$parameter[1]]],
            y = abc_results[["Iteration_2"]]$parameters[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    } else if (method == "abc-mcmc") {
        legend_label <- "ABC-MCMC"
        posterior_df <- data.frame(
            x = abc_results[["Iteration_2"]]$parameters[[parameters_labels$parameter[1]]],
            y = abc_results[["Iteration_2"]]$parameters[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    }
    if (!is.null(lims)) posterior_df <- apply_lims(posterior_df)
    #---Plot joint distribution
    plots <- plots +
        geom_density_2d(data = posterior_df, aes(x = x, y = y, color = legend), linewidth = 3, bins = nBins)
    #---Add label for parameter
    if (new_plot == TRUE) {
        if ("label" %in% colnames(parameters_labels)) {
            plots <- plots + labs(x = eval(parse(text = parameters_labels$label[1])), y = eval(parse(text = parameters_labels$label[2])))
        } else {
            plots <- plots + labs(x = parameters_labels$parameter[1], y = parameters_labels$parameter[2])
        }
    }
    #---Beautify plot
    plots <- plots +
        scale_color_manual(values = color_scheme, name = "", breaks = legend_order) +
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
    file_name <- paste0("comparison-joint-parameters=", parameters_labels$parameter[1], "-vs-", parameters_labels$parameter[2], ".png")
    png(file_name, res = 150, width = 30, height = 31, units = "in", pointsize = 12)
    print(plots)
    dev.off()
    return(plots)
}
