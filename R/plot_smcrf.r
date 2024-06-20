plot_smcrf_marginal <- function(smcrf_results,
                                parameters_truth = NULL,
                                parameters_labels = NULL,
                                statistics_labels = NULL,
                                plot_statistics = FALSE,
                                xlimit = NULL,
                                alpha = 0.3,
                                plot_hist = FALSE) {
    nIterations <- smcrf_results[["nIterations"]]
    if (is.null(parameters_labels)) parameters_labels <- smcrf_results[["parameters_labels"]]
    if (is.null(statistics_labels)) statistics_labels <- smcrf_results[["statistics_labels"]]
    #---Set up color scheme for plotting
    color_scheme <- c(
        "True Posterior Distribution" = "black", "Data statistic" = "black",
        "Prior Distribution" = "gray",
        setNames(rev(rainbow(nIterations)), paste0("Iter. ", 1:nIterations))
    )
    #---Set up legend order for plotting
    legend_order <- c("True Posterior Distribution", "Data statistic", "Prior Distribution", paste0("Iter. ", 1:nIterations))
    #---Plot marginal distributions for each parameter
    for (i in 1:length(parameters_labels$parameter)) {
        parameter_id <- parameters_labels$parameter[i]
        #   Begin plot
        p <- ggplot()
        #   Plot True Posterior distribution (if provided)
        if (!is.null(parameters_truth)) {
            true_posterior_df <- data.frame(value = parameters_truth[[parameter_id]], legend = "True Posterior Distribution")
            if (nrow(parameters_truth) == 1) {
                p <- p + geom_vline(data = true_posterior_df, aes(xintercept = value), linetype = "solid", linewidth = 5)
            } else {
                if (plot_hist) {
                    p <- p + geom_histogram(data = true_posterior_df, aes(x = value, y = ..density.., fill = legend, color = legend), alpha = 1)
                } else {
                    p <- p + geom_density(data = true_posterior_df, aes(x = value, fill = legend, color = legend), alpha = 1, linewidth = 2)
                }
            }
        }
        #   Plot prior distribution
        # prior_df <- data.frame(value = smcrf_results[["Iteration_1"]]$parameters[[parameter_id]], legend = "Prior Distribution")
        prior_df <- data.frame(value = smcrf_results[["Iteration_1"]]$parameters_unperturbed[[parameter_id]], legend = "Prior Distribution")
        if (plot_hist) {
            p <- p + geom_histogram(data = prior_df, aes(x = value, y = ..density.., fill = legend, color = legend), alpha = alpha)
        } else {
            p <- p + geom_density(data = prior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
        }
        #   Plot posterior distribution for each iteration
        for (iteration in 1:nIterations) {
            posterior_df <- data.frame(
                # value = smcrf_results[[paste0("Iteration_", iteration + 1)]]$parameters[[parameter_id]],
                value = smcrf_results[[paste0("Iteration_", iteration + 1)]]$parameters_unperturbed[[parameter_id]],
                legend = paste0("Iter. ", iteration)
            )
            if (plot_hist) {
                p <- p + geom_histogram(data = posterior_df, aes(x = value, y = ..density.., fill = legend, color = legend), alpha = alpha)
            } else {
                p <- p + geom_density(data = posterior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
            }
            if (!is.null(xlimit)) {
                p <- p +
                    xlim(xlimit$min[which(xlimit$parameter == parameter_id)], xlimit$max[which(xlimit$parameter == parameter_id)])
            }
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
            if (plot_hist) {
                p <- p + geom_histogram(data = prior_df, aes(x = value, y = ..density.., fill = legend, color = legend), alpha = alpha)
            } else {
                p <- p + geom_density(data = prior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
            }

            #   Plot posterior distribution for each iteration
            for (iteration in 1:nIterations) {
                posterior_df <- data.frame(
                    value = smcrf_results[[paste0("Iteration_", iteration + 1)]]$statistics[[statistic_id]],
                    legend = paste0("Iter. ", iteration)
                )
                if (plot_hist) {
                    p <- p + geom_histogram(data = posterior_df, aes(x = value, y = ..density.., fill = legend, color = legend), alpha = alpha)
                } else {
                    p <- p + geom_density(data = posterior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
                }
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
                                  plot_hist = FALSE,
                                  plot_hist_point = FALSE,
                                  alpha = 0.3,
                                  plot_prior = FALSE) {
    if (is.null(parameters_labels)) parameters_labels <- abc_results[["parameters_labels"]]
    if (is.null(statistics_labels)) statistics_labels <- abc_results[["statistics_labels"]]
    method <- abc_results[["method"]]
    #---Set up color scheme for plotting
    color_scheme <- c(
        "Prior Distribution" = "gray",
        "True Posterior Distribution" = "black",
        "ABC-rejection" = "forestgreen",
        "ABC-RF" = "magenta4",
        "DRF" = "royalblue2",
        "MCMC" = "goldenrod2",
        "ABC-MCMC" = "goldenrod2",
        "ABC-SMC" = "goldenrod2",
        "SMC-RF for single parameters" = "salmon",
        "SMC-RF for multiple parameters" = "salmon"
    )
    #---Set up legend order for plotting
    legend_order <- c(
        "Prior Distribution",
        "True Posterior Distribution",
        "ABC-rejection",
        "ABC-MCMC",
        "ABC-SMC",
        "ABC-RF",
        "DRF",
        "MCMC",
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
    #---Plot True Posterior Distribution parameter distributions (if provided)
    if (!is.null(parameters_truth)) {
        for (parameter_id in parameters_labels$parameter) {
            if (any(grepl("density", colnames(parameters_truth)))) {
                true_posterior_df <- data.frame(value = parameters_truth[[parameter_id]], density = parameters_truth[["density"]], legend = "True Posterior Distribution")
                if (plot_hist || plot_hist_point) {
                    plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                        geom_histogram(data = true_posterior_df, aes(x = value, weight = density, y = ..density.., fill = legend, color = legend), alpha = 0.5)
                } else {
                    plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                        geom_line(data = true_posterior_df, aes(x = value, y = density, color = legend), alpha = 0.5, linetype = "solid", linewidth = 5) +
                        geom_area(alpha = 0.5, fill = legend)
                }
            } else {
                true_posterior_df <- data.frame(value = parameters_truth[[parameter_id]], legend = "True Posterior Distribution")
                if (nrow(parameters_truth) == 1) {
                    plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                        geom_vline(data = true_posterior_df, aes(xintercept = value), linetype = "solid", linewidth = 5)
                } else {
                    if (plot_hist || plot_hist_point) {
                        plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                            geom_histogram(data = true_posterior_df, aes(x = value, y = ..density.., fill = legend, color = legend), alpha = 0.5)
                    } else {
                        plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                            geom_density(data = true_posterior_df, aes(x = value, fill = legend, color = legend), alpha = 0.5, linewidth = 2)
                    }
                }
            }
        }
    }
    #---Plot prior distributions (if provided)
    if (plot_prior) {
        for (parameter_id in parameters_labels$parameter) {
            prior_df <- data.frame(value = abc_results[["Iteration_1"]]$parameters[[parameter_id]], legend = "Prior Distribution")
            if (plot_hist) {
                plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                    geom_histogram(data = prior_df, aes(x = value, y = ..density.., fill = legend, color = legend), alpha = alpha)
            } else {
                plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                    geom_density(data = prior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
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
        # parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters
        parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed
        if (!is.null(sample_num)) {
            # parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters[1:sample_num, , drop = FALSE]
            parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed[1:sample_num, , drop = FALSE]
        }
        if (plot_statistics) statistics_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$statistics
    } else if (method == "smcrf-multi-param") {
        nIterations <- abc_results[["nIterations"]]
        if (nIterations == 1) {
            legend_label <- "DRF"
        } else {
            legend_label <- "SMC-RF for multiple parameters"
        }
        # parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters
        parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed
        if (plot_statistics) statistics_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$statistics
    } else if (method == "abc-rejection") {
        legend_label <- "ABC-rejection"
        # parameters_values <- abc_results[["Iteration_2"]]$parameters
        parameters_values <- abc_results[["Iteration_2"]]$parameters_unperturbed
        if (plot_statistics) statistics_values <- abc_results[["Iteration_2"]]$statistics
    } else if (method == "abc-smc") {
        legend_label <- "ABC-SMC"
        # parameters_values <- abc_results[["Iteration_2"]]$parameters
        parameters_values <- abc_results[["Iteration_2"]]$parameters_unperturbed
        if (plot_statistics) statistics_values <- abc_results[["Iteration_2"]]$statistics
    } else if (method == "abc-mcmc") {
        legend_label <- "ABC-MCMC"
        # parameters_values <- abc_results[["Iteration_2"]]$parameters
        parameters_values <- abc_results[["Iteration_2"]]$parameters_unperturbed
        if (plot_statistics) statistics_values <- abc_results[["Iteration_2"]]$statistics
    } else if (method == "mcmc") {
        legend_label <- "MCMC"
        parameters_values <- abc_results$selected_params
        if (plot_statistics) statistics_values <- abc_results$statistics
    }
    #---Plot marginal distributions for each parameter
    for (parameter_id in parameters_labels$parameter) {
        posterior_df <- data.frame(
            value = parameters_values[[parameter_id]],
            legend = legend_label
        )
        if (plot_hist) {
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                geom_histogram(data = posterior_df, aes(x = value, y = ..density.., fill = legend, color = legend), alpha = alpha)
        } else if (plot_hist_point) {
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                stat_bin(data = posterior_df, aes(x = value, y = ..density.., fill = legend, color = legend), geom = "point", bins = 30, position = "identity", size = 10) +
                stat_bin(data = posterior_df, aes(x = value, y = ..density.., fill = legend, color = legend), geom = "line", bins = 30, position = "identity", size = 3, show.legend = FALSE)
        } else {
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                geom_density(data = posterior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
        }
        if (!is.null(xlimit)) {
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                xlim(xlimit$min[which(xlimit$parameter == parameter_id)], xlimit$max[which(xlimit$parameter == parameter_id)])
        }
    }
    if (plot_statistics) {
        for (statistic_id in statistics_labels$ID) {
            posterior_df <- data.frame(
                value = statistics_values[[statistic_id]],
                legend = legend_label
            )
            if (plot_hist) {
                plots$statistics[[statistic_id]] <- plots$statistics[[statistic_id]] +
                    geom_histogram(data = posterior_df, aes(x = value, y = ..density.., fill = legend, color = legend), alpha = alpha)
            } else if (plot_hist_point) {
                plots$statistics[[statistic_id]] <- plots$statistics[[statistic_id]] +
                    stat_bin(data = posterior_df, aes(x = value, y = ..density.., color = legend), geom = "point", bins = 30, position = "identity", size = 10) +
                    stat_bin(data = posterior_df, aes(x = value, y = ..density.., color = legend), geom = "line", bins = 30, position = "identity", size = 3)
            } else {
                plots$statistics[[statistic_id]] <- plots$statistics[[statistic_id]] +
                    geom_density(data = posterior_df, aes(x = value, fill = legend, color = legend), alpha = alpha, linewidth = 2)
            }
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
    legend_order <- c("True Posterior Distribution", "Prior Distribution", paste0("Iter. ", 1:nIterations))
    #---Function to apply limits to data (if provided)
    # apply_lims <- function(df) {
    #     df <- df %>%
    #         filter(
    #             x >= lims$min[which(lims$ID == parameters_labels$parameter[1])],
    #             x <= lims$max[which(lims$ID == parameters_labels$parameter[1])],
    #             y >= lims$min[which(lims$ID == parameters_labels$parameter[2])],
    #             y <= lims$max[which(lims$ID == parameters_labels$parameter[2])]
    #         )
    # }
    #---Begin plot
    p <- ggplot()
    #---Plot True Posterior distribution (if provided)
    if (!is.null(parameters_truth)) {
        true_posterior_df <- data.frame(
            x = parameters_truth[[parameters_labels$parameter[1]]],
            y = parameters_truth[[parameters_labels$parameter[2]]]
        )
        # if (!is.null(lims)) true_posterior_df <- apply_lims(true_posterior_df)
        p <- p + geom_density_2d_filled(data = true_posterior_df, aes(x = x, y = y), show.legend = FALSE)
    }
    #---Plot prior distribution
    # prior_df <- data.frame(
    #     x = smcrf_results[["Iteration_1"]]$parameters[[parameters_labels$parameter[1]]],
    #     y = smcrf_results[["Iteration_1"]]$parameters[[parameters_labels$parameter[2]]],
    #     legend = "Prior Distribution"
    # )
    prior_df <- data.frame(
        x = smcrf_results[["Iteration_1"]]$parameters_unperturbed[[parameters_labels$parameter[1]]],
        y = smcrf_results[["Iteration_1"]]$parameters_unperturbed[[parameters_labels$parameter[2]]],
        legend = "Prior Distribution"
    )
    # if (!is.null(lims)) prior_df <- apply_lims(prior_df)
    p <- p + geom_density_2d(data = prior_df, aes(x = x, y = y, color = legend), linewidth = 3, bins = nBins)
    #---Plot posterior distribution for each iteration
    for (iteration in 1:nIterations) {
        # posterior_df <- data.frame(
        #     x = smcrf_results[[paste0("Iteration_", iteration + 1)]]$parameters[[parameters_labels$parameter[1]]],
        #     y = smcrf_results[[paste0("Iteration_", iteration + 1)]]$parameters[[parameters_labels$parameter[2]]],
        #     legend = paste0("Iter. ", iteration)
        # )
        posterior_df <- data.frame(
            x = smcrf_results[[paste0("Iteration_", iteration + 1)]]$parameters_unperturbed[[parameters_labels$parameter[1]]],
            y = smcrf_results[[paste0("Iteration_", iteration + 1)]]$parameters_unperturbed[[parameters_labels$parameter[2]]],
            legend = paste0("Iter. ", iteration)
        )
        # if (!is.null(lims)) posterior_df <- apply_lims(posterior_df)
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
    if (!is.null(lims)) {
        p <- p +
            xlim(c(lims$min[which(lims$parameter == parameters_labels$parameter[1])], lims$max[which(lims$parameter == parameters_labels$parameter[1])])) +
            ylim(c(lims$min[which(lims$parameter == parameters_labels$parameter[2])], lims$max[which(lims$parameter == parameters_labels$parameter[2])]))
    }
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
    library(RColorBrewer)
    method <- abc_results[["method"]]

    if (is.null(parameters_labels)) parameters_labels <- abc_results[["parameters_labels"]]
    if (nrow(parameters_labels) != 2) stop("ERROR: plot_compare_joint only works for two parameters. Please check parameters_labels.")
    #---Set up color scheme for plotting
    # color_scheme <- c(
    #     "Prior Distribution" = "gray",
    #     "True Posterior Distribution" = "black",
    #     "ABC-rejection" = "forestgreen",
    #     "ABC-RF" = "royalblue2",
    #     "DRF" = "royalblue2",
    #     "MCMC" = "goldenrod2",
    #     "ABC-MCMC" = "goldenrod2",
    #     "ABC-SMC" = "magenta4",
    #     "SMC-RF for single parameters" = "salmon",
    #     "SMC-RF for multiple parameters" = "salmon"
    # )

    color_scheme <- c(
        "ABC-rejection" = "forestgreen",
        "ABC-RF" = "magenta4",
        # "DRF" = "burlywood2",
        "DRF" = "cyan1",
        "MCMC" = "khaki",
        "ABC-MCMC" = "khaki",
        "ABC-SMC" = "goldenrod2",
        "SMC-RF for single parameters" = "firebrick1",
        # "SMC-RF for multiple parameters" = "plum1"
        "SMC-RF for multiple parameters" = "firebrick1"
    )
    #---Set up legend order for plotting
    legend_order <- c(
        "ABC-rejection",
        "MCMC",
        "ABC-MCMC",
        "ABC-SMC",
        "ABC-RF",
        "DRF",
        "SMC-RF for single parameters",
        "SMC-RF for multiple parameters"
    )
    my_palette <- colorRampPalette(c("#170756ad", "#0b50a4", "#0beac8", "#ffff0f"))

    #---Function to apply limits to data (if provided)
    # apply_lims <- function(df) {
    #     df <- df %>%
    #         filter(
    #             x >= lims$min[which(lims$ID == parameters_labels$parameter[1])],
    #             x <= lims$max[which(lims$ID == parameters_labels$parameter[1])],
    #             y >= lims$min[which(lims$ID == parameters_labels$parameter[2])],
    #             y <= lims$max[which(lims$ID == parameters_labels$parameter[2])]
    #         )
    # }
    #---Begin plots
    if (is.null(plots)) {
        plots <- ggplot()
        new_plot <- TRUE
    } else {
        new_plot <- FALSE
    }
    #---Plot True Posterior parameter distributions (if provided)
    if (!is.null(parameters_truth)) {
        true_posterior_df <- data.frame(
            x = parameters_truth[[parameters_labels$parameter[1]]],
            y = parameters_truth[[parameters_labels$parameter[2]]]
        )
        # if (!is.null(lims)) true_posterior_df <- apply_lims(true_posterior_df)
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
        # posterior_df <- data.frame(
        #     x = abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters[[parameters_labels$parameter[1]]],
        #     y = abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters[[parameters_labels$parameter[2]]],
        #     legend = legend_label
        # )
        posterior_df <- data.frame(
            x = abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed[[parameters_labels$parameter[1]]],
            y = abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    } else if (method == "smcrf-multi-param") {
        nIterations <- abc_results[["nIterations"]]
        if (nIterations == 1) {
            legend_label <- "DRF"
        } else {
            legend_label <- "SMC-RF for multiple parameters"
        }
        # posterior_df <- data.frame(
        #     x = abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters[[parameters_labels$parameter[1]]],
        #     y = abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters[[parameters_labels$parameter[2]]],
        #     legend = legend_label
        # )
        posterior_df <- data.frame(
            x = abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed[[parameters_labels$parameter[1]]],
            y = abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    } else if (method == "abc-rejection") {
        legend_label <- "ABC-rejection"
        # posterior_df <- data.frame(
        #     x = abc_results[["Iteration_2"]]$parameters[[parameters_labels$parameter[1]]],
        #     y = abc_results[["Iteration_2"]]$parameters[[parameters_labels$parameter[2]]],
        #     legend = legend_label
        # )
        posterior_df <- data.frame(
            x = abc_results[["Iteration_2"]]$parameters_unperturbed[[parameters_labels$parameter[1]]],
            y = abc_results[["Iteration_2"]]$parameters_unperturbed[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    } else if (method == "abc-smc") {
        legend_label <- "ABC-SMC"
        # posterior_df <- data.frame(
        #     x = abc_results[["Iteration_2"]]$parameters[[parameters_labels$parameter[1]]],
        #     y = abc_results[["Iteration_2"]]$parameters[[parameters_labels$parameter[2]]],
        #     legend = legend_label
        # )
        posterior_df <- data.frame(
            x = abc_results[["Iteration_2"]]$parameters_unperturbed[[parameters_labels$parameter[1]]],
            y = abc_results[["Iteration_2"]]$parameters_unperturbed[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    } else if (method == "abc-mcmc") {
        legend_label <- "ABC-MCMC"
        # posterior_df <- data.frame(
        #     x = abc_results[["Iteration_2"]]$parameters[[parameters_labels$parameter[1]]],
        #     y = abc_results[["Iteration_2"]]$parameters[[parameters_labels$parameter[2]]],
        #     legend = legend_label
        # )
        posterior_df <- data.frame(
            x = abc_results[["Iteration_2"]]$parameters_unperturbed[[parameters_labels$parameter[1]]],
            y = abc_results[["Iteration_2"]]$parameters_unperturbed[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    } else if (method == "mcmc") {
        legend_label <- "MCMC"
        posterior_df <- data.frame(
            x = abc_results$selected_params[[parameters_labels$parameter[1]]],
            y = abc_results$selected_params[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    } else if (method == "true-joint") {
        legend_label <- "TRUE-JOINT"
        posterior_df <- data.frame(
            x = abc_results$parameters_truth[[parameters_labels$parameter[1]]],
            y = abc_results$parameters_truth[[parameters_labels$parameter[2]]],
            legend = legend_label
        )
    }
    # if (!is.null(lims)) {
    #     posterior_df <- apply_lims(posterior_df)
    # }
    #---Plot joint distribution
    if (method == "true-joint") {
        plots <- plots +
            geom_density_2d_filled(data = posterior_df, aes(x = x, y = y), show.legend = FALSE)
    }
    # else if (method == "mcmc") {
    #     plots <- plots +
    #         geom_density_2d_filled(data = posterior_df, aes(x = x, y = y), show.legend = FALSE)
    # }
    else {
        plots <- plots +
            geom_density_2d(data = posterior_df, aes(x = x, y = y, color = legend), linewidth = 3, bins = 3)
    }
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
    if (!is.null(lims)) {
        plots <- plots +
            xlim(c(lims$min[which(lims$parameter == parameters_labels$parameter[1])], lims$max[which(lims$parameter == parameters_labels$parameter[1])])) +
            ylim(c(lims$min[which(lims$parameter == parameters_labels$parameter[2])], lims$max[which(lims$parameter == parameters_labels$parameter[2])]))
    }
    #---Print joint distribution plot
    file_name <- paste0("comparison-joint-parameters=", parameters_labels$parameter[1], "-vs-", parameters_labels$parameter[2], ".png")
    png(file_name, res = 150, width = 30, height = 31, units = "in", pointsize = 12)
    print(plots)
    dev.off()
    return(plots)
}

plot_compare_qqplot <- function(plots = NULL,
                                abc_results,
                                parameters_truth = NULL,
                                parameters_labels = NULL,
                                statistics_labels = NULL,
                                sample_num = NULL,
                                plot_statistics = FALSE,
                                plot_hist = FALSE,
                                lims = NULL,
                                alpha = 0.3,
                                plot_prior = FALSE) {
    if (is.null(parameters_labels)) parameters_labels <- abc_results[["parameters_labels"]]
    method <- abc_results[["method"]]
    #---Set up color scheme for plotting
    color_scheme <- c(
        "ABC-rejection" = "forestgreen",
        "ABC-RF" = "magenta4",
        "DRF" = "royalblue2",
        "MCMC" = "goldenrod2",
        "ABC-MCMC" = "goldenrod2",
        "ABC-SMC" = "goldenrod2",
        "SMC-RF for single parameters" = "salmon",
        "SMC-RF for multiple parameters" = "salmon"
    )
    #---Set up legend order for plotting
    legend_order <- c(
        "ABC-rejection",
        "ABC-MCMC",
        "ABC-SMC",
        "ABC-RF",
        "DRF",
        "MCMC",
        "SMC-RF for single parameters",
        "SMC-RF for multiple parameters"
    )
    #---Begin plots
    if (is.null(plots)) {
        plots <- list()
        for (parameter_id in parameters_labels$parameter) {
            plots$parameters[[parameter_id]] <- ggplot()
        }
        new_plot <- TRUE
    } else {
        new_plot <- FALSE
    }
    #---Extract final posterior distributions
    if (method == "smcrf-single-param") {
        nIterations <- abc_results[["nIterations"]]
        if (nIterations == 1) {
            legend_label <- "ABC-RF"
        } else {
            legend_label <- "SMC-RF for single parameters"
        }
        # parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters
        parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed
        if (!is.null(sample_num)) {
            # parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters[1:sample_num, , drop = FALSE]
            parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed[1:sample_num, , drop = FALSE]
        }
    } else if (method == "smcrf-multi-param") {
        nIterations <- abc_results[["nIterations"]]
        if (nIterations == 1) {
            legend_label <- "DRF"
        } else {
            legend_label <- "SMC-RF for multiple parameters"
        }
        # parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters
        parameters_values <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed
    } else if (method == "abc-rejection") {
        legend_label <- "ABC-rejection"
        # parameters_values <- abc_results[["Iteration_2"]]$parameters
        parameters_values <- abc_results[["Iteration_2"]]$parameters_unperturbed
    } else if (method == "abc-smc") {
        legend_label <- "ABC-SMC"
        # parameters_values <- abc_results[["Iteration_2"]]$parameters
        parameters_values <- abc_results[["Iteration_2"]]$parameters_unperturbed
    } else if (method == "abc-mcmc") {
        legend_label <- "ABC-MCMC"
        # parameters_values <- abc_results[["Iteration_2"]]$parameters
        parameters_values <- abc_results[["Iteration_2"]]$parameters_unperturbed
    } else if (method == "mcmc") {
        legend_label <- "MCMC"
        parameters_values <- abc_results$selected_params
    }

    if (new_plot) {
        for (parameter_id in parameters_labels$parameter) {
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                geom_abline(intercept = 0, slope = 1, color = "black", size = 1.5)
        }
    }

    for (parameter_id in parameters_labels$parameter) {
        true_para <- qqplot(x = parameters_truth[[parameter_id]], y = parameters_values[[parameter_id]], plot.it = FALSE)$x
        fitted_para <- qqplot(x = parameters_truth[[parameter_id]], y = parameters_values[[parameter_id]], plot.it = FALSE)$y
        qq_data <- data.frame(True_parameter = true_para, Fitted_parameter = fitted_para, legend = legend_label)
        plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
            geom_point(data = qq_data, aes(x = True_parameter, y = Fitted_parameter, color = legend), size = 10)
    }
    if (new_plot == TRUE) {
        if ("label" %in% colnames(parameters_labels)) {
            for (parameter_id in parameters_labels$parameter) {
                label_for_id <- parameters_labels$label[parameters_labels$parameter == parameter_id]
                label_string <- gsub("expression\\((.*)\\)", "\\1", label_for_id)
                label_expression <- parse(text = label_string)[[1]]
                plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                    labs(
                        x = bquote("Likelihood-based " ~ .(label_expression)),
                        y = bquote("Inferred " ~ .(label_expression)),
                    )
            }
        } else {
            for (parameter_id in parameters_labels$parameter) {
                plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                    labs(
                        x = paste0(
                            "Likelihood ", parameter_id
                        ),
                        y = paste0(
                            "Inferred ", parameter_id
                        )
                    )
            }
        }
    }
    #---Beautify plots
    if (new_plot) {
        for (parameter_id in parameters_labels$parameter) {
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                # scale_fill_manual(values = color_scheme, name = "", breaks = legend_order) +
                scale_color_manual(values = color_scheme, name = "", breaks = legend_order) +
                guides(color = guide_legend(override.aes = list(size = 10))) +
                # guides(fill = guide_legend(nrow = 1, keywidth = 2.5, keyheight = 1)) +
                theme(
                    text = element_text(size = 50),
                    panel.background = element_rect(fill = "white", colour = "white"),
                    panel.grid.major = element_line(colour = "white"),
                    panel.grid.minor = element_line(colour = "white"),
                    legend.key.size = unit(10, "points"),
                    legend.position = "top",
                    legend.justification = c(0, 0.5)
                )
        }
    }
    if (!is.null(lims)) {
        for (parameter_id in parameters_labels$parameter) {
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                xlim(c(lims$min[which(lims$parameter == parameter_id)], lims$max[which(lims$parameter == parameter_id)])) +
                ylim(c(lims$min[which(lims$parameter == parameter_id)], lims$max[which(lims$parameter == parameter_id)]))
        }
    }
    #---Print marginal distribution plots
    for (parameter_id in parameters_labels$parameter) {
        file_name <- paste0("qq-plot-compare-parameter=", parameter_id, ".png")
        png(file_name, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
        print(plots$parameters[[parameter_id]])
        dev.off()
    }
    return(plots)
}
