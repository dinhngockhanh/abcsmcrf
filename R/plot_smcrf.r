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
            p <- p + geom_density(data = true_posterior_df, aes(x = value, fill = legend, color = legend), linewidth = 2)
        }
        #   Plot prior distribution
        prior_df <- data.frame(value = smcrf_results[["Iteration_1"]]$parameters[[parameter_id]], legend = "Prior Distribution")
        p <- p + geom_density(data = prior_df, aes(x = value, fill = legend, color = legend), alpha = 0.2, linewidth = 2)
        #   Plot posterior distribution for each iteration
        for (iteration in 1:nIterations) {
            posterior_df <- data.frame(
                value = smcrf_results[[paste0("Iteration_", iteration)]]$parameters[[parameter_id]],
                weight = smcrf_results[[paste0("Iteration_", iteration)]]$weights[[parameter_id]],
                legend = paste0("Iter. ", iteration)
            )
            p <- p + geom_density(data = posterior_df, aes(x = value, weight = weight, fill = legend, color = legend), alpha = alpha, linewidth = 2)
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
        file_name <- paste0(smcrf_results[["method"]], "-marginal-", parameter_id, ".png")
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
            data_df <- data.frame(value = smcrf_results$statistics_target[[statistic_id]])
            p <- p + geom_vline(data = data_df, aes(xintercept = value), linetype = "solid", linewidth = 5)
            #   Plot prior distribution
            prior_df <- data.frame(value = smcrf_results[["Iteration_1"]]$statistics[[statistic_id]], legend = "Prior Distribution")
            p <- p + geom_density(data = prior_df, aes(x = value, fill = legend, color = legend), alpha = 0.2, linewidth = 2)
            #   Plot posterior distribution for each iteration
            for (iteration in 2:nIterations) {
                posterior_df <- data.frame(
                    value = smcrf_results[[paste0("Iteration_", iteration)]]$statistics[[statistic_id]],
                    legend = paste0("Iter. ", iteration - 1)
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
            file_name <- paste0(smcrf_results[["method"]], "-marginal-", statistic_id, ".png")
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
                                  plot_statistics = FALSE,
                                  alpha = 0.3) {
    if (is.null(parameters_labels)) parameters_labels <- abc_results[["parameters_labels"]]
    if (is.null(statistics_labels)) statistics_labels <- abc_results[["statistics_labels"]]
    method <- abc_results[["method"]]
    #---Set up color scheme for plotting
    color_scheme <- c(
        "True Posterior" = "black",
        "SMC-RF for single parameters" = "salmon",
        "SMC-RF for multiple parameters" = "royalblue4"
    )
    #---Set up legend order for plotting
    legend_order <- c(
        "True Posterior",
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
            plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
                geom_density(data = true_posterior_df, aes(x = value, fill = legend, color = legend), linewidth = 2)
        }
    }

    #---Extract final posterior distributions
    if (method == "smcrf-single-param") {
        legend_label <- "SMC-RF for single parameters"
        nIterations <- abc_results[["nIterations"]]
        parameters_values <- abc_results[[paste0("Iteration_", nIterations)]]$parameters
        parameters_weights <- abc_results[[paste0("Iteration_", nIterations)]]$weights
        if (plot_statistics) statistics_values <- abc_results[[paste0("Iteration_", nIterations)]]$statistics
    } else if (method == "smcrf-multi-param") {
        legend_label <- "SMC-RF for multiple parameters"
        nIterations <- abc_results[["nIterations"]]
        parameters_values <- abc_results[[paste0("Iteration_", nIterations)]]$parameters
        parameters_weights <- abc_results[[paste0("Iteration_", nIterations)]]$weights
        if (plot_statistics) statistics_values <- abc_results[[paste0("Iteration_", nIterations)]]$statistics
    }
    #---Plot marginal distributions for each parameter
    for (parameter_id in parameters_labels$parameter) {
        posterior_df <- data.frame(
            value = parameters_values[[parameter_id]],
            weight = parameters_weights[[parameter_id]],
            legend = legend_label
        )
        plots$parameters[[parameter_id]] <- plots$parameters[[parameter_id]] +
            geom_density(data = posterior_df, aes(x = value, weight = weight, fill = legend, color = legend), alpha = alpha, linewidth = 2)
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
    if ("label" %in% colnames(parameters_labels) & new_plot == TRUE) {
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
        file_name <- paste0("Comparison-marginal-", parameter_id, ".png")
        png(file_name, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
        print(plots$parameters[[parameter_id]])
        dev.off()
    }
    if (plot_statistics) {
        for (statistic_id in statistics_labels$ID) {
            file_name <- paste0("Comparison-marginal-", statistic_id, ".png")
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
    for (iteration in 2:nIterations) {
        posterior_df <- data.frame(
            x = smcrf_results[[paste0("Iteration_", iteration)]]$parameters[[parameters_labels$parameter[1]]],
            y = smcrf_results[[paste0("Iteration_", iteration)]]$parameters[[parameters_labels$parameter[2]]],
            legend = paste0("Iter. ", iteration - 1)
        )
        if (!is.null(lims)) posterior_df <- apply_lims(posterior_df)
        p <- p + geom_density_2d(data = posterior_df, aes(x = x, y = y, color = legend), linewidth = 2, bins = nBins)
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
    file_name <- paste0(smcrf_results[["method"]], "-joint-", parameters_labels$parameter[1], "-", parameters_labels$parameter[2], ".png")
    png(file_name, res = 150, width = 30, height = 31, units = "in", pointsize = 12)
    print(p)
    dev.off()
}
