abc_rejection <- function(statistics_target,
                          model,
                          parameters_initial,
                          ...) {
    library(abc)
    parameters_ids <- colnames(parameters_initial)
    cat(paste0("\n\nABC-REJECTION FOR PARAMETERS: ...\n"))
    prior_parameters <- data.frame(parameters_initial)
    colnames(prior_parameters) <- parameters_ids
    #---Simulate statistics
    cat("Making model simulations...\n")
    reference <- model(parameters = prior_parameters)
    statistics <- data.frame(reference[, colnames(reference)[!colnames(reference) %in% parameters_ids]])
    colnames(statistics) <- colnames(reference)[!colnames(reference) %in% parameters_ids]
    rej_model <- abc(
        target = statistics_target,
        param = prior_parameters,
        sumstat = statistics,
        tol = .1,
        method = "rejection"
    )
    post_parameters <- data.frame(rej_model$unadj.values)
    colnames(post_parameters) <- parameters_ids
    #---Save ABC-REJ results
    ABCREJ <- list()
    ABCREJ[["prior_parameters"]] <- prior_parameters
    ABCREJ[["reference"]] <- reference
    ABCREJ[["statistics"]] <- statistics
    ABCREJ[["rej_model"]] <- rej_model
    ABCREJ[["post_parameters"]] <- post_parameters
    ABCREJ[["method"]] <- "ABC_rejection"
    ABCREJ[["statistics_target"]] <- statistics_target
    ABCREJ[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
    ABCREJ[["statistics_labels"]] <- data.frame(ID = colnames(reference)[!colnames(reference) %in% parameters_ids])
    return(ABCREJ)
}

plot_abc_marginal <- function(abc_results,
                              parameters_truth = NULL,
                              parameters_labels = NULL,
                              statistics_labels = NULL,
                              plot_statistics = FALSE,
                              alpha = 0.3) {
    if (is.null(parameters_labels)) parameters_labels <- abc_results[["parameters_labels"]]
    if (is.null(statistics_labels)) statistics_labels <- abc_results[["statistics_labels"]]
    #---Set up color scheme for plotting
    color_scheme <- c(
        "True Posterior" = "black", "Data statistic" = "black",
        "Prior Distribution" = "gray",
        "Posterior Distribution" = "red"
    )
    legend_order <- c("True Posterior", "Data statistic", "Prior Distribution", "Posterior Distribution")
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
        prior_df <- data.frame(value = abc_results$prior_parameters[[parameter_id]], legend = "Prior Distribution")
        p <- p + geom_density(data = prior_df, aes(x = value, fill = legend, color = legend), alpha = 0.2, linewidth = 2)
        #   Plot posterior distribution
        post_df <- data.frame(value = abc_results$post_parameters[[parameter_id]], legend = "Posterior Distribution")
        p <- p + geom_density(data = post_df, aes(x = value, fill = legend, color = legend), alpha = 0.2, linewidth = 2)
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
        file_name <- paste0(abc_results[["method"]], "-marginal-", parameter_id, ".png")
        png(file_name, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
        print(p)
        dev.off()
    }
}
