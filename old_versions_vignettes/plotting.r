densityPlot_df <- function(object,
                           obs,
                           training,
                           add = TRUE,
                           main = "Posterior density",
                           color_prior,
                           chosen_para = NULL,
                           color_posterior,
                           protocol = "",
                           color_vline,
                           log = "",
                           xlim = NULL,
                           ylim = NULL,
                           xlab = NULL,
                           ylab = NULL,
                           paral = FALSE,
                           cutbound = FALSE,
                           lower = NULL,
                           upper = NULL,
                           ncores = if (paral) max(detectCores() - 1, 1) else 1, ...) {
    findweights <- getFromNamespace("findweights", "abcrf")
    ### Checking arguments
    if (!inherits(object, "regAbcrf")) {
        stop("object not of class regAbcrf")
    }

    if (!inherits(training, "data.frame")) {
        stop("training needs to be a data.frame object")
    }

    if (!inherits(obs, "data.frame")) {
        stop("obs needs to be a data.frame object")
    }
    if (nrow(obs) == 0L || is.null(nrow(obs))) {
        stop("no data in obs")
    }
    if (nrow(training) == 0L || is.null(nrow(training))) {
        stop("no simulation in the training reference table (response, sumstat)")
    }

    if ((!is.logical(add)) || (length(add) != 1L)) {
        stop("add should be TRUE or FALSE")
    }
    if ((!is.logical(paral)) || (length(paral) != 1L)) {
        stop("paral should be TRUE or FALSE")
    }
    if (is.na(ncores)) {
        warning("Unable to automatically detect the number of CPU cores, \n1 CPU core will be used or please specify ncores.")
        ncores <- 1
    }

    if (!is.character(log)) {
        stop("log needs to be a character string")
    }
    x <- obs
    if (!is.null(x)) {
        if (is.vector(x)) {
            x <- matrix(x, ncol = 1)
        }
        if (nrow(x) == 0) {
            stop("obs has 0 rows")
        }
        if (any(is.na(x))) {
            stop("missing values in obs")
        }
    }

    # resp and sumsta recover

    mf <- match.call(expand.dots = FALSE)
    mf <- mf[1]
    mf$formula <- object$formula


    mf$data <- training


    training <- mf$data

    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    resp <- model.response(mf)

    obj <- object$model.rf
    inbag <- matrix(unlist(obj$inbag.counts, use.names = FALSE), ncol = obj$num.trees, byrow = FALSE)

    obj[["origNodes"]] <- predict(obj, training, predict.all = TRUE, num.threads = ncores)$predictions
    obj[["origObs"]] <- model.response(mf)

    #####################

    origObs <- obj$origObs
    origNodes <- obj$origNodes

    nodes <- predict(obj, x, predict.all = TRUE, num.threads = ncores)$predictions
    if (is.null(dim(nodes))) nodes <- matrix(nodes, nrow = 1)
    ntree <- obj$num.trees
    nobs <- object$model.rf$num.samples
    nnew <- nrow(x)

    weights <- findweights(origNodes, nodes, inbag, as.integer(nobs), as.integer(nnew), as.integer(ntree)) # cpp function call
    weights.std <- weights / ntree

    priorDensity <- density(resp)

    if (add) {
        rangex <- range(priorDensity$x)
        rangey <- range(priorDensity$y)

        for (i in 1:nnew) {
            postDensity <- density(resp, weights = weights.std[, i], ...)
            rangex <- range(rangex, postDensity$x)
            rangey <- range(rangey, postDensity$y)
        }

        # plot(priorDensity$x, priorDensity$y,
        #     type = "l", main = main, log = log,
        #     xlim = if (is.null(xlim)) rangex else xlim,
        #     ylim = if (is.null(ylim)) rangey else ylim,
        #     xlab = xlab, ylab = ylab, col = "grey"
        # )
        # for (i in 1:nnew) {
        #     postDensity <- density(resp, weights = weights.std[, i], ...)
        #     points(postDensity$x, postDensity$y, type = "l")
        # }
    } else {
        for (i in 1:nnew) {
            postDensity <- density(resp, weights = weights.std[, i], ...)

            # plot(postDensity$x, postDensity$y,
            #     type = "l", main = main, log = log,
            #     xlim = if (is.null(xlim)) range(postDensity$x, priorDensity$x) else xlim,
            #     ylim = if (is.null(ylim)) range(postDensity$y, priorDensity$y) else ylim,
            #     xlab = xlab, ylab = ylab
            # )
            # points(priorDensity$x, priorDensity$y, type = "l", col = "grey")
            # if (nnew > 1 && i < nnew) readline("Press <ENTER> to Continue")
        }
    }



    if (protocol == "TSG") {
        resp <- 1 / resp
    }

    if (cutbound == TRUE) {
        dist_prior <- density(resp, weights = rep(1 / length(resp), length(resp)), from = lower, to = upper)
        dist_posterior <- density(resp, weights = weights.std[, i], from = lower, to = upper)
    } else if (cutbound == FALSE) {
        dist_prior <- density(resp, weights = rep(1 / length(resp), length(resp)))
        dist_posterior <- density(resp, weights = weights.std[, i])
    }
    df_plot_prior <- data.frame(x = dist_prior$x, y = dist_prior$y)
    # df_plot_posterior <- data.frame(x = dist_posterior$x, y = dist_posterior$y)

    df_plot <- data.frame(x = dist_prior$x, y_prior = dist_prior$y, y_posterior = dist_posterior$y)

    # df_plot <- data.frame(dist_raw = resp)
    # df_plot$weight_prior <- 1 / nrow(df_plot)
    # df_plot$weight_posterior <- weights.std[, i]

    return(df_plot)
}

plot_ABC_inference <- function(object,
                               obs,
                               training,
                               add = TRUE,
                               main = "Posterior density",
                               protocol,
                               color_prior = "lightblue",
                               color_posterior = "darkblue",
                               highlight_values = NULL,
                               highlight_colors = NULL,
                               highlight_linetype = NULL,
                               log = "",
                               xlim = NULL,
                               ylim = NULL,
                               xlab = NULL,
                               ylab = NULL,
                               paral = FALSE,
                               fontsize = 50,
                               plot_ABC_prior_as_uniform = FALSE,
                               cutbound = FALSE,
                               para_lower_bound = NULL,
                               para_upper_bound = NULL,
                               ncores = if (paral) max(detectCores() - 1, 1) else 1, ...) {
    df_plot <- densityPlot_df(
        object = object,
        obs = obs,
        training = training,
        add = add,
        main = main,
        color_prior = color_prior,
        chosen_para = chosen_para,
        color_posterior = color_posterior,
        protocol = protocol,
        color_vline = color_vline,
        log = log,
        xlim = xlim,
        ylim = ylim,
        xlab = xlab,
        ylab = ylab,
        paral = paral,
        ncores = ncores,
        cutbound = cutbound,
        lower = para_lower_bound,
        upper = para_upper_bound
    )

    if (plot_ABC_prior_as_uniform) {
        p_plot <- ggplot(df_plot) +
            annotate("rect",
                xmin = para_lower_bound, xmax = para_upper_bound,
                ymin = 0, ymax = 1 / (para_upper_bound - para_lower_bound),
                color = color_prior, fill = color_prior, alpha = 0.3
            )
    } else {
        p_plot <- ggplot(df_plot) +
            geom_area(aes(x = x, y = y_prior),
                color = color_prior, fill = color_prior, alpha = 0.3
            )
    }
    p_plot <- p_plot +
        geom_area(aes(x = x, y = y_posterior), color = color_posterior, fill = color_posterior, alpha = 0.3) +
        # geom_density(aes(x = dist_raw, kernel = "gaussian", weight = weight_prior), color = color_prior, fill = color_prior, alpha = 0.3) +
        # geom_density(aes(x = dist_raw, kernel = "gaussian", weight = weight_posterior), color = color_posterior, fill = color_posterior, alpha = 0.3) +
        xlab("") +
        ylab("") +
        ggtitle(main) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = fontsize)) +
        theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0))
    if (!is.null(highlight_values)) {
        p_plot <- p_plot +
            geom_vline(
                xintercept = highlight_values,
                colour = highlight_colors,
                linetype = highlight_linetype,
                size = 2
            )
    }
    return(p_plot)
}

#' @export
plot_statistics_correlation <- function(filename) {
    #---------------------------------------------------Load the reference table
    library("ggcorrplot")
    load(file = filename)
    #-------------------------------------------------------------Get parameters
    parameters <- ABC_input$sim_param
    # Parameter names appeared in the plot
    param_names <- c(
        "log10(prob_misseg)",
        "Sel.rate(Chrom 1)",
        "Sel.rate(Chrom 2)",
        "Sel.rate(Chrom 3)",
        "Sel.rate(Chrom 4)",
        "Sel.rate(Chrom 5)",
        "Sel.rate(Chrom 6)",
        "Sel.rate(Chrom 7)",
        "Sel.rate(Chrom 8)",
        "Sel.rate(Chrom 9)",
        "Sel.rate(Chrom 10)",
        "Sel.rate(Chrom 11)",
        "Sel.rate(Chrom 12)",
        "Sel.rate(Chrom 13)",
        "Sel.rate(Chrom 14)",
        "Sel.rate(Chrom 15)",
        "Sel.rate(Chrom 16)",
        "Sel.rate(Chrom 17)",
        "Sel.rate(Chrom 18)",
        "Sel.rate(Chrom 19)",
        "Sel.rate(Chrom 20)",
        "Sel.rate(Chrom 21)",
        "Sel.rate(Chrom 22)"
    )
    colnames(parameters) <- param_names
    #---------------------------------------------------Get statistics
    # list of statistics used by misseg rate
    misseg_stat_names <- c(
        #-----bulk
        "Wasserstein_dist_bulk_genome",
        "Wasserstein_dist_sc_genome",
        "Mean_Shannon_genome",
        "Mean_Clonal_misseg_count_bulk_genome",
        "Mean_Clonal_misseg_count_sc_genome",
        "Mean_Subclonal_misseg_count_sc_genome",
        "Mean_Clonal_armmisseg_count_sc_genome",
        "Mean_Subclonal_armmisseg_count_sc_genome",
        #-----sc_subclonal CN
        "Var_Shannon_genome",
        "Var_Clonal_misseg_count_bulk_genome",
        "Var_Clonal_misseg_count_sc_genome",
        "Var_Subclonal_misseg_count_sc_genome",
        "Var_Clonal_armmisseg_count_sc_genome",
        "Var_Subclonal_armmisseg_count_sc_genome",
        #-----sc_Phylo Stats
        "Mean_Cherries_genome",
        "Mean_Pitchforks_genome",
        "Mean_IL_number_genome",
        "Mean_AvgLadder_genome",
        "Var_Cherries_genome",
        "Var_Pitchforks_genome",
        "Var_IL_number_genome",
        "Var_AvgLadder_genome",
        # balance
        "Mean_Stairs_genome",
        "Mean_Colless_genome",
        "Mean_Sackin_genome",
        "Mean_B2_genome",
        "Mean_MaxDepth_genome",
        "Var_Stairs_genome",
        "Var_Colless_genome",
        "Var_Sackin_genome",
        "Var_B2_genome",
        "Var_MaxDepth_genome"
    )

    # List of statistics used by selection rates
    sel_stat_names <- c(
        #-----bulk
        "Wasserstein_dist_bulk_chr",
        "Wasserstein_dist_sc_chr",
        "Mean_Shannon_chr",
        "Mean_Clonal_misseg_count_bulk_chr",
        #-----sc_subclonal CN
        "Mean_Clonal_misseg_count_sc_chr",
        "Mean_Subclonal_misseg_count_sc_chr",
        "Mean_Clonal_armmisseg_count_sc_chr",
        "Mean_Subclonal_armmisseg_count_sc_chr",
        "Var_Shannon_chr",
        "Var_Clonal_misseg_count_bulk_chr",
        "Var_Clonal_misseg_count_sc_chr",
        "Var_Subclonal_misseg_count_sc_chr",
        "Var_Clonal_armmisseg_count_sc_chr",
        "Var_Subclonal_armmisseg_count_sc_chr",
        #-----sc_Phylo Stats
        "Mean_Cherries_genome",
        "Mean_Pitchforks_genome",
        "Mean_IL_number_genome",
        "Mean_AvgLadder_genome",
        "Var_Cherries_genome",
        "Var_Pitchforks_genome",
        "Var_IL_number_genome",
        "Var_AvgLadder_genome",
        # balance
        "Mean_Stairs_genome",
        "Mean_Colless_genome",
        "Mean_Sackin_genome",
        "Mean_B2_genome",
        "Mean_MaxDepth_genome",
        "Var_Stairs_genome",
        "Var_Colless_genome",
        "Var_Sackin_genome",
        "Var_B2_genome",
        "Var_MaxDepth_genome"
    )

    # List of statistics in the reference table
    list_target_stats <- c(
        #-----bulk
        "Wasserstein_dist_bulk_genome",
        "Mean_Clonal_misseg_count_bulk_genome",
        "Var_Clonal_misseg_count_bulk_genome",
        "Wasserstein_dist_bulk_chr",
        "Mean_Clonal_misseg_count_bulk_chr",
        "Var_Clonal_misseg_count_bulk_chr",
        #-----sc_subclonal CN
        "Wasserstein_dist_sc_genome",
        "Mean_Shannon_genome",
        "Mean_Clonal_misseg_count_sc_genome",
        "Mean_Subclonal_misseg_count_sc_genome",
        "Mean_Clonal_armmisseg_count_sc_genome",
        "Mean_Subclonal_armmisseg_count_sc_genome",
        "Var_Shannon_genome",
        "Var_Clonal_misseg_count_sc_genome",
        "Var_Subclonal_misseg_count_sc_genome",
        "Var_Clonal_armmisseg_count_sc_genome",
        "Var_Subclonal_armmisseg_count_sc_genome",
        "Wasserstein_dist_sc_chr",
        "Mean_Shannon_chr",
        "Mean_Clonal_misseg_count_sc_chr",
        "Mean_Subclonal_misseg_count_sc_chr",
        "Mean_Clonal_armmisseg_count_sc_chr",
        "Mean_Subclonal_armmisseg_count_sc_chr",
        "Var_Shannon_chr",
        "Var_Clonal_misseg_count_sc_chr",
        "Var_Subclonal_misseg_count_sc_chr",
        "Var_Clonal_armmisseg_count_sc_chr",
        "Var_Subclonal_armmisseg_count_sc_chr",
        #-----sc_Phylo Stats
        "Mean_Cherries_genome",
        "Mean_Pitchforks_genome",
        "Mean_IL_number_genome",
        "Mean_AvgLadder_genome",
        "Var_Cherries_genome",
        "Var_Pitchforks_genome",
        "Var_IL_number_genome",
        "Var_AvgLadder_genome",
        # balance
        "Mean_Stairs_genome",
        "Mean_Colless_genome",
        "Mean_Sackin_genome",
        "Mean_B2_genome",
        "Mean_MaxDepth_genome",
        "Var_Stairs_genome",
        "Var_Colless_genome",
        "Var_Sackin_genome",
        "Var_B2_genome",
        "Var_MaxDepth_genome"
    )

    # Name of statistics appeared on the plot
    stat_names <- c(
        "Bulk CN distance",
        "Single-cell CN distance",
        "Mean(sc Shannon index)",
        "Mean(misseg. count in bulk)",
        "Mean(clonal misseg. count in sc)",
        "Mean(subclonal misseg. count in sc)",
        "NA",
        "NA",
        "Var(sc Shannon index)",
        "Var(misseg. count in bulk)",
        "Var(clonal misseg. count in sc)",
        "Var(subclonal misseg. count in sc)",
        "NA",
        "NA",
        "Mean(cherry count)",
        "Mean(pitchfork count)",
        "Mean(IL number)",
        "Mean(average ladder)",
        "Var(cherry count)",
        "Var(pitchfork count)",
        "Var(IL number)",
        "Var(average ladder)",
        "Mean(stairs)",
        "Mean(Colless index)",
        "Mean(Sackin index)",
        "Mean(B2 index)",
        "Mean(max depth)",
        "Var(stairs)",
        "Var(Colless index)",
        "Var(Sackin index)",
        "Var(B2 index)",
        "Var(max depth)"
    )

    #-----------------=--------------------------Get correlation matrix (misseg)
    statistics_misseg <- data.frame(matrix(ncol = 0, nrow = 100000))
    for (i in 1:length(misseg_stat_names)) {
        statistics_misseg[, i] <- ABC_input$sim_stat[[which(grepl(misseg_stat_names[i], list_target_stats))]]
    }
    prob_misseg <- data.frame(ABC_input$sim_param[, 1])
    corr_mtx_misseg <- cor(y = prob_misseg, x = statistics_misseg)
    #-----------------------------------------Get correlation matrix (selection)
    corr_mtx_sel <- NULL
    for (i in 2:length(param_names)) {
        statistics_sel <- data.frame(matrix(ncol = 0, nrow = 100000))
        for (j in 1:length(sel_stat_names)) {
            if (grepl("genome", sel_stat_names[j])) {
                statistics_sel[, j] <- ABC_input$sim_stat[[which(grepl(sel_stat_names[j], list_target_stats))]]
            } else if (grepl("chr", sel_stat_names[j])) {
                statistics_sel[, j] <- ABC_input$sim_stat[[which(grepl(sel_stat_names[j], list_target_stats))]][, (i - 1)]
            }
        }
        corr_mtx_sel <- cbind(corr_mtx_sel, cor(y = parameters[, i], x = statistics_sel))
    }
    corr_mtx <- cbind(corr_mtx_misseg, corr_mtx_sel)
    rownames(corr_mtx) <- stat_names
    colnames(corr_mtx) <- param_names
    #------------------------------------------------------Plot correlation plot
    plot <- ggcorrplot(
        corr_mtx,
        method = "circle",
        ggtheme = ggplot2::theme_bw,
        legend.title = "Correlation"
    ) +
        theme(axis.text.x = element_text(colour = c(
            "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38", "#00BA38",
            "#00BA38", "#00BA38", "#00BA38",
            "#619CFF", "#619CFF", "#619CFF", "#619CFF",
            "#619CFF", "#619CFF", "#619CFF", "#619CFF",
            "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D",
            "#F8766D", "#F8766D", "#F8766D", "#F8766D", "#F8766D"
        ), angle = 45), plot.margin = margin(t = 0.2, r = 0.2, b = 0.2, l = 0.2, unit = "in"))
    ggsave(file = "correlation_plot.png", width = 12, height = 10, units = "in", plot = plot, dpi = 300, limitsize = TRUE)
    return(plot)
}

#' @export
plot_ABC_correlation <- function(inference_result = parameters_inferred,
                                 plot_name = "ABC_correlation_plot.jpg",
                                 value_x = "Ground_truth",
                                 value_y = NULL,
                                 title_x = "",
                                 title_y = "",
                                 error_x = NULL,
                                 error_y = NULL,
                                 title_plot = "",
                                 color_data = "red",
                                 fontsize = 50,
                                 plot_diagonal = FALSE,
                                 plot_Error = FALSE) {
    library(ggplot2)
    library(ehaGoF)
    if (!is.null(error_y)) {
        parameters_inferred$value_y_min <- parameters_inferred[[value_y]] - parameters_inferred[[error_y]]
        parameters_inferred$value_y_max <- parameters_inferred[[value_y]] + parameters_inferred[[error_y]]
        y_min <- min(parameters_inferred$value_y_min)
        y_max <- max(parameters_inferred$value_y_max)
    } else {
        y_min <- min(parameters_inferred[[value_y]])
        y_max <- max(parameters_inferred[[value_y]])
    }
    if (!is.null(error_x)) {
        parameters_inferred$value_x_min <- parameters_inferred[[value_x]] - parameters_inferred[[error_x]]
        parameters_inferred$value_x_max <- parameters_inferred[[value_x]] + parameters_inferred[[error_x]]
        x_min <- min(parameters_inferred$value_x_min)
        x_max <- max(parameters_inferred$value_x_max)
    } else {
        x_min <- min(parameters_inferred[[value_x]])
        x_max <- max(parameters_inferred[[value_x]])
    }
    corr_plot <- ggplot(parameters_inferred, mapping = aes_string(x = value_x, y = value_y)) +
        geom_point(colour = color_data, size = 10) +
        xlim(min(x_min, y_min), max(x_max, y_max)) +
        ylim(min(x_min, y_min), max(x_max, y_max)) +
        xlab(title_x) +
        ylab(title_y) +
        ggtitle(title_plot) +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = fontsize)) +
        theme(aspect.ratio = 1)
    if (!is.null(error_x)) {
        corr_plot <- corr_plot +
            geom_errorbar(aes(ymin = value_x_min, ymax = value_x_max), width = 0, size = 1, colour = color_data)
    }
    if (!is.null(error_y)) {
        corr_plot <- corr_plot +
            geom_errorbar(aes(ymin = value_y_min, ymax = value_y_max), width = 0, size = 1, colour = color_data)
    }
    if (plot_diagonal) {
        corr_plot <- corr_plot + geom_abline(intercept = 0)
    }
    if (plot_Error) {
        Error <- round(compute_error(
            results = parameters_inferred,
            ID_actual = value_x, ID_predicted = value_y
        ), digits = 3)
        text_x <- min(x_min, y_min) + 0.05 * (max(x_max, y_max) - min(x_min, y_min))
        text_y <- min(x_min, y_min) + 0.95 * (max(x_max, y_max) - min(x_min, y_min))
        corr_plot <- corr_plot +
            annotate(
                "text",
                x = text_x, y = text_y,
                label = paste("RMSE=", Error),
                size = unit(fontsize / 2, "pt"), color = color_data, hjust = 0
            )
    }
    jpeg(plot_name, width = 1500, height = 1500)
    print(corr_plot)
    dev.off()
}

#' @export
plot_stats_corr <- function(true_stat, combined_df) {
    library(ggplot2)
    library(grid)
    library(gridExtra)
    layout <- matrix(NA, nrow = 3, ncol = 5)
    gs <- list()
    id <- 0
    id_tmp <- 0
    true_stat <- true_stat[true_stat$simulation %in% combined_df$Simulation, ]
    variables <- c(
        "data.sc.target.genome.variable.event_count.type.clonal.event.missegregation",
        "data.sc.target.genome.variable.event_count.type.subclonal.event.missegregation",
        "data.sc.target.genome.variable.cherries",
        "data.sc.target.genome.variable.pitchforks",
        "data.sc.target.genome.variable.IL_number",
        "data.sc.target.genome.variable.avgLadder",
        "data.sc.target.genome.variable.stairs",
        "data.sc.target.genome.variable.colless",
        "data.sc.target.genome.variable.sackin",
        "data.sc.target.genome.variable.B2",
        "data.sc.target.genome.variable.maxDepth"
    )
    for (var in variables) {
        id <- id + 1
        if (var == "data.sc.target.genome.variable.cherries") {
            id <- 6
        } else if (var == "data.sc.target.genome.variable.stairs") {
            id <- 11
        }
        print(id)
        col <- id %% 5
        if (col == 0) {
            col <- 5
        }
        print(col)
        row <- ceiling(id / 5)
        print(row)
        id_tmp <- id_tmp + 1
        layout[row, col] <- id_tmp
        x <- true_stat[, var]
        y <- combined_df[, var]
        log10_misseg_rate <- combined_df[, "log10_misseg_rate"]
        if (length(x) != length(y)) {
            stop("x and y must be of the same length")
        }
        data_for_plot <- data.frame(x = x, y = y, misseg = log10_misseg_rate)
        if (var == "data.sc.target.genome.variable.cherries") {
            title <- "cherry count"
            title_color <- "#619CFF"
        } else if (var == "data.sc.target.genome.variable.pitchforks") {
            title <- "pitchfork count"
            title_color <- "#619CFF"
        } else if (var == "data.sc.target.genome.variable.IL_number") {
            title <- "IL number"
            title_color <- "#619CFF"
        } else if (var == "data.sc.target.genome.variable.avgLadder") {
            title <- "average ladder"
            title_color <- "#619CFF"
        } else if (var == "data.sc.target.genome.variable.stairs") {
            title <- "stairs"
            title_color <- "#F8766D"
        } else if (var == "data.sc.target.genome.variable.colless") {
            title <- "Colless index"
            title_color <- "#F8766D"
        } else if (var == "data.sc.target.genome.variable.sackin") {
            title <- "Sackin index"
            title_color <- "#F8766D"
        } else if (var == "data.sc.target.genome.variable.B2") {
            title <- "B2 index"
            title_color <- "#F8766D"
        } else if (var == "data.sc.target.genome.variable.maxDepth") {
            title <- "max depth"
            title_color <- "#F8766D"
        } else if (var == "data.sc.target.genome.variable.event_count.type.clonal.event.missegregation") {
            title <- "clonal misseg. count"
            title_color <- "#00BA38"
        } else if (var == "data.sc.target.genome.variable.event_count.type.subclonal.event.missegregation") {
            title <- "subclonal misseg. count"
            title_color <- "#00BA38"
        }
        if (var == "data.sc.target.genome.variable.event_count.type.subclonal.event.missegregation") {
            labels <- waiver()
        } else {
            labels <- NULL
        }
        gs[[id_tmp]] <- ggplot(data_for_plot, aes(x = x, y = y, color = misseg)) +
            geom_point(size = 10, alpha = 0.3) +
            scale_color_gradient(limits = c(-5, -3), low = "#dc1dbf", high = "#00fff7", labels = labels) +
            labs(color = NULL) +
            theme(
                legend.key.size = unit(1, "cm"),
                legend.key.height = unit(2, "cm"), # change legend key height
                legend.key.width = unit(1, "cm") # change legend key width
            ) +
            xlab("Ground truth") +
            ylab("Medicc2") +
            xlim(min(min(x), min(y)), max(max(x), max(y))) +
            ylim(min(min(x), min(y)), max(max(x), max(y))) +
            ggtitle(title) +
            theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
            theme(text = element_text(size = 30)) +
            theme(plot.title = element_text(color = title_color)) +
            theme(aspect.ratio = 1) +
            geom_abline(intercept = 0, slope = 1)
        filename <- "MEDICC_correlation_plot.jpg"
        jpeg(filename, width = 2500, height = 1500)
        p <- grid.arrange(grobs = gs, layout_matrix = layout)
        print(p)
        dev.off()
    }
}

#' @export
plot_RF <- function(df) {
    library(ggplot2)
    library(grid)
    library(gridExtra)
    plot <- ggplot(df, aes(x = log10_misseg_rate, y = TreeDistance)) +
        geom_point(colour = "#a8a0a0", size = 10, alpha = 0.5) +
        xlab("log10(prob_misseg)") +
        ylab("Generalized RF distance") +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = 50))
    jpeg("MEDICC_RF_plot.jpg", width = 2000, height = 1000)
    print(plot)
    dev.off()
}
