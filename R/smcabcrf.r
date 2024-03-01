smcabcrf <- function(target,
                     model,
                     n_samples_per_parameter_set,
                     nNoise,
                     perturb,
                     #  range,
                     parameters_initial,
                     nIter,
                     nParticles,
                     parallel,
                     ...) {
    library(abcrf)
    if (length(nParticles) < nIter) nParticles[(length(nParticles) + 1):nIter] <- nParticles[length(nParticles)]
    parameters_id <- colnames(parameters_initial)
    SMCRF <- list()
    for (iteration in 1:nIter) {
        #   Sample parameters for this round of iteration
        if (iteration == 1) {
            parameters <- parameters_initial[1:nParticles[iteration], ]
        } else {
            parameters_new <- data.frame(matrix(NA, nrow = nParticles[iteration], ncol = 0))
            for (parameter_id in parameters_id) {
                parameters_new[, parameter_id] <-
                    sample(parameters[, parameter_id], size = nParticles[iteration], prob = ABCRF_weights[, parameter_id], replace = T)
            }
            parameters <- perturb(parameters_new)
        }
        #   Simulate statistics
        reference <- model(parameters = parameters, n_samples_per_parameter_set = n_samples_per_parameter_set, nNoise = nNoise)
        #   Run ABCRF for each parameter
        ABCRF_weights <- data.frame(matrix(NA, nrow = nParticles[iteration], ncol = 0))
        for (parameter_id in parameters_id) {
            mini_reference <- reference[, c(parameter_id, colnames(reference)[!colnames(reference) %in% parameters_id])]
            RFmodel <- regAbcrf(
                formula = as.formula(paste0(parameter_id, " ~ .")),
                data = mini_reference,
                paral = parallel
            )
            posterior_gamma_RF <- predict(
                object = RFmodel,
                obs = target,
                training = mini_reference,
                paral = parallel,
                rf.weights = T
            )
            ABCRF_weights[, parameter_id] <- posterior_gamma_RF$weights
        }
        #   Save SMC-RF results from this iteration
        SMCRF_iteration <- list()
        SMCRF_iteration$parameters <- parameters
        SMCRF_iteration$statistics <- reference[, -c(1:length(parameters_id))]
        SMCRF_iteration$weights <- ABCRF_weights
        SMCRF_iteration$rf_model <- RFmodel
        SMCRF_iteration$rf_predict <- posterior_gamma_RF
        SMCRF[[iteration]] <- SMCRF_iteration
    }
    SMCRF[["method"]] <- "smcrf-single-param"
    SMCRF[["nIter"]] <- nIter
    SMCRF[["nParticles"]] <- nParticles
    SMCRF[["target"]] <- target
    return(SMCRF)
}

plotting_smcrf <- function(
    parameters_truth = NULL,
    parameters_initial = NULL,
    parameters_id,
    outputdata,
    Plot_stats = FALSE) {
    nIter <- outputdata[["nIter"]]
    fixed_colors <- c("True Posterior" = "black", "Prior Distribution" = "gray")
    iter_colors <- rev(rainbow(nIter))
    iter_names <- paste0("Iteration ", 1:nIter)
    color_scheme <- c(fixed_colors, setNames(iter_colors, iter_names))
    legend_order <- c("True Posterior", "Prior Distribution", iter_names)
    for (parameter_id in parameters_id) {
        p_para <- ggplot() +
            theme(
                text = element_text(size = 30),
                panel.background = element_rect(fill = "white", colour = "white"),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white")
            )
        if (!is.null(parameters_truth)) {
            if (parameter_id == "theta") {
                vec_t <- seq(1, 20, 0.1)
                const <- integrate(function(t) t^target$K * gamma(t) / gamma(n_samples_per_parameter_set + t), lower = 1, upper = 20)
                vec_pdf <- vec_t^target$K * gamma(vec_t) / (const$value * gamma(100 + vec_t))
                p_para <- p_para +
                    geom_area(aes(x = vec_t, y = vec_pdf, fill = "True Posterior", color = "True Posterior"), size = 4)
            } else {
                true_posterior_df <- data.frame(value = parameters_truth[[parameter_id]], legend = "True Posterior")
                p_para <- p_para +
                    geom_density(data = true_posterior_df, aes(x = value, fill = legend, color = legend), size = 2)
            }
        }
        if (!is.null(parameters_initial)) {
            prior_df <- data.frame(value = parameters_initial[[parameter_id]], legend = "Prior Distribution")
            p_para <- p_para +
                geom_density(data = prior_df, aes(x = value, fill = legend, color = legend), alpha = 0.2, size = 2)
        }

        for (iteration in 1:nIter) {
            if (Plot_stats) {
                if (iteration == 1) {
                    iter_df <- data.frame(value = outputdata[[iteration]]$statistics[[parameter_id]], legend = "Prior Distribution")
                } else {
                    iter_df <- data.frame(value = outputdata[[iteration]]$statistics[[parameter_id]], legend = paste0("Iteration ", iteration - 1))
                }
                p_para <- p_para +
                    geom_density(
                        data = iter_df,
                        aes(x = value, fill = legend, color = legend), alpha = 0.2, size = 2
                    ) + xlab(parameter_id)
            } else {
                para_tmp <- outputdata[[iteration]]$parameters[[parameter_id]]
                if (outputdata[["method"]] == "smcrf-single-param") {
                    weights_tmp <- outputdata[[iteration]]$weights[[parameter_id]]
                } else if (outputdata[["method"]] == "smcrf-multi-param") {
                    weights_tmp <- outputdata[[iteration]]$weights[, 1]
                }
                iter_df <- data.frame(value = para_tmp, weight = weights_tmp, legend = paste0("Iteration ", iteration))
                p_para <- p_para +
                    geom_density(
                        data = iter_df,
                        aes(x = value, weight = weight, fill = legend, color = legend), alpha = 0.2, size = 2
                    ) + xlab(parameter_id)
            }
        }
        if (Plot_stats) {
            p_para <- p_para + geom_vline(
                xintercept = outputdata$target[[parameter_id]],
                colour = "black",
                linetype = "solid",
                size = 2
            )
        }
        p_para <- p_para +
            scale_fill_manual(values = color_scheme, name = "", breaks = legend_order) +
            scale_color_manual(values = color_scheme, name = "", breaks = legend_order) +
            theme(legend.position = c(0, 1), legend.justification = c(0, 0.5)) +
            guides(fill = guide_legend(nrow = 1, keywidth = 1, keyheight = 1))

        file_name <- paste0(outputdata[["method"]], "_iteration_plots_", parameter_id, ".png")
        png(file_name, res = 150, width = 30, height = 15, units = "in", pointsize = 12)
        print(p_para)
        dev.off()
    }
}
