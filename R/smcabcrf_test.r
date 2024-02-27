smcabcrf_fitting <- function(target,
                             model,
                             N,
                             n_samples_per_parameter_set, # Number of samples per parameter set
                             nNoise,
                             perturb,
                             parameters_initial,
                             parameters_truth,
                             nIter,
                             nParticles,
                             parallel,
                             ...) {
    library(abcrf)
    library(matrixStats)
    library(data.table)
    ########################################################  TEST - END
    if (length(nParticles) < nIter) nParticles[(length(nParticles) + 1):nIter] <- nParticles[length(nParticles)]
    parameters_id <- colnames(parameters_initial)
    # ===================================================SAVE OUTPUT CSV
    list_parameters_output <- data.frame(matrix(ncol = 6, nrow = 0))
    colnames(list_parameters_output) <- c("Iteration", "Variable", "Mean", "Median", "Mode", "Sd")
    #####################################################
    # Estimate density for theta1 and theta2
    density_theta1 <- density(parameters_truth$theta1)
    density_theta2 <- density(parameters_truth$theta2)
    # Calculate the weighted mean
    weighted_mean_theta1 <- weighted.mean(density_theta1$x, w = density_theta1$y)
    weighted_mean_theta2 <- weighted.mean(density_theta2$x, w = density_theta2$y)
    weighted_median_theta1 <- weightedMedian(density_theta1$x, w = density_theta1$y)
    weighted_median_theta2 <- weightedMedian(density_theta2$x, w = density_theta2$y)
    mode_theta1 <- density_theta1$x[which.max(density_theta1$y)]
    mode_theta2 <- density_theta2$x[which.max(density_theta2$y)]
    weighted_sd_theta1 <- weightedSd(density_theta1$x, density_theta1$y)
    weighted_sd_theta2 <- weightedSd(density_theta2$x, density_theta2$y)
    list_parameters_output[nrow(list_parameters_output) + 1, ] <- c(
        "True_posteriors", "theta1", weighted_mean_theta1, weighted_median_theta1,
        mode_theta1, weighted_sd_theta1
    )
    list_parameters_output[nrow(list_parameters_output) + 1, ] <- c(
        "True_posteriors", "theta2", weighted_mean_theta2, weighted_median_theta2,
        mode_theta2, weighted_sd_theta2
    )
    #####################################################
    ################################################################
    filename <- paste0("ABCSMC_RF_OUTPUT.rda")
    ABCSMC_RF_output <- list()
    ABCSMC_RF_output$RFmodel <- list()
    ABCSMC_RF_output$posterior_gamma_RF <- list()
    ABCSMC_RF_output$priors <- parameters_initial[1:nParticles[1], ]
    ABCSMC_RF_output$truth <- parameters_truth
    for (parameter_id in parameters_id) {
        ABCSMC_RF_output$posterior_gamma_RF[[parameter_id]] <- list()
        ABCSMC_RF_output$RFmodel[[parameter_id]] <- list()
    }
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
        reference <- model(
            parameters = parameters,
            n_samples_per_parameter_set = n_samples_per_parameter_set,
            nNoise = nNoise
        )
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
            post_mean <- weighted.mean(mini_reference[, parameter_id], posterior_gamma_RF$weights)
            post_sd <- weightedSd(mini_reference[, parameter_id], posterior_gamma_RF$weights)
            post_median <- weightedMedian(mini_reference[, parameter_id], posterior_gamma_RF$weights)
            post_mode <- mini_reference[, parameter_id][which(posterior_gamma_RF$weights == max(posterior_gamma_RF$weights))]
            # ===========================================SAVE OUTPUT CSV
            list_parameters_output[nrow(list_parameters_output) + 1, ] <- c(
                paste0("Iteration_", iteration), parameter_id, post_mean, post_median,
                post_mode, post_sd
            )
            ####################################################################
            cat("Posterior mean: ", post_mean, "\n")
            cat("Posterior median: ", post_median, "\n")
            cat("Posterior mode: ", post_mode, "\n")
            cat("Posterior sd: ", post_sd, "\n")
            ABCRF_weights[, parameter_id] <- posterior_gamma_RF$weights
            if (iteration == 1) {
                ABCSMC_RF_output$RFmodel[[parameter_id]] <- list(RFmodel)
                ABCSMC_RF_output$posterior_gamma_RF[[parameter_id]] <- list(posterior_gamma_RF)
            } else {
                ABCSMC_RF_output$RFmodel[[parameter_id]] <- c(ABCSMC_RF_output$RFmodel[[parameter_id]], list(RFmodel))
                ABCSMC_RF_output$posterior_gamma_RF[[parameter_id]] <- c(ABCSMC_RF_output$posterior_gamma_RF[[parameter_id]], list(posterior_gamma_RF))
            }
        }
        if (iteration == 1) {
            ABCSMC_RF_output$params <- reference[, c(1:ncol(parameters_initial))]
            ABCSMC_RF_output$statistics <- reference[, -c(1:ncol(parameters_initial))]
            ABCSMC_RF_output$ABCRF_weights <- ABCRF_weights
        } else {
            ABCSMC_RF_output$params <- rbind(ABCSMC_RF_output$params, reference[, c(1:ncol(parameters_initial))])
            ABCSMC_RF_output$statistics <- rbind(ABCSMC_RF_output$statistics, reference[, -c(1:ncol(parameters_initial))])
            ABCSMC_RF_output$ABCRF_weights <- rbind(ABCSMC_RF_output$ABCRF_weights, ABCRF_weights)
        }
        ##################################################  TEST - BEGIN
        print(paste0("Iter-", iteration))
        w <- ABCRF_weights
        colnames(w) <- paste0("weights_", colnames(w))
        w$iteration <- iteration
    }
    filename <- "ABCSMC_RF_output.rda"
    save(ABCSMC_RF_output, file = filename)
    write.csv(list_parameters_output, paste0("ABCSMC_RF_para_output.csv"))
    return(ABCSMC_RF_output)
}

plotting_smcrf <- function(
    parameters_truth = NULL,
    parameters_initial = NULL,
    parameters_id,
    outputdata,
    nParticles,
    color_scheme) {
    for (parameter_id in parameters_id) {
        color_scheme <- c(
            "True Posterior" = "black",
            "Prior Distribution" = "gray",
            "Iter-1" = "purple",
            "Iter-2" = "blue",
            "Iter-3" = "cyan",
            "Iter-4" = "green",
            "Iter-5" = "yellow",
            "Iter-6" = "orange",
            "Iter-7" = "red"
        )
        p_para <- ggplot() +
            theme(
                text = element_text(size = 20),
                panel.background = element_rect(fill = "white", colour = "white"),
                panel.grid.major = element_line(colour = "white"),
                panel.grid.minor = element_line(colour = "white")
            )
        if (!is.null(parameters_truth)) {
            p_para <- p_para +
                geom_density(data = data.frame(parameters_truth[[parameter_id]]), aes(x = parameters_truth[[parameter_id]]), fill = color_scheme["True Posterior"], color = color_scheme["True Posterior"], size = 2)
        }

        if (!is.null(parameters_initial)) {
            prior_tmp <- data.frame(parameters_initial[[parameter_id]])
            names(prior_tmp) <- parameter_id
            p_para <- p_para +
                geom_density(data = prior_tmp, aes(x = prior_tmp[[parameter_id]]), fill = color_scheme["Prior Distribution"], color = color_scheme["Prior Distribution"], alpha = 0.2, size = 2)
        }
        for (iteration in 1:length(nParticles)) {
            para_tmp <- outputdata$params[(sum(nParticles[0:(iteration - 1)]) + 1):sum(nParticles[0:(iteration)]), ][[parameter_id]]
            weights_tmp <- outputdata$ABCRF_weights[(sum(nParticles[0:(iteration - 1)]) + 1):sum(nParticles[0:(iteration)]), ][[parameter_id]]
            data_tmp <- data.frame(para_tmp, weight = weights_tmp, iter = paste0("Iter-", iteration))
            p_para <- p_para +
                geom_density(
                    data = data_tmp,
                    aes(x = para_tmp, weight = weights_tmp, fill = iter, color = iter), alpha = 0.2, size = 2
                )
        }
        color_scheme <- c(
            "True Posterior" = "black",
            "Prior Distribution" = "gray",
            "Iter-1" = "purple",
            "Iter-2" = "blue",
            "Iter-3" = "cyan",
            "Iter-4" = "green",
            "Iter-5" = "yellow",
            "Iter-6" = "orange",
            "Iter-7" = "red"
        )
        p_para <- p_para +
            scale_fill_manual(values = color_scheme, name = "Legend") +
            scale_color_manual(values = color_scheme, name = "Legend")


        file_name <- paste0("Iteration_plots_", parameter_id, ".png")
        png(file_name, res = 150, width = 17, height = 17, units = "in", pointsize = 12)
        print(p_para)
        dev.off()
    }
}
