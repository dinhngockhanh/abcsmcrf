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
    # =====================================Plot True Marginal Posteriors
    p_theta1 <- ggplot() +
        geom_density(data = data.frame(parameters_truth$theta1), aes(x = parameters_truth$theta1, fill = "True Posterior", color = "True Posterior"), size = 2) +
        theme(
            text = element_text(size = 20),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white")
        )
    p_theta2 <- ggplot() +
        geom_density(data = data.frame(parameters_truth$theta2), aes(x = parameters_truth$theta2, fill = "True Posterior", color = "True Posterior"), size = 2) +
        theme(
            text = element_text(size = 20),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white")
        )
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
    for (iteration in 1:nIter) {
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
        ##################################################  TEST - BEGIN
        if (iteration == 1) {
            p_theta1 <- p_theta1 + geom_density(data = parameters, aes(x = theta1, fill = "Prior Distribution", color = "Prior Distribution"), alpha = 0.2, size = 2)
            p_theta2 <- p_theta2 + geom_density(data = parameters, aes(x = theta2, fill = "Prior Distribution", color = "Prior Distribution"), alpha = 0.2, size = 2)
        }
        ####################################################  TEST - END
        #   Simulate statistics
        reference <- model(
            parameters = parameters,
            n_samples_per_parameter_set = n_samples_per_parameter_set,
            nNoise = nNoise
        )
        filename <- paste0("ABCSMC_RF_Input_Iter_", iteration, ".rda")
        ABCSMC_RF_input <- list()
        ABCSMC_RF_input$params <- as.list(reference[c(1:ncol(parameters_initial))])
        ABCSMC_RF_input$statistics <- as.list(reference[-c(1:ncol(parameters_initial))])
        save(ABCSMC_RF_input, file = filename)
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
            ############################################
            df_dist <- densityPlot_df(
                object = RFmodel, obs = target, training = mini_reference
            )

            dist_output <- list()
            dist_output$x <- df_dist$x
            dist_output$y_prior <- df_dist$y_prior
            dist_output$y_posterior <- df_dist$y_posterior
            filename <- paste0("dist_output_iter_", iteration, "_", parameter_id, ".rda")
            save(dist_output, file = filename)
            post_mean <- weighted.mean(df_dist$x, df_dist$y_posterior)
            post_sd <- weightedSd(df_dist$x, df_dist$y_posterior)
            post_median <- weightedMedian(df_dist$x, df_dist$y_posterior)
            post_mode <- df_dist$x[which(df_dist$y_posterior == max(df_dist$y_posterior))]
            ############################################
            test_mean <- weighted.mean(mini_reference[, parameter_id], posterior_gamma_RF$weights)
            test_sd <- weightedSd(mini_reference[, parameter_id], posterior_gamma_RF$weights)
            test_median <- weightedMedian(mini_reference[, parameter_id], posterior_gamma_RF$weights)
            test_mode <- mini_reference[, parameter_id][which(posterior_gamma_RF$weights == max(posterior_gamma_RF$weights))]
            # ===========================================SAVE OUTPUT CSV
            list_parameters_output[nrow(list_parameters_output) + 1, ] <- c(
                paste0("Iteration_", iteration), parameter_id, post_mean, post_median,
                post_mode, post_sd
            )
            ####################################################################
            cat("Posterior mean: ", test_mean, ";", post_mean, "\n")
            cat("Posterior median: ", test_median, ";", post_median, "\n")
            cat("Posterior mode: ", test_mode, ";", post_mode, "\n")
            cat("Posterior sd: ", test_sd, ";", post_sd, "\n")
            ABCRF_weights[, parameter_id] <- posterior_gamma_RF$weights
        }
        ##################################################  TEST - BEGIN
        print(paste0("Iter-", iteration))
        w <- ABCRF_weights
        colnames(w) <- paste0("weights_", colnames(w))
        w$iteration <- iteration
        p_theta1 <- p_theta1 +
            geom_density(
                data = cbind(parameters, w),
                aes(x = theta1, weight = weights_theta1, fill = paste0("Iter-", iteration), color = paste0("Iter-", iteration)), alpha = 0.2, size = 2
            ) +
            scale_fill_manual(values = color_scheme, name = "Iteration") +
            scale_color_manual(values = color_scheme, name = "Iteration")
        p_theta2 <- p_theta2 +
            geom_density(
                data = cbind(parameters, w),
                aes(x = theta2, weight = weights_theta2, fill = paste0("Iter-", iteration), color = paste0("Iter-", iteration)), alpha = 0.2, size = 2
            ) +
            scale_fill_manual(values = color_scheme, name = "Iteration") +
            scale_color_manual(values = color_scheme, name = "Iteration")
        ####################################################  TEST - END
    }
    write.csv(list_parameters_output, paste0("ABCSMC_RF_para_output.csv"))
    ######################################################  TEST - BEGIN
    png("TEST.png", res = 150, width = 17, height = 17, units = "in", pointsize = 12)
    grid.arrange(p_theta1, p_theta2, ncol = 1)
    dev.off()
    ########################################################  TEST - END
}



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
