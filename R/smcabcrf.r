smcabcrf <- function(target,
                     model,
                     perturb,
                     parameters_initial,
                     nIter,
                     nParticles,
                     parallel,
                     ...) {
    library(abcrf)
    if (length(nParticles) < nIter) nParticles[(length(nParticles) + 1):nIter] <- nParticles[length(nParticles)]
    parameters_id <- colnames(parameters_initial)
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
        reference <- model(parameters)
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
    }
}

smcabcrf_test <- function(target,
                          model,
                          perturb,
                          parameters_initial,
                          nIter,
                          nParticles,
                          parallel,
                          ...) {
    library(abcrf)
    ######################################################  TEST - BEGIN
    N <- 100000
    n <- 20
    #   True marginal posteriors for theta_1
    s_2 <- target[, "variance"] * (n - 1)
    ybar <- target[, "expectation"]
    mean_theta1 <- n / (n + 1) * ybar
    scale_theta1 <- sqrt((2 * (3 + s_2 / 2 + n * ybar^2 / (2 * n + 2))) / ((n + 1) * (n + 8)))
    t_deviate <- rt(N, df = n + 8)
    theta1_true <- mean_theta1 + t_deviate * scale_theta1
    p_theta1 <- ggplot() +
        geom_density(data = data.frame(theta1_true = theta1_true), aes(x = theta1_true, fill = "True Posterior", color = "True Posterior"), size = 2) +
        theme(
            text = element_text(size = 20),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white")
        )
    #   True marginal posterior for theta_2
    shape_theta2 <- n / 2 + 4
    scale_theta2 <- 0.5 * (s_2 + 6 + n * ybar^2 / (n + 1))
    theta2_true <- rinvgamma(N, shape = shape_theta2, scale = 1 / scale_theta2)
    p_theta2 <- ggplot() +
        geom_density(data = data.frame(theta2_true = theta2_true), aes(x = theta2_true, fill = "True Posterior", color = "True Posterior"), size = 2) +
        theme(
            text = element_text(size = 20),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white")
        )
    ########################################################  TEST - END
    if (length(nParticles) < nIter) nParticles[(length(nParticles) + 1):nIter] <- nParticles[length(nParticles)]
    parameters_id <- colnames(parameters_initial)
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
        ######################################################  TEST - BEGIN
        if (iteration == 1) {
            p_theta1 <- p_theta1 + geom_density(data = parameters, aes(x = theta1, fill = "Prior Distribution", color = "Prior Distribution"), alpha = 0.2, size = 2)
            p_theta2 <- p_theta2 + geom_density(data = parameters, aes(x = theta2, fill = "Prior Distribution", color = "Prior Distribution"), alpha = 0.2, size = 2)
        }
        ########################################################  TEST - END
        #   Simulate statistics
        reference <- model(parameters)
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
        ######################################################  TEST - BEGIN
        print(paste0("Iter-", iteration))
        w <- ABCRF_weights
        colnames(w) <- paste0("weights_", colnames(w))
        w$iteration <- iteration
        p_theta1 <- p_theta1 +
            geom_density(
                data = cbind(parameters, w),
                aes(x = theta1, weight = weights_theta1, fill = paste0("Iter-", iteration), color = paste0("Iter-", iteration)), alpha = 0.2, size = 2
            )
        p_theta2 <- p_theta2 +
            geom_density(
                data = cbind(parameters, w),
                aes(x = theta2, weight = weights_theta2, fill = paste0("Iter-", iteration), color = paste0("Iter-", iteration)), alpha = 0.2, size = 2
            )
        ########################################################  TEST - END
    }
    ######################################################  TEST - BEGIN
    p_theta1 <- p_theta1 +
        scale_fill_manual(values = c(
            "True Posterior" = "black",
            "Prior Distribution" = "gray",
            "Iter-1" = "purple",
            "Iter-2" = "blue",
            "Iter-3" = "cyan",
            "Iter-4" = "green",
            "Iter-5" = "yellow",
            "Iter-6" = "orange",
            "Iter-7" = "red"
        ), name = "Legend") +
        scale_colour_manual(values = c(
            "True Posterior" = "black",
            "Prior Distribution" = "gray",
            "Iter-1" = "purple",
            "Iter-2" = "blue",
            "Iter-3" = "cyan",
            "Iter-4" = "green",
            "Iter-5" = "yellow",
            "Iter-6" = "orange",
            "Iter-7" = "red"
        ), name = "Iteration")
    p_theta2 <- p_theta2 +
        scale_fill_manual(values = c(
            "True Posterior" = "black",
            "Prior Distribution" = "gray",
            "Iter-1" = "purple",
            "Iter-2" = "blue",
            "Iter-3" = "cyan",
            "Iter-4" = "green",
            "Iter-5" = "yellow",
            "Iter-6" = "orange",
            "Iter-7" = "red"
        ), name = "Distribution")
    scale_colour_manual(values = c(
        "True Posterior" = "black",
        "Prior Distribution" = "gray",
        "Iter-1" = "purple",
        "Iter-2" = "blue",
        "Iter-3" = "cyan",
        "Iter-4" = "green",
        "Iter-5" = "yellow",
        "Iter-6" = "orange",
        "Iter-7" = "red"
    ), name = "Distribution")
    png("TEST.png", res = 150, width = 17, height = 17, units = "in", pointsize = 12)
    grid.arrange(p_theta1, p_theta2, ncol = 1)
    dev.off()
    ########################################################  TEST - END
}

smcabcrf_test_2 <- function(target,
                            model,
                            perturb,
                            parameters_initial,
                            nIter,
                            nParticles,
                            parallel,
                            parameters_truth,
                            ...) {
    library(abcrf)
    ######################################################  TEST - BEGIN
    vec_t <- seq(1, 20, 0.1)

    n <- 100



    const <- integrate(function(t) t^target$K * gamma(t) / gamma(n + t), lower = 1, upper = 20)
    vec_pdf <- vec_t^target$K * gamma(vec_t) / (const$value * gamma(100 + vec_t))



    p_theta <- ggplot() +
        geom_area(aes(x = vec_t, y = vec_pdf, fill = "Ground truth", color = "Ground truth"), size = 4) +
        geom_vline(aes(xintercept = 1, color = "Bounds"), size = 1) +
        geom_vline(aes(xintercept = 20, color = "Bounds"), size = 1) +
        theme(
            text = element_text(size = 20),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white")
        )
    p_K <- ggplot() +
        theme(
            text = element_text(size = 20),
            panel.background = element_rect(fill = "white", colour = "white"),
            panel.grid.major = element_line(colour = "white"),
            panel.grid.minor = element_line(colour = "white")
        )
    ########################################################  TEST - END
    if (length(nParticles) < nIter) nParticles[(length(nParticles) + 1):nIter] <- nParticles[length(nParticles)]
    parameters_id <- colnames(parameters_initial)
    for (iteration in 1:nIter) {
        cat(paste0("Iteration ", iteration, "...\n"))
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
        ######################################################  TEST - BEGIN
        if (iteration == 1) {
            p_theta <- p_theta + geom_density(data = parameters, aes(x = theta, fill = "Prior Distribution", color = "Prior Distribution"), alpha = 0.2, size = 2)
        }
        ########################################################  TEST - END
        #   Simulate statistics
        reference <- model(parameters)
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
        ######################################################  TEST - BEGIN
        w <- ABCRF_weights
        colnames(w) <- paste0("weights_", colnames(w))
        w$iteration <- iteration
        p_theta <- p_theta +
            geom_density(
                data = cbind(parameters, w),
                aes(x = theta, weight = weights_theta, fill = paste0("Iter-", iteration), color = paste0("Iter-", iteration)), alpha = 0.2, size = 2
            )
        p_K <- p_K +
            geom_density(
                data = reference,
                aes(x = K, fill = paste0("Iter-", iteration), color = paste0("Iter-", iteration)), alpha = 0.2, size = 2
            )
        ########################################################  TEST - END
    }
    ######################################################  TEST - BEGIN
    # cat(paste0("\n\n\n\n\nChecking original ABC-RF with training set of same size...\n"))
    # parameters <- parameters_initial[1:sum(nParticles[1:nIter]), ]
    # reference <- model(parameters)
    # cat(paste0("\n\n\n\n\n"))
    # ABCRF_weights <- data.frame(matrix(NA, nrow = sum(nParticles[1:nIter]), ncol = 0))
    # for (parameter_id in parameters_id) {
    #     mini_reference <- reference[, c(parameter_id, colnames(reference)[!colnames(reference) %in% parameters_id])]
    #     RFmodel <- regAbcrf(
    #         formula = as.formula(paste0(parameter_id, " ~ .")),
    #         data = mini_reference,
    #         paral = parallel
    #     )
    #     posterior_gamma_RF <- predict(
    #         object = RFmodel,
    #         obs = target,
    #         training = mini_reference,
    #         paral = parallel,
    #         rf.weights = T
    #     )
    #     ABCRF_weights[, parameter_id] <- posterior_gamma_RF$weights
    # }
    # w <- ABCRF_weights
    # colnames(w) <- paste0("weights_", colnames(w))
    # w$iteration <- iteration
    # p_theta <- p_theta +
    #     geom_density(
    #         data = cbind(parameters, w),
    #         aes(x = theta, weight = weights_theta, fill = "ABC-RF with same training size", color = "ABC-RF with same training size"), alpha = 0.2, size = 2
    #     )





    p_theta <- p_theta +
        scale_fill_manual(values = c(
            "Ground truth" = "black",
            "Prior Distribution" = "gray",
            "Iter-1" = "purple",
            "Iter-2" = "blue",
            "Iter-3" = "cyan",
            "Iter-4" = "green",
            "Iter-5" = "yellow",
            "Iter-6" = "orange",
            "Iter-7" = "red",
            "ABC-RF with same training size" = "brown"
        ), name = "Legend") +
        scale_colour_manual(values = c(
            "Ground truth" = "black",
            "Prior Distribution" = "gray",
            "Iter-1" = "purple",
            "Iter-2" = "blue",
            "Iter-3" = "cyan",
            "Iter-4" = "green",
            "Iter-5" = "yellow",
            "Iter-6" = "orange",
            "Iter-7" = "red",
            "ABC-RF with same training size" = "brown"
        ), name = "Legend") +
        guides(fill = "none")
    p_K <- p_K +
        scale_fill_manual(values = c(
            "Ground truth" = "black",
            "Prior Distribution" = "gray",
            "Iter-1" = "purple",
            "Iter-2" = "blue",
            "Iter-3" = "cyan",
            "Iter-4" = "green",
            "Iter-5" = "yellow",
            "Iter-6" = "orange",
            "Iter-7" = "red",
            "ABC-RF with same training size" = "brown"
        ), name = "Legend") +
        scale_colour_manual(values = c(
            "Ground truth" = "black",
            "Prior Distribution" = "gray",
            "Iter-1" = "purple",
            "Iter-2" = "blue",
            "Iter-3" = "cyan",
            "Iter-4" = "green",
            "Iter-5" = "yellow",
            "Iter-6" = "orange",
            "Iter-7" = "red",
            "ABC-RF with same training size" = "brown"
        ), name = "Legend") +
        guides(fill = "none")
    png("TEST.png", res = 150, width = 16, height = 8, units = "in", pointsize = 12)
    grid.arrange(p_theta, ncol = 1)
    dev.off()
    png("TEST_stats.png", res = 150, width = 16, height = 8, units = "in", pointsize = 12)
    grid.arrange(p_K, ncol = 1)
    dev.off()
    ########################################################  TEST - END
}
