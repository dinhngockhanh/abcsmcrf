smcdrf_test_2 <- function(target,
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
        para_id <- colnames(parameters_initial)
        for (iteration in 1:nIter) {
            cat(paste0("Iteration ", iteration, "...\n"))
            #   Sample parameters for this round of iteration
            if (iteration == 1) {
            parameters <- parameters_initial[1:nParticles[iteration], ]
            } else {
                # parameters_new <- sample(parameters, size = nParticles[iteration], prob = DRF_weights, replace = T)
                sample_new <- sample(c(1:1000), size = nParticles[iteration], prob = DRF_weights, replace = T)
                parameters_new <- parameters[sample_new,]
                parameters <- perturb(parameters_new)
                rownames(parameters)<-NULL
            }
            ######################################################  TEST - BEGIN
            if (iteration == 1) {
                p_theta <- p_theta + geom_density(data = parameters, aes(x = theta, fill = "Prior Distribution", color = "Prior Distribution"), alpha = 0.2, size = 2)
            }
            ########################################################  TEST - END
            #   Simulate statistics
            reference <- model(parameters)
            #   Run DRF for all parameters
            Xdrf <- reference[, colnames(reference)[!colnames(reference) %in% para_id]]
            Ydrf <- reference[, para_id]
            drfmodel <- drf(Xdrf, Ydrf, mtry = 3, splitting.rule = "CART")
            def_pred <- predict(drfmodel, target)
            DRF_weights <- as.vector(get_sample_weights(drfmodel, target))
            ######################################################  TEST - BEGIN
            w <- DRF_weights
            # colnames(w) <- paste0("weights_", colnames(w))
            # colnames(w) <- c("weights")
            # w$iteration <- iteration
            p_theta <- p_theta +
                geom_density(
                    data = cbind(parameters, w),
                    aes(x = theta, weight = w, fill = paste0("Iter-", iteration), color = paste0("Iter-", iteration)), alpha = 0.2, size = 2
                )
            p_K <- p_K +
                geom_density(
                    data = reference,
                    aes(x = K, fill = paste0("Iter-", iteration), color = paste0("Iter-", iteration)), alpha = 0.2, size = 2
                )
            ########################################################  TEST - END
        }


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
        png("TEST_drf.png", res = 150, width = 16, height = 8, units = "in", pointsize = 12)
        grid.arrange(p_theta, ncol = 1)
        dev.off()
        png("TEST_stats_drf.png", res = 150, width = 16, height = 8, units = "in", pointsize = 12)
        grid.arrange(p_K, ncol = 1)
        dev.off()
        ########################################################  TEST - END
    }


# joint plot

# ground truth
# overlap each iteration
# if method is multi take the first col in weights eachsample is 
# birth death strong correlation   joint should be very pernuced
# birth rate and death rate should be correlated
# use both version for each model
# an potion that specfies which two parameters

# single sample from each columns  vertically and put together


plot_joint <- function(outputdata,
                       para_id){
    library(MASS)
    #   Plot multidimensional contours
    nIter <- outputdata[["nIter"]]
    nParticles <- outputdata[["nParticles"]]
    method <- outputdata[["method"]] 
    for (iter in 1:nIter){
        if (method == "smcrf-multi-param"){
            weights <- outputdata[[iter]]$weights[,1]
            nParticle <- nParticles[iter]
            indx.sim <- sample(1:nParticle, size = 1000, replace = T, prob = weights)
            Ysim <- outputdata[[iter]]$parameters[indx.sim, ]
            para1 <- Ysim[, para_id[1]]
            para2 <- Ysim[, para_id[2]]
            png(paste0("smcdrf_joint_", para_id[1], "_", para_id[2], "_iter_", iter, ".png"))
        }else{
            weights.1 <- outputdata[[iter]]$weights[,para_id[1]]
            weights.2 <- outputdata[[iter]]$weights[,para_id[2]]
            indx.sim.1 <- sample(1:nParticle, size = 1000, replace = T, prob = weights.1)
            indx.sim.2 <- sample(1:nParticle, size = 1000, replace = T, prob = weights.2)
            para1 <- outputdata[[iter]]$parameters[indx.sim.1, para_id[1]]
            para2 <- outputdata[[iter]]$parameters[indx.sim.2, para_id[2]]
            png(paste0("smcabcrf_joint_", para_id[1], "_", para_id[2], "_iter_", iter, ".png"))
        }
        z <- kde2d(para1, para2, n = 50)
        plot(para1, para2, pch = 19)
        filled.contour(z, lwd = 2, add = TRUE, col = hcl.colors(10, "Spectral"))
        dev.off()
    }
}

png(paste0("smcdrf_joint_", para_id[1], "_", para_id[2], ".png"))
for (iter in 1:nIter) {
    weights <- outputdata[[iter]]$weights[, 1]
    nParticle <- nParticles[iter]
    indx.sim <- sample(1:nParticle, size = 1000, replace = TRUE, prob = weights)
    Ysim <- outputdata[[iter]]$parameters[indx.sim, ]
    para1 <- Ysim[, para_id[1]]
    para2 <- Ysim[, para_id[2]]
    data <- data.frame(x = para1, y = para2)
    colnames(data) <- para_id
    # z <- kde2d(para1, para2, n = 50)
    if (iter == 1) {
        p <- ggplot(data, aes_string(x = para_id[1], y = para_id[2])) +
            # geom_point() +
            # geom_density_2d(alpha = 0.3, color = "black") +
            geom_density_2d_filled(alpha = 0.3) +
            # scale_fill_gradient(low = "blue", high = "red") + 
            labs(x = para_id[1], y = para_id[2])
    } else {
        p <- p + geom_density_2d_filled(data = data, aes_string(x = para_id[1], y = para_id[2]), alpha = 0.3)
    }
}
grid.arrange(p, ncol = 1) 
dev.off()
