plot_joint <- function(outputdata,
                       para_id) {
    library(MASS)
    #   Plot multidimensional contours
    nIter <- outputdata[["nIter"]]
    nParticles <- outputdata[["nParticles"]]
    method <- outputdata[["method"]]
    for (iter in 1:nIter) {
        if (method == "smcrf-multi-param") {
            weights <- outputdata[[iter]]$weights[, 1]
            nParticle <- nParticles[iter]
            indx.sim <- sample(1:nParticle, size = 1000, replace = T, prob = weights)
            Ysim <- outputdata[[iter]]$parameters[indx.sim, ]
            para1 <- Ysim[, para_id[1]]
            para2 <- Ysim[, para_id[2]]
            png(paste0("smcdrf_joint_", para_id[1], "_", para_id[2], "_iter_", iter, ".png"))
        } else {
            weights.1 <- outputdata[[iter]]$weights[, para_id[1]]
            weights.2 <- outputdata[[iter]]$weights[, para_id[2]]
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
