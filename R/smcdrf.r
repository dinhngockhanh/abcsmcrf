smcdrf <- function(target,
                   model,
                   n_samples_per_parameter_set,
                   nNoise,
                   perturb,
                   range,
                   parameters_initial,
                   nIter,
                   nParticles,
                   parallel,
                   splitting.rule = "CART",
                   ...) {
    library(drf)
    if (length(nParticles) < nIter) nParticles[(length(nParticles) + 1):nIter] <- nParticles[length(nParticles)]
    parameters_id <- colnames(parameters_initial)
    SMCDRF <- list()
    for (iteration in 1:nIter) {
        #   Sample parameters for this round of iteration
        if (iteration == 1) {
            parameters <- parameters_initial[1:nParticles[iteration], ]
        } else {
            # sample_new <- sample(c(1:nrow(parameters)), size = nParticles[iteration], prob = weights, replace = T)
            # parameters_new <- parameters[sample_new, ]
            # parameters <- perturb(parameters_new)
            # rownames(parameters) <- NULL
            parameters_new <- data.frame(matrix(NA, nrow = nParticles[iteration], ncol = 0))
            for (parameter_id in parameters_id) {
                para_min <- range$min[which(range$parameter == parameter_id)]
                para_max <- range$max[which(range$parameter == parameter_id)]
                # parameters_new[, parameter_id] <-
                #     sample(parameters[, parameter_id], size = nParticles[iteration], prob = weights, replace = T)
                for (row in 1:nParticles[iteration]) {
                    repeat{
                        sampled_para <- sample(parameters[, parameter_id], size = 1, prob = weights, replace = T)
                        perturb_para <- perturb(sampled_para)
                        if (perturb_para >= para_min & perturb_para <= para_max) {
                            parameters_new[row, parameter_id] <- perturb_para
                            break
                        }
                    }
                }
            }
            parameters <- parameters_new
        }
        #   Simulate statistics
        reference <- model(parameters = parameters, n_samples_per_parameter_set = n_samples_per_parameter_set, nNoise = nNoise)
        #   Run DRF for all parameters
        Xdrf <- reference[, colnames(reference)[!colnames(reference) %in% parameters_id]]
        Ydrf <- reference[, parameters_id]
        drfmodel <- drf(Xdrf, Ydrf, splitting.rule = splitting.rule)
        drf_pred <- predict(drfmodel, target)
        weights <- as.vector(get_sample_weights(drfmodel, target))
        DRF_weights <- matrix(rep(weights, length(parameters_id)), ncol = length(parameters_id))
        #   Save SMC-RF results from this iteration
        SMCDRF_iteration <- list()
        SMCDRF_iteration$parameters <- Ydrf
        SMCDRF_iteration$statistics <- Xdrf
        SMCDRF_iteration$weights <- DRF_weights
        SMCDRF_iteration$rf_model <- drfmodel
        SMCDRF_iteration$rf_predict <- drf_pred
        SMCDRF[[iteration]] <- SMCDRF_iteration
    }
    SMCDRF[["method"]] <- "smcrf-multi-param"
    SMCDRF[["nIter"]] <- nIter
    SMCDRF[["nParticles"]] <- nParticles
    SMCDRF[["target"]] <- target
    return(SMCDRF)
}
