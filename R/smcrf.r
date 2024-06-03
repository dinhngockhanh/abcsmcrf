smcrf <- function(method = "smcrf-single-param",
                  statistics_target,
                  model,
                  perturb,
                  range = NULL,
                  parameters_initial,
                  nParticles,
                  parallel,
                  n_cores = NULL,
                  ...) {
    if (method == "smcrf-single-param") {
        return(smcrf_single_param(
            statistics_target = statistics_target,
            model = model,
            perturb = perturb,
            range = range,
            parameters_initial = parameters_initial,
            nParticles = nParticles,
            parallel = parallel,
            n_cores = n_cores,
            ...
        ))
    } else if (method == "smcrf-multi-param") {
        return(smcrf_multi_param(
            statistics_target = statistics_target,
            model = model,
            perturb = perturb,
            range = range,
            parameters_initial = parameters_initial,
            nParticles = nParticles,
            parallel = parallel,
            n_cores = n_cores,
            ...
        ))
    } else {
        stop("Invalid method for smcrf")
    }
}

# smcrf_single_param_old <- function(statistics_target,
#                                model,
#                                perturb,
#                                range = NULL,
#                                parameters_initial,
#                                nParticles,
#                                parallel,
#                                n_cores = NULL,
#                                ...) {
#     library(abcrf)
#     nIterations <- length(nParticles)
#     parameters_ids <- colnames(parameters_initial)
#     SMCRF <- list()
#     nSimulations <<- 0
#     cat("\n\n==============================================================================================================================================\n")
#     for (iteration in 1:(nIterations + 1)) {
#         if (iteration == (nIterations + 1)) {
#             cat(paste0("\n\n++++++++++++++++++++++++\nSimulation count: ", nSimulations, "\n++++++++++++++++++++++++\n\n\n"))
#             cat("\n\nSIMULATING STATISTICS FROM FINAL POSTERIOR DISTRIBUTION...\n")
#         } else {
#             cat(paste0("SMC-RF FOR SINGLE PARAMETERS: iteration ", iteration, "...\n"))
#         }
#         #---Sample prior parameters for this round of iteration...
#         if (iteration == 1) {
#             #   ... For iteration 1: sample from initial parameters
#             parameters <- data.frame(parameters_initial[1:nParticles[iteration], ])
#             colnames(parameters) <- parameters_ids
#         } else {
#             #   ... For later iterations:
#             ifelse(iteration == (nIterations + 1), nrow <- nParticles[nIterations], nrow <- nParticles[iteration])
#             parameters_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
#             colnames(parameters_next) <- parameters_ids
#             for (parameter_id in parameters_ids) {
#                 invalid_indices <- 1:nrow
#                 while (length(invalid_indices) > 0) {
#                     #   Sample parameters from previous posterior distribution
#                     parameter_replace <- data.frame(sample(parameters[, parameter_id], size = length(invalid_indices), prob = ABCRF_weights[, parameter_id], replace = TRUE))
#                     colnames(parameter_replace) <- parameter_id
#                     #   Perturb parameters
#                     if (iteration < (nIterations + 1)) parameter_replace <- perturb(parameters = parameter_replace)
#                     #   Check if parameters are within range, otherwise redo the failed parameters
#                     parameters_next[[parameter_id]][invalid_indices] <- parameter_replace[[parameter_id]]
#                     if (is.null(range)) {
#                         invalid_indices <- which(is.na(parameters_next[, parameter_id]))
#                     } else {
#                         invalid_indices <- which(
#                             is.na(parameters_next[, parameter_id]) |
#                                 parameters_next[, parameter_id] < range$min[which(range$parameter == parameter_id)] |
#                                 parameters_next[, parameter_id] > range$max[which(range$parameter == parameter_id)]
#                         )
#                     }
#                 }
#             }
#             parameters <- parameters_next
#         }
#         #---Simulate statistics
#         cat("Making model simulations...\n")
#         reference <- model(parameters = parameters)
#         statistics <- data.frame(reference[, colnames(reference)[!colnames(reference) %in% parameters_ids]])
#         colnames(statistics) <- colnames(reference)[!colnames(reference) %in% parameters_ids]
#         if (iteration == (nIterations + 1)) {
#             SMCRF_iteration <- list()
#             SMCRF_iteration$reference <- reference
#             SMCRF_iteration$parameters <- parameters
#             SMCRF_iteration$statistics <- statistics
#             SMCRF[[paste0("Iteration_", iteration)]] <- SMCRF_iteration
#             break
#         }
#         #---Run ABCRF for each parameter
#         cat("Performing Random Forest prediction...\n")
#         ABCRF_weights <- data.frame(matrix(NA, nrow = nParticles[iteration], ncol = 0))
#         for (parameter_id in parameters_ids) {
#             mini_reference <- reference[, c(parameter_id, colnames(reference)[!colnames(reference) %in% parameters_ids])]
#             RFmodel <- regAbcrf(
#                 formula = as.formula(paste0(parameter_id, " ~ .")),
#                 data = mini_reference,
#                 paral = parallel
#             )
#             posterior_gamma_RF <- predict(
#                 object = RFmodel,
#                 obs = statistics_target,
#                 training = mini_reference,
#                 paral = parallel,
#                 rf.weights = T
#             )
#             ABCRF_weights[, parameter_id] <- posterior_gamma_RF$weights
#         }
#         cat("\n\n")
#         #---Save SMC-RF results from this iteration
#         SMCRF_iteration <- list()
#         SMCRF_iteration$reference <- reference
#         SMCRF_iteration$parameters <- parameters
#         SMCRF_iteration$statistics <- statistics
#         SMCRF_iteration$weights <- ABCRF_weights
#         SMCRF_iteration$rf_model <- RFmodel
#         SMCRF_iteration$rf_predict <- posterior_gamma_RF
#         SMCRF[[paste0("Iteration_", iteration)]] <- SMCRF_iteration
#     }
#     SMCRF[["method"]] <- "smcrf-single-param"
#     SMCRF[["nIterations"]] <- nIterations
#     SMCRF[["nParticles"]] <- nParticles
#     SMCRF[["statistics_target"]] <- statistics_target
#     SMCRF[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
#     SMCRF[["statistics_labels"]] <- data.frame(ID = colnames(reference)[!colnames(reference) %in% parameters_ids])
#     return(SMCRF)
# }

# smcrf_multi_param_old <- function(statistics_target,
#                               model,
#                               perturb,
#                               range = NULL,
#                               parameters_initial,
#                               nParticles,
#                               parallel,
#                               n_cores = NULL,
#                               splitting.rule = "CART",
#                               ...) {
#     library(drf)
#     nIterations <- length(nParticles)
#     parameters_ids <- colnames(parameters_initial)
#     SMCDRF <- list()
#     nSimulations <<- 0
#     for (iteration in 1:(nIterations + 1)) {
#         if (iteration == (nIterations + 1)) {
#             cat(paste0("\n\n++++++++++++++++++++++++\nSimulation count: ", nSimulations, "\n++++++++++++++++++++++++\n\n\n"))
#             cat("\n\nSIMULATING STATISTICS FROM FINAL POSTERIOR DISTRIBUTION...\n")
#         } else {
#             cat(paste0("SMC-RF FOR MULTIPLE PARAMETERS: iteration ", iteration, "...\n"))
#         }
#         #---Sample prior parameters for this round of iteration...
#         if (iteration == 1) {
#             #   ... For iteration 1: sample from initial parameters
#             parameters <- data.frame(parameters_initial[1:nParticles[iteration], ])
#             colnames(parameters) <- parameters_ids
#         } else {
#             #   ... For later iterations:
#             ifelse(iteration == (nIterations + 1), nrow <- nParticles[nIterations], nrow <- nParticles[iteration])
#             parameters_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
#             colnames(parameters_next) <- parameters_ids
#             invalid_indices <- 1:nrow
#             while (length(invalid_indices) > 0) {
#                 #   Sample parameters from previous posterior distribution
#                 parameter_replace <- data.frame(parameters[sample(nrow(parameters), size = length(invalid_indices), prob = DRF_weights[, 1], replace = T), ])
#                 colnames(parameter_replace) <- parameters_ids
#                 #   Perturb parameters
#                 if (iteration < (nIterations + 1)) parameter_replace <- perturb(parameters = parameter_replace)
#                 #   Check if parameters are within range, otherwise redo the failed parameters
#                 parameters_next[invalid_indices, ] <- parameter_replace
#                 if (is.null(range)) {
#                     invalid_indices <- which(apply(parameters_next, 1, function(x) any(is.na(x))))
#                 } else {
#                     invalid_indices <- c()
#                     for (parameter_id in parameters_ids) {
#                         invalid_indices <- union(invalid_indices, which(
#                             is.na(parameters_next[, parameter_id]) |
#                                 parameters_next[, parameter_id] < range$min[which(range$parameter == parameter_id)] |
#                                 parameters_next[, parameter_id] > range$max[which(range$parameter == parameter_id)]
#                         ))
#                     }
#                 }
#             }
#             parameters <- parameters_next
#         }
#         #---Simulate statistics
#         cat("Making model simulations...\n")
#         reference <- model(parameters = parameters)
#         statistics <- data.frame(reference[, colnames(reference)[!colnames(reference) %in% parameters_ids]])
#         colnames(statistics) <- colnames(reference)[!colnames(reference) %in% parameters_ids]
#         if (iteration == (nIterations + 1)) {
#             SMCDRF_iteration <- list()
#             SMCDRF_iteration$reference <- reference
#             SMCDRF_iteration$parameters <- parameters
#             SMCDRF_iteration$statistics <- statistics
#             SMCDRF[[paste0("Iteration_", iteration)]] <- SMCDRF_iteration
#             break
#         }
#         #---Run DRF for all parameter
#         cat("Performing Random Forest prediction...\n")
#         Xdrf <- statistics
#         Ydrf <- reference[, parameters_ids]
#         print("==============================")
#         print(min(ceiling(sqrt(ncol(Xdrf)) + 20), ncol(Xdrf)))
#         print("==============================")
#         # drfmodel <- drf(Xdrf, Ydrf, splitting.rule = splitting.rule)
#         drfmodel <- drf(Xdrf, Ydrf, mtry = 5, splitting.rule = splitting.rule)
#         def_pred <- predict(drfmodel, statistics_target)
#         DRF_weights <- as.vector(get_sample_weights(drfmodel, statistics_target))
#         DRF_weights <- data.frame(matrix(rep(DRF_weights, length(parameters_ids)), ncol = length(parameters_ids)))
#         colnames(DRF_weights) <- parameters_ids
#         #---Save SMC-RF results from this iteration
#         SMCDRF_iteration <- list()
#         SMCDRF_iteration$reference <- reference
#         SMCDRF_iteration$parameters <- parameters
#         SMCDRF_iteration$statistics <- statistics
#         SMCDRF_iteration$weights <- DRF_weights
#         SMCDRF_iteration$rf_model <- drfmodel
#         SMCDRF_iteration$rf_predict <- def_pred
#         SMCDRF[[paste0("Iteration_", iteration)]] <- SMCDRF_iteration
#     }
#     SMCDRF[["method"]] <- "smcrf-multi-param"
#     SMCDRF[["nIterations"]] <- nIterations
#     SMCDRF[["nParticles"]] <- nParticles
#     SMCDRF[["statistics_target"]] <- statistics_target
#     SMCDRF[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
#     SMCDRF[["statistics_labels"]] <- data.frame(ID = colnames(reference)[!colnames(reference) %in% parameters_ids])
#     return(SMCDRF)
# }

smcrf_single_param <- function(statistics_target,
                               model,
                               perturb,
                               range = NULL,
                               parameters_initial,
                               nParticles,
                               parallel,
                               n_cores = NULL,
                               ...) {
    library(abcrf)
    nIterations <- length(nParticles)
    parameters_ids <- colnames(parameters_initial)
    SMCRF <- list()
    nSimulations <<- 0
    cat("\n\n==============================================================================================================================================\n")
    for (iteration in 1:(nIterations + 1)) {
        if (iteration == (nIterations + 1)) {
            cat(paste0("\n\n++++++++++++++++++++++++\nSimulation count: ", nSimulations, "\n++++++++++++++++++++++++\n\n\n"))
            cat("\n\nSIMULATING STATISTICS FROM FINAL POSTERIOR DISTRIBUTION...\n")
        } else {
            cat(paste0("SMC-RF FOR SINGLE PARAMETERS: iteration ", iteration, "...\n"))
        }
        #---Sample prior parameters for this round of iteration...
        if (iteration == 1) {
            #   ... For iteration 1: sample from initial parameters
            parameters <- data.frame(parameters_initial[1:nParticles[iteration], ])
            colnames(parameters) <- parameters_ids
        } else {
            #   ... For later iterations:
            ifelse(iteration == (nIterations + 1), nrow <- nParticles[nIterations], nrow <- nParticles[iteration])
            parameters_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
            colnames(parameters_next) <- parameters_ids
            for (parameter_id in parameters_ids) {
                invalid_indices <- 1:nrow
                while (length(invalid_indices) > 0) {
                    #   Sample parameters from previous posterior distribution
                    parameter_replace <- data.frame(sample(parameters[, parameter_id], size = length(invalid_indices), prob = ABCRF_weights[, parameter_id], replace = TRUE))

                    colnames(parameter_replace) <- parameter_id
                    #   Perturb parameters
                    if (iteration < (nIterations + 1)) parameter_replace <- perturb(parameters = parameter_replace)
                    #   Check if parameters are within range, otherwise redo the failed parameters
                    parameters_next[[parameter_id]][invalid_indices] <- parameter_replace[[parameter_id]]
                    if (is.null(range)) {
                        invalid_indices <- which(is.na(parameters_next[, parameter_id]))
                    } else {
                        invalid_indices <- which(
                            is.na(parameters_next[, parameter_id]) |
                                parameters_next[, parameter_id] < range$min[which(range$parameter == parameter_id)] |
                                parameters_next[, parameter_id] > range$max[which(range$parameter == parameter_id)]
                        )
                    }
                }
            }
            parameters <- parameters_next
        }
        #---Simulate statistics
        cat("Making model simulations...\n")
        reference <- model(parameters = parameters)
        # statistics <- data.frame(reference[, colnames(reference)[!colnames(reference) %in% parameters_ids]])
        # colnames(statistics) <- colnames(reference)[!colnames(reference) %in% parameters_ids]
        # ################################################
        ################################################
        invalid_rows <- rowSums(is.na(reference)) == ncol(reference)
        while (any(invalid_rows)) {
            print("++++++")
            print(sum(invalid_rows))
            print(reference[invalid_rows, ])
            print("++++++")
            if (iteration == 1) {
                parameters <- data.frame(parameters_initial[nrow(parameters) + 1:(nrow(parameters) + sum(invalid_rows)), ])
                colnames(parameters) <- parameters_ids
                reference[invalid_rows, ] <- model(parameters = parameters)
                invalid_rows <- rowSums(is.na(reference)) == ncol(reference)
            } else {
                #   ... For later iterations:
                # ifelse(iteration == (nIterations + 1), nrow <- nParticles[nIterations], nrow <- nParticles[iteration])
                parameters_next <- data.frame(matrix(NA, nrow = sum(invalid_rows), ncol = length(parameters_ids)))
                colnames(parameters_next) <- parameters_ids
                parameters_tmp <- parameters
                for (parameter_id in parameters_ids) {
                    # print(sample(nrow(parameters[, parameter_id]), size = sum(invalid_rows), prob = ABCRF_weights[, parameter_id], replace = T))
                    # parameters_tmp <- parameters[, parameter_id]
                    # data.frame(parameters[sample(length(parameters[, parameter_id]), size = sum(invalid_rows), prob = ABCRF_weights[, parameter_id], replace = T), parameter_id])
                    # parameters_next <- data.frame(matrix(NA, nrow = sum(invalid_rows), ncol = length(parameters_ids)))
                    # colnames(parameters_next) <- parameters_id
                    invalid_indices <- 1:sum(invalid_rows)
                    while (length(invalid_indices) > 0) {
                        #   Sample parameters from previous posterior distribution
                        parameter_replace <- data.frame(parameters_tmp[sample(nrow(parameters_tmp), size = length(invalid_indices), prob = ABCRF_weights[, parameter_id], replace = T), parameter_id])
                        colnames(parameter_replace) <- parameter_id
                        #   Perturb parameters
                        if (iteration < (nIterations + 1)) parameter_replace <- perturb(parameters = parameter_replace)
                        #   Check if parameters are within range, otherwise redo the failed parameters
                        parameters_next[[parameter_id]][invalid_indices] <- parameter_replace[[parameter_id]]
                        if (is.null(range)) {
                            invalid_indices <- which(is.na(parameters_next[, parameter_id]))
                        } else {
                            invalid_indices <- which(
                                is.na(parameters_next[, parameter_id]) |
                                    parameters_next[, parameter_id] < range$min[which(range$parameter == parameter_id)] |
                                    parameters_next[, parameter_id] > range$max[which(range$parameter == parameter_id)]
                            )
                        }
                    }
                }
                parameters_tmp <- parameters_next
                colnames(parameters_tmp) <- parameters_ids
                reference[invalid_rows, ] <- model(parameters = parameters_tmp)
                invalid_rows <- rowSums(is.na(reference)) == ncol(reference)
            }
        }
        statistics <- data.frame(reference[, colnames(reference)[!colnames(reference) %in% parameters_ids]])
        parameters <- data.frame(reference[, colnames(reference)[colnames(reference) %in% parameters_ids]])
        colnames(parameters) <- parameters_ids
        colnames(statistics) <- colnames(reference)[!colnames(reference) %in% parameters_ids]
        ################################################
        if (iteration == (nIterations + 1)) {
            SMCRF_iteration <- list()
            SMCRF_iteration$reference <- reference
            SMCRF_iteration$parameters <- parameters
            SMCRF_iteration$statistics <- statistics
            SMCRF[[paste0("Iteration_", iteration)]] <- SMCRF_iteration
            break
        }
        #---Run ABCRF for each parameter
        cat("Performing Random Forest prediction...\n")
        ABCRF_weights <- data.frame(matrix(NA, nrow = nParticles[iteration], ncol = 0))
        RFmodels <- list()
        posterior_gamma_RFs <- list()
        for (parameter_id in parameters_ids) {
            mini_reference <- reference[, c(parameter_id, colnames(reference)[!colnames(reference) %in% parameters_ids])]
            colnames(mini_reference)[1] <- "para"
            f <- as.formula("para ~.")
            RFmodel <- regAbcrf(
                formula = f,
                data = mini_reference,
                paral = parallel,
                ...
            )
            posterior_gamma_RF <- predict(
                object = RFmodel,
                obs = statistics_target,
                training = mini_reference,
                paral = parallel,
                rf.weights = T
            )
            ABCRF_weights[, parameter_id] <- posterior_gamma_RF$weights
            RFmodels[[parameter_id]] <- RFmodel
            posterior_gamma_RFs[[parameter_id]] <- posterior_gamma_RF
        }
        cat("\n\n")
        #---Save SMC-RF results from this iteration
        SMCRF_iteration <- list()
        SMCRF_iteration$reference <- reference
        SMCRF_iteration$parameters <- parameters
        SMCRF_iteration$statistics <- statistics
        SMCRF_iteration$weights <- ABCRF_weights
        SMCRF_iteration$rf_model <- RFmodels
        SMCRF_iteration$rf_predict <- posterior_gamma_RFs
        SMCRF[[paste0("Iteration_", iteration)]] <- SMCRF_iteration
    }
    SMCRF[["method"]] <- "smcrf-single-param"
    SMCRF[["nIterations"]] <- nIterations
    SMCRF[["nParticles"]] <- nParticles
    SMCRF[["statistics_target"]] <- statistics_target
    SMCRF[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
    SMCRF[["statistics_labels"]] <- data.frame(ID = colnames(reference)[!colnames(reference) %in% parameters_ids])
    return(SMCRF)
}

smcrf_multi_param <- function(statistics_target,
                              model,
                              perturb,
                              range = NULL,
                              parameters_initial,
                              nParticles,
                              parallel,
                              n_cores = NULL,
                              splitting.rule = "CART",
                              ...) {
    library(drf)
    nIterations <- length(nParticles)
    parameters_ids <- colnames(parameters_initial)
    SMCDRF <- list()
    nSimulations <<- 0
    for (iteration in 1:(nIterations + 1)) {
        if (iteration == (nIterations + 1)) {
            cat(paste0("\n\n++++++++++++++++++++++++\nSimulation count: ", nSimulations, "\n++++++++++++++++++++++++\n\n\n"))
            cat("\n\nSIMULATING STATISTICS FROM FINAL POSTERIOR DISTRIBUTION...\n")
        } else {
            cat(paste0("SMC-RF FOR MULTIPLE PARAMETERS: iteration ", iteration, "...\n"))
        }
        #---Sample prior parameters for this round of iteration...
        if (iteration == 1) {
            #   ... For iteration 1: sample from initial parameters

            parameters <- data.frame(parameters_initial[1:nParticles[iteration], ])
            colnames(parameters) <- parameters_ids
        } else {
            #   ... For later iterations:
            ifelse(iteration == (nIterations + 1), nrow <- nParticles[nIterations], nrow <- nParticles[iteration])
            parameters_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
            colnames(parameters_next) <- parameters_ids
            invalid_indices <- 1:nrow
            while (length(invalid_indices) > 0) {
                #   Sample parameters from previous posterior distribution
                parameter_replace <- data.frame(parameters[sample(nrow(parameters), size = length(invalid_indices), prob = DRF_weights[, 1], replace = T), ])
                colnames(parameter_replace) <- parameters_ids
                #   Perturb parameters
                if (iteration < (nIterations + 1)) parameter_replace <- perturb(parameters = parameter_replace)
                #   Check if parameters are within range, otherwise redo the failed parameters
                parameters_next[invalid_indices, ] <- parameter_replace
                if (is.null(range)) {
                    invalid_indices <- which(apply(parameters_next, 1, function(x) any(is.na(x))))
                } else {
                    invalid_indices <- c()
                    for (parameter_id in parameters_ids) {
                        invalid_indices <- union(invalid_indices, which(
                            is.na(parameters_next[, parameter_id]) |
                                parameters_next[, parameter_id] < range$min[which(range$parameter == parameter_id)] |
                                parameters_next[, parameter_id] > range$max[which(range$parameter == parameter_id)]
                        ))
                    }
                }
            }
            parameters <- parameters_next
        }
        #---Simulate statistics
        cat("Making model simulations...\n")
        reference <- model(parameters = parameters)
        # statistics <- data.frame(reference[, colnames(reference)[!colnames(reference) %in% parameters_ids]])
        ################################################
        ################################################
        invalid_rows <- rowSums(is.na(reference)) == ncol(reference)

        while (any(invalid_rows)) {
            print("++++++")
            print(sum(invalid_rows))
            print(reference[invalid_rows, ])
            print("++++++")
            if (iteration == 1) {
                parameters <- data.frame(parameters_initial[(nrow(parameters) + 1):(nrow(parameters) + sum(invalid_rows)), ])
                colnames(parameters) <- parameters_ids
                print("=======")
                reference[invalid_rows, ] <- model(parameters = parameters)
                invalid_rows <- rowSums(is.na(reference)) == ncol(reference)
            } else {
                # parameters_tmp <- data.frame(parameters[sample(nrow(parameters), size = sum(invalid_rows), prob = DRF_weights[, 1], replace = T), ])
                parameters_tmp <- parameters
                #   ... For later iterations:
                parameters_next <- data.frame(matrix(NA, nrow = sum(invalid_rows), ncol = length(parameters_ids)))
                colnames(parameters_next) <- parameters_ids
                invalid_indices <- 1:sum(invalid_rows)
                while (length(invalid_indices) > 0) {
                    print("###################################")
                    print(reference[invalid_rows, ])
                    print("###################################")
                    #   Sample parameters from previous posterior distribution
                    parameter_replace <- data.frame(parameters_tmp[sample(nrow(parameters_tmp), size = length(invalid_indices), prob = DRF_weights[, 1], replace = T), ])
                    colnames(parameter_replace) <- parameters_ids
                    #   Perturb parameters
                    if (iteration < (nIterations + 1)) parameter_replace <- perturb(parameters = parameter_replace)
                    #   Check if parameters are within range, otherwise redo the failed parameters
                    parameters_next[invalid_indices, ] <- parameter_replace
                    if (is.null(range)) {
                        invalid_indices <- which(apply(parameters_next, 1, function(x) any(is.na(x))))
                    } else {
                        invalid_indices <- c()
                        for (parameter_id in parameters_ids) {
                            invalid_indices <- union(invalid_indices, which(
                                is.na(parameters_next[, parameter_id]) |
                                    parameters_next[, parameter_id] < range$min[which(range$parameter == parameter_id)] |
                                    parameters_next[, parameter_id] > range$max[which(range$parameter == parameter_id)]
                            ))
                        }
                    }
                }
                parameters_tmp <- parameters_next
                colnames(parameters_tmp) <- parameters_ids
                reference[invalid_rows, ] <- model(parameters = parameters_tmp)
                invalid_rows <- rowSums(is.na(reference)) == ncol(reference)
            }
        }
        statistics <- data.frame(reference[, colnames(reference)[!colnames(reference) %in% parameters_ids]])
        parameters <- data.frame(reference[, colnames(reference)[colnames(reference) %in% parameters_ids]])
        colnames(parameters) <- parameters_ids
        ################################################
        ################################################
        colnames(statistics) <- colnames(reference)[!colnames(reference) %in% parameters_ids]
        if (iteration == (nIterations + 1)) {
            SMCDRF_iteration <- list()
            SMCDRF_iteration$reference <- reference
            SMCDRF_iteration$parameters <- parameters
            SMCDRF_iteration$statistics <- statistics
            SMCDRF[[paste0("Iteration_", iteration)]] <- SMCDRF_iteration
            break
        }
        #---Run DRF for all parameter
        cat("Performing Random Forest prediction...\n")
        Xdrf <- statistics
        Ydrf <- reference[, parameters_ids]
        drfmodel <- drf(Xdrf, Ydrf, splitting.rule = splitting.rule, ...)
        # drfmodel <- drf(Xdrf, Ydrf, mtry = 5, splitting.rule = splitting.rule)
        def_pred <- predict(drfmodel, statistics_target)
        DRF_weights <- as.vector(get_sample_weights(drfmodel, statistics_target))
        DRF_weights <- data.frame(matrix(rep(DRF_weights, length(parameters_ids)), ncol = length(parameters_ids)))
        colnames(DRF_weights) <- parameters_ids
        # ---Save SMC-RF results from this iteration
        SMCDRF_iteration <- list()
        SMCDRF_iteration$reference <- reference
        SMCDRF_iteration$parameters <- parameters
        SMCDRF_iteration$statistics <- statistics
        SMCDRF_iteration$weights <- DRF_weights
        SMCDRF_iteration$rf_model <- drfmodel
        SMCDRF_iteration$rf_predict <- def_pred
        SMCDRF[[paste0("Iteration_", iteration)]] <- SMCDRF_iteration
    }
    SMCDRF[["method"]] <- "smcrf-multi-param"
    SMCDRF[["nIterations"]] <- nIterations
    SMCDRF[["nParticles"]] <- nParticles
    SMCDRF[["statistics_target"]] <- statistics_target
    SMCDRF[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
    SMCDRF[["statistics_labels"]] <- data.frame(ID = colnames(reference)[!colnames(reference) %in% parameters_ids])
    return(SMCDRF)
}
