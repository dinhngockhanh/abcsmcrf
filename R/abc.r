#' @export
abc_rejection <- function(statistics_target,
                          model,
                          parameters_labels,
                          prior_distributions,
                          nParticles,
                          tolerance_quantile,
                          ...) {
    library(EasyABC)
    parameters_ids <- parameters_labels$parameter
    #---Run ABC-rejection
    ABC_REJECTION <- list()
    nSimulations <<- 0
    cat("\n\n==============================================================================================================================================\n")
    cat("ABC-REJECTION...\n")
    model_wrap <- function(x) {
        parameters <- data.frame(matrix(x, nrow = 1))
        colnames(parameters) <- parameters_ids
        as.numeric(model(parameters, parallel = FALSE)[-c(1:ncol(parameters))])
    }
    stat <- as.numeric(statistics_target[1, ])
    cat("Performing ABC-rejection...\n")
    rej_model <- ABC_rejection(
        model = model_wrap, prior = prior_distributions, nb_simul = nParticles,
        summary_stat_target = stat,
        tol = tolerance_quantile,
        ...
    )
    ABC_REJECTION[["Iteration_1"]]$parameters <- data.frame(rej_model$param)
    colnames(ABC_REJECTION[["Iteration_1"]]$parameters) <- parameters_ids
    ABC_REJECTION[["Iteration_1"]]$statistics <- data.frame(rej_model$stats)
    colnames(ABC_REJECTION[["Iteration_1"]]$statistics) <- colnames(statistics_target)[!colnames(statistics_target) %in% parameters_ids]
    ABC_REJECTION[["Iteration_1"]]$weights <- data.frame(matrix(rep(rej_model$weights, length(parameters_ids)), ncol = length(parameters_ids)))
    colnames(ABC_REJECTION[["Iteration_1"]]$weights) <- parameters_ids
    ABC_REJECTION[["Iteration_1"]]$reference <- cbind(ABC_REJECTION[["Iteration_1"]]$parameters, ABC_REJECTION[["Iteration_1"]]$statistics)
    ABC_REJECTION$rej_model <- rej_model
    cat(paste0("\n\n\n\n++++++++++++++++++++++++\nSimulation count: ", nSimulations, "\n++++++++++++++++++++++++\n\n\n\n\n"))
    #---Produce statistics from parameters sampled from the posterior
    cat("Simulating statistics from posterior distribution...\n")
    parameters <- ABC_REJECTION[["Iteration_1"]]$parameters
    weights <- ABC_REJECTION[["Iteration_1"]]$weights
    parameters <- data.frame(parameters[sample(nrow(parameters), size = round(nParticles * tolerance_quantile), prob = weights[, 1], replace = T), ])
    colnames(parameters) <- parameters_ids
    reference <- model(parameters = parameters)
    statistics <- data.frame(reference[, colnames(reference)[!colnames(reference) %in% parameters_ids]])
    colnames(statistics) <- colnames(reference)[!colnames(reference) %in% parameters_ids]
    ABC_REJECTION[["Iteration_2"]]$reference <- reference
    ABC_REJECTION[["Iteration_2"]]$parameters <- parameters
    ABC_REJECTION[["Iteration_2"]]$statistics <- statistics
    ABC_REJECTION[["method"]] <- "abc-rejection"
    ABC_REJECTION[["nParticles"]] <- nParticles
    ABC_REJECTION[["tolerance_quantile"]] <- tolerance_quantile
    ABC_REJECTION[["statistics_target"]] <- statistics_target
    ABC_REJECTION[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
    ABC_REJECTION[["statistics_labels"]] <- data.frame(ID = colnames(statistics_target)[!colnames(statistics_target) %in% parameters_ids])
    return(ABC_REJECTION)
}

#' @export
abc_mcmc <- function(statistics_target,
                     model,
                     parameters_labels,
                     prior_distributions,
                     nParticles,
                     method = "Marjoram",
                     ...) {
    library(EasyABC)
    parameters_ids <- parameters_labels$parameter
    #---Run ABC-MCMC
    AMC_MCMC <- list()
    nSimulations <<- 0
    cat("\n\n==============================================================================================================================================\n")
    cat("ABC-MCMC...\n")
    model_wrap <- function(x) {
        parameters <- data.frame(matrix(x, nrow = 1))
        colnames(parameters) <- parameters_ids
        as.numeric(model(parameters, parallel = FALSE)[-c(1:ncol(parameters))])
    }
    stat <- as.numeric(statistics_target[1, ])
    cat("Performing ABC-MCMC...\n")
    mcmc_model <- ABC_mcmc(
        model = model_wrap, prior = prior_distributions, n_rec = nParticles,
        summary_stat_target = stat,
        method = method,
        ...
    )
    AMC_MCMC[["Iteration_1"]]$parameters <- data.frame(mcmc_model$param)
    colnames(AMC_MCMC[["Iteration_1"]]$parameters) <- parameters_ids
    AMC_MCMC[["Iteration_1"]]$statistics <- data.frame(mcmc_model$stats)
    colnames(AMC_MCMC[["Iteration_1"]]$statistics) <- colnames(statistics_target)[!colnames(statistics_target) %in% parameters_ids]
    AMC_MCMC[["Iteration_1"]]$weights <- data.frame(matrix(1 / nrow(mcmc_model$param), nrow = nrow(mcmc_model$param), ncol = ncol(mcmc_model$param)))
    colnames(AMC_MCMC[["Iteration_1"]]$weights) <- parameters_ids
    AMC_MCMC[["Iteration_1"]]$reference <- cbind(AMC_MCMC[["Iteration_1"]]$parameters, AMC_MCMC[["Iteration_1"]]$statistics)
    AMC_MCMC$mcmc_model <- mcmc_model
    cat(paste0("\n\n\n\n++++++++++++++++++++++++\nSimulation count: ", nSimulations, "\n++++++++++++++++++++++++\n\n\n\n\n"))
    #---Produce statistics from parameters sampled from the posterior
    cat("Simulating statistics from posterior distribution...\n")
    parameters <- AMC_MCMC[["Iteration_1"]]$parameters
    weights <- AMC_MCMC[["Iteration_1"]]$weights
    parameters <- data.frame(parameters[sample(nrow(parameters), size = nParticles, prob = weights[, 1], replace = T), ])
    colnames(parameters) <- parameters_ids
    reference <- model(parameters = parameters)
    statistics <- data.frame(reference[, colnames(reference)[!colnames(reference) %in% parameters_ids]])
    colnames(statistics) <- colnames(reference)[!colnames(reference) %in% parameters_ids]
    AMC_MCMC[["Iteration_2"]]$reference <- reference
    AMC_MCMC[["Iteration_2"]]$parameters <- parameters
    AMC_MCMC[["Iteration_2"]]$statistics <- statistics
    AMC_MCMC[["method"]] <- "abc-mcmc"
    AMC_MCMC[["submethod"]] <- method
    AMC_MCMC[["statistics_target"]] <- statistics_target
    AMC_MCMC[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
    AMC_MCMC[["statistics_labels"]] <- data.frame(ID = colnames(statistics_target)[!colnames(statistics_target) %in% parameters_ids])
    return(AMC_MCMC)
}

#' @export
abc_smc <- function(statistics_target,
                    model,
                    parameters_labels,
                    prior_distributions,
                    nParticles,
                    method = "Beaumont",
                    tolerance,
                    ...) {
    library(EasyABC)
    parameters_ids <- parameters_labels$parameter
    #---Run ABC-MCMC
    AMC_SMC <- list()
    nSimulations <<- 0
    cat("\n\n==============================================================================================================================================\n")
    cat("ABC-SMC...\n")
    model_wrap <- function(x) {
        parameters <- data.frame(matrix(x, nrow = 1))
        colnames(parameters) <- parameters_ids
        as.numeric(model(parameters, parallel = FALSE)[-c(1:ncol(parameters))])
    }
    stat <- as.numeric(statistics_target[1, ])
    cat("Performing ABC-SMC...\n")
    smc_model <- ABC_sequential(
        model = model_wrap, prior = prior_distributions, nb_simul = nParticles,
        summary_stat_target = stat,
        method = method, tolerance_tab = tolerance,
        ...
    )
    AMC_SMC[["Iteration_1"]]$parameters <- data.frame(smc_model$param)
    colnames(AMC_SMC[["Iteration_1"]]$parameters) <- parameters_ids
    AMC_SMC[["Iteration_1"]]$statistics <- data.frame(smc_model$stats)
    colnames(AMC_SMC[["Iteration_1"]]$statistics) <- colnames(statistics_target)[!colnames(statistics_target) %in% parameters_ids]
    AMC_SMC[["Iteration_1"]]$weights <- data.frame(matrix(1 / nrow(smc_model$param), nrow = nrow(smc_model$param), ncol = ncol(smc_model$param)))
    colnames(AMC_SMC[["Iteration_1"]]$weights) <- parameters_ids
    AMC_SMC[["Iteration_1"]]$reference <- cbind(AMC_SMC[["Iteration_1"]]$parameters, AMC_SMC[["Iteration_1"]]$statistics)
    AMC_SMC$smc_model <- smc_model
    cat(paste0("\n\n\n\n++++++++++++++++++++++++\nSimulation count: ", nSimulations, "\n++++++++++++++++++++++++\n\n\n\n\n"))
    #---Produce statistics from parameters sampled from the posterior
    cat("Simulating statistics from posterior distribution...\n")
    parameters <- AMC_SMC[["Iteration_1"]]$parameters
    weights <- AMC_SMC[["Iteration_1"]]$weights
    parameters <- data.frame(parameters[sample(nrow(parameters), size = nParticles, prob = weights[, 1], replace = T), ])
    colnames(parameters) <- parameters_ids
    reference <- model(parameters = parameters)
    statistics <- data.frame(reference[, colnames(reference)[!colnames(reference) %in% parameters_ids]])
    colnames(statistics) <- colnames(reference)[!colnames(reference) %in% parameters_ids]
    AMC_SMC[["Iteration_2"]]$reference <- reference
    AMC_SMC[["Iteration_2"]]$parameters <- parameters
    AMC_SMC[["Iteration_2"]]$statistics <- statistics
    AMC_SMC[["method"]] <- "abc-smc"
    AMC_SMC[["submethod"]] <- method
    AMC_SMC[["statistics_target"]] <- statistics_target
    AMC_SMC[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
    AMC_SMC[["statistics_labels"]] <- data.frame(ID = colnames(statistics_target)[!colnames(statistics_target) %in% parameters_ids])
    return(AMC_SMC)
}
