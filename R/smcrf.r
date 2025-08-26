#' Approximate Bayesian Computation sequential Monte Carlo via random forests
#'
#' \code{\link{smcrf}} uses random forests to find the posterior distribution(s) for one or more parameters in a model.
#' It implements the sequential Monte Carlo framework, where each iteration
#' uses either ABC-RF (functions \code{regAbcrf} and \code{predict} in R package \code{abcrf})
#' or ABC-DRF (functions \code{drf} and \code{predict} in R package \code{drf}) to update the posterior distribution(s).
#'
#' @param method Random forest method to implement in each iteration (\code{"smcrf-single-param"} by default).
#' method = \code{"smcrf-single-param"}: implements ABC-RF for each parameter and results in their marginal posterior distributions.
#' method = \code{"smcrf-multi-param"}: implements ABC-DRF for all parameters and results in the joint posterior distribution.
#' @param statistics_target A dataframe containing statistics from data.
#' Column names are the statistics IDs.
#' \code{\link{smcrf}} only supports one row of statistics.
#' If there are multiple observations, we recommend applying \code{\link{smcrf}} to each row individually.
#' @param statistics_selection A dataframe indicating selection of statistics for fitting individual parameters (only works for method \code{"smcrf-single-param"}; \code{NULL} by default).
#' Each column's name matches one statistic ID, and each row's name matches one parameter ID.
#' The value is 1 if the statistic is used for the parameter, 0 otherwise.
#' @param smcrf_results An existing ABC-SMC-RF result.
#' If provided, smcrf will continue ABC-SMC-RF from the last iteration of the previous run.
#' @param model Model for the statistics.
#' The function must take two inputs: a dataframe \code{parameters} and logic variable \code{parallel}.
#' The model must output a reference table, where each row contains parameters for each simulation and corresponding statistics.
#' The column names of the reference table must match the parameter and statistics IDs.
#' @param rprior Function to generate particles from the prior distribution.
#' The function must take one input: \code{Nparameters}, the number of particles to generate.
#' The output is a dataframe where column names match parameter IDs,
#' and each row contains one parameter set.
#' @param dprior Function to compute the prior density.
#' The function must take two inputs: \code{parameters} and \code{parameter_id}.
#' The dataframe \code{parameters} contains parameter sets in each row, with column names as parameter IDs.
#' The \code{parameter_id} is either \code{"all"} or one of the parameter IDs.
#' The output is a vector of prior probabilities corresponding to rows in \code{parameters},
#' either for the parameter indicated by \code{parameter_id} or jointly for all parameters (if \code{parameter_id} = \code{"all"}).
#' @param perturbation Perturbation method for the parameters.
#' \code{\link{smcrf}} supports \code{perturbation} = \code{"Gaussian"} (default) or \code{"Uniform"}.
#' @param perturbation_parameters A dataframe containing the parameters for the perturbation.
#' Each row corresponds to one iteration, and each column corresponds to one parameter (the column
#' names must match parameter_ids).
#' The values are the normal distribution variances (for \code{perturbation} = \code{"Gaussian"}) or ranges (for \code{perturbation} = \code{"Uniform"}).
#' @param nParticles A vector of particle counts.
#' Each entry indicates the number of simulations (e.g. particles) in the corresponding iteration.
#' @param final_sample A logic variable (\code{TRUE} by default).
#' If \code{final_sample} = \code{TRUE}, parameters will be sampled from the final iteration and saved in iteration \code{length(nParticles)+1}.
#' @param model_redo_if_NA A logic variable (\code{FALSE} by default).
#' If \code{model_redo_if_NA} = \code{TRUE}, the particles where \code{model} returns \code{NA} will be simulated again.
#' @param verbose A logic variable (\code{TRUE} by default).
#' If \code{verbose} = \code{TRUE}, progress messages will be printed during the fitting process.
#' @param parallel A logic variable (\code{FALSE} by default).
#' If \code{parallel} = \code{TRUE}, the ABC-RF functions will be computed in parallel.
#' @param save_model A logic variable (\code{FALSE} by default).
#' If \code{save_model} = \code{TRUE}, the random forest will be saved in \code{\link{smcrf}}'s output.
#' @param save_rds A logic variable (\code{FALSE} by default).
#' If \code{save_rds} = \code{TRUE}, the ABC-SMC-RF results will be saved in an rds file.
#' @param filename_rds A string (\code{"ABCSMCDRF.rds"} by default).
#' If \code{save_rds} = \code{TRUE}, the output from ABC-SMC-(D)RF will be saved in a file with this name.
#' @param ... Additional arguments to be passed to \code{abcrf} or \code{drf}.
#' @return An object \code{smcrf_results} containing the results of the inference.
#' If the posterior distributions have not converged to a satisfactory level,
#' the user may continue with \code{smcrf(smcrf_results = smcrf_results, ...)},
#' in which case ABC-SMC-(D)RF will continue iterating from the last run in \code{smcrf_results}.
#' @export
#' @examples
#' library(abcsmcrf)
#' #--------------------------------------------------------------------
#' #---------------------------ABC-SMC-RF for a model with one parameter
#' #--------------------------------------------------------------------
#' #    Data to be fitted consists of two statistics s1 and s2
#' statistics_target <- data.frame(s1 = 0, s2 = 2)
#' #    We define a parametrized model for the statistics
#' model <- function(parameters) {
#'     statistics <- data.frame(
#'         s1 = parameters$theta - 1 + runif(nrow(parameters), -0.1, 0.1),
#'         s2 = parameters$theta + 1 + runif(nrow(parameters), -0.1, 0.1)
#'     )
#'     cbind(parameters, statistics)
#' }
#' #    and the perturbation parameters for a uniform perturbation
#' perturbation_parameters <- data.frame(
#'     theta = rep(0.1, 2) # vector length is equal to number of ABC-SMC-(D)RF iterations
#' )
#' #    We then define rprior and dprior
#' rprior <- function(Nparameters) {
#'     theta <- runif(Nparameters, -10, 10)
#'     return(data.frame(theta = theta))
#' }
#' dprior <- function(parameters, parameter_id = "theta") {
#'     return(rep(1 / 20, nrow(parameters)))
#' }
#' #    Finally, we run ABC-SMC-RF with 2 iterations, each with 1000 particles
#' smcrf_results <- smcrf(
#'     method = "smcrf-single-param",
#'     statistics_target = statistics_target,
#'     model = model,
#'     rprior = rprior,
#'     dprior = dprior,
#'     perturbation = "Uniform",
#'     perturbation_parameters = perturbation_parameters,
#'     nParticles = c(1000, 1000),
#' )
#' #    Now we examine the posterior distribution of theta
#' posterior_iteration <- paste0("Iteration_", (smcrf_results$nIterations + 1))
#' posterior_theta <- smcrf_results[[posterior_iteration]]$parameters$theta
#' #    We look at the posterior mean of theta
#' theta_mean <- mean(posterior_theta)
#' print(theta_mean)
#' #    Notice that the mean is close to 1, the true value of theta.
#' #    We can also look at the posterior variance of theta
#' theta_var <- var(posterior_theta)
#' print(theta_var)
#' #    We can also continue the ABC-SMC-RF run if the posterior convergence is not satisfactory
#' smcrf_results <- smcrf(
#'     method = "smcrf-single-param",
#'     smcrf_results = smcrf_results,
#'     model = model,
#'     rprior = rprior,
#'     dprior = dprior,
#'     perturbation = "Uniform",
#'     perturbation_parameters = perturbation_parameters,
#'     nParticles = c(1000, 1000),
#' )
#' #    We can look again at the posterior mean and variance of theta
#' posterior_iteration <- paste0("Iteration_", (smcrf_results$nIterations + 1))
#' posterior_theta <- smcrf_results[[posterior_iteration]]$parameters$theta
#' theta_mean <- mean(posterior_theta)
#' print(theta_mean)
#' theta_var <- var(posterior_theta)
#' print(theta_var)
#' #    and notice whether there is any improvement in the posterior distribution
#' #    We can continue the runs of ABC-SMC-(D)RF similarly for the examples below
#' #--------------------------------------------------------------------
#' #---------------------ABC-SMC-RF for a model with multiple parameters
#' #--------------------------------------------------------------------
#' #    Data to be fitted consists of two statistics s1 and s2
#' statistics_target <- data.frame(s1 = 4, s2 = 4)
#' #    We then define a parametrized model for the statistics
#' model <- function(parameters) {
#'     statistics <- data.frame(
#'         s1 = parameters$mu + parameters$theta + runif(nrow(parameters), -0.1, 0.1),
#'         s2 = parameters$mu * parameters$theta + runif(nrow(parameters), -0.1, 0.1)
#'     )
#'     cbind(parameters, statistics)
#' }
#' #    and the perturbation parameters
#' perturbation_parameters <- data.frame(
#'     theta = rep(0.1, 3),
#'     mu = rep(0.1, 3)
#' )
#' #    We define the rprior and dprior functions
#' rprior <- function(Nparameters) {
#'     theta <- runif(Nparameters, -10, 10)
#'     mu <- runif(Nparameters, -10, 10)
#'     return(data.frame(theta = theta, mu = mu))
#' }
#' dprior <- function(parameters, parameter_id = "all") {
#'     probs <- rep(1, nrow(parameters))
#'     if (parameter_id %in% c("all", "theta")) {
#'         probs <- probs * dunif(parameters[["theta"]], -10, 10)
#'     }
#'     if (parameter_id %in% c("all", "mu")) {
#'         probs <- probs * dunif(parameters[["mu"]], -10, 10)
#'     }
#'     return(probs)
#' }
#' #    Finally, we run ABC-SMC-RF with 3 iterations, each with 1000 particles
#' smcrf_results <- smcrf(
#'     method = "smcrf-single-param",
#'     statistics_target = statistics_target,
#'     model = model,
#'     rprior = rprior,
#'     dprior = dprior,
#'     perturbation = "Uniform",
#'     perturbation_parameters = perturbation_parameters,
#'     nParticles = c(1000, 1000, 1000),
#' )
#' #    Now we examine the posterior distribution of each parameter
#' posterior_iteration <- paste0("Iteration_", (smcrf_results$nIterations + 1))
#' posterior_params <- smcrf_results[[posterior_iteration]]$parameters
#' posterior_means <- colMeans(posterior_params)
#' posterior_vars <- var(posterior_params)
#' print(posterior_means)
#' print(posterior_vars)
#' #--------------------------------------------------------------------
#' #--------------------------------ABC-SMC-DRF for a multivariate model
#' #--------------------------------------------------------------------
#' #    Data to be fitted consists of two statistics s1 and s2
#' statistics_target <- data.frame(s1 = 9, s2 = 18)
#' #    We then define a parametrized model for the statistics
#' model <- function(parameters) {
#'     statistics <- data.frame(
#'         s1 = parameters$mu + parameters$theta + runif(nrow(parameters), -0.1, 0.1),
#'         s2 = parameters$mu * parameters$theta + runif(nrow(parameters), -0.1, 0.1)
#'     )
#'     cbind(parameters, statistics)
#' }
#' #    and the perturbation parameters
#' perturbation_parameters <- data.frame(
#'     theta = rep(0.1, 3),
#'     mu = rep(0.1, 3)
#' )
#' #    We define the rprior and dprior functions
#' rprior <- function(Nparameters) {
#'     theta <- runif(Nparameters, -10, 10)
#'     mu <- runif(Nparameters, -10, 10)
#'     return(data.frame(theta = theta, mu = mu))
#' }
#' dprior <- function(parameters, parameter_id = "all") {
#'     probs <- rep(1, nrow(parameters))
#'     if (parameter_id %in% c("all", "theta")) {
#'         probs <- probs * dunif(parameters[["theta"]], -10, 10)
#'     }
#'     if (parameter_id %in% c("all", "mu")) {
#'         probs <- probs * dunif(parameters[["mu"]], -10, 10)
#'     }
#'     return(probs)
#' }
#' #    Finally, we run ABC-SMC-DRF with 3 iterations, each with 1000 particles
#' smcrf_results <- smcrf(
#'     method = "smcrf-multi-param",
#'     statistics_target = statistics_target,
#'     model = model,
#'     rprior = rprior,
#'     dprior = dprior,
#'     perturbation = "Uniform",
#'     perturbation_parameters = perturbation_parameters,
#'     nParticles = c(1000, 1000, 1000),
#' )
#' #    Now we examine the posterior distribution of each parameter
#' posterior_iteration <- paste0("Iteration_", (smcrf_results$nIterations + 1))
#' posterior_params <- smcrf_results[[posterior_iteration]]$parameters
#' posterior_means <- colMeans(posterior_params)
#' posterior_vars <- var(posterior_params)
#' print(posterior_means)
#' print(posterior_vars)
smcrf <- function(method = "smcrf-single-param",
                  statistics_target = NULL,
                  statistics_selection = NULL,
                  smcrf_results = NULL,
                  model,
                  rprior,
                  dprior,
                  perturbation = "Gaussian",
                  perturbation_parameters = NULL,
                  nParticles,
                  final_sample = TRUE,
                  model_redo_if_NA = FALSE,
                  verbose = TRUE,
                  parallel = FALSE,
                  save_model = TRUE,
                  save_rds = FALSE,
                  filename_rds = "ABCSMCDRF.rds",
                  ...) {
    suppressPackageStartupMessages(library(matrixStats))
    suppressPackageStartupMessages(library(Hmisc))
    suppressPackageStartupMessages(library(crayon))
    if (perturbation == "Uniform") {
        if (is.null(perturbation_parameters)) {
            stop("Ranges must be specified for uniform perturbations")
        }
    } else if (perturbation != "Gaussian") {
        stop("Invalid perturbation method")
    }
    if (method == "smcrf-multi-param" & !is.null(statistics_selection)) {
        stop("statistics_selection is only available for method 'smcrf-single-param'")
    }
    if (method == "smcrf-single-param") {
        return(smcrf_single_param(
            statistics_target = statistics_target,
            statistics_selection = statistics_selection,
            model = model,
            rprior = rprior,
            dprior = dprior,
            perturbation = perturbation,
            perturbation_parameters = perturbation_parameters,
            nParticles = nParticles,
            final_sample = final_sample,
            model_redo_if_NA = model_redo_if_NA,
            verbose = verbose,
            parallel = parallel,
            save_model = save_model,
            save_rds = save_rds,
            filename_rds = filename_rds,
            smcrf_single_param_results = smcrf_results,
            ...
        ))
    } else if (method == "smcrf-multi-param") {
        if (!is.null(statistics_selection)) stop("statistics_selection is only available for method 'smcrf-single-param'")
        return(smcrf_multi_param(
            statistics_target = statistics_target,
            model = model,
            rprior = rprior,
            dprior = dprior,
            perturbation = perturbation,
            perturbation_parameters = perturbation_parameters,
            nParticles = nParticles,
            final_sample = final_sample,
            model_redo_if_NA = model_redo_if_NA,
            verbose = verbose,
            parallel = parallel,
            save_model = save_model,
            save_rds = save_rds,
            filename_rds = filename_rds,
            smcrf_multi_param_results = smcrf_results,
            ...
        ))
    } else {
        stop("Invalid method for smcrf")
    }
}

smcrf_single_param <- function(statistics_target,
                               statistics_selection = NULL,
                               model,
                               rprior,
                               dprior,
                               perturbation,
                               perturbation_parameters,
                               nParticles,
                               final_sample,
                               model_redo_if_NA,
                               verbose,
                               parallel,
                               save_model = TRUE,
                               save_rds = FALSE,
                               filename_rds,
                               smcrf_single_param_results,
                               ...) {
    suppressPackageStartupMessages(library(abcrf))
    nSimulations <<- 0
    #---Obtain information from previous SMC-RF
    if (!is.null(smcrf_single_param_results)) {
        #   ... If continuing from a previous ABC-SMC-DRF chain:
        SMCRF <- smcrf_single_param_results
        parameters_ids <- colnames(SMCRF[["Iteration_1"]]$parameters)
        statistics_ids <- colnames(SMCRF[["Iteration_1"]]$statistics)
        nIterations <- length(nParticles) + SMCRF[["nIterations"]]
        iteration_start <- 1 + SMCRF[["nIterations"]]
        if (is.null(statistics_target)) statistics_target <- SMCRF[["statistics_target"]]
        nParticles <- c(SMCRF[["nParticles"]], nParticles)
        SMCRF[[paste0("Iteration_", iteration_start)]] <- c()
        ABCRF_weights <- SMCRF[[paste0("Iteration_", iteration_start - 1)]]$weights
        parameters <- SMCRF[[paste0("Iteration_", iteration_start - 1)]]$parameters
    } else {
        #   ... If starting a new ABC-SMC-DRF chain:
        SMCRF <- list()
        parameters_ids <- colnames(rprior(Nparameters = 1))
        statistics_ids <- colnames(statistics_target)
        nIterations <- length(nParticles)
        iteration_start <- 1
    }
    iteration_end <- ifelse(final_sample, iteration_start + nIterations, iteration_start + nIterations - 1)
    SMCRF[["method"]] <- "smcrf-single-param"
    SMCRF[["nIterations"]] <- nIterations
    SMCRF[["nParticles"]] <- nParticles
    SMCRF[["statistics_target"]] <- statistics_target
    SMCRF[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
    SMCRF[["statistics_labels"]] <- data.frame(ID = colnames(statistics_target))
    if (!verbose) {
        pb <- txtProgressBar(
            min = 1,
            max = length(nParticles) + 1,
            style = 3,
            width = 50,
            char = "+"
        )
    }
    for (iteration in iteration_start:iteration_end) {
        if (verbose) {
            if (iteration == (nIterations + 1)) {
                cat(bold(red("ABC-SMC-DRF FOR SINGLE PARAMETERS:")), paste0(bold(yellow("final posterior distribution", "\n"))))
            } else {
                cat(bold(red("ABC-SMC-DRF FOR SINGLE PARAMETERS:")), paste0(bold(yellow("iteration", iteration, "\n"))))
            }
        } else {
            setTxtProgressBar(pb, iteration)
        }
        #---Compute Beaumont variances from the previous iteration
        if (iteration > 1) {
            tmp <- parameters[sample(nrow(parameters), size = 10000, prob = ABCRF_weights[, 1], replace = T), ]
            if (!is.data.frame(tmp)) {
                Beaumont_variances <- as.data.frame(max(var(tmp), 1e-10))
                colnames(Beaumont_variances) <- parameters_ids
            } else {
                Beaumont_variances <- pmax(sapply(tmp, var), 1e-10)
            }
        }
        #---Create training set
        if (verbose) cat(blue("Sampling parameters and computing model simulations...\n"))
        nrow <- ifelse(iteration == (nIterations + 1), nParticles[nIterations], nParticles[iteration])
        invalid_indices <- 1:nrow
        parameters_unperturbed <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
        colnames(parameters_unperturbed) <- parameters_ids
        parameters_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
        colnames(parameters_next) <- parameters_ids
        reference_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids) + ncol(statistics_target)))
        colnames(reference_next) <- colnames(model(parameters = rprior(Nparameters = 1)))
        while (length(invalid_indices) > 0) {
            #   Generate new particles...
            if (iteration == 1) {
                #   ... for iteration 1:
                #   Sample from the prior distribution
                parameters_unperturbed[invalid_indices, ] <- rprior(Nparameters = length(invalid_indices))
                parameters_next[invalid_indices, ] <- parameters_unperturbed[invalid_indices, ]
            } else {
                #   ... for later iterations:
                #   Sample from the previous posterior distribution
                parameters_previous <- parameters
                weights_previous <- ABCRF_weights
                parameter_replace <- data.frame(matrix(NA, nrow = length(invalid_indices), ncol = length(parameters_ids)))
                colnames(parameter_replace) <- parameters_ids
                for (parameter_id in parameters_ids) {
                    parameter_replace[, parameter_id] <- parameters_previous[sample(nrow(parameters_previous), size = length(invalid_indices), prob = weights_previous[, parameter_id], replace = T), parameter_id]
                }
                parameters_unperturbed[invalid_indices, ] <- parameter_replace
                #   Perturb parameters
                if (iteration < (nIterations + 1)) {
                    if (perturbation == "Gaussian") {
                        if (is.null(perturbation_parameters)) {
                            Gaussian_variances <- Beaumont_variances
                        } else {
                            Gaussian_variances <- perturbation_parameters[iteration - 1, ]
                        }
                        for (parameter_id in parameters_ids) {
                            parameter_replace[[parameter_id]] <- rnorm(
                                n = nrow(parameter_replace),
                                mean = parameter_replace[[parameter_id]],
                                sd = sqrt(Gaussian_variances[[parameter_id]])
                            )
                        }
                    } else if (perturbation == "Uniform") {
                        Uniform_ranges <- perturbation_parameters[iteration - 1, ]
                        for (parameter_id in parameters_ids) {
                            parameter_replace[[parameter_id]] <- runif(
                                n = nrow(parameter_replace),
                                min = parameter_replace[[parameter_id]] - Uniform_ranges[[parameter_id]],
                                max = parameter_replace[[parameter_id]] + Uniform_ranges[[parameter_id]]
                            )
                        }
                    }
                }
                parameters_next[invalid_indices, ] <- parameter_replace
            }
            invalid_indices_next <- which(dprior(parameters_next, parameter_id = "all") <= 0)
            #   Generate statistics for the new particles
            ids <- setdiff(invalid_indices, invalid_indices_next)
            tmp <- parameters_next[ids, ]
            if (!is.data.frame(tmp)) {
                tmp <- as.data.frame(tmp, ncol = length(parameters_ids))
                colnames(tmp) <- parameters_ids
            }
            if (length(ids) > 0) reference_next[ids, ] <- model(parameters = tmp)
            #   Find particles that need to be regenerated
            invalid_indices <- invalid_indices_next
            if (model_redo_if_NA) invalid_indices <- unique(c(invalid_indices, which(apply(reference_next, 1, function(x) any(is.na(x))))))
        }
        reference <- reference_next
        parameters <- reference[, parameters_ids]
        if (!is.data.frame(parameters)) {
            parameters <- as.data.frame(parameters, ncol = length(parameters_ids))
            colnames(parameters) <- parameters_ids
        }
        statistics <- reference[, statistics_ids]
        if (!is.data.frame(statistics)) {
            statistics <- as.data.frame(statistics, ncol = length(statistics_ids))
            colnames(statistics) <- statistics_ids
        }
        #   Finish the last iteration
        if (iteration == (nIterations + 1)) {
            SMCRF_iteration <- list()
            SMCRF_iteration$reference <- reference
            SMCRF_iteration$parameters <- parameters
            SMCRF_iteration$parameters_unperturbed <- parameters_unperturbed
            SMCRF_iteration$statistics <- statistics
            SMCRF[[paste0("Iteration_", iteration)]] <- SMCRF_iteration
            if (save_rds == TRUE) {
                saveRDS(SMCRF, file = filename_rds)
            }
            if (!verbose) cat("\n")
            return(SMCRF)
        }
        #---Run ABCRF for all parameters
        if (verbose) cat(blue("Performing Random Forest prediction...\n"))
        ABCRF_weights <- data.frame(matrix(NA, nrow = nParticles[iteration], ncol = 0))
        if (save_model == TRUE) {
            RFmodels <- list()
            posterior_gamma_RFs <- list()
        }
        for (parameter_id in parameters_ids) {
            #   Select statistics for the parameter
            if (is.null(statistics_selection)) {
                mini_reference <- reference[, c(parameter_id, colnames(reference)[!colnames(reference) %in% parameters_ids])]
            } else {
                mini_reference <- reference[, c(parameter_id, colnames(statistics_selection)[which(statistics_selection[rownames(statistics_selection) == parameter_id, ] == 1)])]
            }
            colnames(mini_reference)[1] <- "para"
            RFmodel <- regAbcrf(
                formula = as.formula("para ~."),
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
            if (save_model == TRUE) {
                RFmodels[[parameter_id]] <- RFmodel
                posterior_gamma_RFs[[parameter_id]] <- posterior_gamma_RF
            }
        }
        #---Modify ABCRF weights
        if (verbose) cat(blue("Recalibrating Random Forest weights...\n"))
        if (iteration > 1) {
            for (parameter_id in parameters_ids) {
                if (perturbation == "Gaussian") {
                    weight_modifiers <- w_modifiers_normal(
                        parameters = parameters,
                        parameters_previous = parameters_previous,
                        weights_previous = weights_previous[, parameter_id],
                        perturb_variances = Gaussian_variances,
                        dprior = dprior,
                        parameter_id = parameter_id
                    )
                } else if (perturbation == "Uniform") {
                    weight_modifiers <- w_modifiers_uniform(
                        parameters = parameters,
                        parameters_previous = parameters_previous,
                        weights_previous = weights_previous[, parameter_id],
                        perturb_ranges = Uniform_ranges,
                        dprior = dprior,
                        parameter_id = parameter_id
                    )
                    weight_modifiers <- rep(1, length(weight_modifiers))
                }
                ABCRF_weights[, parameter_id] <- as.vector(ABCRF_weights[, parameter_id]) * weight_modifiers
                ABCRF_weights[, parameter_id] <- ABCRF_weights[, parameter_id] / sum(ABCRF_weights[, parameter_id])
            }
        }
        #---Save SMC-RF results from this iteration
        SMCRF_iteration <- list()
        SMCRF_iteration$reference <- reference
        SMCRF_iteration$parameters <- parameters
        SMCRF_iteration$parameters_unperturbed <- parameters_unperturbed
        SMCRF_iteration$statistics <- statistics
        SMCRF_iteration$weights <- ABCRF_weights
        if (save_model == TRUE) {
            SMCRF_iteration$rf_model <- RFmodels
            SMCRF_iteration$rf_predict <- posterior_gamma_RFs
        }
        SMCRF[[paste0("Iteration_", iteration)]] <- SMCRF_iteration
        if (save_rds == TRUE) {
            saveRDS(SMCRF, file = filename_rds)
        }
    }
    if (save_rds == TRUE) {
        saveRDS(SMCRF, file = filename_rds)
    }
    if (!verbose) cat("\n")
    return(SMCRF)
}

smcrf_multi_param <- function(statistics_target,
                              model,
                              rprior,
                              dprior,
                              perturbation,
                              perturbation_parameters,
                              nParticles,
                              final_sample,
                              model_redo_if_NA,
                              verbose,
                              parallel,
                              save_model,
                              save_rds,
                              filename_rds,
                              splitting.rule = "CART",
                              smcrf_multi_param_results,
                              ...) {
    suppressPackageStartupMessages(library(drf))
    nSimulations <<- 0
    #---Obtain information from previous SMC-DRF
    if (!is.null(smcrf_multi_param_results)) {
        #   ... If continuing from a previous ABC-SMC-DRF chain:
        SMCDRF <- smcrf_multi_param_results
        parameters_ids <- colnames(SMCDRF[["Iteration_1"]]$parameters)
        statistics_ids <- colnames(SMCDRF[["Iteration_1"]]$statistics)
        nIterations <- length(nParticles) + SMCDRF[["nIterations"]]
        iteration_start <- 1 + SMCDRF[["nIterations"]]
        if (is.null(statistics_target)) statistics_target <- SMCDRF[["statistics_target"]]
        nParticles <- c(SMCDRF[["nParticles"]], nParticles)
        SMCDRF[[paste0("Iteration_", iteration_start)]] <- c()
        DRF_weights <- SMCDRF[[paste0("Iteration_", iteration_start - 1)]]$weights
        parameters <- SMCDRF[[paste0("Iteration_", iteration_start - 1)]]$parameters
    } else {
        #   ... If starting a new ABC-SMC-DRF chain:
        SMCDRF <- list()
        parameters_ids <- colnames(rprior(Nparameters = 1))
        statistics_ids <- colnames(statistics_target)
        nIterations <- length(nParticles)
        iteration_start <- 1
    }
    iteration_end <- ifelse(final_sample, iteration_start + nIterations, iteration_start + nIterations - 1)
    SMCDRF[["method"]] <- "smcrf-multi-param"
    SMCDRF[["nIterations"]] <- nIterations
    SMCDRF[["nParticles"]] <- nParticles
    SMCDRF[["statistics_target"]] <- statistics_target
    SMCDRF[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
    SMCDRF[["statistics_labels"]] <- data.frame(ID = statistics_ids)
    if (!verbose) {
        pb <- txtProgressBar(
            min = 1,
            max = length(nParticles) + 1,
            style = 3,
            width = 50,
            char = "+"
        )
    }
    for (iteration in iteration_start:iteration_end) {
        if (verbose) {
            if (iteration == (nIterations + 1)) {
                cat(bold(red("ABC-SMC-DRF FOR MULTIPLE PARAMETERS:")), paste0(bold(yellow("final posterior distribution", "\n"))))
            } else {
                cat(bold(red("ABC-SMC-DRF FOR MULTIPLE PARAMETERS:")), paste0(bold(yellow("iteration", iteration, "\n"))))
            }
        } else {
            setTxtProgressBar(pb, iteration)
        }
        #---Compute Beaumont variances from the previous iteration
        if (iteration > 1) {
            tmp <- parameters[sample(nrow(parameters), size = 10000, prob = DRF_weights[, 1], replace = T), ]
            if (!is.data.frame(tmp)) {
                Beaumont_variances <- as.data.frame(max(var(tmp), 1e-10))
                colnames(Beaumont_variances) <- parameters_ids
            } else {
                Beaumont_variances <- pmax(sapply(tmp, var), 1e-10)
            }
        }
        #---Create training set
        if (verbose) cat(blue("Sampling parameters and computing model simulations...\n"))
        nrow <- ifelse(iteration == (nIterations + 1), nParticles[nIterations], nParticles[iteration])
        invalid_indices <- 1:nrow
        parameters_unperturbed <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
        colnames(parameters_unperturbed) <- parameters_ids
        parameters_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids)))
        colnames(parameters_next) <- parameters_ids
        reference_next <- data.frame(matrix(NA, nrow = nrow, ncol = length(parameters_ids) + ncol(statistics_target)))
        colnames(reference_next) <- colnames(model(parameters = rprior(Nparameters = 1)))
        while (length(invalid_indices) > 0) {
            #   Generate new particles...
            if (iteration == 1) {
                #   ... for iteration 1:
                #   Sample from the prior distribution
                parameters_unperturbed[invalid_indices, ] <- rprior(Nparameters = length(invalid_indices))
                parameters_next[invalid_indices, ] <- parameters_unperturbed[invalid_indices, ]
            } else {
                #   ... for later iterations:
                #   Sample from the previous posterior distribution
                parameters_previous <- parameters
                weights_previous <- DRF_weights
                parameter_replace <- parameters[sample(nrow(parameters), size = length(invalid_indices), prob = DRF_weights[, 1], replace = T), ]
                parameters_unperturbed[invalid_indices, ] <- parameter_replace
                #   Perturb parameters
                if (iteration < (nIterations + 1)) {
                    if (perturbation == "Gaussian") {
                        if (is.null(perturbation_parameters)) {
                            Gaussian_variances <- Beaumont_variances
                        } else {
                            Gaussian_variances <- perturbation_parameters[iteration - 1, ]
                        }
                        for (parameter_id in parameters_ids) {
                            parameter_replace[[parameter_id]] <- rnorm(
                                n = nrow(parameter_replace),
                                mean = parameter_replace[[parameter_id]],
                                sd = sqrt(Gaussian_variances[[parameter_id]])
                            )
                        }
                    } else if (perturbation == "Uniform") {
                        Uniform_ranges <- perturbation_parameters[iteration - 1, ]
                        for (parameter_id in parameters_ids) {
                            parameter_replace[[parameter_id]] <- runif(
                                n = nrow(parameter_replace),
                                min = parameter_replace[[parameter_id]] - Uniform_ranges[[parameter_id]],
                                max = parameter_replace[[parameter_id]] + Uniform_ranges[[parameter_id]]
                            )
                        }
                    }
                }
                parameters_next[invalid_indices, ] <- parameter_replace
            }
            invalid_indices_next <- which(dprior(parameters_next, parameter_id = "all") <= 0)
            #   Generate statistics for the new particles
            ids <- setdiff(invalid_indices, invalid_indices_next)
            tmp <- parameters_next[ids, ]
            if (!is.data.frame(tmp)) {
                tmp <- as.data.frame(tmp, ncol = length(parameters_ids))
                colnames(tmp) <- parameters_ids
            }
            if (length(ids) > 0) reference_next[ids, ] <- model(parameters = tmp)
            #   Find particles that need to be regenerated
            invalid_indices <- invalid_indices_next
            if (model_redo_if_NA) invalid_indices <- unique(c(invalid_indices, which(apply(reference_next, 1, function(x) any(is.na(x))))))
        }
        reference <- reference_next
        parameters <- reference[, parameters_ids]
        if (!is.data.frame(parameters)) {
            parameters <- as.data.frame(parameters, ncol = length(parameters_ids))
            colnames(parameters) <- parameters_ids
        }
        statistics <- reference[, statistics_ids]
        if (!is.data.frame(statistics)) {
            statistics <- as.data.frame(statistics, ncol = length(statistics_ids))
            colnames(statistics) <- statistics_ids
        }
        #   Finish the last iteration
        if (iteration == (nIterations + 1)) {
            SMCDRF_iteration <- list()
            SMCDRF_iteration$reference <- reference
            SMCDRF_iteration$parameters <- parameters
            SMCDRF_iteration$parameters_unperturbed <- parameters_unperturbed
            SMCDRF_iteration$statistics <- statistics
            SMCDRF[[paste0("Iteration_", iteration)]] <- SMCDRF_iteration
            if (save_rds == TRUE) {
                saveRDS(SMCDRF, file = filename_rds)
            }
            if (!verbose) cat("\n")
            return(SMCDRF)
        }
        #---Run DRF for all parameters
        if (verbose) cat(blue("Performing Random Forest prediction...\n"))
        drfmodel <- drf(
            X = statistics,
            Y = parameters,
            splitting.rule = splitting.rule,
            ...
        )
        def_pred <- predict(
            object = drfmodel,
            newdata = statistics_target
        )
        DRF_weights <- as.vector(get_sample_weights(drfmodel, statistics_target))
        #---Modify DRF weights
        if (verbose) cat(blue("Recalibrating Random Forest weights...\n"))
        if (iteration > 1) {
            if (perturbation == "Gaussian") {
                weight_modifiers <- w_modifiers_normal(
                    parameters = parameters,
                    parameters_previous = parameters_previous,
                    weights_previous = weights_previous[, 1],
                    perturb_variances = Gaussian_variances,
                    dprior = dprior,
                    parameter_id = "all"
                )
            } else if (perturbation == "Uniform") {
                weight_modifiers <- w_modifiers_uniform(
                    parameters = parameters,
                    parameters_previous = parameters_previous,
                    weights_previous = weights_previous[, 1],
                    perturb_ranges = Uniform_ranges,
                    dprior = dprior,
                    parameter_id = "all"
                )
                weight_modifiers <- rep(1, length(weight_modifiers))
            }
            DRF_weights <- DRF_weights * weight_modifiers
            DRF_weights <- DRF_weights / sum(DRF_weights)
        }
        DRF_weights <- data.frame(matrix(rep(DRF_weights, length(parameters_ids)), ncol = length(parameters_ids)))
        colnames(DRF_weights) <- parameters_ids
        #---Save SMC-DRF results from this iteration
        SMCDRF_iteration <- list()
        SMCDRF_iteration$reference <- reference
        SMCDRF_iteration$parameters <- parameters
        SMCDRF_iteration$parameters_unperturbed <- parameters_unperturbed
        SMCDRF_iteration$statistics <- statistics
        SMCDRF_iteration$weights <- DRF_weights
        if (save_model == TRUE) {
            SMCDRF_iteration$rf_model <- drfmodel
            SMCDRF_iteration$rf_predict <- def_pred
        }
        SMCDRF[[paste0("Iteration_", iteration)]] <- SMCDRF_iteration
        if (save_rds == TRUE) {
            saveRDS(SMCDRF, file = filename_rds)
        }
    }
    if (save_rds == TRUE) {
        saveRDS(SMCDRF, file = filename_rds)
    }
    if (!verbose) cat("\n")
    return(SMCDRF)
}

w_modifiers_normal <- function(parameters,
                               parameters_previous,
                               weights_previous,
                               perturb_variances,
                               dprior,
                               parameter_id = "all") {
    n_parameters_previous <- dim(parameters_previous)[1]
    n_parameters <- dim(parameters)[1]
    if (parameter_id == "all") {
        parameter_ids <- colnames(parameters)
    } else {
        parameter_ids <- parameter_id
    }
    #---Compute denominators for weight recalibration
    weights <- array(0, n_parameters)
    for (i in 1:n_parameters_previous) {
        tab_temp <- array(weights_previous[i], n_parameters)
        for (j in parameter_ids) {
            tab_temp <- tab_temp *
                exp(
                    -0.5 * (parameters[[j]] - parameters_previous[[j]][i]) *
                        (parameters[[j]] - parameters_previous[[j]][i]) /
                        perturb_variances[[j]]
                )
        }
        weights <- weights + tab_temp
    }
    #---Compute numerators for weight recalibration
    weights_previous_prior <- dprior(parameters, parameter_id = parameter_id)
    #---Compute weights for new particles
    weights <- weights_previous_prior / weights
    weights <- weights / sum(weights)
    return(weights)
}

w_modifiers_uniform <- function(parameters,
                                parameters_previous,
                                weights_previous,
                                perturb_ranges,
                                dprior,
                                parameter_id = "all") {
    n_parameters_previous <- dim(parameters_previous)[1]
    n_parameters <- dim(parameters)[1]
    if (parameter_id == "all") {
        parameter_ids <- colnames(parameters)
    } else {
        parameter_ids <- parameter_id
    }
    #---Compute denominators for weight recalibration
    weights <- array(0, n_parameters)
    for (i in 1:n_parameters_previous) {
        tab_temp <- array(weights_previous[i], n_parameters)
        for (j in parameter_ids) {
            tab_temp[which(parameters[[j]] < parameters_previous[[j]][i] - perturb_ranges[j] |
                parameters[[j]] > parameters_previous[[j]][i] + perturb_ranges[j])] <- 0
        }
        weights <- weights + tab_temp
    }
    #---Compute numerators for weight recalibration
    weights_previous_prior <- dprior(parameters, parameter_id = parameter_id)
    #---Compute weights for new particles
    weights <- weights_previous_prior / weights
    weights <- weights / sum(weights)
    return(weights)
}
