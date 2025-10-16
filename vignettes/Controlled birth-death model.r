library(abcsmcrf)
library(LaplacesDemon)
library(truncnorm)
set.seed(1)
# ===================================Function to compute estimation of m
estimated_mean <- function(abc_results) {
    nIterations <- abc_results[["nIterations"]]
    ref_df <- abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed
    nb_rows <- nrow(ref_df)
    m <- numeric(nb_rows)
    for (i in seq_len(nb_rows)) {
        kappa_i <- ref_df$kappa[i]
        js <- 0:kappa_i
        prob_js <- unlist(ref_df[i, paste0("prob_", js)])
        m[i] <- sum(js * prob_js)
    }
    abc_results[[paste0("Iteration_", nIterations + 1)]]$parameters_unperturbed$m <- m
    return(abc_results)
}
# =====================================Model for the birth-death process
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters, parallel = TRUE) {
    nTimes <- 10
    K_max <- 7
    if (exists("nSimulations")) nSimulations <<- nSimulations + nrow(parameters)
    #   Make simulations & compute summary statistics (population sizes at each time point)
    if (parallel) {
        library(parallel)
        library(pbapply)
        library(data.table)
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cl, varlist = c("BD_model"))
        stats <- pblapply(
            cl = cl, X = 1:nrow(parameters),
            FUN = function(i) BD_model(parameters$kappa[i], as.numeric(parameters[i, c(paste0("prob_", 0:K_max))]), parameters$gamma[i], nTimes)
        )
        stopCluster(cl)
        stats <- rbindlist(stats)
        class(stats) <- "data.frame"
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            stats <- rbind(stats, BD_model(parameters$kappa[i], as.numeric(parameters[i, c(paste0("prob_", 0:K_max))]), parameters$gamma[i], nTimes))
        }
    }
    data <- stats
    return(data)
}
BD_model <- function(kappa, prob_X, gamma, nTimes) {
    #   This simulates a controlled branching process
    #   Input:  lambda = birth rate
    #           mu = death rate
    kappa <- floor(kappa)
    phi_sim <- rep(0, nTimes)
    Zsim <- rep(0, nTimes + 1)
    #   Now generate observations on phi and Z
    zv <- 1
    Zsim[1] <- zv
    for (j in seq(2, length(Zsim))) {
        if (Zsim[j - 1] == 0) {
            size_sim <- 0
        } else {
            size_sim <- Zsim[j - 1] + floor(log(Zsim[j - 1]))
        }
        phi_sim[j - 1] <- rbinom(1, size = size_sim, prob = gamma)
        if (phi_sim[j - 1] == 0) {
            Zsim[j] <- 0
        } else {
            if (phi_sim[j - 1] < .Machine$integer.max) {
                prob_X_new <- rep(0, length(prob_X))
                prob_X_new[1:(kappa + 1)] <- prob_X[1:(kappa + 1)]
                X <- sample(0:(length(prob_X) - 1), size = phi_sim[j - 1], replace = TRUE, prob = prob_X_new)
                Zsim[j] <- sum(X)
            } else {
                Zsim[j] <- 0
            }
        }
    }
    if (Zsim[length(Zsim)] == 0) {
        stats <- data.frame(matrix(rep(NA, (length(Zsim) + 3 + length(prob_X))), nrow = 1, ncol = (length(Zsim) + 3 + length(prob_X))))
    } else {
        stats <- data.frame(matrix(c(kappa, prob_X[seq_along(prob_X)], gamma, Zsim[1:length(Zsim)], phi_sim[length(Zsim) - 1]), nrow = 1))
    }
    colnames(stats) <- c("kappa", paste0("prob_", 0:(length(prob_X) - 1)), "gamma", paste0("Z_", 0:(length(Zsim) - 1)), paste0("Phi_", length(Zsim) - 2))
    return(stats)
}
# =====================================================Target statistics
statistics_target <- data.frame(matrix(c(
    1, 4, 12, 30, 84, 249, 728, 2148, 6165, 17883, 51412, 14281
), nrow = 1))
colnames(statistics_target) <- c(paste0("Z_", 0:10), "Phi_9")
K_max <- 7
rho <- 0.9
kappa <- 4
gamma <- 0.8
prob_X <- rep(0, K_max + 1)
x <- 0:kappa
prob_X[1:(kappa + 1)] <- dbinom(x, size = kappa, prob = rho)
parameters_target <- data.frame(
    m = 3.6,
    gamma = 0.8,
    kappa = 4
)
# ====================================================Prior distribution
dprior <- function(parameters, parameter_id = "all") {
    probs <- rep(1, nrow(parameters))
    if (parameter_id %in% c("all", "gamma")) {
        probs <- probs * dbeta(parameters[["gamma"]], shape1 = 1, shape2 = 1)
    }
    if (parameter_id %in% c("all", "kappa")) {
        invalid_kappa <- (parameters[["kappa"]] <= 0) | (parameters[["kappa"]] > (K_max + 1))
        probs[invalid_kappa] <- 0
        valid_idx <- !invalid_kappa
        probs[valid_idx] <- probs[valid_idx] * 1 / (K_max - 1)
    }
    if (parameter_id %in% c("all", paste0("prob_", 0:K_max))) {
        for (i in 1:nrow(parameters)) {
            kappa_i <- floor(parameters$kappa[i])
            condition_kappa <- (kappa_i >= 0) && (kappa_i < K_max + 1)
            if (condition_kappa) {
                x <- as.numeric(parameters[i, paste0("prob_", 0:kappa_i)])
                if (all(x >= 0)) {
                    probs[i] <- probs[i] * ddirichlet(x, alpha = rep(1, kappa_i + 1))
                } else {
                    probs[i] <- 0
                }
            } else {
                probs[i] <- 0
            }
        }
    }
    return(probs)
}
rprior <- function(Nparameters) {
    gamma <- rbeta(Nparameters, shape1 = 1, shape2 = 1)
    kappa <- sample(2:K_max, size = Nparameters, replace = TRUE)
    probabilities <- matrix(0, nrow = Nparameters, ncol = K_max + 1)
    for (i in 1:Nparameters) {
        prior_sim <- rdirichlet(1, alpha = rep(1, kappa[i] + 1))
        probabilities[i, 1:(kappa[i] + 1)] <- prior_sim
    }
    probabilities <- as.data.frame(probabilities)
    colnames(probabilities) <- paste0("prob_", 0:K_max)
    parameter_df <- data.frame(gamma = gamma, kappa = kappa)
    return(cbind(parameter_df, probabilities))
}
# =============================================Perturbation distribution
dperturb <- function(parameters, parameters_previous, parameters_previous_sampled, iteration, parameter_id = "all") {
    probs <- rep(1, nrow(parameters))
    if (parameter_id %in% c("all", "gamma")) {
        probs <- probs * dunif(
            parameters[["gamma"]],
            min = pmax(0, parameters_previous[["gamma"]] - 0.05),
            max = pmin(1, parameters_previous[["gamma"]] + 0.05)
        )
    }
    if (parameter_id %in% c("all", "kappa")) {
        probs <- probs * dunif(
            parameters[["kappa"]],
            min = pmax(0, parameters_previous[["kappa"]] - 3),
            max = pmin(K_max + 1, parameters_previous[["kappa"]] + 3)
        )
    }
    for (i in 0:K_max) {
        if (parameter_id %in% c("all", paste0("prob_", i))) {
            probs <- probs * dunif(
                parameters[[paste0("prob_", i)]],
                min = pmax(0, parameters_previous[[paste0("prob_", i)]] - 0.05),
                max = pmin(1, parameters_previous[[paste0("prob_", i)]] + 0.05)
            )
        }
    }
    return(probs)
}
rperturb <- function(parameters_unperturbed, parameters_previous_sampled, iteration) {
    parameters_perturbed <- parameters_unperturbed
    parameters_perturbed[["gamma"]] <- runif(
        n = nrow(parameters_perturbed),
        min = pmax(0, parameters_perturbed[["gamma"]] - 0.05),
        max = pmin(1, parameters_perturbed[["gamma"]] + 0.05)
    )
    parameters_perturbed[["kappa"]] <- runif(
        n = nrow(parameters_perturbed),
        min = pmax(0, parameters_perturbed[["kappa"]] - 3),
        max = pmin(K_max + 1, parameters_perturbed[["kappa"]] + 3)
    )
    for (i in 0:K_max) {
        parameters_perturbed[[paste0("prob_", i)]] <- runif(
            n = nrow(parameters_perturbed),
            min = pmax(0, parameters_perturbed[[paste0("prob_", i)]] - 0.05),
            max = pmin(K_max + 1, parameters_perturbed[[paste0("prob_", i)]] + 0.05)
        )
    }
    return(parameters_perturbed)
}
# ====================================Labels for parameters in the plots
parameters_labels <- data.frame(
    parameter = c("m", "gamma", "kappa", paste0("prob_", 0:K_max)),
    label = c(deparse(expression(m)), deparse(expression(gamma)), deparse(expression(kappa)), paste0("prob_", 0:K_max))
)
# ===================================================================DRF
#---Run DRF
drf_results <- smcrf(
    method = "smcrf-multi-param",
    statistics_target = statistics_target,
    model = model,
    rprior = rprior,
    dprior = dprior,
    nParticles = rep(20000, 1),
    num.trees = 2500,
    save_model = FALSE,
    splitting.rule = "FourierMMD",
    parallel = TRUE,
    model_redo_if_NA = TRUE,
    honesty = FALSE
)
#---Compute m
drf_results <- estimated_mean(drf_results)
#---Plot posterior joint distributions against other methods
plots_joint <- plot_compare_joint(
    abc_results = drf_results,
    parameters_labels = data.frame(
        parameter = c("m", "gamma"),
        label = c(deparse(expression(m)), deparse(expression(gamma)))
    ),
    nBins = 9
)
#---Plot qqplots against other methods
plots_compare_qqplot <- plot_compare_qqplot(
    abc_results = drf_results,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels,
)
#---Plot marginal distributions against other methods
plots_marginal <- plot_compare_marginal(
    abc_results = drf_results,
    parameters_labels = parameters_labels,
    parameters_truth = parameters_target,
    plot_hist = TRUE
)
# ====================================================================RF
#---Run ABC-RF
abcrf_results <- smcrf(
    method = "smcrf-single-param",
    statistics_target = statistics_target,
    model = model,
    rprior = rprior,
    dprior = dprior,
    ntree = 2500,
    nParticles = rep(20000, 1),
    save_model = FALSE,
    parallel = TRUE,
    model_redo_if_NA = TRUE
)
#---Compute m
abcrf_results <- estimated_mean(abcrf_results)
#---Plot qqplots against other methods
plots_compare_qqplot <- plot_compare_qqplot(
    plots = plots_compare_qqplot,
    abc_results = abcrf_results,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels
)
#---Plot marginal distributions against other methods
plots_marginal <- plot_compare_marginal(
    plots = plots_marginal,
    abc_results = abcrf_results,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)
# ================================================================SMC-RF
#---Run ABC-SMC-RF
abcrf_results <- smcrf(
    method = "smcrf-single-param",
    statistics_target = statistics_target,
    model = model,
    rprior = rprior,
    dprior = dprior,
    rperturb = rperturb,
    dperturb = dperturb,
    ntree = 2500,
    nParticles = rep(5000, 4),
    save_model = FALSE,
    parallel = TRUE,
    model_redo_if_NA = TRUE
)
#---Compute m
abcrf_results <- estimated_mean(abcrf_results)
#---Plot marginal distributions against other methods
plots_marginal <- plot_compare_marginal(
    abc_results = abcrf_results,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)
#---Plot qqplots against other methods
plots_compare_qqplot <- plot_compare_qqplot(
    plots = plots_compare_qqplot,
    abc_results = abcrf_results,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels
)
# ===============================================================SMC-DRF
#---Run SMC-DRF
smcrf_results_multi_param <- smcrf(
    method = "smcrf-multi-param",
    statistics_target = statistics_target,
    model = model,
    rprior = rprior,
    dprior = dprior,
    rperturb = rperturb,
    dperturb = dperturb,
    nParticles = rep(5000, 4),
    save_model = FALSE,
    num.trees = 2500,
    splitting.rule = "FourierMMD",
    model_redo_if_NA = TRUE,
    parallel = TRUE
)
#---Compute m
smcrf_results_multi_param <- estimated_mean(smcrf_results_multi_param)
#---Plot posterior joint distributions against other methods
plots_joint <- plot_compare_joint(
    plots = plots_joint,
    abc_results = smcrf_results_multi_param,
    parameters_labels = data.frame(
        parameter = c("m", "gamma"),
        label = c(deparse(expression(m)), deparse(expression(gamma)))
    )
)
#---Plot qqplots
plots_compare_qqplot <- plot_compare_qqplot(
    plots = plots_compare_qqplot,
    abc_results = smcrf_results_multi_param,
    lims = data.frame(
        parameter = c("kappa", paste0("prob_", 0:(length(prob_X) - 1)), "gamma"),
        min = c(0, rep(0, length(prob_X)), 0),
        max = c(10, rep(1, length(prob_X)), 1)
    ),
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels
)
#---Plot marginal distributions against other methods
plots_marginal <- plot_compare_marginal(
    plots = plots_marginal,
    abc_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)
