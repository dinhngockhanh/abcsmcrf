library(abcsmcrf)
# =====================================Model for the predator–prey process
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
# nSimulations <- 0
model <- function(parameters, parallel = TRUE) {
    nNoise <- 0
    if (exists("nSimulations")) {
        nSimulations <<- nSimulations + nrow(parameters)
    }
    # if (exists("nSimulations")) nSimulations <<- nSimulations + nrow(parameters)
    #   Make simulations & compute summary statistics (population sizes at each time point)
    if (parallel) {
        library(parallel)
        library(pbapply)
        library(data.table)
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cl, varlist = c("LV_model"))
        stats <- parLapply(
            cl = cl, 1:nrow(parameters),
            function(i) LV_model(parameters$a[i], parameters$b[i])
        )
        stopCluster(cl)
        stats <- rbindlist(stats)
        class(stats) <- "data.frame"
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            stats <- rbind(stats, LV_model(parameters$a[i], parameters$b[i]))
        }
    }
    #   Add noise statistics
    noise <- matrix(runif(nrow(parameters) * nNoise), nrow(parameters), nNoise)
    data <- cbind(stats, noise)
    #   Add column names
    data <- data.frame(data)
    if (nNoise > 0) {
        colnames(data) <- c(
            colnames(stats),
            paste0("noise_", c(1:nNoise))
        )
    }
    return(data)
}
LV_model <- function(a, b, parallel = FALSE) {
    nRuns <- 1
    library(deSolve)
    LVmatrix <- function(Time, State, parms) {
        with(as.list(c(State, parms)), {
            dx <- x * (alpha - beta * y)
            dy <- -y * (gamma - delta * x)
            return(list(c(dx, dy)))
        })
    }
    parms <- c(alpha = a, beta = 1, gamma = 1, delta = b)
    State <- c(x = 1, y = 0.5)
    Time <- seq(0, 15, by = 0.01)
    Xsim <- numeric(length = 8)
    Ysim <- numeric(length = 8)
    X <- NA
    Y <- NA
    for (i in 1:nRuns) {
        out <- as.data.frame(ode(func = LVmatrix, y = State, parms = parms, times = Time))
        tsim <- c(1.1, 2.4, 3.9, 5.6, 7.5, 9.6, 11.9, 14.4)
        index <- c(111, 241, 391, 561, 751, 961, 1191, 1441)
        X <- out$x[index]
        Y <- out$y[index]
        if (anyNA(X)) {
            Xend <- as.double(tail(out, n = 1)["x"])
            X[is.na(X)] <- Xend
            Yend <- as.double(tail(out, n = 1)["y"])
            Y[is.na(Y)] <- Yend
        }
        Xsim <- Xsim + X / nRuns
        Ysim <- Ysim + Y / nRuns
    }
    stats <- data.frame(matrix(c(a, b, Xsim, Ysim), nrow = 1))
    colnames(stats) <- c("a", "b", paste0("X_", 1:8), paste0("Y_", 1:8))
    return(stats)
}
# =====================================================Target statistics
parameters_truth <- data.frame(
    a = 1,
    b = 1
)
Xt <- c(1.87, 0.65, 0.22, 0.31, 1.64, 1.15, 0.24, 2.91)
Yt <- c(0.49, 2.62, 1.54, 0.02, 1.14, 1.68, 1.07, 0.88)
statistics_target <- data.frame(matrix(c(Xt, Yt), nrow = 1))
colnames(statistics_target) <- c(paste0("X_", 1:8), paste0("Y_", 1:8))
# ====================================================Prior distribution
dprior <- function(parameters, parameter_id = "all") {
    probs <- rep(1, nrow(parameters))
    if (parameter_id %in% c("all", "a")) {
        probs <- probs * dunif(parameters[["a"]], min = -10, max = 10)
    }
    if (parameter_id %in% c("all", "b")) {
        probs <- probs * dunif(parameters[["b"]], min = -10, max = 10)
    }
    return(probs)
}
rprior <- function(Nparameters) {
    parameters <- data.frame(
        a = runif(Nparameters, min = -10, max = 10),
        b = runif(Nparameters, min = -10, max = 10)
    )
    return(parameters)
}
# ====================================Labels for parameters in the plots
parameters_labels <- data.frame(
    parameter = c("a", "b"),
    label = c(deparse(expression(a)), deparse(expression(b)))
)
# ========================================SMC-RF for multiple parameters
#---Run SMC-RF for multiple parameters
smcrf_results_multi_param <- smcrf(
    method = "smcrf-multi-param",
    statistics_target = statistics_target,
    model = model,
    rprior = rprior,
    dprior = dprior,
    nParticles = rep(5000, 4),
    perturbation = "Uniform",
    perturbation_parameters = data.frame(a = rep(0.1, 3), b = rep(0.1, 3)),
    parallel = TRUE
)
#---Plot marginal distributions
plot_smcrf_marginal(
    smcrf_results = smcrf_results_multi_param,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)
#---Plot posterior marginal distributions against other methods
plots <- plot_compare_marginal(
    abc_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels,
    parameters_truth = parameters_truth,
    plot_statistics = FALSE,
    plot_hist = TRUE,
    plot_prior = TRUE
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
    parallel = TRUE
)
#---Plot posterior marginal distributions against other methods
plots <- plot_compare_marginal(
    plots = plots,
    abc_results = drf_results,
    parameters_labels = parameters_labels,
    parameters_truth = parameters_truth,
    plot_statistics = FALSE,
    plot_hist = TRUE
)
# ===============================================================ABC-SMC
# ---Run ABC-SMC
abc_smc_results <- abc_smc(
    statistics_target = statistics_target,
    model = model,
    parameters_labels = parameters_labels,
    prior_distributions = list(c("unif", -10, 10), c("unif", -10, 10)),
    nParticles = 1000,
    method = "Beaumont",
    tolerance = c(30.0, 16.0, 6.0, 5.0, 4.3),
    progress_bar = TRUE
)
#---Plot posterior marginal distributions against other methods
plots <- plot_compare_marginal(
    plots = plots,
    abc_results = abc_smc_results,
    parameters_labels = parameters_labels,
    parameters_truth = parameters_truth,
    plot_statistics = FALSE,
    plot_hist = TRUE
)
