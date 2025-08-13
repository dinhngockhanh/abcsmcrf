library(abcsmcrf)
set.seed(1)
# =====================================Model for the birth-death process
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters, parallel = TRUE) {
    nTimes <- 25
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
            FUN = function(i) BD_model(parameters$lambda_minus_mu[i], parameters$mu_over_lambda[i], nTimes)
        )
        stopCluster(cl)
        stats <- rbindlist(stats)
        class(stats) <- "data.frame"
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            stats <- rbind(stats, BD_model(parameters$lambda_minus_mu[i], parameters$mu_over_lambda[i], nTimes))
        }
    }
    data <- data.frame(stats)
    return(data)
}
BD_model <- function(lambda_minus_mu, mu_over_lambda, nTimes) {
    #   This simulates a birth-death process
    #   Input:  lambda = birth rate
    #           mu = death rate
    mu <- (lambda_minus_mu * mu_over_lambda) / (1 - mu_over_lambda)
    lambda <- mu + lambda_minus_mu
    times <- seq(1, nTimes) / nTimes
    npt <- length(times)
    alpha <- rep(0, npt)
    beta <- rep(0, npt)
    tdiff <- c(times[1], diff(times))
    Fsim <- rep(0, npt)
    Zsim <- rep(0, npt)
    elm <- exp((lambda - mu) * tdiff)
    alpha <- mu * (elm - 1) / (lambda * elm - mu)
    beta <- lambda * alpha / mu
    #   Now generate observations on F and Z.
    #   We condition on F > 0 at each stage, i.e. Z_25 > 0
    zv <- 10
    Fsim[1] <- rbinom(1, size = zv, prob = 1 - alpha[1])
    if (Fsim[1] == 0) {
        Zsim[1] <- 0
    } else {
        Zsim[1] <- rnbinom(1, size = Fsim[1], prob = 1 - beta[1]) + Fsim[1]
    }
    for (j in seq(2, npt)) {
        Fsim[j] <- rbinom(1, size = Zsim[j - 1], prob = 1 - alpha[j])
        if (Fsim[j] == 0) {
            Zsim[j] <- 0
        } else {
            Zsim[j] <- rnbinom(1, size = Fsim[j], prob = 1 - beta[j]) + Fsim[j]
        }
    }
    if (Zsim[npt] == 0) {
        stats <- data.frame(matrix(rep(NA, (npt + 2)), nrow = 1, ncol = npt + 2))
    } else {
        stats <- data.frame(matrix(c(lambda_minus_mu, mu_over_lambda, Zsim[1:npt]), nrow = 1))
    }
    colnames(stats) <- c("lambda_minus_mu", "mu_over_lambda", paste0("Z_", 1:npt))
    return(stats)
}
# =====================================================Target statistics
statistics_target <- data.frame(matrix(c(
    11, 13, 16, 17, 17, 21, 23, 23, 28, 32, 36, 48, 61, 68, 80, 96, 123, 135, 161, 173, 187, 197, 221, 253, 289
), nrow = 1))
colnames(statistics_target) <- c(paste0("Z_", 1:25))
parameters_target <- data.frame(
    lambda_minus_mu = 3,
    mu_over_lambda = 0.5
)
# ====================================Labels for parameters in the plots
parameters_labels <- data.frame(
    parameter = c("lambda_minus_mu", "mu_over_lambda"),
    label = c(
        deparse(expression(lambda - mu)),
        deparse(expression(mu / lambda))
    )
)
# =====================================True joint density for parameters
#---Define the negative log-likelihood function
nLL <- function(param) {
    #   Changes for reparametrization
    lambda_minus_mu <- param[1]
    mu_over_lambda <- param[2]
    if (mu_over_lambda == 1) {
        return(-log(0))
    } else {
        mu <- (lambda_minus_mu * mu_over_lambda) / (1 - mu_over_lambda)
        lambda <- mu + lambda_minus_mu
        zvals <- c(10, unlist(statistics_target))
        zvals <- unname(zvals)
        times <- seq(1, 25) / 25
        npt <- length(times)
        alpha <- rep(0, npt)
        beta <- rep(0, npt)
        tdiff <- c(times[1], diff(times))
        #   Compute alpha and beta vectors
        if (lambda == mu) {
            alpha <- lambda * tdiff / (1 + lambda * tdiff)
            beta <- alpha
        } else {
            elm <- exp((lambda - mu) * tdiff)
            alpha <- mu * (elm - 1) / (lambda * elm - mu)
            beta <- lambda * alpha / mu
        }
        totLL <- 0.0
        #   Now calculate log-likelihood from the decomposition
        for (j in seq(1, npt)) {
            n <- zvals[j]
            m <- zvals[j + 1]
            alph <- alpha[j]
            bet <- beta[j]
            mn <- min(m, n)
            tot <- 0
            for (l in 1:mn) {
                tot <- tot + dbinom(l, size = n, prob = 1 - alph) *
                    dnbinom(m - l, size = l, prob = 1 - bet)
            }
            totLL <- totLL - log(tot)
        }
        return(-totLL)
    }
}
#---Define the prior distribution
shape <- 3
rate <- 1
shape1 <- 5
shape2 <- 5
prior <- function(param) {
    lambda_minus_mu <- param[1]
    mu_over_lambda <- param[2]
    if ((lambda_minus_mu < 0) || (mu_over_lambda > 1)) {
        return(log(0))
    } else {
        return(dgamma(lambda_minus_mu, shape = shape, rate = rate, log = T) + dbeta(mu_over_lambda, shape1 = shape1, shape2 = shape2))
    }
}
#---Compute the posterior distribution
posterior <- function(param) {
    return(nLL(param) + prior(param))
}
lambdas_minus_mu <- seq(0.1, 20, 0.1)
mus_over_lambda <- seq(0.1, 0.99, 0.001)
vec_lambda <- c()
vec_mu <- c()
vec_nLL <- c()
cnt <- 0
for (i in 1:length(lambdas_minus_mu)) {
    for (j in 1:length(mus_over_lambda)) {
        cnt <- cnt + 1
        lambda <- lambdas_minus_mu[i]
        mu <- mus_over_lambda[j]
        vec_lambda <- c(vec_lambda, lambda)
        vec_mu <- c(vec_mu, mu)
        vec_nLL <- c(vec_nLL, posterior(c(lambda, mu)))
    }
}
max_likelihood <- max(vec_nLL)
lln_weights <- exp(vec_nLL - max_likelihood)
lln_weights <- lln_weights / sum(sort(lln_weights))
df_nLL <- data.frame(lambda_minus_mu = vec_lambda, mu_over_lambda = vec_mu, likelihood = vec_nLL)
sampled_parameters <- df_nLL[sample(nrow(df_nLL), 2000, prob = lln_weights, replace = T), ]
parameters_truth <- data.frame(
    lambda_minus_mu = sampled_parameters$lambda_minus_mu + runif(nrow(sampled_parameters), min = -0.05, max = 0.05),
    mu_over_lambda = sampled_parameters$mu_over_lambda + runif(nrow(sampled_parameters), min = -0.05, max = 0.05)
)
true_joint <- c()
true_joint$likelihood <- df_nLL
true_joint$method <- "true-joint"
true_joint$parameters_truth <- parameters_truth
#---Plot posterior joint distributions against other methods
plots_joint <- plot_compare_joint(
    abc_results = true_joint,
    parameters_labels = parameters_labels
)
# ====================================================Prior distribution
dprior <- function(parameters, parameter_id = "all") {
    probs <- rep(1, nrow(parameters))
    if (parameter_id %in% c("all", "lambda_minus_mu")) {
        probs <- probs * dgamma(parameters[["lambda_minus_mu"]], shape = shape, rate = rate)
    }
    if (parameter_id %in% c("all", "mu_over_lambda")) {
        probs <- probs * dbeta(parameters[["mu_over_lambda"]], shape1 = shape1, shape2 = shape2)
    }
    return(probs)
}
rprior <- function(Nparameters) {
    lambda_minus_mu <- rgamma(Nparameters, shape = shape, rate = rate)
    mu_over_lambda <- rbeta(Nparameters, shape1 = shape1, shape2 = shape2)
    return(data.frame(lambda_minus_mu = lambda_minus_mu, mu_over_lambda = mu_over_lambda))
}
# ============================================================ABC-SMC-RF
#---Run ABC-RF
abcrf_results <- smcrf(
    method = "smcrf-single-param",
    statistics_target = statistics_target,
    model = model,
    rprior = rprior,
    dprior = dprior,
    perturbation = "Uniform",
    perturbation_parameters = data.frame(
        lambda_minus_mu = rep(1, 4),
        mu_over_lambda = rep(0.05, 4)
    ),
    nParticles = rep(5000, 4),
    ntree = 2500,
    save_model = FALSE,
    parallel = TRUE,
    model_redo_if_NA = TRUE
)
#---Plot qqplots against other methods
plots_compare_qqplot <- plot_compare_qqplot(
    abc_results = abcrf_results,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels
)
#---Plot marginal distributions against other methods
plots_marginal <- plot_compare_marginal(
    abc_results = abcrf_results,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)
# ===========================================================ABC-SMC-DRF
smcrf_results_multi_param <- smcrf(
    method = "smcrf-multi-param",
    statistics_target = statistics_target,
    model = model,
    rprior = rprior,
    dprior = dprior,
    perturbation = "Uniform",
    perturbation_parameters = data.frame(
        lambda_minus_mu = rep(1, 4),
        mu_over_lambda = rep(0.05, 4)
    ),
    nParticles = rep(5000, 4),
    save_model = FALSE,
    num.trees = 2500,
    splitting.rule = "FourierMMD",
    model_redo_if_NA = TRUE,
    parallel = TRUE
)
#---Plot posterior joint distributions against other methods
plots_joint <- plot_compare_joint(
    plots = plots_joint,
    lims = data.frame(
        parameter = c("lambda_minus_mu", "mu_over_lambda"),
        min = c(1, 0),
        max = c(5, 1)
    ),
    abc_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels,
    nBins = 8
)
#---Plot marginal distributions against other methods
plots_marginal <- plot_compare_marginal(
    plots = plots_marginal,
    abc_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)
#---Plot qqplots against other methods
plots_compare_qqplot <- plot_compare_qqplot(
    plots = plots_compare_qqplot,
    abc_results = smcrf_results_multi_param,
    lims = data.frame(
        parameter = c("lambda_minus_mu", "mu_over_lambda"),
        min = c(1, 0),
        max = c(5, 1)
    ),
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels
)
