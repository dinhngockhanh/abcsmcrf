# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
# R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/Birth_death/test_0411(try new model)/perturb=+-2;mtry=default;num.tree=1000;5000*4/test_on_targets"
# R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/Birth_death/test_0411(try new model)/perturb=+-2;mtry=default;num.tree=1000;5000*4/test/1010"
R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/Birth_death/test_0411(try new model)/perturb=+-2;mtry=default;num.tree=1000;5000*4/test_on_targets/previous one/try_final_result"
R_libPaths <- ""
R_libPaths_extra <- "/Users/xiangzijin/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zhihan - Macbook
# R_workplace <- "/Users/lexie/Documents/DNA/SMC-RF/vignettes/birth_death/test"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/lexie/Documents/DNA/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - HPC
# R_workplace <- getwd()
# R_libPaths <- "/burg/iicd/users/zx2406/rpackages"
# R_libPaths_extra <- "/burg/iicd/users/zx2406/R_smcrf"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(ggplot2)
library(gridExtra)
library(grid)
library(invgamma)
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)
# set.seed(1)
# =====================================Model for the birth-death process
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters, parallel = TRUE) {
    nTimes <- 25
    nNoise <- 0
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
            FUN = function(i) BD_model(parameters$lambda[i], parameters$mu[i], nTimes)
        )
        stopCluster(cl)
        stats <- rbindlist(stats)
        class(stats) <- "data.frame"
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            stats <- rbind(stats, BD_model(parameters$lambda[i], parameters$mu[i], nTimes))
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
BD_model <- function(lambda, mu, nTimes) {
    #   This simulates a birth-death process
    #   Input:  lambda = birth rate
    #           mu = death rate
    times <- seq(1, nTimes) / nTimes
    npt <- length(times)
    alpha <- rep(0, npt)
    beta <- rep(0, npt)
    tdiff <- c(times[1], diff(times))
    Fsim <- rep(0, npt)
    Zsim <- rep(0, npt)
    # compute alpha and beta vectors
    # if (lambda == mu) {
    #     alpha <- lambda * tdiff / (1 + lambda * tdiff)
    #     beta <- alpha
    # } else {
    elm <- exp((lambda - mu) * tdiff)
    alpha <- mu * (elm - 1) / (lambda * elm - mu)
    beta <- lambda * alpha / mu
    # }
    # Now generate observations on F and Z. We assume F > 0 at each stage,
    # so we are conditioning on Z > 0 at the end
    zv <- 10 # start from single individual
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
    # stats <- data.frame(matrix(c(lambda, mu, Zsim[1:npt]), nrow = 1))
    if (Zsim[npt] == 0) {
        stats <- data.frame(matrix(rep(NA, (npt + 2)), nrow = 1, ncol = npt + 2))
    } else {
        stats <- data.frame(matrix(c(lambda, mu, Zsim[1:npt]), nrow = 1))
    }
    colnames(stats) <- c("lambda", "mu", paste0("Z_", 1:npt))
    return(stats)
}
# =====================================================Target statistics
# 19.8740943334997 15.763969453983
# 4.44521791301668 2.23759585060179
# 17.7337687229738 13.2654156023636
# 13.7864806316793 4.14433028083295
# 16.7352986009791 1.64104933850467


# lambda = 4, mu = 1
# lambda = 10, mu = 7
# lambda = 16, mu = 13
# lambda = 16, mu = 3
# lambda = 16, mu = 10

parameters_target <- data.frame(
    lambda = 6,
    mu = 3
    # lambda = runif(1, 0, 20),
    # mu = runif(1, 0, 20)
    # lambda = 13.79,
    # mu = 4.14
)
statistics_target <- data.frame(matrix(c(
    11, 13, 16, 17, 17, 21, 23, 23, 28, 32, 36, 48, 61, 68, 80, 96, 123, 135, 161, 173, 187, 197, 221, 253, 289
), nrow = 1))
# colnames(statistics_target) <- c(paste0("Z_", 1:25))

# statistics_target <- model(parameters = parameters_target, parallel = FALSE)[-c(1:ncol(parameters_target))]
colnames(statistics_target) <- c(paste0("Z_", 1:25))
target <- c()
target$parameters_target <- parameters_target
target$statistics_target <- statistics_target
save(target, file = "target_data.rda")
# load("/Users/xiangzijin/Documents/ABC_SMCRF/Birth_death/test_0411(try new model)/perturb=+-2;mtry=default;num.tree=1000;5000*4/test_on_targets/target_data.rda")
# parameters_target <- target$parameters_target
# statistics_target <- target$statistics_target
# ======================================Model for parameter perturbation
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters, after perturbation
perturb <- function(parameters) {
    for (i in 1:ncol(parameters)) parameters[[i]] <- parameters[[i]] + runif(nrow(parameters), min = -2, max = 2)
    return(parameters)
}
# ======================================Define ranges for the parameters
range <- data.frame(
    parameter = c("lambda", "mu"),
    min = c(0, 0),
    max = c(20, 20)
    # max = c(10, 10)
)
# ========================================Initial guesses for parameters
# ====================================(sampled from prior distributions)
lambda <- runif(100000, 0, 20)
mu <- runif(100000, 0, 20)
# lambda <- runif(100000, 0, 10)
# mu <- runif(100000, 0, 10)
parameters_initial <- data.frame(
    lambda = lambda,
    mu = mu
)
# test <- model(parameters_initial)

# print(head(test))
# ====================================Labels for parameters in the plots
parameters_labels <- data.frame(
    parameter = c("lambda", "mu"),
    label = c(deparse(expression(lambda)), deparse(expression(mu)))
)
# # =====================================True joint density for parameters
# #---Define the negative log-likelihood function
# nLL <- function(param) { # this version edited for MCMC example
#     lambda <- param[1]
#     mu <- param[2]
#     zvals <- c(10, unlist(statistics_target))
#     zvals <- unname(zvals)
#     times <- seq(1, 25) / 25
#     npt <- length(times)
#     alpha <- rep(0, npt)
#     beta <- rep(0, npt)
#     tdiff <- c(times[1], diff(times))
#     # compute alpha and beta vectors
#     if (lambda == mu) {
#         alpha <- lambda * tdiff / (1 + lambda * tdiff)
#         beta <- alpha
#     } else {
#         elm <- exp((lambda - mu) * tdiff)
#         alpha <- mu * (elm - 1) / (lambda * elm - mu)
#         beta <- lambda * alpha / mu
#     }
#     totLL <- 0.0
#     # Now calculate log-likelihood from the decomposition
#     for (j in seq(1, npt)) {
#         n <- zvals[j]
#         m <- zvals[j + 1]
#         alph <- alpha[j]
#         bet <- beta[j]
#         mn <- min(m, n)
#         tot <- 0
#         for (l in 1:mn) {
#             tot <- tot + dbinom(l, size = n, prob = 1 - alph) *
dnbinom(m - l, size = l, prob = 1 - bet)
#         }
#         totLL <- totLL - log(tot) # negative log-likelihood
#     }
#     return(-totLL) # return log-likelihood
# }

# # lambdas <- seq(0.05, 10, 0.05)
# # mus <- seq(0.05, 10, 0.05)
# lambdas <- seq(0.1, 20, 0.1)
# mus <- seq(0.1, 20, 0.1)
# vec_lambda <- c()
# vec_mu <- c()
# vec_nLL <- c()
# cnt <- 0
# for (i in 1:length(lambdas)) {
#     for (j in 1:length(mus)) {
#         cnt <- cnt + 1
#         print(cnt)
#         print("=====")
#         print(length(lambdas) * length(mus))
#         lambda <- lambdas[i]
#         mu <- mus[j]
#         vec_lambda <- c(vec_lambda, lambda)
#         vec_mu <- c(vec_mu, mu)
#         vec_nLL <- c(vec_nLL, nLL(c(lambda, mu)))
#     }
# }
# max_likelihood <- max(vec_nLL)
# lln_weights <- exp(vec_nLL - max_likelihood)
# lln_weights <- lln_weights / sum(sort(lln_weights))
# df_nLL <- data.frame(lambda = vec_lambda, mu = vec_mu, likelihood = vec_nLL)
# sampled_parameters <- df_nLL[sample(nrow(df_nLL), 2000, prob = lln_weights, replace = T), ]
# parameters_truth <- data.frame(
#     lambda = sampled_parameters$lambda + runif(nrow(sampled_parameters), min = -0.05, max = 0.05),
#     mu = sampled_parameters$mu + runif(nrow(sampled_parameters), min = -0.05, max = 0.05)
#     # lambda = sampled_parameters$lambda + runif(nrow(sampled_parameters), min = -0.025, max = 0.25),
#     # mu = sampled_parameters$mu + runif(nrow(sampled_parameters), min = -0.025, max = 0.025)
# )
# true_joint <- c()
# true_joint$likelihood <- df_nLL
# true_joint$method <- "true-joint"
# true_joint$parameters_truth <- parameters_truth
# save(true_joint, file = "likelihood_test.rda")
# load("likelihood_test.rda")
#-----------------------------------------------------------------------
#---Load likelihood
# df_nLL <- load("/Users/xiangzijin/Documents/ABC_SMCRF/Birth_death/test_0411(try new model)/perturb=+-2;mtry=default;num.tree=1000;5000*4/test/1010/likelihood_1010.rda")
df_nLL <- load("/Users/xiangzijin/Documents/ABC_SMCRF/Birth_death/test_0411(try new model)/perturb=+-2;mtry=default;num.tree=1000;5000*4/final results for paper/likelihood.rda")
# df_nLL <- load("/Users/xiangzijin/Documents/ABC_SMCRF/Birth_death/test_0411(try new model)/perturb=+-2;mtry=default;num.tree=1000;5000*4/test_on_targets/likelihood_test.rda")
parameters_truth <- true_joint$parameters_truth
#---Plot posterior joint distributions against other methods
plots_joint <- plot_compare_joint(
    abc_results = true_joint,
    parameters_labels = parameters_labels
)
# =================================================================MCMC
# proposalfunction <- function(param) {
#     # prior for lambda is unif (a,b)
#     # prior for mu is unif (c,d)
#     a <- 0
#     b <- 20 # changed it from 10 to 20
#     c <- 0
#     d <- 20 # changed it from 10 to 20
#     gridl <- 20 # also needs changing in qxylogdens
#     # update lambda
#     span1 <- (b - a) / gridl
#     minlambda <- max(a, param[1] - span1)
#     maxlambda <- min(b, param[1] + span1)
#     newlambda <- runif(1, min = minlambda, max = maxlambda)
#     # print(c(span1,minlambda,maxlambda))
#     # update mu
#     span2 <- (d - c) / gridl
#     minmu <- max(c, param[2] - span2)
#     maxmu <- min(d, param[2] + span2)
#     newmu <- runif(1, min = minmu, max = maxmu)
#     # print(c(span2,minmu,maxmu))
#     return(c(newlambda, newmu))
# }

# qxylogdens <- function(x) {
#     # computes log-density of q(x to y) for uniform updates
#     qval <- vector(length = 2)
#     # prior for lambda is unif (a,b)
#     # prior for mu is unif (c,d)
#     a <- 0
#     b <- 20 # changed it from 10 to 20
#     c <- 0
#     d <- 20 # changed it from 10 to 20
#     gridl <- 20 # gridl must be at least 2; also change in proposalfunction
#     # update lambda
#     span1 <- (b - a) / gridl
#     if (a + span1 <= x[1] & x[1] <= b - span1) {
#         qval[1] <- 1 / (2 * span1)
#     }
#     if (x[1] <= a + span1) {
#         qval[1] <- 1 / (span1 + x[1] - a)
#     }
#     if (x[1] > b - span1) {
#         qval[1] <- 1 / (b - x[1] + span1)
#     }
#     # update mu
#     span2 <- (d - c) / gridl
#     if (c + span2 <= x[2] & x[2] <= d - span2) {
#         qval[2] <- 1 / (2 * span2)
#     }
#     if (x[2] <= c + span2) {
#         qval[2] <- 1 / (span2 + x[2] - a)
#     }
#     if (x[2] > b - span2) {
#         qval[2] <- 1 / (b - x[2] + span2)
#     }
#     return(log(qval))
# }


# run_metropolis_MCMC <- function(startvalue, iterations) {
#     # set.seed(301)
#     chain <- array(dim = c(iterations + 1, 2))
#     chain[1, ] <- startvalue
#     for (i in 1:iterations) {
#         print(i)
#         proposal <- proposalfunction(chain[i, ])
#         qratio <- exp(sum(qxylogdens(proposal)) - sum(qxylogdens(chain[i, ])))
#         probab <- exp(posterior(proposal) - posterior(chain[i, ])) * qratio
#         if (runif(1) < probab) {
#             chain[i + 1, ] <- proposal
#         } else {
#             chain[i + 1, ] <- chain[i, ]
#         }
#     }
#     return(chain)
# }

# # Prior distribution
# prior <- function(param) {
#     lambdaprior <- dunif(param[1], min = 0, max = 20, log = T)
#     muprior <- dunif(param[2], min = 0, max = 20, log = T)
#     return(lambdaprior + muprior)
# }

# posterior <- function(param) {
#     return(nLL(param) + prior(param))
# }

# startvalue <- c(parameters_target$lambda, parameters_target$mu)
# run10A <- run_metropolis_MCMC(startvalue, 100000)
# burnIn <- 2500
# acceptance <- 1 - mean(duplicated(run10A[-(1:burnIn), ]))
# vec <- seq(2501, 100000, 50)
# library(MASS)
# png(paste0("MCMC_lambda.png"))
# truehist(run10A[vec, 1], xlim = c(0, 20), xlab = "Value of lambda")
# dev.off()
# png(paste0("MCMC_mu.png"))
# truehist(run10A[vec, 2], xlim = c(0, 20), xlab = "Value of mu")
# dev.off()
# png(paste0("MCMC_contour.png"))
# x <- run10A[vec, 1]
# y <- run10A[vec, 2]
# twodplot <- kde2d(x, y, n = 1950, lims = c(range(x), range(y)))
# contour(twodplot$x, twodplot$y, twodplot$z, xlim = c(0, 20), ylim = c(0, 20), xlab = "lambda", ylab = "mu")
# dev.off()
# MCMC <- list()
# MCMC$start_value <- parameters_target
# MCMC$full_chain <- run10A
# MCMC$burnIn <- burnIn
# MCMC$acceptance <- acceptance
# MCMC$time_points <- vec
# selected_lambda <- run10A[vec, 1]
# selected_mu <- run10A[vec, 2]
# selected_params <- data.frame(lambda = selected_lambda, mu = selected_mu)
# MCMC$selected_params <- selected_params
# MCMC$fitted_modes <- data.frame(lambda = twodplot$x[which(twodplot$z == max(twodplot$z), arr.ind = T)[[1]]], mu = twodplot$y[which(twodplot$z == max(twodplot$z), arr.ind = T)[[2]]])
# MCMC[["method"]] <- "mcmc"
# reference <- model(selected_params)
# parameters_ids <- parameters_labels$parameter
# statistics <- data.frame(reference[, colnames(reference)[!colnames(reference) %in% parameters_ids]])
# colnames(statistics) <- colnames(reference)[!colnames(reference) %in% parameters_ids]
# MCMC$statistics <- statistics
# MCMC[["statistics_target"]] <- statistics_target
# MCMC[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
# MCMC[["statistics_labels"]] <- data.frame(ID = colnames(statistics_target)[!colnames(statistics_target) %in% parameters_ids])
# save(MCMC, file = "MCMC_test.rda")
# load("MCMC_test.rda")
MCMC_20_results <- load("/Users/xiangzijin/Documents/ABC_SMCRF/Birth_death/MCMC_2020.rda")
# # MCMC_10_results <- load("/Users/xiangzijin/Documents/ABC_SMCRF/Birth_death/MCMC_1010.rda")
# #---Plot marginal distributions compare
plots_marginal <- plot_compare_marginal(
    # plots = plots_marginal,
    abc_results = MCMC,
    parameters_labels = parameters_labels,
    parameters_truth = parameters_truth,
    plot_hist = TRUE
    # plot_hist_point = TRUE
)
plots_compare_qqplot <- plot_compare_qqplot(
    # plots = plot_compare_qqplot,
    abc_results = MCMC,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels
)
#---Plot posterior joint distributions against other methods
plots_joint <- plot_compare_joint(
    abc_results = MCMC,
    plots = plots_joint,
    parameters_labels = parameters_labels
)
# ===================================================================DRF
#---Run ABC-DRF
# abcrf_results <- smcrf(
#     method = "smcrf-multi-param",
#     statistics_target = statistics_target,
#     parameters_initial = parameters_initial,
#     model = model,
#     perturb = perturb,
#     range = range,
#     nParticles = rep(20000, 1),
#     num.trees = 1000,
#     parallel = TRUE
# )
# save(abcrf_results, file = "drf.rda")
load("drf.rda")
# results <- load("/Users/xiangzijin/Documents/ABC_SMCRF/Birth_death/test_0411(try new model)/perturb=+-2;mtry=default;num.tree=1000;5000*4/drf.rda")
# load("/Users/xiangzijin/Documents/ABC_SMCRF/Birth_death/test_0411(try new model)/perturb=+-2;mtry=default;num.tree=1000;5000*4/test/1010/drf.rda")
# ---Plot posterior joint distributions against other methods
plots_joint <- plot_compare_joint(
    plots = plots_joint,
    abc_results = abcrf_results,
    parameters_labels = parameters_labels
)
plots_compare_qqplot <- plot_compare_qqplot(
    plots = plots_compare_qqplot,
    abc_results = abcrf_results,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels
)
#---Plot marginal distributions compare
plots_marginal <- plot_compare_marginal(
    plots = plots_marginal,
    abc_results = abcrf_results,
    # parameters_truth = parameters_truth,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
    # plot_hist_point = TRUE
)
# ===================================================================RF
#---Run ABC-RF
# abcrf_results <- smcrf(
#     method = "smcrf-single-param",
#     statistics_target = statistics_target,
#     parameters_initial = parameters_initial,
#     model = model,
#     perturb = perturb,
#     range = range,
#     ntree = 2500,
#     nParticles = rep(20000, 1),
#     parallel = TRUE
# )
# save(abcrf_results, file = "abc-rf.rda")
load("abc-rf.rda")
#---Plot posterior joint distributions against other methods
# plots_joint <- plot_compare_joint(
#     plots = plots_joint,
#     abc_results = abcrf_results,
#     parameters_labels = parameters_labels
# )
#---Plot marginal distributions compare
# plots_marginal <- plot_compare_marginal(
#     plots = plots_marginal,
#     abc_results = abcrf_results,
#     parameters_labels = parameters_labels,
#     # plot_hist = TRUE
#     plot_hist_point = TRUE
# )
# plots_compare_qqplot <- plot_compare_qqplot(
#     plots = plots_compare_qqplot,
#     abc_results = abcrf_results,
#     parameters_truth = parameters_truth,
#     parameters_labels = parameters_labels
# )
# ========================================SMC-RF for multiple parameters
#---Run SMC-RF for multiple parameters
# smcrf_results_multi_param <- smcrf(
#     method = "smcrf-multi-param",
#     statistics_target = statistics_target,
#     parameters_initial = parameters_initial,
#     model = model,
#     perturb = perturb,
#     range = range,
#     nParticles = rep(5000, 4),
#     num.trees = 2000,
#     parallel = TRUE
# )
# save(smcrf_results_multi_param, file = "smc-drf.rda")
# smcrf_results_multi_param <- smcrf(
#     method = "smcrf-multi-param",
#     statistics_target = statistics_target,
#     parameters_initial = parameters_initial,
#     model = model,
#     perturb = perturb,
#     range = range,
#     nParticles = rep(4000, 5),
#     num.trees = 1500,
#     parallel = TRUE
# )
# save(smcrf_results_multi_param, file = "smc-drf.rda")
load("smc-drf.rda")
# smcdrf_results <- load("/Users/xiangzijin/Documents/ABC_SMCRF/Birth_death/test_0411(try new model)/perturb=+-2;mtry=default;num.tree=1000;5000*4/test/1010/smc-drf.rda")
# smcdrf_results <- load("/Users/xiangzijin/Documents/ABC_SMCRF/Birth_death/test_0411(try new model)/perturb=+-2;mtry=default;num.tree=1000;5000*4/final results for paper/likelihood.rda/smc-drf.rda")
#---Plot joint distributions
plot_smcrf_joint(
    smcrf_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels
)
#---Plot posterior joint distributions against other methods
plots_joint <- plot_compare_joint(
    plots = plots_joint,
    abc_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels
)
#---Plot marginal distributions compare
plots_marginal <- plot_compare_marginal(
    plots = plots_marginal,
    abc_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
    # plot_hist_point = TRUE
)
plots_compare_qqplot <- plot_compare_qqplot(
    plots = plots_compare_qqplot,
    abc_results = smcrf_results_multi_param,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels
)
plot_smcrf_marginal(
    smcrf_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels,
    plot_hist = TRUE
)
# ========================================SMC-RF for single parameters
#---Run SMC-RF for single parameters
# smcrf_results_single_param <- smcrf(
#     method = "smcrf-single-param",
#     statistics_target = statistics_target,
#     parameters_initial = parameters_initial,
#     model = model,
#     perturb = perturb,
#     range = range,
#     nParticles = rep(2000, 5),
#     parallel = TRUE
# )
# save(smcrf_results_single_param, file = "smc-rf.rda")
# #---Plot joint distributions
# plot_smcrf_joint(
#     smcrf_results = smcrf_results_single_param,
#     parameters_labels = parameters_labels
# )
# #---Plot posterior joint distributions against other methods
# plots_joint <- plot_compare_joint(
#     plots = plots_joint,
#     abc_results = smcrf_results_single_param,
#     parameters_labels = parameters_labels
# )
# #---Plot marginal distributions compare
# plots_marginal <- plot_compare_marginal(
#     plots = plots_marginal,
#     abc_results = smcrf_results_single_param,
#     parameters_labels = parameters_labels,
#     parameters_truth = parameters_truth,
#     plot_hist = TRUE
# )
# plot_smcrf_marginal(
#     smcrf_results = smcrf_results_single_param,
#     parameters_labels = parameters_labels,
#     parameters_truth = parameters_truth,
#     plot_hist = TRUE
# )
# ========================================QQ-plots
# truth_lambda <- parameters_truth$lambda
# truth_mu <- parameters_truth$mu
# # smcrf_lambda <- smcrf_results_multi_param$Iteration_5$parameters$lambda
# # smcrf_mu <- smcrf_results_multi_param$Iteration_5$parameters$mu
# drf_lambda <- abcrf_results$Iteration_2$parameters$lambda
# # drf_mu <- abcrf_results$Iteration_2$parameters$mu
# mcmc_lambda <- MCMC$selected_params$lambda
# mcmc_mu <- MCMC$selected_params$mu
# qqplot(mcmc_lambda, truth_lambda)
# abline(0, 1)
# qqplot(mcmc_mu, truth_mu)
# abline(0, 1)
# qqplot(truth_lambda, smcrf_lambda)
# abline(0, 1)
# qqplot(truth_lambda, drf_lambda)
# abline(0, 1)
# qqplot(truth_mu, smcrf_mu)
# abline(0, 1)
# qqplot(truth_mu, drf_mu)
# abline(0, 1)
# qqplot(mcmc_lambda, smcrf_lambda)
# abline(0, 1)
# qqplot(mcmc_lambda, drf_lambda)
# abline(0, 1)
# qqplot(mcmc_mu, smcrf_mu)
# abline(0, 1)
# qqplot(mcmc_mu, drf_mu)
# abline(0, 1)

######## Metropolis algorithm ################
# proposalfunction <- function(param) {
#     # prior for lambda is unif (a,b)
#     # prior for mu is unif (c,d)
#     a <- 0
#     b <- 20 # changed it from 10 to 20
#     c <- 0
#     d <- 20 # changed it from 10 to 20
#     gridl <- 20 # also needs changing in qxylogdens
#     # update lambda
#     span1 <- (b - a) / gridl
#     minlambda <- max(a, param[1] - span1)
#     maxlambda <- min(b, param[1] + span1)
#     newlambda <- runif(1, min = minlambda, max = maxlambda)
#     # print(c(span1,minlambda,maxlambda))

#     # update mu
#     span2 <- (d - c) / gridl
#     minmu <- max(c, param[2] - span2)
#     maxmu <- min(d, param[2] + span2)
#     newmu <- runif(1, min = minmu, max = maxmu)
#     # print(c(span2,minmu,maxmu))
#     return(c(newlambda, newmu))
# }

# qxylogdens <- function(x) {
#     # computes log-density of q(x to y) for uniform updates
#     qval <- vector(length = 2)
#     # prior for lambda is unif (a,b)
#     # prior for mu is unif (c,d)
#     a <- 0
#     b <- 20 # changed it from 10 to 20
#     c <- 0
#     d <- 20 # changed it from 10 to 20
#     gridl <- 20 # gridl must be at least 2; also change in proposalfunction

#     # update lambda
#     span1 <- (b - a) / gridl
#     if (a + span1 <= x[1] & x[1] <= b - span1) {
#         qval[1] <- 1 / (2 * span1)
#     }
#     if (x[1] <= a + span1) {
#         qval[1] <- 1 / (span1 + x[1] - a)
#     }
#     if (x[1] > b - span1) {
#         qval[1] <- 1 / (b - x[1] + span1)
#     }

#     # update mu
#     span2 <- (d - c) / gridl
#     if (c + span2 <= x[2] & x[2] <= d - span2) {
#         qval[2] <- 1 / (2 * span2)
#     }
#     if (x[2] <= c + span2) {
#         qval[2] <- 1 / (span2 + x[2] - a)
#     }
#     if (x[2] > b - span2) {
#         qval[2] <- 1 / (b - x[2] + span2)
#     }

#     return(log(qval))
# }


# run_metropolis_MCMC <- function(startvalue, iterations) {
#     # set.seed(301)
#     chain <- array(dim = c(iterations + 1, 2))
#     chain[1, ] <- startvalue
#     for (i in 1:iterations) {
#         proposal <- proposalfunction(chain[i, ])
#         qratio <- exp(sum(qxylogdens(proposal)) - sum(qxylogdens(chain[i, ])))
#         probab <- exp(posterior(proposal) - posterior(chain[i, ])) * qratio
#         if (runif(1) < probab) {
#             chain[i + 1, ] <- proposal
#         } else {
#             chain[i + 1, ] <- chain[i, ]
#         }
#     }
#     return(chain)
# }

# # Prior distribution
# prior <- function(param) {
#     lambdaprior <- dunif(param[1], min = 0, max = 20, log = T)
#     muprior <- dunif(param[2], min = 0, max = 20, log = T)
#     return(lambdaprior + muprior)
# }

# posterior <- function(param) {
#     return(nLL(param) + prior(param))
# }

# # ***************************
# # Load necessary libraries
# library(ggplot2)
# library(RColorBrewer)

# # Define the negative log-likelihood function
# nLL <- function(param) { # this version edited for MCMC example
#     lambda <- param[1]
#     mu <- param[2]

#     # note: zvals vector includes Z_0 as first value
#     # nLL is - log-likelihood

#     # set the data (fix the input later!)
#     # avals <-c(8,11,11,15,14,17,19,20,21,27,30,35,42,54,67,72,93,111,119,146,157,163,173,202,230)
#     zvals <- c(10, 11, 13, 16, 17, 17, 21, 23, 23, 28, 32, 36, 48, 61, 68, 80, 96, 123, 135, 161, 173, 187, 197, 221, 253, 289)
#     times <- seq(1, 25) / 25
#     # zvals <- c(10,48,289)
#     # times <- c(0.5,1)
#     # zvals <- c(10,289)
#     # times <- c(1.0)

#     npt <- length(times)
#     alpha <- rep(0, npt)
#     beta <- rep(0, npt)
#     tdiff <- c(times[1], diff(times))

#     # compute alpha and beta vectors

#     if (lambda == mu) {
#         alpha <- lambda * tdiff / (1 + lambda * tdiff)
#         beta <- alpha
#     } else {
#         elm <- exp((lambda - mu) * tdiff)
#         alpha <- mu * (elm - 1) / (lambda * elm - mu)
#         beta <- lambda * alpha / mu
#     }
#     # print(alpha)
#     # print(beta)

#     totLL <- 0.0

#     # Now calculate log-likelihood from the decomposition

#     for (j in seq(1, npt)) {
#         n <- zvals[j]
#         m <- zvals[j + 1]
#         alph <- alpha[j]
#         bet <- beta[j]

#         mn <- min(m, n)
#         tot <- 0
#         for (l in 1:mn) {
#             tot <- tot + dbinom(l, size = n, prob = 1 - alph) * dnbinom(m - l, size = l, prob = 1 - bet)
#         }
#         # print(log(tot))
#         totLL <- totLL - log(tot) # negative log-likelihood
#     }
#     # print(totLL)
#     return(-totLL) # return log-likelihood
# }

# # Define parameter ranges for lambda and mu
# lambdas <- seq(0.1, 20, 0.1)
# mus <- seq(0.1, 20, 0.1)

# # Create a grid of lambda and mu values
# vec_lambda <- rep(lambdas, each = length(mus))
# vec_mu <- rep(mus, times = length(lambdas))
# vec_nLL <- numeric(length(vec_lambda))

# # Calculate nLL for each combination of lambda and mu
# for (i in 1:length(vec_lambda)) {
#     vec_nLL[i] <- nLL(c(vec_lambda[i], vec_mu[i]))
# }

# # Convert negative log-likelihood to likelihood and scale it to 0-10
# min_likelihood <- min(vec_nLL)
# max_likelihood <- max(vec_nLL)
# scaled_loglikelihood <- ((vec_nLL - min_likelihood) / (max_likelihood - min_likelihood))^1000

# # Prepare the data frame for plotting
# df_nLL <- data.frame(lambda = vec_lambda, mu = vec_mu, likelihood = scaled_loglikelihood)
# p <- ggplot() +
#     geom_tile(data = df_nLL, aes(x = lambda, y = mu, fill = likelihood), show.legend = FALSE) +
#     scale_fill_gradientn(colors = brewer.pal(10, "RdYlBu")) +
#     theme_minimal() +
#     labs(x = "Lambda", y = "Mu", fill = "Scaled Likelihood", title = "Heatmap of Scaled Log-Likelihood")
# png("test.png")
# # Print the plot
# print(p)
# dev.off()


# +++++++++++++++
# startvalue <- c(6, 3)
# run10A <- run_metropolis_MCMC(startvalue, 100000)
# burnIn <- 2500
# acceptance <- 1 - mean(duplicated(run10A[-(1:burnIn), ]))
# vec <- seq(2501, 100000, 50)

# png(paste0("MCMC_lambda.png"))
# truehist(run10A[vec, 1], xlim = c(0, 20), xlab = "Value of lambda")
# dev.off()
# png(paste0("MCMC_mu.png"))
# truehist(run10A[vec, 2], xlim = c(0, 20), xlab = "Value of mu")
# dev.off()
# png(paste0("MCMC_contour.png"))
# x <- run10A[vec, 1]
# y <- run10A[vec, 2]
# twodplot <- kde2d(x, y, n = 1950, lims = c(range(x), range(y)))
# contour(twodplot$x, twodplot$y, twodplot$z, xlim = c(2, 10), ylim = c(0, 8), xlab = "lambda", ylab = "mu")
# dev.off()
# which(twodplot$z == max(twodplot$z), arr.ind = T)
# twodplot$x[336]
# twodplot$y[307]
# which(twodplot$z == max(twodplot$z), arr.ind = T)
# MCMC <- list()
# MCMC$start_value <- data.frame(
#     lambda = 6,
#     mu = 3
# )
# MCMC$full_chain <- full_chain
# MCMC$burnIn <- burnIn
# MCMC$acceptance <- acceptance
# MCMC$time_points <- vec
# selected_lambda <- run10A[vec, 1]
# selected_mu <- run10A[vec, 2]
# selected_params <- data.frame(lambda = selected_lambda, mu = selected_mu)
# MCMC$selected_params <- selected_params
# MCMC$fitted_modes <- data.frame(lambda = twodplot$x[614], mu = twodplot$y[535])
# MCMC[["method"]] <- "mcmc"
# reference <- model(selected_params)
# parameters_ids <- parameters_labels$parameter
# statistics <- data.frame(reference[, colnames(reference)[!colnames(reference) %in% parameters_ids]])
# colnames(statistics) <- colnames(reference)[!colnames(reference) %in% parameters_ids]
# MCMC$statistics <- statistics
# MCMC[["statistics_target"]] <- statistics_target
# MCMC[["parameters_labels"]] <- data.frame(parameter = parameters_ids)
# MCMC[["statistics_labels"]] <- data.frame(ID = colnames(statistics_target)[!colnames(statistics_target) %in% parameters_ids])
# save(MCMC, file = "MCMC.rda")
# Sample data for demonstration
