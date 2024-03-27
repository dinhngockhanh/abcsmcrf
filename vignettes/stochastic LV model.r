# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zhihan - Macbook
R_workplace <- "/Users/lexie/Documents/DNA/SMC-RF/vignettes/LV model"
R_libPaths <- ""
R_libPaths_extra <- "/Users/lexie/Documents/DNA/SMC-RF/R"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)
library(ggplot2)
library(gridExtra)
library(grid)
library(invgamma)
library(SimBIID)
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)

set.seed(1)

# ======================================LV
library(GillespieSSA)
parms <- c(c1=10, c2=0.01, c3=10)
x0 <- c(X=1000,Y=1000)
a <- c(
    "c1*X",
    "c2*X*Y",
    "c3*Y")
nu <- matrix(c(+1, -1, 0, 0, +1, -1), nrow = 2, byrow = TRUE)
out <- ssa(x0,a,nu,parms,tf=1,method=ssa.d(),simName="Predator-prey model")
ssa.plot(out)

# =====================================Model for the predator–prey process
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters, parallel = TRUE) {
    nNoise <- 0
    #   Make simulations & compute summary statistics (population sizes at each time point)
    if (parallel) {
        library(parallel)
        library(pbapply)
        library(data.table)
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cl, varlist = c("LV_model"))
        stats <- pblapply(
            cl = cl, X = 1:nrow(parameters),
            FUN = function(i) LV_model(parameters$c1[i], parameters$c2[i], parameters$c3[i])
        )
        stopCluster(cl)
        stats <- rbindlist(stats)
        class(stats) <- "data.frame"
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            stats <- rbind(stats, LV_model(parameters$c1[i], parameters$c2[i], parameters$c3[i]))
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
LV_model <- function(c1, c2, c3) {
    #   This simulates a predator–prey process
    #   Input: 
    # 
    nRuns <- 3
    library(GillespieSSA)          
    parms <- c(c1=c2, c2=c2, c3=c3)
    x0 <- c(X=1000,Y=1000)
    a <- c(
        "c1*X",
        "c2*X*Y",
        "c3*Y")
    nu <- matrix(c(+1, -1, 0, 0, +1, -1), nrow = 2, byrow = TRUE)
    Xsim <- numeric(length = 19)
    Ysim <- numeric(length = 19)
    for (i in nRuns){
        out <- ssa(x0,a,nu,parms,tf=1,method=ssa.d(),simName="Predator-prey model")
        nSteps <- out$stat$nSteps
        sample_id <- seq(1, nSteps + 2, by = nSteps %/% 20)[2:20]
        X <- out$data[sample_id,"X"]
        Y <- out$data[sample_id,"Y"]
        Xsim <- (Xsim + X)/nRuns
        Ysim <- (Ysim + Y)/nRuns
    }
    stats <- data.frame(matrix(c(c1, c2, c3, Xsim, Ysim), nrow = 1))
    colnames(stats) <- c("c1", "c2", "c3", paste0("X_", 1:19), paste0("Y_", 1:19))
    return(stats)
}
# =====================================================Target statistics
parameters_truth <- data.frame(
    c1 = 10,
    c2 = 0.01,
    c3 = 10
)
statistics_target <- model(parameters = parameters_truth, parallel = FALSE)[-c(1:ncol(parameters_truth))]
# ======================================Model for parameter perturbation
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters, after perturbation
perturb <- function(parameters) {
    parameters[[1]] <- parameters[[1]] + runif(nrow(parameters), min = 0, max = 1.0)
    parameters[[2]] <- parameters[[2]] + runif(nrow(parameters), min = 0, max = 0.0025)
    parameters[[3]] <- parameters[[3]] + runif(nrow(parameters), min = 0, max = 1.0)
    return(parameters)
}
# ======================================Define ranges for the parameters
range <- data.frame(
    parameter = c("c1", "c2", "c3"),
    min = c(0, 0, 0),
    max = c(10, 0.01, 10)
)
# ========================================Initial guesses for parameters
# ====================================(sampled from prior distributions)
c1 <- runif(100, 0, 28)
c2 <- runif(100, 0, 0.04)
c3 <- runif(100, 0, 28)
parameters_initial <- data.frame(
    c1 = c1,
    c2 = c2,
    c3 = c3
)
# ====================================Labels for parameters in the plots
parameters_labels <- data.frame(
    parameter = c("c1", "c2", "c3"),
    label = c(deparse(expression(c[1])), deparse(expression(c[2])), deparse(expression(c[3])))
)
# ========================================SMC-RF for multiple parameters
#---Run SMC-RF for multiple parameters
smcrf_results_multi_param <- smcrf(
    method = "smcrf-multi-param",
    statistics_target = statistics_target,
    parameters_initial = parameters_initial,
    model = model,
    perturb = perturb,
    range = range,
    nParticles = rep(100, 5),
    parallel = TRUE
)
#---Plot marginal distributions
plot_smcrf_marginal(
    smcrf_results = smcrf_results_multi_param,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels
)
#---Plot joint distributions
plot_smcrf_joint(
    smcrf_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels[1:2,]
)
plot_smcrf_joint(
    smcrf_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels[2:3,]
)
plot_smcrf_joint(
    smcrf_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels[c(1,3),]
)
#---Plot posterior joint distributions against other methods
plots.1 <- plot_compare_joint(
    abc_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels[1:2,]
)
plots.2 <- plot_compare_joint(
    abc_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels[2:3,]
)
plots.3 <- plot_compare_joint(
    abc_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels[c(1,3),]
)
# ========================================DRF
#---Run DRF
abcrf_results <- smcrf(
    method = "smcrf-multi-param",
    statistics_target = statistics_target,
    parameters_initial = parameters_initial,
    model = model,
    perturb = perturb,
    range = range,
    nParticles = rep(100, 1),
    parallel = TRUE
)
#---Plot posterior joint distributions against other methods
plots.1 <- plot_compare_joint(
    plots = plots.1,
    abc_results = abcrf_results,
    parameters_labels = parameters_labels[1:2,]
)
plots.2 <- plot_compare_joint(
    abc_results = abcrf_results,
    parameters_labels = parameters_labels[2:3,]
)
plots.3 <- plot_compare_joint(
    abc_results = abcrf_results,
    parameters_labels = parameters_labels[c(1,3),]
)
#---Plot marginal distributions
plot_smcrf_marginal(
    smcrf_results = abcrf_results,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels
)
# ===============================================================ABC-SMC
#---Run ABC-SMC
# abc_smc_results <- abc_smc(
#     statistics_target = statistics_target,
#     model = model,
#     parameters_labels = parameters_labels,
#     prior_distributions = list(c("unif", 0, 28),c("unif", 0, 0.04),c("unif", 0, 28)), 
#     nParticles = 10, method = "Beaumont", progress_bar = TRUE,
#     tolerance = c(4000, 2900, 2000, 1900, 1800),
#     dist_weights = rep(1, ncol(statistics_target))
# )
#---Plot posterior marginal distributions against other methods
# plots <- plot_compare_marginal(
#     plots = plots,
#     abc_results = abc_smc_results,
#     parameters_labels = parameters_labels,
#     plot_statistics = TRUE
# )