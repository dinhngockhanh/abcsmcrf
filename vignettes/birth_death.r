# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zhihan - Macbook
R_workplace <- "/Users/lexie/Documents/DNA/SMC-RF/Results for birth death model"
R_libPaths <- ""
R_libPaths_extra <- "/Users/lexie/Documents/DNA/SMC-RF/R"
# =======================================SET UP FOLDER PATHS & LIBRARIES
.libPaths(R_libPaths)

library(ggplot2)
library(gridExtra)
library(grid)

setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)

# ======================================================================
nNoise <- 20
N <- 1000
set.seed(1)
BD_model <- function(lambda, mu, n_times) {
    #   This simulates a birth-death process
    #   Input:  lambda = birth rate
    #           mu = death rate
    times <- seq(1, n_times) / n_times
    npt <- length(times)
    alpha <- rep(0, npt)
    beta <- rep(0, npt)
    tdiff <- c(times[1], diff(times))
    Fsim <- rep(0, npt)
    Zsim <- rep(0, npt)
    # compute alpha and beta vectors
    if (lambda == mu) {
        alpha <- lambda * tdiff / (1 + lambda * tdiff)
        beta <- alpha
    } else {
        elm <- exp((lambda - mu) * tdiff)
        alpha <- mu * (elm - 1) / (lambda * elm - mu)
        beta <- lambda * alpha / mu
    }
    # Now generate observations on F and Z. We assume F > 0 at each stage,
    # so we are conditioning on Z > 0 at the end
    zv <- 1 # start from single individual
    Fsim[1] <- sample(1:zv, size = 1, prob = dbinom(1:zv, size = zv, 1 - alpha[1]))
    Zsim[1] <- rnbinom(1, size = Fsim[1], prob = 1 - beta[1]) + Fsim[1]
    for (j in seq(2, npt)) {
        zv <- Zsim[j - 1]
        Fsim[j] <- sample(1:zv, size = 1, prob = dbinom(1:zv, size = zv, 1 - alpha[j]))
        Zsim[j] <- rnbinom(1, size = Fsim[j], prob = 1 - beta[j]) + Fsim[j]
    }

    stats <- data.frame(matrix(c(lambda, mu, Zsim[1:npt]), nrow = 1))
    colnames(stats) <- c("lambda", "mu", paste0("Z_", 1:npt))
    # stats <- data.frame(matrix(c(lambda, mu, Fsim[1:npt], Zsim[1:npt]), nrow = 1))
    # colnames(stats) <- c("lambda", "mu", paste0("F_", 1:npt), paste0("Z_", 1:npt))

    return(stats)
}
#---Model for the birth-death process
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters,
                  n_samples_per_parameter_set,
                  nNoise) {
    library(parallel)
    library(pbapply)
    library(data.table)
    lambda <- parameters$lambda
    mu <- parameters$mu

    n_times <- 50

    nNoise <- nNoise # Number of noise variables

    cl <- makePSOCKcluster(detectCores() - 1)
    BD_model <<- BD_model
    clusterExport(cl, varlist = c("BD_model"))
    stats <- pblapply(cl = cl, X = 1:nrow(parameters), FUN = function(i) BD_model(lambda[i], mu[i], n_times))
    stats <- rbindlist(stats)
    class(stats) <- "data.frame"
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
#---Model for parameter perturbation
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters, after perturbation
perturb <- function(parameters) {
    parameters$lambda <- pmin(pmax(parameters$lambda + runif(nrow(parameters), min = -1, max = 1), 0.0001), 15)
    parameters$mu <- pmin(pmax(parameters$mu + runif(nrow(parameters), min = -1, max = 1), 0.0001), 15)
    return(parameters)
}
#---Target statistics
lambda <- runif(1, 0, 15)
mu <- runif(1, 0, lambda)
target <- model(data.frame(lambda = lambda, mu = mu), 
                n_samples_per_parameter_set=1, 
                nNoise=nNoise)[-c(1:2)]
cat(paste0("\n\n\n\n\nTRUE VALUE FOR THETA = ", lambda, "\n"))
cat(paste0("TRUE VALUE FOR MU = ", mu, "\n"))
cat(paste0("STATISTICS:\n"))
print(target)
cat("\n\n\n\n\n")
#---Initial guesses for parameters (sampled from prior distributions)
lambda <- runif(10000, 0, 15)
mu <- runif(10000, 0, lambda)
parameters_initial <- data.frame(lambda = lambda, mu = mu)



# #---Run SMC-ABCRF
# smcabcrf_test_3(
#     target = target,
#     model = model,
#     perturb = perturb,
#     parameters_initial = parameters_initial,
#     nIter = 7, # Number of iterations
#     nParticles = rep(1000, 7), # Number of particles for each iteration
#     # ntree = 2000,
#     parallel = T,
#     parameters_truth = data.frame(theta = theta, beta = beta)
# )

# =========================================================Run SMC-DRF
drf_output <- smcdrf(
    target = target,
    model = model,
    n_samples_per_parameter_set = 1,
    nNoise = nNoise,
    perturb = perturb,
    parameters_initial = parameters_initial,
    nIter = 7, # Number of iterations
    nParticles = rep(1000, 7), # Number of particles for each iteration
    # ntree = 2000,
    parallel = T
)

# filename <- "ABCSMC_DRF_output.rda"
# save(drf_output, file = filename)
# =========================================================Run SMC-ABCRF
rf_output <- smcabcrf(
    target = target,
    model = model,
    n_samples_per_parameter_set = 1,
    nNoise = nNoise,
    perturb = perturb,
    parameters_initial = parameters_initial,
    nIter = 7, # Number of iterations
    nParticles = rep(N, 7), # Number of particles for each iteration
    # ntree = 2000,
    parallel = T
)
# filename <- "ABCSMC_RF_output.rda"
# save(rf_output, file = filename)
# =========================================Plot SMC-ABCRF for Parameters
# plotting_smcrf(
#     parameters_truth = parameters_truth,
#     parameters_initial = parameters_initial,
#     parameters_id = colnames(parameters_initial),
#     outputdata = drf_output
# )
# plotting_smcrf(
#     parameters_truth = parameters_truth,
#     parameters_initial = parameters_initial,
#     parameters_id = colnames(parameters_initial),
#     outputdata = rf_output
# )
# ==============================================Plot SMC-ABCRF for Stats
plotting_smcrf(
    parameters_id = names(target),
    outputdata = drf_output,
    Plot_stats = TRUE
)
# plotting_smcrf(
#     parameters_id = names(target),
#     outputdata = rf_output,
#     Plot_stats = TRUE
# )

plot_joint(
    para_id=c("lambda","mu"), 
    outputdata=drf_output
)
plot_joint(
    para_id=c("lambda","mu"), 
    outputdata=rf_output
)