# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/0305_test/afs"
R_libPaths <- ""
R_libPaths_extra <- "/Users/xiangzijin/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zhihan - Macbook
# R_workplace <- "/Users/lexie/Documents/DNA/SMC-RF/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/lexie/Documents/DNA/SMC-RF/R"
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



set.seed(1)



# =========================Model for the Allele Frequency Spectrum (AFS)
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters) {
    library(parallel)
    library(pbapply)
    library(data.table)
    nSamples <- 100
    nNoise <- 0
    #   Make simulations & compute summary statistics (allele count)
    cl <- makePSOCKcluster(detectCores() - 1)
    clusterExport(cl, varlist = c("AFS_model"))
    stats <- pblapply(
        cl = cl, X = 1:nrow(parameters),
        FUN = function(i) {
            AFS_model(theta = parameters$theta[i], beta = 0, model_type = 1, n = nSamples)
        }
    )
    stopCluster(cl)
    stats <- rbindlist(stats)
    class(stats) <- "data.frame"
    #   Add noise statistics
    noise <- matrix(runif(nrow(parameters) * nNoise), nrow(parameters), nNoise)
    #   Add column names
    data <- data.frame(cbind(stats, noise))
    if (nNoise > 0) {
        colnames(data) <- c(
            colnames(stats),
            paste0("noise_", c(1:nNoise))
        )
    }
    return(data)
}
AFS_model <- function(theta, beta, model_type, n) {
    #   This simulates haplotype and site information
    #   Model_type = 1 for regular coalescent (beta = 0)
    #              = 2 for exponential growth coalescent model
    #              = 3 for pure death process with rate 1
    #---Set up arrays and initialisation
    hist <- list() # list of sets of common ancestor labels
    muts <- list() # list of mutation lists
    muts[[2 * n - 1]] <- 0
    mutindx <- 0 # label of last mutation to arise
    #---Generate coalescence times and related vectors
    ctimes <- rep(0, n) # coalescence times set to 0
    stimes <- rep(0, n) # jump times set to 0
    if (model_type != 3) {
        #   Generate ordinary coalescence times
        for (j in n:2) {
            ctimes[j] <- rexp(1, rate = j * (j - 1) / 2)
        }
        if (model_type == 2) {
            #   Calculate times from GT eqn (2.7)
            #   First, compute standard coalescent jump times
            stimes[n] <- ctimes[n]
            for (l in seq(n - 1, 2)) {
                stimes[l] <- ctimes[l] + stimes[l + 1]
            }
            #   Now compute exp growth jump times from GT (2.7)
            for (j in seq(2, n)) {
                stimes[j] <- log(1 + beta * stimes[j]) / beta
            }
            #   Then compute exp growth coal times by differencing
            ctimes[n] <- stimes[n]
            for (l in seq(n - 1, 2)) {
                ctimes[l] <- stimes[l] - stimes[l + 1]
            }
        }
    } else {
        for (j in n:2) {
            ctimes[j] <- rexp(1, rate = j)
        }
    }
    #---Simulate tree topology
    hist[[1]] <- 1:n
    for (j in n:2) { # number of things to choose from
        current <- hist[[n + 1 - j]] # current ancestor labels
        indx <- sample(1:j, 2, replace = F) # which two to merge
        current[indx[1]] <- 2 * n + 1 - j
        current[indx[2]] <- 2 * n + 1 - j
        hist[[n + 2 - j]] <- unique(current)
    }
    #---Pour mutations down the tree
    for (j in n:2) {
        parent <- setdiff(hist[[j]], hist[[j - 1]]) # parent label
        children <- setdiff(hist[[j - 1]], hist[[j]]) # children's labels
        for (l in children) {
            s <- rpois(1, theta * ctimes[n + 2 - j] / 2)
            if (s == 0) {
                muts[[l]] <- muts[[parent]]
            }
            if (s > 0) {
                mutstoadd <- seq(mutindx + 1, mutindx + s)
                mutindx <- mutindx + s # new last mutation label
                muts[[l]] <- union(muts[[parent]], mutstoadd)
            }
        }
        notmerged <- setdiff(hist[[j - 1]], children) # Set of labels w/ no merges at this level
        for (l in notmerged) {
            s <- rpois(1, theta * ctimes[n + 2 - j] / 2) # generate mutations on branch label l
            if (s > 0) {
                mutstoadd <- seq(mutindx + 1, mutindx + s)
                mutindx <- mutindx + s
                muts[[l]] <- union(muts[[l]], mutstoadd)
            }
        }
    }
    #---Find the haplotype distribution
    alleles <- table(sapply(muts[1:n], paste, collapse = " "))
    dimnames(alleles) <- NULL
    #---Find allele frequency spectrum
    afs <- tabulate(alleles, nbins = n)
    nalleles <- sum(afs)
    ncnt <- sum(1:n * afs)
    #---Find site frequency spectrum
    b1 <- unlist(muts[1:n])
    b2 <- tabulate(b1)
    sfs <- tabulate(b2, nbins = n - 1)
    sval <- sum(sfs)
    #---Combine table of statistics:
    #---number of alleles, the homozygoity statistic and the first sqrt(n) frequencies in the ESF
    ss <- 0
    for (i in 1:n) ss <- ss + afs[i] * (i / n)^2
    lvec <- floor(sqrt(n))
    if (model_type == 1) {
        stats <- data.frame(matrix(c(theta, nalleles), nrow = 1))
        colnames(stats) <- c("theta", "K")
    } else {
        stats <- data.frame(matrix(c(theta, beta, nalleles), nrow = 1))
        colnames(stats) <- c("theta", "beta", "K")
    }
    return(stats)
}
# =====================================================Target statistics
theta <- runif(1, 1, 20)
parameters_truth <- data.frame(
    theta = theta
)
statistics_target <- model(parameters = parameters_truth)[-c(1:ncol(parameters_truth))]
# ======================================Model for parameter perturbation
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters, after perturbation
perturb <- function(parameters) {
    for (i in 1:ncol(parameters)) parameters[[i]] <- parameters[[i]] + runif(nrow(parameters), min = -0.5, max = 0.5)
    return(parameters)
}
# ======================================Define ranges for the parameters
range <- data.frame(
    parameter = c("theta"),
    min = c(1),
    max = c(20)
)
# ========================================Initial guesses for parameters
# ====================================(sampled from prior distributions)
theta <- runif(1000, 1, 20)
parameters_initial <- data.frame(
    theta = theta
)
# ====================================Labels for parameters in the plots
parameters_labels <- data.frame(
    parameter = c("theta"),
    label = c(deparse(expression(theta)))
)
# ==========================================SMC-RF for single parameters
#---Run SMC-RF for single parameters
smcrf_results_single_param <- smcrf(
    method = "smcrf-single-param",
    statistics_target = statistics_target,
    parameters_initial = parameters_initial,
    model = model,
    perturb = perturb,
    range = range,
    nParticles = rep(1000, 7),
    parallel = TRUE
)
#---Plot marginal distributions
plot_smcrf_marginal(
    smcrf_results = smcrf_results_single_param,
    parameters_labels = parameters_labels,
    plot_statistics = TRUE
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
    nParticles = rep(1000, 7),
    parallel = TRUE
)
#---Plot marginal distributions
plot_smcrf_marginal(
    smcrf_results = smcrf_results_multi_param,
    parameters_labels = parameters_labels,
    plot_statistics = TRUE
)

# =========================================ABC-REJ for single parameters
# ========================================Initial guesses for parameters
# ====================================(sampled from prior distributions)
theta <- runif(7000, 1, 20)
parameters_initial <- data.frame(
    theta = theta
)
#---Run SMC-RF for single parameters
abcrej_results <- abc_rejection(
    statistics_target = statistics_target,
    model = model,
    parameters_initial = parameters_initial
)
#---Plot marginal distributions
plot_abc_marginal(
    abc_results = abcrej_results,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels
)
