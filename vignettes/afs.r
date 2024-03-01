# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Khanh - Macbook
# R_workplace <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/vignettes"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/0224_test/afs"
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
# =============================================SET UP INITIAL PARAMETERS
N <- 10000
n_samples_per_parameter_set <- 100
nNoise <- 10 # Number of noise variables
# ======================================================================
set.seed(1)
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
        # cat("\n")
        # cat("Parent label: ",parent,"\n")
        # cat("mutation labels: ",muts[[parent]],"\n")
        children <- setdiff(hist[[j - 1]], hist[[j]]) # children's labels
        # cat("Child labels: ",children,"\n")
        for (l in children) {
            s <- rpois(1, theta * ctimes[n + 2 - j] / 2)
            # cat("label, number of mutations on branch: ",c(l,s),"\n")
            if (s == 0) {
                muts[[l]] <- muts[[parent]]
                # cat("mutation labels: ",muts[[l]],"\n")
            }
            if (s > 0) {
                # cat("mutindx: ",mutindx,"\n")
                mutstoadd <- seq(mutindx + 1, mutindx + s)
                # cat("Mutation labels to add: ",mutstoadd,"\n")
                mutindx <- mutindx + s # new last mutation label
                muts[[l]] <- union(muts[[parent]], mutstoadd)
                # cat("mutation labels: ",muts[[l]],"\n")
            }
        }
        notmerged <- setdiff(hist[[j - 1]], children) # Set of labels w/ no merges at this level
        for (l in notmerged) {
            s <- rpois(1, theta * ctimes[n + 2 - j] / 2) # generate mutations on branch label l
            # cat("label, number of mutations on branch: ",c(l,s),"\n")
            # if(s == 0) cat("mutation labels: ",muts[[l]],"\n")
            if (s > 0) {
                mutstoadd <- seq(mutindx + 1, mutindx + s)
                mutindx <- mutindx + s
                muts[[l]] <- union(muts[[l]], mutstoadd)
                # cat("mutation labels: ",muts[[l]],"\n")
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
    # stats <- data.frame(matrix(c(theta, beta, nalleles, ss, sval, afs[1:lvec]), nrow = 1))
    # colnames(stats) <- c("theta", "beta", "K", "homzy", "S", paste0("AFS_", 1:lvec))
    stats <- data.frame(matrix(c(theta, beta, nalleles), nrow = 1))
    colnames(stats) <- c("theta", "beta", "K")
    return(stats)
}
#---Model for the Allele Frequency Spectrum (AFS)
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
    thetas <- parameters$theta
    betas <- parameters$beta
    model_type <- 1
    cl <- makePSOCKcluster(detectCores() - 1)
    AFS_model <<- AFS_model
    n_samples_per_parameter_set <<- n_samples_per_parameter_set
    clusterExport(cl, varlist = c("AFS_model", "n_samples_per_parameter_set"))
    stats <- pblapply(cl = cl, X = 1:nrow(parameters), FUN = function(i) AFS_model(thetas[i], betas[i], model_type, n_samples_per_parameter_set))
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
# ======================================Model for parameter perturbation
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters, after perturbation
# perturb <- function(parameters) {
#     parameters$theta <- pmin(pmax(parameters$theta + runif(nrow(parameters), min = -1, max = 1), 1), 20)
#     parameters$beta <- pmin(pmax(parameters$beta + runif(nrow(parameters), min = -1, max = 1), 1), 20)
#     return(parameters)
# }
perturb <- function(parameter) {
    # parameter <- parameter + runif(1, min = -0.5, max = 0.5)
    parameter <- parameter + runif(1, min = -1, max = 1)
    # parameters$theta2 <- pmax(parameters$theta2 + runif(nrow(parameters), min = -0.5, max = 0.5), 0)
    return(parameter)
}

range <- data.frame(
    parameter = c("theta", "beta"),
    min = c(1, 0),
    max = c(20, 1)
)

# =====================================================Target statistics
target_theta <- runif(1, 1, 10)
target_beta <- 0
parameters_target <- data.frame(
    theta = target_theta,
    beta = target_beta
)
target <- model(
    parameters = parameters_target,
    n_samples_per_parameter_set = n_samples_per_parameter_set,
    nNoise = nNoise
)[-c(1:ncol(parameters_target))]
parameters_truth <- data.frame(
    theta = target_theta,
    beta = target_beta
)
cat(paste0("\n\n\n\n\nTRUE VALUE FOR THETA = ", target_theta, "\n"))
cat(paste0("STATISTICS:\n"))
print(target)
cat("\n\n\n\n\n")
# ========================================Initial guesses for parameters
# ====================================(sampled from prior distributions)
theta <- runif(N, 1, 20)
beta <- runif(N, 0, 1)
parameters_initial <- data.frame(
    theta = theta,
    beta = beta
)
# ========================================================Run SMC-ABCDRF
drf_output <- smcdrf(
    target = target,
    model = model,
    n_samples_per_parameter_set = n_samples_per_parameter_set,
    nNoise = nNoise,
    perturb = perturb,
    range = range,
    parameters_initial = parameters_initial,
    nIter = 7, # Number of iterations
    nParticles = rep(1000, 7), # Number of particles for each iteration
    # ntree = 2000,
    parallel = T,
)
filename <- "ABCSMC_DRF_output.rda"
save(drf_output, file = filename)
# =========================================================Run SMC-ABCRF
rf_output <- smcabcrf(
    target = target,
    model = model,
    n_samples_per_parameter_set = n_samples_per_parameter_set,
    nNoise = nNoise,
    perturb = perturb,
    range = range,
    parameters_initial = parameters_initial,
    nIter = 7, # Number of iterations
    nParticles = rep(1000, 7), # Number of particles for each iteration
    # ntree = 2000,
    parallel = T
)
filename <- "ABCSMC_RF_output.rda"
save(rf_output, file = filename)
# =========================================Plot SMC-ABCRF for Parameters
plotting_smcrf(
    parameters_truth = parameters_truth,
    parameters_initial = parameters_initial,
    parameters_id = colnames(parameters_initial),
    outputdata = drf_output
)
plotting_smcrf(
    parameters_truth = parameters_truth,
    parameters_initial = parameters_initial,
    parameters_id = colnames(parameters_initial),
    outputdata = rf_output
)
# ==============================================Plot SMC-ABCRF for Stats
plotting_smcrf(
    parameters_id = names(target),
    outputdata = drf_output,
    Plot_stats = TRUE
)
plotting_smcrf(
    parameters_id = names(target),
    outputdata = rf_output,
    Plot_stats = TRUE
)
