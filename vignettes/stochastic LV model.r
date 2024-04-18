# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zhihan - Macbook
R_workplace <- "/Users/lexie/Documents/DNA/SMC-RF/vignettes/LV model"
R_libPaths <- ""
R_libPaths_extra <- "/Users/lexie/Documents/DNA/SMC-RF/R"
# R_libPaths_extra <- "/Users/lexie/Documents/DNA/EasyABC-master/R"
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


# ======================================LV
# library(GillespieSSA)
# c1<-10
# c2<-0.01
# c3<-10
# parms <- c(c1=10, c2=0.01, c3=10)
# # parms <- c(c1=5.763955,c2=  0.03841533,c3= 14.46528194  )
# x0 <- c(X=1000,Y=1000)
# a <- c("c1*X","c2*X*Y","c3*Y")
# nu <- matrix(c(+1, -1, 0, 0, +1, -1), nrow = 2, byrow = TRUE)
# method <- ssa.btl()
# out <- ssa(x0,a,nu,parms,tf=19,method=method,
#             verbose = TRUE,consoleInterval=1,
#             censusInterval=1,
#             maxWallTime = 30,
#             simName="Lotka predator-prey model")
# # status <- out$stats$terminationStatus
# # title <- paste0('ssa_', out$args$method$name, 's.png')
# # png(title, res = 150, width = 30, height = 15, units = "in",)
# ssa.plot(out)
# # dev.off()


library(GillespieSSA2)
sim_name <- "Lotka Predator-Prey model"
params <- c(c1 = 10, c2 = .01, c3 = 10)
final_time <- 2
initial_state <- c(Y1 = 1000, Y2 = 1000)
reactions <- list(
  reaction("c1 * Y1",      c(Y1 = +1)),
  reaction("c2 * Y1 * Y2", c(Y1 = -1, Y2 = +1)),
  reaction("c3 * Y2",      c(Y2 = -1))
)
out <- GillespieSSA2::ssa(
  initial_state = initial_state,
  reactions = reactions,
  params = params,
  final_time = final_time,
  method = ssa_exact(),
  census_interval = .001,
  verbose = TRUE,
  sim_name = sim_name
) 
plot_ssa(out)


nSSA<-0
nParam<-0
nSimulations<-0

# =====================================Model for the predator–prey process
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters, parallel = FALSE) {
    nNoise <- 0
    if (exists("nSimulations")) {
        nSimulations <<- nSimulations + nrow(parameters)
        print(nSimulations)
        }
    # if (exists("nSimulations")) nSimulations <<- nSimulations + nrow(parameters)
    #   Make simulations & compute summary statistics (population sizes at each time point)
    if (parallel) {
        library(parallel)
        library(pbapply)
        library(data.table)
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cl, varlist = c("LV_model"))
        print('parallel')
        stats <- pblapply(
            cl = cl, X = 1:nrow(parameters),
            FUN = function(i) LV_model(parameters$c1[i], parameters$c2[i], parameters$c3[i])
        )
        stopCluster(cl)
        print(stats)
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
    # print(data)
    return(data)
}
LV_model <- function(c1, c2, c3,parallel=FALSE) {
    #   This simulates a predator–prey process
    #   Input: 
    nRuns <- 3
    library(GillespieSSA)          
    parms <- c(c1=c1, c2=c2, c3=c3)
    x0 <- c(X=1000,Y=1000)
    a <- c("c1*X", "c2*X*Y", "c3*Y")
    nu <- matrix(c(+1, -1, 0, 0, +1, -1), nrow = 2, byrow = TRUE)
    Xsim <- numeric(length = 19)
    Ysim <- numeric(length = 19)
    print('=====================')
    print(parms)
    # print('start')
    # if (exists("nParam")) {
    #     nParam <<- nParam + 1
    # }  
    for (i in 1:nRuns){
        #=================================
        if (exists("nSSA")) {
                nSSA <<- nSSA + 1
                print(nSSA)
            } 
        out <- ssa(x0,a,nu,parms,tf=19,
            method=ssa.btl(),
            verbose = TRUE, 
            consoleInterval=1, censusInterval = 1,
            maxWallTime = 30,
            simName="Predator-prey model"
        )
        X <- out$data[-c(1),"X"]
        Y <- out$data[-c(1),"Y"]
        ll <- length(X)
        if (ll < 19){
            X_final <-X[ll]
            Y_final <-Y[ll]
            X <- c(X, rep(X_final, 19 - ll))
            Y <- c(Y, rep(Y_final, 19 - ll))
        }
        Xsim <- Xsim + X/nRuns
        Ysim <- Ysim + Y/nRuns
    }
    stats <- data.frame(matrix(c(c1, c2, c3, Xsim, Ysim), nrow = 1))
    colnames(stats) <- c("c1", "c2", "c3", paste0("X_", 1:19), paste0("Y_", 1:19))
    print(stats)
    return(stats)
}
# =====================================================Target statistics
parameters_truth <- data.frame(
    c1 = 10,
    c2 = 0.01,
    c3 = 10
)
# statistics_target <- model(parameters = parameters_truth, parallel = FALSE)[-c(1:ncol(parameters_truth))]
Xt <- c(1026.6, 1007.9, 998.5, 876.6, 1064.7, 815.7, 1047.0, 780.0, 1156.1, 631.1, 1724.0, 564.0, 1424.8, 377.2, 2028.0, 812.1, 903.0, 1575.3, 466.8)
Yt <- c(1084.7, 1053.9, 1050.9, 894.4, 1258.5, 640.8, 1686.8, 529.8, 1937.6, 449.0, 1177.3, 560.1, 1371.0, 815.5, 637.6, 1338.8, 287.2, 1631.3, 914.6)
statistics_target <- data.frame(matrix(c(Xt, Yt), nrow = 1))
colnames(statistics_target) <- c(paste0("X_", 1:19), paste0("Y_", 1:19))
# ======================================Model for parameter perturbation
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters, after perturbation
perturb <- function(parameters) {
    parameters[[1]] <- parameters[[1]] + runif(nrow(parameters), min = -1.0, max = 1.0)
    parameters[[2]] <- parameters[[2]] + runif(nrow(parameters), min = -0.0025, max = 0.0025)
    parameters[[3]] <- parameters[[3]] + runif(nrow(parameters), min = -1.0, max = 1.0)
    return(parameters)
}
# ======================================Define ranges for the parameters
range <- data.frame(
    parameter = c("c1", "c2", "c3"),
    min = c(0, 0, 0),
    max = c(28, 0.04, 28)
)
# ========================================Initial guesses for parameters
# ====================================(sampled from prior distributions)
c1 <- runif(100000, 0, 28)
c2 <- runif(100000, 0, 0.04)
c3 <- runif(100000, 0, 28)
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
# # ========================================SMC-RF for multiple parameters
# #---Run SMC-RF for multiple parameters
# smcrf_results_multi_param <- smcrf(
#     method = "smcrf-multi-param",
#     statistics_target = statistics_target,
#     parameters_initial = parameters_initial,
#     model = model,
#     perturb = perturb,
#     range = range,
#     nParticles = rep(10, 5),
#     parallel = TRUE
# )
# #---Plot marginal distributions
# plot_smcrf_marginal(
#     smcrf_results = smcrf_results_multi_param,
#     parameters_truth = parameters_truth,
#     parameters_labels = parameters_labels
# )
# #---Plot joint distributions
# for (i in 1:length(parameters_labels$parameter)){
#     plot_smcrf_joint(
#         smcrf_results = smcrf_results_multi_param,
#         parameters_labels = parameters_labels[-c(i),]
#     )
# }
# #---Plot posterior joint distributions against other methods
# plots_joint <- list()
# for (i in 1:length(parameters_labels$parameter)) {
#     plots_joint[[i]] <- plot_compare_joint(
#         abc_results = smcrf_results_multi_param,
#         parameters_labels = parameters_labels[-c(i),]
#     )
# }
# #---Plot posterior marginal distributions against other methods
# plots <- plot_compare_marginal(
#     abc_results = smcrf_results_multi_param,
#     parameters_labels = parameters_labels,
#     parameters_truth = parameters_truth,
#     plot_statistics = FALSE
# )
# # ========================================DRF
# #---Run DRF
# abcrf_results <- smcrf(
#     method = "smcrf-multi-param",
#     statistics_target = statistics_target,
#     parameters_initial = parameters_initial,
#     model = model,
#     perturb = perturb,
#     range = range,
#     nParticles = rep(5000, 1),
#     parallel = TRUE
# )
# #---Plot posterior joint distributions against other methods
# for (i in 1:length(parameters_labels$parameter)) {
#     plots_joint[[i]] <- plot_compare_joint(
#         plots = plots_joint[[i]],
#         abc_results = abcrf_results,
#         parameters_labels = parameters_labels[-c(i),]
#     )
# }

# # =========================================================ABC-rejection
# #---Run ABC-rejection
source('/Users/lexie/Documents/DNA/SMC-RF/vignettes/easyabc-r/EasyABC-internal.R')
abc_rej_results <- abc_rejection(
    statistics_target = statistics_target,
    model = model,
    parameters_labels = parameters_labels,
    prior_distributions = list(c("unif", 0, 28),c("unif", 0, 0.04),c("unif", 0, 28)), 
    nParticles = 10, tolerance_quantile = 0.1                                                , progress_bar = TRUE,
    # n_cluster=5, use_seed=TRUE,
    verbose=TRUE
)
# #---Plot posterior marginal distributions against other methods
# plots <- plot_compare_marginal(
#     # plots = plots,
#     abc_results = abc_rej_results,
#     parameters_labels = parameters_labels,
#     parameters_truth = parameters_truth,
#     plot_statistics = FALSE
# )
# #---Plot posterior joint distributions against other methods
# plots_joint <- list()
# for (i in 1:length(parameters_labels$parameter)) {
#     plots_joint[[i]] <- plot_compare_joint(
#         # plots = plots_joint[[i]],
#         abc_results = abc_rej_results,
#         parameters_labels = parameters_labels[-c(i),]
#     )
# }
# # ===============================================================ABC-SMC
# parameters_test <- do.call(rbind, replicate(50, parameters_truth, simplify = FALSE))
# statistics_test <- model(parameters = parameters_test,parallel=TRUE)[-c(1:ncol(parameters_test))]
# distance_matrix <- as.matrix(dist(statistics_test, method = "euclidean"))^2
# tolerance_min <- sum(distance_matrix) / (nrow(distance_matrix) * (ncol(distance_matrix) - 1))
# 14322574905 24223775063
# #---Run ABC-SMC
# abc_smc_results <- abc_smc(
#     statistics_target = statistics_target,
#     model = model,
#     parameters_labels = parameters_labels,
#     prior_distributions = list(c("unif", 0, 28),c("unif", 0, 0.04),c("unif", 0, 28)), 
#     nParticles = 100, method = "Beaumont", 
#     tolerance = c(4000^2, 3500^2, 3200^2), #, 3000^2, 2900^2
#     # tolerance = c(4000, 2900, 2000, 1900, 1800),
#     progress_bar = TRUE,
#     # dist_weights = rep(1/ncol(statistics_target), ncol(statistics_target)),
#     # n_cluster = 10,use_seed=TRUE
# )
# # # ~40-50min
# #---Plot posterior marginal distributions against other methods
# plots <- plot_compare_marginal(
#     # plots = plots,
#     abc_results = abc_smc_results,
#     parameters_labels = parameters_labels,
#     parameters_truth = parameters_truth,
#     plot_statistics = FALSE,
#     plot_hist=TRUE
# )
# #---Plot posterior joint distributions against other methods
# plots_joint <- list()
# for (i in 1:length(parameters_labels$parameter)) {
#     plots_joint[[i]] <- plot_compare_joint(
#         # plots = plots_joint[[i]],
#         abc_results = abc_smc_results,
#         parameters_labels = parameters_labels[-c(i),],
#     )
# }
