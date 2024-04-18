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
library(GillespieSSA)
# a<-1
# b<-1
# parms <- c(a=1,b=1)
# x0 <- c(X=1,Y=0.5)
# a <- c("a*X","X*Y","b*X*Y","Y")
# nu <- matrix(c(+1,-1,0,0,0,0,+1,-1), nrow = 2, byrow = TRUE)
# method <- ssa.etl()
# out <- ssa(x0,a,nu,parms,tf=15,method=method,
#             verbose = TRUE,consoleInterval=1,
#             censusInterval=1,
#             # maxWallTime = 30,
#             simName="Lotka predator-prey model")
# status <- out$stats$terminationStatus
# title <- paste0('ssa_', out$args$method$name, 's.png')
# png(title, res = 150, width = 30, height = 15, units = "in",)
# ssa.plot(out)
# dev.off()

nSSA<-0
nParam<-0
nSimulations<-0


library(deSolve)
LVmatrix <- function (Time, State, parms) {
    with(as.list(c(State, parms)), {
        dx = x*(alpha - beta*y)
        dy = -y*(gamma - delta*x)
        return(list(c(dx, dy)))
    })
}
parms <- c(alpha = 4, beta = 1, gamma = 1, delta = -1)
State <- c(x = 1, y = 0.5)
Time <- seq(0, 15, by = 0.01)
out <- ode(func = LVmatrix, y = State, parms = parms, times = Time)
matplot(out[,-1], type = "l", xlab = "time", ylab = "population")
legend("topright", c("y", "x"), lty = c(1,2), col = c(1,2), box.lwd = 0)
ylim(c(0,4))
# points(index,X, pch=19, col="red")
# points(index,Y, pch=19, col="#2c4081")



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
        # print('parallel')
        stats <- pblapply(
            cl = cl, X = 1:nrow(parameters),
            FUN = function(i) LV_model(parameters$a[i], parameters$b[i])
        )
        stopCluster(cl)
        # print(stats)
        stats <- rbindlist(stats)
        class(stats) <- "data.frame"
    } else {
        stats <- c()
        # print('=====')
        # print(parameters)
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
    # print(data)
    return(data)
}
LV_model <- function(a,b,parallel=FALSE) {
    #   This simulates a predator–prey process
    #   Input: 
    nRuns <- 1 
    library(deSolve)
    LVmatrix <- function (Time, State, parms) {
            with(as.list(c(State, parms)), {
                dx = x*(alpha - beta*y)
                dy = -y*(gamma - delta*x)
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
    # print('=====================')
    # print(parms) 
    for (i in 1:nRuns){
        out <- as.data.frame(ode(func = LVmatrix, y = State, parms = parms, times = Time)) 
        tsim <- c(1.1, 2.4, 3.9, 5.6, 7.5, 9.6, 11.9, 14.4)
        index <- c(111,241,391,561,751,961,1191,1441)
        X <- out$x[index]
        Y <- out$y[index]
        if(anyNA(X)){
            Xend <- as.double(tail(out,n=1)['x'])
            X[is.na(X)] <- Xend
            Yend <- as.double(tail(out,n=1)['y'])
            Y[is.na(Y)] <- Yend
        }
        # X <- X + rnorm(length(X),mean=0,sd=0.5)
        # Y <- Y + rnorm(length(Y),mean=0,sd=0.5)
        Xsim <- Xsim + X/nRuns
        Ysim <- Ysim + Y/nRuns
    }
    stats <- data.frame(matrix(c(a, b, Xsim, Ysim), nrow = 1))
    colnames(stats) <- c("a", "b", paste0("X_", 1:8), paste0("Y_", 1:8))
    # print(stats)
    return(stats)
}
# =====================================================Target statistics
parameters_truth <- data.frame(
    a = 1,
    b = 1
)
# statistics_target <- model(parameters = parameters_truth, parallel = FALSE)[-c(1:ncol(parameters_truth))]
Xt <- c(1.87, 0.65, 0.22, 0.31, 1.64, 1.15, 0.24, 2.91)
Yt <- c(0.49, 2.62, 1.54, 0.02, 1.14, 1.68, 1.07, 0.88)
statistics_target <- data.frame(matrix(c(Xt, Yt), nrow = 1))
colnames(statistics_target) <- c(paste0("X_", 1:8), paste0("Y_", 1:8))
# ======================================Model for parameter perturbation
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters, after perturbation
perturb <- function(parameters) {
    parameters[[1]] <- parameters[[1]] + runif(nrow(parameters), min=-0.1,max=0.1)
    parameters[[2]] <- parameters[[2]] + runif(nrow(parameters), min=-0.1,max=0.1)
    return(parameters)
}
# ======================================Define ranges for the parameters
range <- data.frame(
    parameter = c("a", "b"),
    min = c(-10, -10),
    max = c(10, 10)
)
# ========================================Initial guesses for parameters
# ====================================(sampled from prior distributions)
a <- runif(100000, -10, 10)
b <- runif(100000, -10, 10)
parameters_initial <- data.frame(
    a = a,
    b = b
)
# ====================================Labels for parameters in the plots
parameters_labels <- data.frame(
    parameter = c("a", "b"),
    label = c(deparse(expression(a)), deparse(expression(b)))
)
# ========================================SMC-RF for multiple parameters
#---Run SMC-RF for multiple parameters
smcrf_results_multi_param_test <- smcrf(
    method = "smcrf-multi-param-test",
    statistics_target = statistics_target,
    parameters_initial = parameters_initial,
    model = model,
    perturb = perturb,
    range = range,
    nParticles = rep(2000, 4),
    parallel = TRUE
)
# smcrf_results_multi_param <- smcrf(
#     method = "smcrf-multi-param",
#     statistics_target = statistics_target,
#     parameters_initial = parameters_initial,
#     model = model,
#     perturb = perturb,
#     range = range,
#     nParticles = rep(2000, 5),
#     parallel = TRUE
# )
# save(smcrf_results_multi_param, file='smcrf_results_multi_param.rda')
save(smcrf_results_multi_param_test, file='smcrf_results_multi_param_test.rda')
load("smcrf_results_multi_param.rda")
#---Plot marginal distributions
plot_smcrf_marginal(
    smcrf_results = smcrf_results_multi_param_test,
    parameters_truth = parameters_truth,
    parameters_labels = parameters_labels,
)
#---Plot posterior marginal distributions against other methods
plots <- plot_compare_marginal(
    # plots = plots,
    abc_results = smcrf_results_multi_param_test,
    parameters_labels = parameters_labels,
    parameters_truth = parameters_truth,
    plot_statistics = FALSE,
    plot_hist=TRUE,
    plot_prior=TRUE
)
# ========================================DRF
#---Run DRF
drf_results_test <- smcrf(
    method = "smcrf-multi-param-test",
    statistics_target = statistics_target,
    parameters_initial = parameters_initial,
    model = model,
    perturb = perturb,
    range = range,
    nParticles = rep(8000, 1),
    parallel = TRUE
)
save(drf_results_test, file='drf_results_test.rda')
# save(drf_results, file='drf_results.rda')
load("drf_results_test.rda")
#---Plot posterior marginal distributions against other methods
plots <- plot_compare_marginal(
    plots = plots,
    abc_results = drf_results_test,
    parameters_labels = parameters_labels,
    parameters_truth = parameters_truth,
    plot_statistics = FALSE,
    plot_hist=TRUE,
)
# ===============================================================ABC-SMC
# #---Run ABC-SMC
# abc_smc_results <- abc_smc(
#     statistics_target = statistics_target,
#     model = model,
#     parameters_labels = parameters_labels,
#     prior_distributions = list(c("unif", -10, 10),c("unif", -10,10)), 
#     nParticles = 1000, method = "Beaumont", 
#     tolerance = c(30.0, 16.0, 6.0, 5.0, 4.3),
#     progress_bar = TRUE,
#     # dist_weights = rep(1/ncol(statistics_target), ncol(statistics_target)),
#     # n_cluster = 10,use_seed=TRUE
# )
# save(abc_smc_results, file='abc_smc_results.rda')
# Simulation count:  ~20min 55536
load('abc_smc_results.rda')
#---Plot posterior marginal distributions against other methods
plots <- plot_compare_marginal(
    plots = plots,
    abc_results = abc_smc_results,
    parameters_labels = parameters_labels,
    parameters_truth = parameters_truth,
    plot_statistics = FALSE,
    plot_hist=TRUE,
)