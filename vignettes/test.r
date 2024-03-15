# R_workplace <- getwd()
# R_libPaths_extra <- "/Users/dinhngockhanh/Library/CloudStorage/GoogleDrive-knd2127@columbia.edu/My Drive/RESEARCH AND EVERYTHING/Projects/GITHUB/SMC-RF/EasyABC"
# setwd(R_libPaths_extra)
# files_sources <- list.files(pattern = "\\.[rR]$")
# sapply(files_sources, source)
# setwd(R_workplace)

.check_prior_test <- function(nb_parameters, prior_test) {
    m <- gregexpr("[xX]([0-9])+", prior_test)
    result <- regmatches(prior_test, m)
    for (i in 1:length(result[[1]])) {
        parameter_index <- substring(result[[1]][i], 2)
        # print(result)
        print(parameter_index)
        print(nb_parameters)
        print(parameter_index > nb_parameters)
        if (parameter_index > nb_parameters) {
            stop(paste("Parameter out of range: ", result[[1]][i], sep = ""))
        }
    }
    TRUE
}
print(.check_prior_test(9, "X1 > 2; X2 < 4"))

# library(EasyABC)

# set.seed(1)

# toy_model2 <- function(x) {
#     c(x[1] + x[2] + rnorm(1, 0, 0.1), x[1] * x[2] + rnorm(1, 0, 0.1))
# }
# toy_prior2 <- list(c("unif", 0, 1), c("normal", 1, 2))
# n <- 100
# sum_stat_obs2 <- c(1.5, 0.5)
# #-----------------------------------------------------------------------
# tolerance <- 0.2
# ABC_rej <- ABC_rejection(
#     model = toy_model2, prior = toy_prior2, summary_stat_target = sum_stat_obs2,
#     nb_simul = n,
#     tol = tolerance
# )
# #-----------------------------------------------------------------------
# ABC_Marjoram <- ABC_mcmc(
#     model = toy_model2, prior = toy_prior2, summary_stat_target = sum_stat_obs2,
#     n_rec = n,
#     method = "Marjoram"
# )
# #-----------------------------------------------------------------------
# tolerance <- c(1.5, 0.5)
# ABC_Beaumont <- ABC_sequential(
#     model = toy_model2, prior = toy_prior2, summary_stat_target = sum_stat_obs2,
#     nb_simul = n,
#     method = "Beaumont", tolerance_tab = tolerance
# )
