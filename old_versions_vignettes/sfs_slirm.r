# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/0326_sfs"
R_libPaths <- ""
R_libPaths_extra <- "/Users/xiangzijin/SMC-RF/R"
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macmini
# R_workplace <- "/Users/khanhngocdinh/Documents/Zijin/0322_sfs_experiments"
# R_libPaths <- ""
# R_libPaths_extra <- "/Users/khanhngocdinh/Documents/Zijin/SMC-RF/R"
# =======================================SET UP FOLDER PATHS & LIBRARIES
# devtools::install_github("rdinnager/slimr", force = TRUE)
# slim_setup()
# usethis::edit_r_environ() SLIM_PATH="/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/slimr"
# Sys.setenv(SLIM_PATH = "/Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/library/slimr")
.libPaths(R_libPaths)
library(ggplot2)
library(gridExtra)
library(grid)
library(invgamma)
setwd(R_libPaths_extra)
files_sources <- list.files(pattern = "\\.[rR]$")
sapply(files_sources, source)
setwd(R_workplace)

# ===========================Model for the Site Frequency Spectrum (SFS)
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters & statistics, each row contains statistics for one set of parameters:
#           first columns = input parameters
#           next columns = summary statistics
model <- function(parameters, parallel = TRUE) {
    nNoise <- 0
    if (exists("nSimulations")) nSimulations <<- nSimulations + nrow(parameters)
    #   Make simulations & compute summary statistics (allele count)
    if (parallel) {
        library(parallel)
        library(pbapply)
        library(data.table)
        cl <- makePSOCKcluster(detectCores() - 1)
        clusterExport(cl, varlist = c("SFS_model"))
        stats <- pblapply(
            cl = cl, X = 1:nrow(parameters),
            FUN = function(i) {
                mut_rate <<- parameters$theta[i]
                SFS_model(mut_rate = mut_rate, generation_num = 100)
            }
        )
        stopCluster(cl)
        stats <- rbindlist(stats)
        class(stats) <- "data.frame"
    } else {
        stats <- c()
        for (i in 1:nrow(parameters)) {
            mut_rate <<- parameters$theta[i]
            stats <- rbind(stats, SFS_model(mut_rate = mut_rate, generation_num = 100))
        }
    }
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

SFS_model <- function(mut_rate, generation_num = 100) {
    library(slimr)
    library(adegenet)
    generation_num <- 100
    mut_rate <- mut_rate
    script_1 <- slim_script(
        slim_block(initialize(), {
            mut_rate <- r_inline(mut_rate)
            initializeMutationRate(mut_rate)
            ## m1 mutation type: neutral
            initializeMutationType("m1", 0.5, "e", 0.0)
            ## Mean selection coefficient for Exponential Distribution
            ## M1 is the neutral mutation type

            ## g1 genomic element type: uses m1 for all mutations
            initializeGenomicElementType("g1", m1, 1.0)

            ## uniform chromosome of length 100 kb
            initializeGenomicElement(g1, 0, 99999)
            ## uniform recombination along the chromosome
            initializeRecombinationRate(0)
        }),
        slim_block(1, early(), {
            ## create a population of 33 individuals
            sim.addSubpop("p1", 100)
        }),
        slim_block(1, 100, late(), {
            # create an output
            r_output_full("out", do_every = 10)
        }),
        slim_block(100, late(), {
            sim.simulationFinished()
        })
    )
    # ===================================================RUN SLIM SCRIPT
    slim_result <- slim_run(script_1, show_output = TRUE, capture_output = TRUE)
    mutlist <- slim_extract_full(slim_result$output_data, type = "genomes")[slim_extract_full(slim_result$output_data, type = "genomes")["generation"] == generation_num, ]["mut_list"]
    mutlist_all <- unlist(mutlist)
    mutlist_all[mutlist_all == 0] <- max(mutlist_all) + 1
    sfs <- tabulate(tabulate(mutlist_all))
    mutation_count <- sum(sfs)
    # ====================================================GET STATISTICS
    # stats <- data.frame(matrix(c(mut_rate, mutation_count, sfs[1:20]), nrow = 1))
    stats <- data.frame(matrix(c(mut_rate, mutation_count), nrow = 1))
    # colnames(stats) <- c("theta", "Mutation_count_S", paste0("SFS_", 1:20))
    colnames(stats) <- c("theta", "Mutation_count_S")
    return(stats)
}
# =====================================================Target statistics
set.seed(1)
theta <- runif(1, 1e-5, 1e-3)
parameters_ground_truth <- data.frame(
    theta = theta
)
statistics_target <- model(parameters = parameters_ground_truth, parallel = FALSE)[-c(1:ncol(parameters_ground_truth))]
print(statistics_target)
# ======================================Model for parameter perturbation
#   Input:  data frame of parameters, each row is one set of parameters
#   Output: data frame of parameters, after perturbation
perturb <- function(parameters) {
    for (i in 1:ncol(parameters)) parameters[[i]] <- parameters[[i]] + runif(nrow(parameters), min = -1, max = 1)
    return(parameters)
}
# ======================================Define ranges for the parameters
range <- data.frame(
    parameter = c("theta"),
    min = c(1e-5),
    max = c(1e-3)
)
# ========================================Initial guesses for parameters
# ====================================(sampled from prior distributions)
theta <- runif(10000, 1e-5, 1e-3)
parameters_initial <- data.frame(
    theta = theta
)
# ====================================Labels for parameters in the plots
parameters_labels <- data.frame(
    parameter = c("theta"),
    label = c(deparse(expression(theta)))
)
# ================================================================ABC-RF
#---Run ABC-RF
abcrf_results <- smcrf(
    method = "smcrf-single-param",
    statistics_target = statistics_target,
    parameters_initial = parameters_initial,
    model = model,
    perturb = perturb,
    range = range,
    nParticles = rep(500, 1),
    parallel = TRUE
)
#---Plot posterior marginal distributions against other methods
plots <- plot_compare_marginal(
    plots = plots,
    abc_results = abcrf_results,
    parameters_truth = parameters_ground_truth,
    parameters_labels = parameters_labels,
    plot_statistics = TRUE
)

save(abcrf_results, file = "ABCRF_results.rda")
# ========================================
#   Plot the out-of-bag estimates (equivalent to cross-validation)
abcrf_results[["Iteration_1"]][["rf_model"]][["model.rf"]]$predictions
abcrf_results$Iteration_1$rf_model$model.rf$predictions
abcrf_results$Iteration_1$parameters$theta
png(paste0("NEUTRAL_abcrf_theta_out_of_bag.png"))
plot(abcrf_results$Iteration_1$parameters$theta,
    abcrf_results$Iteration_1$rf_model$model.rf$predictions,
    xlab = "True value",
    ylab = "Out-of-bag estimate"
) + abline(a = 0, b = 1, col = "red")
dev.off()
#   Can the error be lowered by increasing the number of trees?
library(abcrf)
oob_error <- err.regAbcrf(abcrf_results$Iteration_1$rf_model, training = abcrf_results$Iteration_1$reference, paral = T)
png(paste0("NEUTRAL_abcrf_theta_error_by_ntree.png"))
plot(oob_error[, "ntree"], oob_error[, "oob_mse"], type = "l", xlab = "Number of trees", ylab = "Out-of-bag MSE")
dev.off()
#   Variance Importance of each statistic in inferring gamma
png(paste0("NEUTRAL_abcrf_theta_variance_importance.png"), width = 1500, height = 800, res = 150)
n.var <- min(30, length(abcrf_results$Iteration_1$rf_model$model.rf$variable.importance))
imp <- abcrf_results$Iteration_1$rf_model$model.rf$variable.importance
names(imp) <- colnames(statistics_target)
ord <- rev(order(imp, decreasing = TRUE)[1:n.var])
xmin <- 0
xlim <- c(xmin, max(imp) + 1)
dotchart(imp[ord], pch = 19, xlab = "Variable Importance", ylab = "", xlim = xlim, main = NULL, bg = "white", cex = 0.7)
dev.off()
# =========================================================ABC-rejection
#---Run ABC-rejection
# print("here")
# abc_rej_results <- abc_rejection(
#     statistics_target = statistics_target,
#     model = model,
#     parameters_labels = parameters_labels,
#     prior_distributions = list(c("unif", 1e-5, 1e-3)),
#     nParticles = 300, tolerance_quantile = 0.1, progress_bar = TRUE
# )
# #---Plot posterior marginal distributions against other methods
# plots <- plot_compare_marginal(
#     # plots = plots,
#     abc_results = abc_rej_results,
#     parameters_labels = parameters_labels,
#     parameters_truth = parameters_ground_truth,
#     plot_statistics = TRUE
# )
# ========================================
# parameters <- data.frame(parameters_initial[1:100, ])
# colnames(parameters) <- colnames(parameters_initial)
# reference <- model(parameters = parameters)
# library(abcrf)
# model_rf <- regAbcrf(
#     formula = as.formula(paste0("theta", " ~ .")),
#     data = reference,
#     paral = TRUE
# )
# posterior_gamma_RF <- predict(
#     object = model_rf,
#     obs = statistics_target,
#     training = reference,
#     paral = TRUE,
#     rf.weights = T
# )
#   Variance Importance of each statistic in inferring gamma
# png(paste0("NEUTRAL_abcrf_test_variance_importance.png"), width = 1500, height = 800, res = 150)
# n.var <- min(30, length(model_rf$model.rf$variable.importance))
# imp <- model_rf$model.rf$variable.importance
# names(imp) <- colnames(statistics_target)
# ord <- rev(order(imp, decreasing = TRUE)[1:n.var])
# xmin <- 0
# xlim <- c(xmin, max(imp) + 1)
# dotchart(imp[ord], pch = 19, xlab = "Variable Importance", ylab = "", xlim = xlim, main = NULL, bg = "white", cex = 0.7)
# dev.off()
# imp
# groups <- c()
# for (i in 1:length(names(imp))) {
#     if (any(grepl(names(imp)[i], copy_number_stats, fixed = TRUE))) {
#         group_idx <- "copy_number_stats"
#     } else if (any(grepl(names(imp)[i], phylo_balance_stats, fixed = TRUE))) {
#         group_idx <- "phylo_balance_stats"
#     } else if (any(grepl(names(imp)[i], phylo_tip_stats, fixed = TRUE))) {
#         group_idx <- "phylo_tip_stats"
#     }
#     groups <- c(groups, group_idx)
# }

# group_colors <- c("copy_number_stats" = "#00BA38", "phylo_balance_stats" = "#F8766D", "phylo_tip_stats" = "#619CFF")

# imp
# ord <- rev(order(imp, decreasing = TRUE)[1:n.var])
# # colors <- group_colors[groups[ord]]
# xmin <- 0
# xlim <- c(xmin, max(imp) + 1)
# dotchart(imp[ord], pch = 19, xlab = "Variable Importance", ylab = "", xlim = xlim, main = NULL, bg = "white", cex = 0.7)
# # legend("bottomright", legend = names(group_colors), fill = group_colors, bty = "n", cex = 0.7)
# dev.off()
# # =====================================================Target statistics
# mut_rates <- runif(100, 1e-6, 2e-5)
# parameters_initial <- data.frame(
#     mut_rate = mut_rates
# )
# # parameters_ground_truth$mut_rate[9]

# set.seed(1)
# # statistics_target <- model(parameters = parameters_initial, parallel = FALSE)
# statistics_target <- model(parameters = parameters_initial, parallel = TRUE)
# print(statistics_target)


# library(slimr)
# library(adegenet)
# mut_rate <- 1e-5
# generation_num <- 100
# pop_size <- 100
# script_1 <- slim_script(
#     slim_block(initialize(), {
#         mut_rate <- r_inline(mut_rate)
#         initializeMutationRate(mut_rate)
#         ## m1 mutation type: neutral
#         initializeMutationType("m1", 0.5, "e", 0.0)
#         ## Mean selection coefficient for Exponential Distribution
#         ## M1 is the neutral mutation type

#         ## g1 genomic element type: uses m1 for all mutations
#         initializeGenomicElementType("g1", m1, 1.0)

#         ## uniform chromosome of length 100 kb
#         initializeGenomicElement(g1, 0, 99999)
#         ## uniform recombination along the chromosome
#         initializeRecombinationRate(0)
#     }),
#     slim_block(1, early(), {
#         ## create a population of 33 individuals
#         pop_size <- r_inline(pop_size)
#         sim.addSubpop("p1", pop_size)
#     }),
#     slim_block(1, 100, late(), {
#         # create an output
#         r_output_full("out", do_every = 10)
#     }),
#     slim_block(100, late(), {
#         sim.simulationFinished()
#     })
# )

# slim_result <- slim_run(script_1, show_output = TRUE, capture_output = TRUE)
# mutlist <- slim_extract_full(slim_result$output_data, type = "genomes")[slim_extract_full(slim_result$output_data, type = "genomes")["generation"] == generation_num, ]["mut_list"]
# View(mutlist)
# mutlist_all <- unlist(mutlist)
# tabulate(mutlist_all)
# tabulate(tabulate(mutlist_all))
# sfs <- sort(tabulate(mutlist_all), decreasing = TRUE)
#     mutation_count <- sum(sfs)
# table(mutlist_all)
# tail(order(tabulate(mutlist_all)), 50)
# sort(tabulate(mutlist_all),decreasing = TRUE)
# max_numbers

#  tabulate(c(-2, 0, 2, 3, 3, 5))
# tabulate(unlist(mutlist))
# mutlist_all
# order(tabulate(mutlist_all[1:911]), decreasing = TRUE)
# # mutlist <- (slim_extract_full(sr$output_data, type = "genomes")[slim_extract_full(sr$output_data, type = "genomes")["generation"] == 100, ]["mut_list"])


# max(mutlist_all)
# mutlist_all
# mutlist_all[mutlist_all == 0] <- max(mutlist_all) + 1
# max(mutlist_all)


# order(tabulate(unlist(a)), decreasing = TRUE)
# # we want to run SLiM using 5 local cores
# # please note on Linux and Windows you would like to use multicore as this is faster as it is a 'real' fork,
# # but multisession works on all systems
# plan(multisession, workers = 5)
# sr <- slim_run(script_2, parallel = TRUE, show_output = TRUE, capture_output = TRUE)
# View(sr$output_data)

# length(slim_extract_full(sr$output_data, type = "genomes")[slim_extract_full(sr$output_data, type = "genomes")["generation"] == 100, ]["mut_list"])

# slim_extract_full(sr$output_data, type = "mutations")[slim_extract_full(sr$output_data, type = "mutations")["generation"] == 100, ]["mut_id"]
# # gls <- slim_extract_genlight(sr, by = c("generation"))
# # gls

# # View(slim_extract_full(sr$output_data, type = "full_individual"))
# View(slim_extract_full(sr$output_data, type = "genomes")[slim_extract_full(sr$output_data, type = "genomes")["generation"] == 100, ])
# # slim_extract_full(sr$output_data, type = "genomes")
# # the genlight object for generation 10
# # gls$genlight[[10]]
# # a <- unlist(slim_extract_full(sr$output_data, type = "genomes")[slim_extract_full(sr$output_data, type = "genomes")["generation"] == 100, ]["mut_list"])
# # a[1]
# # View(a)
# a <- slim_extract_full(sr$output_data, type = "genomes")[slim_extract_full(sr$output_data, type = "genomes")["generation"] == 100, ]["mut_list"]
# (frequency(table(unlist(a)))))
# hist(data.frame(table(unlist(a))["Freq"]))
# hist(data.frame(table(unlist(a)))["Freq"])
# table(unlist(a))
# data.frame(table(unlist(a)))
# order(data.frame(table(unlist(a)))["Freq"], decreasing = TRUE)[1:10]
# data.frame(table(unlist(a)))[order(data.frame(table(unlist(a)))["Freq"], decreasing = TRUE)[1:10], ]
# # number of loci at generation 10
# # nLoc(gls$genlight[[1]])

# # nloci <- unlist(lapply(gls$genlight, nLoc))
# # barplot(nloci,
# # col = rainbow(length(gls$genlight)), names.arg = 1:10, xlab = "generation",
# # ylab = "# loci"
# # )
# ?order
# # initializeGenomicElementType()


# ########################


# # library(adegenet)
# # slim_script(
# #     slim_block(initialize(), {
# #         initializeMutationRate(1e-5)
# #         ## m1 mutation type: neutral
# #         initializeMutationType("m1", 0.5, "f", 0.0)
# #         ## g1 genomic element type: uses m1 for all mutations
# #         initializeGenomicElementType("g1", m1, 1.0)
# #         ## uniform chromosome of length 100 kb
# #         initializeGenomicElement(g1, 0, 99999)
# #         ## uniform recombination along the chromosome
# #         initializeRecombinationRate(0e-8)
# #     }),
# #     slim_block(1, early(), {
# #         ## create a population of 33 individuals
# #         sim.addSubpop("p1", 33)
# #     }),
# #     slim_block(1, 100, late(), {
# #         # create an output
# #         r_output(p1.genomes.output(), "p1", do_every = 10)
# #     }),
# #     slim_block(100, {
# #         r_output(p1.outputVCFSample(sampleSize = 10), name = "VCF")
# #         sim.simulationFinished()
# #     }),
# #     slim_block(100, late(), {
# #         sim.simulationFinished()
# #     })
# # ) -> script_2

# # library(future)
# # # we want to run SLiM using 5 local cores
# # # please note on Linux and Windows you would like to use multicore as this is faster as it is a 'real' fork,
# # # but multisession works on all systems
# # plan(multisession, workers = 5)
# # sr <- slim_run(script_2, parallel = TRUE, show_output = TRUE, capture_output = TRUE)
# # sr$output_data
# # gls <- slim_extract_genlight(sr, by = c("generation"))
# # gls


# # sr
# # slim_extract_output_data(sr)
# # # the genlight object for generation 10
# # gls$genlight[[1]]

# # # number of loci at generation 10
# # nLoc(gls$genlight[[1]])

# # nloci <- unlist(lapply(gls$genlight, nLoc))
# # barplot(nloci,
# #     col = rainbow(length(gls$genlight)), names.arg = 1:10, xlab = "generation",
# #     ylab = "# loci"
# # )

# # initializeGenomicElementType()
