# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Zijin - Macbook
R_workplace <- "/Users/xiangzijin/Documents/ABC_SMCRF/0321_sfs"
R_libPaths <- ""
R_libPaths_extra <- "/Users/xiangzijin/SMC-RF/R"
# =======================================SET UP FOLDER PATHS & LIBRARIES
# devtools::install_github("rdinnager/slimr", force = TRUE)
# slim_setup()
# usethis::edit_r_environ() SLIM_PATH="/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/slimr"
# Sys.setenv(SLIM_PATH = "/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/slimr")
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
                # mut_rate <- 1e-6
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

SFS_model <- function(mut_rate = c(1e-5), generation_num = 100) {
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
    stats <- data.frame(matrix(c(mut_rate, mutation_count, sfs[1:20]), nrow = 1))
    colnames(stats) <- c("theta", "Mutation_count_S", paste0("SFS_", 1:20))
    return(stats)
}
# =====================================================Target statistics
theta <- runif(1, 1e-6, 1e-4)
parameters_ground_truth <- data.frame(
    theta = theta
)
set.seed(1)
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
    min = c(1e-6),
    max = c(1e-4)
)
# ========================================Initial guesses for parameters
# ====================================(sampled from prior distributions)
theta <- runif(100000, 1e-6, 1e-4)
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
    abc_results = abcrf_results,
    parameters_truth = parameters_ground_truth,
    parameters_labels = parameters_labels,
    plot_statistics = TRUE
)

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
