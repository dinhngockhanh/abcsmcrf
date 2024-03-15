model_dlp <- function(parameters) {
    library(parallel)
    library(pbapply)
    library(data.table)
    nNoise <- 0
    parameter_IDs <- colnames(parameters)
    #   Make simulations & compute summary statistics (allele count)
    cl <- makePSOCKcluster(detectCores() - 1)
    clusterExport(cl, varlist = c(
        "get_statistics_dlp", "func_stat", "assign_paras", "model_variables", "simulator_full_program", "get_each_clonal_CN_profiles",
        "list_parameters", "get_arm_CN_profiles", "cn_table", "get_each_statistics", "find_clonal_ancestry", "find_event_count",
        "get_statistics", "ground_truth_cn_data_sc", "ground_truth_cn_data_bulk", "cohort_distance", "cn_distance", "sample_distance"
    ))
    stats <- pblapply(
        cl = cl, X = 1:nrow(parameters),
        FUN = function(i) {
            get_statistics_dlp(
                parameter_IDs = parameter_IDs,
                simulation_parameters = data.frame(parameters[i, ]),
                model_variables = model_variables,
                ground_truth_cn_data_sc = ground_truth_cn_data_sc,
                ground_truth_cn_data_bulk = ground_truth_cn_data_bulk
            )
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
    } else {
        colnames(data) <- colnames(stats)
    }
    return(data)
}








get_statistics_dlp <- function(parameter_IDs,
                               simulation_parameters,
                               model_variables,
                               ground_truth_cn_data_sc,
                               ground_truth_cn_data_bulk) {
    #---------------------------------List of parameter IDs to be fitted
    n_simulations_sc <- 3
    n_simulations_bulk <- 4
    arm_level <- TRUE

    list_targets_library <- c(
        #---Bulk DNA: CN
        "data=bulk;target=genome;statistic=dist;variable=average_CN;metric=euclidean",
        "data=bulk;target=genome;statistic=mean;representative_CN=average_CN;variable=event_count;type=total;event=missegregation",
        "data=bulk;target=genome;statistic=var;representative_CN=average_CN;variable=event_count;type=total;event=missegregation",
        "data=bulk;target=genome;statistic=mean;representative_CN=average_CN;variable=event_count;type=total;event=chromosome-arm-missegregation",
        "data=bulk;target=genome;statistic=var;representative_CN=average_CN;variable=event_count;type=total;event=chromosome-arm-missegregation",
        #---Single-cell DNA: subclonal CN
        "data=sc;target=genome;statistic=dist;variable=clonal_CN;metric=euclidean",
        "data=sc;target=genome;statistic=mean;variable=shannon",
        "data=sc;target=genome;statistic=mean;variable=event_count;type=clonal;event=missegregation",
        "data=sc;target=genome;statistic=mean;variable=event_count;type=subclonal;event=missegregation",
        "data=sc;target=genome;statistic=mean;variable=event_count;type=clonal;event=chromosome-arm-missegregation",
        "data=sc;target=genome;statistic=mean;variable=event_count;type=subclonal;event=chromosome-arm-missegregation",
        "data=sc;target=genome;statistic=var;variable=shannon",
        "data=sc;target=genome;statistic=var;variable=event_count;type=clonal;event=missegregation",
        "data=sc;target=genome;statistic=var;variable=event_count;type=subclonal;event=missegregation",
        "data=sc;target=genome;statistic=var;variable=event_count;type=clonal;event=chromosome-arm-missegregation",
        "data=sc;target=genome;statistic=var;variable=event_count;type=subclonal;event=chromosome-arm-missegregation",
        #---Single-cell DNA: phylo stats for balance
        "data=sc;target=genome;statistic=mean;variable=stairs", # proportion of subtrees that are imbalanced
        "data=sc;target=genome;statistic=mean;variable=colless", # balance index of phylogeny tree
        "data=sc;target=genome;statistic=mean;variable=sackin", # balance index of phylogeny tree
        "data=sc;target=genome;statistic=mean;variable=B2", # balance index of phylogeny tree
        "data=sc;target=genome;statistic=mean;variable=maxDepth" # height of phylogeny tree
    )
    # colnames(simulation_parameters) <- parameter_IDs
    simulation_parameters <- as.matrix(simulation_parameters)
    # ===============GET NEW MODEL VARIABLES BY ASSIGNING NEW PARAMETERS

    model_variables <- assign_paras(model_variables, parameter_IDs, parameters = simulation_parameters)

    # ==========GET TABLE OF CHROMOSOME LENGTHS AND CENTROMERE LOCATIONS
    cn_table <- model_variables$cn_info
    cn_bin_length <- as.numeric(model_variables$general_variables$Value[which(model_variables$general_variables$Variable == "size_CN_block_DNA")])
    cn_table$Length <- cn_table$Bin_count * cn_bin_length
    cn_table$Centromere <- cn_table$Centromere_location * cn_bin_length

    #   Make single-cell simulations and get relevant statistics
    simulations_sc <- simulator_full_program(
        model = model_variables, model_prefix = "", n_simulations = n_simulations_sc,
        stage_final = 3,
        save_simulation = FALSE, report_progress = TRUE,
        lite_memory = TRUE,
        output_variables = c(
            "evolution_origin",
            "evolution_genotype_changes",
            "sample_clone_ID",
            "sample_genotype_unique",
            "sample_genotype_unique_profile",
            "phylogeny_clustering_truth"
        )
    )
    simulations_clonal_CN_sc <- get_each_clonal_CN_profiles(simulations_sc, arm_level, cn_table)
    list_targets_sc <- list_targets_library[grepl("data=sc", list_targets_library)]
    simulations_statistics_sc <- get_each_statistics(simulations_sc, simulations_clonal_CN_sc, list_targets_sc)
    #   Make bulk simulations and get relevant statistics
    simulations_bulk <- simulator_full_program(
        model = model_variables, model_prefix = "", n_simulations = n_simulations_bulk,
        stage_final = 2,
        save_simulation = FALSE, report_progress = TRUE,
        lite_memory = TRUE,
        output_variables = c(
            "evolution_origin",
            "evolution_genotype_changes",
            "sample_clone_ID",
            "sample_genotype_unique",
            "sample_genotype_unique_profile"
        )
    )
    simulations_clonal_CN_bulk <- get_each_clonal_CN_profiles(simulations_bulk, arm_level, cn_table, bulk = TRUE)
    list_targets_bulk <- list_targets_library[grepl("data=bulk", list_targets_library)]
    simulations_statistics_bulk <- get_each_statistics(simulations_bulk, simulations_clonal_CN_bulk, list_targets_bulk)

    list_stats <- get_statistics(
        simulations_statistics_sc = simulations_statistics_sc,
        simulations_statistics_bulk = simulations_statistics_bulk,
        cn_data_sc = ground_truth_cn_data_sc,
        cn_data_bulk = ground_truth_cn_data_bulk,
        list_targets = list_targets_library,
        arm_level = arm_level,
        cn_table = cn_table
    )
    colnames(simulation_parameters) <- parameter_IDs
    stats <- cbind(simulation_parameters, data.frame(list_stats$statistics))
    return(stats)
}
