compute_error <- function(results, ID_actual, ID_predicted) {
    library(Metrics)
    error <- rmse(results[[ID_actual]], results[[ID_predicted]])
    # library(ehaGoF)
    # error <- gofRRMSE(results[[ID_actual]], results[[ID_predicted]], dgt = 3)
    return(error)
}

#---Function to get arm-level CN profiles
#---from bin-level CN profiles
get_arm_CN_profiles <- function(clonal_CN_profiles) {
    clonal_CN_profiles_arm <- vector("list", length(clonal_CN_profiles))
    for (sample in 1:length(clonal_CN_profiles)) {
        clonal_CN_profiles_arm[[sample]] <- vector("list", length(clonal_CN_profiles[[sample]]))
        for (clone in 1:length(clonal_CN_profiles[[sample]])) {
            clone_CN_profile <- clonal_CN_profiles[[sample]][[clone]]
            clone_CN_profile_arm <- data.frame(matrix(ncol = (ncol(clone_CN_profile) + 1), nrow = 0))
            colnames(clone_CN_profile_arm) <- c(colnames(clone_CN_profile), "arm")
            for (chrom in cn_table$Chromosome) {
                for (arm in c("p", "q")) {
                    if (arm == "p") {
                        start <- 1
                        end <- cn_table[cn_table$Chromosome == chrom, ]$Centromere
                    } else if (arm == "q") {
                        start <- 1 + cn_table[cn_table$Chromosome == chrom, ]$Centromere
                        end <- cn_table[cn_table$Chromosome == chrom, ]$Length
                    }
                    locs <- which(
                        clone_CN_profile$chr == chrom & clone_CN_profile$start >= start & clone_CN_profile$end <= end
                    )
                    Min <- as.numeric(names(sort(summary(as.factor(clone_CN_profile[locs, ]$Min)), decreasing = T)[1]))
                    Maj <- as.numeric(names(sort(summary(as.factor(clone_CN_profile[locs, ]$Maj)), decreasing = T)[1]))
                    tmp <- clone_CN_profile[locs, ]
                    state <- Min + Maj
                    copy <- state
                    clone_CN_profile_arm[nrow(clone_CN_profile_arm) + 1, ] <- 0
                    clone_CN_profile_arm$chr[nrow(clone_CN_profile_arm)] <- chrom
                    clone_CN_profile_arm$start[nrow(clone_CN_profile_arm)] <- start
                    clone_CN_profile_arm$end[nrow(clone_CN_profile_arm)] <- end
                    clone_CN_profile_arm$Min[nrow(clone_CN_profile_arm)] <- Min
                    clone_CN_profile_arm$Maj[nrow(clone_CN_profile_arm)] <- Maj
                    clone_CN_profile_arm$state[nrow(clone_CN_profile_arm)] <- state
                    clone_CN_profile_arm$copy[nrow(clone_CN_profile_arm)] <- copy
                    clone_CN_profile_arm$arm[nrow(clone_CN_profile_arm)] <- arm
                }
            }
            clonal_CN_profiles_arm[[sample]][[clone]] <- clone_CN_profile_arm
        }
    }
    return(clonal_CN_profiles_arm)
}

get_each_clonal_CN_profiles <- function(simulations,
                                        arm_level = FALSE,
                                        cn_table = NULL,
                                        bulk = FALSE) {
    #-----------------------------------------Find bin-level CN profiles
    clonal_CN_profiles_all_sims <- list()
    for (i in 1:length(simulations)) {
        simulation <- simulations[[i]]
        clonal_CN_profiles <- simulation$sample$sample_genotype_unique_profile
        sample_genotype_unique <- simulation$sample$sample_genotype_unique
        clonal_CN_populations <- rep(0, length(clonal_CN_profiles))
        for (j in 1:length(clonal_CN_profiles)) {
            clonal_CN_populations[j] <- sum(simulation$sample$sample_clone_ID == sample_genotype_unique[j])
        }
        clonal_CN_profiles_all_sims[["variable=clonal_CN_profiles"]][[i]] <- clonal_CN_profiles
        clonal_CN_profiles_all_sims[["variable=clonal_CN_populations"]][[i]] <- clonal_CN_populations
    }
    #----------------------Convert to arm-level CN profiles if requested
    if (arm_level) {
        clonal_CN_profiles_all_sims[["variable=clonal_CN_profiles"]] <-
            get_arm_CN_profiles(clonal_CN_profiles_all_sims[["variable=clonal_CN_profiles"]])
    }
    #--------------------------------Convert to bulk format if necessary
    if (bulk) {
        clonal_CN_profiles_all_sims_new <- list()
        for (i in 1:length(simulations)) {
            clonal_populations <- clonal_CN_profiles_all_sims[["variable=clonal_CN_populations"]][[i]]
            #   Find maximal CN profile (most popular among clones)
            loc_max <- which.max(clonal_populations)
            clonal_CN_profiles_all_sims_new[["variable=maximal_CN_profile"]][[i]] <- clonal_CN_profiles_all_sims[["variable=clonal_CN_profiles"]][[i]][[loc_max]]
            #   Find average CN profile (among all clones)
            for (clone in 1:length(clonal_populations)) {
                if (clone == 1) {
                    average_CN_profile <- clonal_CN_profiles_all_sims[["variable=clonal_CN_profiles"]][[i]][[1]]
                    average_CN_profile$Min <- clonal_populations[1] * average_CN_profile$Min
                    average_CN_profile$Maj <- clonal_populations[1] * average_CN_profile$Maj
                } else {
                    average_CN_profile$Min <- average_CN_profile$Min + clonal_populations[clone] * clonal_CN_profiles_all_sims[["variable=clonal_CN_profiles"]][[i]][[clone]]$Min
                    average_CN_profile$Maj <- average_CN_profile$Maj + clonal_populations[clone] * clonal_CN_profiles_all_sims[["variable=clonal_CN_profiles"]][[i]][[clone]]$Maj
                }
            }
            average_CN_profile$Min <- round(average_CN_profile$Min / sum(clonal_populations))
            average_CN_profile$Maj <- round(average_CN_profile$Maj / sum(clonal_populations))
            average_CN_profile$state <- average_CN_profile$Min + average_CN_profile$Maj
            average_CN_profile$copy <- average_CN_profile$state
            clonal_CN_profiles_all_sims_new[["variable=average_CN_profile"]][[i]] <- average_CN_profile
        }
        clonal_CN_profiles_all_sims <- clonal_CN_profiles_all_sims_new
    }
    #-----------------------------------Return the clonal CN information
    return(clonal_CN_profiles_all_sims)
}

#---------------------------------Get statistics for each simulation
get_each_statistics <- function(simulations, simulations_clonal_CN, list_targets) {
    library(vegan)
    library(matrixStats)
    library(transport)
    library(ape)
    library(phyloTop)
    library(treebalance)
    simulations_statistics <- list()
    #----------------------Initialize matrices for simulation statistics
    for (stat in list_targets) {
        stat_details <- strsplit(stat, ";")[[1]]
        stat_target <- strsplit(stat_details[grep("target=", stat_details)], "=")[[1]][2]
        stat_ID <- paste(stat_details[!grepl("statistic=", stat_details)], collapse = ";")
        if (stat_target == "genome") {
            simulations_statistics[[stat_ID]] <- matrix(0, nrow = length(simulations), ncol = 1)
        } else if (stat_target == "chromosome") {
            stat_chromosome_ID <- strsplit(strsplit(stat_details[grep("chromosome=", stat_details)], "=")[[1]][2], ",")[[1]]
            simulations_statistics[[stat_ID]] <- matrix(0, nrow = length(simulations), ncol = length(stat_chromosome_ID))
        }
    }
    #--------------------------------------Compute simulation statistics
    for (i in 1:length(simulations)) {
        simulation <- simulations[[i]]
        #   Find clonal ancestors for single-cell or bulk average-CN-based event counts
        if (any(grepl("variable=event_count", list_targets) & grepl("data=sc", list_targets)) |
            any(grepl("variable=event_count", list_targets) & grepl("data=bulk", list_targets) & grepl("representative_CN=average_CN", list_targets))) {
            evolution_origin <- simulation$clonal_evolution$evolution_origin
            sample_genotype_unique <- simulation$sample$sample_genotype_unique
            evolution_genotype_changes <- simulation$clonal_evolution$evolution_genotype_changes
            subclonal_ancestry <- vector("list", length(sample_genotype_unique))
            for (j in 1:length(sample_genotype_unique)) {
                ancestry <- sample_genotype_unique[j]
                while (ancestry[1] != 0) ancestry <- c(evolution_origin[ancestry[1]], ancestry)
                subclonal_ancestry[[j]] <- ancestry
            }
            clonal_ancestry <- find_clonal_ancestry(subclonal_ancestry)
        }
        #   Find clonal ancestors for bulk maximal-CN-based event counts
        if (any(grepl("variable=event_count", list_targets) & grepl("data=bulk", list_targets) & grepl("representative_CN=maximal_CN", list_targets))) {
            #   ....
            #   ....
            #   ....
        }
        #   Get requested statistics for this simulation
        for (stat in list_targets) {
            stat_details <- strsplit(stat, ";")[[1]]
            stat_target <- strsplit(stat_details[grep("target=", stat_details)], "=")[[1]][2]
            stat_ID <- paste(stat_details[!grepl("statistic=", stat_details)], collapse = ";")
            stat_variable <- strsplit(stat_details[grep("variable=", stat_details)], "=")[[1]][2]
            if (stat_target == "chromosome") {
                stat_chromosome_ID <- strsplit(strsplit(stat_details[grep("chromosome=", stat_details)], "=")[[1]][2], ",")[[1]]
                if (stat_variable == "shannon") {
                    #   Extract CN for each chromosome from each unique clone
                    ls_chrom_profiles <- vector("list", length(stat_chromosome_ID))
                    for (j in 1:length(simulation$sample$sample_genotype_unique)) {
                        genome_profile <- simulation$sample$sample_genotype_unique_profile[[j]]
                        for (k in 1:length(stat_chromosome_ID)) {
                            vec_CN <- paste(genome_profile$copy[genome_profile$chr == stat_chromosome_ID[k]], collapse = "")
                            ls_chrom_profiles[[k]][j] <- vec_CN
                        }
                    }
                    #   Get Shannon index for each chromosome
                    diversity_by_chromosome <- rep(0, length(stat_chromosome_ID))
                    for (j in 1:length(stat_chromosome_ID)) {
                        tmp_frequency_table <- cbind(simulation$sample$sample_genotype_unique, match(ls_chrom_profiles[[j]], unique(ls_chrom_profiles[[j]])))
                        tmp_sample_clone_ID <- simulation$sample$sample_clone_ID
                        for (k in 1:nrow(tmp_frequency_table)) {
                            tmp_sample_clone_ID[which(simulation$sample$sample_clone_ID == tmp_frequency_table[k, 1])] <- tmp_frequency_table[k, 2]
                        }
                        simulations_statistics[[stat_ID]][i, j] <- diversity(table(tmp_sample_clone_ID), index = "shannon")
                    }
                } else if (stat_variable == "event_count") {
                    #   Get source of data
                    stat_data <- strsplit(stat_details[grep("data=", stat_details)], "=")[[1]][2]
                    #   Get count of clonal/subclonal events of a given type
                    clonal_type <- strsplit(stat_details[grep("type=", stat_details)], "=")[[1]][2]
                    event_type <- strsplit(stat_details[grep("event=", stat_details)], "=")[[1]][2]
                    if (stat_data == "sc") {
                        #   Get event count from single-cell samples
                        if (clonal_type == "clonal") {
                            for (j in 1:length(stat_chromosome_ID)) {
                                simulations_statistics[[stat_ID]][i, j] <- find_event_count(clonal_ancestry, event_type, evolution_genotype_changes, by_chromosome = stat_chromosome_ID[j])
                            }
                        } else if (clonal_type == "subclonal") {
                            sample_genotype_event_counts <- rep(0, length(sample_genotype_unique))
                            for (j in 1:length(stat_chromosome_ID)) {
                                for (k in 1:length(sample_genotype_unique)) {
                                    sample_genotype_event_counts[k] <- find_event_count(subclonal_ancestry[[k]], event_type, evolution_genotype_changes, by_chromosome = stat_chromosome_ID[j])
                                }
                                simulations_statistics[[stat_ID]][i, j] <-
                                    sum(sample_genotype_event_counts * table(simulation$sample$sample_clone_ID)) /
                                    length(simulation$sample$sample_clone_ID) -
                                    find_event_count(clonal_ancestry, event_type, evolution_genotype_changes, by_chromosome = stat_chromosome_ID[j])
                            }
                        } else {
                            stop(paste0("Error: Unknown clonal type: ", stat))
                        }
                    } else if (stat_data == "bulk") {
                        #   Get event count from bulk samples
                        representative_CN_type <- strsplit(stat_details[grep("representative_CN=", stat_details)], "=")[[1]][2]
                        if (representative_CN_type == "average_CN") {
                            sample_genotype_event_counts <- rep(0, length(sample_genotype_unique))
                            for (j in 1:length(stat_chromosome_ID)) {
                                for (k in 1:length(sample_genotype_unique)) {
                                    sample_genotype_event_counts[k] <- find_event_count(subclonal_ancestry[[k]], event_type, evolution_genotype_changes, by_chromosome = stat_chromosome_ID[j])
                                }
                                simulations_statistics[[stat_ID]][i, j] <-
                                    sum(sample_genotype_event_counts * table(simulation$sample$sample_clone_ID)) / length(simulation$sample$sample_clone_ID)
                            }
                        } else if (representative_CN_type == "maximal_CN") {
                            simpleError("Error: event count from maximal CN in bulk data not implemented yet")
                        }
                    }
                } else if (stat_variable == "clonal_CN") {
                    #   Get clonal CN profiles and their populations
                    simulations_statistics[["variable=clonal_CN_profiles"]][[i]] <-
                        simulations_clonal_CN[["variable=clonal_CN_profiles"]][[i]]
                    simulations_statistics[["variable=clonal_CN_populations"]][[i]] <-
                        simulations_clonal_CN[["variable=clonal_CN_populations"]][[i]]
                } else if (stat_variable == "maximal_CN") {
                    #   Get maximal CN profile
                    simulations_statistics[["variable=maximal_CN_profile"]][[i]] <-
                        simulations_clonal_CN[["variable=maximal_CN_profile"]][[i]]
                } else if (stat_variable == "average_CN") {
                    #   Get average CN profile
                    simulations_statistics[["variable=average_CN_profile"]][[i]] <-
                        simulations_clonal_CN[["variable=average_CN_profile"]][[i]]
                } else {
                    stop(paste0("Error: Unknown statistic: ", stat))
                }
            } else if (stat_target == "genome") {
                if (stat_variable == "shannon") {
                    #   Get Shannon index
                    simulations_statistics[[stat_ID]][i] <- diversity(table(simulation$sample$sample_clone_ID), index = "shannon")
                } else if (stat_variable == "event_count") {
                    #   Get source of data
                    stat_data <- strsplit(stat_details[grep("data=", stat_details)], "=")[[1]][2]
                    #   Get count of clonal/subclonal events of a given type
                    clonal_type <- strsplit(stat_details[grep("type=", stat_details)], "=")[[1]][2]
                    event_type <- strsplit(stat_details[grep("event=", stat_details)], "=")[[1]][2]
                    if (stat_data == "sc") {
                        #   Get event count from single-cell samples
                        if (clonal_type == "clonal") {
                            simulations_statistics[[stat_ID]][i] <- find_event_count(clonal_ancestry, event_type, evolution_genotype_changes)
                        } else if (clonal_type == "subclonal") {
                            sample_genotype_event_counts <- rep(0, length(sample_genotype_unique))
                            for (j in 1:length(sample_genotype_unique)) {
                                sample_genotype_event_counts[j] <- find_event_count(subclonal_ancestry[[j]], event_type, evolution_genotype_changes)
                            }
                            simulations_statistics[[stat_ID]][i] <-
                                sum(sample_genotype_event_counts * table(simulation$sample$sample_clone_ID)) /
                                length(simulation$sample$sample_clone_ID) -
                                find_event_count(clonal_ancestry, event_type, evolution_genotype_changes)
                        } else {
                            stop(paste0("Error: Unknown clonal type: ", stat))
                        }
                    } else if (stat_data == "bulk") {
                        #   Get event count from bulk samples
                        representative_CN_type <- strsplit(stat_details[grep("representative_CN=", stat_details)], "=")[[1]][2]
                        if (representative_CN_type == "average_CN") {
                            sample_genotype_event_counts <- rep(0, length(sample_genotype_unique))
                            for (j in 1:length(sample_genotype_unique)) {
                                sample_genotype_event_counts[j] <- find_event_count(subclonal_ancestry[[j]], event_type, evolution_genotype_changes)
                            }
                            simulations_statistics[[stat_ID]][i] <-
                                sum(sample_genotype_event_counts * table(simulation$sample$sample_clone_ID)) / length(simulation$sample$sample_clone_ID)
                        } else if (representative_CN_type == "maximal_CN") {
                            simpleError("Error: event count from maximal CN in bulk data not implemented yet")
                        }
                    }
                } else if (stat_variable == "cherries") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get number of cherries
                    simulations_statistics[[stat_ID]][i] <- cherries(tree, normalise = TRUE)
                } else if (stat_variable == "pitchforks") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get number of pitchforks
                    simulations_statistics[[stat_ID]][i] <- pitchforks(tree, normalise = TRUE)
                } else if (stat_variable == "colless") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get colless index
                    simulations_statistics[[stat_ID]][i] <- colless.phylo(tree, normalise = TRUE)
                } else if (stat_variable == "sackin") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get sackin index
                    simulations_statistics[[stat_ID]][i] <- sackin.phylo(tree, normalise = TRUE)
                } else if (stat_variable == "IL_number") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get IL_number
                    simulations_statistics[[stat_ID]][i] <- ILnumber(tree, normalise = TRUE)
                } else if (stat_variable == "avgLadder") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get avgLadder
                    simulations_statistics[[stat_ID]][i] <- avgLadder(tree, normalise = TRUE)
                } else if (stat_variable == "maxDepth") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get maxDepth
                    simulations_statistics[[stat_ID]][i] <- maxDepth(tree)
                } else if (stat_variable == "stairs") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get stairness
                    simulations_statistics[[stat_ID]][i] <- stairs(tree)[1]
                } else if (stat_variable == "B2") {
                    #   Get cell phylogeny tree
                    tree <- simulation$sample_phylogeny$phylogeny_clustering_truth$tree
                    #   Get B2Index
                    simulations_statistics[[stat_ID]][i] <- B2I(tree, logbase = 2)
                } else if (stat_variable == "clonal_CN") {
                    #   Get clonal CN profiles and their populations
                    simulations_statistics[["variable=clonal_CN_profiles"]][[i]] <-
                        simulations_clonal_CN[["variable=clonal_CN_profiles"]][[i]]
                    simulations_statistics[["variable=clonal_CN_populations"]][[i]] <-
                        simulations_clonal_CN[["variable=clonal_CN_populations"]][[i]]
                } else if (stat_variable == "maximal_CN") {
                    #   Get maximal CN profile
                    simulations_statistics[["variable=maximal_CN_profile"]][[i]] <-
                        simulations_clonal_CN[["variable=maximal_CN_profile"]][[i]]
                } else if (stat_variable == "average_CN") {
                    #   Get average CN profile
                    simulations_statistics[["variable=average_CN_profile"]][[i]] <-
                        simulations_clonal_CN[["variable=average_CN_profile"]][[i]]
                } else {
                    stop(paste0("Error: Unknown statistic: ", stat))
                }
            }
        }
    }
    return(simulations_statistics)
}

#----------------------------------------------------Local functions
#---Function to find shared ancestral clones between subclones
find_clonal_ancestry <- function(list_subclonal_ancestry) {
    if (length(list_subclonal_ancestry) == 0) {
        clonal_ancestry <- c()
    } else if (length(list_subclonal_ancestry) == 1) {
        clonal_ancestry <- list_subclonal_ancestry[[1]]
    } else {
        clonal_ancestry <- list_subclonal_ancestry[[1]]
        for (i in 2:length(list_subclonal_ancestry)) {
            clonal_ancestry <- intersect(clonal_ancestry, list_subclonal_ancestry[[i]])
        }
    }
    return(clonal_ancestry)
}
#---Function to find the number of events
#---given a list of clones and event type
find_event_count <- function(clone_ancestry, event_type, evolution_genotype_changes, by_chromosome = NULL) {
    event_count <- 0
    for (clone in clone_ancestry) {
        if (clone <= 0) next
        if (length(evolution_genotype_changes[[clone]]) == 0) next
        for (j in 1:length(evolution_genotype_changes[[clone]])) {
            if (evolution_genotype_changes[[clone]][[j]][1] == event_type) {
                if (is.null(by_chromosome)) {
                    event_count <- event_count + 1
                } else {
                    if (evolution_genotype_changes[[clone]][[j]][2] == by_chromosome) {
                        event_count <- event_count + 1
                    }
                }
            }
        }
    }
    return(event_count)
}
#---Function to find distance between two CN profiles
cn_distance <- function(from_CN_profile,
                        to_CN_profile,
                        metric) {
    if (metric == "euclidean") {
        distance <- sqrt(sum((from_CN_profile - to_CN_profile)^2))
    } else if (metric == "wcnd") {
        distance <- 0
    } else if (metric == "correlation") {
        library(energy)
        distance <- dcor(from_CN_profile, to_CN_profile)
    } else if (metric == "mahalanobis") {
        library(MASS)
        distance <- sqrt(t(from_CN_profile - to_CN_profile) %*% ginv(cov(from_CN_profile %*% t(from_CN_profile))) %*% (from_CN_profile - to_CN_profile))
    }
    return(distance)
}

#---Function to find distance between two samples
#---based on clonal population and CN profiles
sample_distance <- function(from_clonal_CN_populations,
                            from_clonal_CN_profiles,
                            to_clonal_CN_populations,
                            to_clonal_CN_profiles,
                            metric) {
    library(transport)
    #   Compute distance between any pair of CN profiles in the samples
    n_sample_from <- length(from_clonal_CN_populations)
    n_sample_to <- length(to_clonal_CN_populations)
    cost_matrix <- matrix(nrow = n_sample_from, ncol = n_sample_to)
    for (i in 1:n_sample_from) {
        for (j in 1:n_sample_to) {
            from_CN_profile <- from_clonal_CN_profiles[[i]]$copy
            to_CN_profile <- to_clonal_CN_profiles[[j]]$copy
            cost_matrix[i, j] <- cn_distance(from_CN_profile, to_CN_profile, metric)
        }
    }
    #   Normalize clonal populations
    from_clonal_CN_populations <- from_clonal_CN_populations / sum(from_clonal_CN_populations)
    to_clonal_CN_populations <- to_clonal_CN_populations / sum(to_clonal_CN_populations)
    #   Compute the Wasserstein distance between the two samples
    wasserstein_dist <- wasserstein(
        from_clonal_CN_populations,
        to_clonal_CN_populations,
        costm = cost_matrix, p = 1
    )
    return(wasserstein_dist)
}
#---Function to find distance between two sample cohorts
#---based on clonal population and CN profiles
cohort_distance <- function(cohort_from, cohort_to, metric, bulk_CN_input = "", bulk = FALSE, by_chromosome = NULL) {
    library(transport)
    #-----------------------------Define cost matrix between two cohorts
    if (is.null(by_chromosome)) {
        if (bulk) {
            if (bulk_CN_input == "average") {
                n_cohort_from <- length(cohort_from[["variable=average_CN_profile"]])
                n_cohort_to <- length(cohort_to[["variable=average_CN_profile"]])
            } else if (bulk_CN_input == "maximal") {
                n_cohort_from <- length(cohort_from[["variable=maximal_CN_profile"]])
                n_cohort_to <- length(cohort_to[["variable=maximal_CN_profile"]])
            }
            sample_distance_mtx <- matrix(0, nrow = n_cohort_from, ncol = n_cohort_to)
            for (i in 1:n_cohort_from) {
                for (j in 1:n_cohort_to) {
                    if (bulk_CN_input == "average") {
                        from_CN_profile <- cohort_from[["variable=average_CN_profile"]][[i]]$copy
                        to_CN_profile <- cohort_to[["variable=average_CN_profile"]][[j]]$copy
                    } else if (bulk_CN_input == "maximal") {
                        from_CN_profile <- cohort_from[["variable=maximal_CN_profile"]][[i]]$copy
                        to_CN_profile <- cohort_to[["variable=maximal_CN_profile"]][[j]]$copy
                    }
                    sample_distance_mtx[i, j] <- cn_distance(from_CN_profile, to_CN_profile, metric)
                }
            }
        } else {
            #   Compute distance between any pair of samples in the cohorts
            n_cohort_from <- length(cohort_from[["variable=clonal_CN_populations"]])
            n_cohort_to <- length(cohort_to[["variable=clonal_CN_populations"]])
            sample_distance_mtx <- matrix(0, nrow = n_cohort_from, ncol = n_cohort_to)
            for (i in 1:n_cohort_from) {
                for (j in 1:n_cohort_to) {
                    from_clonal_CN_populations <- cohort_from[["variable=clonal_CN_populations"]][[i]]
                    from_clonal_CN_profiles <- cohort_from[["variable=clonal_CN_profiles"]][[i]]
                    to_clonal_CN_populations <- cohort_to[["variable=clonal_CN_populations"]][[j]]
                    to_clonal_CN_profiles <- cohort_to[["variable=clonal_CN_profiles"]][[j]]
                    sample_distance_mtx[i, j] <- sample_distance(
                        from_clonal_CN_populations = from_clonal_CN_populations,
                        from_clonal_CN_profiles = from_clonal_CN_profiles,
                        to_clonal_CN_populations = to_clonal_CN_populations,
                        to_clonal_CN_profiles = to_clonal_CN_profiles,
                        metric = metric
                    )
                }
            }
        }
        #---------------------------------------Normalize clonal populations
        dist_from <- rep(1 / n_cohort_from, n_cohort_from)
        dist_to <- rep(1 / n_cohort_to, n_cohort_to)
        #-----------Compute the Wasserstein distance between the two cohorts
        wasserstein_dist <- wasserstein(
            dist_from,
            dist_to,
            costm = sample_distance_mtx, p = 1
        )
    } else {
        if (bulk) {
            if (bulk_CN_input == "average") {
                n_cohort_from <- length(cohort_from[["variable=average_CN_profile"]])
                n_cohort_to <- length(cohort_to[["variable=average_CN_profile"]])
            } else if (bulk_CN_input == "maximal") {
                n_cohort_from <- length(cohort_from[["variable=maximal_CN_profile"]])
                n_cohort_to <- length(cohort_to[["variable=maximal_CN_profile"]])
            }
            sample_distance_mtx <- matrix(0, nrow = n_cohort_from, ncol = n_cohort_to)
            for (i in 1:n_cohort_from) {
                for (j in 1:n_cohort_to) {
                    if (bulk_CN_input == "average") {
                        from_CN_profile_genome <- cohort_from[["variable=average_CN_profile"]][[i]]
                        to_CN_profile_genome <- cohort_to[["variable=average_CN_profile"]][[j]]
                        from_CN_profile <- from_CN_profile_genome[which(from_CN_profile_genome$chr == by_chromosome), ]$copy
                        to_CN_profile <- to_CN_profile_genome[which(to_CN_profile_genome$chr == by_chromosome), ]$copy
                    } else if (bulk_CN_input == "maximal") {
                        from_CN_profile_genome <- cohort_from[["variable=maximal_CN_profile"]][[i]]
                        to_CN_profile_genome <- cohort_to[["variable=maximal_CN_profile"]][[j]]
                        from_CN_profile <- from_CN_profile_genome[which(from_CN_profile_genome$chr == by_chromosome), ]$copy
                        to_CN_profile <- to_CN_profile_genome[which(to_CN_profile_genome$chr == by_chromosome), ]$copy
                    }
                    sample_distance_mtx[i, j] <- cn_distance(from_CN_profile, to_CN_profile, metric)
                }
            }
        } else {
            #   Compute distance between any pair of samples in the cohorts
            n_cohort_from <- length(cohort_from[["variable=clonal_CN_populations"]])
            n_cohort_to <- length(cohort_to[["variable=clonal_CN_populations"]])
            sample_distance_mtx <- matrix(0, nrow = n_cohort_from, ncol = n_cohort_to)
            for (i in 1:n_cohort_from) {
                for (j in 1:n_cohort_to) {
                    from_clonal_CN_populations <- cohort_from[["variable=clonal_CN_populations"]][[i]]
                    from_clonal_CN_profiles_genome <- cohort_from[["variable=clonal_CN_profiles"]][[i]]
                    from_clonal_CN_profiles <- list()
                    for (k in 1:length(from_clonal_CN_profiles_genome)) {
                        from_clonal_CN_profiles[[k]] <- from_clonal_CN_profiles_genome[[k]][which(from_clonal_CN_profiles_genome[[k]]$chr == by_chromosome), ]
                    }
                    to_clonal_CN_populations <- cohort_to[["variable=clonal_CN_populations"]][[j]]
                    to_clonal_CN_profiles_genome <- cohort_to[["variable=clonal_CN_profiles"]][[j]]
                    to_clonal_CN_profiles <- list()
                    for (k in 1:length(to_clonal_CN_profiles_genome)) {
                        to_clonal_CN_profiles[[k]] <- to_clonal_CN_profiles_genome[[k]][which(to_clonal_CN_profiles_genome[[k]]$chr == by_chromosome), ]
                    }
                    sample_distance_mtx[i, j] <- sample_distance(
                        from_clonal_CN_populations = from_clonal_CN_populations,
                        from_clonal_CN_profiles = from_clonal_CN_profiles,
                        to_clonal_CN_populations = to_clonal_CN_populations,
                        to_clonal_CN_profiles = to_clonal_CN_profiles,
                        metric = metric
                    )
                }
            }
        }
        #---------------------------------------Normalize clonal populations
        dist_from <- rep(1 / n_cohort_from, n_cohort_from)
        dist_to <- rep(1 / n_cohort_to, n_cohort_to)
        #-----------Compute the Wasserstein distance between the two cohorts
        wasserstein_dist <- wasserstein(
            dist_from,
            dist_to,
            costm = sample_distance_mtx, p = 1
        )
    }
    return(wasserstein_dist)
}

#' @export
get_statistics <- function(list_targets,
                           n_simulations_sc = NULL,
                           n_simulations_bulk = NULL,
                           simulations_sc = NULL,
                           simulations_bulk = NULL,
                           simulations_statistics_sc = NULL,
                           simulations_statistics_bulk = NULL,
                           cn_data_sc = NULL,
                           cn_data_bulk = NULL,
                           arm_level = FALSE,
                           cn_table = NULL,
                           save_sample_statistics = FALSE) {
    #-------------Get clonal CN profiles for all single-cell simulations
    if (is.null(simulations_statistics_sc)) {
        if (!is.null(n_simulations_sc)) {
            simulations_sc <- simulations_sc[1:n_simulations_sc]
        }
        if (any(grepl("variable=clonal_CN", list_targets))) {
            simulations_clonal_CN_sc <- get_each_clonal_CN_profiles(simulations_sc, arm_level, cn_table)
        }
        list_targets_sc <- list_targets[grepl("data=sc", list_targets)]
        simulations_statistics_sc <- get_each_statistics(simulations_sc, simulations_clonal_CN_sc, list_targets_sc)
    } else {
        if (!is.null(n_simulations_sc)) {
            statistics <- names(simulations_statistics_sc)
            for (statistic in statistics) {
                simulations_statistics_sc[[statistic]] <- simulations_statistics_sc[[statistic]][1:n_simulations_sc, ]
            }
        }
    }
    #------------Get representative CN profiles for all bulk simulations
    if (is.null(simulations_statistics_bulk)) {
        if (!is.null(n_simulations_bulk)) {
            simulations_bulk <- simulations_bulk[1:n_simulations_bulk]
        }
        if ((any(grepl("variable=average_CN", list_targets))) | (any(grepl("variable=maximal_CN", list_targets)))) {
            simulations_clonal_CN_bulk <- get_each_clonal_CN_profiles(simulations_bulk, arm_level, cn_table, bulk = TRUE)
        }
        list_targets_bulk <- list_targets[grepl("data=bulk", list_targets)]
        simulations_statistics_bulk <- get_each_statistics(simulations_bulk, simulations_clonal_CN_bulk, list_targets_bulk)
    } else {
        if (!is.null(n_simulations_bulk)) {
            statistics <- names(simulations_statistics_bulk)
            for (statistic in statistics) {
                simulations_statistics_bulk[[statistic]] <- simulations_statistics_bulk[[statistic]][1:n_simulations_bulk, ]
            }
        }
    }
    #-----------------------------------------Get statistics for fitting
    # statistics <- vector("list", length(list_targets))
    statistics <- list()
    for (stat in list_targets) {
        stat_details <- strsplit(stat, ";")[[1]]
        stat_ID <- paste(stat_details[!grepl("statistic=", stat_details)], collapse = ";")
        stat_data <- strsplit(stat_details[grep("data=", stat_details)], "=")[[1]][2]
        stat_variable <- strsplit(stat_details[grep("variable=", stat_details)], "=")[[1]][2]
        stat_type <- strsplit(stat_details[grepl("statistic=", stat_details)], "=")[[1]][2]
        stat_target <- strsplit(stat_details[grep("target=", stat_details)], "=")[[1]][2]
        if (stat_type == "mean") {
            if (stat_data == "sc") {
                statistics[[stat]] <- apply(simulations_statistics_sc[[stat_ID]], 2, mean)
            } else if (stat_data == "bulk") {
                statistics[[stat]] <- apply(simulations_statistics_bulk[[stat_ID]], 2, mean)
            }
        } else if (stat_type == "var") {
            if (stat_data == "sc") {
                statistics[[stat]] <- apply(simulations_statistics_sc[[stat_ID]], 2, var)
            } else if (stat_data == "bulk") {
                statistics[[stat]] <- apply(simulations_statistics_bulk[[stat_ID]], 2, var)
            }
        } else if (stat_type == "dist") {
            if (stat_target == "chromosome") {
                stat_chromosome_ID <- strsplit(strsplit(stat_details[grep("chromosome=", stat_details)], "=")[[1]][2], ",")[[1]]
                if (stat_variable == "clonal_CN" & stat_data == "sc") {
                    stat_metric <- strsplit(stat_details[grepl("metric=", stat_details)], "=")[[1]][2]
                    if (is.null(cn_data_sc)) cn_data_sc <- simulations_statistics_sc
                    for (j in 1:length(stat_chromosome_ID)) {
                        statistics[[stat]][j] <- cohort_distance(
                            cohort_from = simulations_statistics_sc,
                            cohort_to = cn_data_sc,
                            metric = stat_metric,
                            by_chromosome = stat_chromosome_ID[j]
                        )
                    }
                } else if (stat_variable == "average_CN" & stat_data == "bulk") {
                    stat_metric <- strsplit(stat_details[grepl("metric=", stat_details)], "=")[[1]][2]
                    if (is.null(cn_data_bulk)) cn_data_bulk <- simulations_statistics_bulk
                    for (j in 1:length(stat_chromosome_ID)) {
                        statistics[[stat]][j] <- cohort_distance(
                            cohort_from = simulations_statistics_bulk,
                            cohort_to = cn_data_bulk,
                            metric = stat_metric,
                            bulk_CN_input = "average",
                            bulk = TRUE,
                            by_chromosome = stat_chromosome_ID[j]
                        )
                    }
                } else if (stat_variable == "maximal_CN" & stat_data == "bulk") {
                    stat_metric <- strsplit(stat_details[grepl("metric=", stat_details)], "=")[[1]][2]
                    if (is.null(cn_data_bulk)) cn_data_bulk <- simulations_statistics_bulk
                    for (j in 1:length(stat_chromosome_ID)) {
                        statistics[[stat]][j] <- cohort_distance(
                            cohort_from = simulations_statistics_bulk,
                            cohort_to = cn_data_bulk,
                            metric = stat_metric,
                            bulk_CN_input = "maximal",
                            bulk = TRUE,
                            by_chromosome = stat_chromosome_ID[j]
                        )
                    }
                }
            } else {
                if (stat_variable == "clonal_CN" & stat_data == "sc") {
                    stat_metric <- strsplit(stat_details[grepl("metric=", stat_details)], "=")[[1]][2]
                    if (is.null(cn_data_sc)) cn_data_sc <- simulations_statistics_sc
                    statistics[[stat]] <- cohort_distance(
                        cohort_from = simulations_statistics_sc,
                        cohort_to = cn_data_sc,
                        metric = stat_metric
                    )
                } else if (stat_variable == "average_CN" & stat_data == "bulk") {
                    stat_metric <- strsplit(stat_details[grepl("metric=", stat_details)], "=")[[1]][2]
                    if (is.null(cn_data_bulk)) cn_data_bulk <- simulations_statistics_bulk
                    statistics[[stat]] <- cohort_distance(
                        cohort_from = simulations_statistics_bulk,
                        cohort_to = cn_data_bulk,
                        metric = stat_metric,
                        bulk_CN_input = "average",
                        bulk = TRUE
                    )
                } else if (stat_variable == "maximal_CN" & stat_data == "bulk") {
                    stat_metric <- strsplit(stat_details[grepl("metric=", stat_details)], "=")[[1]][2]
                    if (is.null(cn_data_bulk)) cn_data_bulk <- simulations_statistics_bulk
                    statistics[[stat]] <- cohort_distance(
                        cohort_from = simulations_statistics_bulk,
                        cohort_to = cn_data_bulk,
                        metric = stat_metric,
                        bulk_CN_input = "maximal",
                        bulk = TRUE
                    )
                }
            }
        } else {
            stop(paste0("Error: Unknown statistic type: ", stat))
        }
    }
    #-----------------------------------------------------Prepare output
    output <- list()
    output$statistics <- statistics
    if (save_sample_statistics) {
        output$simulations_statistics_sc <- simulations_statistics_sc
        output$simulations_statistics_bulk <- simulations_statistics_bulk
    }
    return(output)
}

#---Function to assign parameters to proper positions
assign_paras <- function(model_variables, parameter_IDs, parameters) {
    for (i in 1:length(parameter_IDs)) {
        parameter_ID_input <- parameter_IDs[i]
        parameter_value_input <- parameters[i]
        #   Prepare values for operation on parameter
        if (grepl(":", parameter_ID_input)) {
            parameter_ID <- sub(".*:", "", parameter_ID_input)
            parameter_operator <- sub(":.*", "", parameter_ID_input)
            parameter_operator <- paste0(parameter_operator, "(parameter_value_input)")
        } else {
            parameter_ID <- parameter_ID_input
            parameter_operator <- "parameter_value_input"
        }
        parameter_value <- eval(parse(text = parameter_operator))
        #   Input parameter
        if (parameter_ID %in% model_variables$general_variables$Variable) {
            model_variables$general_variables$Value[which(model_variables$general_variables$Variable == parameter_ID)] <- parameter_value
        } else if (parameter_ID %in% model_variables$chromosome_arm_library$Arm_ID) {
            model_variables$chromosome_arm_library$s_rate[which(model_variables$chromosome_arm_library$Arm_ID == parameter_ID)] <- parameter_value
        }
    }
    return(model_variables)
}
#---Objective function for ABC fitting
func_stat <- function(parameters,
                      parameter_IDs,
                      model_variables,
                      list_targets,
                      arm_level,
                      n_simulations_sc,
                      n_simulations_bulk) {
    #   Assign parameters in model variables
    model_variables <- assign_paras(model_variables, parameter_IDs, parameters)
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
    list_targets_sc <- list_targets[grepl("data=sc", list_targets)]
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
    list_targets_bulk <- list_targets[grepl("data=bulk", list_targets)]
    simulations_statistics_bulk <- get_each_statistics(simulations_bulk, simulations_clonal_CN_bulk, list_targets_bulk)
    #   Make output list
    output <- list()
    output$parameter_IDs <- parameter_IDs
    output$parameters <- parameters
    output$stats_sc <- simulations_statistics_sc
    output$stats_bulk <- simulations_statistics_bulk
    return(output)
}

func_ABC <- function(parameters,
                     parameter_IDs,
                     model_variables,
                     list_targets,
                     cn_data_sc = NULL,
                     cn_data_bulk = NULL,
                     arm_level = FALSE,
                     save_sample_statistics = FALSE,
                     n_simulations_sc,
                     n_simulations_bulk) {
    #   Assign parameters in model variables
    model_variables <- assign_paras(model_variables, parameter_IDs, parameters)
    #   Make single-cell simulations
    SIMS_sc <- simulator_full_program(
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
    #   Make bulk simulations
    SIMS_bulk <- simulator_full_program(
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
    #   Get statistics from simulations
    stat <- get_statistics(
        simulations_sc = SIMS_sc,
        simulations_bulk = SIMS_bulk,
        cn_data_sc = cn_data_sc,
        cn_data_bulk = cn_data_bulk,
        list_targets = list_targets,
        arm_level = arm_level,
        cn_table = cn_table,
        save_sample_statistics = save_sample_statistics
    )
    return(stat)
}

#' @export
library_simulations <- function(library_name,
                                model_variables,
                                list_parameters,
                                list_targets_library,
                                ABC_simcount_start = 0,
                                ABC_simcount = 10000,
                                arm_level = FALSE,
                                cn_table = NULL,
                                n_simulations_sc,
                                n_simulations_bulk,
                                ##############################################
                                n_cores = NULL) {
    library(parallel)
    library(pbapply)
    library(abcrf)
    library(grid)
    library(gridExtra)
    library(ggplot2)
    library(signals)
    if (is.null(n_cores)) n_cores <- max(detectCores() - 1, 1)
    #---------------------------------List of parameter IDs to be fitted
    parameter_IDs <- list_parameters$Variable
    # =============================================CREATE REFERENCE DATA
    #---------------------------------------Simulate table of parameters
    sim_param <- matrix(0, nrow = ABC_simcount, ncol = nrow(list_parameters))
    for (col in 1:ncol(sim_param)) {
        sim_param[, col] <- runif(
            ABC_simcount,
            min = as.numeric(list_parameters$Lower_bound[col]),
            max = as.numeric(list_parameters$Upper_bound[col])
        )
    }
    #-----------------------------------------------Make reference table
    start_time <- Sys.time()
    #   Configure parallel pool
    cl <- makePSOCKcluster(n_cores)
    cat("Creating reference table for ABC...\n")
    model_name <<- model_name
    sim_param <<- sim_param
    parameter_IDs <<- parameter_IDs
    model_variables <<- model_variables
    func_stat <<- func_stat
    assign_paras <<- assign_paras
    list_targets_library <<- list_targets_library
    cn_table <<- cn_table
    arm_level <<- arm_level
    n_simulations_sc <<- n_simulations_sc
    n_simulations_bulk <<- n_simulations_bulk
    clusterExport(cl, varlist = c(
        "list_targets_library", "sim_param", "parameter_IDs", "model_name", "model_variables", "cn_table", "arm_level",
        "func_stat", "n_simulations_sc", "n_simulations_bulk", "assign_paras", "get_statistics", "get_each_clonal_CN_profiles", "get_cn_profile", "get_arm_CN_profiles",
        "find_clonal_ancestry", "find_event_count", "cn_distance", "sample_distance", "cohort_distance", "get_each_statistics",
        "vec_CN_block_no", "vec_centromeres",
        "BUILD_driver_library", "simulator_full_program", "one_simulation",
        "SIMULATOR_VARIABLES_for_simulation",
        "SIMULATOR_FULL_PHASE_1_main", "SIMULATOR_FULL_PHASE_1_clonal_population_cleaning",
        "SIMULATOR_FULL_PHASE_1_CN_chrom_arm_missegregation", "SIMULATOR_FULL_PHASE_1_CN_cnloh_interstitial", "SIMULATOR_FULL_PHASE_1_CN_cnloh_terminal", "SIMULATOR_FULL_PHASE_1_CN_focal_amplification", "SIMULATOR_FULL_PHASE_1_CN_focal_deletion", "SIMULATOR_FULL_PHASE_1_CN_missegregation", "SIMULATOR_FULL_PHASE_1_CN_whole_genome_duplication", "SIMULATOR_FULL_PHASE_1_drivers",
        "SIMULATOR_FULL_PHASE_1_genotype_cleaning", "SIMULATOR_FULL_PHASE_1_genotype_comparison", "SIMULATOR_FULL_PHASE_1_genotype_initiation", "SIMULATOR_FULL_PHASE_1_genotype_update", "SIMULATOR_FULL_PHASE_1_selection_rate",
        "SIMULATOR_FULL_PHASE_2_main",
        "SIMULATOR_FULL_PHASE_3_main",
        "p2_cn_profiles_long", "p2_cn_profiles_wide", "p2_readcount_model"
    ))
    e <- new.env()
    e$libs <- .libPaths()
    clusterExport(cl, "libs", envir = e)
    clusterEvalQ(cl, .libPaths(libs))
    #   Create simulated results in parallel
    pbo <- pboptions(type = "txt")
    dir.create(library_name)
    sim_results_list <- pblapply(cl = cl, X = (ABC_simcount_start + 1):(ABC_simcount_start + ABC_simcount), FUN = function(iteration) {
        simulation_parameters <- sim_param[iteration - ABC_simcount_start, ]
        simulation_statistics <- func_stat(
            parameters = simulation_parameters, parameter_IDs = parameter_IDs, model_variables = model_variables, list_targets = list_targets_library, arm_level = arm_level, n_simulations_sc = n_simulations_sc, n_simulations_bulk = n_simulations_bulk
        )
        filename <- paste0(library_name, "/", library_name, "_ABC_simulation_statistics_", iteration, ".rda")
        save(simulation_statistics, file = filename)
        return(simulation_statistics)
    })
    stopCluster(cl)
    end_time <- Sys.time()
    print(end_time - start_time)
}

#' @export
library_statistics <- function(library_name,
                               library_sensitivity_name = NULL,
                               model_variables,
                               list_parameters,
                               list_targets_library,
                               ABC_simcount_start = NULL,
                               ABC_simcount,
                               compute_parallel = TRUE,
                               n_cores = NULL,
                               cn_data_sc = NULL,
                               cn_data_bulk = NULL,
                               cn_table = NULL,
                               arm_level = FALSE) {
    library(parallel)
    library(pbapply)
    library(abcrf)
    library(grid)
    library(gridExtra)
    library(ggplot2)
    library(signals)
    library(data.table)
    library(matrixStats)
    library(combinat)
    # ================FIND STATISTICS FOR EACH SIMULATION IN THE LIBRARY
    if (is.null(ABC_simcount_start)) {
        ABC_simcount_start_real <- 1
        ABC_simcount_end_real <- ABC_simcount
    } else {
        ABC_simcount_start_real <- ABC_simcount_start + 1
        ABC_simcount_end_real <- ABC_simcount_start + ABC_simcount
    }
    if (is.null(library_sensitivity_name)) {
        library_sensitivity_name <- library_name
    }
    #   Get statistics for each simulation w.r.t. data
    if (compute_parallel == FALSE) {
        sim_output_list <- vector("list", ABC_simcount)
        for (iteration in ABC_simcount_start_real:ABC_simcount_end_real) {
            print(paste0("iteration: ", iteration))
            start_time <- Sys.time()
            filename <- paste0(library_name, "/", library_name, "_ABC_simulation_statistics_", iteration, ".rda")
            load(filename)
            stat <- get_statistics(
                simulations_statistics_sc = simulation_statistics$stats_sc,
                simulations_statistics_bulk = simulation_statistics$stats_bulk,
                cn_data_sc = cn_data_sc,
                cn_data_bulk = cn_data_bulk,
                list_targets = list_targets_library,
                arm_level = arm_level,
                cn_table = cn_table
            )
            output <- list()
            output$parameters <- simulation_statistics$parameters
            output$stat <- stat
            sim_output_list[[iteration - ABC_simcount_start]] <- output
            end_time <- Sys.time()
            print(end_time - start_time)
        }
    } else {
        if (is.null(n_cores)) {
            n_cores <- max(detectCores() - 1, 1)
        }
        cl <- makePSOCKcluster(n_cores)
        cat("\nGetting statistics...\n")
        get_statistics <<- get_statistics
        ABC_simcount <<- ABC_simcount
        list_targets_library <<- list_targets_library
        library_name <<- library_name
        cohort_distance <<- cohort_distance
        cn_data_sc <<- cn_data_sc
        cn_data_bulk <<- cn_data_bulk
        arm_level <<- arm_level
        cn_table <<- cn_table
        clusterExport(cl, varlist = c(
            "ABC_simcount", "get_statistics", "list_targets_library", "library_name", "cn_data_sc",
            "cn_data_bulk", "arm_level", "cn_table", "cohort_distance", "cn_distance", "sample_distance"
        ))
        e <- new.env()
        e$libs <- .libPaths()
        clusterExport(cl, "libs", envir = e)
        clusterEvalQ(cl, .libPaths(libs))
        #   Create simulated results in parallel
        pbo <- pboptions(type = "txt")
        sim_output_list <- pblapply(cl = cl, X = ABC_simcount_start_real:ABC_simcount_end_real, FUN = function(iteration) {
            filename <- paste0(library_name, "/", library_name, "_ABC_simulation_statistics_", iteration, ".rda")
            load(filename)
            stat <- get_statistics(
                simulations_statistics_sc = simulation_statistics$stats_sc,
                simulations_statistics_bulk = simulation_statistics$stats_bulk,
                cn_data_sc = cn_data_sc,
                cn_data_bulk = cn_data_bulk,
                list_targets = list_targets_library,
                arm_level = arm_level,
                cn_table = cn_table
            )
            output <- list()
            output$parameters <- simulation_statistics$parameters
            output$stat <- stat
            return(output)
        })
        stopCluster(cl)
    }
    #   Group simulated statistics into one table
    sim_param <- NULL
    sim_stat <- list()
    for (row in 1:ABC_simcount) {
        sim_param <- rbind(sim_param, sim_output_list[[row]]$parameters)
        for (stat in list_targets_library) {
            sim_stat[[stat]] <- rbind(sim_stat[[stat]], sim_output_list[[row]]$stat$statistics[[stat]])
        }
    }
    sim_param <- data.frame(sim_param)
    colnames(sim_param) <- list_parameters$Variable
    # ======================SAVE RESULTS FROM GET STATISTICS TO LIBRARY
    dir.create(library_sensitivity_name)
    ABC_input <- list()
    ABC_input$model_variables <- model_variables
    ABC_input$sim_param <- sim_param
    ABC_input$sim_stat <- sim_stat
    if (is.null(ABC_simcount_start)) {
        filename <- paste0(library_sensitivity_name, "/", library_sensitivity_name, "_ABC_input.rda")
    } else {
        filename <- paste0(library_sensitivity_name, "/", library_sensitivity_name, "_ABC_input_", (ABC_simcount_start + ABC_simcount) / ABC_simcount, ".rda")
    }
    save(ABC_input, file = filename)
}

#' @export
fitting_parameters <- function(library_name,
                               copynumber_DATA,
                               parameters_truth = NULL,
                               rda_name = NULL,
                               list_parameters,
                               list_targets_by_parameter,
                               n_cores = NULL,
                               plot_ABC_prior_as_uniform = FALSE) {
    library(parallel)
    library(pbapply)
    library(abcrf)
    library(grid)
    library(gridExtra)
    library(ggplot2)
    library(signals)
    library(data.table)
    library(matrixStats)
    library(combinat)
    if (is.null(n_cores)) {
        n_cores <- max(detectCores() - 1, 1)
    }
    # ======================================================LOAD LIBRARY
    if (is.null(rda_name)) {
        rda_name <- library_name
    }
    filename <- paste0(rda_name, "_ABC_input.rda")
    load(filename)

    model_variables <- ABC_input$model_variables
    sim_param <- ABC_input$sim_param
    sim_stat <- ABC_input$sim_stat
    ABC_input <- list()
    # ====================================FITTING WITH ABC RANDOM FOREST
    #---List of target statistics
    list_target_statistics <- colnames(list_targets_by_parameter)[-1]
    #---Fit each parameter with ABC-rf
    list_parameters_output <- list_parameters
    list_parameters_output$Mean <- 0
    list_parameters_output$Median <- 0
    list_parameters_output$Mode <- 0
    if (!is.null(parameters_truth)) {
        list_parameters_output$Ground_truth <- 0
        for (i in 1:length(parameters_truth$Variable)) {
            list_parameters_output$Ground_truth[which(list_parameters_output$Variable == parameters_truth$Variable[i])] <- parameters_truth$Value[i]
        }
    }
    nrow_plot <- 6
    layout <- matrix(NA, nrow = nrow_plot, ncol = ceiling(nrow(list_parameters) / nrow_plot))
    gs <- list()
    id <- 0
    for (para in 1:nrow(list_parameters)) {
        start_time <- Sys.time()
        para_ID <- list_parameters$Variable[para]
        para_type <- list_parameters$Type[para]
        cat(paste("\nABC for parameter ", para_ID, " [", para, "/", nrow(list_parameters), "]", "\n", sep = ""))
        #   Prepare matrices of prepared statistics library & data observation
        flags_targets <- list_targets_by_parameter[which(list_targets_by_parameter$Variable == para_ID), 2:ncol(list_targets_by_parameter)]
        list_targets <- list_target_statistics[which(flags_targets == 1)]
        mini_data <- NULL
        mini_obs <- NULL
        for (i in 1:length(list_targets)) {
            stat <- list_targets[i]
            stat_details <- strsplit(stat, ";")[[1]]
            stat_target <- strsplit(stat_details[grep("target=", stat_details)], "=")[[1]][2]
            if (stat_target == "genome") {
                mini_data <- cbind(mini_data, sim_stat[[stat]])
                mini_obs <- cbind(mini_obs, copynumber_DATA$statistics[[stat]])
            } else if (stat_target == "chromosome") {
                selected_chrom <- list_parameters$Chromosome[para]
                stat_chromosome_ID <- strsplit(strsplit(stat_details[grep("chromosome=", stat_details)], "=")[[1]][2], ",")[[1]]
                selected_column <- which(stat_chromosome_ID == selected_chrom)
                mini_data <- cbind(mini_data, sim_stat[[stat]][, selected_column])
                mini_obs <- cbind(mini_obs, copynumber_DATA$statistics[[stat]][selected_column])
            }
        }
        mini_data <- data.frame(mini_data)
        colnames(mini_data) <- paste0("stat_", 1:ncol(mini_data))
        mini_obs <- data.frame(matrix(mini_obs, nrow = 1))
        colnames(mini_obs) <- paste0("stat_", 1:ncol(mini_obs))
        #   Prepare library of parameters for this parameter
        data_rf <- cbind(sim_param[para_ID], mini_data)
        #   Train the random forest
        colnames(data_rf)[1] <- "para"
        f <- as.formula("para ~.")
        model_rf <- regAbcrf(
            formula = f, data_rf,
            paral = TRUE, ncores = n_cores,
            # ntree = ntree,
            # sampsize = nrow(data_rf),
            # save.memory = TRUE
        )
        #   Predict posterior distribution based on found random forest
        post_rf <- predict(model_rf, mini_obs, data_rf, paral = TRUE, ncores = n_cores)
        #   Choose best values from posterior distribution
        df_dist <- densityPlot_df(
            object = model_rf, obs = mini_obs, training = data_rf
        )
        dist_output <- list()
        dist_output$x <- df_dist$x
        dist_output$y_prior <- df_dist$y_prior
        dist_output$y_posterior <- df_dist$y_posterior
        filename <- paste0(library_name, "_dist_output_", para_ID, ".rda")
        save(dist_output, file = filename)

        post_mean <- weighted.mean(df_dist$x, df_dist$y_posterior)
        post_sd <- weightedSd(df_dist$x, df_dist$y_posterior)
        post_median <- weightedMedian(df_dist$x, df_dist$y_posterior)
        post_mode <- df_dist$x[which(df_dist$y_posterior == max(df_dist$y_posterior))]
        if (!is.null(parameters_truth)) cat("True value: ", parameters_truth$Value[which(parameters_truth$Variable == para_ID)], "\n")
        cat("Posterior mean: ", post_mean, "\n")
        cat("Posterior median: ", post_median, "\n")
        cat("Posterior mode: ", post_mode, "\n")
        cat("Posterior sd: ", post_sd, "\n")
        #   Save the parameter output table
        list_parameters_output$Mean[which(list_parameters_output$Variable == para_ID)] <- post_mean
        list_parameters_output$Median[which(list_parameters_output$Variable == para_ID)] <- post_median
        list_parameters_output$Mode[which(list_parameters_output$Variable == para_ID)] <- post_mode
        list_parameters_output$Sd[which(list_parameters_output$Variable == para_ID)] <- post_sd
        #   Save results for fitting this parameter
        ABC_output <- list()
        ABC_output$para_ID <- para_ID
        ABC_output$post_mode <- post_mode
        ABC_output$post_mean <- post_mean
        ABC_output$post_median <- post_median
        ABC_output$post_sd <- post_sd
        filename <- paste0(library_name, "_ABC_output_", para_ID, ".rda")
        save(ABC_output, file = filename)
        #   Plot the prior, posterior and chosen best parameter for all variables
        true_para <- NULL
        if (!is.null(parameters_truth)) {
            true_para <- parameters_truth$Value[which(parameters_truth$Variable == para_ID)]
        }
        id <- id + 1
        row <- id %% nrow_plot
        if (row == 0) row <- nrow_plot
        col <- ceiling(id / nrow_plot)
        layout[row, col] <- id
        gs[[id]] <- plot_ABC_inference(
            model_rf, mini_obs, data_rf,
            protocol = "arm",
            ###
            highlight_values = c(true_para, post_mode, post_mean, post_median),
            highlight_colors = c("black", "#d03333", "#3939ae", "#3da53d"),
            highlight_linetype = c("solid", "dashed", "dashed", "dashed"),
            ###
            fontsize = 20,
            main = list_parameters$Title[para],
            plot_ABC_prior_as_uniform = plot_ABC_prior_as_uniform,
            cutbound = TRUE,
            para_lower_bound = as.numeric(list_parameters$Lower_bound[para]),
            para_upper_bound = as.numeric(list_parameters$Upper_bound[para])
        )
        end_time <- Sys.time()
        cat(paste0("Best parameter: ", post_mode, "\n"))
        print(end_time - start_time)
        #   Clear memory
        model_rf <- c()
        data_rf <- c()
        post_rf <- c()
        post_mode <- c()
    }
    #   Output the csv of parameter values
    write.csv(list_parameters_output, paste0(library_name, "_para_output.csv"))
    #   Plot the prior, posterior and chosen best parameter for all variables
    filename <- paste0(library_name, "_ABC_all.jpeg")
    jpeg(filename, width = 3000, height = 1500)
    p <- grid.arrange(grobs = gs, layout_matrix = layout)
    print(p)
    dev.off()
}

#' @export
sensitivity_library_statistics <- function(library_name,
                                           library_sensitivity_name,
                                           model_variables,
                                           sensitivity_parameter,
                                           sensitivity_values,
                                           list_parameters,
                                           list_targets_library,
                                           ABC_simcount_start = NULL,
                                           ABC_simcount,
                                           compute_parallel = TRUE,
                                           n_cores = NULL,
                                           cn_data_sc = NULL,
                                           cn_data_bulk = NULL,
                                           cn_table = NULL,
                                           arm_level = FALSE) {
    if (sensitivity_parameter == "ABC_simcount") {
        #   Load large rda with all simulations
        #   Save smaller rda with requested simulation counts
        ABC_input_mini <- list()
        for (ABC_simcount in sensitivity_values) {
            load(paste0(library_name, "_ABC_input.rda"))
            ABC_input_mini$model_variables <- ABC_input$model_variables
            ABC_input_mini$sim_param <- ABC_input$sim_param[1:ABC_simcount, ]
            stat_IDs <- names(ABC_input$sim_stat)
            ABC_input_mini$sim_stat <- list()
            for (stat in stat_IDs) {
                ABC_input_mini$sim_stat[[stat]] <- ABC_input$sim_stat[[stat]][1:ABC_simcount, ]
            }
            ABC_input <- ABC_input_mini
            ABC_input_mini <- list()
            filename <- paste0(library_sensitivity_name, "_", ABC_simcount, "_ABC_input.rda")
            save(ABC_input, file = filename)
            ABC_input <- list()
        }
    } else if (sensitivity_parameter == "N_data_sc") {
        N_data_bulk <- as.numeric(names(ground_truth_cn_data_bulk_whole))
        for (N_data_sc in sensitivity_values) {
            cn_data_sc <- ground_truth_cn_data_sc_whole[[as.character(N_data_sc)]]
            cn_data_bulk <- ground_truth_cn_data_bulk_whole[[1]]
            each_library_sensitivity_name <- paste0(library_sensitivity_name, "_", as.character(N_data_sc))
            library_statistics(
                library_name = library_name,
                library_sensitivity_name = each_library_sensitivity_name,
                model_variables = model_variables,
                list_parameters = list_parameters,
                list_targets_library = list_targets_library,
                ABC_simcount_start = ABC_simcount_start,
                ABC_simcount = ABC_simcount,
                compute_parallel = compute_parallel,
                n_cores = n_cores,
                cn_data_sc = cn_data_sc,
                cn_data_bulk = cn_data_bulk,
                cn_table = cn_table,
                arm_level = arm_level
            )
        }
    } else if (sensitivity_parameter == "N_data_bulk") {
        N_data_sc <- as.numeric(names(ground_truth_cn_data_sc_whole))
        for (N_data_bulk in sensitivity_values) {
            cn_data_bulk <- ground_truth_cn_data_bulk_whole[[as.character(N_data_bulk)]]
            cn_data_sc <- ground_truth_cn_data_sc_whole[[1]]
            each_library_sensitivity_name <- paste0(library_sensitivity_name, "_", as.character(N_data_bulk))
            library_statistics(
                library_name = library_name,
                library_sensitivity_name = each_library_sensitivity_name,
                model_variables = model_variables,
                list_parameters = list_parameters,
                list_targets_library = list_targets_library,
                ABC_simcount_start = ABC_simcount_start,
                ABC_simcount = ABC_simcount,
                compute_parallel = compute_parallel,
                n_cores = n_cores,
                cn_data_sc = cn_data_sc,
                cn_data_bulk = cn_data_bulk,
                cn_table = cn_table,
                arm_level = arm_level
            )
        }
    }
}

#' @export
sensitivity_fitting_and_plotting <- function(library_name,
                                             library_sensitivity_name,
                                             sensitivity_title,
                                             sensitivity_parameter,
                                             sensitivity_values,
                                             Error_targets = c("CNA_probability", "Selection_rate"),
                                             Error_metrics = c("Variance", "Standard deviation", "RMSE"),
                                             #  Error_titles = c("All parameters", "Prob(misseg)", "Sel. rates"),
                                             copynumber_DATA = NULL,
                                             parameters_truth = NULL,
                                             list_parameters,
                                             list_targets_by_parameter = NULL,
                                             n_cores = NULL,
                                             plot_ABC_prior_as_uniform = FALSE,
                                             fontsize = 50,
                                             fitting = FALSE) {
    library(ggplot2)
    library(tidyr)
    #   Compute matrix of Error
    list_Error <- list()
    for (Error_target in Error_targets) {
        list_Error[[Error_target]] <- data.frame(Value = sensitivity_values)
        for (Error_metric in Error_metrics) {
            list_Error[[Error_target]][[Error_metric]] <- NA
        }
    }

    for (sensitivity_value in sensitivity_values) {
        if (sensitivity_parameter == "ABC_simcount") {
            copynumber_DATA_to_use <- copynumber_DATA
        } else {
            copynumber_DATA_to_use <- copynumber_DATA[[as.character(sensitivity_value)]]
        }
        #   ABC fitting for each sensitivity value
        library_name_mini <- paste0(library_sensitivity_name, "_", sensitivity_value)
        if (fitting == TRUE) {
            fitting_parameters(
                library_name = library_name_mini,
                copynumber_DATA = copynumber_DATA_to_use,
                parameters_truth = parameters_truth,
                list_parameters = list_parameters,
                list_targets_by_parameter = list_targets_by_parameter,
                n_cores = n_cores,
                plot_ABC_prior_as_uniform = plot_ABC_prior_as_uniform
            )
        }
        #   Input the csv of parameter values
        filename <- paste0(library_name_mini, "_para_output.csv")
        list_parameters_output_mini <- read.csv(filename, header = TRUE)
        # print(list_parameters_output_mini)
        #   Compute Error
        for (j in 1:length(Error_targets)) {
            for (k in 1:length(Error_metrics)) {
                if (Error_metrics[k] == "Standard deviation") {
                    list_Error[[Error_targets[j]]][[Error_metrics[k]]][which(list_Error[[Error_targets[j]]]$Value == sensitivity_value)] <-
                        mean(list_parameters_output_mini$Sd[which(list_parameters_output_mini$Type == Error_targets[j])])
                    print(list_Error[[Error_targets[j]]][[paste0(Error_targets[j], "_", Error_metrics[k])]])
                } else if (Error_metrics[k] == "RMSE") {
                    list_Error[[Error_targets[j]]][[Error_metrics[k]]][which(list_Error[[Error_targets[j]]]$Value == sensitivity_value)] <- compute_error(
                        results = list_parameters_output_mini[which(list_parameters_output_mini$Type == Error_targets[j]), ],
                        ID_actual = "Ground_truth",
                        ID_predicted = "Mean"
                    )
                } else if (Error_metrics[k] == "Variance") {
                    list_Error[[Error_targets[j]]][[Error_metrics[k]]][which(list_Error[[Error_targets[j]]]$Value == sensitivity_value)] <-
                        (mean(list_parameters_output_mini$Sd[which(list_parameters_output_mini$Type == Error_targets[j])]))^2
                }
            }
        }
    }

    # if (sensitivity_parameter == "ABC_simcount") {
    #     list_Error <- data.frame(cbind(
    #         Value = list_Error$CNA_probability$Value,
    #         CNA_RMSE = list_Error$CNA_probability$RMSE,
    #         Sel_RMSE = list_Error$Selection_rate$RMSE
    #     ))
    #     print(list_Error)
    #     cols <- c("c1" = "#ff00ff", "c2" = "#3399ff")
    #     p <- ggplot(list_Error, aes(Value)) +
    #         xlab("Simulation count in ABC library") +
    #         ylab("") +
    #         geom_line(aes(y = Sel_RMSE, color = "Sel.rate"), size = 1) +
    #         geom_point(aes(y = Sel_RMSE, color = "Sel.rate"), size = 10) +
    #         geom_line(aes(y = CNA_RMSE, color = "log10(prob_misseg)"), size = 1) +
    #         geom_point(aes(y = CNA_RMSE, color = "log10(prob_misseg)"), size = 10) +
    #         theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
    #         theme(text = element_text(size = fontsize)) +
    #         scale_color_manual(values = c("Sel.rate" = "#ff00ff", "log10(prob_misseg)" = "#3399ff")) +
    #         theme(legend.position = "top", legend.justification = "left", legend.title = element_blank(), legend.text = element_text(size = fontsize))
    #     plot_name <- paste0(library_sensitivity_name, ".jpeg")
    #     jpeg(plot_name, width = 2000, height = 1000)
    #     print(p)
    #     dev.off()
    # } else {
    #   Plot Error with respect to sensitivity parameter
    for (i in 1:length(Error_targets)) {
        Error_target <- Error_targets[i]
        # Define a custom color palette
        custom_colors <- c("Variance" = "blue", "Standard deviation" = "#0084ff", "RMSE" = "#ff8c00")
        # Reshape the DataFrame into long format
        df_long <- pivot_longer(list_Error[[Error_target]], cols = -Value, names_to = "Variable", values_to = "Value2")
        # Create the ggplot and map the colors using scale_color_manual
        p <- ggplot(df_long, aes(x = Value, y = Value2, color = Variable)) +
            xlab(sensitivity_title) +
            ylab("") +
            theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
            theme(text = element_text(size = fontsize)) +
            geom_line(size = 1) +
            geom_point(size = 10) +
            scale_color_manual(values = custom_colors) +
            theme(legend.position = "top", legend.justification = "left", legend.title = element_blank(), legend.text = element_text(size = fontsize))
        plot_name <- paste0(library_sensitivity_name, Error_target, ".jpeg")
        jpeg(plot_name, width = 2000, height = 1000)
        print(p)
        dev.off()
    }
    # }
}

#' @export
statistics_fitting_and_plotting <- function(library_name,
                                            library_statistics_name,
                                            statistics_title,
                                            statistics_values,
                                            statistics_IDs,
                                            Error_targets = c("CNA_probability", "Selection_rate"),
                                            Error_IDs = c("CNA probability", "Selection rate"),
                                            #  Error_metrics = c("Var", "Sd", "RMSE"),
                                            #  Error_titles = c("All parameters", "Prob(misseg)", "Sel. rates"),
                                            copynumber_DATA = NULL,
                                            parameters_truth = NULL,
                                            list_parameters,
                                            list_targets_library = NULL,
                                            list_targets_misseg = NULL,
                                            list_targets_selection = NULL,
                                            n_cores = NULL,
                                            plot_ABC_prior_as_uniform = FALSE,
                                            fontsize = 50,
                                            fitting = FALSE,
                                            plot_name) {
    library(ggplot2)
    library(tidyr)
    list_Error <- data.frame(matrix(ncol = 5, nrow = 0))
    colnames(list_Error) <- c("Error_Targets", "Statistics_values", "Error", "Error_Titles", "Target_Titles")
    for (stat_value in statistics_values) {
        library_name_mini <- paste0(library_statistics_name, "_", stat_value)
        if (fitting == TRUE) {
            list_targets <- data.frame(matrix(0, ncol = (length(list_targets_library) + 1), nrow = length(list_parameters$Variable)))
            colnames(list_targets) <- c("Variable", list_targets_library)
            list_targets[, 1] <- list_parameters$Variable
            for (row in 1:nrow(list_targets)) {
                if (any(grep("Green", stat_value))) {
                    list_targets[row, which(colnames(list_targets) %in%
                        list_targets_library[grep(
                            "variable=(average_CN|clonal_CN|event_count|shannon)",
                            list_targets_library
                        )])] <- 1
                } else if (any(grep("Blue", stat_value))) {
                    list_targets[row, which(colnames(list_targets) %in%
                        list_targets_library[grep(
                            "variable=(cherries|pitchforks|IL_number|avgLadder)",
                            list_targets_library
                        )])] <- 1
                } else if (any(grep("Red", stat_value))) {
                    list_targets[row, which(colnames(list_targets) %in%
                        list_targets_library[grep(
                            "variable=(stairs|colless|sackin|B2)",
                            list_targets_library
                        )])] <- 1
                } else if (stat_value == "All") {
                    list_targets[row, which(colnames(list_targets) %in%
                        list_targets_library)] <- 1
                }
            }
            fitting_parameters(
                rda_name = library_name,
                library_name = library_name_mini,
                copynumber_DATA = copynumber_DATA,
                parameters_truth = parameters_truth,
                list_parameters = list_parameters,
                list_targets_by_parameter = list_targets,
                n_cores = n_cores,
                plot_ABC_prior_as_uniform = plot_ABC_prior_as_uniform
            )
        }
        #   Input the csv of parameter values
        filename <- paste0(library_name_mini, "_para_output.csv")
        list_parameters_output_mini <- read.csv(filename, header = TRUE)
        #   Compute Error
        for (Error_target in Error_targets) {
            list_Error[nrow(list_Error) + 1, ] <- c(
                Error_target, stat_value, round(compute_error(
                    results = list_parameters_output_mini[which(list_parameters_output_mini$Type == Error_target), ],
                    ID_actual = "Ground_truth",
                    ID_predicted = "Mean"
                ), digits = 3),
                statistics_IDs[which(statistics_values == stat_value)],
                Error_IDs[which(Error_targets == Error_target)]
            )
        }
    }

    list_Error$Error <- as.numeric(list_Error$Error)

    list_Error$Error_Targets <- factor(list_Error$Error_Targets, levels = Error_targets)
    list_Error$Error_Titles <- factor(list_Error$Error_Titles, levels = statistics_IDs)

    list_Error <<- list_Error
    print(list_Error)
    p <- ggplot(data = list_Error, aes(x = Error_Titles, y = Error, fill = Target_Titles)) +
        geom_bar(stat = "identity", position = position_dodge(), width = 0.6) +
        scale_fill_manual(values = c("#774d28", "#70c972"))
    # scale_fill_manual(values = c("" = "#dd4751", "Selection_rate" = "darkblue", alpha = 0.3))
    p <- p +
        xlab(statistics_title) +
        ylab("RMSE") +
        # scale_y_continuous() +
        theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
        theme(text = element_text(size = fontsize)) +
        theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
        theme(legend.position = "top", legend.justification = "left", legend.title = element_blank(), legend.text = element_text(size = fontsize)) +
        jpeg(plot_name, width = 1500, height = 1500)
    print(p)
    dev.off()
}
