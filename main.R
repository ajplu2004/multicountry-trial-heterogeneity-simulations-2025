#load in files
source("simulate trial.R")

#####################run simulations
run_main_simulation <- function(
    total_sample_size = 964,
    max_sites = 5,
    iterations = 10,
    generation_method = "distribution",  # or "bootstrap"
    treatment_proportion = 0.5 
) {
  results_list <- list()
  
  for (num_sites in 2:max_sites) {
    cat("Running simulation with", num_sites, "sites...\n")
    
    # Run full set of simulations for this site count
    trial_results <- simulate_multiple_trials(
      num_sites = num_sites,
      specification = "random_selection",
      iterations = iterations,
      total_sample_size = total_sample_size,
      site_counts_per_site = "same",
      treatment_proportions = treatment_proportion,
      generation_method = generation_method
    )
    
    # Store as a row with two columns: num_sites, simulation_data
    results_list[[num_sites]] <- data.frame(
      num_sites = num_sites,
      stringsAsFactors = FALSE
    )
    results_list[[num_sites]]$simulation_data <- list(trial_results)
  }
  
  # Combine all rows into a single dataframe
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  
  return(results_df)
}

###################running and saving the simulation########################
results_df <- run_main_simulation(total_sample_size = 1032,
                                  max_sites = 25,
                                  iterations = 2,
                                  generation_method = "distribution",  # or "bootstrap"
                                  treatment_proportion = 0.5 
                                  )

colnames(results_df)
colnames(results_df$simulation_data[[1]]) #output of site number
results_df$simulation_data[[1]]$z_test[1] #output of the first set number's first repeat of simulation
saveRDS(results_df, file = "results_df_pw_naive_15_85.rds")
results_df <- readRDS("results_df_pw.rds")



##############plotting#####################
library(ggplot2)
library(dplyr)
library(tidyr)

plot_test_power_by_site <- function(results_df) {
  # Extract test column names to summarize
  test_names <- c(
    "z_test",
    "naive_t_test",
    "fix_effect_model_test",
    "mixed_effect_model_test",
    "mantel_haenszel_test"
  )
  
  # Prepare summary dataframe: one row per (site, test)
  summary_df <- lapply(seq_len(nrow(results_df)), function(i) {
    site_num <- results_df$num_sites[i]
    sim_data <- results_df$simulation_data[[i]]
    
    means <- sapply(test_names, function(test) {
      mean(sim_data[[test]], na.rm = TRUE)
    })
    
    data.frame(num_sites = site_num, t(means))
  }) %>% bind_rows()
  
  # Convert to long format for ggplot
  plot_data <- summary_df %>%
    pivot_longer(
      cols = all_of(test_names),
      names_to = "test_type",
      values_to = "mean_result"
    )
  
  # Plot using ggplot2
  ggplot(plot_data, aes(x = num_sites, y = mean_result, color = test_type)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    labs(
      title = "Mean Power of Each Test by Number of Sites",
      x = "Number of Sites",
      y = "Proportion of Simulations Passed (Power)",
      color = "Test Type"
    ) +  ylim(0.75, 1)+
    theme_minimal() +
    theme(text = element_text(size = 12))
}

plot_test_power_by_site(results_df)

