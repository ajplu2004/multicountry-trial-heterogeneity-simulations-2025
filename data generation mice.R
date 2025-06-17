# ------------------------------------
# 0. PACKAGES & GLOBAL SETTINGS
# ------------------------------------
library(dplyr)
library(mice)
library(readr)

syn_vars <- c("cal_agey", "sex", "vap", "hpd_admreason", "comorbidities_CCI")

model_out_dir <- "site_cart_models"
if (!dir.exists(model_out_dir)) dir.create(model_out_dir)

# ------------------------------------
# 1. Learn CART models for each site
# ------------------------------------
learn_site_cart_models <- function(data,
                                   site_col = "site",
                                   vars     = syn_vars,
                                   out_dir  = model_out_dir,
                                   m_sets   = 10000,
                                   maxit    = 1) {
  stopifnot(all(c(site_col, vars) %in% names(data)))
  
  unique_sites <- unique(data[[site_col]])
  
  for (s in unique_sites) {
    message("Fitting CART for site: ", s)
    
    site_data <- data %>%
      filter(.data[[site_col]] == s) %>%
      select(all_of(vars)) %>%
      mutate(across(where(is.character), as.factor)) %>%
      mutate(vap = as.factor(vap))
    
    where_mat <- make.where(site_data, "all")
    method_vec <- rep("cart", length(vars))
    names(method_vec) <- names(site_data)
    
    cart_obj <- mice(site_data,
                     m         = m_sets,
                     maxit     = maxit,
                     method    = method_vec,
                     where     = where_mat,
                     printFlag = FALSE)
    
    safe_name <- gsub("[^A-Za-z0-9]", "_", s)
    saveRDS(cart_obj, file = file.path(out_dir,
                                       paste0("cartmodel_", safe_name, ".rds")))
  }
  
  invisible(TRUE)
}


df <- read_csv("final_cleaned_data.csv")

# Step 1: Save site-wise preprocessed data
#learn_site_cart_models(df)


# ------------------------------------
# 2. Generate synthetic data using rotating random m=100 imputations
# ------------------------------------
generate_random_dataset_distribution <- function(sites,
                                                 site_counts,
                                                 site_col  = "site",
                                                 vars      = syn_vars,
                                                 model_dir = model_out_dir) {
  stopifnot(length(sites) == length(site_counts))
  out_rows <- vector("list", length(sites))
  
  for (i in seq_along(sites)) {
    s         <- sites[i]
    n_needed  <- site_counts[i]
    safe_name <- gsub("[^A-Za-z0-9]", "_", s)
    model_f   <- file.path(model_dir, paste0("cartmodel_", safe_name, ".rds"))
    
    if (!file.exists(model_f))
      stop("No CART model saved for site ", s)
    
    cart_obj <- readRDS(model_f)
    m_total  <- cart_obj$m
    block_size <- 10 #change here if we want
    
    reps <- ceiling(n_needed / block_size)
    syn_blocks <- vector("list", reps)
    
    for (j in seq_len(reps)) {
      k <- sample(1:m_total, 1, replace = TRUE)# randomly pick an imputation
      block_data <- complete(cart_obj, action = k)
      
      block_data <- block_data %>%
        mutate(across(where(is.factor), as.character))
      
      if ("vap" %in% names(block_data)) {
        block_data$vap <- as.numeric(block_data$vap)
      }
      
      block_data[[site_col]] <- s
      syn_blocks[[j]] <- block_data %>%
        slice_sample(n = block_size, replace = block_size > nrow(block_data))
    }
    
    syn_data <- bind_rows(syn_blocks) %>%
      slice_head(n = n_needed)
    
    if (s == "SL 001") {
      syn_data$hpd_admreason <- "ONC"
    }
    
    out_rows[[i]] <- syn_data
  }
  
  bind_rows(out_rows)
}



# need 14 400 synthetic rows for TH 001
synthetic_TH <- generate_random_dataset_distribution(
  sites       = "TH 001",
  site_counts = 1000)

# peek
head(synthetic_TH)




library(ggplot2)

# Filter the original data to match TH 001
original_TH <- df %>% filter(site == "TH 001")

# Plot density curves
ggplot() +
  geom_density(data = original_TH,
               aes(x = comorbidities_CCI, color = "Original"),
               linetype = "solid", size = 1.2, na.rm = TRUE) +
  geom_density(data = synthetic_TH,
               aes(x = comorbidities_CCI, color = "Synthetic"),
               linetype = "dashed", size = 1.2, na.rm = TRUE) +
  labs(title = "Density of comorbidities_CCI: Real vs Synthetic (TH 001), 1000 synthetic data",
       x = "comorbidities_CCI",
       y = "Density",
       color = "Data Source") +   
  theme_minimal() +
  scale_color_manual(values = c("Original" = "blue", "Synthetic" = "red")) +
  theme(legend.position = "top")

