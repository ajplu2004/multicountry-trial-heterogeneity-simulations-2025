library(dplyr)
library(readr)

df <- read_csv("final_cleaned_data.csv")

learn_site_distributions <- function(df) {
  sites <- unique(df$site)
  site_distributions <- list()
  
  for (site in sites) {
    site_data <- df %>% filter(site == !!site)
    
    # --- Binary variables ---
    p_sex_M <- mean(site_data$sex == "M", na.rm = TRUE)
    p_vap_1 <- mean(site_data$vap == 1, na.rm = TRUE)
    
    # --- Categorical variable ---
    admin_dist <- prop.table(table(site_data$hpd_admreason))
    admin_probs <- as.list(admin_dist)
    
    # --- Truncated Normal for age (> 0) ---
    age_data <- site_data$cal_agey[site_data$cal_agey >= 0]
    age_mean <- mean(age_data, na.rm = TRUE)
    age_sd <- sd(age_data, na.rm = TRUE)
    
    # --- Zero-Inflated Poisson for CCI ---
    cci_data <- site_data$comorbidities_CCI
    p_zero <- mean(cci_data == 0, na.rm = TRUE)
    lambda <- mean(cci_data[cci_data > 0], na.rm = TRUE)
    
    # Save all parameters into a list for the site
    site_distributions[[site]] <- list(
      p_sex_M = p_sex_M,
      p_vap_1 = p_vap_1,
      admin_probs = admin_probs,
      age_mean = age_mean,
      age_sd = age_sd,
      cci_p_zero = p_zero,
      cci_lambda = lambda
    )
  }
  
  return(site_distributions)
}

# Learn site distributions
#site_params <- learn_site_distributions(df)

# Save to file
#saveRDS(site_params, file = "site_params.rds")

site_params <- readRDS("site_params.rds")

generate_random_dataset_distribution <- function(sites, site_counts) {
  if (length(sites) != length(site_counts)) {
    stop("Length of sites and site_counts must be equal.")
  }
  
  # Truncated normal sampler (> 0)
  rtruncnorm_pos <- function(n, mean, sd) {
  x <- rnorm(n, mean, sd)
  while (any(x < 0)) {
    x[x < 0] <- rnorm(sum(x < 0), mean, sd)
  }
  return(x)
}

  
  all_rows <- list()
  
  for (i in seq_along(sites)) {
    site <- sites[i]
    n <- site_counts[i]
    
    # Retrieve parameters for this site
    params <- site_params[[site]]
    if (is.null(params)) stop(paste("Missing site parameters for", site))
    
    # Generate variables
    cal_agey <- round(rtruncnorm_pos(n, params$age_mean, params$age_sd))
    sex <- ifelse(runif(n) < params$p_sex_M, "M", "F")
    vap <- rbinom(n, 1, params$p_vap_1)
    
    adm_choices <- names(params$admin_probs)
    adm_probs <- unlist(params$admin_probs)
    hpd_admreason <- sample(adm_choices, size = n, replace = TRUE, prob = adm_probs)
    
    is_zero <- rbinom(n, 1, params$cci_p_zero)
    comorbidities_CCI <- ifelse(is_zero == 1, 0, rpois(n, lambda = params$cci_lambda))
    
    site_df <- data.frame(
      cal_agey = cal_agey,
      sex = sex,
      vap = vap,
      hpd_admreason = hpd_admreason,
      comorbidities_CCI = comorbidities_CCI,
      site = site,
      stringsAsFactors = FALSE
    )
    
    all_rows[[i]] <- site_df
  }
  
  return(do.call(rbind, all_rows))
}


# Example: simulate 5 from BD 001, 3 from JP 001
#sites <- c("TH 001")
#site_counts <- c(14400)

# Assume site_params already exists
#synthetic_df <- generate_synthetic_patients(sites, site_counts)
#print(synthetic_df)



# Base plotting
#plot(
#  density(df$cal_agey[df$site == "TH 001"], na.rm = TRUE),
#  main = "Density of cal_agey: Real vs Synthetic (BD 001)",
#  xlab = "cal_agey",
#  ylab = "Density",
#  col = "blue",
#  lwd = 2,
#  ylim = c(0, 0.035)
#)

# Add synthetic line
#lines(
#  density(synthetic_df$cal_agey[synthetic_df$site == "TH 001"], na.rm = TRUE),
#  col = "red",
#  lwd = 2,
#  lty = 2
#)

# Add legend
#legend("topright", legend = c("Real BD 001", "Synthetic BD 001"),
#       col = c("blue", "red"), lty = c(1, 2), lwd = 2)

