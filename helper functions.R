valid_sites <- c(
  "BD 001", "BN 001", "CN 001", "CN 002", "HK 000", "ID 004", "ID 005", "ID 006",
  "IN 004", "IN 007", "JP 001", "KH 002", "MG 001", "MG 000", "MY 001", "MY 002",
  "MY 003", "MY 004", "MY 005", "NP 001", "PH 001", "PH 002", "PK 001", "PK 002",
  "SG 002", "PK 003", "SG 001", "SL 001", "TH 001", "TH 002", "TH 004", "TH 005",
  "TL 001", "TW 001", "VN 001", "VN 004", "VN 005"
)



select_sites <- function(num_sites, specification = "random_selection") {
  
  if (specification == "random_selection") {
    #print("inside")
    if (num_sites > length(valid_sites)) {
      #print("error")
      stop("Error: Number of sites requested exceeds number of available sites.")
    }
    return(sample(valid_sites, num_sites, replace = FALSE))
  }
  
  if (specification == "different_country") {
    # Extract country codes (first 2 characters)
    country_codes <- substr(valid_sites, 1, 2)
    country_site_map <- split(valid_sites, country_codes)
    #print(length(country_site_map))
    if (num_sites > length(country_site_map)) {
      stop("Error: Number of sites requested exceeds number of distinct countries.")
    }
    
    selected_countries <- sample(names(country_site_map), num_sites)
    selected_sites <- sapply(selected_countries, function(cty) {
      sample(country_site_map[[cty]], 1)
    })
    
    return(unname(selected_sites))
  }
  
  stop("Error: Unknown specification. Use 'random_selection' or 'different_country'.")
}

#select_sites(10, "different_country")
df <- read.csv("final_cleaned_data.csv")
site_table <- table(df$site)

generate_site_counts_per_site <- function(site_list, total_sample) {  
  # Initialize the vector to store counts from the table
  site_counts <- numeric(length(site_list))
  # Retrieve the count for each site using a for loop
  for (i in seq_along(site_list)) {
    site_name <- site_list[i]
    if (site_name %in% names(site_table)) {
      site_counts[i] <- site_table[site_name]
    } else {
      stop(paste("Site name", site_name, "not found in site_table."))
    }
  }
  #print(site_counts)
  
  # Calculate the total number of patients across the specified sites
  total_site_count <- sum(site_counts)
  #print(total_site_count)
  
  # Calculate the allocated samples per site and round up
  allocated_samples <- ceiling(total_sample * site_counts / total_site_count)
  
  return(allocated_samples)
}


sites_to_allocate <- c("VN 001", "SG 001", "VN 005")
total_sample_size <- 300

allocated_counts <- generate_site_counts_per_site(sites_to_allocate, total_sample_size)
print(allocated_counts)


#change this function to test different distributions of mean
generate_trial_effect <- function(num_sites) {
  #mu <- 0.0
  #sd <- sqrt(0.0033333)
  #Draw site-specific treatment effects from a normal distribution
  #return (rnorm(n = num_sites, mean = mu, sd = sd))
  #return (generate_left_skewed_effect(num_sites))
  #return(runif(num_sites, min = 0, max = 0.2))
  return(rep(0.1,num_sites))
  
  #shape <- 3
  #scale <- 0.033333
  #return(rgamma(num_sites, shape = shape, scale = scale))
  
}

#test
#hist(generate_trial_effect(1000))


generate_left_skewed_effect <- function(n, mu = 0.1, shape = 0.8, rate = 20) {
  effects <- numeric(n)
  i <- 1
  while (i <= n) {
    x <- rgamma(1, shape = shape, rate = rate)
    val <- mu - x
    if (val >= 0) {
      effects[i] <- val
      i <- i + 1
    }
    # else: discard and resample
  }
  return(effects)
}



calculate_condition_number <- function(df) {
  # Create the design matrix from the model formula
  X <- model.matrix(death ~ experiment + sex + cal_agey + vap + 
                        comorbidities_CCI + site, data = df)[, -1]  # Remove intercept
  
  # Compute the singular values of X
  svd_vals <- svd(X)$d
  
  # Calculate the condition number (ratio of largest to smallest singular value)
  condition_number <- max(svd_vals) / min(svd_vals)
  
  cat("Condition Number of the Design Matrix:", round(condition_number, 2), "\n")
  
  # General rule of thumb interpretation
  if (condition_number < 100) {
    cat("✅ Multicollinearity is low.\n")
  } else if (condition_number < 1000) {
    cat("⚠️ Moderate multicollinearity.\n")
  } else {
    cat("🚨 High multicollinearity! Consider variable reduction or regularization.\n")
  }
  
  return(condition_number)
}

mean_mortality_by_site <- function(df) {
  site_means <- aggregate(mortality ~ site, data = df, FUN = mean)
  mean(site_means$mortality)
}


calculate_total_sample_size  <- function(
    p_control,           # Control group mortality
    p_experiment,        # Experimental group mortality
    tau_sq,              # Between-site variance in treatment effect (τ²)
    J,                   # Number of sites
    alpha = 0.05,        # Type I error rate
    power = 0.9          # Power (1 - beta)
) {
  # Standard normal quantiles
  z_alpha <- qnorm(1 - alpha / 2)
  z_beta <- qnorm(power)
  
  # Treatment effect
  delta <- abs(p_control - p_experiment)
  
  # Variance term with between-site variance adjustment
  var_combined <- p_control * (1 - p_control) +
    p_experiment * (1 - p_experiment) +
    tau_sq / J
  print(combined <- p_control * (1 - p_control) +
          p_experiment * (1 - p_experiment))
  # Sample size calculation
  n <- ((z_alpha + z_beta)^2 * var_combined) / (delta^2)
  print((z_alpha + z_beta))
  return(ceiling(2*n))  # Round up to next whole number
}


calculate_total_sample_size(
  p_control = 0.41,
  p_experiment = 0.31,
  tau_sq = 0,
  J = 1
)



calculate_total_sample_size  <- function(
    p_control,           # Control group mortality
    p_experiment,        # Experimental group mortality
    sigma_sq,              # Between-site variance in treatment effect (τ²)
    s,                   # Number of sites
    alpha = 0.05,        # Type I error rate
    power = 0.9          # Power (1 - beta)
) {
  # Standard normal quantiles
  z_alpha <- qnorm(1 - alpha / 2)
  z_beta <- qnorm(power)
  z_term <- (z_alpha + z_beta)
  p_term <- p_control * (1 - p_control) + p_experiment * (1 - p_experiment)
  
  # Treatment effect
  delta <- abs(p_control - p_experiment)
  
  # Sample size calculation
  numerator <- p_term * (z_term^2)
  denominator <- delta - ((z_term*sqrt(sigma_sq))/sqrt(s))
  print(p_term)
  print(z_term)
  print(denominator)
  n <- numerator / (denominator^2)
  return(ceiling(2*n))  # Round up to next whole number
}



calculate_total_sample_size <- function(
    p_control,           # Control group mortality (p0)
    p_experiment,        # Experimental group mortality (p1 = p0 - mu)
    sigma_sq,            # Between-site variance in treatment effect (tau^2)
    s,                   # Number of sites
    alpha = 0.05,        # Type I error rate
    power = 0.9          # Desired power (1 - beta)
) {
  # Calculate z-scores for alpha and beta
  z_alpha <- qnorm(1 - alpha / 2)
  z_beta <- qnorm(power)
  z_term <- z_alpha + z_beta
  
  # Treatment effect
  mu <- p_control - p_experiment
  
  # V_within formula
  v_within <- 2 * p_control * (1 - p_control) + 
    (1 - 2 * p_control) * mu - 
    mu^2 - sigma_sq
  
  # Check for invalid conditions
  print("v within is")
  print(v_within)
  if (abs(mu) <= z_term * sqrt(sigma_sq / s)) {
    stop("Effect size too small compared to site heterogeneity; denominator becomes non-positive.")
  }
  
  # Sample size formula
  n <- (z_term^2 * 2 * v_within) / (abs(mu) - z_term * sqrt(sigma_sq / s))^2
  
  return(n)
}


calculate_total_sample_size(
  p_control = 0.414,
  p_experiment = 0.1,
  sigma_sq = 0,
  s = 1,
  alpha = 0.05,
  power = 0.9
)

p1 <- 0.414
p2 <- 0.314
p1 * (1 - p1) + p2 * (1 - p2)

2*p1*(1-p1) -0.1+2*p1*0.1-0.01

2*p1*(1-p1) -0.1*(1-2*p1)- (0.01)











compute_g_difference <- function(mean_p, var_p, delta, S = 10000) {
  # Function g_delta(p)
  g_delta <- function(p, delta) {
    p * (1 - p) + (p - delta) * (1 - p + delta)
  }
  
  # Generate site-level p_s with given mean and variance
  sd_p <- sqrt(var_p)
  p_s <- rnorm(S, mean = mean_p, sd = sd_p)
  
  # Clip p_s to [0, 1] since it's a probability
  p_s <- pmin(pmax(p_s, 0), 1)
  
  # Compute site-level g values
  g_vals <- g_delta(p_s, delta)
  g_bar <- g_delta(mean(p_s), delta)  # empirically compute g_delta(mean of p_s)
  
  # Return the difference
  return(mean(g_vals) - g_bar)
}

compute_g_difference(mean_p = 0.4152047, var_p = 0.1662231, delta = 0.1) 

  
  
  
# Example baseline linear predictor (fixed part): eta = -1.2
logit_sd_to_prob_sd <- function(logit_sd, baseline_prob, n = 10000) {
  eta <- qlogis(baseline_prob)  # logit(p)
  u <- rnorm(n, mean = 0, sd = logit_sd)
  p <- plogis(eta + u)
  return(sd(p))  # std dev of mortality probabilities across sites
}

# Your case:
logit_sd <- 0.7582
baseline_prob <- 0.4152047

logit_sd_to_prob_sd(logit_sd, baseline_prob)

