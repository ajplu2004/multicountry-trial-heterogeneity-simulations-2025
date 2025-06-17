library(lme4)

perform_mixed_effect_model_test <- function(df, alpha = 0.05) {
  # Ensure correct data types
  df$death <- as.numeric(df$death)
  df$experiment <- as.numeric(df$experiment)
  df$site <- as.factor(df$site)
  
  # Fit random-intercepts logistic regression
  model <- glmer(death ~ experiment + (1 | site), data = df, family = binomial())
  #print(model)
  # Extract fixed effect for treatment
  coef_summary <- summary(model)$coefficients
  treat_coef <- coef_summary["experiment", "Estimate"]
  treat_pval <- coef_summary["experiment", "Pr(>|z|)"]
  
  # Check if effect is negative and p-value < 0.5
  result <- (treat_coef < 0) & ((treat_pval/2) < alpha)
  
  #return(list(
  #  passed = result,
  #  coef = treat_coef,
  #  p_value = treat_pval,
  #  model_summary = coef_summary
  #))
  return(as.numeric(result))
}

#outcome <- perform_mixed_effect_model_test(sim_data)
