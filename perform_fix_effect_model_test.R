
perform_fix_effect_model_test <- function(df, alpha = 0.05) {
  # Ensure variables are of correct type
  df$death <- as.numeric(df$death)
  df$experiment <- as.numeric(df$experiment)
  df$site <- as.factor(df$site)
  
  # Fit fixed-effects logistic regression
  model <- glm(death ~ experiment + site, data = df, family = binomial())
  
  # Extract treatment coefficient and p-value
  coef_summary <- summary(model)$coefficients
  treat_coef <- coef_summary["experiment", "Estimate"]
  treat_pval <- coef_summary["experiment", "Pr(>|z|)"]
  
  # Check if effect is negative and p-value < alpha
  result <- (treat_coef < 0) & ((treat_pval/2) < alpha)
  #print(result)
  #return(list(
  #  passed = result,
  #  coef = treat_coef,
  #  p_value = treat_pval,
  #  model_summary = coef_summary
  #))
  return (as.numeric(result))
}

#outcome <- perform_fix_effect_model_test(sim_data)
