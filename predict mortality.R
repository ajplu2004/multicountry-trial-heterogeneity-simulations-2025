
predict_mortality_from_data <- function(df, model_dir = "models") {
  predictions <- numeric(nrow(df))
  sites <- unique(df$site)
  
  for (site in sites) {
    #cat("Processing site:", site, "\n")
    
    # Get all rows for this site
    rows <- which(df$site == site)
    site_data <- df[rows, ]
    
    # Build file path
    safe_site <- gsub(" ", "_", site)
    model_path <- file.path(model_dir, paste0("fixed_effect_model_", safe_site, ".rds"))
    
    # Load model
    model <- readRDS(model_path)
    
    # Predict
    pred <- predict(model, newdata = site_data, type = "response")
    predictions[rows] <- pred
  }
  
  return(predictions)
}
#  # Create a small dummy dataframe
#  df_test <- data.frame(
#    cal_agey = c(67, 54, 81),
#    sex = c("F", "M", "F"),
#    vap = c(0, 1, 1),
#    hpd_admreason = c("GIT", "GIT", "OTH"),
#    comorbidities_CCI = c(1, 2, 0),
#    site = c("BD 001", "BD 001", "HK 000")
#  )
  
  # Run prediction
#  preds <- predict_mortality_from_data(df_test)
  
  # Show result
#  print(preds)


