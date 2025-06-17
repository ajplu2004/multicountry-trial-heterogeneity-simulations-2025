# Load libraries
library(readr)
library(dplyr)

#read in data
df <- read_csv("final_cleaned_data.csv")
full_predictors <- c("cal_agey", "sex", "vap", "hpd_admreason", "comorbidities_CCI")

# Create output folder if it doesn't exist
if (!dir.exists("models")) {
  dir.create("models")
}

# Get list of unique sites
sites <- unique(df$site)

# Loop through each site
for (site in unique(df$site)) {
  #print(site)
  safe_site <- gsub(" ", "_", site)
  model_name <- paste0("fixed_effect_model_", safe_site)
  model_path <- file.path("models", paste0(model_name, ".rds"))
  
  site_data <- df %>% filter(site == !!site)
  
  # Ensure factors are set
  site_data <- site_data %>%
    mutate(
      sex = factor(sex),
      hpd_admreason = factor(hpd_admreason)
    )
  
  # Remove predictors with no variation
  valid_predictors <- full_predictors[sapply(site_data[, full_predictors], function(x) length(unique(x)) > 1)]
  
  # Try decreasing subsets of predictors until model fits
  formula_text <- paste("mortality ~", paste(valid_predictors, collapse = " + "))
  #print(formula_text)
  model <- tryCatch({
    suppressWarnings(glm(as.formula(formula_text),
                         data = site_data,
                         family = binomial(link = "logit"),
                         control = list(maxit = 100)))
  }, error = function(e) NULL)
  saveRDS(model, model_path)
  assign(model_name, model, envir = .GlobalEnv)
}


##############test function########################
# Define dummy patient input
new_patient <- list(
  cal_agey = 107,
  sex = "F",
  vap = 1,
  hpd_admreason = "OTH",
  comorbidities_CCI = 1
)

# Convert to data frame
new_data <- as.data.frame(new_patient, stringsAsFactors = FALSE)

# Get all model file names from models directory
model_files <- list.files("models", pattern = "\\.rds$", full.names = TRUE)

# Loop over all models
for (file in model_files) {
  model_name <- basename(file)
  print(model_name)
  # Read model
  model <- readRDS(file)
  print(predict(model, new_patient, type = "response"))
}


