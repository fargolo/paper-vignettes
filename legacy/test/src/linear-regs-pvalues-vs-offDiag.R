# Get a list of CSV files in your directory
# setwd("outputs")
csv_files <- list.files(pattern = ".csv")

# Create an empty dataframe to store results
results_df <- data.frame(File = character(), Label = character(), 
                         Intercept = numeric(), Coeff = numeric(),
                         PValue = numeric(), 
                         AdjustedRSquare = numeric(), 
                         stringsAsFactors = FALSE)

# Iterate through each CSV file
for (file in csv_files) {
  # Read the CSV file into a dataframe
  df <- read.csv(file)
  
  # Get unique labels in the dataframe
  unique_labels <- unique(df$Label)
  
  # Iterate through each label and perform linear regression
  for (label in unique_labels) {
    # Subset the dataframe for the current label
    subset_df <- subset(df, Label == label)
    
    # Perform linear regression for "windows" vs. "average"
    lm_result <- lm(Average ~ Window, data = subset_df)
    
    # Extract relevant information from the model summary
    intercept <- coef(lm_result)[1]
    lin_coeff <- coef(lm_result)[2]
    p_value <- summary(lm_result)$coefficients[2, 4]
    adjusted_r_square <- summary(lm_result)$adj.r.squared
    
    # Store the results in the dataframe
    results_df <- rbind(results_df, data.frame(File = file, Label = label, 
                                               Intercept = intercept, 
                                               Coeff = lin_coeff, 
                                               PValue = p_value, 
                                               AdjustedRSquare = adjusted_r_square))
  }
}

# Write the results dataframe to a CSV file
write.csv(results_df, "linear_regression_results.csv", row.names = FALSE)
