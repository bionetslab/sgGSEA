# Load the xCell library
library(xCell)

# Set the directory path
directory_path <- "/data_nfs/je30bery/hiwi/sgGSEA/data/for_R_xcell"

# Get a list of all CSV files in the directory
expression_files <- list.files(path = directory_path, pattern = "\\.csv$", full.names = TRUE)

# Loop through each expression file
for (expr_file in expression_files) {
  # Read the expression matrix
  exprMatrix <- read.table(expr_file, header = TRUE, row.names = 1, as.is = TRUE, sep = ",")

  # Perform xCell analysis
  result <- xCellAnalysis(exprMatrix)

  # Create a unique result file name based on the input expression file
  result_file <- file.path(directory_path, paste0("result_", tools::file_path_sans_ext(basename(expr_file)), ".csv"))

  # Write the result to a CSV file
  write.csv(result, result_file, row.names = TRUE)

  # Print a message indicating the completion for each file
  cat("Analysis completed for:", expr_file, "\n")
}

# Print a final message indicating the completion of the script
cat("All xCell analyses completed for the expression files in the directory.\n")
