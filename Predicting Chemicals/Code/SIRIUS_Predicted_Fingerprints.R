library(readr)
library(dplyr)
library(tibble)
library(writexl)

# Set path to SIRIUS project folder
sirius_folder <- "C:/Users/91962/Predicting Chemicals/validation compound/"

# List compound folders
compound_folders <- list.dirs(sirius_folder, recursive = FALSE)

# Initialize result list
results <- list()

# Loop over each compound folder
for (folder in compound_folders) {
  fingerprint_file <- file.path(folder, "fingerid", "fingerprint.csi.fingerprint")
  
  if (file.exists(fingerprint_file)) {
    lines <- read_lines(fingerprint_file)
    
    compound_name <- basename(folder)
    fp_vector <- lines[grepl("^Un", lines)]  # Only keep fingerprint lines
    
    # Extract names and values
    values <- sapply(fp_vector, function(line) {
      parts <- unlist(strsplit(line, "="))
      if (length(parts) == 2) {
        val <- as.numeric(trimws(parts[2]))
        ifelse(val > 0.5, 1, 0)
      } else {
        NA
      }
    })
    
    names(values) <- sapply(fp_vector, function(line) {
      parts <- unlist(strsplit(line, "="))
      trimws(parts[1])
    })
    
    results[[compound_name]] <- values
  }
}

# Convert to data frame
df <- bind_rows(results, .id = "Compound")

# Write to Excel
write.csv(df, "SIRIUS_Predicted_Fingerprints.csv")
