library(readxl)
library(writexl)
library(dplyr)

# Load the Excel files
#  Update these paths to your actual file locations
chemical_file <- "C:/Users/91962/All_319/chemical_compounds.xlsx"
d2_file <- "C:/Users/91962/df2_final.xlsx"

chemical_data <- read_xlsx(chemical_file)
d2_data <- read_xlsx(d2_file)

# Check and clean up column names
# Make sure both files have a "Preferred.Name" column
print(names(chemical_data))
print(names(d2_data))

# Standardize column for matching
# Rename to a common column for joining
chemical_data <- chemical_data %>% rename(Name = `PREFERRED.NAME`)
d2_data <- d2_data %>% rename(Name = `PREFERRED.NAME`)

# Trim any spaces just in case
chemical_data$Name <- trimws(chemical_data$Name)
d2_data$Name <- trimws(d2_data$Name)

# Remove matching compounds from d2_data
filtered_data <- d2_data %>%
  filter(!(Name %in% chemical_data$Name))

#new Excel file
output_dir <- "C:/Users/91962"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Save the filtered Excel file to that path
write_xlsx(filtered_data, file.path(output_dir, "Filtered_compounds.xlsx"))

cat("Done! Saved filtered dataset to:", file.path(output_dir, "Filtered_compounds.xlsx"), "\n")


------------------------------------------------------------------------------------------

# Set your folder paths
folder_a <- "C:/Users/91962/100/chemical_compounds.xlsx"
folder_b <- "C:/Users/91962/100/Filtered_compounds.xlsx"

#List file names (filenames only, not full paths)
files_a <- list.files(folder_a)
files_b <- list.files(folder_b)

#Find files unique to each folder
only_in_a <- setdiff(files_a, files_b)
only_in_b <- setdiff(files_b, files_a)

#results
cat(" Files only in Folder A:\n")
print(only_in_a)

cat("\n Files only in Folder B:\n")
print(only_in_b)

----------------------------------------------------------------------------
#With set  
data <- read_excel("C:/Users/91962/All_319/Filtered_compounds.xlsx")

# Set seed for reproducibility
set.seed(123)

#a new column 'Set' with 80% Train, 20% Test
data$Set <- ifelse(runif(nrow(data)) < 0.8, "Train", "Test")

write_xlsx(data, "C:/Users/91962/All_319/Filtered_compounds_with_Set1.xlsx")