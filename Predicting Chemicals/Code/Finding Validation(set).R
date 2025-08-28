# Load libraries
library(readxl)
library(dplyr)
library(writexl)

# Set file paths
file_a_path <- "C:/Users/91962/Predicting Chemicals/D1.csv"  # Large file
file_b_path <- "C:/Users/91962/Predicting Chemicals/Validation.xlsx"  # Smaller file

# Read both Excel files
df_a <- read.csv(file_a_path)
df_b <- read_excel(file_b_path)

# Select only DTXSID and HIT.CALL from File A
df_a_subset <- df_a %>% select(DTXSID, HIT.CALL)

# Join HIT.CALL from File A to File B based on DTXSID
df_b_updated <- df_b %>% left_join(df_a_subset, by = "DTXSID")

# Save the updated file B
write_xlsx(df_b_updated, "C:/Users/91962/Predicting Chemicals/Validation1.xlsx")

cat(" HIT.CALL column added to File B and saved successfully.\n")
