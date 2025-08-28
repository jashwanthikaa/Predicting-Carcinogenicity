#Preprocessing Removing N/A and just keeping oraganic atoms (from molecular formula) 

library(readxl)
library(writexl)
library(dplyr)
library(stringr)
library(tidyr)

# Load your Excel file
df <- read.csv("C:/Users/91962/Predicting Chemicals/df2_final.csv")

# Remove rows with any NA values
df <- df %>% drop_na()

# Define allowed atoms
allowed_atoms <- c("C", "H", "O", "N", "S")

# Function to check MOLECULAR_FORMULA for allowed atoms
only_allowed_atoms <- function(formula) {
  atoms <- str_extract_all(formula, "[A-Z][a-z]*")[[1]]
  all(atoms %in% allowed_atoms)
}

# Filter rows based on allowed atoms in molecular formula
df <- df %>% filter(sapply(MOLECULAR_FORMULA, only_allowed_atoms))

# Remove columns that are named like atom symbols but are NOT C, H, O, N, S
atom_col_pattern <- "^[A-Z][a-z]?$"
atom_cols <- names(df)[str_detect(names(df), atom_col_pattern)]
drop_cols <- setdiff(atom_cols, allowed_atoms)
df <- df %>% select(-all_of(drop_cols))

# Save cleaned file
write.csv(df, "C:/Users/91962/Predicting Chemicals//Filtered_compounds_cleaned.csv")

cat(" Done! Cleaned file saved to Filtered_compounds_cleaned.csv\n")
