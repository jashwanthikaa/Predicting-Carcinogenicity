#Chemical Compound Separation

# Load isotope data
data(isotopes)

# SET FOLDERS 
input_dir <- "C:/Users/91962/Predicting Chemicals/All_319/"
output_dir <- "C:/Users/91962/All_319/Compound_Details"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

files <- list.files(input_dir, pattern = "\\.txt$", full.names = TRUE)
compound_map <- list()

# READ FILES AND GROUP BY COMPOUND
for (file in files) {
  lines <- readLines(file)
  
  # Extract core details
  compound_line <- grep("^CH\\$NAME:", lines, value = TRUE)[1]
  compound <- str_trim(str_remove(compound_line, "CH\\$NAME:"))
  
  ion_line <- grep("^MS\\$FOCUSED_ION: PRECURSOR_TYPE", lines, value = TRUE)
  ionization <- str_trim(str_remove(ion_line, "MS\\$FOCUSED_ION: PRECURSOR_TYPE"))
  
  pm_line <- grep("^MS\\$FOCUSED_ION: PRECURSOR_M/Z", lines, value = TRUE)
  parentmass <- str_trim(str_remove(pm_line, "MS\\$FOCUSED_ION: PRECURSOR_M/Z"))
  
  formula_line <- grep("^CH\\$FORMULA:", lines, value = TRUE)
  formula <- if (length(formula_line) > 0) str_trim(str_remove(formula_line, "CH\\$FORMULA:")) else NA
  
  ce_line <- grep("^AC\\$MASS_SPECTROMETRY: COLLISION_ENERGY", lines, value = TRUE)
  ce_value <- str_extract(ce_line, "\\d+")
  ce_title <- paste0("Collision ", ce_value)
  
  # Extract peak data
  peak_start <- grep("^PK\\$PEAK: m/z int. rel.int.", lines)
  if (length(peak_start) == 0) next
  
  peak_lines <- lines[(peak_start + 1):length(lines)]
  peak_lines <- peak_lines[peak_lines != "" & !grepl("^//", peak_lines)]
  
  key <- paste(compound, ionization, parentmass, sep = "||")
  
  if (!key %in% names(compound_map)) {
    # Generate MS1 using enviPat
    ms1_block <- c()
    if (!is.na(formula)) {
      try({
        formula_list <- list(formula)
        result <- isopattern(isotopes, formula_list, threshold = 0.1, charge = 0, emass = 0.00054858, algo = 1)
        ms1_block <- apply(result[[1]], 1, function(row) paste0(sprintf("%.5f", as.numeric(row[1])), " ", round(as.numeric(row[2]), 1)))
      }, silent = TRUE)
    }
    
    compound_map[[key]] <- list(
      compound = compound,
      ionization = ionization,
      parentmass = parentmass,
      formula = formula,
      ms1 = ms1_block,
      collisions = list()
    )
  }
  
  compound_map[[key]]$collisions[[ce_title]] <- peak_lines
}

# WRITE OUTPUT FILES (.ms) WITH CLEANED FORMATTING
for (entry in compound_map) {
  compound <- entry$compound
  ionization <- entry$ionization
  parentmass <- entry$parentmass
  ms1_block <- entry$ms1
  collisions <- entry$collisions
  
  output_lines <- c(
    paste0("Compound: ", compound),
    paste0("Ionization: ", ionization),
    paste0("Parent Mass: ", parentmass),
    ""
  )
  
  # Clean MS1 block
  if (length(ms1_block) > 0) {
    ms1_cleaned <- lapply(ms1_block, function(line) {
      parts <- str_split(str_trim(line), "\\s+")[[1]]
      if (length(parts) >= 2) paste(parts[1], parts[2]) else NULL
    })
    output_lines <- c(output_lines, "MS1", unlist(ms1_cleaned), "")
  }
  
  # Clean collision peaks (remove 3rd value)
  for (ce in sort(names(collisions))) {
    output_lines <- c(output_lines, ce)
    
    cleaned_peaks <- lapply(collisions[[ce]], function(line) {
      parts <- str_split(str_trim(line), "\\s+")[[1]]
      if (length(parts) >= 2) paste(parts[1], parts[2]) else NULL
    })
    
    output_lines <- c(output_lines, unlist(cleaned_peaks), "")
  }
  
  # Save as .ms file
  filename <- paste0(gsub("[^a-zA-Z0-9_]", "_", compound), ".ms")
  writeLines(output_lines, con = file.path(output_dir, filename))
}

cat("All .ms files written successfully to:", output_dir, "\n")

--------------------------------------------------------------------

# Store compound to accession + formula mapping
accession_map <- list()

for (file in files) {
  lines <- readLines(file)
  
  # Compound name
  compound_line <- grep("^CH\\$NAME:", lines, value = TRUE)[1]
  compound <- str_trim(str_remove(compound_line, "CH\\$NAME:"))
  
  # Accession
  accession_line <- grep("^ACCESSION:", lines, value = TRUE)
  accession <- if (length(accession_line) > 0) str_trim(str_remove(accession_line, "ACCESSION:")) else "NA"
  
  # Chemical formula
  formula_line <- grep("^CH\\$FORMULA:", lines, value = TRUE)
  formula <- if (length(formula_line) > 0) str_trim(str_remove(formula_line, "CH\\$FORMULA:")) else "NA"
  
  # Store as a list per compound
  accession_map[[compound]] <- list(
    accession = accession,
    formula = formula
  )
}

# Write to text file (Name: Accession - Formula)
accession_output <- sapply(names(accession_map), function(name) {
  acc <- accession_map[[name]]$accession
  form <- accession_map[[name]]$formula
  paste0(name, ": ", acc, " | Formula: ", form)
})

# Save to file
writeLines(accession_output, file.path(output_dir, "accession_lists.txt"))



--------------------------------------------------------------------------

  library(writexl)

# List of combined compound names
combined_compounds <- sapply(compound_map, function(x) x$compound)

# Create empty data frame with desired columns
combined_df <- data.frame(
  Compound = character(),
  CASRN = character(),
  Formula = character(),
  Mass = character(),
  stringsAsFactors = FALSE
)

# Loop through files and match only the combined compounds
for (file in files) {
  lines <- readLines(file)
  
  # Get compound name
  name_line <- grep("^CH\\$NAME:", lines, value = TRUE)[1]
  compound <- str_trim(str_remove(name_line, "CH\\$NAME:"))
  
  # Only include if it's in the combined list
  if (!(compound %in% combined_compounds)) next
  
  # CASRN
  cas_line <- grep("^CH\\$LINK: CAS", lines, value = TRUE)
  casrn <- if (length(cas_line) > 0) str_trim(str_remove(cas_line, "CH\\$LINK: CAS")) else NA
  
  # Formula
  formula_line <- grep("^CH\\$FORMULA:", lines, value = TRUE)
  formula <- if (length(formula_line) > 0) str_trim(str_remove(formula_line, "CH\\$FORMULA:")) else NA
  
  # Mass
  mass_line <- grep("^CH\\$EXACT_MASS:", lines, value = TRUE)
  mass <- if (length(mass_line) > 0) str_trim(str_remove(mass_line, "CH\\$EXACT_MASS:")) else NA
  
  # Add to dataframe
  combined_df <- rbind(combined_df, data.frame(
    Compound = compound,
    CASRN = casrn,
    Formula = formula,
    Mass = mass,
    stringsAsFactors = FALSE
  ))
}

# Remove duplicates (based on Compound column)
combined_df <- combined_df[!duplicated(combined_df$Compound), ]

# Ensure output directory exists
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Save Excel output
write_xlsx(combined_df, file.path(output_dir, "chemical_compounds.xlsx"))

cat("Saved combined compound list (deduplicated) to Excel.\n")
