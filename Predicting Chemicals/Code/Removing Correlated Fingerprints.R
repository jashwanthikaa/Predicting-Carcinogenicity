
library(readxl)
library(writexl)
library(caret)

df <- read_excel("C:/Users/91962/Predicting Chemicals/All_319/Filtered_compounds_with_Set1.xlsx")

#Remove non-fingerprint (non-numeric) columns
fingerprints <- subset(df, select = -c(
  DTXSID, CASRN, Name, MOLECULAR_FORMULA, 
  MONOISOTOPIC.MASS, AVERAGE_MASS, SMILES, 
  ToxCast.Active, ToxCast.Total, X..ToxCast.Active, Set
))

#Convert all columns to numeric
fingerprints[] <- lapply(fingerprints, function(x) as.numeric(as.character(x)))

#Remove columns with NA values
fingerprints <- fingerprints[, colSums(is.na(fingerprints)) == 0]

# Report initial feature count
cat("Number of features before NZV removal:", ncol(fingerprints), "\n")

# Remove near-zero variance features
nzv <- nearZeroVar(fingerprints)
filtereddf <- fingerprints[, -nzv]

#Report after NZV
cat("Number of features after NZV removal:", ncol(filtereddf), "\n")
cat("NZV features removed:", length(nzv), "\n")

#Compute correlation matrix
descrCor <- cor(filtereddf, use = "pairwise.complete.obs")

#Loop over correlation cutoffs and report retained features
cutoffs <- c(0.90, 0.95, 0.96, 0.97, 0.98)

for (cutoff in cutoffs) {
  highlyCorDescr <- findCorrelation(descrCor, cutoff = cutoff)
  filteredcor <- filtereddf[, -highlyCorDescr]
  
  cat("\n=== Correlation cutoff:", cutoff, "===\n")
  cat("Highly correlated features removed:", length(highlyCorDescr), "\n")
  cat("Remaining features:", ncol(filteredcor), "\n")
  
  file_name <- paste0("C:/Users/91962/All_319/filtered_fingerprints_cutoff_", cutoff, ".xlsx")
  write_xlsx(filteredcor, file_name)
}
