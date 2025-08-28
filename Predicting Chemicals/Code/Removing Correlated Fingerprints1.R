
library(readxl)
library(writexl)
library(caret)

df <- read_excel("C:/Users/91962/Predicting Chemicals/All_319/Filtered_compounds_with_Set1.xlsx")

#Remove non-fingerprint columns
fingerprints <- subset(df, select = -c(
  DTXSID, CASRN, Name, MOLECULAR_FORMULA, 
  MONOISOTOPIC.MASS, AVERAGE_MASS, SMILES, 
  ToxCast.Active, ToxCast.Total, X..ToxCast.Active
))

#Convert all fingerprint columns to numeric
fingerprints[] <- lapply(fingerprints, function(x) as.numeric(as.character(x)))

#Remove columns with any missing values (NA)
fingerprints <- fingerprints[, colSums(is.na(fingerprints)) == 0]

#Print original number of fingerprint features
original_total <- ncol(fingerprints)
cat("Original number of fingerprint features (after NA removal):", original_total, "\n")

#Remove near-zero variance predictors
nzv <- nearZeroVar(fingerprints)
filtereddf <- fingerprints[, -nzv]

#Report number of features after NZV removal
after_nzv_total <- ncol(filtereddf)
cat("Number of features after NZV removal:", after_nzv_total, "\n")
cat("Total NZV features removed:", length(nzv), "\n")

#Create correlation matrix
descrCor <- cor(filtereddf, use = "pairwise.complete.obs")

#Try multiple correlation cutoffs and save results
cutoffs <- c(0.90, 0.95, 0.96, 0.97, 0.98)

for (cutoff in cutoffs) {
  highlyCorDescr <- findCorrelation(descrCor, cutoff = cutoff)
  filteredcor <- filtereddf[, -highlyCorDescr]
  
  #summary for this cutoff
  cat("\n===== Correlation Cutoff:", cutoff, "=====\n")
  cat("Number of highly correlated features removed:", length(highlyCorDescr), "\n")
  cat("Number of features remaining:", ncol(filteredcor), "\n")
  
  output_filename <- paste0("C:/Users/91962/All_319/filtered_fingerprints_cutoff_", cutoff, ".xlsx")
  write_xlsx(filteredcor, output_filename)
}

