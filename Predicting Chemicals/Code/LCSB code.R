casrn_file <- "C:/Users/91962/Both AB.txt"
source_folder <- "C:/Users/91962/LCSB/"
output_folder <- "C:/Users/91962/Detail"


# Read CASRNs from the text file
casrns <- readLines(casrn_file)
casrns <- trimws(casrns)
casrns <- casrns[casrns != ""]  # Remove empty lines

# List all files in source folder
files <- list.files(source_folder, full.names = TRUE)

# Loop through each file and check if it contains any CASRN
for (file_path in files) {
  file_content <- tryCatch(readLines(file_path, warn = FALSE), error = function(e) return(NULL))
  
  if (!is.null(file_content)) {
    # Check if any CASRN is present in the content
    found <- any(sapply(casrns, function(casrn) any(grepl(casrn, file_content, fixed = TRUE))))
    
    if (found) {
      file.copy(file_path, file.path(output_folder, basename(file_path)))
      cat("Copied:", basename(file_path), "\n")
    }
  }
}
