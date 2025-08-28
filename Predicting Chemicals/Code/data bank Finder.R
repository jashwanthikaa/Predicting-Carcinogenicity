#Data Bank Finder
casrn_file <- "C:/Users/91962/Predicting-Toxicity/Both AB (1).txt"
source_folder <- "C:/Users/91962/Predicting-Toxicity/LCSB"
output_folder <- "C:/Users/91962/Predicting-Toxicity/Detail Data"

if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

casrns <- trimws(readLines(casrn_file, warn = FALSE))
casrns <- casrns[casrns != ""]  

files <- list.files(source_folder, pattern = "\\.txt$", full.names = TRUE)

matched_files <- character()

suppressWarnings({
  for (file_path in files) {
    content <- tryCatch(readLines(file_path, warn = FALSE), error = function(e) NULL)
    
    if (!is.null(content)) {
      for (cas in casrns) {
        if (any(grepl(paste0("\\b", cas, "\\b"), content))) {
          copied <- file.copy(file_path, file.path(output_folder, basename(file_path)))
          if (copied) {
            matched_files <- c(matched_files, basename(file_path))
          }
          break  
        }
      }
    }
  }
})