setwd("C:/Users/91962/Predicting Chemicals/LCSB/")

library(readr)

getcas <- function(txtfile){
  readfile <- read_lines(txtfile)
  cas = readfile[grep('CAS.*', readfile)]
  print(cas)
}


filelist = list.files(pattern = ".*.txt", full.names = TRUE)

caslist = getcas(filelist)

df <- data.frame(caslist)

library(writexl)
write_xlsx(df, "Cas numbers.xlsx")
