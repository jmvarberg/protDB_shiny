#github action to run database refresh and publish to shinyapps.io
library(janitor)
library(dplyr)
library(tidyr)
library(data.table)
library(curl)
library(rsconnect)
source("./scripts/NCBI_FTP_reference_table.R")
source("./scripts/UniProt_FTP_reference_table.R")

#push to shinyapps.io
rsconnect::deployApp(appName = "protoDB_dev", appFiles = c("app.r", "./databases/ncbi_reference_table.fst", "./databases/uniprot_ids_reference_table.fst", "README.md"), account = "jmvarberg", server = "shinyapps.io") 
