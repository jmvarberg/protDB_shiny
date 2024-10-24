#github action to run database refresh and publish to shinyapps.io
options(repos = BiocManager::repositories())
library(janitor)
library(dplyr)
library(tidyr)
library(data.table)
library(curl)
library(rsconnect)
source("./scripts/NCBI_FTP_reference_table.R")
source("./scripts/UniProt_FTP_reference_table.R")

#push to shinyapps.io
rsconnect::deployApp(appName = "protDB_generator", 
                     account = "jmvarberg", 
                     server = "shinyapps.io") 
