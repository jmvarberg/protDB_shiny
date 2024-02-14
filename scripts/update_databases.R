#github action to run database refresh and publish to shinyapps.io
library(janitor)
library(dplyr)
library(tidyr)
library(data.table)
library(curl)
source("./scripts/NCBI_FTP_reference_table.R")
source("./scripts/UniProt_FTP_reference_table.R")

