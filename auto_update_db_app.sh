#!/bin/bash

# Load modules
ml simr
ml R/4.3.1

# Navigate to app directory
cd ~/ShinyApps/protDB_shiny

# Run the R script to update the databases and to publish the app to shinyapps.io
Rscript ./scripts/update_databases.R

# Git operations to push updated database Rds objects to remote
git add .
git commit -m "Automated commit with results from script.R"
git push origin main