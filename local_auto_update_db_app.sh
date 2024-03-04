#!/bin/bash

# Run the R script to update the databases and to publish the app to shinyapps.io
Rscript ./scripts/update_databases.R

# Git operations to push updated database objects to remote GitHub repository (https://github.com/jmvarberg/protDB_shiny)
git add .
git commit -m "Automated monthly database update"
git push origin main