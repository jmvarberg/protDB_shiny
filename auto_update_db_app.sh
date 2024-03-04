#!/bin/bash
#SBATCH --job-name=protDB_generator
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH --time=1:00:00

# Load modules
ml simr
ml R/4.2.3

# Navigate to app directory
cd ~/ShinyApps/protDB_shiny

# Run the R script to update the databases and to publish the app to shinyapps.io
Rscript ./scripts/update_databases.R

# Git operations to push updated database objects to remote GitHub repository (https://github.com/jmvarberg/protDB_shiny)
git add .
git commit -m "Automated monthly database update"
git push origin main