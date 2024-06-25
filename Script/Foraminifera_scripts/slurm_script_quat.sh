#!/bin/bash

# Account name:
#SBATCH --account=nn8105k

# Project name:
#SBATCH --job-name=Hodell


# Time limit:
#SBATCH --time=7-00:0:0

#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=52


module --quiet purge  # Reset the modules to the system default
module load R/4.2.1-foss-2022a
module list


R --vanilla -f "$SCRIPT"

