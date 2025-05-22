#!/bin/bash
#SBATCH --job-name=job_name               # Job name (set this)
#SBATCH --output=logs/job_%j.out          # STDOUT log file
#SBATCH --error=logs/job_%j.err           # STDERR log file
#SBATCH --ntasks=1                        # Number of tasks
#SBATCH --cpus-per-task=4                 # Cores per task (adjust as needed)
#SBATCH --mem=16G                         # Memory (adjust as needed)
#SBATCH --time=02:00:00                   # Time limit (HH:MM:SS)
#SBATCH --partition=standard              # Partition name (adjust to your system)

# Load Conda (adjust if using a module system instead)
# Uncomment the appropriate one for your system

# Option 1: Conda via system path
# source /apps/ohpc/pub/apps/conda/3/etc/profile.d/conda.sh

# Option 2: Conda via module (example)
# module load anaconda/2023.03

# Activate r-env environment
conda activate r-env

# Optional sanity check
echo "Running R version:"
which Rscript
Rscript --version

echo "Running from: $(pwd)"
ls -l

# Run the R script (pass as variable or relative path)
R_SCRIPT="standard-workflow/main_workflow.r"
Rscript -e "options(future.globals.maxSize=200*1024^3); source('$R_SCRIPT')"
