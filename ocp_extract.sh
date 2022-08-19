#!/bin/bash
#SBATCH --nodes=1
#SBATCH --error=4c_error.log
#SBATCH --output=4c_output.log
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --job-name=exocp
#SBATCH --mem=100G
#SBATCH --partition=short

source activate ocp-models
python ocp_extract.py