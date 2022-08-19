#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --time=7:00:00
#SBATCH --job-name=exocp
#SBATCH --mem=100GB
#SBATCH --ntasks=1
#SBATCH --error=error.log
#SBATCH --output=output.log

source activate ocp-models
python ocp_extract_10_100k.py